#!/bin/bash

JobCMD="$0 $@ ### `date` `hostname` `pwd`"
echo "### $0 $@" ; echo ""
FirstDir=$PWD

MyCFG=""
MyDataset=""
MyListFile=""
MyMaxFiles="1"
MyJSON=""
MyCMSARG=""
MyURL=""
MyAddOutFiles=""
myCondorName=""
MyCMDFiles=""
AutoSubmit="0"
HostSubmit=""
MyCondorOutSRM="NULL"
ReSubmitStrings=""
AddSRCs=""
MyCheckDir=""
MyXMLDirs=""

MySRMHOST=root://cms-xrdr.sdfarm.kr:1094
MySRMUSERDIR=/xrd/store/user/$USER

function dirString { dataStr=${1//\//_D_}; echo ${dataStr/_D_/}; }

function stringInFile { numDup=0; 
   for str in `cat $1`; do bname1=`basename $str`; dname1=`dirname $str`; bname2=`basename $2`; dname2=`dirname $2`;
      if [ "${bname1}" == "${bname2}" ] && [ "${dname1}" == "${dname2}" ]; then ((numDup++)); fi
   done
   echo $numDup 
}

function removeDupLineInFile {
   file=$1 ; tmpFile="${file}.temp${RANDOM}" ; rm -rf $tmpFile ; touch $tmpFile
   for str in `cat $file`; do
      if [ "`stringInFile $tmpFile $str`" == "0" ]; then 
			echo "$str" >> $tmpFile
      else
         echo " @ DupLine $str, ignored it"
      fi
   done
   mv -f $tmpFile $file
}

function PrintSet() {
	echo "##################"
	echo "@python $MyCFG"
	if [ $MyCMSARG ]; then echo "@cmsRunARG $MyCMSARG"; fi
	if [ $MyDataset ]; then echo "@dataset $MyDataset"; fi
	if [ $MyListFile ]; then echo "@list $MyListFile"; fi
	if [ $MyJSON ]; then echo @"json $MyJSON"; fi
	echo "@MaxFiles $MyMaxFiles"
	echo "@CondorName $myCondorName"
	echo "@CurrentDir $PWD"
	echo "@CondorWork $myCondorWorkDir"
	echo "##################"
	echo ""
}

function PrintUsage() {
	echo "usage $0 -help --> PrintUsage" 
	echo "usageToSubmit $0 -py=cmsRun.py -dataset=DATASET" 
	echo "usageToSubmit $0 -py=cmsRun.py -list=listfile"
	echo "usageToSubmit $0 -py=cmsRun.py -dataset=DATASET -nfiles=1 -json=JSONfile -output=hist.root  -arg=\"isMC=True isData=False\" -url=root://cms-xrdr.sdfarm.kr:1094//xrd/"
	echo "usageToCheck  $0 -check=condorRunDirectory"
	echo "usageToCheck  $0 -makejson=condorRunDirectory"
	exit
}

function argvPar {
	if [ ! $1 ]; then PrintUsage; fi
	doCheck=0
   for argv in $@
   do
      arg1=`echo $argv | cut -d"=" -f 1`
      arg2=`echo $argv | cut -d"=" -f 2-`
		if [ "${arg1}" == "-h" ]; then PrintUsage; fi
		if [ "${arg1}" == "-help" ]; then PrintUsage; fi
      if [ "${arg1}" == "-py"      ]; then MyCFG=${arg2}      ; continue; fi
      if [ "${arg1}" == "-list"    ]; then MyListFile=${arg2} ; continue; fi
      if [ "${arg1}" == "-json"    ]; then MyJSON=${arg2}     ; continue; fi
      if [ "${arg1}" == "-nfiles"  ]; then MyMaxFiles=${arg2} ; continue; fi
      if [ "${arg1}" == "-dataset" ]; then MyDataset=${arg2}  ; continue; fi
      if [ "${arg1}" == "-arg"     ]; then MyCMSARG=${arg2}   ; continue; fi
      if [ "${arg1}" == "-url"     ]; then MyURL=${arg2}      ; continue; fi
      if [ "${arg1}" == "-output"  ]; then MyAddOutFiles="$MyAddOutFiles ${arg2}"      ; continue; fi
      if [ "${arg1}" == "-submit"  ]; then AutoSubmit="1"; HostSubmit=${arg2}    ; continue; fi
      if [ "${arg1}" == "-srmdir"  ]; then MyCondorOutSRM="${arg2}" ; continue; fi
      if [ "${arg1}" == "-CMDFiles"  ]; then MyCMDFiles="${arg2}"     ; continue; fi
      if [ "${arg1}" == "-addSRC"  ]; then AddSRCs="${AddSRCs} ${arg2}" ; continue; fi
      if [ "${arg1}" == "-check"   ]; then doCheck=1; MyCheckDir="${arg2}" ; continue; fi
      if [ "${arg1}" == "-makejson"  ]; then doMakeJson=1; MyXMLDirs="${MyXMLDirs} ${arg2}" ; continue; fi
      echo -e "\e[31mUnknownOption $argv Ignored it! \e[0m" 
   done
	if [ "${doCheck}" == "1" ]; then
		if [ "${MyCheckDir}" == "-check" ]; then 
			MyCheckDir=`find . -maxdepth 1 -type d -name "condor_*" | xargs -n 1 basename`
		fi
	fi

	if [ "`echo ${MyXMLDirs} | grep makejson | wc -l`" == "1" ]; then 
		MyXMLDirs=`find . -maxdepth 1 -type d -name "condor_*" | xargs -n 1 basename`
	fi

	if [ ! -f $MyCFG ]; then echo "Error NotFound $MyCFG"; exit; fi
	if [ $MyDataset ]; then 
		if [ $MyListFile ]; then
			echo "$MyDataset or $MyListFile"; exit;
		fi
	fi
	if [ $MyListFile ] && [ ! -f $MyListFile ]; then echo "Error NotFound $MyListFile"; exit; fi
	if [ $MyJSON ] && [ ! -f $MyJSON ]; then echo "Error NotFound $MyJSON"; exit; fi

	DataString=`dirString ${MyDataset}${MyListFile}`
	myCondorName="condor_${DataString}_`date +"%y%m%d_%H%M%S"`"
	myCondorWorkDir=${FirstDir}/${myCondorName}
}

function MakeDataList {
   echo "# Make List for ${MyDataset}${MyListFile}"
	export SSL_CERT_DIR='/etc/grid-security/certificates'
   if [ $MyListFile ]; then
      ThisGURL=""
      cp `readlink -f $MyListFile` ${myCondorWorkDir}/.inputfiles.das
      cp ${myCondorWorkDir}/.inputfiles.das ${myCondorWorkDir}/inputfiles.das
      chmod -w ${myCondorWorkDir}/.inputfiles.das
      echo -n "   @ "
      echo "`cat ${myCondorWorkDir}/inputfiles.das | wc -l` files"
      echo -n "   @ "
      head -n 1 ${myCondorWorkDir}/inputfiles.das
      if [ "`cat ${myCondorWorkDir}/inputfiles.das | wc -l`" == "0" ]; then echo "Not Found Files in $Dataset" >> $myCondorWorkDir/.jobConfigError; fi
   else
      if [ "${MyDataset:(-4)}" == "USER" ]; then instance="instance=prod/phys03"; fi
		DataString=`dirString ${MyDataset}`
      das_client --query="file dataset=${MyDataset} ${instance} system=dbs3 | grep file.name, file.nevents" --limit=0 | grep ".root" | sed "s/[\[\'\"\,]//g" >  ${myCondorWorkDir}/.inputfiles.das
      sites=`das_client --query="site dataset=${MyDataset} ${instance}" --limit=0`
      #sites=${sites//'N/A'/}
      sitesStr=""
      for site in $sites
		do 
			if [ "${site:0:2}" == "T1" ] || [ "${site:0:2}" == "T2" ] || [ "${site:0:2}" == "T3" ]; then sitesStr="$sitesStr $site"; fi
		done
      isKISTI=`echo $sites | grep -E "T3_KR_KISTI|cms-se.sdfarm.kr" | wc -l`
      URL="root://cms-xrd-global.cern.ch/"
      if [ "${isKISTI}" != "0" ]; then URL="root://cms-xrdr.sdfarm.kr:1094//xrd/"; fi
      ThisGURL=$URL
		if [ ! $MyURL ]; then MyURL=${ThisGURL} ; fi
      cp -r ${myCondorWorkDir}/.inputfiles.das ${myCondorWorkDir}/inputfiles.das
      chmod -w ${myCondorWorkDir}/.inputfiles.das
      echo -n "   @ " >> ${myCondorWorkDir}/.das${DataString}.summary
      echo "`cat ${myCondorWorkDir}/inputfiles.das | wc -l` files `cat ${myCondorWorkDir}/inputfiles.das | awk '{NTotalEvents += $2} END{print NTotalEvents}'` events in $sitesStr URL=$URL" >> ${myCondorWorkDir}/.das${DataString}.summary
      echo -n "   @ " >> ${myCondorWorkDir}/.das${DataString}.summary
      head -n 1 ${myCondorWorkDir}/inputfiles.das >> ${myCondorWorkDir}/.das${DataString}.summary
      if [ "`cat ${myCondorWorkDir}/inputfiles.das | wc -l`" == "0" ]; then echo "Not Found Files in $Dataset" >> $myCondorWorkDir/.jobConfigError; fi
		echo "   @ SITES: $sitesStr" >> ${myCondorWorkDir}/.das${DataString}.summary
		cat ${myCondorWorkDir}/.das${DataString}.summary
		cat ${myCondorWorkDir}/inputfiles.das >> ${myCondorWorkDir}/.das${DataString}.summary
   fi
   if [ "`cat ${myCondorWorkDir}/inputfiles.das | wc -l`" == "0" ]; then echo "   Error Not Found inputfiles in ${Dataset} !!! Goto exit"; exit; fi
}

function MakeCMSSW {
	echo "### Make PSet for Condor"
	CondorSrc="$myCondorWorkDir/input/src"
	if [ ! -d $CondorSrc ]; then mkdir -p $CondorSrc ; fi
	MakePSetFile="${CondorSrc}/makePSet"
	CMSSWPKL="${CondorSrc}/CMSSW.pkl"
	CMSSWCFG="${CondorSrc}/CMSSW.py"
	CMSSWOUTFILE="${CondorSrc}/CMSSW.outname"
#from FWCore.ParameterSet.VarParsing import VarParsing
cat << EOF >  ${MakePSetFile}
#!/usr/bin/env python
import FWCore.ParameterSet.Config as cms
from optparse import OptionParser
import sys, os, pickle, imp
parser = OptionParser()
(options,args) = parser.parse_args()
filename = args[0]
handle = open(filename,'r')
cfo = imp.load_source("pycfg",filename, handle)
cmsProcess = cfo.process
if hasattr(cmsProcess, 'options'):
	cmsProcess.options.wantSummary = cms.untracked.bool(True)
else:
	cmsProcess.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
pklFileName = "${CMSSWPKL}"
pklFile = open(pklFileName,"wb")
pickle.dump(cmsProcess, pklFile)
pklFile.close()
outNameFile = open('${CMSSWOUTFILE}','w')
if hasattr(cmsProcess, 'TFileService'):
	outNameFile.write(cmsProcess.TFileService.fileName.value()+" \\n")
for modName in cmsProcess.outputModules_():
	outNameFile.write(getattr(cmsProcess, modName).fileName.value()+" \\n")
outNameFile.close()
EOF
	python ${MakePSetFile} ${MyCFG} ${MyCMSARG} >& ${MakePSetFile}.log
	ThisJobStatus=$?
	if [ "$ThisJobStatus" != "0" ]; then echo "   @ Error MakePSet "; cat ${MakePSetFile}.log; exit; fi
#	rm -rf $MakePSetFile ${MakePSetFile}.log
cat << EOF > $CMSSWCFG
import FWCore.ParameterSet.Config as cms
import pickle, os
process = pickle.load(open('CMSSW.pkl', 'rb'))
#process.options.wantSummary = cms.untracked.bool(True) # make some problem for unschedules 
f = open('ThisInputFiles.list','r')
tempInFiles = f.readlines()
f.close()
InFiles = [(InFile.replace('\\n','')) for InFile in tempInFiles]
process.source.fileNames = InFiles
for file in InFiles:
   print 'cmsRunInFile:', file
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
EOF
	if [ "$MyJSON" != "" ]; then
		bJSON=`basename $MyJSON`
		cp -r `readlink -f $MyJSON` ${CondorSrc}/$bJSON
cat << EOF >> $CMSSWCFG
import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '${bJSON}').getVLuminosityBlockRange()
EOF
	fi
   if [ "${MyAddOutFiles}" != "" ]; then
      for addout in ${MyAddOutFiles//,/ }; do echo "$addout" >> ${CMSSWOUTFILE}; done
   fi
   removeDupLineInFile ${CMSSWOUTFILE}
   numOutFile=`cat ${CMSSWOUTFILE} | wc -l`
   echo -n "   @ Found OutFileName $numOutFile " 
   for out_name_file in `cat ${CMSSWOUTFILE}`; do OutROOTNames="${out_name_file} ${OutROOTNames}"; done
   echo $OutROOTNames

   echo "### MakeInput"
   cd $myCondorWorkDir/input
	if [ "${MyCMDFiles}" != "" ]; then 
		for MyCMDFile in ${MyCMDFiles//,/ }
		do
			echo "### ADD $MyCMDFile in $CMSSW_BASE/src/" 
			cp -r `readlink -f $FirstDir/$MyCMDFile` ${CondorSrc}
		done
	fi
   ln -s $CMSSW_BASE/lib .
   ln -s $CMSSW_BASE/biglib .
	ln -s $CMSSW_BASE/external .
###################
	cDir=$PWD
	inDataDirs=`find $CMSSW_BASE/src -type d -name "data"`
	for inDataDir in $inDataDirs
	do
   	size=`du $inDataDir | tail -n 1 | awk '{print $1}'`
   	if [ $size -gt 1000000 ]; then echo "### $inDataDir is bigger than 1GB"; fi
		echo -n "   @cpying data dir "
		du -h $inDataDir | tail -n 1
	   nameDataDir="${inDataDir/"${CMSSW_BASE}/src/"/}"
		nameDataDir=`dirname $nameDataDir`
		cd src
		mkdir -p $nameDataDir
		cd $nameDataDir
		ln -s $inDataDir .
		cd $cDir
	done
	if [ "${AddSRCs}" != "" ]; then
		cd src
		for AddSRC in $AddSRCs; do ln -s $CMSSW_BASE/src/$AddSRC .; done
		cd ../
	fi
##################
   tar zcfh ../input.tgz *
   rm -rf ${myCondorWorkDir}/input
   cd $FirstDir
}

function MakeJob {
	echo "$JobCMD" > ${myCondorWorkDir}/.jobCMD
   echo "# Make Condor Config "
   mkdir -p ${myCondorWorkDir}/condorLog
   totalNFile=`cat $myCondorWorkDir/inputfiles.das | wc -l`
   totalNFileT=`expr $totalNFile + $MyMaxFiles - 1`
   nJob=`expr $totalNFileT / $MyMaxFiles`
   echo "   @ NJobs $nJob for $totalNFile / $MyMaxFiles"
	SRMOutDir="NULL"
	if [ "${MyCondorOutSRM}" != "NULL" ]; then SRMOutDir="${MyCondorOutSRM}/`basename ${myCondorWorkDir}`"; fi
cat << EOF > ${myCondorWorkDir}/job.jdl
executable = \$ENV(PWD)/runCMSRUN.sh
universe = vanilla
output   = condorLog/condorLog_\$(Cluster)_\$(Process).log
error    = condorLog/condorLog_\$(Cluster)_\$(Process).log
log      = /dev/null
use_x509userproxy = True
should_transfer_files = yes
initialdir = \$ENV(PWD)
transfer_input_files = input.tgz, inputfiles.das
when_to_transfer_output = ON_EXIT
transfer_output_files = condorOut, cmsRunLog 
environment = "CONDOR_ID=\$(Cluster)_\$(Process)  SRMOUTDIR=${SRMOutDir}"
requirements = (machine != "wn3003.sdfarm.kr") && (machine != "wn3004.sdfarm.kr") && (machine != "wn3005.sdfarm.kr")
#requirements = (machine != "wn3011.sdfarm.kr") && (machine != "wn3018.sdfarm.kr")  #condor_status
#arguments = ${MyMaxFiles} root://cms-xrd-global.cern.ch/
#arguments = ${MyMaxFiles} root://cms-xrdr.sdfarm.kr:1094//xrd/
arguments = ${MyMaxFiles} ${MyURL}
queue $nJob 
EOF
	cp ${myCondorWorkDir}/job.jdl ${myCondorWorkDir}/.job.jdl
	tar zchf ${myCondorName}.tgz ${myCondorName}
	tarFile=${myCondorName}/.${myCondorName}.tgz
	mv ${myCondorName}.tgz $tarFile
	if [ -f $myCondorWorkDir/.jobConfigError ]; then
	   echo ""
	   echo -e "### Condor Setup Error "
	   cat $myCondorWorkDir/.jobConfigError
	else
	   echo "ToSubmit \$> cd ${myCondorName}; condor_submit job.jdl; cd -"
	fi
	if [ "$AutoSubmit" == "1" ]; then
		host=`echo "${HostSubmit}:" | cut -f 1 -d":"`
		hostDir=`echo "${HostSubmit}:" | cut -f 2 -d":"`
		echo "SubmitTo $host $hostDir"
		if [ "$host" == "ui02" ] || [ "$host" == "ui02.sdfarm.kr" ]; then
			host=ui02.sdfarm.kr
			if [ "$hostDir" == "" ]; then hostDir=/hcp/data/data02/jelee/rCondor ; fi
			if [ "${hostDir:0:1}" != "/" ]; then hostDir=/hcp/data/data02/jelee/rCondor/${hostDir} ; fi
			bTarFile=`basename ${myCondorName}/.${myCondorName}.tgz`
			name=${bTarFile/./}
			name=${name/.tgz/}
			ssh -p 4280 $host "if [ ! -d $hostDir ]; then mkdir -p $hostDir ; fi"
			isOK=`ssh -p 4280 $host "if [ -d $hostDir ]; then echo OK; fi"`
			if [ "$isOK" == "OK" ]; then
				echo "#RemoteSubmit $tarFile ${host}:${hostDir}/"
				scp -P 4280  $tarFile ${host}:${hostDir}/
				ssh -p 4280 ${host} "cd ${hostDir}/; tar zxvf $bTarFile; mv $bTarFile $name; cd $name; condor_submit job.jdl; exit;" | tee ${myCondorName}/.jobClusterID_ui02
			fi
		elif [ "$host" == "ui10" ] || [ "${host}" == "ui10.sdfarm.kr" ]; then
			host=134.75.124.121
			if [ "$hostDir" == "" ]; then hostDir=/cms/home/jelee/scratch/rCondor ; fi
			if [ "${hostDir:0:1}" != "/" ]; then hostDir=/cms/home/jelee/scratch/rCondor/${hostDir} ; fi
			bTarFile=`basename ${myCondorName}/.${myCondorName}.tgz`
			name=${bTarFile/./}
			name=${name/.tgz/}
			ssh -p 4280 $host "if [ ! -d $hostDir ]; then mkdir -p $hostDir ; fi"
			isOK=`ssh -p 4280 $host "if [ -d $hostDir ]; then echo OK; fi"`
			if [ "$isOK" == "OK" ]; then
				echo "#RemoteSubmit $tarFile ${host}:${hostDir}/"
				scp -P 4280 $tarFile ${host}:${hostDir}/
				ssh -p 4280 ${host} "cd ${hostDir}/; tar zxvf $bTarFile; mv $bTarFile $name; cd $name; condor_submit job.jdl; exit;" | tee ${myCondorName}/.jobClusterID_ui10
			fi
		elif [ "$HostSubmit" == "-submit" ]; then
			cd ${myCondorName}
			condor_submit job.jdl | tee .jobClusterID_local
			cd -
		else  
			echo "### Error I don't know $HostSubmit"
		fi
	else
		echo "NotSubmitYet" > ${myCondorName}/.jobClusterID_NotYet
	fi
}

function Check() {
dir=$1
if [ -f $dir/DoneCheck ]; then cat $dir/DoneCheck; return ; fi

if [ ! -d $dir ]; then echo "#Error NotFound $dir directory"; exit;  fi

nJob=`grep "^queue"  $dir/.job.jdl  | awk '{print $2}'`
nLog=`ls -al $dir/condorLog/condorLog*.log | wc -l`
nOut=`find $dir/condorOut  -maxdepth 1 -type f | wc -l`

echo "#nJobs $nJob nLogs $nLog nOut $nOut"

idx=0
nOK=0
nNot=0
rm -rf $dir/OKtemp
for logFile in `ls -1vt $dir/condorLog/condorLog*.log`
do
	bName=`basename $logFile`
	bName=${bName/condorLog_/}
	bName=${bName/.log/}
	((idx++))
	ok=`tail -n 1000 $logFile | grep "#FinalCondorRunResult cmsRunStatus 0 DiffnFiles 0 DiffEvtDAS_JSON 0 DiffEvtJSONGOOD_cmsRun 0" | wc -l`
	if [ "$ok" == "1" ]; then
		((nOK++)); 
		tail -n 1000 $logFile | grep "^ThisRunFileOK" | awk '{print $2}' >> $dir/OKtemp
	else 
		((nNot++));
		if [ ! -d ${dir}/fail/condorLog ]; then mkdir -p ${dir}/fail/condorLog ; fi
		if [ ! -d ${dir}/fail/condorOut ]; then mkdir -p ${dir}/fail/condorOut ; fi
		if [ ! -d ${dir}/fail/cmsRunLog ]; then mkdir -p ${dir}/fail/cmsRunLog ; fi
		mv $logFile  ${dir}/fail/condorLog/
		mv $dir/condorOut/*_${bName}.* ${dir}/fail/condorOut/
		mv $dir/cmsRunLog/*_${bName}.* ${dir}/fail/cmsRunLog/
	fi
	if [ "$ok" == "0" ]; then
		echo -e "\e[1;31mIdx $idx/$nJob nOK $nOK nFail $nNot This $ok $logFile \e[0m"
	else
		echo "Idx $idx/$nJob nOK $nOK nFail $nNot This $ok $logFile"
	fi
done

rm -rf $dir/failed_inputfiles.das
rm -rf $dir/succes_inputfiles.das
idx=0
nFiles=`cat $dir/.inputfiles.das | wc -l`
div=`expr $nFiles / 99`
if [ $nFiles -le 100 ]; then div="1"; fi
echo "Checking InputFile"
while read line
do
	((idx++))
#	if [ "`expr $idx % $div`" == "0" ]; then ((per100++)); echo "   @CheckingInputFile $per100 % "; fi
	fileName=`echo $line | awk '{print $1}'`
	isOK=`grep "$fileName" $dir/OKtemp | wc -l`
	if [ "$isOK" == "1" ]; then
		echo $line >> $dir/succes_inputfiles.das
	else
		echo $line >> $dir/failed_inputfiles.das
	fi
done < $dir/.inputfiles.das  

nFileTota=`cat $dir/.inputfiles.das | wc -l`
nEvtsTota=`cat $dir/.inputfiles.das | awk '{sum+=$2} END{print sum}'`
nFileGood=`cat $dir/succes_inputfiles.das | wc -l`
nEvtsGood=`cat $dir/succes_inputfiles.das | awk '{sum+=$2} END{print sum}'`
nFilesBad=0
nEvtssBad=0
if [ -f $dir/failed_inputfiles.das ]; then
	nFilesBad=`cat $dir/failed_inputfiles.das | wc -l`
	nEvtssBad=`cat $dir/failed_inputfiles.das | awk '{sum+=$2} END{print sum}'`
fi

nChFile=`expr $nFileTota - $nFileGood - $nFilesBad`
nChEvts=`expr $nEvtsTota - $nEvtsGood - $nEvtssBad`
echo "#FinalCheck $nChFile $nChEvts nFile: total $nFileTota ok $nFileGood bad $nFilesBad nEvents: total $nEvtsTota ok $nEvtsGood bad $nEvtssBad for $dir" > $dir/check.log
rm -rf $dir/OKtemp 
rm -rf $dir/succes_inputfiles.das

nDAS=`cat $dir/.inputfiles.das | wc -l`
nLog=`ls $dir/condorLog/*.log  | wc -l`
nOut=`ls $dir/condorOut/*.root | wc -l`
nCMS=`ls $dir/cmsRunLog/*.log  | wc -l`
nXML=`ls $dir/cmsRunLog/*.xml  | wc -l`

echo "#CondorFile DAS $nDAS Log $nLog Out $nOut CMS $nCMS XML $nXML xml $dir" >> $dir/check.log
cat $dir/check.log
if [ -f $dir/failed_inputfiles.das ]; then
	mv $dir/failed_inputfiles.das $dir/inputfiles.das
	sed 's/queue/#queue/g' $dir/job.jdl > $dir/job.jdlRe
	sed 's/arguments/#arguments/g' $dir/job.jdlRe > $dir/job.jdlRe1
	rm -rf $dir/job.jdlRe
	grep "^arguments" $dir/job.jdl | awk '{$3="1"; print}' >> $dir/job.jdlRe1
	echo "queue $nFilesBad" >> $dir/job.jdlRe1
	echo "### "
	mv -f $dir/job.jdlRe1 $dir/job.jdl
	tail -n 1 $dir/job.jdl
	echo -e "\e[1;31m cd $dir ; condor_submit job.jdl; cd - # $nFilesBad  ToReSubmit \e[0m"
	ReSubmitStrings="${ReSubmitStrings} cd $dir ; condor_submit job.jdl; cd - # $nFilesBad\n"	
	if [ "$AutoSubmit" == "1" ]; then
		cd $dir ; condor_submit job.jdl; cd -
	fi
fi
nFinalDAS=`cat $dir/.inputfiles.das | wc -l`
nFinalCMS=`tail -n 1000 $dir/condorLog/condorLog_*.log | grep ThisRunFileOK | wc -l`
if [ "$nFinalCMS" == "$nFinalDAS" ]; then 
	echo -e "\e[1;32m#FinalDone $nFinalDAS $nFinalCMS files in $dir \e[0m"; 
	echo "#FinalDone $nFinalDAS $nFinalCMS files in $dir" > $dir/DoneCheck
fi

}

function MakeRun() {
cat << EOF > ${myCondorWorkDir}/runCMSRUN.sh
#!/bin/bash

echo "### Condor Running \`whoami\`@\`hostname\`:\`pwd\` CONDOR_ID=\${CONDOR_ID} \`date\` "
export MyCondorCluster=\`echo \${CONDOR_ID} | cut -d_ -f1\`
export MyCondorISection=\`echo \${CONDOR_ID} | cut -d_ -f2\`
export MyCondorMaxFile=\$1
export MyCondorInputURL="\$2"

echo "MyCondorCluster    \$MyCondorCluster  " 
echo "MyCondorISection   \$MyCondorISection " 
echo "MyCondorMaxFile    \$MyCondorMaxFile  " 
echo "MyCondorInputURL   \$MyCondorInputURL "

function JsonEvents() {
   JSON=\$1
   rootFile=\$2
   xrdURL="root://cms-xrd-global.cern.ch/"
   if [ \$3 ]; then xrdURL=\$3 ; fi
   
   bName=\`basename \$rootFile\`
   bDir=".tempGR_\${bName}"
   mkdir -p \${bDir}
   txt1="\${bDir}/file1"
   txt2="\${bDir}/file2"
   
   echo "edmFileUtil --eventsInLumis \${rootFile} >& \${txt1}"
   edmFileUtil --eventsInLumis \${rootFile} >& \${txt1}
   nTotal=\`grep "runs," \${txt1} | grep root | grep lumis | grep events | awk '{print \$6}'\`
   nLumi=\`grep "runs," \${txt1} | grep root | grep lumis | grep events | awk '{print \$4}'\`
   grep -A \$nLumi "Run           Lumi       # Events" \${txt1} | tail -n \${nLumi} > \${txt2}
   
   numGood=0
   numBad=0
   numTotal=0
   while read CMD; 
   do
      a=\`echo \$CMD | awk '{print \$1}'\`
      b=\`echo \$CMD | awk '{print \$2}'\`
      c=\`echo \$CMD | awk '{print \$3}'\`
      echo "{\"\${a}\":[[\${b},\${b}]]}" > \${bDir}/lumi_\${a}_\${b}.txt
      isOK=\`compareJSON.py --and \${bDir}/lumi_\${a}_\${b}.txt \${JSON} | grep "\$a" | wc -l\`
      numTotal=\`expr \$numTotal + \$c\`
      if [ "\$isOK" == "1" ]; then
         numGood=\`expr \$numGood + \$c\`
      else
         numBad=\`expr \$numBad + \$c\`
      fi
      echo "   GoodLumiCheck \$isOK \$rootFile \$a \$b \$c \$1"
      rm -rf \${bDir}/lumi_\${a}_\${b}.txt
   done < \${txt2}
   
   echo "#edmFileJSONresult \$rootFile Good \$numGood Bad \$numBad Sum \$numTotal == Total \$nTotal \$1"
   rm -rf \$bDir
}

FirstDir=\$PWD
condorOutDir=\${FirstDir}/condorOut
cmsRunLogDir=\${FirstDir}/cmsRunLog
mkdir -p \${condorOutDir}
mkdir -p \${cmsRunLogDir}

echo "@ CMSSW setup"
export SCRAM_ARCH=$SCRAM_ARCH
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source \${VO_CMS_SW_DIR}/cmsset_default.sh
scramv1 project CMSSW $CMSSW_VERSION
cd $CMSSW_VERSION
tar zxf \${FirstDir}/input.tgz
cd src
eval \`scramv1 runtime -sh\`

echo "@ Check InputFiles \$MyCondorISection \$MyCondorMaxFile"
fileIdx1=\`expr \$MyCondorISection \* \$MyCondorMaxFile + 1\`
fileIdx2=\`expr \$fileIdx1 + \$MyCondorMaxFile - 1\`
fileIdx=0
rm -rf ThisInputFiles.das
while read line 
do
   ((fileIdx++))
   if [ \$fileIdx -le \$fileIdx2 ] && [ \$fileIdx -ge \$fileIdx1 ]; then
      echo \$line >> ThisInputFiles.das
      echo "@CondorThisInputFile \$line"
   fi
done < \$FirstDir/inputfiles.das
cat ThisInputFiles.das | awk '{print "'\${MyCondorInputURL}'"\$1}' > ThisInputFiles.list
echo "###################"
cat ThisInputFiles.list
echo "###################"
NumberOfInputFiles=\`cat ThisInputFiles.list | wc -l\`
ThisSumNEventsInDAS=\`cat ThisInputFiles.das | awk '{sum+=\$2} END{print sum}'\`
echo "   @ NumberOfInputFiles= \${NumberOfInputFiles} files \${ThisSumNEventsInDAS} events"
echo ""

echo "@ OutNameSet"
for name in \`cat CMSSW.outname\`
do
   echo "OutName= \${name}"
   thisOutDirName=\`dirname \${name}\`
   if [ "\${thisOutDirName}" != "." ]; then
      mkdir -p \${thisOutDirName}
      mkdir -p \${condorOutDir}/\${thisOutDirName}
   fi
done 
echo ""

echo "@ Processing nEvent Check"
IsLHE=0
MyJSON=${MyJSON}
if [ \$MyJSON ]; then
	bMyJSON=\`basename \$MyJSON\`
	echo "# GoodLumiCheck \$bMyJSON"
	for file in \`cat ThisInputFiles.list | awk '{print \$1}'\`
	do
		JsonEvents \${bMyJSON} \${file} 2>&1 
	done | tee goodLumi_\${CONDOR_ID}.log
	nJSONfile=\`grep "#edmFileJSONresult" goodLumi_\${CONDOR_ID}.log | wc -l\`
	nJSONevtGoo=\`grep "#edmFileJSONresult" goodLumi_\${CONDOR_ID}.log | awk '{sum+=\$4} END{print sum}'\`
	nJSONevtBad=\`grep "#edmFileJSONresult" goodLumi_\${CONDOR_ID}.log | awk '{sum+=\$6} END{print sum}'\`
	nJSONevtTot=\`grep "#edmFileJSONresult" goodLumi_\${CONDOR_ID}.log | awk '{sum+=\$11} END{print sum}'\`
	echo "@finalEMDJSON \$nJSONfile files \$nJSONevtTot events \$nJSONevtGoo Good \$nJSONevtBad Bad "
	echo ""
	if [ "\`cat ThisInputFiles.das | wc -l\`" != "\$nJSONfile" ]; then
		echo "#Error nFile ThisInputFiles.das and goodLumi_\${CONDOR_ID}.log"
		echo "#ThisInputFiles.das"
		cat ThisInputFiles.das
		echo "#GoodLumi file"
		cat goodLumi_\${CONDOR_ID}.log
		echo ""
		exit
	fi
else
   rm -rf edmEvents_\${CONDOR_ID}.log
   while read line
   do
      f1=\`echo \$line | awk '{print \$1}'\`
      f2=\`echo \$line | awk '{print \$2}'\`
      if [ \$f2 ]; then 
   		echo "\$f1 \$f2" >> edmEvents_\${CONDOR_ID}.log 
   	fi
      if [ ! \$f2 ]; then 
          for i in 1 2 3 4 5 6 6 6 6 6
          do
              rm -rf tempEDM.txt
              if [ "\${f1:(-5)}" == ".root" ]; then
						if [ "\${f1:0:6}" == "/store" ]; then
	              		#echo "    @edmEventSize -v root://cms-xrd-global.cern.ch/\$f1 try \$i"
   	               #edmEventSize -v root://cms-xrd-global.cern.ch/\$f1 | grep "^File" > tempEDM.txt
	              		echo "    @edmEventSize -v \${MyCondorInputURL}/\$f1 try \$i"
   	               edmEventSize -v \${MyCondorInputURL}//\$f1 | grep "^File" > tempEDM.txt
						else
	              		echo "    @edmEventSize -v \$f1 try \$i"
   	               edmEventSize -v \$f1 | grep "^File" > tempEDM.txt
						fi
						cat tempEDM.txt
              elif [ "\${f1:(-4)}" == ".lhe" ]; then
						IsLHE=1
						tempLHEName=\$f1
                  tempLHEName=\${tempLHEName/file:/}
                  nLHEevents=\`grep "<event>" \$tempLHEName | wc -l\`
                  echo "\$tempLHEName LHE File \$nLHEevents" > tempEDM.txt
						cat tempEDM.txt
					else
  						tempLHEName=\$f1
                  tempLHEName=\${tempLHEName/file:/}
  	               echo "\$tempLHEName Unkonwn File 0" > tempEDM.txt
						cat tempEDM.txt
					fi
              if [ "\`grep File tempEDM.txt | wc -l\`" == "1" ]; then 
                  f2=\`grep "File" tempEDM.txt | awk '{print \$4}'\`
                  echo "\$f1 \$f2" >> edmEvents_\${CONDOR_ID}.log
                  break
					else
						sleep "\${i}0"
              fi
          done
      fi
   done < ThisInputFiles.das
	if [ "\`cat ThisInputFiles.das | wc -l\`" != "\`cat edmEvents_\${CONDOR_ID}.log | wc -l\`" ]; then 
		echo "Error edmFileCheck"
		echo "### ThisInputFiles.das"
		cat ThisInputFiles.das
		echo "### edmEvents"
		cat edmEvents_\${CONDOR_ID}.log
		echo ""
		exit
	fi
   nJSONfile=\`cat edmEvents_\${CONDOR_ID}.log | wc -l\`
   nJSONevtGoo=\`cat edmEvents_\${CONDOR_ID}.log | awk '{sum+=\$2} END{print sum}'\`
   nJSONevtBad=0
   nJSONevtTot=\${nJSONevtGoo}
	echo "@finalEDMJSON \$nJSONfile files \$nJSONevtTot events \$nJSONevtGoo Good \$nJSONevtBad Bad "
	echo ""
fi
if [ "\${ThisSumNEventsInDAS}" == "0" ]; then ThisSumNEventsInDAS=\${nJSONevtTot} ; fi

echo "# cmsRun"
echo "@@@ cmsRun Start \`date\`"
if [ "\${IsLHE}" == "1" ]; then 
	firstRunNumber=\`expr \$MyCondorISection + 1\`
	echo "#### Set FirstRun Number is \$firstRunNumber for MC Prod"
	echo "process.source.firstRun = cms.untracked.uint32(\${firstRunNumber})" >> CMSSW.py
fi
if [ -f beforeCMD ]; then echo "### Start beforeCMD" ; source beforeCMD;  echo "### END beforeCMD"; fi
cmsRun -j cmsRunLog_\${CONDOR_ID}.xml  CMSSW.py 2>&1 | tee cmsRunLog_\${CONDOR_ID}.log
cmsRunStatus=\$?
echo "@@@ cmsRun End \${cmsRunStatus}  \`date\`"
echo ""
if [ -f afterCMD ]; then echo "### Start afterCMD" ; source afterCMD;  echo "### END afterCMD"; fi

bRootFileNames=""
for outRootFile in \`cat CMSSW.outname\`
do
	if [ "\${outRootFile:0:5}" == "file:" ]; then outRootFile=\${outRootFile/file:/}; fi ### ycyang
   thisOutDirName=\`dirname \$outRootFile\`
   bRootFile=\`basename \${outRootFile}\`
	thisFileNameName="\${bRootFile%.*}"
	thisFileNameExt="\${bRootFile##*.}"
	mv \${outRootFile}   \${FirstDir}/condorOut/\${thisOutDirName}/\${thisFileNameName}_\${CONDOR_ID}.\${thisFileNameExt}
done
echo ""


echo "# Check cmsRun Log"
if [ ! -f cmsRunLog_\${CONDOR_ID}.log ]; then echo "NotFound cmsRunLog_\${CONDOR_ID}.log" ; exit; fi
if [ ! -f cmsRunLog_\${CONDOR_ID}.xml ]; then echo "NotFound cmsRunLog_\${CONDOR_ID}.xml" ; exit; fi
mv cmsRunLog_\${CONDOR_ID}.log cmsRunLog_\${CONDOR_ID}.xml  \${cmsRunLogDir}/

echo "# Final check"
nTriEventTotal=\`grep "^TrigReport Events total = " \${cmsRunLogDir}/cmsRunLog_\${CONDOR_ID}.log | grep passed | grep failed | awk '{print \$5}'\`
nTriEventPass=\`grep "^TrigReport Events total = " \${cmsRunLogDir}/cmsRunLog_\${CONDOR_ID}.log | grep passed | grep failed | awk '{print \$8}'\`
nTriEventFail=\`grep "^TrigReport Events total = " \${cmsRunLogDir}/cmsRunLog_\${CONDOR_ID}.log | grep passed | grep failed | awk '{print \$11}'\`
echo "@finalTriReport \$nTriEventTotal total \$nTriEventPass pass \$nTriEventFail fail"

echo "@finalDAS \$NumberOfInputFiles \$ThisSumNEventsInDAS"

#nRunEvents=\`grep "^Begin processing " \${cmsRunLogDir}/cmsRunLog_\${CONDOR_ID}.log | wc -l\`
nRunEvents=\${nTriEventTotal}
echo "@finalCmsRun \$nRunEvents"

checkNFile=\`expr \$NumberOfInputFiles - \$nJSONfile\`
checkNEvt1=\`expr \$ThisSumNEventsInDAS - \$nJSONevtGoo - \$nJSONevtBad\`
checkNEvt2=\`expr \$nJSONevtGoo - \$nRunEvents\`
echo "#FinalCondorRunResult cmsRunStatus \$cmsRunStatus DiffnFiles \$checkNFile DiffEvtDAS_JSON \$checkNEvt1 DiffEvtJSONGOOD_cmsRun \$checkNEvt2"
echo "##########"
fileBit="-9"
if [ "\${checkNFile}\${checkNEvt1}\${checkNEvt2}" == "000" ]; then
	fileBit="OK"
else 
	fileBit="BAD"
fi

while read line
do
	echo "ThisRunFile\${fileBit} \$line"
done < ThisInputFiles.das 
echo ""

if [ "\${cmsRunStatus}" == "0" ] && [ "\${SRMOUTDIR:0:4}" != "NULL" ]; then
   SRMHOST=${MySRMHOST}
   SRMUSERDIR=${MySRMUSERDIR}
   echo "### OutPut To \${SRMHOST}/\${SRMUSERDIR}/\${SRMOUTDIR}"
   ln -s \$FirstDir/condorOut condorOutSRM
   for file in \`find condorOutSRM/ -type f\`
   do
      outName=\${file/condorOutSRM\//}
      echo "\$file -> \$SRMHOST/\$SRMUSERDIR/\$SRMOUTDIR/condorOut/\$outName"
      xrdcp \$file \$SRMHOST/\$SRMUSERDIR/\$SRMOUTDIR/condorOut/\$outName
      echo "\$SRMHOST/\$SRMUSERDIR/\$SRMOUTDIR/condorOut/\$outName" > \${file}.srm
		rm -rf \${file}
   done
   rm -rf condorOutSRM
   ln -s \$FirstDir/cmsRunLog cmsRunLogSRM
   for file in \`find cmsRunLogSRM/ -type f\`
   do
      outName=\${file/cmsRunLogSRM\//}
      echo "\$file -> \$SRMHOST/\$SRMUSERDIR/\$SRMOUTDIR/cmsRunLog/\$outName"
      xrdcp \$file \$SRMHOST/\$SRMUSERDIR/\$SRMOUTDIR/cmsRunLog/\$outName
      echo "\$SRMHOST/\$SRMUSERDIR/\$SRMOUTDIR/cmsRunLog/\$outName" > \${file}.srm
		rm -rf \${file}
   done
   rm -rf cmsRunLogSRM
	echo "### Done SRM"
fi
echo "### condorOut ##############################"
ls -al \$FirstDir/condorOut
echo "### cmsRunLog ##############################"
ls -al \$FirstDir/cmsRunLog
echo "#################################"
echo "Done bye bye"
EOF
	chmod +x ${myCondorWorkDir}/runCMSRUN.sh
}

function MakeLocalTest() {
cat << EOF > ${myCondorWorkDir}/.localTest.sh
#!/bin/bash 

idx=0
if [ \$1 ]; then idx=\$1; fi
thisdir=\`echo \$(cd \$(dirname "\$0") && pwd -P)\`
cd \$thisdir
mkdir -p localTest
cp \$thisdir/input.tgz \$thisdir/inputfiles.das \$thisdir/runCMSRUN.sh localTest/
cd localTest
echo "PWD \$PWD"
export CONDOR_ID=localTest_0
#./runCMSRUN.sh 1 root://cms-xrd-global.cern.ch/ 2>&1 | tee condorLog_\${CONDOR_ID}.log
./runCMSRUN.sh 1 $MyURL 2>&1 | tee condorLog_\${CONDOR_ID}.log
ls -al \$PWD
EOF
	chmod +x ${myCondorWorkDir}/.localTest.sh
}

argvPar $@
if [ "${MyCheckDir}" != "" ]; then 
	for doCheckDir in $MyCheckDir; do  Check $doCheckDir;  echo ""; done
	ReSubmitStringsP=$(printf "$ReSubmitStrings")
	echo -e "\e[1;31m"
	printf '%s\n' "${ReSubmitStringsP}"
 	echo -e "\e[0m"
	exit
fi

if [ "${MyXMLDirs}" != "" ]; then
	for MyXMLDir in $MyXMLDirs
	do
		thisDir=`readlink -e $MyXMLDir`
		name=`basename $thisDir`
		echo "Making... $name"
		fjr2json.py --output ${MyXMLDir}/${name}.json ${MyXMLDir}/cmsRunLog/*.xml
		ls -al ${MyXMLDir}/${name}.json
	done
	exit
fi


PrintSet
if [ "$CMSSW_BASE" == "" ]; then echo "NotFound CMSSW; Setup CMSSW"; exit; fi
mkdir -p ${myCondorWorkDir}
MakeDataList 
MakeCMSSW
MakeRun
MakeLocalTest
MakeJob




