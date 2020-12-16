#!/bin/bash
echo "$HOSTNAME"


TopDir=`pwd`

export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source ${VO_CMS_SW_DIR}/cmsset_default.sh
DIR=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU
cd $DIR
eval `scramv1 runtime -sh`

cd $TopDir
mcexeFile=$DIR/Test_trkIso_mc.exe
dataexeFile=$DIR/Test_trkIso_data.exe
export LD_LIBRARY_PATH=${CMSSW_BASE}/src/MakeNpKNU/ProdNpKNU/src/:${LD_LIBRARY_PATH}
echo $exeFile

if [ $1 == "0" ]; then
  $mcexeFiles 0     16.146  0.0029432   20000     /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WprimeToMuNu_M-3800_PUMoriond17.root ; fi
if [ $1 == "1" ]; then
  $mcexeFiles 0     16.146  0.27238     20000     /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WprimeToMuNu_M-1800_PUMoriond17.root ; fi
if [ $1 == "3" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunB_reMiniaod_tot.root ; fi 
if [ $1 == "4" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunC_reMiniaod_tot.root ; fi 
if [ $1 == "5" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunD_reMiniaod_tot.root ; fi 
if [ $1 == "6" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunE_reMiniaod_tot.root ; fi 
if [ $1 == "7" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunF_reMiniaod_tot.root ; fi 
if [ $1 == "8" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunG_reMiniaod_tot.root ; fi 
if [ $1 == "9" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH2_reMiniaod_tot.root; fi
if [ $1 == "10" ]; then
   $dataexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH3_reMiniaod.root    ; fi


 
   #$exeFile  2000  16.146  0.000556  997582  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-2000_PUMoriond17.root ; fi
   #$exeFile  0     16.146   0.0013     94799     /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_1400_2300_PUMoriond17.root ; fi
