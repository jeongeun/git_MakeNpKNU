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
EBexe=$DIR/MakeEff_EB_18Nov.exe
EEexe=$DIR/MakeEff_EE_18Nov.exe

export LD_LIBRARY_PATH=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ProdNpKNU/src/:${LD_LIBRARY_PATH}

outName=$2

if [ ! -d condorOut ]; then mkdir condorOut; fi
if [ $1 == "0" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016B-03Feb2017.root condorOut/hEff_Run2016B_reMiniaod_EE_pho.root ; fi
if [ $1 == "1" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016C-03Feb2017.root condorOut/hEff_Run2016C_reMiniaod_EE_pho.root ; fi
if [ $1 == "2" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016D-03Feb2017.root condorOut/hEff_Run2016D_reMiniaod_EE_pho.root ; fi
if [ $1 == "3" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016E-03Feb2017.root condorOut/hEff_Run2016E_reMiniaod_EE_pho.root ; fi
if [ $1 == "4" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016F-03Feb2017.root condorOut/hEff_Run2016F_reMiniaod_EE_pho.root ; fi
if [ $1 == "5" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016G-03Feb2017.root condorOut/hEff_Run2016G_reMiniaod_EE_pho.root ; fi
if [ $1 == "6" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016H-03Feb2017_ver2.root condorOut/hEff_Run2016Hv2_reMiniaod_EE_pho.root ; fi
if [ $1 == "7" ]; then
  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016H-03Feb2017_ver3.root condorOut/hEff_Run2016Hv3_reMiniaod_EE_pho.root ; fi

###if [ $1 == "0" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016B-03Feb2017.root condorOut/hEff_Run2016B_reMiniaod_EB_pho.root ; fi
###if [ $1 == "1" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016B-03Feb2017.root condorOut/hEff_Run2016B_reMiniaod_EE_pho.root ; fi
###if [ $1 == "2" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016C-03Feb2017.root condorOut/hEff_Run2016C_reMiniaod_EB_pho.root ; fi
###if [ $1 == "3" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016C-03Feb2017.root condorOut/hEff_Run2016C_reMiniaod_EE_pho.root ; fi
###if [ $1 == "4" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016D-03Feb2017.root condorOut/hEff_Run2016D_reMiniaod_EB_pho.root ; fi
###if [ $1 == "5" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016D-03Feb2017.root condorOut/hEff_Run2016D_reMiniaod_EE_pho.root ; fi
###if [ $1 == "6" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016E-03Feb2017.root condorOut/hEff_Run2016E_reMiniaod_EB_pho.root ; fi
###if [ $1 == "7" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016E-03Feb2017.root condorOut/hEff_Run2016E_reMiniaod_EE_pho.root ; fi
###if [ $1 == "8" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016F-03Feb2017.root condorOut/hEff_Run2016F_reMiniaod_EB_pho.root ; fi
###if [ $1 == "9" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016F-03Feb2017.root condorOut/hEff_Run2016F_reMiniaod_EE_pho.root ; fi
###if [ $1 == "10" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016G-03Feb2017.root condorOut/hEff_Run2016G_reMiniaod_EB_pho.root ; fi
###if [ $1 == "11" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016G-03Feb2017.root condorOut/hEff_Run2016G_reMiniaod_EE_pho.root ; fi
###if [ $1 == "12" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016H-03Feb2017_ver2.root condorOut/hEff_Run2016Hv2_reMiniaod_EB_pho.root ; fi
###if [ $1 == "13" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016H-03Feb2017_ver2.root condorOut/hEff_Run2016Hv2_reMiniaod_EE_pho.root ; fi
###if [ $1 == "0" ]; then
###  $EBexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016H-03Feb2017_ver3.root condorOut/hEff_Run2016Hv3_reMiniaod_EB_pho.root ; fi
###if [ $1 == "1" ]; then
###  $EEexe 1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/ntupleData/ntNpKNUdata_SingleMuon_Run2016H-03Feb2017_ver3.root condorOut/hEff_Run2016Hv3_reMiniaod_EE_pho.root ; fi

