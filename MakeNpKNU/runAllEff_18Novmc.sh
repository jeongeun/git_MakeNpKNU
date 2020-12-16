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
$EEexe  0 35.9  61526.7  57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_April.root   condorOut/hEff_WJets_mad_EE_pho.root  ; fi
if [ $1 == "1" ]; then
$EEexe  0 35.9  1627.45   29503700   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-100To200_PUMoriond17_April.root   condorOut/hEff_WJets_HT100To200_EE_pho.root  ; fi
if [ $1 == "2" ]; then
$EEexe  0 35.9  86.89745    8068473  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17_April_0.root   condorOut/hEff_TT_EE_0.root  ; fi
if [ $1 == "3" ]; then
$EEexe  0 35.9  744.86  69160868  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17_April_1.root   condorOut/hEff_TT_EE_1.root  ; fi
if [ $1 == "4" ]; then
$EEexe  0 35.9  59.181    5796237    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-400To600_PUMoriond17_April.root  condorOut/hEff_WJets_HT400To600_EE_pho.root  ; fi
if [ $1 == "5" ]; then
$EEexe  0 35.9  435.2369  19914590   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-200To400_PUMoriond17_April.root  condorOut/hEff_WJets_HT200To400_EE_pho.root  ; fi
if [ $1 == "6" ]; then
$EEexe  0 35.9    6.6562   6200954   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-800To1200_PUMoriond17_April.root  condorOut/hEff_WJets_HT800To1200_EE_pho.root  ; fi
if [ $1 == "7" ]; then
$EEexe  0 35.9    1.60809  6627909   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-1200To2500_PUMoriond17_April.root  condorOut/hEff_WJets_HT1200To2500_EE_pho.root  ; fi
if [ $1 == "8" ]; then
$EEexe  0 35.9   0.03891359  2384260  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-2500ToInf_PUMoriond17_April.root  condorOut/hEff_WJets_HT2500ToInf_EE_pho.root  ; fi

###if [ $1 == "0" ]; then
###$EBexe  0 35.9  61526.7  57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_April.root   condorOut/hEff_WJets_mad_EB_pho.root  ; fi
###if [ $1 == "1" ]; then
###$EEexe  0 35.9  61526.7  57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_April.root   condorOut/hEff_WJets_mad_EE_pho.root  ; fi
###if [ $1 == "2" ]; then
###$EBexe  0 35.9  1627.45   29503700  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-100To200_PUMoriond17_April.root   condorOut/hEff_WJets_HT100To200_EB_pho.root  ; fi
###if [ $1 == "3" ]; then
###$EEexe  0 35.9  1627.45   29503700   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-100To200_PUMoriond17_April.root   condorOut/hEff_WJets_HT100To200_EE_pho.root  ; fi
###if [ $1 == "4" ]; then
###$EBexe  0 35.9  86.89745    8068473  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17_April_0.root   condorOut/hEff_TT_EB_0.root  ; fi
###if [ $1 == "5" ]; then
###$EEexe  0 35.9  86.89745    8068473  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17_April_0.root   condorOut/hEff_TT_EE_0.root  ; fi
###if [ $1 == "6" ]; then
###$EBexe  0 35.9  744.86   69160868  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17_April_1.root   condorOut/hEff_TT_EB_1.root  ; fi
###if [ $1 == "7" ]; then
###$EEexe  0 35.9  744.86  69160868  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17_April_1.root   condorOut/hEff_TT_EE_1.root  ; fi
###if [ $1 == "8" ]; then
###$EBexe  0 35.9  59.181    5796237    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-400To600_PUMoriond17_April.root  condorOut/hEff_WJets_HT400To600_EB_pho.root  ; fi
###if [ $1 == "9" ]; then
###$EEexe  0 35.9  59.181    5796237    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-400To600_PUMoriond17_April.root  condorOut/hEff_WJets_HT400To600_EE_pho.root  ; fi
###if [ $1 == "10" ]; then
###$EBexe  0 35.9  435.2369  19914590   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-200To400_PUMoriond17_April.root  condorOut/hEff_WJets_HT200To400_EB_pho.root  ; fi
###if [ $1 == "11" ]; then
###$EEexe  0 35.9  435.2369  19914590   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-200To400_PUMoriond17_April.root  condorOut/hEff_WJets_HT200To400_EE_pho.root  ; fi
###if [ $1 == "12" ]; then
###$EBexe  0 35.9    6.6562   6200954   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-800To1200_PUMoriond17_April.root  condorOut/hEff_WJets_HT800To1200_EB_pho.root  ; fi
###if [ $1 == "13" ]; then
###$EEexe  0 35.9    6.6562   6200954   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-800To1200_PUMoriond17_April.root  condorOut/hEff_WJets_HT800To1200_EE_pho.root  ; fi
###if [ $1 == "14" ]; then
###$EBexe  0 35.9    1.60809  6627909   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-1200To2500_PUMoriond17_April.root  condorOut/hEff_WJets_HT1200To2500_EB_pho.root  ; fi
###if [ $1 == "15" ]; then
###$EEexe  0 35.9    1.60809  6627909   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-1200To2500_PUMoriond17_April.root  condorOut/hEff_WJets_HT1200To2500_EE_pho.root  ; fi
###if [ $1 == "16" ]; then
###$EBexe  0 35.9    14.5805  14908339   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-600To800_PUMoriond17_April.root  condorOut/hEff_WJets_HT600To800_EB_pho.root  ; fi
###if [ $1 == "17" ]; then
###$EEexe  0 35.9    14.5805  14908339   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-600To800_PUMoriond17_April.root  condorOut/hEff_WJets_HT600To800_EE_pho.root  ; fi
###if [ $1 == "18" ]; then
###$EBexe  0 35.9   0.03891359  2384260  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-2500ToInf_PUMoriond17_April.root  condorOut/hEff_WJets_HT2500ToInf_EB_pho.root  ; fi
###if [ $1 == "19" ]; then
###$EEexe  0 35.9   0.03891359  2384260  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-2500ToInf_PUMoriond17_April.root  condorOut/hEff_WJets_HT2500ToInf_EE_pho.root  ; fi
###
