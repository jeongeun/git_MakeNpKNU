#!/bin/bash

echo "$HOSTNAME"

TopDir=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU
export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source ${VO_CMS_SW_DIR}/cmsset_default.sh
DIR=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana
cd $DIR
eval `scramv1 runtime -sh`
cd $TopDir
#BFexeFile=$DIR/MakeWpHist_mc_17_wjet_BF.exe
#GHexeFile=$DIR/MakeWpHist_mc_17_wjet_GH.exe
BFMBexeFile=$TopDir/MakeWpHist_MC_BF_MB.exe
GHMBexeFile=$TopDir/MakeWpHist_MC_GH_MB.exe
BFMEexeFile=$TopDir/MakeWpHist_MC_BF_ME.exe
GHMEexeFile=$TopDir/MakeWpHist_MC_GH_ME.exe

export LD_LIBRARY_PATH=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ProdNpKNU/src/:${LD_LIBRARY_PATH}

outName=$2
#if [ $1 == "6" ]; then
#$exeFile  9    16.1  61526.7     24120319  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_amcatnloFXFX_PUMoriond17.root    condorOut/hMC_WJets_amc_GH_GenTest.root  ; fi 

if [ ! -d condorOut ]; then mkdir condorOut; fi
if [ $1 == "0" ]; then
$GHMEexeFile  100  16.1  163.15      1985397   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-100_PUMoriond17.root    condorOut/hMC_W_M100_GH_ME.root  ; fi 
if [ $1 == "1" ]; then
$GHMEexeFile  200  16.1  6.236       996128    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-200_PUMoriond17.root    condorOut/hMC_W_M200_GH_ME.root  ; fi 
if [ $1 == "2" ]; then
$GHMEexeFile  500  16.1  0.2138      997511    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-500_PUMoriond17.root    condorOut/hMC_W_M500_GH_ME.root  ; fi 
if [ $1 == "3" ]; then
$GHMEexeFile  1000 16.1  0.01281     1000000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-1000_PUMoriond17.root   condorOut/hMC_W_M1000_GH_ME.root ; fi 
if [ $1 == "4" ]; then
$GHMEexeFile  2000 16.1  0.000556    997582    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-2000_PUMoriond17.root   condorOut/hMC_W_M2000_GH_ME.root ; fi 
if [ $1 == "5" ]; then
$GHMEexeFile  3000 16.1  0.00002904  997236    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-3000_PUMoriond17.root   condorOut/hMC_W_M3000_GH_ME.root ; fi 
if [ $1 == "6" ]; then
$GHMEexeFile  9    16.1  61526.7     57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_April.root   condorOut/hMC_WJets_mad_GH_ME.root  ; fi 
if [ $1 == "7" ]; then
$GHMEexeFile  92   16.1  435.2369  19914590   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-200To400_PUMoriond17_April.root  condorOut/hMC_WJets_HT200To400_GH_ME.root  ; fi 
if [ $1 == "8" ]; then
$GHMEexeFile  94   16.1  59.181    5796237    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-400To600_PUMoriond17_April.root  condorOut/hMC_WJets_HT400To600_GH_ME.root  ; fi 
if [ $1 == "9" ]; then
$GHMEexeFile  96   16.1  14.5805  14908339   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-600To800_PUMoriond17_April.root  condorOut/hMC_WJets_HT600To800_GH_ME.root  ; fi 
if [ $1 == "10" ]; then
$GHMEexeFile  98   16.1  6.6562   6200954   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-800To1200_PUMoriond17_April.root  condorOut/hMC_WJets_HT800To1200_GH_ME.root  ; fi 
if [ $1 == "11" ]; then
$GHMEexeFile  912  16.1  1.60809  6627909   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-1200To2500_PUMoriond17_April.root  condorOut/hMC_WJets_HT1200-2500_GH_ME.root  ; fi 
if [ $1 == "12" ]; then
$GHMEexeFile  925   16.1 0.03891359  2384260  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_HT-2500ToInf_PUMoriond17_April.root  condorOut/hMC_WJets_HT2500-Inf_GH_ME.root  ; fi 
