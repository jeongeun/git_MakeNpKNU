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
MEexeFile=$DIR/MakeWpHist_Data_ME.exe
MBexeFile=$DIR/MakeWpHist_Data_MB.exe

#MBhighexe=$DIR/FindHighMT_data_MB.exe
#MEhighexe=$DIR/FindHighMT_data_ME.exe
#runHexe=$DIR/runH2Test.exe
#runH2exe=$DIR/runH_Test.exe
export LD_LIBRARY_PATH=${CMSSW_BASE}/src/MakeNpKNU/ProdNpKNU/src/:${LD_LIBRARY_PATH}

outName=$2

if [ ! -d condorOut ]; then mkdir condorOut; fi
#
#if [ $1 == "0" ]; then
#$runH2exe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH2_reMiniaod_tot.root ; fi
#   $MBhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunD_1.root condorOut/hmt_ADD_RunD1_MB.root ; fi
#if [ $1 == "1" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunD_1.root condorOut/hmt_ADD_RunD1_ME.root ; fi
#if [ $1 == "2" ]; then
#   $MBhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunD_2.root condorOut/hmt_ADD_RunD2_MB.root ; fi
#if [ $1 == "3" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunD_2.root condorOut/hmt_ADD_RunD2_ME.root ; fi
#if [ $1 == "4" ]; then
#   $MBhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunF.root condorOut/hmt_ADD_RunF_MB.root ; fi
#if [ $1 == "5" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunF.root condorOut/hmt_ADD_RunF_ME.root ; fi
#if [ $1 == "6" ]; then
#   $MBhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunG.root condorOut/hmt_ADD_RunG_MB.root ; fi
#if [ $1 == "7" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_highMt_RunG.root condorOut/hmt_ADD_RunG_ME.root ; fi
if [ $1 == "0" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunB_reMiniaod_tot.root   condorOut/hData_RunB_reMiniaod_MB.root  ; fi
if [ $1 == "1" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunC_reMiniaod_tot.root   condorOut/hData_RunC_reMiniaod_MB.root  ; fi
if [ $1 == "2" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunD_reMiniaod_tot.root   condorOut/hData_RunD_reMiniaod_MB.root  ; fi
if [ $1 == "3" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunE_reMiniaod_tot.root   condorOut/hData_RunE_reMiniaod_MB.root  ; fi
if [ $1 == "4" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunF_reMiniaod_tot.root   condorOut/hData_RunF_reMiniaod_MB.root  ; fi
if [ $1 == "5" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunG_reMiniaod_tot.root   condorOut/hData_RunG_reMiniaod_MB.root  ; fi
if [ $1 == "6" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH2_reMiniaod_tot.root  condorOut/hData_RunH2_reMiniaod_MB.root ; fi
if [ $1 == "7" ]; then
$MBexeFile  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH3_reMiniaod.root      condorOut/hData_RunH3_reMiniaod_MB.root ; fi
#if [ $1 == "8" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunB_reMiniaod_tot.root   condorOut/hmt_RunB_reMiniaod_ME.root  ; fi
#if [ $1 == "9" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunC_reMiniaod_tot.root   condorOut/hmt_RunC_reMiniaod_ME.root  ; fi
#if [ $1 == "10" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunD_reMiniaod_tot.root   condorOut/hmt_RunD_reMiniaod_ME.root  ; fi
#if [ $1 == "11" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunE_reMiniaod_tot.root   condorOut/hmt_RunE_reMiniaod_ME.root  ; fi
#if [ $1 == "12" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunF_reMiniaod_tot.root   condorOut/hmt_RunF_reMiniaod_ME.root  ; fi
#if [ $1 == "13" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunG_reMiniaod_tot.root   condorOut/hmt_RunG_reMiniaod_ME.root  ; fi
#if [ $1 == "14" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH2_reMiniaod_tot.root  condorOut/hmt_RunH2_reMiniaod_ME.root ; fi
#if [ $1 == "15" ]; then
#   $MEhighexe  1 1 1 1 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/ntuples/NtreMiniaod/ntNpKNUdata_RunH3_reMiniaod.root      condorOut/hmt_RunH3_reMiniaod_ME.root ; fi
