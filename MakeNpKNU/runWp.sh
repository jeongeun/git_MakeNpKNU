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
exeb=$DIR/MakeWpHist_datab_17.exe
exe=$DIR/MakeWpHist_data_16.exe

export LD_LIBRARY_PATH=${DIR}/ProdNpKNU/src/:${LD_LIBRARY_PATH}

outName=$2

if [ ! -d condorOut ]; then mkdir condorOut; fi
if [ $1 == "0" ]; then
  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunD_v1_1.root  condorOut/hEle_RunD_v1_2016_1.root ; fi
if [ $1 == "1" ]; then                                                               
  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunD_v1_2.root  condorOut/hEle_RunD_v1_2016_2.root ; fi
if [ $1 == "2" ]; then                                                               
  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunD_v1_3.root  condorOut/hEle_RunD_v1_2016_3.root ; fi
if [ $1 == "3" ]; then                                                               
  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunD_v1_4.root  condorOut/hEle_RunD_v1_2016_4.root ; fi
if [ $1 == "4" ]; then                                                               
  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunD_v1_5.root  condorOut/hEle_RunD_v1_2016_5.root ; fi
if [ $1 == "5" ]; then                                                               
  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunD_v1_6.root  condorOut/hEle_RunD_v1_2016_6.root ; fi
#if [ $1 == "6" ]; then                                                               
#  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunG_v1_7.root  condorOut/hEle_RunG_v1_2016_7.root ; fi
#if [ $1 == "7" ]; then                                                                
#  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunG_v1_8.root  condorOut/hEle_RunG_v1_2016_8.root ; fi
##if [ $1 == "8" ]; then                                                                
##  $exe 1 1 1 1 ${DIR}/ntuple2016/SE2016/ntNpKNUdata2016_RunF_v2.root  condorOut/hEle_RunF_v2.root ; fi
##
