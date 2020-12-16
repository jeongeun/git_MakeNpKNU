#!/bin/bash

echo "$HOSTNAME"

TopDir=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU
#`pwd`
export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source ${VO_CMS_SW_DIR}/cmsset_default.sh
DIR=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana
cd $DIR
eval `scramv1 runtime -sh`
cd $TopDir
BFMBexeFile=$TopDir/MakeWpHist_MC_BF_MB_tot.exe
GHMBexeFile=$TopDir/MakeWpHist_MC_GH_MB_tot.exe
BFMEexeFile=$TopDir/MakeWpHist_MC_BF_ME_tot.exe
GHMEexeFile=$TopDir/MakeWpHist_MC_GH_ME_tot.exe

export LD_LIBRARY_PATH=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ProdNpKNU/src/:${LD_LIBRARY_PATH}

outName=$2

if [ ! -d condorOut ]; then mkdir condorOut; fi

###if [ $1 == "0" ]; then
###$GHMBexeFile  100  16.1  163.15      1985397   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-100_PUMoriond17.root    condorOut/hMC_W_M100_GH_MB.root  ; fi 
###if [ $1 == "1" ]; then
###$GHMBexeFile  200  16.1  6.236       996128    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-200_PUMoriond17.root    condorOut/hMC_W_M200_GH_MB.root  ; fi 
###if [ $1 == "2" ]; then
###$GHMBexeFile  500  16.1  0.2138      997511    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-500_PUMoriond17.root    condorOut/hMC_W_M500_GH_MB.root  ; fi 
###if [ $1 == "3" ]; then
###$GHMBexeFile  1000 16.1  0.01281     1000000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-1000_PUMoriond17.root   condorOut/hMC_W_M1000_GH_MB.root ; fi 
###if [ $1 == "4" ]; then
###$GHMBexeFile  2000 16.1  0.000556    997582    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-2000_PUMoriond17.root   condorOut/hMC_W_M2000_GH_MB.root ; fi 
###if [ $1 == "5" ]; then
###$GHMBexeFile  3000 16.1  0.00002904  997236    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-3000_PUMoriond17.root   condorOut/hMC_W_M3000_GH_MB.root ; fi 
###if [ $1 == "6" ]; then
###$GHMBexeFile  99   16.1  61526.7     24120319  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_amcatnloFXFX_PUMoriond17.root    condorOut/hMC_WJets_amc_GH_MB.root  ; fi 
###if [ $1 == "7" ]; then
###$GHMBexeFile  99   16.1  61526.7     57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root   condorOut/hMC_WJets_mad_GH_MB.root  ; fi 
#if [ $1 == "0" ]; then
#$GHMBexeFile  0    16.1  831.76      77229341  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17.root                         condorOut/hMC_TT_GH_MB.root         ; fi 
#if [ $1 == "1" ]; then
#$GHMBexeFile  0    16.1  3.36        1000000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_s-channel_PUMoriond17.root               condorOut/hMC_ST_sch_GH_MB.root     ; fi 
#if [ $1 == "2" ]; then
#$GHMBexeFile  0    16.1  44.33       5993676   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_t-channel_top_PUMoriond17.root           condorOut/hMC_ST_tch_t_GH_MB.root   ; fi 
#if [ $1 == "3" ]; then
#$GHMBexeFile  0    16.1  26.38       3928063   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_t-channel_antitop_PUMoriond17.root       condorOut/hMC_ST_tch_at_GH_MB.root  ; fi 
#if [ $1 == "4" ]; then
#$GHMBexeFile  0    16.1  35.85        992024   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_tW_top_PUMoriond17.root                 condorOut/hMC_ST_tWch_t_GH_MB.root  ; fi
#if [ $1 == "5" ]; then
#$GHMBexeFile  0    16.1  35.85        998276   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_tW_antitop_PUMoriond17.root             condorOut/hMC_ST_tWch_at_GH_MB.root ; fi
#if [ $1 == "6" ]; then
#$GHMBexeFile  0    16.1  1975         2977600  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_50_120_PUMoriond17.root          condorOut/hMC_Z_M_50_GH_MB.root     ; fi
#if [ $1 == "15" ]; then
#$GHMBexeFile  0    16.1  19.3         100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_120_200_PUMoriond17.root         condorOut/hMC_Z_M_120_GH_MB.root    ; fi
#if [ $1 == "16" ]; then
#$GHMBexeFile  0    16.1  2.73         100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_200_400_PUMoriond17.root         condorOut/hMC_Z_M_200_GH_MB.root    ; fi
#if [ $1 == "17" ]; then
#$GHMBexeFile  0    16.1  0.241        98400    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_400_800_PUMoriond17.root         condorOut/hMC_Z_M_400_GH_MB.root    ; fi
#if [ $1 == "18" ]; then
#$GHMBexeFile  0    16.1  0.0168       100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_800_1400_PUMoriond17.root        condorOut/hMC_Z_M_800_GH_MB.root    ; fi
#if [ $1 == "19" ]; then
#$GHMBexeFile  0    16.1   0.0013      94799    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_1400_2300_PUMoriond17.root       condorOut/hMC_Z_M_1400_GH_MB.root   ; fi        
#if [ $1 == "20" ]; then
#$GHMBexeFile  0    16.1  0.0008948    100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_2300_3500_PUMoriond17.root       condorOut/hMC_Z_M_2300_GH_MB.root   ; fi
#if [ $1 == "21" ]; then
#$GHMBexeFile  0    16.1  0.00004135   100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_3500_4500_PUMoriond17.root       condorOut/hMC_Z_M_3500_GH_MB.root   ; fi
#if [ $1 == "22" ]; then
#$GHMBexeFile  0    16.1  0.000000456  100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_4500_6000_PUMoriond17.root       condorOut/hMC_Z_M_4500_GH_MB.root   ; fi
#if [ $1 == "23" ]; then
#$GHMBexeFile  0    16.1  0.0000000206  100000  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_6000_Inf_PUMoriond17.root       condorOut/hMC_Z_M_6000_GH_MB.root   ; fi  
#if [ $1 == "24" ]; then
#$GHMBexeFile  0    16.1  47.13         1000000 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WZ_PUMoriond17.root                       condorOut/hMC_WZ_GH_MB.root         ; fi 
#if [ $1 == "25" ]; then
#$GHMBexeFile  0    16.1  23.23         990064  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZZ_PUMoriond17.root                       condorOut/hMC_ZZ_GH_MB.root         ; fi 
#if [ $1 == "26" ]; then
#$GHMBexeFile  0    16.1  118.7         994012  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WW_PUMoriond17.root                       condorOut/hMC_WW_GH_MB.root         ; fi 
#if [ $1 == "27" ]; then
#$GHMBexeFile  0    16.1  0.27238       20000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WprimeToMuNu_M-1800_PUMoriond17.root      condorOut/hMC_Wprime_1800_GH_MB.root; fi  
#if [ $1 == "28" ]; then
#$GHMBexeFile  0    16.1  0.0029432     20000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WprimeToMuNu_M-3800_PUMoriond17.root      condorOut/hMC_Wprime_3800_GH_MB.root; fi  



###if [ $1 == "29" ]; then
###$GHMEexeFile  100  16.1  163.15      1985397   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-100_PUMoriond17.root    condorOut/hMC_W_M100_GH_ME.root  ; fi 
###if [ $1 == "30" ]; then
###$GHMEexeFile  200  16.1  6.236       996128    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-200_PUMoriond17.root    condorOut/hMC_W_M200_GH_ME.root  ; fi 
###if [ $1 == "31" ]; then
###$GHMEexeFile  500  16.1  0.2138      997511    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-500_PUMoriond17.root    condorOut/hMC_W_M500_GH_ME.root  ; fi 
###if [ $1 == "32" ]; then
###$GHMEexeFile  1000 16.1  0.01281     1000000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-1000_PUMoriond17.root   condorOut/hMC_W_M1000_GH_ME.root ; fi 
###if [ $1 == "33" ]; then
###$GHMEexeFile  2000 16.1  0.000556    997582    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-2000_PUMoriond17.root   condorOut/hMC_W_M2000_GH_ME.root ; fi 
###if [ $1 == "34" ]; then
###$GHMEexeFile  3000 16.1  0.00002904  997236    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WToMuNu_M-3000_PUMoriond17.root   condorOut/hMC_W_M3000_GH_ME.root ; fi 
###if [ $1 == "35" ]; then
###$GHMEexeFile  99   16.1  61526.7     24120319  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_amcatnloFXFX_PUMoriond17.root    condorOut/hMC_WJets_amc_GH_ME.root  ; fi 
###if [ $1 == "36" ]; then
###$GHMEexeFile  99   16.1  61526.7     57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root   condorOut/hMC_WJets_mad_GH_ME.root  ; fi 
###if [ $1 == "37" ]; then
###$GHMEexeFile  0    16.1  831.76      77229341  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_TT_PUMoriond17.root                         condorOut/hMC_TT_GH_ME.root         ; fi 
###if [ $1 == "38" ]; then
###$GHMEexeFile  0    16.1  3.36        1000000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_s-channel_PUMoriond17.root               condorOut/hMC_ST_sch_GH_ME.root     ; fi 
###if [ $1 == "39" ]; then
###$GHMEexeFile  0    16.1  44.33       5993676   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_t-channel_top_PUMoriond17.root           condorOut/hMC_ST_tch_t_GH_ME.root   ; fi 
###if [ $1 == "40" ]; then
###$GHMEexeFile  0    16.1  26.38       3928063   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_t-channel_antitop_PUMoriond17.root       condorOut/hMC_ST_tch_at_GH_ME.root  ; fi 
###if [ $1 == "41" ]; then
###$GHMEexeFile  0    16.1  35.85        992024   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_tW_top_PUMoriond17.root                 condorOut/hMC_ST_tWch_t_GH_ME.root  ; fi
###if [ $1 == "42" ]; then
###$GHMEexeFile  0    16.1  35.85        998276   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ST_tW_antitop_PUMoriond17.root             condorOut/hMC_ST_tWch_at_GH_ME.root ; fi
###if [ $1 == "43" ]; then
###$GHMEexeFile  0    16.1  1975         2977600  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_50_120_PUMoriond17.root          condorOut/hMC_Z_M_50_GH_ME.root     ; fi
###if [ $1 == "44" ]; then
###$GHMEexeFile  0    16.1  19.3         100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_120_200_PUMoriond17.root         condorOut/hMC_Z_M_120_GH_ME.root    ; fi
###if [ $1 == "45" ]; then
###$GHMEexeFile  0    16.1  2.73         100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_200_400_PUMoriond17.root         condorOut/hMC_Z_M_200_GH_ME.root    ; fi
###if [ $1 == "46" ]; then
###$GHMEexeFile  0    16.1  0.241        98400    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_400_800_PUMoriond17.root         condorOut/hMC_Z_M_400_GH_ME.root    ; fi
###if [ $1 == "47" ]; then
###$GHMEexeFile  0    16.1  0.0168       100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_800_1400_PUMoriond17.root        condorOut/hMC_Z_M_800_GH_ME.root    ; fi
###if [ $1 == "48" ]; then
###$GHMEexeFile  0    16.1   0.0013      94799    /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_1400_2300_PUMoriond17.root       condorOut/hMC_Z_M_1400_GH_ME.root   ; fi        
###if [ $1 == "49" ]; then
###$GHMEexeFile  0    16.1  0.0008948    100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_2300_3500_PUMoriond17.root       condorOut/hMC_Z_M_2300_GH_ME.root   ; fi
###if [ $1 == "50" ]; then
###$GHMEexeFile  0    16.1  0.00004135   100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_3500_4500_PUMoriond17.root       condorOut/hMC_Z_M_3500_GH_ME.root   ; fi
###if [ $1 == "51" ]; then
###$GHMEexeFile  0    16.1  0.000000456  100000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_4500_6000_PUMoriond17.root       condorOut/hMC_Z_M_4500_GH_ME.root   ; fi
###if [ $1 == "52" ]; then
###$GHMEexeFile  0    16.1  0.0000000206  100000  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZToMuMu_M_6000_Inf_PUMoriond17.root       condorOut/hMC_Z_M_6000_GH_ME.root   ; fi  
###if [ $1 == "53" ]; then
###$GHMEexeFile  0    16.1  47.13         1000000 /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WZ_PUMoriond17.root                       condorOut/hMC_WZ_GH_ME.root         ; fi 
###if [ $1 == "54" ]; then
###$GHMEexeFile  0    16.1  23.23         990064  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_ZZ_PUMoriond17.root                       condorOut/hMC_ZZ_GH_ME.root         ; fi 
###if [ $1 == "55" ]; then
###$GHMEexeFile  0    16.1  118.7         994012  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WW_PUMoriond17.root                       condorOut/hMC_WW_GH_ME.root         ; fi 
if [ $1 == "0" ]; then
$GHMEexeFile  0    16.1  0.27238       20000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WprimeToMuNu_M-1800_PUMoriond17.root      condorOut/hMC_Wprime_1800_GH_ME.root; fi  
if [ $1 == "1" ]; then
$GHMEexeFile  0    16.1  0.0029432     20000   /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WprimeToMuNu_M-3800_PUMoriond17.root      condorOut/hMC_Wprime_3800_GH_ME.root; fi  

