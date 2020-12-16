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
BFexeFile=$DIR/MakeWpHist_mc_17_BtoF.exe
GHexeFile=$DIR/MakeWpHist_mc_17_GtoH.exe
export LD_LIBRARY_PATH=/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ProdNpKNU/src/:${LD_LIBRARY_PATH}
outName=$2
if [ ! -d condorOut ]; then mkdir condorOut; fi

if [ $1 == "0" ]; then
$BFexeFile  99   19.7  61526.7     57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root   condorOut/hMC_WJets_mad_BF.root  ; fi 
if [ $1 == "1" ]; then
$GHexeFile  99   16.1  61526.7     57026058  /hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_25/src/MakeNpKNU/ana/ntuples/ntNpKNUmc_WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root   condorOut/hMC_WJets_mad_GH.root  ; fi 


#
####BACKUP
#####------------------------------------------
##./MakeHist_BF_mc.exe 100  20.236  163.15*1.0      1985397   ./ntuples/ntNpKNUmc_WToMuNu_M-100_PUMoriond17.root            hMCera1_W_M-100.root    >& lera1_W_100.log &
##./MakeHist_BF_mc.exe 200  20.236  6.236*1.1       996128    ./ntuples/ntNpKNUmc_WToMuNu_M-200_PUMoriond17.root            hMCera1_W_M-200.root    >& lera1_W_200.log &
##./MakeHist_BF_mc.exe 500  20.236  0.2138*1.3      997511    ./ntuples/ntNpKNUmc_WToMuNu_M-500_PUMoriond17.root            hMCera1_W_M-500.root    >& lera1_W_500.log &
##./MakeHist_BF_mc.exe 1000 20.236  0.01281*1.1     1000000   ./ntuples/ntNpKNUmc_WToMuNu_M-1000_PUMoriond17.root           hMCera1_W_M-1000.root   >& lera1_W_1000.log &
##./MakeHist_BF_mc.exe 2000 20.236  0.000556*1.2    997582    ./ntuples/ntNpKNUmc_WToMuNu_M-2000_PUMoriond17.root           hMCera1_W_M-2000.root   >& lera1_W_2000.log &
##./MakeHist_BF_mc.exe 3000 20.236  0.00002904*0.85 997236    ./ntuples/ntNpKNUmc_WToMuNu_M-3000_PUMoriond17.root           hMCera1_W_M-3000.root   >& lera1_W_3000.log &
#####-------------------------------------------
##./MakeHist_BF_mc.exe 99   20.236  61526.7         24120319  ./ntuples/ntNpKNUmc_WJetsToLNu_amcatnloFXFX_PUMoriond17.root  hMCera1_WJets_amc.root  >& lera1_WJet.log &
####------------------------------------------
##./MakeHist_BF_mc.exe 0    20.236  831.76          77229341  ./ntuples/ntNpKNUmc_TT_PUMoriond17.root        hMCera1_TT.root         >& lera1_TT.log &
##./MakeHist_BF_mc.exe 0    20.236  3.36            1000000   ./ntuples/ntNpKNUmc_ST_s-channel_PUMoriond17.root             hMCera1_ST_sch.root     >& lera1_ST_sch.log &
##./MakeHist_BF_mc.exe 0    20.236  44.33           5993676   ./ntuples/ntNpKNUmc_ST_t-channel_top_PUMoriond17.root         hMCera1_ST_tch_t.root   >& lera1_ST_tch_t.log &
##./MakeHist_BF_mc.exe 0    20.236  26.38           3928063   ./ntuples/ntNpKNUmc_ST_t-channel_antitop_PUMoriond17.root     hMCera1_ST_tch_at.root  >& lera1_ST_tch_at.log &
##./MakeHist_BF_mc.exe 0    20.236  35.85            992024    ./ntuples/ntNpKNUmc_ST_tW_top_PUMoriond17.root                hMCera1_ST_tWch_t.root  >& lera1_ST_tWch_t.log &
##./MakeHist_BF_mc.exe 0    20.236  35.85            998276    ./ntuples/ntNpKNUmc_ST_tW_antitop_PUMoriond17.root            hMCera1_ST_tWch_at.root >& lera1_ST_tWch_at.log &
####-----------------------------------------
##./MakeHist_BF_mc.exe 0    20.236  1975            2977600   ./ntuples/ntNpKNUmc_ZToMuMu_M_50_120_PUMoriond17.root         hMCera1_Z_M_50.root     >& lera1_Z_50.log & 
##./MakeHist_BF_mc.exe 0    20.236  19.3            100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_120_200_PUMoriond17.root        hMCera1_Z_M_120.root    >& lera1_Z_120.log &    
##./MakeHist_BF_mc.exe 0    20.236  2.73            100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_200_400_PUMoriond17.root        hMCera1_Z_M_200.root    >& lera1_Z_200.log &    
##./MakeHist_BF_mc.exe 0    20.236  0.241           98400     ./ntuples/ntNpKNUmc_ZToMuMu_M_400_800_PUMoriond17.root        hMCera1_Z_M_400.root    >& lera1_Z_400.log &    
##./MakeHist_BF_mc.exe 0    20.236  0.0168          100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_800_1400_PUMoriond17.root       hMCera1_Z_M_800.root    >& lera1_Z_800.log &    
##./MakeHist_BF_mc.exe 0    20.236  0.0008948       100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_2300_3500_PUMoriond17.root      hMCera1_Z_M_2300.root   >& lera1_Z_2300.log &    
##./MakeHist_BF_mc.exe 0    20.236  0.00004135      100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_3500_4500_PUMoriond17.root      hMCera1_Z_M_3500.root   >& lera1_Z_3500.log &    
##./MakeHist_BF_mc.exe 0    20.236  0.000000456     100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_4500_6000_PUMoriond17.root      hMCera1_Z_M_4500.root   >& lera1_Z_4500.log &    
##./MakeHist_BF_mc.exe 0    20.236  0.0000000206    100000    ./ntuples/ntNpKNUmc_ZToMuMu_M_6000_Inf_PUMoriond17.root       hMCera1_Z_M_6000.root   >& lera1_Z_6000.log &    
####------------------------------------------
##./MakeHist_BF_mc.exe 0    20.236  47.13           1000000   ./ntuples/ntNpKNUmc_WZ_PUMoriond17.root                       hMCera1_WZ.root         >& lera1_WZ.log &
##./MakeHist_BF_mc.exe 0    20.236  16.523          990064    ./ntuples/ntNpKNUmc_ZZ_PUMoriond17.root                       hMCera1_ZZ.root         >& lera1_ZZ.log &
##./MakeHist_BF_mc.exe 0    20.236   118.7          994012    ./ntuples/ntNpKNUmc_WW_PUMoriond17.root                       hMCera1_WW.root         >& lera1_WW.log &
####------------------------------------------
##./MakeHist_BF_mc.exe 0    20.236  12.178          1999000   ./ntuples/ntNpKNUmc_WWTo2L2Nu_PUMoriond17.root                hMCera1_WWTo2L2Nu.root  >& lera1_WW2L2Nu.log &
##./MakeHist_BF_mc.exe 0    20.236  0.1322          200000    ./ntuples/ntNpKNUmc_WWTo2L2Nu_Mll_200To600_PUMoriond17.root   hMCera1_WWTo2L2Nu_Mll_200.root  >& lera1_WWTo2L2Nu_200.log &
##./MakeHist_BF_mc.exe 0    20.236  0.005404        200000    ./ntuples/ntNpKNUmc_WWTo2L2Nu_Mll_600To1200_PUMoriond17.root  hMCera1_WWTo2L2Nu_Mll_600.root  >& lera1_WWTo2L2Nu_600.log &
##./MakeHist_BF_mc.exe 0    20.236  0.00033931      200000    ./ntuples/ntNpKNUmc_WWTo2L2Nu_Mll_1200To2500_PUMoriond17.root hMCera1_WWTo2L2Nu_Mll_1200.root >& lera1_WWTo2L2Nu_1200.log &
##./MakeHist_BF_mc.exe 0    20.236  0.0000051484    38969     ./ntuples/ntNpKNUmc_WWTo2L2Nu_Mll_2500ToInf_PUMoriond17.root  hMCera1_WWTo2L2Nu_Mll_2500.root >& lera1_WWTo2L2Nu_2500.log &
####------------------------------------------
##./MakeHist_BF_mc.exe 0    20.236  59.813          20000     ./ntuples/ntNpKNUmc_WprimeToMuNu_M-2400_PUMoriond17.root      hMCera1_Wprime_M-2400.root      >& lera1_Wprime2400.log &
##./MakeHist_BF_mc.exe 0    20.236  1.3926          20000     ./ntuples/ntNpKNUmc_WprimeToMuNu_M-4200_PUMoriond17.root      hMCera1_Wprime_M-4200.root      >& lera1_Wprime4200.log &
##./MakeHist_BF_mc.exe 0    20.236  0.27238          20000     ./ntuples/ntNpKNUmc_WprimeToMuNu_M-1800_PUMoriond17.root      hMCera1_Wprime_M-1800.root      >& lera1_Wprime1800.log &
##./MakeHist_BF_mc.exe 0    20.236  0.0029432        20000     ./ntuples/ntNpKNUmc_WprimeToMuNu_M-3800_PUMoriond17.root      hMCera1_Wprime_M-3800.root      >& lera1_Wprime3800.log &
