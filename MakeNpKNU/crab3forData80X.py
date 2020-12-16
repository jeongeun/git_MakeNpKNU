from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'ntReminiaod_RunH_PR_v2'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'NtMaker_Reminiaod'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName ='/hcp/data/data02/jelee/Wprime13TeV/Moriond17/CMSSW_8_0_26_patch2/src/MakeNpKNU/makeNpKNU_80X_reMiniAOD_runH.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['ntNpKNUdata.root','histFileData.root']
config.JobType.maxMemoryMB = 2500
config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' #1115
config.Data.runRange ='280919-284044' 
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/jelee/'

config.section_('User')
config.section_("Site")
config.Site.storageSite = 'T2_KR_KNU'
