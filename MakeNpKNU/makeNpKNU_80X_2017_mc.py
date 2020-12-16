import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = dict()
varOptions = VarParsing('analysis')
varOptions.register("isMC", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.parseArguments()

process = cms.Process("MakeNpKNU")

### Standard modules Loading
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

### Message Logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# SetOption -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
	allowUnscheduled = cms.untracked.bool(True),
	wantSummary = cms.untracked.bool(True),
        SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#configurable options ============================================
runOnData=False
#usePrivateSQlite=True #use external JECs (sqlite file)
#useHFCandidates=False #create an additionnal NoHF slimmed MET collection if the option is set to false
#redoPuppi=True # rebuild puppiMETi
#DoReclusterMET=True #re-cluster and get the proper uncertainties
#redoPuppiMET=False # rebuild puppiMET
#METTagging=True
DEBUG=False
# ================================================================

## DEBUG Tracer
if (DEBUG):
        process.Tracer = cms.Service("Tracer")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if (varOptions.isMC):
	process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'##reMiniaodJEC Regression included# '80X_mcRun2_asymptotic_2016_TrancheIV_v6'##'80X_mcRun2_asymptotic_2016_TrancheIV_v7' ##For Moriond17 MC (2017/Jan) #EGMRegression(2017)
        inputFileName = '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root'
        #inputFileName = 'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/WToMuNu_M-1000_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/40000/021B26DA-5630-E611-AA8E-90B11C2AA16C.root'
	ntFileName = "ntNpKNUmc.root"
	histFileName = 'histFileMC.root'
else:
	process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' #B-G(Feb03) #80X_dataRun2_2016SeptRepro_v6'  #'80X_dataRun2_Prompt_ICHEP16JEC_v0'#EGMRegression(2017)
#	process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v16' #H(Feb03) #80X_dataRun2_2016SeptRepro_v6'  #'80X_dataRun2_Prompt_ICHEP16JEC_v0'#EGMRegression(2017)
        inputFileName = '/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/FEF25E85-82EC-E611-A8F0-0CC47A4D7670.root'
        #inputFileName = '/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/086A6153-A6EA-E611-B985-001E67A3F92F.root'
        #inputFileName = 'root://cms-xrdr.sdfarm.kr:1094//store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/FEF25E85-82EC-E611-A8F0-0CC47A4D7670.root'
	#inputFileName = 'file:/hcp/data/data01/ycyang/store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/A2C0F697-B19C-E611-A4D8-F04DA275BFF2.root'
	#inputFileName = '/store/data/Run2016G/SingleMuon/MINIAOD/PromptReco-v1/000/278/820/00000/0667AC34-2464-E611-84CE-02163E011979.root'
	#inputFileName = 'root://cms-xrd-global.cern.ch//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/A2C0F697-B19C-E611-A4D8-F04DA275BFF2.root'
	#inputFileName = 'root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/A8816A41-1E1A-E611-907B-02163E0128EB.root'
	ntFileName = 'ntNpKNUdata.root'
	histFileName = 'histFileData.root'

print "GlobalTag:", process.GlobalTag.globaltag
print "inFileName:", inputFileName
print "ntFileName:", ntFileName
print "histFileName:", histFileName

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFileName))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.TFileService = cms.Service("TFileService",
	closeFileFast = cms.untracked.bool(True),	
	fileName = cms.string(ntFileName) )

###### LOAD DATABASE
#from CondCore.CondDB.CondDB_cfi import *
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename='./etc/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()
#############

######################################
# 2017  EGM Smear/Scale, Regression
######################################
#https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression#TrainingBoundary (17/2/27)
#https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer
print "EGM smearing"
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                  calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(8675389),
                                                      engineName = cms.untracked.string('TRandom3'),
                                                      ),
                                                   )
process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')

#process.selectedElectrons = cms.EDFilter("PATElectronSelector",
#    src = cms.InputTag("slimmedElectrons"),
#    cut = cms.string("pt > 5 && abs(eta)<2.5")
#)
# Path and EndPath definitions
#process.EGMRegression = cms.Path(process.regressionApplication)
process.EGMSmearerElectrons = cms.Path(process.calibratedPatElectrons)
process.calibratedPatElectrons.isMC = cms.bool(True) #varOptions.isMC
#process.calibratedPatElectrons.electrons = "selectedElectrons"


# Set the lines below to True or False depending if you are correcting the scale (data) or smearing the resolution (MC) 
#process.EGMSmearerElectrons.isMC = cms.untracked.bool(False)  #varOptions.isMC # cms.bool(False)

###https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression (17/01/15)
###To read from the DB, first tell your code where to find the regression.
print "HEEPV70"

#####################
# For EGM ID (New HEEPV70)
#####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupVIDElectronSelection
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupAllVIDIdsInModule
from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat
# turn on VID producer, indicate data format  to be eleProducer = "PatElectronSelectorByValueMap"
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
#	'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
	'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#	'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff'
]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.heepIDVarValueMaps.elesMiniAOD         = cms.InputTag("slimmedElectrons")
print "MET Filter"

######################
# Including METFilter 
##################### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
#process.load("RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff")
#process.badGlobalMuonTagger.muons= cms.InputTag("slimmedMuons")
#process.badGlobalMuonTagger.vtx = cms.InputTag("offlineSlimmedPrimaryVertices")

process.load("RecoMET.METFilters.BadPFMuonFilter_cfi")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates= cms.InputTag("packedPFCandidates")
#process.BadPFMuonFilter.debug = cms.bool(True)
#process.BadPFMuonFilter.taggingMode=False,
		
process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadChargedCandidateFilter.debug = cms.bool(True)
#process.BadChargedCandidateFilter.taggingMode=False,

###########################
#  Recompute Corrected MET 
###########################
#https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py#L136
#default configuration for miniAOD reprocessing, change the isData flag to run on data for a full met computation, remove the pfCandColl input
#https://github.com/MiT-HEP/NeroProducer/blob/master/Nero/test/testNero.py
###############################################
print "MET Correction"
#  2017 re-miniaod Recompute Corrected MET 
###############################################
#https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#MET_Recipes (17/3/3)
# How to re-correct MET based on bad muons on the fly for MC (Release: CMSSW_8_0_X, X>=26_patch1) 
### Following lines are for default MET for Type1 corrections.
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#
## If you only want to re-correct for JEC and get the proper uncertainties for the default MET
#runMetCorAndUncFromMiniAOD(process,
#                       isData= False # (or False),
#                      )
#
## Now you are creating the bad muon corrected MET
#process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
#process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
#process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
#
#from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation
#
#muonRecoMitigation(
#                   process = process,
#                   pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
#                   runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
#                   selection="", #You can use a custom selection for your bad muons. Leave empty if you would like to use the bad muon recipe definition.
#                   muonCollection="", #The muon collection name where your custom selection will be applied to. Leave empty if you would like to use the bad muon recipe definition.
#                   cleanCollName="cleanMuonsPFCandidates", #output pf candidate collection ame
#                   cleaningScheme="computeAllApplyClone", #Options are: "all", "computeAllApplyBad","computeAllApplyClone". Decides which (or both) bad muon collections to be used for MET cleaning coming from the bad muon recipe.
#                   postfix="" #Use if you would like to add a post fix to your muon / pf collections
#                   )
#
#runMetCorAndUncFromMiniAOD(process,
#                          isData=False,
#                          pfCandColl="cleanMuonsPFCandidates",
#                          recoMetFromPFCs=True,
#                          postfix="MuClean"
#                          )
##del process.slimmedMETsMuEGClean.caloMET
#
## Repeated calls to runMetCorAnd.. overwrites the MET source of patCaloMet
#process.patCaloMet.metSource = 'metrawCaloPuppi'
#
#process.mucorMET = cms.Sequence(                     
#    process.badGlobalMuonTaggerMAOD *
#    process.cloneGlobalMuonTaggerMAOD *
#    #process.badMuons * # If you are using cleaning mode "all", uncomment this line
#    process.cleanMuonsPFCandidates *
#    process.fullPatMetSequenceMuClean
#    )
#
#

###How to re-correct MET based on e/gamma gain switch correction on the fly for Re-Miniaod Data (Release: CMSSW_8_0_X, X>=26_patch1) 

### Following lines are for default MET for Type1 corrections.
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#
## If you only want to re-correct for JEC and get the proper uncertainties for the default MET
#runMetCorAndUncFromMiniAOD(process,
#                       isData= False #(or False),
#                      )
#
## Now you are creating the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
#from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
#corMETFromMuonAndEG(process,
#                 pfCandCollection="", #not needed                                                                                                                                                                                                                                                                                                                      
#                 electronCollection="slimmedElectronsBeforeGSFix",
#                 photonCollection="slimmedPhotonsBeforeGSFix",
#                 corElectronCollection="slimmedElectrons",
#                 corPhotonCollection="slimmedPhotons",
#                 allMETEGCorrected=True,
#                 muCorrection=False,
#                 eGCorrection=True,
#                 runOnMiniAOD=True,
#                 postfix="MuEGClean"
#                 )
#process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
#process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
#process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
#process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
#del process.slimmedMETsMuEGClean.caloMET
#
## If you are running in the scheduled mode:
#process.egcorrMET = cms.Sequence(
#   process.cleanedPhotonsMuEGClean+process.cleanedCorPhotonsMuEGClean+
#   process.matchedPhotonsMuEGClean + process.matchedElectronsMuEGClean +
#   process.corMETPhotonMuEGClean+process.corMETElectronMuEGClean+
#   process.patPFMetT1MuEGClean+process.patPFMetRawMuEGClean+
#   process.patPFMetT1SmearMuEGClean+process.patPFMetT1TxyMuEGClean+
#   process.patPFMetTxyMuEGClean+process.patPFMetT1JetEnUpMuEGClean+
#   process.patPFMetT1JetResUpMuEGClean+process.patPFMetT1SmearJetResUpMuEGClean+
#   process.patPFMetT1ElectronEnUpMuEGClean+process.patPFMetT1PhotonEnUpMuEGClean+
#   process.patPFMetT1MuonEnUpMuEGClean+process.patPFMetT1TauEnUpMuEGClean+
#   process.patPFMetT1UnclusteredEnUpMuEGClean+process.patPFMetT1JetEnDownMuEGClean+
#   process.patPFMetT1JetResDownMuEGClean+process.patPFMetT1SmearJetResDownMuEGClean+
#   process.patPFMetT1ElectronEnDownMuEGClean+process.patPFMetT1PhotonEnDownMuEGClean+
#   process.patPFMetT1MuonEnDownMuEGClean+process.patPFMetT1TauEnDownMuEGClean+
#   process.patPFMetT1UnclusteredEnDownMuEGClean+process.slimmedMETsMuEGClean)
#

#start MakeNpKNU	
from MakeNpKNU.ProdNpKNU.ProdNpElectron_cfi import *
process.ProdNpElectron = ProdNpElectron.clone()
process.ProdNpElectron.DoMakeHist = cms.untracked.bool(True)
process.ProdNpElectron.electronIDTags = cms.VInputTag(
	"egmGsfElectronIDs:heepElectronID-HEEPV70",
)
# ProdNpMuon DoMakeHist 
from MakeNpKNU.ProdNpKNU.ProdNpMuon_cfi import *
process.ProdNpMuon = ProdNpMuon.clone()
process.ProdNpMuon.DoMakeHist = cms.untracked.bool(False)

# ProdNpJet 
from MakeNpKNU.ProdNpKNU.ProdNpJet_cfi import *
process.ProdNpJet = ProdNpJet.clone()

# ProdNpMET DoMakeHist
from MakeNpKNU.ProdNpKNU.ProdNpMET_cfi import *
process.ProdNpMET = ProdNpMET.clone()
process.ProdNpMET.DoMakeHist = cms.untracked.bool(False)
process.ProdNpMET.metscleanmu=cms.InputTag('slimmedMETsMuClean','','MakeNpKNU')
process.ProdNpMET.metCorrToken      = cms.InputTag('slimmedMETs','','MakeNpKNU')
process.ProdNpMET.metPUPPICorrToken = cms.InputTag('slimmedMETsPuppi','','MakeNpKNU')
BadMuonFilter = cms.InputTag("BadPFMuonFilter","")
BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter","")
# ProdNpGenerator 
from MakeNpKNU.ProdNpKNU.ProdNpGenerator_cfi import *
process.ProdNpGenerator = ProdNpGenerator.clone()

# ProdNpMisc(Pileup, Vertex, Trigger)
from MakeNpKNU.ProdNpKNU.ProdNpMisc_cfi import *
process.ProdNpMisc = ProdNpMisc.clone()
process.ProdNpMisc.triggerIdentifiers = cms.vstring(['HLT_Ele*','HLT_Mu*','HLT_TkMu*','HLT_ECALHT*','HLT_Photon175*','HLT_PFMET*','HLT_MET*'])
process.ProdNpMisc.DoMakeHist = cms.untracked.bool(False)

process.MakeNpKNU = cms.EDAnalyzer('MakeNpKNU',
     electronToken     = cms.InputTag('ProdNpElectron','Electron'),
     muonToken         = cms.InputTag('ProdNpMuon','Muon'),
 #    photonToken       = cms.InputTag('ProdNpPhoton','Photon'),
     jetToken          = cms.InputTag('ProdNpJet','Jet'),
     metToken          = cms.InputTag('ProdNpMET','MET'),
     vertexToken       = cms.InputTag('ProdNpMisc','Vertex'),
     triggerToken      = cms.InputTag('ProdNpMisc','Trigger'),
     triggerObjectToken= cms.InputTag('ProdNpMisc','TriggerObject'),
     pileupToken       = cms.InputTag('ProdNpMisc','Pileup'),
     genParticleToken  = cms.InputTag('ProdNpGenerator','GenParticle'),
     genInfoToken      = cms.InputTag('ProdNpGenerator','GenInfo'),
     HistFileName      = cms.string(histFileName),
)

#DEBUG
#print "Process=",process, process.__dict__.keys()
#------------------------------------------------------



process.p = cms.Path(
# process.fullPatMetSequence ## no puppi
#* process.egcorrMET 
#* process.mucorMET
#* process.BadPFMuonFilter 
#* process.BadChargedCandidateFilter 
#* process.regressionApplication 
#* process.calibratedPatElectrons 
#* process.selectedElectrons
#* process.EGMRegression
#* process.EGMSmearerElectrons
 process.egmGsfElectronIDSequence #* process.egmPhotonIDSequence
#* process.puppiMETSequence ##puppi
#* process.fullPatMetSequencePuppi ##puppi
* process.ProdNpElectron
* process.ProdNpMuon
#* process.ProdNpPhoton
* process.ProdNpJet
* process.ProdNpMET
* process.ProdNpGenerator
* process.ProdNpMisc
* process.MakeNpKNU
)
#####################

# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string('edmMakeNpKNU.root')
#     ,outputCommands = cms.untracked.vstring('keep *')
# )
# process.e = cms.EndPath(process.out)


