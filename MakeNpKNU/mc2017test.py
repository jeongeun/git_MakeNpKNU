import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = dict()
varOptions = VarParsing('analysis')
varOptions.register("isMC", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.parseArguments()

process = cms.Process("MakeNpKNU")

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

### Message Logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
	allowUnscheduled = cms.untracked.bool(True),
	wantSummary = cms.untracked.bool(True),
        SkipEvent = cms.untracked.vstring('ProductNotFound')
)

runOnData=False
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

if (varOptions.isMC):
	process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'##reMiniaodJEC Regression included# '80X_mcRun2_asymptotic_2016_TrancheIV_v6'##'80X_mcRun2_asymptotic_2016_TrancheIV_v7' ##For Moriond17 MC (2017/Jan) #EGMRegression(2017)
        inputFileName = '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root'
	ntFileName = "ntNpKNUmc.root"
	histFileName = 'histFileMC.root'
else:
	process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7' #B-G(Feb03) #80X_dataRun2_2016SeptRepro_v6'  #'80X_dataRun2_Prompt_ICHEP16JEC_v0'#EGMRegression(2017)
        inputFileName = '/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/FEF25E85-82EC-E611-A8F0-0CC47A4D7670.root'
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

process.EGMSmearerElectrons = cms.Path(process.calibratedPatElectrons)
process.calibratedPatElectrons.isMC = cms.bool(True) #varOptions.isMC
print "HEEPV70"

from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupVIDElectronSelection
from PhysicsTools.SelectorUtils.tools.vid_id_tools import setupAllVIDIdsInModule
from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
	'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.heepIDVarValueMaps.elesMiniAOD         = cms.InputTag("slimmedElectrons")
print "MET Filter"

process.load("RecoMET.METFilters.BadPFMuonFilter_cfi")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates= cms.InputTag("packedPFCandidates")
process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
print "MET Correction"
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

process.p = cms.Path(
# process.fullPatMetSequence ## no puppi
#* process.BadPFMuonFilter 
#* process.BadChargedCandidateFilter 
#* process.regressionApplication 
#* process.calibratedPatElectrons 
 process.egmGsfElectronIDSequence #* process.egmPhotonIDSequence
* process.ProdNpElectron
* process.MakeNpKNU
)
