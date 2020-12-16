import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = dict()
varOptions = VarParsing('analysis')
varOptions.register("isMC", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "isMC" )
varOptions.parseArguments()

process = cms.Process("MakeNpKNU")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

if (varOptions.isMC):
	process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
	inputFileName = 'root://cms-xrd-global.cern.ch//store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/002ABFCA-A0B9-E511-B9BA-0CC47A57CD6A.root'
	ntFileName = "ntNpKNUmc.root"
	histFileName = 'histFileMC.root'
else:
	process.GlobalTag.globaltag = '76X_dataRun2_v15'
	inputFileName = 'root://cms-xrd-global.cern.ch//store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/50000/9C88480E-79A6-E511-92D4-0CC47A4D769E.root'
	ntFileName = 'ntNpKNUdata.root'
	histFileName = 'histFileData.root'

print "GlobalTag:", process.GlobalTag.globaltag
print "inFileName:", inputFileName
print "ntFileName:", ntFileName
print "histFileName:", histFileName

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFileName))
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(ntFileName) )

####################
# For Electron ID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff'
]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# For Photon ID.
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, dataFormat)
my_photonID_modules = [
'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff',
]
for idmod in my_photonID_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#start MakeNpKNU	
from MakeNpKNU.ProdNpKNU.ProdNpElectron_cfi import *
process.ProdNpElectron = ProdNpElectron.clone()
process.ProdNpElectron.DoMakeHist = cms.untracked.bool(True)
process.ProdNpElectron.electronIDTags = cms.VInputTag(
"egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto",
"egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose",
"egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium",
"egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight",
"egmGsfElectronIDs:heepElectronID-HEEPV60",
)
process.ProdNpElectron.mvaMediumIdMap   =  cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90")
process.ProdNpElectron.mvaTightIdMap    =  cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80")
process.ProdNpElectron.mvaValuesMap     =  cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values")
process.ProdNpElectron.mvaCategoriesMap =  cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories")

from MakeNpKNU.ProdNpKNU.ProdNpPhoton_cfi import *
process.ProdNpPhoton = ProdNpPhoton.clone()
process.ProdNpPhoton.DoMakeHist = cms.untracked.bool(True)
process.ProdNpPhoton.photonIDTags = cms.VInputTag(
"egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose",
"egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium",
"egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight",
)
process.ProdNpPhoton.mvaMediumIdMap         = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90")
process.ProdNpPhoton.mvaMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90")
process.ProdNpPhoton.mvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values")
process.ProdNpPhoton.mvaCategoriesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Categories")

from MakeNpKNU.ProdNpKNU.ProdNpMuon_cfi import *
process.ProdNpMuon = ProdNpMuon.clone()

from MakeNpKNU.ProdNpKNU.ProdNpMuon_cfi import *
process.ProdNpMuon = ProdNpMuon.clone()

from MakeNpKNU.ProdNpKNU.ProdNpJet_cfi import *
process.ProdNpJet = ProdNpJet.clone()

from MakeNpKNU.ProdNpKNU.ProdNpMET_cfi import *
process.ProdNpMET = ProdNpMET.clone()

from MakeNpKNU.ProdNpKNU.ProdNpGenerator_cfi import *
process.ProdNpGenerator = ProdNpGenerator.clone()

from MakeNpKNU.ProdNpKNU.ProdNpMisc_cfi import *
process.ProdNpMisc = ProdNpMisc.clone()
process.ProdNpMisc.triggerIdentifiers = cms.vstring(['HLT_Ele*','HLT_Mu*'])

process.MakeNpKNU = cms.EDAnalyzer('MakeNpKNU',
   electronToken = cms.InputTag('ProdNpElectron','Electron'),
   muonToken = cms.InputTag('ProdNpMuon','Muon'),
   photonToken = cms.InputTag('ProdNpPhoton','Photon'),
   jetToken = cms.InputTag('ProdNpJet','Jet'),
   metToken = cms.InputTag('ProdNpMET','MET'),
   vertexToken = cms.InputTag('ProdNpMisc','Vertex'),
   triggerToken = cms.InputTag('ProdNpMisc','Trigger'),
   triggerObjectToken = cms.InputTag('ProdNpMisc','TriggerObject'),
   pileupToken = cms.InputTag('ProdNpMisc','Pileup'),
   genParticleToken = cms.InputTag('ProdNpGenerator','GenParticle'),
   genInfoToken = cms.InputTag('ProdNpGenerator','GenInfo'),
   HistFileName = cms.string(histFileName),
)

process.p = cms.Path(
  process.egmGsfElectronIDSequence
* process.egmPhotonIDSequence
* process.ProdNpElectron
* process.ProdNpMuon
* process.ProdNpPhoton
* process.ProdNpMET
* process.ProdNpJet
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


