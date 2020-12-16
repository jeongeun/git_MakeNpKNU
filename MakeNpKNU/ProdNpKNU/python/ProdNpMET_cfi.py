import FWCore.ParameterSet.Config as cms

ProdNpMET = cms.EDProducer('ProdNpMET',
         metToken = cms.InputTag("slimmedMETs","","PAT"),
         metEGCleanToken = cms.InputTag("slimmedMETsEGClean","","PAT"),
         metMuEGCleanToken = cms.InputTag("slimmedMETsMuEGClean","","PAT"),
         metMuEGCleanCorrToken = cms.InputTag("slimmedMETsMuEGClean","","MakeNpKNU"),
         metUncorrectedToken = cms.InputTag("slimmedMETsUncorrected"),
         metPUPPIToken   = cms.InputTag("slimmedMETsPuppi"),
	 BadPFMuonFilter = cms.InputTag("BadPFMuonFilter",""),
	 BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
	 #badGlobalMuonTagger = cms.InputTag("BadGlobalMuonTagger",""),
         #badGlobalMuonFilter = cms.InputTag("badGlobalMuonTagger","bad"),
         #duplicateMuonFilter = cms.InputTag("cloneGlobalMuonTagger","bad"),
	 metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),
	 metFilterBits_data = cms.InputTag("TriggerResults", "", "RECO"),
	 particleFlowEGammaGSFixed = cms.InputTag("particleFlowEGammaGSFixed", "dupECALClusters", "PAT"),
         ecalMultiAndGSGlobalRecHitEB = cms.InputTag("ecalMultiAndGSGlobalRecHitEB","hitsNotReplaced", "PAT"),
         MinPtCut = cms.untracked.double(800.0),
  	 PrintNum = cms.untracked.bool(True),
)

