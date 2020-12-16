import FWCore.ParameterSet.Config as cms

ProdNpMisc = cms.EDProducer('ProdNpMisc',
	DoPileup = cms.untracked.bool(True),
        pileupToken = cms.InputTag("slimmedAddPileupInfo"),
	DoVertex = cms.untracked.bool(True),
        vertexToken = cms.InputTag("offlineSlimmedPrimaryVertices"),
	DoTrigger = cms.untracked.bool(True),
        triggerBitsToken = cms.InputTag("TriggerResults","","HLT"),
        triggerObjectsToken = cms.InputTag("selectedPatTrigger"),
        triggerPrescalesToken = cms.InputTag("patTrigger"),
        triggerIdentifiers = cms.vstring(['HLT_Ele*','HLT_Mu*','HLT_TkMu*']),
        PrintTriggerResult = cms.untracked.bool(False),
)

