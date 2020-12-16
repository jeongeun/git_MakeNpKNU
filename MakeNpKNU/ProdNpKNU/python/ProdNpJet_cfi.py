import FWCore.ParameterSet.Config as cms

ProdNpJet = cms.EDProducer('ProdNpJet',
   jets = cms.InputTag("slimmedJets"),
	MinPtCut = cms.untracked.double(10.0),
)

