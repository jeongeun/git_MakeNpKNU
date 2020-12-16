import FWCore.ParameterSet.Config as cms

ProdNpGenerator = cms.EDProducer('ProdNpGenerator',
	DoGenParticle = cms.untracked.bool(True),
	genParticles = cms.InputTag('prunedGenParticles'),
	PrintParticle = cms.untracked.bool(False),
	PrintGenParticleNum = cms.untracked.bool(False),
	DoGenInfo    = cms.untracked.bool(True),
        genInfoToken = cms.InputTag('generator'),
)

