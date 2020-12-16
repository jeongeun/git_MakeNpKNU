import FWCore.ParameterSet.Config as cms

ProdNpPhoton = cms.EDProducer('ProdNpPhoton',
        photons = cms.InputTag("slimmedPhotons"),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
        phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
	genParticles = cms.InputTag("prunedGenParticles"),
        electrons = cms.InputTag("slimmedElectrons"),
        cutBasedIdTokenTight  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
        cutBasedIdTokenMedium = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
        cutBasedIdTokenLoose  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
	mvaMediumIdMap         = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"),
	mvaMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"),
	mvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
	mvaCategoriesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Categories"),
	MinPtCut = cms.untracked.double(10.0),
	PrintNum = cms.untracked.bool(False),
)

