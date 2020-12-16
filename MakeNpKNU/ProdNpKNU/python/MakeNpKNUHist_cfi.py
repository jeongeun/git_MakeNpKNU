import FWCore.ParameterSet.Config as cms

MakeNpKNUHist =  cms.EDAnalyzer('MakeNpKNUHist', 
        pileup = cms.InputTag("slimmedAddPileupInfo"),
)

