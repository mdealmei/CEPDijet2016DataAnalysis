import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
        'file:ggTo2Jet_13TeV_noPU_TuneCUETP8M1_Realistic50ns13TeVCollision.root'
    )
)

# Adding Good Primary Vertex
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

''' 
# Adding Tracks Associator with Vertex Collection
process.ak4JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
       tracks = cms.InputTag("generalTracks"),
       jets = cms.InputTag("ak4PFJets"),
       coneSize = cms.double(0.4)
)
'''
process.demo = cms.EDAnalyzer('AnalyzerPPS'
		, genParticle = cms.InputTag('genParticles')
		, RunWithWeightGen = cms.bool(False)
		, genjets = cms.InputTag("ak4GenJets")
		, vertex      = cms.InputTag('goodOfflinePrimaryVertices')
		, track = cms.InputTag('generalTracks')
		, jets      = cms.InputTag('ak4PFJetsCHS')
		, PFCandidates = cms.InputTag("particleFlow")
		, muons      = cms.InputTag('muons')
		, electrons      = cms.InputTag('gedGsfElectrons')
		, ppsGen   = cms.InputTag('ppssim:PPSGen')
		, ppsSim   = cms.InputTag('ppssim:PPSSim')
		, ppsReco   = cms.InputTag('ppssim:PPSReco')
		, pTPFThresholdCharged = cms.double(0.1)
                , energyPFThresholdBar = cms.double(1.5)
                , energyPFThresholdEnd = cms.double(3.5)
                , energyPFThresholdHF = cms.double(4.0)
		, EBeam = cms.double(13000)
)

process.TFileService = cms.Service("TFileService",
                fileName = cms.string('AnalyzerPPS.root')
)


#process.p = cms.Path(process.goodOfflinePrimaryVertices*process.ak4JetTracksAssociatorAtVertex*process.demo)
process.p = cms.Path(process.goodOfflinePrimaryVertices*process.demo)
