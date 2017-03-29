import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#        'file:samples/miniAOD_GluGluTo2Jets_noPU.root'
	 'file:samples/miniAOD_GluGluTo2Jets_noPU_20k.root'
    )

)


''' 
# Adding Tracks Associator with Vertex Collection
process.ak4JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
       tracks = cms.InputTag("generalTracks"),
       jets = cms.InputTag("ak4PFJets"),
       coneSize = cms.double(0.4)
)
'''

# PPS Vertex Filter
process.ppsVertexFilter = cms.EDFilter( 'ppsVertexFilter'
        , ppsReco               = cms.InputTag('ppssim:PPSReco')
        , vertex                = cms.InputTag('offlineSlimmedPrimaryVertices')
        , ppsZ_resolution       = cms.double(0.3)
)

# PPS Golden Vertex
process.goldenVertex = cms.EDFilter('goldenVertex'
        , ppsReco               = cms.InputTag('ppssim:PPSReco')
        , vertex                = cms.InputTag('offlineSlimmedPrimaryVertices')
)

# FatJets Filter
process.fatjets = cms.EDFilter('fatjets'
	, jets  = cms.InputTag('slimmedJets')
	, wideJetDeltaR = cms.double(1.1)
)

process.demo = cms.EDAnalyzer('AnalyzerPPS'
		, genParticle = cms.InputTag('genParticles')
		, RunWithWeightGen = cms.bool(False)
		, genjets = cms.InputTag('slimmedGenJets')
		, vertex      = cms.InputTag('offlineSlimmedPrimaryVertices')
		, jets      = cms.InputTag('slimmedJets')
		, widejetsTag = cms.untracked.InputTag("fatjets","widejets")
		, ppsGen   = cms.InputTag('ppssim:PPSGen')
		, ppsSim   = cms.InputTag('ppssim:PPSSim')
		, ppsReco   = cms.InputTag('ppssim:PPSReco')
		, EBeam = cms.double(13000)
		#, PPSres = cms.double(0.3)
)


process.ntuplesOut = cms.OutputModule(
     "PoolOutputModule",
     outputCommands = cms.untracked.vstring(
     "keep *_fatjets_*_*"       
     )
)


process.TFileService = cms.Service("TFileService",
                fileName = cms.string('AnalyzerPPS_noPU.root')
)




#process.p = cms.Path(process.goldenVertex*process.demo)
#process.p = cms.Path(process.demo)
#process.p = cms.Path(process.ppsVertexFilter*process.demo)
#process.p = cms.Path(process.fatjets)
#process.p = cms.Path(process.goldenVertex)
process.p = cms.Path(process.goldenVertex*process.fatjets*process.demo)


#process.ntpoutpath = cms.EndPath(process.ntuplesOut)


