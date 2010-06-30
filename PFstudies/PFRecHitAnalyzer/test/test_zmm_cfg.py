import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START3X_V26::All"

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PFstudies.PFRecHitAnalyzer.pfrechitanalyzer_cff")
#process.load("RecoMET.METProducers.TCMET_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/users/fgolf/devel/zmm.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('zmm.root')
)

#process.out.outputCommands = cms.untracked.vstring( 'drop *' )
#process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Analyzer*_*_OWNPARTICLES*'))
 
process.p = cms.Path(process.tcmetSequence)
1
process.e = cms.EndPath(process.out)
