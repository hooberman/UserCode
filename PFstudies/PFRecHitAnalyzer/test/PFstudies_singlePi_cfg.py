import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES2")

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
    'file:/tas05/disk00/benhoob/tcmetTestFiles/singePi_0_2/reco_100_1_gc4.root'
    )
                            )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/tas05/disk00/benhoob/tcmetTestFiles/output/PFstudies_singlePi_0_2.root')
                               )

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Analyzer*_*_OWNPARTICLES2*'))
 
#process.p = cms.Path(process.tcmetSequence)
process.p = cms.Path(process.MYpfRecHitSequence)

process.e = cms.EndPath(process.out)
