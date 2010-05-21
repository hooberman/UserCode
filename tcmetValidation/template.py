import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("RecoMET.Configuration.CaloTowersOptForMET_cff")

process.load("RecoMET.Configuration.RecoMET_cff")

process.load("RecoMET.Configuration.RecoHTMET_cff")

process.load("RecoMET.Configuration.RecoGenMET_cff")

process.load("RecoMET.Configuration.GenMETParticles_cff")

process.load("RecoMET.Configuration.RecoPFMET_cff")

process.load("RecoJets.Configuration.CaloTowersRec_cff")

process.load("Validation.RecoMET.CaloMET_cff")

process.load("Validation.RecoMET.GenMET_cff")

process.load("Validation.RecoMET.HTMET_cff")

process.load("Validation.RecoMET.GenMETFromGenJets_cff")
process.load("RecoMET/Configuration/RecoMET_BeamHaloId_cff")
process.DQMStore = cms.Service("DQMStore")
process.load("DQMOffline/JetMET/BeamHaloAnalyzer_cfi")

process.load("DQMOffline.JetMET.caloTowers_cff")
process.towerSchemeBAnalyzer.FineBinning = cms.untracked.bool(True)
process.towerSchemeBAnalyzer.FolderName =  cms.untracked.string("RecoMETV/MET_CaloTowers/SchemeB")
process.towerOptAnalyzer.FineBinning = cms.untracked.bool(True)
process.towerOptAnalyzer.FolderName =  cms.untracked.string("RecoMETV/MET_CaloTowers/Optimized")

process.load("DQMOffline.JetMET.RecHits_cff")
process.ECALAnalyzer.FineBinning = cms.untracked.bool(True)
process.ECALAnalyzer.FolderName =  cms.untracked.string("RecoMETV/MET_ECAL/data")
process.HCALAnalyzer.FineBinning = cms.untracked.bool(True)
process.HCALAnalyzer.FolderName =  cms.untracked.string("RecoMETV/MET_HCAL/data")

process.load("Validation.RecoMET.PFMET_cff")

process.load("Validation.RecoMET.TCMET_cff")

process.load("Validation.RecoMET.MuonCorrectedCaloMET_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = cms.string("%(MYGLOBALTAG)s")   

process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")

process.DQMStore = cms.Service("DQMStore")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    %(MYFILENAMES)s

    )
                            

)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(%(MYNEVENTS)s) )

    
## Messages & Convenience
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
reportEvery = cms.untracked.int32(100)     # every 1000th only
#    limit = cms.untracked.int32(10)       # or limit to 10 printouts...
))
process.MessageLogger.statistics.append('cout')
      


process.fileSaver = cms.EDAnalyzer("METFileSaver",
    OutputFile = cms.untracked.string('%(MYROOTFILE)s') )
process.p = cms.Path(process.fileSaver*
                     process.calotoweroptmaker*
                     process.calotoweroptmakerWithHO*
                     process.towerMakerWithHO*
                     process.genJetParticles*
                     process.genMETParticles*
                     process.metreco*
                     process.analyzeRecHits*
                     process.analyzecaloTowers*
                     process.analyzeGenMET*
                     process.analyzeGenMETFromGenJets*
                     process.analyzeHTMET*
                     process.analyzeCaloMET*
                     process.analyzePFMET*
                     process.analyzeTCMET*
                     process.analyzeMuonCorrectedCaloMET*
                     process.AnalyzeBeamHalo
)
process.schedule = cms.Schedule(process.p)


