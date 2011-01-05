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

process.GlobalTag.globaltag = cms.string("GR10_P_V6::All")

process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")

process.DQMStore = cms.Service("DQMStore")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/04359DC0-5177-DF11-82A8-0030486792BA.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/041EA906-5377-DF11-A6F1-003048679000.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/041E9CA9-4E77-DF11-AD33-0030486791AA.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/027C15BB-4E77-DF11-B11B-002618943833.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/00CF8BFF-5177-DF11-ABA1-00261894396F.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/00570C24-5277-DF11-92DC-003048678B38.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/00346AE4-5177-DF11-B2F9-001A92971B28.root',
    '/store/data/Commissioning10/MinimumBias/RECO/GOODCOLL-Jun9thSkim_v1/0032/0019751C-5277-DF11-8C08-0018F3D096E8.root'

    )
                            

)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

    
## Messages & Convenience
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
reportEvery = cms.untracked.int32(100)     # every 1000th only
#    limit = cms.untracked.int32(10)       # or limit to 10 printouts...
))
process.MessageLogger.statistics.append('cout')
      


process.fileSaver = cms.EDAnalyzer("METFileSaver",
    OutputFile = cms.untracked.string('data.root') )
process.p = cms.Path(process.fileSaver*
                     process.calotoweroptmaker*
                     process.calotoweroptmakerWithHO*
                     process.towerMakerWithHO*
                     #process.genJetParticles*
                     #process.genMETParticles*
                     process.metreco*
                     process.analyzeRecHits*
                     process.analyzecaloTowers*
                     #process.analyzeGenMET*
                     #process.analyzeGenMETFromGenJets*
                     process.analyzeHTMET*
                     process.analyzeCaloMET*
                     process.analyzePFMET*
                     process.analyzeTCMET*
                     process.analyzeMuonCorrectedCaloMET*
                     process.AnalyzeBeamHalo
)
process.schedule = cms.Schedule(process.p)


