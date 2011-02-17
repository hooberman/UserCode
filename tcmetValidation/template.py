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

process.load("Validation.RecoMET.TCMET_cfi")
process.tcMetAnalyzer.sample = cms.untracked.string("%(MYSAMPLE)s")
process.tcMetAnalyzer.InputMETLabel = cms.InputTag("%(METCOLLECTION)s")
                      
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
      
from RecoMET.Configuration.RecoMET_cff import metreco as metreco2
from RecoJets.Configuration.CaloTowersES_cfi import *
from RecoMET.METProducers.CaloTowersOpt_cfi import *
from RecoMET.METProducers.CaloMET_cfi import *
from RecoMET.METProducers.HTMET_cfi import *
from RecoMET.METProducers.CaloMETSignif_cfi import *
#from RecoMET.METProducers.TCMET_cfi import *
from RecoMET.METProducers.hcalnoiseinfoproducer_cfi import *
from RecoMET.METProducers.MuonMETValueMapProducer_cff import *
from RecoMET.METProducers.MuonTCMETValueMapProducer_cff import *
from RecoMET.METProducers.MetMuonCorrections_cff import *
from RecoMET.Configuration.RecoMET_BeamHaloId_cff import *
from RecoMET.Configuration.RecoTCMET_cff import *


metreco2 = cms.Sequence(
       met+
       metNoHF+
       metHO+
       metNoHFHO+
       calotoweroptmaker+
       metOpt+
       metOptNoHF+
       calotoweroptmakerWithHO+
       metOptHO+
       metOptNoHFHO+
       htMetKT4+
       htMetKT6+
       htMetIC5+
       htMetAK5+
       htMetAK7+
       muonMETValueMapProducer+
       corMetGlobalMuons+
       muonTCMETValueMapProducer+
       tcMet+
       BeamHaloId
       )


if %(MODMETRECO)s:
    process.metreco = metreco2


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
                     process.tcMetAnalyzer*
                     process.analyzeMuonCorrectedCaloMET*
                     process.AnalyzeBeamHalo
)
process.schedule = cms.Schedule(process.p)


