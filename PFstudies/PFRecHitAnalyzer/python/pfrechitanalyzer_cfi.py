import FWCore.ParameterSet.Config as cms

pfRecHitAnalyzer = cms.EDProducer('PFRecHitAnalyzer',
                                  pf_ecal_tag  = cms.InputTag("particleFlowRecHitECAL"          ), 
                                  pf_hcal_tag  = cms.InputTag("particleFlowRecHitHCAL"          ), 
                                  pf_hfem_tag  = cms.InputTag("particleFlowRecHitHCAL", "HFEM"  ), 
                                  pf_hfhad_tag = cms.InputTag("particleFlowRecHitHCAL", "HFHAD" ),
                                  calomet_tag  = cms.InputTag("met"),
                                  genmet_tag   = cms.InputTag("genMetTrue"),
                                  eb_threshold = cms.double(0),
                                  ee_threshold = cms.double(0),
                                  he_threshold = cms.double(0),
                                  hb_threshold = cms.double(0),
                                  hfe_threshold = cms.double(0),
                                  hfh_threshold = cms.double(0)
                                  )
