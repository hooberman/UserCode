import FWCore.ParameterSet.Config as cms

demo = cms.EDProducer('PixTcMet',
                      metInputTag   = cms.InputTag("met"),
                      trackInputTag = cms.InputTag("pixelTracks"),
                      beamSpotInputTag  = cms.InputTag("offlineBeamSpot"),
                      checkTrackPropagation = cms.bool(False),
                      d0cuta = cms.double(0.015),
                      d0cutb = cms.double(0.5),
                      maxd0cut = cms.double(0.3),
                      chi2_max = cms.double(50),
                      eta_max = cms.double(2.5), 
                      pt_min  = cms.double(0.),
                      pt_max  = cms.double(4.),
                      radius  = cms.double(130.),
                      zdist  = cms.double(314.),
                      corner = cms.double(1.479)
                    
)
