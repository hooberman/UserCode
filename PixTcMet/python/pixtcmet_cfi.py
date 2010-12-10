import FWCore.ParameterSet.Config as cms

demo = cms.EDProducer('PixTcMet',
                      metInputTag   = cms.InputTag("met"),
                      trackInputTag = cms.InputTag("pixelTracks"),
                      beamSpotInputTag  = cms.InputTag("offlineBeamSpot"),
                      checkTrackPropagation  = cms.bool(False)
)
