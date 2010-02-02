import FWCore.ParameterSet.Config as cms
maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles )
readFiles.extend( [

#/MinimumBias/BeamCommissioning09-BSCNOBEAMHALO-Jan23Skim-v1/RAW-RECO RUN 123596 
'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Jan23Skim-v1/0015/06AA87A8-F809-DF11-9F02-0026189438B3.root',

] ); 
