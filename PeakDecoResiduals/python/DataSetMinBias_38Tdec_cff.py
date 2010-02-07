import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    
    #/MinimumBias/BeamCommissioning09-BSCNOBEAMHALO-Dec19thSkim_341_v2/RAW-RECO RUN 124275
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/9C264FD4-2CEE-DE11-9CB1-002618943945.root',
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/56B6C4FB-F5ED-DE11-8763-002618943930.root',
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/001FFE50-B9ED-DE11-9A24-002618943882.root',
    '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0005/04CA81AC-A3ED-DE11-9414-0018F3D09630.root'
    
    ] );


secFiles.extend( [
    ] )


















###/MinimumBias/BeamCommissioning09-BSCNOBEAMHALO-Jan23Skim-v1/RAW-RECO RUN 123596
#'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Jan23Skim-v1/0015/06AA87A8-F809-DF11-9F02-0026189438B3.root',

###/MinimumBias/BeamCommissioning09-BSCNOBEAMHALO-Jan29Skim-v2/RAW-RECO RUN 122294
#'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Jan29Skim-v2/0022/C455E3C4-760E-DF11-8F30-00261894393E.root'

#/MinimumBias/BeamCommissioning09-BSCNOBEAMHALO-Dec19thSkim_341_v2/RAW-RECO 
#'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/F682F559-B9ED-DE11-8057-002618943900.root',
#'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/F6361662-B9ED-DE11-877E-002618943950.root',
