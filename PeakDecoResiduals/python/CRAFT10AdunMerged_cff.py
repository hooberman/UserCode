import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
    
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/763/2A7705D2-CC1F-DF11-BD19-001617C3B778.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/761/EA720BFD-991F-DF11-A5A1-0030487CD180.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/760/16BF0AAD-951F-DF11-BD30-0019DB29C5FC.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/759/FE096FB2-AD1F-DF11-86B8-00304879EDEA.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/759/1CD15997-B21F-DF11-8773-000423D95030.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/758/B2146745-A51F-DF11-B9F3-001D09F24D8A.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/758/28040BF3-AC1F-DF11-9EE3-000423D6AF24.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/757/50333405-701F-DF11-AE29-0019B9F705A3.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/755/869E3186-C81F-DF11-BCD3-0030487CD7EA.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/755/1612DD6C-C41F-DF11-96C6-000423D951D4.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/753/4CBD1A80-561F-DF11-A14C-001D09F242EA.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/752/2C8158BA-5D1F-DF11-AF61-001617E30E28.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/750/48294D45-6F1F-DF11-828E-0030487A18A4.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/749/3AF1608A-391F-DF11-AF17-0030487C5CE2.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/748/462724C6-161F-DF11-82B4-001617E30F48.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/747/F2F79810-EF1E-DF11-A0D6-001D09F232B9.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/746/7649BAF7-291F-DF11-8939-001D09F28EA3.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/746/68A43AB1-2F1F-DF11-AAA3-001617E30CE8.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/745/F4C99359-591F-DF11-8421-0030486733B4.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/744/7C20D115-DC1E-DF11-A580-000423D60FF6.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/743/EE17FF3C-DC1E-DF11-8FC4-000423D98BC4.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/741/F4135F6E-C61E-DF11-9DBC-003048D2BBF0.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/740/448AF101-D51E-DF11-85A6-000423D6006E.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/739/B290573F-C71E-DF11-B308-001D09F252F3.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/736/149EBED0-A01E-DF11-99C1-000423D9890C.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/736/040F4B6F-971E-DF11-B024-000423D6CA42.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/734/68358CF6-7B1E-DF11-87F7-0030487CD7C0.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/734/48D28E03-771E-DF11-82F1-001617C3B77C.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/732/126B5470-781E-DF11-B432-000423D98B08.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/730/E625DA4D-5E1E-DF11-9AF3-001D09F2438A.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/730/9A1F2730-631E-DF11-AF63-001D09F24FEC.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/729/CEDDAC2E-381E-DF11-B959-0019B9F72BAA.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/729/4EDDBA49-3A1E-DF11-AB2F-0030487CD700.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/728/52966407-2F1E-DF11-8F5A-0030487C6062.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/725/80C81D70-161E-DF11-AEF1-001D09F34488.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/719/D4BCA149-1B1E-DF11-849B-000423D999CA.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/719/704B11D6-171E-DF11-AFBF-003048673374.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/718/D05F3048-0D1E-DF11-BDF2-0030487A3232.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/715/B28CE4AA-281E-DF11-97B0-000423D944DC.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/713/DC0FE7AE-151E-DF11-9B2C-003048D3750A.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/713/C841AAF6-191E-DF11-B3E5-000423D99F3E.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/712/22359F9C-0E1E-DF11-8DC2-000423D99658.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/712/1291F03A-141E-DF11-9D1D-000423D98E30.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/711/3A2BBC70-051E-DF11-94B3-0030487C778E.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/711/046B7B84-001E-DF11-8AE2-000423D94494.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/706/8C466A1B-EF1D-DF11-858B-001D09F2424A.root',
    '/store/data/Commissioning10/Cosmics/ALCARECO/v3/000/128/706/387F6DD4-F31D-DF11-B9AE-0019B9F704D6.root'  
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
