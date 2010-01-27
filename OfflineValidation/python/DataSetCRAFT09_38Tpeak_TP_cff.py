#/Cosmics/CRAFT09-TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/RAW-RECO 109011-109624
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles )
readFiles.extend( [
    ###109046
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0012/E4E9BE9C-49BB-DE11-94BD-001A92971B8E.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0008/003531A1-DBB9-DE11-AC25-0030486792B6.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/F2AE14A7-AAB8-DE11-881B-001A92971B3A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/EA75FFF0-ACB8-DE11-9CC6-003048678E2A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/E865D886-AFB8-DE11-AE15-001A92810AC6.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/D2598E3F-D6B8-DE11-A97A-001A92810AAE.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/C23AB601-BCB8-DE11-85E6-001A92810AA0.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/BE748966-C6B8-DE11-9DE2-001A92971B90.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/AA7AC47C-ABB8-DE11-AF93-001731AF68C3.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/A0EDABD4-B7B8-DE11-AEB3-001A92971B1A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/9A6332F2-ADB8-DE11-B04D-0017312B577F.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/4445656E-C3B8-DE11-9CC8-003048679168.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/3E8752C6-C9B8-DE11-9FDC-001A92971B48.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/300B8366-B0B8-DE11-A82D-0018F3D096C8.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/2EB6B16F-B6B8-DE11-95EA-00304867906C.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0003/185D7BD0-B4B8-DE11-B480-003048678B7C.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/EEDA19F9-99B8-DE11-A0F3-003048678B30.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/D410109C-9EB8-DE11-8248-00304867906C.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/CAB3E696-A2B8-DE11-93BE-001A928116F4.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/BC64991D-99B8-DE11-8EC2-001A928116F4.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/A4552C87-9BB8-DE11-AEB1-001A92811724.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/A0A0BD1C-A1B8-DE11-A609-003048678A78.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/8C98D038-9DB8-DE11-9A74-003048678DD6.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/78666380-9FB8-DE11-AFEA-0018F3D0961A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/54C015A9-98B8-DE11-8960-001731A281B1.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/4C803A16-A4B8-DE11-943B-003048678A6A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/44CC7A0A-A9B8-DE11-9D0A-00304867C136.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/402D7B15-A5B8-DE11-84CB-001731AF677B.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/2A057CBA-97B8-DE11-9440-001731AF6651.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0002/10BC2186-A7B8-DE11-9AE5-003048678B30.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/EE772AB5-94B8-DE11-8FF0-0018F3D0961A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/E678FCC8-8DB8-DE11-B5EB-00304867BFF2.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/DE415537-83B8-DE11-B27C-0018F3D0962A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/B87A3ECE-8BB8-DE11-8096-003048678B86.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/B47AEFA2-91B8-DE11-8F7B-001731AF66AB.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/AEDECD5E-8AB8-DE11-B72F-001A92971B7C.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/9A145BB4-85B8-DE11-93EE-001A928116F2.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/96C02A03-88B8-DE11-B892-003048678E92.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/7CEFC746-92B8-DE11-9375-001A92811736.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/705A5146-96B8-DE11-9B2A-001A92810AD0.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/324F0E9D-8FB8-DE11-8480-003048678FE4.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/1ACBD7AC-84B8-DE11-8A3E-0018F3D0965C.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/1A1DAA80-96B8-DE11-A63E-003048679030.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/08319FC6-90B8-DE11-B56E-001A92971B54.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/02315244-93B8-DE11-99CA-001A9281173E.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0001/0212CB0E-88B8-DE11-A0F8-0030486791F2.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/FCA980F2-6BB8-DE11-8308-001731AF692F.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/F4E26D70-7FB8-DE11-B336-001A92810A94.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/F232FD58-6EB8-DE11-A656-001731AF6845.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/D20B877E-74B8-DE11-B202-0018F3D09634.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/BAA16338-6FB8-DE11-B2B5-001A9281174C.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/BA06FCEC-80B8-DE11-8556-00304867915A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/AEEBFD07-7CB8-DE11-89F0-001A9281174A.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/A4EECDD6-7FB8-DE11-B711-001A92810A94.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/668E4528-78B8-DE11-BDAD-001731AF66B9.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/664EBA5A-75B8-DE11-8D6E-003048678BC6.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/628ED5B0-77B8-DE11-9A38-0018F3D09660.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/5C4734D2-76B8-DE11-9B18-003048678BEA.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/4EF9D16F-7EB8-DE11-A2A4-003048678F62.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/46893017-7DB8-DE11-BB7B-0018F3D095FC.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/40FF6A2B-73B8-DE11-AF35-003048679150.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/36BF95AC-7AB8-DE11-8E8A-003048679166.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/1EA182C9-6FB8-DE11-B6C5-0017312B5DC9.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/1E6F9E9D-70B8-DE11-8001-001731AF66B3.root',
    '/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0000/125FF631-73B8-DE11-9B23-0018F3D09612.root'

    ###109032
    #'/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0008/70748AFD-DAB9-DE11-A8AF-003048678AC0.root',
    #'/store/data/CRAFT09/Cosmics/RAW-RECO/TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/0004/7E9F6429-F7B8-DE11-9EB4-003048678BAC.root',
] ); 
