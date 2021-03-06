import FWCore.ParameterSet.Config as cms
maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles )
readFiles.extend( [
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/032/DED339AF-FE7B-DE11-85FF-001D09F2525D.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/035/78EEC1B5-ED7B-DE11-B395-001D09F25109.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/039/32DAD72A-E57B-DE11-B3D6-001D09F24FEC.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/046/E6C9E810-BD7C-DE11-8F94-000423D9997E.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/046/2884F426-BD7C-DE11-9A2C-000423D98E6C.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/049/42383E8C-8B7C-DE11-80FB-000423D985E4.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/134/C615C4D9-E77C-DE11-8EC4-000423D94A04.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/142/DCB0463C-FC7C-DE11-86AE-000423D99E46.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/142/D68D84EC-FA7C-DE11-9C5E-000423D98E54.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/143/C4DD14F6-FF7C-DE11-92B6-001D09F2305C.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/144/7A30E64A-037D-DE11-8BDE-000423D94534.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/145/6AAF2547-037D-DE11-9DCF-000423D98C20.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/146/2E78D37F-267D-DE11-8B5D-000423D952C0.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/459/98BE3D3D-A37D-DE11-A086-001D09F24498.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/467/DA2120A7-937D-DE11-A10F-001D09F27067.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/468/EC419160-DC7D-DE11-853B-0030487A1990.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/470/AC92CD13-F57D-DE11-B46F-001D09F29524.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/472/F64A50E3-EB7D-DE11-8270-001D09F2305C.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/474/44E3D24F-167E-DE11-A335-000423D99F3E.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/490/EC644E95-217E-DE11-9DB7-000423D990CC.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/504/BA5FDDAC-347E-DE11-A653-001D09F24DDA.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/508/1AD008CD-727E-DE11-868C-000423D99614.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/508/B88E0B46-727E-DE11-B5D8-001D09F248FD.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/508/20DEEF4F-727E-DE11-898A-000423D94C68.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/519/B27C8339-B07E-DE11-BC2D-000423D98750.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/519/9CDDD942-B07E-DE11-89C6-001D09F2441B.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/524/CC3C30D1-7180-DE11-BC0F-000423D98C20.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/524/9663A4DA-7180-DE11-997D-000423D6B444.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/524/6E93B5D5-7180-DE11-BD52-000423D98634.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/524/FCE768D6-7180-DE11-84C4-000423D94E70.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/562/D003D58D-3B80-DE11-A38E-001D09F2532F.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/562/2E306B90-3B80-DE11-B434-001D09F24493.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/573/5C90805E-A77F-DE11-B67D-000423D999CA.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/575/C04D06E6-AA7F-DE11-B9BC-000423D99614.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/578/566F4E8A-7E80-DE11-B8BB-001D09F2447F.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/584/2CA9FBCF-C77F-DE11-80E3-000423D99394.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/584/E6D30BD6-C77F-DE11-A407-000423D6AF24.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/586/063A140A-DD7F-DE11-9682-0019B9F705A3.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/603/FC0B8E17-3880-DE11-A34C-001D09F2AF1E.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/606/A8D74481-4080-DE11-97BD-000423D996C8.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/616/80083D48-4180-DE11-A95D-000423D9880C.root',
'/store/data/CRAFT09/Cosmics/ALCARECO/v1/000/109/624/B8FB48A9-8080-DE11-B147-001D09F23D04.root'
] ); 
