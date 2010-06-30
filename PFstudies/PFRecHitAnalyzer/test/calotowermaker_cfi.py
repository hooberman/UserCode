import FWCore.ParameterSet.Config as cms
 
calotowermaker = cms.EDProducer("CaloTowersCreator",
     # Energy threshold for EB 5x5 crystal inclusion [GeV]
     EBSumThreshold = cms.double(0.2),
     # Weighting factor for HF short-fiber readouts
     HF2Weight = cms.double(1.0),
     # Weighting factor for EB   
     EBWeight = cms.double(1.0),
     # Label of HFRecHitCollection to use
     hfInput = cms.InputTag("hfreco"),
     # Energy threshold for EE crystals-in-tower inclusion [GeV]
     EESumThreshold = cms.double(0.45),
     # Energy threshold for HO cell inclusion [GeV]
     HOThreshold0 = cms.double(1.1),
     HOThresholdPlus1 = cms.double(3.5),
     HOThresholdMinus1 = cms.double(3.5),
     HOThresholdPlus2 = cms.double(3.5),
     HOThresholdMinus2 = cms.double(3.5),
     HBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     # Energy threshold for HB cell inclusion [GeV]
     HBThreshold = cms.double(0.7),                            #<<<---HB threshold
     EEWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     # Energy threshold for long-fiber HF readout inclusion [GeV]
     HF1Threshold = cms.double(0.5),                           #<<<---HF long threshold
     HF2Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     HOWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     EEGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     # Weighting factor for HE 10-degree cells   
     HEDWeight = cms.double(1.0),
     # Weighting factor for EE   
     EEWeight = cms.double(1.0),
     # HO on/off flag for tower energy reconstruction
     UseHO = cms.bool(True),
     HBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     # Weighting factor for HE 5-degree cells   
     HESWeight = cms.double(1.0),
     # Weighting factor for HF long-fiber readouts 
     HF1Weight = cms.double(1.0),
     HF2Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     HEDWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     HF1Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     EBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     # Weighting factor for HO 
     HOWeight = cms.double(1.0),
                               
     # Energy threshold for EB crystal inclusion [GeV]
     EBThreshold = cms.double(0.07),                         #<<<---EB threshold
     # Energy threshold for EE crystal inclusion [GeV]
     EEThreshold = cms.double(0.3),                          #<<<---EE threshold
     # Flags specifying if the above thresholds
     # should be applied to Et (UseEtEXTreshold='True') or E ('False')
     # Flags for use of symmetric thresholds: |e|>threshold                          
     UseEtEBTreshold = cms.bool(False),
     UseSymEBTreshold = cms.bool(True),
     UseEtEETreshold = cms.bool(False),
     UseSymEETreshold = cms.bool(True),
 
 
     # Label of HBHERecHitCollection to use
     hbheInput = cms.InputTag("hbhereco"),
     # Global energy threshold on Hcal [GeV]
     HcalThreshold = cms.double(-1000.0),
     # Energy threshold for short-fiber HF readout inclusion [GeV]
     HF2Threshold = cms.double(0.85),                                  #<<<---HF short threshold
 
     # Energy threshold for 5-degree (phi) HE cell inclusion [GeV]
     HESThreshold = cms.double(0.8),                                   #<<<---HE threshold
     HF1Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     # Label of HORecHitCollection to use
     hoInput = cms.InputTag("horeco"),
     HESGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     #
     HESWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
     # Energy threshold for 10-degree (phi) HE cel inclusion [GeV]
     HEDThreshold = cms.double(0.8),                                   #<<<---HE threshold
     # Global energy threshold on tower [GeV]
     EcutTower = cms.double(-1000.0),
     HEDGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     # Label of EcalRecHitCollections to use
     ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEE")),
     # Weighting factor for HB   
     HBWeight = cms.double(1.0),
     HOGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     # Energy dependent weights and energy scale to be used
     EBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
     #
     #
     # momentum assignment
     # Method for momentum reconstruction
     MomConstrMethod = cms.int32(1),                           
     # Depth, fraction of the respective calorimeter [0,1]
     MomHBDepth = cms.double(0.2),
     MomHEDepth = cms.double(0.4),   
     MomEBDepth = cms.double(0.3),
     MomEEDepth = cms.double(0.0),
 
 # parameters for handling of anomalous cells
 # 
     # acceptable severity level
     HcalAcceptSeverityLevel = cms.uint32(9),
     EcalAcceptSeverityLevel = cms.uint32(1),
        
     # use of recovered hits
     UseHcalRecoveredHits = cms.bool(True),
       # The CaloTower code treats recovered cells as a separate category.
       # The flag to use (or not use) them should be explicitly set
       # regardless of the specified severity level in EcalAcceptSeverityLevel.                       
     UseEcalRecoveredHits = cms.bool(False),
                         
                               
 
 # flag to allow/disallow missing inputs
     AllowMissingInputs = cms.bool(False)
     
