import FWCore.ParameterSet.Config as cms

from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
es_prefer_l1GtTriggerMaskTechTrig = cms.ESPrefer("L1GtTriggerMaskTechTrigTrivialProducer","l1GtTriggerMaskTechTrig")

# bit 40+41 no BH  selection

bit4041 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)'))
seqBSC4041Selection = cms.Sequence(bit4041)

#bit 0 selection: BPTX_AND

bptxAnd = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('0'))
seqBPTXAndSelection = cms.Sequence(bptxAnd)


# BPTX OR
 
bptxOr = cms.EDFilter("L1Filter",
                      inputTag = cms.InputTag("gtDigis"),
                      useAODRecord = cms.bool(False),
                      useFinalDecision = cms.bool(False),
                      algorithms = cms.vstring("L1_BptxPlusORMinus")
                      )

seqBPTXOrSelection = cms.Sequence(bptxOr * ~bptxAnd)

bptxPlus = bptxOr.clone(algorithms=cms.vstring("L1_BptxPlus"))

seqBPTXPlusSelection = cms.Sequence(bptxPlus + ~bptxAnd)

bptxMinus = bptxOr.clone(algorithms=cms.vstring("L1_BptxMinus"))

seqBPTXMinusSelection = cms.Sequence(bptxMinus + ~bptxAnd)

# PhysDecl bit

#from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
#physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])
#seqPhysDeclBitSelection = cms.Sequence(physDecl)

from HLTrigger.special.hltPhysicsDeclared_cfi import *
hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'
seqPhysDeclBitSelection = cms.Sequence(hltPhysicsDeclared)



seqPhysCollSelection = cms.Sequence(seqPhysDeclBitSelection + seqBPTXAndSelection)
seqBitSelection = cms.Sequence(seqPhysDeclBitSelection + seqBPTXAndSelection + seqBSC4041Selection)

seqPhysCollSelectionMC = cms.Sequence(seqPhysDeclBitSelection)
seqBitSelectionMC = cms.Sequence(seqPhysDeclBitSelection + seqBSC4041Selection)

seqNonCollidingSelection = cms.Sequence(seqPhysDeclBitSelection + seqBPTXOrSelection)
seqNonCollidingPlusSelection = cms.Sequence(seqPhysDeclBitSelection + seqBPTXPlusSelection)
seqNonCollidingMinusSelection = cms.Sequence(seqPhysDeclBitSelection + seqBPTXMinusSelection)

seqNonCollidingBitSelection = cms.Sequence(seqPhysDeclBitSelection + seqBPTXOrSelection + seqBSC4041Selection)

seqBitSelectionNoPhys = cms.Sequence(seqBPTXAndSelection + seqBSC4041Selection)

seqNonCollidingSelectionNoPhys = cms.Sequence(seqBPTXOrSelection)
seqNonCollidingPlusSelectionNoPhys = cms.Sequence(seqBPTXPlusSelection)
seqNonCollidingMinusSelectionNoPhys = cms.Sequence(seqBPTXMinusSelection)


# there is another prescription to use HLTHighLevel instead of HLTHighLevelDev


