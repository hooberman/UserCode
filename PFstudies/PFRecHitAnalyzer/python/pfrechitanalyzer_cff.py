import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import towerMakerPF
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHCAL_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitPS_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHFEM_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHFHAD_cfi import *
from PFstudies.PFRecHitAnalyzer.pfrechitanalyzer_cfi import *

MYpfRecHitSequence = cms.Sequence(towerMakerPF * particleFlowRecHitECAL * particleFlowRecHitHCAL * particleFlowRecHitPS * particleFlowClusterHFEM * particleFlowClusterHFHAD * pfRecHitAnalyzer)
