import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineValidator")


process.load("Alignment.PeakDecoResiduals.Commissioning10_Apr20Skim_GOODCOLL_v1_cff")

## Maximum number of Events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


#process.source.inputCommands = cms.untracked.vstring('keep *', 'drop *_MEtoEDMConverter_*_*') # hack to get rid of the memory consumption problem in 2_2_X and beond
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(False),
   Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
   fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)


    
## Messages & Convenience
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
reportEvery = cms.untracked.int32(1000) # every 1000th only
#    limit = cms.untracked.int32(10)       # or limit to 10 printouts...
))
process.MessageLogger.statistics.append('cout')
      
#-- Track hit filter
# TrackerTrackHitFilter takes as input the tracks/trajectories coming out from TrackRefitter1
process.load("RecoTracker.FinalTrackSelectors.TrackerTrackHitFilter_cff")
process.TrackerTrackHitFilter.src = 'TrackRefitter1'

#-- Alignment Track Selection
process.load("Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi")
process.AlignmentTrackSelector.src = 'HitFilteredTracks'
process.AlignmentTrackSelector.filter = True



##### For Tracks:
process.AlignmentTrackSelector.applyBasicCuts = True
process.AlignmentTrackSelector.pMin    = 3
process.AlignmentTrackSelector.pMax    = 9999.
process.AlignmentTrackSelector.ptMin   = 0.
process.AlignmentTrackSelector.ptMax   = 9999.
process.AlignmentTrackSelector.etaMin  = -999.
process.AlignmentTrackSelector.etaMax  = 999.
process.AlignmentTrackSelector.nHitMin = 8
process.AlignmentTrackSelector.nHitMin2D = 2
process.AlignmentTrackSelector.chi2nMax = 999.
process.AlignmentTrackSelector.applyMultiplicityFilter = False
process.AlignmentTrackSelector.maxMultiplicity = 1
process.AlignmentTrackSelector.applyNHighestPt = False
process.AlignmentTrackSelector.nHighestPt = 1
process.AlignmentTrackSelector.seedOnlyFrom = 0 
process.AlignmentTrackSelector.applyIsolationCut = False
process.AlignmentTrackSelector.minHitIsolation = 0.8
process.AlignmentTrackSelector.applyChargeCheck = False
process.AlignmentTrackSelector.minHitChargeStrip = 50.

##### For Hits:
process.TrackerTrackHitFilter.useTrajectories= True # this is needed only if you require some selections; but it will work even if you don't ask for them  
process.TrackerTrackHitFilter.minimumHits = 6
process.TrackerTrackHitFilter.commands = cms.vstring("keep PXB","keep PXE","keep TIB","keep TID","keep TOB","keep TEC")
process.TrackerTrackHitFilter.detsToIgnore = [
    # see https://hypernews.cern.ch/HyperNews/CMS/get/tracker-performance/484.html
    # TIB / TID
    #369136710, 369136714, 402668822,
    # TOB
    #436310989, 436310990, 436299301, 436299302,
    # TEC
    #470340521, 470063045, 470063046, 470114669, 470114670, 470161093, 470161094, 470164333, 470164334, 470312005, 470312006, 470312009, 470067405, 470067406, 470128813
]
process.TrackerTrackHitFilter.replaceWithInactiveHits = True
process.TrackerTrackHitFilter.stripAllInvalidHits = False
process.TrackerTrackHitFilter.rejectBadStoNHits = True
process.TrackerTrackHitFilter.StoNcommands = cms.vstring("ALL 14.0")
#process.TrackerTrackHitFilter.StoNcommands = cms.vstring("ALL 8.0")
process.TrackerTrackHitFilter.rejectLowAngleHits= True
process.TrackerTrackHitFilter.TrackAngleCut= 0.17 # in rads, starting from the module surface
process.TrackerTrackHitFilter.usePixelQualityFlag= True


#now we give the TrackCandidate coming out of the TrackerTrackHitFilter to the track producer
import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff
process.HitFilteredTracks = RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff.ctfWithMaterialTracks.clone(
   src = 'TrackerTrackHitFilter',
   NavigationSchool = "",

###    TrajectoryInEvent = True,
    TTRHBuilder = "WithAngleAndTemplate"    
)


## Load and Configure TrackRefitter1
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")


## parameters for TrackRefitter
#process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff")
import RecoTracker.TrackProducer.TrackRefitters_cff
process.TrackRefitter1 = process.TrackRefitter.clone(
  #src = 'ALCARECOTkAlMinBias',
   src = 'generalTracks',
   TrajectoryInEvent = True,
   TTRHBuilder = "WithAngleAndTemplate",
   NavigationSchool = ""
)
process.TrackRefitter2 = process.TrackRefitter1.clone(
#    src = 'HitFilteredTracks')
     src = 'AlignmentTrackSelector'
)

## Get the BeamSpot
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")

## GlobalTag Conditions (if needed)
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_35X_V7A::All"

#use lorentz angle from global tag
  
## Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
 
## Magnetic Field
process.load("Configuration/StandardSequences/MagneticField_38T_cff")


#from CondCore.DBCommon.CondDBSetup_cfi import *
from CalibTracker.Configuration.Common.PoolDBESSource_cfi import poolDBESSource

##include private db object

#include Alignment sqlite file
import CalibTracker.Configuration.Common.PoolDBESSource_cfi
process.trackerAlignment =  CalibTracker.Configuration.Common.PoolDBESSource_cfi.poolDBESSource.clone(
    connect = cms.string('sqlite_file:TrackerAlignment_Feb2010Cosmics_38T.db'),
    timetype = cms.string("runnumber"),
    toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                               tag = cms.string('Alignments')
                               ))
    )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")

#include TrackerAlignmentErrors
from CondCore.DBCommon.CondDBSetup_cfi import *
process.APE = poolDBESSource.clone(
    connect = cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X'),
    timetype = cms.string("runnumber"),
    toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                               tag = cms.string('TrackerIdealGeometryErrors210_mc')
                               ))
    )
process.es_prefer_APE = cms.ESPrefer("PoolDBESSource", "APE")


#include custom sqlite file
process.stripLorentzAngle = cms.ESSource("PoolDBESSource",CondDBSetup,
                                         connect = cms.string('sqlite_file:SiStripLorentzAngle_Deco.db'),  
                                         toGet = cms.VPSet(cms.PSet(record = cms.string('SiStripLorentzAngleRcd'),
                                                                    tag = cms.string('SiStripLorentzAngle_Deco') ))
                                         )
process.es_prefer_stripLorentzAngle = cms.ESPrefer("PoolDBESSource", "stripLorentzAngle")


###--------------------------------------------------------------------------------
###filters from https://twiki.cern.ch/twiki/bin/view/CMS/TrackerAlignmentWithBeams2010

## process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
## process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
## process.L1T1=process.hltLevel1GTSeed.clone()
## process.L1T1.L1TechTriggerSeeding = cms.bool(True)
## process.L1T1.L1SeedsLogicalExpression=cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND (NOT 43)) OR (43 AND (NOT 42)))')

# Filter for HLT PhysicsDeclared

## process.hltHighLevel = cms.EDFilter("HLTHighLevel",
## TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
## HLTPaths = cms.vstring('HLT_PhysicsDeclared'),           # provide list of HLT paths (or patterns) you want
## eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
## andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
## throw = cms.bool(True)    # throw exception on unknown path names 
##  )

# Filter for high purity tracks

import Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi
process.HighPuritySelector = Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi.AlignmentTrackSelector.clone(
applyBasicCuts = True,
filter = True,
src = '.oO[TrackCollection]Oo.',
trackQualities = ["highPurity"]
)

# Back-plane correction for the TOB in deconvolution mode

process.load("RecoLocalTracker.SiStripRecHitConverter.OutOfTime_cff")
process.OutOfTime.TOBlateBP=0.075
process.OutOfTime.TIBlateBP=0.036

# The corrected beamspot is included:

############################################################################
# load corrected beamspot
###########################################################################
## process.myBeamspot = poolDBESSource.clone(
##     connect = cms.string('frontier://PromptProd/CMS_COND_31X_BEAMSPOT'),
##     toGet = cms.VPSet(
##       cms.PSet(
##        record = cms.string('BeamSpotObjectsRcd'),
##         tag = cms.string('BeamSpotObjects_2009_v7_offline')
##        )
##       )
##    )
## process.es_prefer_corrBeamspot = cms.ESPrefer("PoolDBESSource","myBeamspot")

###--------------------------------------------------------------------------------

process.load("Alignment.PeakDecoResiduals.bitSelectionSequence_cff")


## to apply misalignments
#TrackerDigiGeometryESModule.applyAlignment = True
   
### Load and Configure PeakDecoResiduals
process.load("Alignment.PeakDecoResiduals.PeakDecoResiduals_cfi")
process.PeakDecoResiduals.Tracks = 'TrackRefitter2'
process.PeakDecoResiduals.trajectoryInput = 'TrackRefitter2'
process.PeakDecoResiduals.debug = cms.bool(False)
process.PeakDecoResiduals.runOnCosmics = cms.bool(False)
process.PeakDecoResiduals.createTree = cms.bool(False)


#process.TFileService.fileName = '/tmp/benhoob/temp.root'
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("temp.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.p = cms.Path(
    process.seqBitSelection*
    process.offlineBeamSpot*
    process.TrackRefitter1*
    process.TrackerTrackHitFilter*
    process.HitFilteredTracks*
    process.AlignmentTrackSelector*
    process.TrackRefitter2*
    process.PeakDecoResiduals
    )

