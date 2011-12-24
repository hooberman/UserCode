import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineValidator")

## Maximum number of Events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/659/7A4F15D6-2AC0-E011-A7D8-003048F024FE.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/659/166BBCD6-2AC0-E011-BE8B-0030487CD6B4.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/635/E2BCA3E5-06C1-E011-A6E0-003048F1C420.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/630/D8D14B71-E1C0-E011-BFEF-003048F11DE2.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/627/B01AA350-08C0-E011-8A21-BCAEC5329702.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/625/2CD8D5D0-07C0-E011-AE1D-BCAEC5329730.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/622/B66B697A-FDBF-E011-BE79-BCAEC518FF69.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/622/A2155A7B-FDBF-E011-9931-003048673374.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/622/6AD3987A-FDBF-E011-8A83-E0CB4E55365C.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/622/44F3C977-FDBF-E011-8016-BCAEC518FF8F.root',
#'/store/data/Run2011A/Cosmics/RECO/PromptReco-v6/000/172/620/D451FF13-42C0-E011-86C4-BCAEC5329732.root'

    ),
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


#--------------------------------
# filter beam scraping events
#--------------------------------

process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(True),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

#--------------------------------
# require a good vertex
#--------------------------------

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(15), 
                                           maxd0 = cms.double(2) 
                                           )


    
#############################################################
##select only high purity tracks
##has to run first as necessary information
##is only available in initial track selection
##(Quality information is thrown away by the tracker refitters)
##########################################################

import Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi
process.HighPuritySelector = Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi.AlignmentTrackSelector.clone(
    applyBasicCuts = True,
    filter = True,
    #src = 'ALCARECOTkAlMinBias',
    src = 'ctfWithMaterialTracksP5',
    #src = 'generalTracks',
    trackQualities = ["highPurity"]
    )
      
#-- Track hit filter
# TrackerTrackHitFilter takes as input the tracks/trajectories coming out from TrackRefitter1
process.load("RecoTracker.FinalTrackSelectors.TrackerTrackHitFilter_cff")
process.TrackerTrackHitFilter.src = 'TrackRefitter1'

#-- Alignment Track Selection
process.load("Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi")
process.AlignmentTrackSelector.src = 'HitFilteredTracks'
process.AlignmentTrackSelector.filter = True


##### For Tracks:collisions taken in deco mode

####################################################
##include back-plane correction for deco mode
######################################################
#process.load("RecoLocalTracker.SiStripRecHitConverter.OutOfTime_cff.py")
#process.load("RecoLocalTracker.SiStripRecHitConverter.OutOfTime_cff")
#process.OutOfTime.TOBlateBP=0.071
#process.OutOfTime.TIBlateBP=0.036

process.AlignmentTrackSelector.applyBasicCuts = True
# Note that pMin is overridden and set to zero in
# the offlineTemplate0T
process.AlignmentTrackSelector.pMin    = 3
process.AlignmentTrackSelector.pMax    = 9999.
process.AlignmentTrackSelector.ptMin   = 0.65
process.AlignmentTrackSelector.ptMax   = 9999.
process.AlignmentTrackSelector.etaMin  = -999.
process.AlignmentTrackSelector.etaMax  = 999.
process.AlignmentTrackSelector.nHitMin = 7
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
#process.AlignmentTrackSelector.trackQualities = ["highPurity"]
#process.AlignmentTrackSelector.iterativeTrackingSteps = ["iter1","iter2"]




##### For Hits:
process.TrackerTrackHitFilter.useTrajectories= True  # this is needed only if you require some selections; but it will work even if you don't ask for them
process.TrackerTrackHitFilter.minimumHits = 7
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
#process.TrackerTrackHitFilter.StoNcommands = cms.vstring("ALL 14.0")
process.TrackerTrackHitFilter.StoNcommands = cms.vstring("ALL 12.0")
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
#process.TrackRefitter1 = process.TrackRefitter.clone(
process.TrackRefitter1 = process.TrackRefitterP5.clone(
#   src = 'HighPuritySelector',
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
process.GlobalTag.globaltag = "GR_P_V22::All"

#use lorentz angle from global tag
  
## Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
 
## Magnetic Field
process.load("Configuration/StandardSequences/MagneticField_38T_cff")


#from CondCore.DBCommon.CondDBSetup_cfi import *
from CalibTracker.Configuration.Common.PoolDBESSource_cfi import poolDBESSource

##include private db object
#import CalibTracker.Configuration.Common.PoolDBESSource_cfi
#process.trackerAlignment =  CalibTracker.Configuration.Common.PoolDBESSource_cfi.poolDBESSource.clone(
#    connect = cms.string('sqlite_file:TrackerAlignment_Feb2010Cosmics_38T.db'),
#    timetype = cms.string("runnumber"),
#    toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
#                               tag = cms.string('Alignments')
#                               ))
#    )
#process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")


from CondCore.DBCommon.CondDBSetup_cfi import *
process.APE = poolDBESSource.clone(
    connect = cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X'),
    timetype = cms.string("runnumber"),
    toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                               tag = cms.string('TrackerIdealGeometryErrors210_mc')
                               ))
    )
process.es_prefer_APE = cms.ESPrefer("PoolDBESSource", "APE")


## to apply misalignments
#TrackerDigiGeometryESModule.applyAlignment = True
   
### Load and Configure PeakDecoResiduals
process.load("Alignment.PeakDecoResiduals4.PeakDecoResiduals4_cfi")
process.PeakDecoResiduals4.Tracks = 'TrackRefitter2'
process.PeakDecoResiduals4.trajectoryInput = 'TrackRefitter2'
process.PeakDecoResiduals4.debug = cms.bool(False)
process.PeakDecoResiduals4.runOnCosmics = cms.bool(False)
process.PeakDecoResiduals4.createTree = cms.bool(True)


#process.TFileService.fileName = '/tmp/benhoob/temp.root'
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("temp.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


process.p = cms.Path(
#    process.noscraping*
#    process.primaryVertexFilter*
    process.offlineBeamSpot*
#    process.HighPuritySelector*
    process.TrackRefitter1*
    process.TrackerTrackHitFilter*
    process.HitFilteredTracks*
    process.AlignmentTrackSelector*
    process.TrackRefitter2*
    process.PeakDecoResiduals4
)

