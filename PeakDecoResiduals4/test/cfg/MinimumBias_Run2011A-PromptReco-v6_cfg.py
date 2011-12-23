import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineValidator")

## Maximum number of Events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/630/883CFB93-4FC0-E011-BA94-003048D375AA.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/630/6633024F-54C0-E011-B3AB-003048F118D4.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/F40D9540-0DC0-E011-84E7-0030486780AC.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/EC29840A-09C0-E011-9CAD-485B39897227.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/D4C31A0B-09C0-E011-BA10-485B3962633D.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/D4C19A0A-09C0-E011-B038-BCAEC518FF80.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/6019850A-09C0-E011-8487-BCAEC5329720.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/4E489C40-0DC0-E011-83F0-003048D2C0F2.root',
#        '/store/data/Run2011A/MinimumBias/RECO/PromptReco-v6/000/172/620/2E678643-0DC0-E011-959F-BCAEC53296F9.root'
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
    src = 'generalTracks',
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
process.TrackRefitter1 = process.TrackRefitter.clone(
   src = 'HighPuritySelector',
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
    process.noscraping*
    process.primaryVertexFilter*
    process.offlineBeamSpot*
    process.HighPuritySelector*
    process.TrackRefitter1*
    process.TrackerTrackHitFilter*
    process.HitFilteredTracks*
    process.AlignmentTrackSelector*
    process.TrackRefitter2*
    process.PeakDecoResiduals4)

