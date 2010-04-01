
import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineValidator")

process.load("Alignment.PeakDecoResiduals.DataSetMinBias_Mar30_2010_dec_cff")

#process.source.inputCommands = cms.untracked.vstring('keep *', 'drop *_MEtoEDMConverter_*_*') # hack to get rid of the memory consumption problem in 2_2_X and beond
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(False),
   Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
   fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

 ##
 ## Maximum number of Events
 ## 
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#lumisToProcess = cms.untracked.VLuminosityBlockRange(
#    #'124275:3-124275:30',
#    '124230:26-124230:max',
#    '124120:1-124120:max',
#    '124030:1-124030:max',
#    '124027:24-124027:max',
#    '124025:3-124025:13',
#    '124024:2-124024:83',
#    '124023:38-124023:96',
#    '124022:66-124022:179',
#    '124020:12-124020:94',
#    '124009:1-124009:68',
#    '124008:1-124008:max',
#    '124006:1-124006:max',
#    '123908:2-123908:12',
#    '123906:18-123906:28',
#    '123818:2-123818:42',
#    '123815:8-123815:max',
#    '123732:62-123732:109',
#    '123615:70-123615:max',
#    '123596:2-123596:max',
#    '123592:2-123592:12'
#    )

 ##   
 ## Messages & Convenience
 ##
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
process.TrackerTrackHitFilter.useTrajectories= True  # this is needed only if you require some selections; but it will work even if you don't ask for them
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

 ##
 ## Load and Configure TrackRefitter1
 ##

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

#############
# parameters for TrackRefitter
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


 ##
 ## Get the BeamSpot
 ##
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
 
 ##
 ## GlobalTag Conditions (if needed)
 ##
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "GR09_R_34X_V2::All"
#process.GlobalTag.globaltag = "GR09_R_V6A::All"
process.GlobalTag.globaltag = "GR10_P_V4::All"

#use lorentz angle from global tag
  
 ##
 ## Geometry
 ##
process.load("Configuration.StandardSequences.Geometry_cff")
 
 ##
 ## Magnetic Field
 ##
process.load("Configuration/StandardSequences/MagneticField_38T_cff")


#from CondCore.DBCommon.CondDBSetup_cfi import *
from CalibTracker.Configuration.Common.PoolDBESSource_cfi import poolDBESSource
##include private db object
##
import CalibTracker.Configuration.Common.PoolDBESSource_cfi
process.trackerAlignment =  CalibTracker.Configuration.Common.PoolDBESSource_cfi.poolDBESSource.clone(
    #connect = cms.string('sqlite_file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/eaguiloc/Merged_1st_900GeV/alignments_iter15.db'),
    connect = cms.string('sqlite_file:TrackerAlignment_Feb2010Cosmics_38T.db'),
    timetype = cms.string("runnumber"),
    toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                               tag = cms.string('Alignments')
                               ))
    )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")


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
process.load("Alignment.PeakDecoResiduals.PeakDecoResiduals_cfi")
process.PeakDecoResiduals.Tracks = 'TrackRefitter2'
process.PeakDecoResiduals.trajectoryInput = 'TrackRefitter2'
process.PeakDecoResiduals.debug = cms.bool(False)
process.PeakDecoResiduals.runOnCosmics = cms.bool(False)
process.PeakDecoResiduals.createTree = cms.bool(True)


#process.TFileService.fileName = '/tmp/benhoob/temp.root'
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("temp.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.p = cms.Path(process.offlineBeamSpot*
                     process.TrackRefitter1*
                     process.TrackerTrackHitFilter*
                     process.HitFilteredTracks*
                     process.AlignmentTrackSelector*
                     process.TrackRefitter2*
                     process.PeakDecoResiduals)

