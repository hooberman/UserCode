
import FWCore.ParameterSet.Config as cms

process = cms.Process("OfflineValidator") 
   
#process.load("Alignment.OfflineValidation.DataSetCRAFT09_109011_109624_cff")
#process.load("Alignment.OfflineValidation.DataSetCRAFT09_38Tpeak_TP_cff")
#process.load("Alignment.OfflineValidation.DataSetCRAFT09_38Tdec_TP_cff")

#process.source.inputCommands = cms.untracked.vstring('keep *', 'drop *_MEtoEDMConverter_*_*') # hack to get rid of the memory consumption problem in 2_2_X and beond
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(False),
   #Rethrow = cms.untracked.vstring("ProductNotFound"), # make this exception fatal
   fileMode  =  cms.untracked.string('NOMERGE') # no ordering needed, but calls endRun/beginRun etc. at file boundaries
)

### Maximum number of Events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
 )

### Output File Configuration
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('temp.root'),
    #fileName = cms.string('/tmp/eaguiloc/alignmentValidation_090902_112207/AlignmentValidation_HIPPeak.root'),
    closeFileFast = cms.untracked.bool(True)
 )
#process.TFileService.closeFileFast = True

### Messages & Convenience
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger",
#    destinations = cms.untracked.vstring('LOGFILE_Offline_HIPPeak', 
#        'cout')
#)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
reportEvery = cms.untracked.int32(1000) # every 1000th only
#    limit = cms.untracked.int32(10)       # or limit to 10 printouts...
))
process.MessageLogger.statistics.append('cout')

 ## report only every 100th record
 ## process.MessageLogger.cerr.FwkReport.reportEvery = 100

    
### Alignment Track Selection
process.load("Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi")
process.AlignmentTrackSelector.src = 'HitFilteredTracks'
process.AlignmentTrackSelector.filter = True
process.AlignmentTrackSelector.applyBasicCuts = True
process.AlignmentTrackSelector.pMin    = 4.
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
   

####  new FILTER
#-- new track hit filter
# TrackerTrackHitFilter takes as input the tracks/trajectories coming out from TrackRefitter1
process.load("RecoTracker.FinalTrackSelectors.TrackerTrackHitFilter_cff")
process.TrackerTrackHitFilter.src = 'TrackRefitter1'
process.TrackerTrackHitFilter.useTrajectories= True  # this is needed only if you require some selections; but it will work even if you don't ask for them
process.TrackerTrackHitFilter.minimumHits = 8
process.TrackerTrackHitFilter.commands = cms.vstring("keep PXB","keep PXE","keep TIB","keep TID","keep TOB","keep TEC")
process.TrackerTrackHitFilter.detsToIgnore = [
]
process.TrackerTrackHitFilter.replaceWithInactiveHits = True
process.TrackerTrackHitFilter.stripAllInvalidHits = False
process.TrackerTrackHitFilter.rejectBadStoNHits = True
process.TrackerTrackHitFilter.StoNcommands = cms.vstring("ALL 8.0")
#process.TrackerTrackHitFilter.StoNcommands = cms.vstring("ALL 14.0")
process.TrackerTrackHitFilter.rejectLowAngleHits= True
process.TrackerTrackHitFilter.TrackAngleCut= 0.35 # in rads, starting from the module surface
process.TrackerTrackHitFilter.usePixelQualityFlag= True

#now we give the TrackCandidate coming out of the TrackerTrackHitFilter to the track producer
import RecoTracker.TrackProducer.CTFFinalFitWithMaterialP5_cff
process.HitFilteredTracks = RecoTracker.TrackProducer.CTFFinalFitWithMaterialP5_cff.ctfWithMaterialTracksP5.clone(
    src = 'TrackerTrackHitFilter',
###    TrajectoryInEvent = True,
    TTRHBuilder = "WithAngleAndTemplate"    
)

 ##
 ## Load and Configure TrackRefitter1
 ##

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

#############3
# parameters for TrackRefitter
#process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff")
import RecoTracker.TrackProducer.TrackRefitters_cff
process.TrackRefitter1 = process.TrackRefitterP5.clone()
#process.TrackRefitter1 = process.TrackRefitter.clone()
#process.TrackRefitter1.src = 'ALCARECOTkAlCosmicsCTF0T'
process.TrackRefitter1.TrajectoryInEvent = True
process.TrackRefitter1.TTRHBuilder = "WithAngleAndTemplate"

process.TrackRefitter2 = process.TrackRefitter1.clone(
#    src = 'HitFilteredTracks'
    src = 'AlignmentTrackSelector',
    )


### Get the BeamSpot
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
 
### GlobalTag Conditions (if needed)
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "CRAFT09_R_V4::All"
#process.GlobalTag.globaltag = "GR09_31X_V4P::All"
#process.GlobalTag.connect="frontier://FrontierProd/CMS_COND_31X_FROM21X"
#process.GlobalTag.connect="frontier://FrontierProd/CMS_COND_31X_GLOBALTAG"

## LAYERWISE Lorentz Angle ###################

#process.SiStripLorentzAngle = cms.ESSource("PoolDBESSource",
#     BlobStreamerName = 
#cms.untracked.string('TBufferBlobStreamingService'),
#     DBParameters = cms.PSet(
#         messageLevel = cms.untracked.int32(2),
#         authenticationPath = 
#         cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
#     ),
#     timetype = cms.string('runnumber'),
#     toGet = cms.VPSet(cms.PSet(
#         record = cms.string('SiStripLorentzAngleRcd'),
#        tag = cms.string('SiStripLA_CRAFT_layers')
#     )),
#     connect = cms.string('sqlite_file:/afs/cern.ch/user/j/jdraeger/public/LA_object/LA_CRAFT_layers.db')
#)
#process.es_prefer_SiStripLorentzAngle = cms.ESPrefer("PoolDBESSource","SiStripLorentzAngle")
  
### Geometry
process.load("Configuration.StandardSequences.Geometry_cff")
 
### Magnetic Field
process.load("Configuration/StandardSequences/MagneticField_38T_cff")


from CondCore.DBCommon.CondDBSetup_cfi import *
process.ZeroAPE = cms.ESSource("PoolDBESSource",CondDBSetup,
                                        connect = cms.string('frontier://FrontierProd/CMS_COND_31X_FROM21X'),
                                        timetype = cms.string("runnumber"),
                                        toGet = cms.VPSet(
                                                          cms.PSet(record = cms.string('TrackerAlignmentErrorRcd'),
                                                                   tag = cms.string('TrackerIdealGeometryErrors210_mc')
                                                                   ))
                                        )
process.es_prefer_ZeroAPE = cms.ESPrefer("PoolDBESSource", "ZeroAPE")



from CondCore.DBCommon.CondDBSetup_cfi import *
process.trackerAlignment = cms.ESSource("PoolDBESSource",CondDBSetup,
                                        #PEAK
                                        #connect = cms.string('sqlite_file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/PayLoads/CRAFT09/TrackerAlignment_2009_v1_prompt/111462-infty/TrackerAlignment_2009_v1_prompt.db'),
                                        connect = cms.string('sqlite_file:TrackerAlignment_2009_v1_prompt.db'),
                                        #DECO
                                        #connect = cms.string('sqlite_file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/eaguiloc/ALCA_merged_3/alignments.db'),
                                        #connect = cms.string('sqlite_file:/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/bonato/CMSSW_3_2_2_patch2/src/Alignment/HIPAlignmentAlgorithm/HIPv4_Step2_SSTduPXLdets_smallAPE/alignments.db'),
                                        timetype = cms.string("runnumber"),
                                        toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
                                                                   tag = cms.string('Alignments')
                                                                   ))
                                        )
process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")

#process.trackerAlignment = cms.ESSource("PoolDBESSource",
#                                        CondDBSetup,
#                                        timetype = cms.string('runnumber'),
#                                        toGet = cms.VPSet(cms.PSet(record = cms.string('TrackerAlignmentRcd'),
#                                                                   tag = cms.string('TrackerAlignment_2009_v1_prompt')
#                                                                   )),
#                                        connect = cms.string('frontier://FrontierProd/CMS_COND_31X')
#                                        )
#process.es_prefer_trackerAlignment = cms.ESPrefer("PoolDBESSource", "trackerAlignment")


## to apply misalignments
#TrackerDigiGeometryESModule.applyAlignment = True
   
### Load and Configure OfflineValidation
process.load("Alignment.OfflineValidation.TrackerOfflineValidation_cfi")
process.TrackerOfflineValidation.Tracks = 'TrackRefitter2'
process.TrackerOfflineValidation.trajectoryInput = 'TrackRefitter2'
process.TrackerOfflineValidation.moduleLevelHistsTransient = cms.bool(True)
process.TrackerOfflineValidation.localCoorHistosOn = cms.bool(True)
process.TrackerOfflineValidation.bookTH1 = cms.bool(True)
process.TrackerOfflineValidation.bookTH2 = cms.bool(False)
process.TrackerOfflineValidation.debug = cms.bool(False)
process.TrackerOfflineValidation.fillTree = cms.bool(False)

# Normalized X Residuals, normal local coordinates (Strip)
process.TrackerOfflineValidation.TH1NormXResStripModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-3.0), xmax = cms.double(3.0)
)

# X Residuals, normal local coordinates (Strip)                      
process.TrackerOfflineValidation.TH1XResStripModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-0.5), xmax = cms.double(0.5)
)

# Normalized X Residuals, native coordinates (Strip)
process.TrackerOfflineValidation.TH1NormXprimeResStripModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-3.0), xmax = cms.double(3.0)
)

# X Residuals, native coordinates (Strip)
process.TrackerOfflineValidation.TH1XprimeResStripModules = cms.PSet(
#    Nbinx = cms.int32(2000), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    Nbinx = cms.int32(200), xmin = cms.double(-1000), xmax = cms.double(1000)
)

# X Residuals vs Theta, native coordinates (Strip)
process.TrackerOfflineValidation.TH2XprimevsThetaStripModules = cms.PSet(
    #    Nbinx = cms.int32(2000), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    Nbinx = cms.int32(100), xmin = cms.double(-1000), xmax = cms.double(1000)
)

# X Residuals vs Theta, native coordinates (Pixel)
process.TrackerOfflineValidation.TH2XprimevsThetaPixelModules = cms.PSet(
    #    Nbinx = cms.int32(2000), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    Nbinx = cms.int32(1), xmin = cms.double(-1000), xmax = cms.double(1000)
)

# Normalized Y Residuals, native coordinates (Strip -> hardly defined)
process.TrackerOfflineValidation.TH1NormYResStripModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-3.0), xmax = cms.double(3.0)
)
# -> very broad distributions expected                                         
process.TrackerOfflineValidation.TH1YResStripModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-10.0), xmax = cms.double(10.0)
)

# Normalized X residuals normal local coordinates (Pixel)                                        
process.TrackerOfflineValidation.TH1NormXResPixelModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-3.0), xmax = cms.double(3.0)
)
# X residuals normal local coordinates (Pixel)                                        
process.TrackerOfflineValidation.TH1XResPixelModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-0.5), xmax = cms.double(0.5)
)
# Normalized X residuals native coordinates (Pixel)                                        
process.TrackerOfflineValidation.TH1NormXprimeResPixelModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-3.0), xmax = cms.double(3.0)
)
# X residuals native coordinates (Pixel)                                        
process.TrackerOfflineValidation.TH1XprimeResPixelModules = cms.PSet(
#    Nbinx = cms.int32(2000), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    Nbinx = cms.int32(1), xmin = cms.double(-0.1), xmax = cms.double(0.1)
)                                        
# Normalized Y residuals native coordinates (Pixel)                                         
process.TrackerOfflineValidation.TH1NormYResPixelModules = cms.PSet(
    Nbinx = cms.int32(1), xmin = cms.double(-3.0), xmax = cms.double(3.0)
)
# Y residuals native coordinates (Pixel)                                         
process.TrackerOfflineValidation.TH1YResPixelModules = cms.PSet(
#    Nbinx = cms.int32(2000), xmin = cms.double(-0.5), xmax = cms.double(0.5)
    Nbinx = cms.int32(1), xmin = cms.double(-0.1), xmax = cms.double(0.1)
)

###Remove pixel template read-in from DB
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi")
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPETemplateReco_cfi")
process.templates.LoadTemplatesFromDB = False

### PATH
process.p = cms.Path(process.offlineBeamSpot*process.TrackRefitter1*process.TrackerTrackHitFilter*process.HitFilteredTracks
                     *process.AlignmentTrackSelector*process.TrackRefitter2*process.TrackerOfflineValidation)

