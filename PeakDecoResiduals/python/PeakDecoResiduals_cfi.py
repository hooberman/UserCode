import FWCore.ParameterSet.Config as cms

# Put here the modules you want the cfg file to use,
# then include this file in your cfg file.
# i.e. in Validator.cfg replace 'module demo = Validator {} '
# with 'include "anlyzerDir/Validator/data/Validator.cfi" '.
# (Remember that filenames are case sensitive.)
PeakDecoResiduals = cms.EDFilter("PeakDecoResiduals",
                                 Tracks = cms.InputTag("TrackRefitter"),
                                 trajectoryInput           = cms.string('TrackRefitter'),
                                 debug = cms.bool(False),
                                 runOnCosmics = cms.bool(False)
                                 )


