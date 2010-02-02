import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_123615.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124009.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124020.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124022.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124024.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124030.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124120.root',
'rfio:/castor/cern.ch/cms/store/user/emiglior/ALCARECO/08Jan10/TkAlMinBias_124230.root',
] );


secFiles.extend( [
               ] )
