
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "CLHEP/Vector/RotationInterfaces.h" 

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"
#include "RecoTracker/TrackProducer/interface/TrackingRecHitLessFromGlobalPosition.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"  
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"

#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiStripLorentzAngleRcd.h"

//#include "CalibTracker/SiStripCommon/interface/ShallowTrackClustersProducer.h"
#include "CalibTracker/SiStripCommon/interface/ShallowTools.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"

bool debug_=false;

TrackerValidationVariables::TrackerValidationVariables(){}

int TrackerValidationVariables::getHisto(float z){
  if(z<-92.)             return 0;
  if(z>=-92. && z<-73.)  return 1;
  if(z>=-73. && z<-55.)  return 2;
  if(z>=-55. && z<-37.)  return 3;
  if(z>=-37. && z<-18.)  return 4;
  if(z>=-18. && z<-0.)   return 5;
  if(z>=0.   && z<18.)   return 6;
  if(z>=18.  && z<37.)   return 7;
  if(z>=37.  && z<55.)   return 8;
  if(z>=55.  && z<73.)   return 9;
  if(z>=73.  && z<92.)   return 10;
  if(z>=92.)             return 11;
  return 0;
}

TrackerValidationVariables::TrackerValidationVariables(const edm::EventSetup& es, const edm::ParameterSet& iSetup) 
  : conf_(iSetup)
{
  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  es.get<IdealMagneticFieldRecord>().get(magneticField_);
  es.get<SiStripLorentzAngleRcd>().get(SiStripLorentzAngle_);  
}

TrackerValidationVariables::~TrackerValidationVariables() {}

//void TrackerValidationVariables::fillHitQuantities(const edm::Event& iEvent, std::vector<AVHitStruct>& v_avhitout ){
void TrackerValidationVariables::fillHitQuantities(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<AVHitStruct>& v_avhitout ){ 

  edm::Handle<std::vector<Trajectory> > trajCollectionHandle;
  iEvent.getByLabel(conf_.getParameter<std::string>("trajectoryInput"),trajCollectionHandle);
  
  TrajectoryStateCombiner tsoscomb;
  LogDebug("TrackerValidationVariables") << "trajColl->size(): " << trajCollectionHandle->size() ;

  int nit=0;

  //cout<<"trajColl->size() "<<trajCollectionHandle->size()<<endl;


  for(std::vector<Trajectory>::const_iterator it = trajCollectionHandle->begin(), itEnd = trajCollectionHandle->end(); 
      it!=itEnd;++it){

    //cout<<"Trajectory "<<nit<<endl;
    nit++;
    int nitTraj=0;

    const std::vector<TrajectoryMeasurement> &tmColl = it->measurements();
    for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = tmColl.begin(), itTrajEnd = tmColl.end(); 
	itTraj != itTrajEnd; ++itTraj) {

      //cout<<"Trajectory measurement "<<nitTraj<<endl;
      nitTraj++;

      if(! itTraj->updatedState().isValid()) continue;
      
      
      TrajectoryStateOnSurface tsos = tsoscomb( itTraj->forwardPredictedState(), itTraj->backwardPredictedState() );
      TransientTrackingRecHit::ConstRecHitPointer hit    = itTraj->recHit();
      
      if(! hit->isValid() || hit->geographicalId().det() != DetId::Tracker ) continue; 
      
      AVHitStruct hitStruct;
      const DetId& hit_detId = hit->geographicalId();
      uint IntRawDetID = (hit_detId.rawId());	
      uint IntSubDetID = (hit_detId.subdetId());
      
      if(IntSubDetID == 0) continue;
      
      float trackOrientation=1.;
      float globalZofunitlocalY=-100.;      

      //Added by Ben-------------------------------------------------------------------------------
      if(IntSubDetID == StripSubdetector::TIB || IntSubDetID == StripSubdetector::TOB || 
	 IntSubDetID == StripSubdetector::TID || IntSubDetID == StripSubdetector::TEC) {

	const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>( tkGeom_->idToDet( hit->geographicalId() ) );

	LocalVector drift = shallow::drift( theStripDet, *magneticField_, *SiStripLorentzAngle_);    
	float localtheta = (theStripDet->toLocal(tsos.globalDirection())).theta(); 
	float localphi   = (theStripDet->toLocal(tsos.globalDirection())).phi();   
	//cout<<"localtheta "<<localtheta<<" localphi "<<localphi<<endl;

	globalZofunitlocalY = (theStripDet->toGlobal(LocalVector(0,1,0))).z();

	float tanTrackAngle = tan(localtheta)*cos(localphi);
	float tanLorentzAngle = drift.x()/drift.z();
	if(tanTrackAngle<tanLorentzAngle) trackOrientation=-1.;

	//float tracktheta    = atan(tantracktheta);
	//string subdet;
	//if(IntSubDetID == StripSubdetector::TIB) subdet="TIB";
	//if(IntSubDetID == StripSubdetector::TOB) subdet="TOB";
	//if(IntSubDetID == StripSubdetector::TEC) subdet="TEC";
	//if(IntSubDetID == StripSubdetector::TID) subdet="TID";
	//cout<<subdet<<" tanTrackAngle "<<tanTrackAngle<<" tanLorentzAngle "<<tanLorentzAngle<<endl;

	hitStruct.tanTrackAngle   = tanTrackAngle;
	hitStruct.tanLorentzAngle = tanLorentzAngle;
	hitStruct.localtheta      = localtheta;
	hitStruct.localphi        = localphi;
	hitStruct.x               = theStripDet->toGlobal(hit->localPosition()).x();
	hitStruct.y               = theStripDet->toGlobal(hit->localPosition()).y();
	hitStruct.z               = theStripDet->toGlobal(hit->localPosition()).z();

	const SiStripRecHit2D* stripHit = dynamic_cast<const SiStripRecHit2D*> ( hit->hit() );   
	const SiStripCluster *const stripCluster = stripHit->cluster().operator->(); 
	  
	float charge=0;
	int nstrips=0;
	const std::vector<uint8_t>& amplitudes_ =  stripCluster->amplitudes();
	for(size_t i=0; i<amplitudes_.size();i++){
	  if (amplitudes_[i]>0){
	    charge+=amplitudes_[i];
	    nstrips++;
	  }
	  
	}

	hitStruct.charge  = charge;
	hitStruct.nstrips = nstrips;

	/*
	//printout for fireworks
	string det;
	if(IntSubDetID == StripSubdetector::TIB) det="TIB";
	if(IntSubDetID == StripSubdetector::TOB) det="TOB";
	if(IntSubDetID == StripSubdetector::TID) det="TID";
	if(IntSubDetID == StripSubdetector::TEC) det="TEC";
	float xhit   = theStripDet->toGlobal(hit->localPosition()).x();
	float yhit   = theStripDet->toGlobal(hit->localPosition()).y();
	float zhit   = theStripDet->toGlobal(hit->localPosition()).z();
	float rhit=sqrt(xhit+xhit+yhit*yhit);
	cout<<"Det "<<det
	    <<" id "<<IntSubDetID
	    <<" rawid "<<IntRawDetID
	    <<" r "<<rhit
	    <<" z "<<zhit
	    <<" ring "<<this->getHisto(zhit)<<endl;
	cout<<"theta "<<localtheta
	    <<" phi "<<localphi
	    <<" tanTrackAngle "<<tanTrackAngle
	    <<" tanLorentzAngle "<<tanLorentzAngle
	    <<" trackOrientation "<<trackOrientation<<endl;*/
      }
      //-------------------------------------------------------------------------------------------
      

      //first calculate residuals in cartesian coordinates in the local module coordinate system
      
      LocalPoint lPHit = hit->localPosition();
      LocalPoint lPTrk = tsos.localPosition();
      
      //LocalError errHit = hit->localPositionError();
      // adding APE to hitError
      AlgebraicROOTObject<2>::SymMatrix mat = asSMatrix<2>(hit->parametersError());
      LocalError errHit = LocalError( mat(0,0),mat(0,1),mat(1,1) );
      LocalError errTrk = tsos.localError().positionError();
      
      //check for negative error values: track error can have negative value, if matrix inversion fails (very rare case)
      //hit error should always give positive values
      if(errHit.xx()<0. || errHit.yy()<0. || errTrk.xx()<0. || errTrk.yy()<0.){
        edm::LogError("Negative error Value")<<"@SUB=TrackerValidationVariables::fillHitQuantities"
	                                     <<"One of the squared error methods gives negative result"
					     <<"\n\terrHit.xx()\terrHit.yy()\terrTrk.xx()\terrTrk.yy()"
					     <<"\n\t"<<errHit.xx()<<"\t"<<errHit.yy()<<"\t"<<errTrk.xx()<<"\t"<<errTrk.yy();
	continue;
      }
      
      //align::LocalVector res = lPTrk - lPHit;
      align::LocalVector res = lPHit - lPTrk;
      
      float resXErr = std::sqrt( errHit.xx() + errTrk.xx() );
      float resYErr = std::sqrt( errHit.yy() + errTrk.yy() );
      
      hitStruct.resX = res.x();
      hitStruct.resY = res.y();
      hitStruct.resErrX = resXErr;
      hitStruct.resErrY = resYErr;
      
      // now calculate residuals taking global orientation of modules and radial topology in TID/TEC into account
      float resXprime(999.F), resYprime(999.F);
      float resXprimeErr(999.F), resYprimeErr(999.F);
      
      if(hit->detUnit()){ // is it a single physical module?
        const GeomDetUnit& detUnit = *(hit->detUnit());
	float uOrientation(-999.F), vOrientation(-999.F);
	float resXTopol(999.F), resYTopol(999.F);


	const Surface& surface = hit->detUnit()->surface();
	LocalPoint lPModule(0.,0.,0.), lUDirection(1.,0.,0.), lVDirection(0.,1.,0.);
	GlobalPoint gPModule    = surface.toGlobal(lPModule),
	            gUDirection = surface.toGlobal(lUDirection),
	            gVDirection = surface.toGlobal(lVDirection);
	
	if(IntSubDetID == PixelSubdetector::PixelBarrel || IntSubDetID == StripSubdetector::TIB || IntSubDetID == StripSubdetector::TOB) {
	  uOrientation = deltaPhi(gUDirection.phi(),gPModule.phi()) >= 0. ? +1.F : -1.F;
	  vOrientation = gVDirection.z() - gPModule.z() >= 0 ? +1.F : -1.F;
	  resXTopol = res.x();
	  resYTopol = res.y();
	  resXprimeErr = resXErr;
	  resYprimeErr = resYErr;
	  //if(StripSubdetector::TIB || IntSubDetID == StripSubdetector::TOB){
	  //cout<<"Barrel: uOrientation "<<uOrientation<<" vOrientation "<<vOrientation<<" globalZofunitlocalY "<<globalZofunitlocalY<<endl;
	  //}

	} else if (IntSubDetID == PixelSubdetector::PixelEndcap) {
	  uOrientation = gUDirection.perp() - gPModule.perp() >= 0 ? +1.F : -1.F;
	  vOrientation = deltaPhi(gVDirection.phi(),gPModule.phi()) >= 0. ? +1.F : -1.F;
	  resXTopol = res.x();
	  resYTopol = res.y();
	  resXprimeErr = resXErr;
	  resYprimeErr = resYErr;
	} else if (IntSubDetID == StripSubdetector::TID || IntSubDetID == StripSubdetector::TEC) {
	  uOrientation = deltaPhi(gUDirection.phi(),gPModule.phi()) >= 0. ? +1.F : -1.F;
	  vOrientation = gVDirection.perp() - gPModule.perp() >= 0. ? +1.F : -1.F;
	
	  //cout<<"Endcap: uOrientation "<<uOrientation<<" vOrientation "<<vOrientation<<" globalZofunitlocalY "<<globalZofunitlocalY<<endl;
	  	  
	  if(!dynamic_cast<const RadialStripTopology*>(&detUnit.topology()))continue;
	  const RadialStripTopology& topol = dynamic_cast<const RadialStripTopology&>(detUnit.topology());
	  
	  MeasurementPoint measHitPos = topol.measurementPosition(lPHit);
          MeasurementPoint measTrkPos = topol.measurementPosition(lPTrk);
	  
	  MeasurementError measHitErr = topol.measurementError(lPHit,errHit);
          MeasurementError measTrkErr = topol.measurementError(lPTrk,errTrk);
	  
	  if(measHitErr.uu()<0. || measHitErr.vv()<0. || measTrkErr.uu()<0. || measTrkErr.vv()<0.){
	    edm::LogError("Negative error Value")<<"@SUB=TrackerValidationVariables::fillHitQuantities"
	                                         <<"One of the squared error methods gives negative result"
	                                         <<"\n\tmeasHitErr.uu()\tmeasHitErr.vv()\tmeasTrkErr.uu()\tmeasTrkErr.vv()"
		                                 <<"\n\t"<<measHitErr.uu()<<"\t"<<measHitErr.vv()<<"\t"<<measTrkErr.uu()<<"\t"<<measTrkErr.vv();
	    continue;
	  }
	  
	  float localStripLengthHit = topol.localStripLength(lPHit);
	  float localStripLengthTrk = topol.localStripLength(lPTrk);
          float phiHit = topol.stripAngle(measHitPos.x());
	  float phiTrk = topol.stripAngle(measTrkPos.x());
          float r_0 = topol.originToIntersection();
	  
	  
	  resXTopol = (phiTrk-phiHit)*r_0;
	  //resYTopol = measTrkPos.y()*localStripLengthTrk - measHitPos.y()*localStripLengthHit;
	  float cosPhiHit(cos(phiHit)), cosPhiTrk(cos(phiTrk)),
	        sinPhiHit(sin(phiHit)), sinPhiTrk(sin(phiTrk));
	  float l_0 = r_0 - topol.detHeight()/2;
	  resYTopol = measTrkPos.y()*localStripLengthTrk - measHitPos.y()*localStripLengthHit + l_0*(1/cosPhiTrk - 1/cosPhiHit);
	  
	  
	  resXprimeErr = std::sqrt(measHitErr.uu()+measTrkErr.uu())*topol.angularWidth()*r_0;
          //resYprimeErr = std::sqrt(measHitErr.vv()*localStripLengthHit*localStripLengthHit + measTrkErr.vv()*localStripLengthTrk*localStripLengthTrk);
	  float helpSummand = l_0*l_0*topol.angularWidth()*topol.angularWidth()*(sinPhiHit*sinPhiHit/pow(cosPhiHit,4)*measHitErr.uu()
	                                                                       + sinPhiTrk*sinPhiTrk/pow(cosPhiTrk,4)*measTrkErr.uu() );
	  resYprimeErr = std::sqrt(measHitErr.vv()*localStripLengthHit*localStripLengthHit
	                         + measTrkErr.vv()*localStripLengthTrk*localStripLengthTrk + helpSummand );
	  
	} else {
	  edm::LogWarning("TrackerValidationVariables") << "@SUB=TrackerValidationVariables::fillHitQuantities" 
	                                                << "No valid tracker subdetector " << IntSubDetID;
	  continue;
	}
	
	//vOrientation=-1.;
	float wOrientation=uOrientation*vOrientation;
	//cout<<"orientation: u "<<uOrientation<<" v "<<vOrientation<<" w "<<wOrientation<<endl;
	resXprime = resXTopol*trackOrientation;
	resYprime = resYTopol*vOrientation;
	hitStruct.uOrientation = uOrientation;
	hitStruct.vOrientation = vOrientation;
	hitStruct.wOrientation = wOrientation;

	//Read in DT/ECAL/HCAL time
	timeStruct ts;
	setTime(iEvent,ts);
	if(debug_){
	  cout<<"Number of muons "<<ts.nmu<<endl;
	  cout<<"DT time   "<<ts.dttime<<"+/-"<<ts.dttimeerr<<" # DT measurements "<<ts.ndt<<endl;
	  cout<<"ECAL time "<<ts.ecaltime<<"+/-"<<ts.ecaltimeerr<<" for ECAL energy "<<ts.ecalenergy<<endl;
	  cout<<"HCAL time "<<ts.hcaltime<<"+/-"<<ts.hcaltimeerr<<" for ECAL energy "<<ts.hcalenergy<<endl;
	}

	hitStruct.dttime     = ts.dttime;
	hitStruct.dttimeerr  = ts.dttimeerr;
	hitStruct.ndt        = ts.ndt;
	hitStruct.hcaltime   = ts.hcaltime;
	hitStruct.nvalidmu   = ts.nvalidmu;

	float dusign=0.;
	if(IntSubDetID == StripSubdetector::TOB){
	  dusign=1.;
	  const DetId& id = hit->geographicalId();
	  if((id&0x1e000000)==0x1a000000) {
	    int mod =(id&0x0000001c) >> 2;
	    int side = (id&0x00003000) >> 12;
	    if(side==1 && (mod == 5 || mod == 6)) dusign = -1.;
	    if(side==2 && (mod == 1 || mod == 2 || mod == 3 || mod == 4)) dusign = -1.;
	  }
	  //cout<<"u "<<uOrientation<<" v "<<vOrientation<<" w "<<uOrientation*vOrientation<<" dusign "<<dusign<<endl;
	  //float duumcorr = duum * dusign;
	}
	hitStruct.dusign = dusign;

	if(IntSubDetID == StripSubdetector::TOB && debug_){
	  cout<<"-------------------------------------------"<<endl;
	  cout<<"z "<<hitStruct.z<<" zbin "<<this->getHisto(hitStruct.z)<<" vOrientation "<<vOrientation
	      <<" dusign "<<dusign<<" trackOrientation "<<trackOrientation<<" res "<<resXprime<<endl;
	}
	
      }else{ // not a detUnit, so must be a virtual 2D-Module
	//FIXME: at present only for det units residuals are calculated and filled in the hitStruct
	// But in principle this method should also be useable for for the gluedDets (2D modules in TIB, TID, TOB, TEC)
	// In this case, only orientation should be taken into account for primeResiduals, but not the radial topology
      }
      
      hitStruct.resXprime = resXprime;
      hitStruct.resYprime = resYprime;
      hitStruct.resXprimeErr = resXprimeErr;
      hitStruct.resYprimeErr = resYprimeErr;

      hitStruct.rawDetId = IntRawDetID;
      hitStruct.phi = tsos.globalDirection().phi();
      hitStruct.eta = tsos.globalDirection().eta();

      if(IntSubDetID == StripSubdetector::TOB && debug_){
	cout<<"z "<<hitStruct.z<<" zbin "<<this->getHisto(hitStruct.z)<<" vOrientation "<<hitStruct.vOrientation
	    <<" dusign "<<hitStruct.dusign<<" res "<<hitStruct.resXprime<<endl;
      }

      //cout<<"theta "<<tsos.globalDirection().theta()<<" phi "<<tsos.globalDirection().phi()<<endl<<endl;
      
      
      // first try for overlapp residuals
      // based on Code from Keith and Wolfgang
      if(itTraj+1 != itTrajEnd) {
	TransientTrackingRecHit::ConstRecHitPointer hit2 = (itTraj+1)->recHit();
	TrackerAlignableId ali1, ali2;
	if(hit2->isValid() && 
	   ali1.typeAndLayerFromDetId(hit->geographicalId()) == ali2.typeAndLayerFromDetId(hit2->geographicalId())  &&
	   hit2->geographicalId().rawId() !=  SiStripDetId::SiStripDetId(IntRawDetID).partnerDetId()  
	   ) {	    
	  
	  float overlapPath_;
	  TrajectoryStateCombiner combiner_;
	  AnalyticalPropagator propagator(&(*magneticField_));
	  // forward and backward predicted states at module 1
	  TrajectoryStateOnSurface fwdPred1 = (itTraj)->forwardPredictedState();
	  TrajectoryStateOnSurface bwdPred1 = (itTraj)->backwardPredictedState();
	  if ( !fwdPred1.isValid() || !bwdPred1.isValid() )  continue;
	  // backward predicted state at module 2
	  TrajectoryStateOnSurface bwdPred2 = (itTraj+1)->backwardPredictedState();
	  TrajectoryStateOnSurface fwdPred2 = (itTraj+1)->forwardPredictedState();
	  if ( !bwdPred2.isValid() )  continue;
	  // extrapolation bwdPred2 to module 1
	  TrajectoryStateOnSurface bwdPred2At1 = propagator.propagate(bwdPred2,fwdPred1.surface());
	  if ( !bwdPred2At1.isValid() )  continue;
	  // combination with fwdPred1 (ref. state, best estimate without hits 1 and 2)
	  TrajectoryStateOnSurface comb1 = combiner_.combine(fwdPred1,bwdPred2At1);
	  if ( !comb1.isValid() )  continue;
	  
	  //
	  // propagation of reference parameters to module 2
	  //
	  std::pair<TrajectoryStateOnSurface,double> tsosWithS =
	    propagator.propagateWithPath(comb1,bwdPred2.surface());
	  TrajectoryStateOnSurface comb1At2 = tsosWithS.first;
	  
	  // Alternative possibility, not used at present
	  //TrajectoryStateOnSurface comb1At2 = propagator.propagate(comb1,bwdPred2.surface());
	  
	  if ( !comb1At2.isValid() )  continue;
	  overlapPath_ = tsosWithS.second;
	  
	  std::vector<GlobalPoint> predictedPositions;
	  predictedPositions.push_back(comb1.globalPosition());
	  predictedPositions.push_back(comb1At2.globalPosition());
	  
	  GlobalVector diff_pred = predictedPositions[0] - predictedPositions[1];
	  
	  TrajectoryStateOnSurface tsos2 = tsoscomb( (itTraj+1)->forwardPredictedState(), (itTraj+1)->backwardPredictedState() );
	  align::LocalVector res2 = tsos2.localPosition() - hit2->localPosition();
	  //float overlapresidual = res2.x() - res.x();
	  float overlapresidual = diff_pred.x();
	  
	  hitStruct.overlapres = std::make_pair(hit2->geographicalId().rawId(),overlapresidual);
	}
      }
      
      v_avhitout.push_back(hitStruct);
      
    } 
  }  
}

void TrackerValidationVariables::setTime(const edm::Event& iEvent, timeStruct &ts){
  
  edm::Handle<reco::MuonCollection> muH;
  //try {
  // iEvent.getByLabel("muonsWitht0Correction",muH);
  iEvent.getByLabel("muons1Leg",muH);
  //}catch (...) {;}
  const reco::MuonCollection& muonsT0 = *(muH.product()); 
  
  if(debug_) cout << endl<< "Event "<<iEvent.id().event()<<" Number of muons = " << muonsT0.size() << "-----------------------------------------"<<endl;

  //reco::MuonCollection muonsT0filtered;
  float time = -9999.;
  float time_error = -9999.;
  float time_ndof = -9999.;
  float energy_ecal=-9999.;
  float time_ecal = -9999.; 
  float timeerr_ecal = -9999.;
  float energy_hcal=-9999.;
  float time_hcal = -9999.;
  float timeerr_hcal = -9999.;
  int nvalidmu=0;
  int imuon=0;
  //float p=-9999.;

  if(muonsT0.size()>0){
    for (unsigned int i=0; i<muonsT0.size(); i++) {

      //DT time
      reco::MuonTime mt0 = muonsT0[i].time();
      time = mt0.timeAtIpInOut;
      time_error = mt0.timeAtIpInOutErr;
      time_ndof = mt0.nDof;
      if(time_ndof>0){
	nvalidmu++;
	imuon=i;
      }

      //ECAL time
      energy_ecal  = muonsT0[i].calEnergy().em;	
      time_ecal    = muonsT0[i].calEnergy().ecal_time;
      timeerr_ecal = 0.;    //timeerr_ecal = muonsT0[i].calEnergy().ecal_timeError;

     //HCAL time
      energy_hcal  = muonsT0[i].calEnergy().had;	
      time_hcal    = muonsT0[i].calEnergy().hcal_time;
      timeerr_hcal = 0.;    //timeerr_hcal = muonsT0[i].calEnergy().hcal_timeError;
      
      if(debug_){
	cout<<"***** Muon" << i << " *****"<<endl;
	cout<<"DT time   "<<time<<"+/-"<<time_error<<" for "<<time_ndof<<" DT measurements"<<endl;
	cout<<"ECAL time "<<time_ecal<<"+/-"<<timeerr_ecal<<" for ECAL energy "<<energy_ecal<<endl;
	cout<<"HCAL time "<<time_hcal<<"+/-"<<timeerr_hcal<<" for HCAL energy "<<energy_hcal<<endl;
      }
    }

    reco::MuonTime mt0 = muonsT0[imuon].time();
    time = mt0.timeAtIpInOut;
    time_error = mt0.timeAtIpInOutErr;
    time_ndof = mt0.nDof;

    energy_ecal = muonsT0[imuon].calEnergy().em;	
    time_ecal = muonsT0[imuon].calEnergy().ecal_time;
    timeerr_ecal = 0.;     //timeerr_ecal = muonsT0[imuon].calEnergy().ecal_timeError;
    energy_hcal = muonsT0[imuon].calEnergy().had;	
    time_hcal = muonsT0[imuon].calEnergy().hcal_time;
    timeerr_hcal = 0.;     //timeerr_hcal = muonsT0[imuon].calEnergy().hcal_timeError;
    
    //p=muonsT0[imuon].muon_p;
    
  }

  ts.nmu              = muonsT0.size();
  ts.dttime           = time;
  ts.dttimeerr        = time_error;
  ts.ndt              = (int)time_ndof;
  ts.ecaltime         = time_ecal;
  ts.ecaltimeerr      = timeerr_ecal;
  ts.ecalenergy       = energy_ecal;
  ts.hcaltime         = time_hcal;
  ts.hcaltimeerr      = timeerr_hcal;
  ts.hcalenergy       = energy_hcal;
  ts.nvalidmu         = nvalidmu;
  //ts.p                = p;

}

void TrackerValidationVariables::fillTrackQuantities(const edm::Event& iEvent,
					std::vector<AVTrackStruct>& v_avtrackout)
{
  // get track collection from the event
  edm::InputTag TkTag = conf_.getParameter<edm::InputTag>("Tracks");
  edm::Handle<reco::TrackCollection> RecoTracks;
  iEvent.getByLabel(TkTag,RecoTracks);
  LogDebug("TrackerValidationVariables")<<"track collection size "<< RecoTracks->size();
  
  // Put here all track based quantities such as eta, phi, pt,.... 
  int i=0;
  for( reco::TrackCollection::const_iterator RecoTrack = RecoTracks->begin(), RecoTrackEnd = RecoTracks->end();
       RecoTrack !=RecoTrackEnd ; ++i, ++RecoTrack) {
    AVTrackStruct trackStruct;
    trackStruct.pt = RecoTrack->pt();
    trackStruct.ptError = RecoTrack->ptError();
    trackStruct.px = RecoTrack->px();
    trackStruct.py = RecoTrack->py();
    trackStruct.pz = RecoTrack->pz();
    trackStruct.eta = RecoTrack->eta();
    trackStruct.phi = RecoTrack->phi();
    trackStruct.chi2 = RecoTrack->chi2();
    trackStruct.normchi2 = RecoTrack->normalizedChi2();
    GlobalPoint gPoint(RecoTrack->vx(), RecoTrack->vy(), RecoTrack->vz());
    double theLocalMagField = magneticField_->inTesla(gPoint).z();
    trackStruct.kappa = -RecoTrack->charge()*0.002998*theLocalMagField/RecoTrack->pt();
    trackStruct.charge = RecoTrack->charge();
    trackStruct.d0 = RecoTrack->d0();
    trackStruct.dz = RecoTrack->dz();
    v_avtrackout.push_back(trackStruct);
  }

}











/*
void TrackerValidationVariables::setHCALTime(const edm::Event& iEvent, std::vector<float> &x){
 
  x.clear();
  edm::Handle<reco::MuonCollection> muH;
  //try {
  // iEvent.getByLabel("muonsWitht0Correction",muH);
  iEvent.getByLabel("muons1Leg",muH);
  //}catch (...) {;}
  const reco::MuonCollection& muonsT0 = *(muH.product()); 
  
  cout << endl<< "Event "<<iEvent.id().event()<<" Number of muons = " << muonsT0.size() << "-----------------------------------------"<<endl;
  //if(muonsT0.size()<1) return -1000.;  
  //reco::MuonCollection muonsT0filtered;
  //float e3x3_ecal = -9999.;

  float energy_ecal = -9999.; 
 
  float time_ecal = -9999.; 
  float timeerr_ecal = -9999.; 
  float time_hcal = -9999.;
  float timeerr_hcal = -9999.;
 
  if(muonsT0.size()>0){
    for (unsigned int i=0; i<muonsT0.size(); i++) {
      
      bool hasCaloEnergyInfo = muonsT0[i].isEnergyValid();
      if (hasCaloEnergyInfo) {
	energy_ecal = muonsT0[i].calEnergy().em;	
	time_ecal = muonsT0[i].calEnergy().ecal_time;
	timeerr_ecal = muonsT0[i].calEnergy().ecal_timeError;
	
	//e3x3_ecal = muonsT0[i].calEnergy().emS9;
	time_hcal = muonsT0[i].calEnergy().hcal_time;
      }
      
      if(hasCaloEnergyInfo)cout<<"****Muon "<<i<<" ecal time "<<time_ecal<<" ecal energy "<<energy_ecal<<" hcal time "<<time_hcal<<endl;
      else                 cout<<"no CaloEnergyInfo"<<endl;
    }
    
    bool hasCaloEnergyInfo = muonsT0[0].isEnergyValid();
    if (hasCaloEnergyInfo) {
      time_ecal = muonsT0[0].calEnergy().ecal_time;
      energy_ecal = muonsT0[0].calEnergy().em;
      e3x3_ecal = muonsT0[0].calEnergy().emS9;
      time_hcal = muonsT0[0].calEnergy().hcal_time;
    }
  }

  x.push_back((float)muonsT0.size());
  x.push_back(time_ecal);
  x.push_back(energy_ecal);
  x.push_back(e3x3_ecal);
  x.push_back(time_hcal);
   
}
*/
