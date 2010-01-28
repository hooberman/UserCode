// -*- C++ -*-
//
// Package:    TrackerOfflineValidation
// Class:      TrackerOfflineValidation
// 
/**\class TrackerOfflineValidation TrackerOfflineValidation.cc Alignment/Validator/src/TrackerOfflineValidation.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Erik Butz
//         Created:  Tue Dec 11 14:03:05 CET 2007
// $Id: TrackerOfflineValidation.cc,v 1.1 2010/01/27 13:48:18 benhoob Exp $
//
//


// system include files
#include <memory>
#include <map>
#include <sstream>
#include <fstream>
#include <math.h>
#include <utility>
#include <vector>
#include <iostream>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "dout.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"
#include "Alignment/OfflineValidation/interface/TkOffTreeVariables.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/CommonAlignment/interface/AlignableComposite.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Alignment/CommonAlignment/interface/Utilities.h"
#include "Alignment/CommonAlignment/interface/AlignableObjectId.h"
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"

#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//
using namespace std;

class TrackerOfflineValidation : public edm::EDAnalyzer {
public:
  explicit TrackerOfflineValidation(const edm::ParameterSet&);
  ~TrackerOfflineValidation();
  
  enum HistogrammType { XResidual, 
			NormXResidual, 
			XprimeResidual, 
			NormXprimeResidual, 
			YprimeResidual, 
			NormYprimeResidual,
			XprimevsThetaResidual};
  
private:

  
  struct ModuleHistos{
    ModuleHistos() :  
		      ResXprimeHisto(), 
		      ResXprimeHisto2(), 
		      tanTrackAngleHisto(), 
		      tanLorentzAngleHisto(),
		      deltaTanHisto(),
		      ResXprimeOverThetaHisto(), 
		      ResXprimeOverThetaHisto_thick(), 
		      ResXprimeOverThetaHisto_thin(), 
		      ResXprimeOverThetaHisto2(), 
		      ResXprimeVsThetaHisto(),
		      duvstrkangleHisto(),
		      //duvslorangleHisto(),
		      HitEtaHisto(), 
		      HitZHisto(),
		      uOrientationHisto(), 
		      vOrientationHisto(), 
		      wOrientationHisto(), 
		      localphiHisto(), 
		      localthetaHisto(),
		      dttimeHisto(), 
		      dttimeerrHisto(),
		      ndtHisto(),
		      hcaltimeHisto(),
		      duvsdttimeHisto(),
		      dwvsdttimeHisto(),
		      chargeHisto(),
		      nstripsHisto(),
		      chargevsdttimeHisto()
      

    {}

  
    TH1* ResXprimeHisto;
    TH1* ResXprimeHisto2[2];
    TH1* tanTrackAngleHisto;
    TH1* tanLorentzAngleHisto;
    TH1* deltaTanHisto;
    TH1* ResXprimeOverThetaHisto;
    TH1* ResXprimeOverThetaHisto_thick;
    TH1* ResXprimeOverThetaHisto_thin;
    TH1* ResXprimeOverThetaHisto2[2];
    TH2* ResXprimeVsThetaHisto;
    TH2* duvstrkangleHisto;
    TH2* duvslorangleHisto;
    TH1* HitEtaHisto;
  
    TH1* HitZHisto;
    TH1* uOrientationHisto;
    TH1* vOrientationHisto;
    TH1* wOrientationHisto;
  
    TH1* localphiHisto;
    TH1* localthetaHisto;
 
    TH1* dttimeHisto;
    TH1* dttimeerrHisto;
    TH1* ndtHisto;
    TH1* hcaltimeHisto;
    TH2* duvsdttimeHisto;
    TH2* dwvsdttimeHisto;

    TH1* chargeHisto;
    TH1* nstripsHisto;
    TH2* chargevsdttimeHisto;

  };


  // container struct to organize collection of histogramms during endJob
  struct SummaryContainer{
    SummaryContainer() : sumXResiduals_(), 
			 sumXResiduals2_(), 
			 summaryXResiduals_(), 
			 sumTanTrackAngle_(), 
			 sumTanLorentzAngle_(),
			 sumdeltaTan_(),
			 sumResXprimeOverTheta_(), 
			 sumResXprimeOverTheta_thick_(), 
			 sumResXprimeOverTheta_thin_(), 
			 sumResXprimeOverTheta2_(), 
			 sumResXprimeVsTheta_(),
			 sumduvstrkangle_(),
			 sumHitEta_(), 
			 sumHitZ_(),
			 sumuOrientation_(), 
			 sumvOrientation_(), 
			 sumwOrientation_(), 
			 sumlocalphi_(), 
			 sumlocaltheta_(),
			 sumdttime_(),
			 sumdttimeerr_(),
			 sumndt_(),
			 sumhcaltime_(),
			 sumduvsdttime_(),
			 sumdwvsdttime_(),
			 sumcharge_(),
			 sumnstrips_(),
			 sumchargevsdttime_()
    {}
    
    TH1* sumXResiduals_;
    TH1* sumXResiduals2_[2];
    TH1* summaryXResiduals_;
    TH1* sumTanTrackAngle_;
    TH1* sumTanLorentzAngle_;
    TH1* sumdeltaTan_;
    TH1* sumResXprimeOverTheta_;
    TH1* sumResXprimeOverTheta_thick_;
    TH1* sumResXprimeOverTheta_thin_;
    TH1* sumResXprimeOverTheta2_[2];
    TH2* sumResXprimeVsTheta_;
    TH2* sumduvstrkangle_;
    TH1 *sumHitEta_;
    TH1 *sumHitZ_;
    TH1* sumuOrientation_;
    TH1* sumvOrientation_;
    TH1* sumwOrientation_;
    TH1* sumlocalphi_;
    TH1* sumlocaltheta_;
    TH1* sumdttime_;
    TH1* sumdttimeerr_;
    TH1* sumndt_;
    TH1* sumhcaltime_;
    TH2* sumduvsdttime_;
    TH2* sumdwvsdttime_;
    TH1* sumcharge_;
    TH1* sumnstrips_;
    TH2* sumchargevsdttime_;
  };

  // 
  // ------------- private member function -------------
  // 
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void checkBookHists(const edm::EventSetup &setup);



  void bookGlobalHists(TFileDirectory &tfd);
  void bookDirHists(TFileDirectory &tfd, const Alignable& ali, const AlignableObjectId &aliobjid);
  void bookHists(TFileDirectory &tfd, const Alignable& ali, align::StructureType type, int i, 
		 const AlignableObjectId &aliobjid);
 
  void collateSummaryHists( TFileDirectory &tfd, const Alignable& ali, int i, 
			    const AlignableObjectId &aliobjid, 
			    std::vector<TrackerOfflineValidation::SummaryContainer > &v_levelProfiles);
  
  void fillTree(TTree &tree, const std::map<int, TrackerOfflineValidation::ModuleHistos> &moduleHist_, 
		TkOffTreeVariables &treeMem, const TrackerGeometry &tkgeom );
  
  TrackerOfflineValidation::SummaryContainer bookSummaryHists(TFileDirectory &tfd, 
							      const Alignable& ali, 
							      align::StructureType type, int i, 
							      const AlignableObjectId &aliobjid); 

  ModuleHistos& getHistStructFromMap(const DetId& detid); 

  bool isBarrel(uint32_t subDetId);
  bool isEndCap(uint32_t subDetId);
  bool isPixel(uint32_t subDetId);
  bool isDetOrDetUnit(align::StructureType type);

  TH1* bookTH1F(bool isTransient, TFileDirectory& tfd, const char* histName, const char* histTitle, 
		int nBinsX, double lowX, double highX);

  TH2* bookTH2F(bool isTransient, TFileDirectory& tfd, const char* histName, const char* histTitle, 
		int nBinsX, double lowX, double highX,int nBinsY, double lowY, double highY);

  void getBinning(uint32_t subDetId, TrackerOfflineValidation::HistogrammType residualtype, 
		  int &nBinsX, double &lowerBoundX, double &upperBoundX);

  void summarizeBinInContainer(int bin, SummaryContainer &targetContainer, 
			       SummaryContainer &sourceContainer);

  void summarizeBinInContainer(int bin, uint32_t subDetId, SummaryContainer &targetContainer, 
			       ModuleHistos &sourceContainer);

  void setSummaryBin(int bin, TH1* targetHist, TH1* sourceHist);
    
  int getHisto(float z);

  float Fwhm(const TH1* hist) const;
  std::pair<float,float> fitResiduals(TH1 *hist) const; //, float meantmp, float rmstmp);
  float getMedian( const TH1 *hist) const; 
 // From MillePedeAlignmentMonitor: Get Index for Arbitary vector<class> by name
  template <class OBJECT_TYPE>  
  int GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name);

  // ---------- member data ---------------------------


  const edm::ParameterSet parset_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  const TrackerGeometry *bareTkGeomPtr_; // ugly hack to book hists only once, but check 

  
  // parameters from cfg to steer
  bool lCoorHistOn_;
  bool moduleLevelHistsTransient_;
  bool overlappOn_;
  bool stripYResiduals_;
  bool useFwhm_;
  bool useFit_;
  bool useOverflowForRMS_;
  bool bookTH1_; 
  bool bookTH2_;
  bool debug_;
  bool fillTree_;
  //bool removePixel_;
  std::map< std::pair<uint32_t, uint32_t >, TH1*> hOverlappResidual;

  // a vector to keep track which pointers should be deleted at the very end
  std::vector<TH1*> vDeleteObjects_;
  std::vector<TH2*> vDeleteObjectsTH2_;

  // 
  std::vector<TH1*> vTrackHistos_;
  std::vector<TProfile*> vTrackProfiles_;
  std::vector<TH2*> vTrack2DHistos_;
  
  std::map<int,TrackerOfflineValidation::ModuleHistos> mPxbResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mPxeResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTibResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTidResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTobResiduals_;
  std::map<int,TrackerOfflineValidation::ModuleHistos> mTecResiduals_;

  ofstream ofile;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
template <class OBJECT_TYPE>  
int TrackerOfflineValidation::GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name)
{
  int result = 0;
  for (typename std::vector<OBJECT_TYPE*>::const_iterator iter = vec.begin(), iterEnd = vec.end();
       iter != iterEnd; ++iter, ++result) {
    if (*iter && (*iter)->GetName() == name) return result;
  }
  edm::LogError("Alignment") << "@SUB=TrackerOfflineValidation::GetIndex" << " could not find " << name;
  return -1;
}

int TrackerOfflineValidation::getHisto(float z){
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
//
// constructors and destructor
//
TrackerOfflineValidation::TrackerOfflineValidation(const edm::ParameterSet& iConfig)
  : parset_(iConfig), bareTkGeomPtr_(0), lCoorHistOn_(parset_.getParameter<bool>("localCoorHistosOn")),
    moduleLevelHistsTransient_(parset_.getParameter<bool>("moduleLevelHistsTransient")),
    overlappOn_(parset_.getParameter<bool>("overlappOn")), 
    stripYResiduals_(parset_.getParameter<bool>("stripYResiduals")), 
    useFwhm_(parset_.getParameter<bool>("useFwhm")),
    useFit_(parset_.getParameter<bool>("useFit")),
    useOverflowForRMS_(parset_.getParameter<bool>("useOverflowForRMS")),
    bookTH1_(parset_.getParameter<bool>("bookTH1")),
    bookTH2_(parset_.getParameter<bool>("bookTH2")),
    debug_(parset_.getParameter<bool>("debug")),
    fillTree_(parset_.getParameter<bool>("fillTree"))
  //,
  //removePixel_(parset_.getParameter<bool>("removePixel"))
  
{

  //string datasetName                = iConfig.getParameter<std::string>("datasetName");
  string datasetName="/Cosmics/CRAFT09-TrackingPointing-CRAFT09_R_V4_CosmicsSeq_v1/RAW-RECO";
  ofile.open("events.txt",ios::trunc);
  ofile<<datasetName<<endl;

  if(bookTH1_) cout<<"Booking TH1 histos"<<endl;
  else         cout<<"Not Booking TH1 histos"<<endl;
  if(bookTH2_) cout<<"Booking TH2 histos"<<endl;
  else         cout<<"Not Booking TH2 histos"<<endl;
  if(debug_)   cout<<"Print debug statements"<<endl;

 //now do what ever initialization is needed
}


TrackerOfflineValidation::~TrackerOfflineValidation()
{
  if(debug_) cout<<"TrackerOfflineValidation::~TrackerOfflineValidation"<<endl;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  for( std::vector<TH1*>::const_iterator it = vDeleteObjects_.begin(), itEnd = vDeleteObjects_.end();it != itEnd; ++it) delete *it;
  for( std::vector<TH2*>::const_iterator itTH2 = vDeleteObjectsTH2_.begin(), itEndTH2 = vDeleteObjectsTH2_.end();itTH2 != itEndTH2; ++itTH2) delete *itTH2;
}


//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void
TrackerOfflineValidation::checkBookHists(const edm::EventSetup &es)
{

  if(debug_) cout<<"TrackerOfflineValidation::checkBookHists"<<endl;

  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  const TrackerGeometry *newBareTkGeomPtr = &(*tkGeom_);
  if (newBareTkGeomPtr == bareTkGeomPtr_) return; // already booked hists, nothing changed

  if (!bareTkGeomPtr_) { // pointer not yet set: called the first time => book hists
    edm::Service<TFileService> fs;    
    AlignableObjectId aliobjid;
    
    // construct alignable tracker to get access to alignable hierarchy 
    AlignableTracker aliTracker(&(*tkGeom_));
    
    edm::LogInfo("TrackerOfflineValidation") << "There are " << newBareTkGeomPtr->detIds().size()
					     << " detUnits in the Geometry record";
    
    //
    // Book Histogramms for global track quantities
    TFileDirectory trackglobal = fs->mkdir("GlobalTrackVariables");  
    this->bookGlobalHists(trackglobal);
    
    // recursively book histogramms on lowest level
//     this->bookDirHists(static_cast<TFileDirectory&>(*fs), aliTracker, aliobjid);  
    this->bookDirHists(*fs, aliTracker, aliobjid);  
  } else { // histograms booked, but changed TrackerGeometry?
    edm::LogWarning("GeometryChange") << "@SUB=checkBookHists"
				      << "TrackerGeometry changed, but will not re-book hists!";
  }

  bareTkGeomPtr_ = newBareTkGeomPtr;

  if(debug_) cout<<"Exit TrackerOfflineValidation::checkBookHists"<<endl;
}


void 
TrackerOfflineValidation::bookGlobalHists(TFileDirectory &tfd )
{
  if(debug_) cout<<"TrackerOfflineValidation::bookGlobalHists"<<endl;

  vTrackHistos_.push_back(tfd.make<TH1F>("h_tracketa",
					 "Track #eta;#eta_{Track};Number of Tracks",
					 90,-3.,3.));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature",
					 "Curvature #kappa;#kappa_{Track};Number of Tracks",
					 100,-.05,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature_pos",
					 "Curvature |#kappa| Positive Tracks;|#kappa_{pos Track}|;Number of Tracks",
					 100,.0,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_curvature_neg",
					 "Curvature |#kappa| Negative Tracks;|#kappa_{neg Track}|;Number of Tracks",
					 100,.0,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_diff_curvature",
					 "Curvature |#kappa| Tracks Difference;|#kappa_{Track}|;# Pos Tracks - # Neg Tracks",
					 100,.0,.05));
  vTrackHistos_.push_back(tfd.make<TH1F>("h_chi2",
					 "#chi^{2};#chi^{2}_{Track};Number of Tracks",
					 500,-0.01,500.));	       
  vTrackHistos_.push_back(tfd.make<TH1F>("h_normchi2",
					 "#chi^{2}/ndof;#chi^{2}/ndof;Number of Tracks",
					 100,-0.01,10.));     
  vTrackHistos_.push_back(tfd.make<TH1F>("h_pt",
					 "p_{T}^{track};p_{T}^{track} [GeV];Number of Tracks",
					 100,0.,2500));           
  vTrackHistos_.push_back(tfd.make<TH1F>("h_ptResolution",
					 "#delta{p_{T}/p_{T}^{track}};#delta_{p_{T}/p_{T}^{track}};Number of Tracks",
					 100,0.,0.5));           

  vTrackProfiles_.push_back(tfd.make<TProfile>("p_d0_vs_phi",
					       "Transverse Impact Parameter vs. #phi;#phi_{Track};#LT d_{0} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_dz_vs_phi",
					       "Longitudinal Impact Parameter vs. #phi;#phi_{Track};#LT d_{z} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_d0_vs_eta",
					       "Transverse Impact Parameter vs. #eta;#eta_{Track};#LT d_{0} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_dz_vs_eta",
					       "Longitudinal Impact Parameter vs. #eta;#eta_{Track};#LT d_{z} #GT [cm]",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_chi2_vs_phi",
					       "#chi^{2} vs. #phi;#phi_{Track};#LT #chi^{2} #GT",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_normchi2_vs_phi",
					       "#chi^{2}/ndof vs. #phi;#phi_{Track};#LT #chi^{2}/ndof #GT",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_chi2_vs_eta",
					       "#chi^{2} vs. #eta;#eta_{Track};#LT #chi^{2} #GT",
					       100,-3.15,3.15));  
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_normchi2_vs_eta",
					       "#chi^{2}/ndof vs. #eta;#eta_{Track};#LT #chi^{2}/ndof #GT",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_kappa_vs_phi",
					       "#kappa vs. #phi;#phi_{Track};#kappa",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_kappa_vs_eta",
					       "#kappa vs. #eta;#eta_{Track};#kappa",
					       100,-3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_ptResolution_vs_phi",
					       "#delta_{p_{T}}/p_{T}^{track};#phi^{track};#delta_{p_{T}}/p_{T}^{track}",
					       100, -3.15,3.15));
  vTrackProfiles_.push_back(tfd.make<TProfile>("p_ptResolution_vs_eta",
					       "#delta_{p_{T}}/p_{T}^{track};#eta^{track};#delta_{p_{T}}/p_{T}^{track}",
					       100, -3.15,3.15));


  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_d0_vs_phi",
					   "Transverse Impact Parameter vs. #phi;#phi_{Track};d_{0} [cm]",
					   100, -3.15, 3.15, 100,-1.,1.) );
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_dz_vs_phi",
					   "Longitudinal Impact Parameter vs. #phi;#phi_{Track};d_{z} [cm]",
					   100, -3.15, 3.15, 100,-100.,100.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_d0_vs_eta",
					   "Transverse Impact Parameter vs. #eta;#eta_{Track};d_{0} [cm]",
					   100, -3.15, 3.15, 100,-1.,1.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_dz_vs_eta",
					   "Longitudinal Impact Parameter vs. #eta;#eta_{Track};d_{z} [cm]",
					   100, -3.15, 3.15, 100,-100.,100.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_chi2_vs_phi",
					   "#chi^{2} vs. #phi;#phi_{Track};#chi^{2}",
					   100, -3.15, 3.15, 500, 0., 500.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_normchi2_vs_phi",
					   "#chi^{2}/ndof vs. #phi;#phi_{Track};#chi^{2}/ndof",
					   100, -3.15, 3.15, 100, 0., 10.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_chi2_vs_eta",
					   "#chi^{2} vs. #eta;#eta_{Track};#chi^{2}",
					   100, -3.15, 3.15, 500, 0., 500.));  
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_normchi2_vs_eta",
					   "#chi^{2}/ndof vs. #eta;#eta_{Track};#chi^{2}/ndof",
					   100,-3.15,3.15, 100, 0., 10.));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_kappa_vs_phi",
					   "#kappa vs. #phi;#phi_{Track};#kappa",
					   100,-3.15,3.15, 100, .0,.05));
  vTrack2DHistos_.push_back(tfd.make<TH2F>("h2_kappa_vs_eta",
					   "#kappa vs. #eta;#eta_{Track};#kappa",
					   100,-3.15,3.15, 100, .0,.05));
 

}


void
TrackerOfflineValidation::bookDirHists( TFileDirectory &tfd, const Alignable& ali, const AlignableObjectId &aliobjid)
{
  if(debug_) cout<<"TrackerOfflineValidation::bookDirHists"<<endl;

  std::vector<Alignable*> alivec(ali.components());
  for(int i=0, iEnd = ali.components().size();i < iEnd; ++i) {
    std::string structurename  = aliobjid.typeToName((alivec)[i]->alignableObjectId());
    LogDebug("TrackerOfflineValidation") << "StructureName = " << structurename;
    std::stringstream dirname;
    if(structurename=="Strip") structurename="MyStrip";

    dirname << structurename;
   
    // add no suffix counter to Strip and Pixel, just aesthetics
    if (structurename != "MyStrip" && structurename != "Pixel") dirname << "_" << i+1;
    
    if(debug_)dout<<endl;

    if (structurename.find("Endcap",0) != std::string::npos )   {
      TFileDirectory f = tfd.mkdir((dirname.str()).c_str());
      bookHists(f, *(alivec)[i], ali.alignableObjectId() , i, aliobjid);
      if(debug_)dout<<endl;
      bookDirHists( f, *(alivec)[i], aliobjid);
      if(debug_)dout<<endl;
    } else if( !(this->isDetOrDetUnit( (alivec)[i]->alignableObjectId()) )
	      || alivec[i]->components().size() > 1) {      
      TFileDirectory f = tfd.mkdir((dirname.str()).c_str());
      bookHists(tfd, *(alivec)[i], ali.alignableObjectId() , i, aliobjid);
     if(debug_)dout<<endl;
     bookDirHists( f, *(alivec)[i], aliobjid);
     if(debug_)dout<<endl;
    } else {
      if(debug_)dout<<endl;
      bookHists(tfd, *(alivec)[i], ali.alignableObjectId() , i, aliobjid);
      if(debug_)dout<<endl;
    }
  }

  if(debug_) cout<<"Exit TrackerOfflineValidation::bookDirHists"<<endl;
}




void 
TrackerOfflineValidation::bookHists(TFileDirectory &tfd, const Alignable& ali, align::StructureType type, int i, const AlignableObjectId &aliobjid)
{
  if(debug_) cout<<"TrackerOfflineValidation::bookHists"<<endl;

  TrackerAlignableId aliid;
  const DetId id = ali.id();

  // comparing subdetandlayer to subdetIds gives a warning at compile time
  // -> subdetandlayer could also be pair<uint,uint> but this has to be adapted
  // in AlignableObjId 
  std::pair<int,int> subdetandlayer = aliid.typeAndLayerFromDetId(id);

  align::StructureType subtype = align::invalid;
  
  // are we on or just above det, detunit level respectively?
  if (type == align::AlignableDetUnit )subtype = type;
  else if( this->isDetOrDetUnit(ali.alignableObjectId()) ) subtype = ali.alignableObjectId();
  
  // construct histogramm title and name
  std::stringstream histoname, 
    histotitle, 
    normhistoname, 
    normhistotitle, 
    xprimehistoname, 
    xprimehistotitle, 
    xprimehistoname2[2], 
    xprimehistotitle2[2], 
    normxprimehistoname, 
    normxprimehistotitle,
    yprimehistoname, 
    yprimehistotitle, 
    normyprimehistoname, 
    normyprimehistotitle,
    tahistoname, 
    lahistoname, 
    dthistoname,
    dthistotitle,
    tahistotitle, 
    lahistotitle,
    resxprimeoverthetahistoname, 
    resxprimeoverthetathickhistoname, 
    resxprimeoverthetathinhistoname, 
    resxprimeoverthetahistoname2[2], 
    resxprimevsthetahistoname,
    resxprimeoverthetahistotitle, 
    resxprimeoverthetathinhistotitle, 
    resxprimeoverthetathickhistotitle, 
    resxprimeoverthetahistotitle2[2], 
    resxprimevsthetahistotitle,
    uorientationhistoname, 
    vorientationhistoname, 
    worientationhistoname,
    uorientationhistotitle, 
    vorientationhistotitle, 
    worientationhistotitle,
    hitetahistoname, 
    hitetahistotitle, 
    hitzhistoname, 
    hitzhistotitle,
    localthetahistoname, 
    localthetahistotitle,
    localphihistoname, 
    localphihistotitle,
    uzhistoname,
    vzhistoname,
    wzhistoname,
    uzhistotitle,
    vzhistotitle,
    wzhistotitle,
    duzhistoname[12],
    dwzhistoname[12],
    duzhistotitle[12],
    dwzhistotitle[12],
    dttimehistotitle,
    dttimehistoname,
    dttimeerrhistotitle,
    dttimeerrhistoname,
    hcaltimehistotitle,
    hcaltimehistoname,
    duvsdttimehistotitle,
    duvsdttimehistoname,
    dwvsdttimehistotitle,
    dwvsdttimehistoname,
    ndthistoname,
    ndthistotitle,
    duvstrkanglehistotitle,
    duvsloranglehistotitle,
    duvstrkanglehistoname,
    duvsloranglehistoname,
    chargehistoname,
    chargehistotitle,
    nstripshistoname,
    nstripshistotitle,
    chargevsdttimehistoname,
    chargevsdttimehistotitle;

  std::string wheel_or_layer;

  if( this->isEndCap(static_cast<uint32_t>(subdetandlayer.first)) ) wheel_or_layer = "_wheel_";
  else if ( this->isBarrel(static_cast<uint32_t>(subdetandlayer.first)) ) wheel_or_layer = "_layer_";
  else edm::LogWarning("TrackerOfflineValidation") << "@SUB=TrackerOfflineValidation::bookHists" 
						   << "Unknown subdetid: " <<  subdetandlayer.first;     
   
  for(int ihist=0;ihist<12;ihist++){
    duzhistoname[ihist] << "h_duz_"<<ihist<<"_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
    
    dwzhistoname[ihist] << "h_dwz_"<<ihist<<"_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
    
    duzhistotitle[ihist] << "delta(u) for Module "<<ihist<<" "<< id.rawId();
    dwzhistotitle[ihist] << "delta(w) for Module "<<ihist<<" "<< id.rawId();
  }

  for(int ihist=0;ihist<2;ihist++){
    xprimehistoname2[ihist] << "h_xprime_residuals_"<<ihist<<"_subdet_" << subdetandlayer.first 
			    << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  
    resxprimeoverthetahistoname2[ihist] << "h_xprimeovertheta_"<<ihist<<"_subdet_" << subdetandlayer.first 
					<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

    xprimehistotitle2[ihist] << "X' Residual "<<ihist<<" for module " << id.rawId() 
			     << ";(x_{pred} - x_{rec})' [cm]";

    resxprimeoverthetahistotitle2[ihist] << "Residual/#Delta tan(#theta) "<<ihist
					 <<" for module " << id.rawId() << ";(x_{pred} - x_{rec}) [cm]";

  }
  histoname << "h_residuals_subdet_" << subdetandlayer.first 
	    << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  chargehistoname << "h_charge_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  nstripshistoname << "h_nstrips_subdet_" << subdetandlayer.first 
		   << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  chargevsdttimehistoname << "h_chargevsdttime_subdet_" << subdetandlayer.first 
			  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  duvsdttimehistoname << "h_duvsdttime_subdet_" << subdetandlayer.first 
		      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  dwvsdttimehistoname << "h_dwvsdttime_subdet_" << subdetandlayer.first 
		      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  dttimehistoname << "h_dttime_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  dttimehistoname << "h_dttime_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  dttimeerrhistoname << "h_dttimeerr_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  ndthistoname << "h_ndt_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  hcaltimehistoname << "h_hcaltime_subdet_" << subdetandlayer.first 
		    << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  uzhistoname << "h_uz_subdet_" << subdetandlayer.first 
	      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  vzhistoname << "h_vz_subdet_" << subdetandlayer.first 
	      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  wzhistoname << "h_wz_subdet_" << subdetandlayer.first 
	      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  hitetahistoname << "h_hiteta_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  
  hitzhistoname << "h_hitz_subdet_" << subdetandlayer.first 
		<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  
  xprimehistoname << "h_xprime_residuals_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  yprimehistoname << "h_yprime_residuals_subdet_" << subdetandlayer.first 
		  << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  normhistoname << "h_normresiduals_subdet_" << subdetandlayer.first 
		<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  normxprimehistoname << "h_normxprimeresiduals_subdet_" << subdetandlayer.first 
		      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  normyprimehistoname << "h_normyprimeresiduals_subdet_" << subdetandlayer.first 
		      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  tahistoname << "h_tantrackangle_subdet_" << subdetandlayer.first 
	      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  dthistoname << "h_deltatan_subdet_" << subdetandlayer.first 
	      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();
  
  lahistoname << "h_tanlorentzangle_subdet_" << subdetandlayer.first 
	      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  resxprimeoverthetahistoname << "h_xprimeovertheta_subdet_" << subdetandlayer.first 
			      << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  resxprimeoverthetathickhistoname << "h_xprimeoverthetathick_subdet_" << subdetandlayer.first 
				   << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  resxprimeoverthetathinhistoname  << "h_xprimeoverthetathin_subdet_" << subdetandlayer.first 
				   << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  resxprimevsthetahistoname << "h_xprimevstheta_subdet_" << subdetandlayer.first 
			    << wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  duvstrkanglehistoname << "h_duvstrkangle_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  duvsloranglehistoname << "h_duvslorangle_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  uorientationhistoname << "h_uOrientation_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  vorientationhistoname << "h_vOrientation_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  worientationhistoname << "h_wOrientation_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  localthetahistoname << "h_localtheta_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  localphihistoname << "h_localphi_subdet_" << subdetandlayer.first 
			<< wheel_or_layer << subdetandlayer.second << "_module_" << id.rawId();

  histotitle << "Residual for module " << id.rawId() << ";x_{pred} - x_{rec} [cm]";

  chargehistotitle << "cluster charge " << id.rawId() << " charge";
  chargevsdttimehistotitle << "cluster charge vs. dt time" << id.rawId() << "charge:dttime";
  nstripshistotitle << "cluster nstrips " << id.rawId() << " nstrips";

  normhistotitle << "Normalized Residual for module " << id.rawId() << ";x_{pred} - x_{rec}/#sigma";
  xprimehistotitle << "X' Residual for module " << id.rawId() << ";(x_{pred} - x_{rec})' [cm]";

  normxprimehistotitle << "Normalized X' Residual for module " << id.rawId() << ";(x_{pred} - x_{rec})'/#sigma";
  yprimehistotitle << "Y' Residual for module " << id.rawId() << ";(y_{pred} - y_{rec})' [cm]";
  normyprimehistotitle << "Normalized Y' Residual for module " << id.rawId() << ";(y_{pred} - y_{rec})'/#sigma";
  
  tahistotitle << "tan(#theta_{track}) for module " << id.rawId() << ";tan(#theta_{track})";
  lahistotitle << "tan(#theta_{LA}) for module " << id.rawId() << ";tan(#theta_{LA})";
  dthistotitle << "tan(#theta_{track})-tan(#theta_{LA}) for module " << id.rawId() << ";tan(#theta_{track})-tan(#theta_{LA})";

  resxprimeoverthetahistotitle << "Residual/#Delta tan(#theta) for module " << id.rawId() << ";(x_{pred} - x_{rec})/(tan(#theta_{track})-tan(#theta_{LA}) [cm]";

  resxprimeoverthetathickhistotitle << "Residual/#Delta tan(#theta) (think) for module " << id.rawId() << ";(x_{pred} - x_{rec})/(tan(#theta_{track})-tan(#theta_{LA}) [cm]";
  
  resxprimeoverthetathinhistotitle  << "Residual/#Delta tan(#theta) (thin) for module " << id.rawId() << ";(x_{pred} - x_{rec})/(tan(#theta_{track})-tan(#theta_{LA}) [cm]";

  resxprimevsthetahistotitle << "Residual vs. #Delta tan(#theta) for module " << id.rawId() << ";(x_{pred} - x_{rec}) [cm]";
 

  duvsdttimehistotitle << "#Delta_{u} vs. DT Time " << id.rawId() << "; #Delta_{u} vs. DT Time";
  dwvsdttimehistotitle << "#Delta_{w} vs. DT Time " << id.rawId() << "; #Delta_{w} vs. DT Time";

  uorientationhistotitle << "uOrientation " << id.rawId() << "; uOrientation";
  vorientationhistotitle << "vOrientation " << id.rawId() << "; vOrientation";
  worientationhistotitle << "wOrientation " << id.rawId() << "; wOrientation";
  localthetahistotitle << "localtheta " << id.rawId() << "; tan(#theta)";
  localphihistotitle << "localphi " << id.rawId() << "; cos(#phi)";
  dttimehistotitle << "DT time " << id.rawId() << "; time";
  dttimeerrhistotitle << "DT timeError " << id.rawId() << "; timeError";
  ndthistotitle << "nDT " << id.rawId() << "; nDT";
  hcaltimehistotitle << "HCAL time " << id.rawId() << "; time";

  hitetahistotitle << "Hit Eta " << id.rawId() << "; #eta";
  hitzhistotitle << "Hit Z " << id.rawId() << "; z";
  uzhistotitle << "uOrientation vs. Hit Z "<< id.rawId();
  vzhistotitle << "vOrientation vs. Hit Z "<< id.rawId();
  wzhistotitle << "wOrientation vs. Hit Z "<< id.rawId();

  duvstrkanglehistotitle <<"#Delta_{u} vs. tan(#theta_{trk}) "<< id.rawId();
  duvsloranglehistotitle <<"#Delta_{u} vs. tan(#theta_{LA}) "<< id.rawId();

  if( this->isDetOrDetUnit( subtype ) ) {
    ModuleHistos &histStruct = this->getHistStructFromMap(id);
    int nbins = 0;
    double xmin = 0., xmax = 0.;

    //book only TH2 histos
    if(bookTH2_){
      if(debug_) cout<<"Booking "<<resxprimevsthetahistoname.str()<<" isTransient "<<moduleLevelHistsTransient_<<endl;   
      if(debug_) cout<<"Booking TH2 histos"<<endl;
      this->getBinning(id.subdetId(), XprimevsThetaResidual, nbins, xmin, xmax);
      
      histStruct.ResXprimeVsThetaHisto = this->bookTH2F(moduleLevelHistsTransient_, tfd, 
							resxprimevsthetahistoname.str().c_str(),resxprimevsthetahistotitle.str().c_str(),
							40, -2, 2, 40, -200, 200);
 
      
      histStruct.duvstrkangleHisto = this->bookTH2F(moduleLevelHistsTransient_, tfd, 
						    duvstrkanglehistoname.str().c_str(),duvstrkanglehistotitle.str().c_str(),
						    40, -2, 2, 40, -200, 200);

      if(debug_) cout<<"Booked "<<resxprimevsthetahistoname.str()<<endl;
      
      histStruct.duvsdttimeHisto = this->bookTH2F(moduleLevelHistsTransient_, tfd, 
						  duvsdttimehistoname.str().c_str(),duvsdttimehistotitle.str().c_str(),
						  50, -25, 25, 50, -500, 500);

      histStruct.chargevsdttimeHisto = this->bookTH2F(moduleLevelHistsTransient_, tfd, 
						      chargevsdttimehistoname.str().c_str(),chargevsdttimehistotitle.str().c_str(),
						      50, -25, 25, 50, 0, 1000);

      histStruct.dwvsdttimeHisto = this->bookTH2F(moduleLevelHistsTransient_, tfd, 
						  dwvsdttimehistoname.str().c_str(),dwvsdttimehistotitle.str().c_str(),
						  50, -25, 25, 50, -1000, 1000);
   
    }
    //book all TH1 histos
    if(bookTH1_){

      // decide via cfg if hists in local coordinates should be booked
      /*
      if(lCoorHistOn_) {
	this->getBinning(id.subdetId(), XResidual, nbins, xmin, xmax);
	histStruct.ResHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					     histoname.str().c_str(),histotitle.str().c_str(),		     
					     nbins, xmin, xmax);
	//cout<<"Booking ResHisto "<<histoname.str()<<endl;
	this->getBinning(id.subdetId(), NormXResidual, nbins, xmin, xmax);
	histStruct.NormResHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd,
						 normhistoname.str().c_str(),normhistotitle.str().c_str(),
						 nbins, xmin, xmax);
      } 
      */

      histStruct.tanTrackAngleHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						     tahistoname.str().c_str(),tahistotitle.str().c_str(),
						     100, -5, 5);

      histStruct.deltaTanHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						dthistoname.str().c_str(),dthistotitle.str().c_str(),
						100, -5, 5);

      histStruct.HitEtaHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					      hitetahistoname.str().c_str(),hitetahistotitle.str().c_str(),
					      100, -5, 5);

      histStruct.dttimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					      dttimehistoname.str().c_str(),dttimehistotitle.str().c_str(),
					      200, -100, 100);

      histStruct.dttimeerrHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						 dttimeerrhistoname.str().c_str(),dttimeerrhistotitle.str().c_str(),
						 100, 0, 100);
      
      histStruct.ndtHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					   ndthistoname.str().c_str(),ndthistotitle.str().c_str(),
					   100, 0, 100);

      histStruct.chargeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					      chargehistoname.str().c_str(),chargehistotitle.str().c_str(),
					      100, 0, 1000);

      histStruct.nstripsHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					       nstripshistoname.str().c_str(),nstripshistotitle.str().c_str(),
					       10, 0, 10);

      histStruct.hcaltimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						hcaltimehistoname.str().c_str(),hcaltimehistotitle.str().c_str(),
						200, -100, 100);
      /*
      for(int ihist=0;ihist<12;ihist++){
	histStruct.du_zHisto[ihist] = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						     duzhistoname[ihist].str().c_str(),duzhistotitle[ihist].str().c_str(),
						     50, -500, 500);

	histStruct.dw_zHisto[ihist] = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						     dwzhistoname[ihist].str().c_str(),dwzhistotitle[ihist].str().c_str(),
						     50, -1000, 1000);
      }
      */

      /*
      histStruct.HitXHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					    hitxhistoname.str().c_str(),hitxhistotitle.str().c_str(),
					    1000, -100, 100);

      histStruct.HitYHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					    hityhistoname.str().c_str(),hityhistotitle.str().c_str(),
					    1000, -100, 100);
      */

      histStruct.HitZHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
					    hitzhistoname.str().c_str(),hitzhistotitle.str().c_str(),
					    100, -100, 100);
    
      histStruct.tanLorentzAngleHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						       lahistoname.str().c_str(),lahistotitle.str().c_str(),
						       100, -0.5, 0.5);

      histStruct.ResXprimeOverThetaHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
							  resxprimeoverthetahistoname.str().c_str(),resxprimeoverthetahistotitle.str().c_str(),
							  100, -1000, 1000);


      histStruct.ResXprimeOverThetaHisto_thick = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
								resxprimeoverthetathickhistoname.str().c_str(),resxprimeoverthetathickhistotitle.str().c_str(),
								100, -1000, 1000);

      histStruct.ResXprimeOverThetaHisto_thin  = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
								resxprimeoverthetathinhistoname.str().c_str(),resxprimeoverthetathinhistotitle.str().c_str(),
								100, -1000, 1000);
      
      histStruct.uOrientationHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						    uorientationhistoname.str().c_str(),uorientationhistotitle.str().c_str(),
						    3, -1.5, 1.5);

      histStruct.vOrientationHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						    vorientationhistoname.str().c_str(),vorientationhistotitle.str().c_str(),
						    3, -1.5, 1.5);

      histStruct.wOrientationHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						    worientationhistoname.str().c_str(),worientationhistotitle.str().c_str(),
						    3, -1.5, 1.5);

//       histStruct.vminusdusignHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
// 						    vminusdusignhistoname.str().c_str(),vminusdusignhistotitle.str().c_str(),
// 						    5, -2.5, 2.5);

      histStruct.localthetaHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						  localthetahistoname.str().c_str(),localthetahistotitle.str().c_str(),
						  100, -5, 5);

      histStruct.localphiHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						localphihistoname.str().c_str(),localphihistotitle.str().c_str(),
						100, -1, 1);

      this->getBinning(id.subdetId(), XprimeResidual, nbins, xmin, xmax);
      histStruct.ResXprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						 xprimehistoname.str().c_str(),xprimehistotitle.str().c_str(),
						 nbins, xmin, xmax);

      for(int ihist=0;ihist<2;ihist++){
	histStruct.ResXprimeOverThetaHisto2[ihist] = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
								    resxprimeoverthetahistoname2[ihist].str().c_str(),
								    resxprimeoverthetahistotitle2[ihist].str().c_str(),
								    100, -1000, 1000);

	histStruct.ResXprimeHisto2[ihist] = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
							   xprimehistoname2[ihist].str().c_str(),xprimehistotitle2[ihist].str().c_str(),
							   nbins, xmin, xmax);
      }	
      
      /*
      this->getBinning(id.subdetId(), NormXprimeResidual, nbins, xmin, xmax);
      histStruct.NormResXprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						     normxprimehistoname.str().c_str(),normxprimehistotitle.str().c_str(),
						     nbins, xmin, xmax);
      */
      /*
      if( this->isPixel(subdetandlayer.first) || stripYResiduals_ ) {
	this->getBinning(id.subdetId(), YprimeResidual, nbins, xmin, xmax);
	histStruct.ResYprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd,
						   yprimehistoname.str().c_str(),yprimehistotitle.str().c_str(),
						   nbins, xmin, xmax);
	this->getBinning(id.subdetId(), NormYprimeResidual, nbins, xmin, xmax);
	histStruct.NormResYprimeHisto = this->bookTH1F(moduleLevelHistsTransient_, tfd, 
						       normyprimehistoname.str().c_str(),normyprimehistotitle.str().c_str(),
						       nbins, xmin, xmax);
      }
      */
    }
    if(debug_)cout<<"End of bookHists"<<endl;
  }
  if(debug_)cout<<"Exit bookHists"<<endl;
}

TH1* TrackerOfflineValidation::bookTH1F(bool isTransient, TFileDirectory& tfd, const char* histName, const char* histTitle, 
		int nBinsX, double lowX, double highX)
{
  if(isTransient) {
    vDeleteObjects_.push_back(new TH1F(histName, histTitle, nBinsX, lowX, highX));
    return vDeleteObjects_.back(); // return last element of vector
  }
  else
    return tfd.make<TH1F>(histName, histTitle, nBinsX, lowX, highX);


}


TH2* TrackerOfflineValidation::bookTH2F(bool isTransient, TFileDirectory& tfd, const char* histName, const char* histTitle, 
		int nBinsX, double lowX, double highX,int nBinsY, double lowY, double highY)
{
  if(isTransient) {
    vDeleteObjectsTH2_.push_back(new TH2F(histName, histTitle, nBinsX, lowX, highX, nBinsY, lowY, highY));
    return vDeleteObjectsTH2_.back(); // return last element of vector
  }
  else{
    return tfd.make<TH2F>(histName, histTitle, nBinsX, lowX, highX, nBinsY, lowY, highY );
  }
}


bool TrackerOfflineValidation::isBarrel(uint32_t subDetId)
{
  return (subDetId == StripSubdetector::TIB ||
	  subDetId == StripSubdetector::TOB ||
	  subDetId == PixelSubdetector::PixelBarrel );

}

bool TrackerOfflineValidation::isEndCap(uint32_t subDetId)
{
  return ( subDetId == StripSubdetector::TID ||
	   subDetId == StripSubdetector::TEC ||
	   subDetId == PixelSubdetector::PixelEndcap);
}

bool TrackerOfflineValidation::isPixel(uint32_t subDetId)
{
  return (subDetId == PixelSubdetector::PixelBarrel || subDetId == PixelSubdetector::PixelEndcap);
}


bool TrackerOfflineValidation::isDetOrDetUnit(align::StructureType type)
{
  return ( type == align::AlignableDet || type == align::AlignableDetUnit);
}

void 
TrackerOfflineValidation::getBinning(uint32_t subDetId, 
				     TrackerOfflineValidation::HistogrammType residualType, 
				     int &nBinsX, double &lowerBoundX, double &upperBoundX)
{
  // determine if 
  const bool isPixel = this->isPixel(subDetId);
  
  edm::ParameterSet binningPSet;
  
  switch(residualType) 
    {
    case XResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1XResPixelModules");                
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1XResStripModules");                
      break;
    case NormXResidual : 
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1NormXResPixelModules");             
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1NormXResStripModules");                
      break;
    case XprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1XprimeResPixelModules");                
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1XprimeResStripModules");                
      break;
    case NormXprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1NormXprimeResPixelModules");
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1NormXprimeResStripModules");
      break;
    case YprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1YResPixelModules");                
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1YResStripModules");                
      break; 
    case NormYprimeResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH1NormYResPixelModules");             
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH1NormYResStripModules");  
      break;
    case XprimevsThetaResidual :
      if(isPixel) binningPSet = parset_.getParameter<edm::ParameterSet>("TH2XprimevsThetaPixelModules");             
      else binningPSet        = parset_.getParameter<edm::ParameterSet>("TH2XprimevsThetaStripModules");  
      break;
    }
  nBinsX      = binningPSet.getParameter<int32_t>("Nbinx");		       
  lowerBoundX = binningPSet.getParameter<double>("xmin");		       
  upperBoundX = binningPSet.getParameter<double>("xmax");     
  
}

void 
TrackerOfflineValidation::setSummaryBin(int bin, TH1* targetHist, TH1* sourceHist)
{
  if(targetHist && sourceHist) {
    targetHist->SetBinContent(bin, sourceHist->GetMean(1));
    if(useFwhm_) targetHist->SetBinError(bin, Fwhm(sourceHist)/2.);
    else targetHist->SetBinError(bin, sourceHist->GetRMS(1) );
  } else {
    return;
  }

}


void 
TrackerOfflineValidation::summarizeBinInContainer( int bin, SummaryContainer &targetContainer, 
						   SummaryContainer &sourceContainer)
{
  
  
  this->setSummaryBin(bin, targetContainer.summaryXResiduals_, sourceContainer.sumXResiduals_);
  //this->setSummaryBin(bin, targetContainer.summaryNormXResiduals_, sourceContainer.sumNormXResiduals_);
  // If no y-residual hists, just returns:
  //this->setSummaryBin(bin, targetContainer.summaryYResiduals_, sourceContainer.sumYResiduals_);
  //this->setSummaryBin(bin, targetContainer.summaryNormYResiduals_, sourceContainer.sumNormYResiduals_);

}

void 
TrackerOfflineValidation::summarizeBinInContainer( int bin, uint32_t subDetId, 
						   SummaryContainer &targetContainer, 
						   ModuleHistos &sourceContainer)
{

  // takes two summary Containers and sets summaryBins for all histogramms
  this->setSummaryBin(bin, targetContainer.summaryXResiduals_, sourceContainer.ResXprimeHisto);
  /*
  this->setSummaryBin(bin, targetContainer.summaryNormXResiduals_, sourceContainer.NormResXprimeHisto);
  if( this->isPixel(subDetId) || stripYResiduals_ ) {
    this->setSummaryBin(bin, targetContainer.summaryYResiduals_, sourceContainer.ResYprimeHisto);
    this->setSummaryBin(bin, targetContainer.summaryNormYResiduals_, sourceContainer.NormResYprimeHisto);
  }
  */
}




TrackerOfflineValidation::ModuleHistos& 
TrackerOfflineValidation::getHistStructFromMap(const DetId& detid)
{

  // get a struct with histogramms from the respective map
  // if no object exist, the reference is automatically created by the map
  // throw exception if non-tracker id is passed
  uint subdetid = detid.subdetId();
  if(subdetid == PixelSubdetector::PixelBarrel ) {
    return mPxbResiduals_[detid.rawId()];
  } else if (subdetid == PixelSubdetector::PixelEndcap) {
    return mPxeResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TIB) {
    return mTibResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TID) {
    return mTidResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TOB) {
    return mTobResiduals_[detid.rawId()];
  } else if(subdetid  == StripSubdetector::TEC) {
    return mTecResiduals_[detid.rawId()];
  } else {
    throw cms::Exception("Geometry Error")
      << "[TrackerOfflineValidation] Error, tried to get reference for non-tracker subdet " << subdetid 
      << " from detector " << detid.det();
    return mPxbResiduals_[0];
  }
  
}


// ------------ method called to for each event  ------------
void
TrackerOfflineValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(debug_) cout<<"TrackerOfflineValidation::analyze"<<endl;
  int evt_run                       = iEvent.id().run()        ;
  int evt_lumiBlock                 = iEvent.luminosityBlock() ;  
  int evt_event                     = iEvent.id().event()      ;
  
  ofile<<evt_run<<" "<<evt_lumiBlock<<" "<<evt_event<<endl;


  if (useOverflowForRMS_)TH1::StatOverflows(kTRUE);
  this->checkBookHists(iSetup); // check whether hists are are booked and do so if not yet done
  
  if(debug_)dout<<endl;
  //using namespace edm;
  TrackerValidationVariables avalidator_(iSetup,parset_);
  edm::Service<TFileService> fs;
    
  std::vector<TrackerValidationVariables::AVTrackStruct> vTrackstruct;
  avalidator_.fillTrackQuantities(iEvent,vTrackstruct);
  if(debug_)dout<<endl;
  std::vector<TrackerValidationVariables::AVHitStruct> v_hitstruct;
  //avalidator_.fillHitQuantities(iEvent,v_hitstruct);
  avalidator_.fillHitQuantities(iEvent,iSetup,v_hitstruct);
  if(debug_)dout<<endl;
  
  for (std::vector<TrackerValidationVariables::AVTrackStruct>::const_iterator it = vTrackstruct.begin(),
	 itEnd = vTrackstruct.end(); it != itEnd; ++it) {
    
    // Fill 1D track histos
    static const int etaindex = this->GetIndex(vTrackHistos_,"h_tracketa");
    vTrackHistos_[etaindex]->Fill(it->eta);
    static const int kappaindex = this->GetIndex(vTrackHistos_,"h_curvature");
    vTrackHistos_[kappaindex]->Fill(it->kappa);
    static const int kappaposindex = this->GetIndex(vTrackHistos_,"h_curvature_pos");
    if(it->charge > 0)
      vTrackHistos_[kappaposindex]->Fill(fabs(it->kappa));
    static const int kappanegindex = this->GetIndex(vTrackHistos_,"h_curvature_neg");
    if(it->charge < 0)
      vTrackHistos_[kappanegindex]->Fill(fabs(it->kappa));
    static const int normchi2index = this->GetIndex(vTrackHistos_,"h_normchi2");
    vTrackHistos_[normchi2index]->Fill(it->normchi2);
    static const int chi2index = this->GetIndex(vTrackHistos_,"h_chi2");
    vTrackHistos_[chi2index]->Fill(it->chi2);
    static const int ptindex = this->GetIndex(vTrackHistos_,"h_pt");
    vTrackHistos_[ptindex]->Fill(it->pt);
    if(it->ptError != 0.) {
      static const int ptResolutionindex = this->GetIndex(vTrackHistos_,"h_ptResolution");
      vTrackHistos_[ptResolutionindex]->Fill(it->ptError/it->pt);
    }
    // Fill track profiles
    static const int d0phiindex = this->GetIndex(vTrackProfiles_,"p_d0_vs_phi");
    vTrackProfiles_[d0phiindex]->Fill(it->phi,it->d0);
    static const int dzphiindex = this->GetIndex(vTrackProfiles_,"p_dz_vs_phi");
    vTrackProfiles_[dzphiindex]->Fill(it->phi,it->dz);
    static const int d0etaindex = this->GetIndex(vTrackProfiles_,"p_d0_vs_eta");
    vTrackProfiles_[d0etaindex]->Fill(it->eta,it->d0);
    static const int dzetaindex = this->GetIndex(vTrackProfiles_,"p_dz_vs_eta");
    vTrackProfiles_[dzetaindex]->Fill(it->eta,it->dz);
    static const int chiphiindex = this->GetIndex(vTrackProfiles_,"p_chi2_vs_phi");
    vTrackProfiles_[chiphiindex]->Fill(it->phi,it->chi2);
    static const int normchiphiindex = this->GetIndex(vTrackProfiles_,"p_normchi2_vs_phi");
    vTrackProfiles_[normchiphiindex]->Fill(it->phi,it->normchi2);
    static const int chietaindex = this->GetIndex(vTrackProfiles_,"p_chi2_vs_eta");
    vTrackProfiles_[chietaindex]->Fill(it->eta,it->chi2);
    static const int normchietaindex = this->GetIndex(vTrackProfiles_,"p_normchi2_vs_eta");
    vTrackProfiles_[normchietaindex]->Fill(it->eta,it->normchi2);
    static const int kappaphiindex = this->GetIndex(vTrackProfiles_,"p_kappa_vs_phi");
    vTrackProfiles_[kappaphiindex]->Fill(it->phi,it->kappa);
    static const int kappaetaindex = this->GetIndex(vTrackProfiles_,"p_kappa_vs_eta");
    vTrackProfiles_[kappaetaindex]->Fill(it->eta,it->kappa);
    static const int ptResphiindex = this->GetIndex(vTrackProfiles_,"p_ptResolution_vs_phi");
    vTrackProfiles_[ptResphiindex]->Fill(it->phi,it->ptError/it->pt);
    static const int ptResetaindex = this->GetIndex(vTrackProfiles_,"p_ptResolution_vs_eta");
    vTrackProfiles_[ptResetaindex]->Fill(it->eta,it->ptError/it->pt);

    // Fill 2D track histos
    static const int d0phiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_d0_vs_phi");
    vTrack2DHistos_[d0phiindex_2d]->Fill(it->phi,it->d0);
    static const int dzphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_dz_vs_phi");
    vTrack2DHistos_[dzphiindex_2d]->Fill(it->phi,it->dz);
    static const int d0etaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_d0_vs_eta");
    vTrack2DHistos_[d0etaindex_2d]->Fill(it->eta,it->d0);
    static const int dzetaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_dz_vs_eta");
    vTrack2DHistos_[dzetaindex_2d]->Fill(it->eta,it->dz);
    static const int chiphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_chi2_vs_phi");
    vTrack2DHistos_[chiphiindex_2d]->Fill(it->phi,it->chi2);
    static const int normchiphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_normchi2_vs_phi");
    vTrack2DHistos_[normchiphiindex_2d]->Fill(it->phi,it->normchi2);
    static const int chietaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_chi2_vs_eta");
    vTrack2DHistos_[chietaindex_2d]->Fill(it->eta,it->chi2);
    static const int normchietaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_normchi2_vs_eta");
    vTrack2DHistos_[normchietaindex_2d]->Fill(it->eta,it->normchi2);
    static const int kappaphiindex_2d = this->GetIndex(vTrack2DHistos_,"h2_kappa_vs_phi");
    vTrack2DHistos_[kappaphiindex_2d]->Fill(it->phi,it->kappa);
    static const int kappaetaindex_2d = this->GetIndex(vTrack2DHistos_,"h2_kappa_vs_eta");
    vTrack2DHistos_[kappaetaindex_2d]->Fill(it->eta,it->kappa);
     
  } // finish loop over track quantities
  if(debug_)cout<<"Finished loop over track quantities"<<endl;
  
  // hit quantities: residuals, normalized residuals
  for (std::vector<TrackerValidationVariables::AVHitStruct>::const_iterator it = v_hitstruct.begin(),
  	 itEnd = v_hitstruct.end(); it != itEnd; ++it) {

    DetId detid(it->rawDetId);
    ModuleHistos &histStruct = this->getHistStructFromMap(detid);
    
    if(it->nvalidmu==1 && bookTH1_){
      if(it->dttimeerr != -999)   histStruct.dttimeerrHisto->Fill(it->dttimeerr);
      if(it->ndt != -999)         histStruct.ndtHisto->Fill(it->ndt);
    }
    
    if(debug_)cout<<"nvalidmu "<<it->nvalidmu<<" ndt "<<it->ndt<<" dttimeerr "<<it->dttimeerr<<endl;
    
    if(it->nvalidmu==1 && it->ndt>=25 &&  it->dttimeerr<10){
      
      if(bookTH2_){
	if(debug_){
	  cout<<"Fill TH2 "<<it->resXprime<<" tanTrackAngle "<<it->tanTrackAngle<<" tanLorentzAngle "<<it->tanLorentzAngle<<endl;
	  if(it->resXprime == -999.) cout<<"Invalid residual deltaTanTheta "<<it->tanTrackAngle-it->tanLorentzAngle<<endl;
	}
	
	if(it->resXprime != -999. && it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.) {
	  if(fabs(it->tanTrackAngle-it->tanLorentzAngle>0)){
	    histStruct.ResXprimeVsThetaHisto -> Fill(it->tanTrackAngle-it->tanLorentzAngle, 10000*it->resXprime);
	    histStruct.duvstrkangleHisto     -> Fill(it->tanTrackAngle,                     10000*it->resXprime);
	  }
	  if(fabs(it->tanTrackAngle-it->tanLorentzAngle<0)){
	    histStruct.ResXprimeVsThetaHisto -> Fill(it->tanTrackAngle-it->tanLorentzAngle, -10000*it->resXprime);
	    histStruct.duvstrkangleHisto     -> Fill(it->tanTrackAngle,                     -10000*it->resXprime);
	  }
	}
	
	if(it->charge != -999. && it->dttime != -999.){
	  histStruct.chargevsdttimeHisto->Fill(it->dttime,it->charge);
	}

	if(it->resXprime != -999. && it->dttime != -999.){
	  histStruct.duvsdttimeHisto->Fill(it->dttime,10000*it->resXprime);

	  if(it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.){
	    float deltaTanTheta=it->tanTrackAngle-it->tanLorentzAngle;
	    if(deltaTanTheta != 0) histStruct.dwvsdttimeHisto->Fill(it->dttime,10000*it->resXprime/fabs(deltaTanTheta));
	    
	  }
	}
	
	if(debug_) cout<<"Filled TH2"<<endl;
      }
      
      //fill TH1 histos
      if(bookTH1_){
	// fill histos in local coordinates if set in cf
	/*
	  if(lCoorHistOn_) {
	  histStruct.ResHisto->Fill(it->resX);
	  if(it->resErrX != 0) histStruct.NormResHisto->Fill(it->resX/it->resErrX);
	  }
	*/
	//if(it->x != -999.) histStruct.HitXHisto->Fill(it->x);
	//if(it->y != -999.) histStruct.HitYHisto->Fill(it->y);

	if(it->dttime != -999)         histStruct.dttimeHisto->  Fill(it->dttime);
	if(it->hcaltime != -999)       histStruct.hcaltimeHisto->Fill(it->hcaltime);
	if(it->z != -999.)             histStruct.HitZHisto->    Fill(it->z);
	if(it->eta != -999.)           histStruct.HitEtaHisto->  Fill(it->eta);
	if(it->charge != -999.)        histStruct.chargeHisto->  Fill(it->charge);
	if(it->nstrips != -999)        histStruct.nstripsHisto-> Fill(it->nstrips);
	
	if(it->uOrientation != -999.)   histStruct.uOrientationHisto->Fill(it->uOrientation);
	if(it->vOrientation != -999.)   histStruct.vOrientationHisto->Fill(it->vOrientation);
	//if(it->vOrientation != -999. && it->dusign != -999.) histStruct.vminusdusignHisto->Fill(it->vOrientation-it->dusign);
	if(it->wOrientation != -999.)   histStruct.wOrientationHisto->Fill(it->wOrientation);
	if(it->localtheta != -999.)     histStruct.localthetaHisto->Fill(tan(it->localtheta));
	if(it->localphi != -999.)       histStruct.localphiHisto->Fill(cos(it->localphi));
	
	/*
	if(it->z != -999. && it->resXprime!=-999. && it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.){
	  float dw=it->resXprime/fabs(it->tanTrackAngle-it->tanLorentzAngle);
	  histStruct.du_zHisto[this->getHisto(it->z)]->Fill(10000*it->resXprime);
	  histStruct.dw_zHisto[this->getHisto(it->z)]->Fill(10000*dw);
	}
	*/

	if(it->resXprime != -999.) {
	  histStruct.ResXprimeHisto->Fill(10000*it->resXprime);
	  
	  if(it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.){
	    float deltaTanTheta=it->tanTrackAngle-it->tanLorentzAngle;
	    if(deltaTanTheta != 0) histStruct.ResXprimeOverThetaHisto->Fill(10000*it->resXprime/fabs(deltaTanTheta));
	    if(it->ectype == 1)    histStruct.ResXprimeOverThetaHisto_thick -> Fill(10000*it->resXprime/fabs(deltaTanTheta)); 
	    if(it->ectype == 2)    histStruct.ResXprimeOverThetaHisto_thin  -> Fill(10000*it->resXprime/fabs(deltaTanTheta)); 
	    
	    int ih=0;
	    if(deltaTanTheta<0)ih=1;
	    histStruct.ResXprimeHisto2[ih]->Fill(10000*it->resXprime);
	    if(deltaTanTheta != 0) histStruct.ResXprimeOverThetaHisto2[ih]->Fill(10000*it->resXprime/fabs(deltaTanTheta));
	    
	  }
	  /*
	    if(it->resXprimeErr != 0 && it->resXprimeErr != -999 ) {	
	    histStruct.NormResXprimeHisto->Fill(it->resXprime/it->resXprimeErr);
	    } 
	  */
	}
	/*
	  if(it->resYprime != -999.) {
	  if( this->isPixel(detid.subdetId())  || stripYResiduals_ ) {
	  histStruct.ResYprimeHisto->Fill(it->resYprime);
	  if(it->resYprimeErr != 0 && it->resYprimeErr != -999. ) {	
	  histStruct.NormResYprimeHisto->Fill(it->resYprime/it->resYprimeErr);
	  } 
	  }
	  }
	*/
	if(it->tanTrackAngle != -999.)	                                  histStruct.tanTrackAngleHisto->Fill(it->tanTrackAngle);
	if(it->tanLorentzAngle != -999.)	                          histStruct.tanLorentzAngleHisto->Fill(it->tanLorentzAngle);
      	if(it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.)	  histStruct.deltaTanHisto->Fill(it->tanTrackAngle-it->tanLorentzAngle);
	

	if(overlappOn_) {
	  std::pair<uint32_t,uint32_t> tmp_pair(std::make_pair(it->rawDetId, it->overlapres.first));
	  if(it->overlapres.first != 0 ) {
	    if( hOverlappResidual[tmp_pair] ) {
	      hOverlappResidual[tmp_pair]->Fill(it->overlapres.second);
	    } else if( hOverlappResidual[std::make_pair( it->overlapres.first, it->rawDetId) ]) {
	      hOverlappResidual[std::make_pair( it->overlapres.first, it->rawDetId) ]->Fill(it->overlapres.second);
	    } else {
	      TFileDirectory tfd = fs->mkdir("OverlappResiduals");
	      hOverlappResidual[tmp_pair] = tfd.make<TH1F>(Form("hOverlappResidual_%d_%d",tmp_pair.first,tmp_pair.second),
							   "Overlapp Residuals",100,-50,50);
	      hOverlappResidual[tmp_pair]->Fill(it->overlapres.second);
	    }
	  }
	} // end overlappOn
      }
    }
  }
  if(debug_) cout<<"End TrackerOfflineValidation::analyze"<<endl;
  if (useOverflowForRMS_) TH1::StatOverflows(kFALSE);  
}



// ------------ method called once each job just after ending the event loop  ------------
void 
TrackerOfflineValidation::endJob()
{
  if(debug_) cout<<"TrackerOfflineValidation::endJob"<<endl;
  ofile.close();

  AlignableTracker aliTracker(&(*tkGeom_));
  edm::Service<TFileService> fs;   
  AlignableObjectId aliobjid;

  if(debug_) cout<<"make TTree"<<endl;
  TTree *tree = fs->make<TTree>("TkOffVal","TkOffVal");
  TkOffTreeVariables *treeMemPtr = new TkOffTreeVariables;
  // We create branches for all members of 'TkOffTreeVariables' (even if not needed).
  // This works because we have a dictionary for 'TkOffTreeVariables'
  // (see src/classes_def.xml and src/classes.h):
  tree->Branch("TkOffTreeVariables", &treeMemPtr); // address of pointer!

  if(debug_) cout<<"fill TTree"<<endl;
  if(bookTH1_ && fillTree_){
    this->fillTree(*tree, mPxbResiduals_, *treeMemPtr, *tkGeom_);
    this->fillTree(*tree, mPxeResiduals_, *treeMemPtr, *tkGeom_);
    this->fillTree(*tree, mTibResiduals_, *treeMemPtr, *tkGeom_);
    this->fillTree(*tree, mTidResiduals_, *treeMemPtr, *tkGeom_);
    this->fillTree(*tree, mTobResiduals_, *treeMemPtr, *tkGeom_);
    this->fillTree(*tree, mTecResiduals_, *treeMemPtr, *tkGeom_);
  }
  delete treeMemPtr; treeMemPtr = 0;

  if(debug_) cout<<"make vTrackHistos"<<endl;
  static const int kappadiffindex = this->GetIndex(vTrackHistos_,"h_diff_curvature");
  vTrackHistos_[kappadiffindex]->Add(vTrackHistos_[this->GetIndex(vTrackHistos_,"h_curvature_neg")],
				     vTrackHistos_[this->GetIndex(vTrackHistos_,"h_curvature_pos")],-1,1);

  // Collate Information for Subdetectors
  // create summary histogramms recursively

  std::vector<TrackerOfflineValidation::SummaryContainer > vTrackerprofiles;
  if(debug_) cout<<"call collateSummaryHists"<<endl;
  this->collateSummaryHists((*fs),(aliTracker), 0, aliobjid, vTrackerprofiles);
   
}


void
TrackerOfflineValidation::collateSummaryHists( TFileDirectory &tfd, const Alignable& ali, int i, 
					       const AlignableObjectId &aliobjid, 
					       std::vector< TrackerOfflineValidation::SummaryContainer > &v_levelProfiles)
{
  
  if(debug_) cout<<"TrackerOfflineValidation::collateSummaryHists"<<endl;
  std::vector<Alignable*> alivec(ali.components());
  if( this->isDetOrDetUnit((alivec)[0]->alignableObjectId()) ) return;

  for(int iComp=0, iCompEnd = ali.components().size();iComp < iCompEnd; ++iComp) {
    std::vector< TrackerOfflineValidation::SummaryContainer > v_profiles;        
    std::string structurename  = aliobjid.typeToName((alivec)[iComp]->alignableObjectId());
 
    LogDebug("TrackerOfflineValidation") << "StructureName = " << structurename;
    std::stringstream dirname;
    if(structurename=="Strip") structurename="MyStrip";   
//     else{
//       if(removePixel_ && structurename=="Pixel"){
// 	cout<<"TrackerOfflineValidation::collateSummaryHists skip pixel"<<endl;
// 	continue;
//       }
//     }
    dirname << structurename;
    
    // add no suffix counter to strip and pixel -> just aesthetics
    //if (structurename != "Strip" && structurename != "Pixel") dirname << "_" << iComp+1;
    if (structurename != "MyStrip" && structurename != "Pixel") dirname << "_" << iComp+1;
    //cout<<"TrackerOfflineValidation::collateSummaryHists dirname "<<dirname.str()<<endl;    

    if(  !(this->isDetOrDetUnit( (alivec)[iComp]->alignableObjectId()) )
	 || (alivec)[0]->components().size() > 1 ) {
      TFileDirectory f = tfd.mkdir((dirname.str()).c_str());
      this->collateSummaryHists( f, *(alivec)[iComp], i, aliobjid, v_profiles);
      v_levelProfiles.push_back(this->bookSummaryHists(tfd, *(alivec[iComp]), ali.alignableObjectId(), iComp, aliobjid));
      //TH1 *hY = v_levelProfiles[iComp].sumYResiduals_;
      //TH1 *hNormY = v_levelProfiles[iComp].sumNormYResiduals_;
      for(uint n = 0; n < v_profiles.size(); ++n) {
	this->summarizeBinInContainer(n+1, v_levelProfiles[iComp], v_profiles[n] );
	
	if(bookTH2_){
	  if(debug_) cout<<"Add sumResXprimeVsTheta to v_levelProfiles"<<endl;
	  v_levelProfiles[iComp].sumResXprimeVsTheta_->Add(v_profiles[n].sumResXprimeVsTheta_);
	  v_levelProfiles[iComp].sumduvstrkangle_->Add(v_profiles[n].sumduvstrkangle_);
	  //v_levelProfiles[iComp].sumduvslorangle_->Add(v_profiles[n].sumduvslorangle_);
	  //v_levelProfiles[iComp].sumuz_->Add(v_profiles[n].sumuz_);
	  //v_levelProfiles[iComp].sumvz_->Add(v_profiles[n].sumvz_);
	  //v_levelProfiles[iComp].sumwz_->Add(v_profiles[n].sumwz_);
	  v_levelProfiles[iComp].sumduvsdttime_->Add(v_profiles[n].sumduvsdttime_);
	  v_levelProfiles[iComp].sumdwvsdttime_->Add(v_profiles[n].sumdwvsdttime_);
	  v_levelProfiles[iComp].sumchargevsdttime_->Add(v_profiles[n].sumchargevsdttime_);
	}
	if(bookTH1_){
	  v_levelProfiles[iComp].sumXResiduals_->Add(v_profiles[n].sumXResiduals_);
	  for(int ihist=0;ihist<2;ihist++){
	    v_levelProfiles[iComp].sumXResiduals2_[ihist]->Add(v_profiles[n].sumXResiduals2_[ihist]);
	    v_levelProfiles[iComp].sumResXprimeOverTheta2_[ihist]->Add(v_profiles[n].sumResXprimeOverTheta2_[ihist]);
	  }
	  v_levelProfiles[iComp].sumTanTrackAngle_->Add(v_profiles[n].sumTanTrackAngle_);
	  v_levelProfiles[iComp].sumTanLorentzAngle_->Add(v_profiles[n].sumTanLorentzAngle_);
	  v_levelProfiles[iComp].sumdeltaTan_->Add(v_profiles[n].sumdeltaTan_);
	  v_levelProfiles[iComp].sumResXprimeOverTheta_->Add(v_profiles[n].sumResXprimeOverTheta_);
	  v_levelProfiles[iComp].sumResXprimeOverTheta_thick_->Add(v_profiles[n].sumResXprimeOverTheta_thick_);
	  v_levelProfiles[iComp].sumResXprimeOverTheta_thin_ ->Add(v_profiles[n].sumResXprimeOverTheta_thin_);
	  v_levelProfiles[iComp].sumHitEta_->Add(v_profiles[n].sumHitEta_);
	  //v_levelProfiles[iComp].sumHitX_->Add(v_profiles[n].sumHitX_);
	  //v_levelProfiles[iComp].sumHitY_->Add(v_profiles[n].sumHitY_);
	  v_levelProfiles[iComp].sumHitZ_->Add(v_profiles[n].sumHitZ_);
	  //v_levelProfiles[iComp].sumNormXResiduals_->Add(v_profiles[n].sumNormXResiduals_);
	  v_levelProfiles[iComp].sumuOrientation_->Add(v_profiles[n].sumuOrientation_);
	  v_levelProfiles[iComp].sumvOrientation_->Add(v_profiles[n].sumvOrientation_);
	  v_levelProfiles[iComp].sumwOrientation_->Add(v_profiles[n].sumwOrientation_);
	  //v_levelProfiles[iComp].sumvminusdusign_->Add(v_profiles[n].sumvminusdusign_);
	  v_levelProfiles[iComp].sumlocaltheta_->Add(v_profiles[n].sumlocaltheta_);
	  v_levelProfiles[iComp].sumlocalphi_->Add(v_profiles[n].sumlocalphi_);
	  v_levelProfiles[iComp].sumdttime_->Add(v_profiles[n].sumdttime_);
	  v_levelProfiles[iComp].sumdttimeerr_->Add(v_profiles[n].sumdttimeerr_);
	  v_levelProfiles[iComp].sumndt_->Add(v_profiles[n].sumndt_);
	  v_levelProfiles[iComp].sumhcaltime_->Add(v_profiles[n].sumhcaltime_);
	  v_levelProfiles[iComp].sumcharge_->Add(v_profiles[n].sumcharge_);
	  v_levelProfiles[iComp].sumnstrips_->Add(v_profiles[n].sumnstrips_);
	  /*
	  for(int ihist=0;ihist<12;ihist++){
	    v_levelProfiles[iComp].sumdu_z_[ihist]->Add(v_profiles[n].sumdu_z_[ihist]);
	    v_levelProfiles[iComp].sumdw_z_[ihist]->Add(v_profiles[n].sumdw_z_[ihist]);
	  }
	  */
	  //if (hY)     hY->Add(v_profiles[n].sumYResiduals_);         // only if existing
	  //if (hNormY) hNormY->Add(v_profiles[n].sumNormYResiduals_); // dito (pxl, stripYResiduals_)
	}
      }
      //add fit values to stat box
      this->fitResiduals(v_levelProfiles[iComp].sumXResiduals_);
      //this->fitResiduals(v_levelProfiles[iComp].sumNormXResiduals_);
      //if (hY)     this->fitResiduals(hY);     // only if existing (pixel or stripYResiduals_)
      //if (hNormY) this->fitResiduals(hNormY); // dito
    } else {
      // nothing to be done for det or detunits
      continue;
    }
  }
}

TrackerOfflineValidation::SummaryContainer 
TrackerOfflineValidation::bookSummaryHists(TFileDirectory &tfd, const Alignable& ali, 
					   align::StructureType type, int i, 
					   const AlignableObjectId &aliobjid)
{

  if(debug_) cout<<"Start TrackerOfflineValidation::bookSummaryHists"<<endl;
  const uint aliSize = ali.components().size();
  const align::StructureType alitype = ali.alignableObjectId();
  const align::StructureType subtype = ali.components()[0]->alignableObjectId();
  const char *aliTypeName = aliobjid.typeToName(alitype).c_str(); // lifetime of char* OK
  const char *aliSubtypeName = aliobjid.typeToName(subtype).c_str();
  const char *typeName = aliobjid.typeToName(type).c_str();

  const DetId aliDetId = ali.id(); 
  // y residuals only if pixel or specially requested for strip:
  //const bool bookResidY = this->isPixel(aliDetId.subdetId()) || stripYResiduals_;

  SummaryContainer sumContainer;
  
  // Book summary hists with one bin per component, 
  // but special case for Det with two DetUnit that we want to summarize one level up 
  // (e.g. in TOBRods with 12 bins for 6 stereo and 6 rphi DetUnit.)
  //    component of ali is not Det or Det with just one components

  
  const uint subcompSize = ali.components()[0]->components().size();
  if (subtype != align::AlignableDet || subcompSize == 1) { // Det with 1 comp. should not exist anymore...
    const TString title(Form("Summary for substructures in %s %d;%s;",aliTypeName,i,aliSubtypeName));
  
    if(debug_) cout<<"Book summary histos"<<endl;
    sumContainer.summaryXResiduals_ = tfd.make<TH1F>(Form("h_summaryX%s_%d",aliTypeName,i), 
						     title + "#LT #Delta x' #GT",
						     aliSize, 0.5, aliSize+0.5);

    /*
    sumContainer.summaryNormXResiduals_ = tfd.make<TH1F>(Form("h_summaryNormX%s_%d",aliTypeName,i), 
							 title + "#LT #Delta x'/#sigma #GT",
							 aliSize,0.5,aliSize+0.5);
    if (bookResidY) {
      sumContainer.summaryYResiduals_ = tfd.make<TH1F>(Form("h_summaryY%s_%d",aliTypeName,i), 
						       title + "#LT #Delta y' #GT",
						       aliSize, 0.5, aliSize+0.5);
      sumContainer.summaryNormYResiduals_ = tfd.make<TH1F>(Form("h_summaryNormY%s_%d",aliTypeName,i), 
							   title + "#LT #Delta y'/#sigma #GT",
							   aliSize,0.5,aliSize+0.5);
    }
    */
  } else if (subtype == align::AlignableDet && subcompSize > 1) { // fixed: was aliSize before
    if (subcompSize != 2) { // strange... expect only 2 DetUnits in DS layers
      // this 2 is hardcoded factor 2 in binning below and also assummed later on
      edm::LogError("Alignment") << "@SUB=bookSummaryHists"
				 << "Det with " << subcompSize << " components";
    }
    // title contains x-title
    const TString title(Form("Summary for substructures in %s %d;%s;", aliTypeName, i,
			     aliobjid.typeToName(ali.components()[0]->components()[0]->alignableObjectId()).c_str()));
    sumContainer.summaryXResiduals_ 
      = tfd.make<TH1F>(Form("h_summaryX%s_%d", aliTypeName, i), 
		       title + "#LT #Delta x' #GT", (2*aliSize), 0.5, 2*aliSize+0.5);
    /*
    sumContainer.summaryNormXResiduals_ 
      = tfd.make<TH1F>(Form("h_summaryNormX%s_%d", aliTypeName, i), 
		       title + "#LT #Delta x'/#sigma #GT", (2*aliSize), 0.5, 2*aliSize+0.5);

    if (bookResidY) {
      sumContainer.summaryYResiduals_ 
	= tfd.make<TH1F>(Form("h_summaryY%s_%d", aliTypeName, i), 
			 title + "#LT #Delta y' #GT", (2*aliSize), 0.5, 2*aliSize+0.5);
      sumContainer.summaryNormYResiduals_ 
	= tfd.make<TH1F>(Form("h_summaryNormY%s_%d", aliTypeName, i), 
			 title + "#LT #Delta y'/#sigma #GT", (2*aliSize), 0.5, 2*aliSize+0.5);
    }
    */
  } else {
    edm::LogError("TrackerOfflineValidation") << "@SUB=TrackerOfflineValidation::bookSummaryHists" 
					      << "No summary histogramm for hierarchy level " 
					      << aliTypeName << " in subdet " << aliDetId.subdetId();
  }
  


  // Now book hists that just sum up the residual histograms from lower levels.
  // Axis title is copied from lowest level module of structure.
  // Should be safe that y-hists are only touched if non-null pointers...
  int nbins = 0;
  double xmin = 0., xmax = 0.;
  const TString sumTitle(Form("Residual for %s %d in %s;", aliTypeName, i, typeName));
  const TString taTitle(Form("tan(#theta_{track}) for %s %d in %s;", aliTypeName, i, typeName));
  const TString laTitle(Form("tan(#theta_{LA}) for %s %d in %s;", aliTypeName, i, typeName));
  const TString dtTitle(Form("tan(#theta_{track})-tan(#theta_{LA}) for %s %d in %s;", aliTypeName, i, typeName));
  const TString uoTitle(Form("uOrientation for %s %d in %s;", aliTypeName, i, typeName));
  const TString voTitle(Form("vOrientation for %s %d in %s;", aliTypeName, i, typeName));
  const TString woTitle(Form("wOrientation for %s %d in %s;", aliTypeName, i, typeName));
  const TString uzTitle(Form("uOrientation vs. Hit Z for %s %d in %s;", aliTypeName, i, typeName));
  const TString vzTitle(Form("vOrientation vs. Hit Z for %s %d in %s;", aliTypeName, i, typeName));
  const TString wzTitle(Form("wOrientation vs. Hit Z for %s %d in %s;", aliTypeName, i, typeName));
  //const TString vdTitle(Form("vOrientation-dusign for %s %d in %s;", aliTypeName, i, typeName));
  const TString ltTitle(Form("tan(localtheta) for %s %d in %s;", aliTypeName, i, typeName));
  const TString lpTitle(Form("cos(localphi) for %s %d in %s;", aliTypeName, i, typeName));
  const TString xTitle(Form("Hit X for %s %d in %s;", aliTypeName, i, typeName));
  const TString yTitle(Form("Hit Y for %s %d in %s;", aliTypeName, i, typeName));
  const TString zTitle(Form("Hit Z for %s %d in %s;", aliTypeName, i, typeName));
  const TString dttimeTitle(Form("DT Time for %s %d in %s;", aliTypeName, i, typeName));
  const TString dudtTitle(Form("#Delta_{u} vs. DT Time for %s %d in %s;", aliTypeName, i, typeName));
  const TString dwdtTitle(Form("#Delta_{w} vs. DT Time for %s %d in %s;", aliTypeName, i, typeName));
  const TString dttimeerrTitle(Form("DT Time Error for %s %d in %s;", aliTypeName, i, typeName));
  const TString ndtTitle(Form("nDT for %s %d in %s;", aliTypeName, i, typeName));
  const TString hcaltimeTitle(Form("HCAL Time for %s %d in %s;", aliTypeName, i, typeName));
  const TString etaTitle(Form("Hit Eta for %s %d in %s;", aliTypeName, i, typeName));
  const TString resxprimeoverthetaTitle(Form("x Residual Over #Delta tan(#theta) for %s %d in %s;", aliTypeName, i, typeName));
  const TString resxprimeoverthetaThickTitle(Form("x Residual Over #Delta tan(#theta) (thick) for %s %d in %s;", aliTypeName, i, typeName));
  const TString resxprimeoverthetaThinTitle(Form("x Residual Over #Delta tan(#theta) (thin) for %s %d in %s;", aliTypeName, i, typeName));
  const TString resxprimevsthetaTitle(Form("x Residual vs. #Delta tan(#theta) for %s %d in %s;", aliTypeName, i, typeName));
  const TString dutrkTitle(Form("#Delta_{u} vs. tan(#theta_{trk}) for %s %d in %s;", aliTypeName, i, typeName));
  const TString dulorTitle(Form("#Delta_{u} vs. tan(#theta_{LA}) for %s %d in %s;", aliTypeName, i, typeName));

  const TString chargeTitle(Form("Cluster Charge for %s %d in %s;", aliTypeName, i, typeName));
  const TString chargedttimeTitle(Form("Cluster Charge vs. DT Time for %s %d in %s;", aliTypeName, i, typeName));
  const TString nstripsTitle(Form("Cluster NStrips for %s %d in %s;", aliTypeName, i, typeName));
 
  const ModuleHistos &xTitHists = this->getHistStructFromMap(aliDetId); // for x-axis titles

  if(bookTH2_){
    this->getBinning(aliDetId.subdetId(), XprimevsThetaResidual, nbins, xmin, xmax);
    if(debug_) cout<<"nbins "<<nbins<<" xmin "<<xmin<<" xmax "<<xmax<<endl;
    if(debug_) cout<<"Booking summary hist "<<resxprimevsthetaTitle<<endl;

    sumContainer.sumResXprimeVsTheta_ = tfd.make<TH2F>(Form("h_resXprimeVsTheta_%s_%d", aliTypeName, i),
						       resxprimevsthetaTitle,
						       40, -2, 2, 40, -200, 200);

    
    sumContainer.sumduvstrkangle_ = tfd.make<TH2F>(Form("h_duvstrkangle_%s_%d", aliTypeName, i),
						   dutrkTitle,
						   40, -2, 2, 40, -200, 200);
    
//       sumContainer.sumduvslorangle_ = tfd.make<TH2F>(Form("h_duvslorangle_%s_%d", aliTypeName, i),
// 						     dulorTitle,
// 						     50, -0.5, 0.5, nbins, xmin, xmax);
    
    

    sumContainer.sumduvsdttime_ = tfd.make<TH2F>(Form("h_duvsdttime_%s_%d", aliTypeName, i),
						 dudtTitle,
						 50, -25, 25, 50, -500, 500);
    
    sumContainer.sumdwvsdttime_ = tfd.make<TH2F>(Form("h_dwvsdttime_%s_%d", aliTypeName, i),
						 dwdtTitle,
						 50, -25, 25, 50, -1000, 1000);

    sumContainer.sumchargevsdttime_ = tfd.make<TH2F>(Form("h_chargevsdttime_%s_%d", aliTypeName, i),
						     chargedttimeTitle,
						     50, -25, 25, 50, 0, 1000);
    
    /*
    sumContainer.sumuz_ = tfd.make<TH2F>(Form("h_uz_%s_%d", aliTypeName, i),
					 uzTitle,
					 12, -110, 110, 3, -1.5, 1.5);

    sumContainer.sumvz_ = tfd.make<TH2F>(Form("h_vz_%s_%d", aliTypeName, i),
					 vzTitle,
					 12, -110, 110, 3, -1.5, 1.5);

    sumContainer.sumwz_ = tfd.make<TH2F>(Form("h_wz_%s_%d", aliTypeName, i),
					 wzTitle,
					 12, -110, 110, 3, -1.5, 1.5);

    */
    if(debug_) cout<<"Booked summary hist"<<endl;
  }

  if(bookTH1_){
    if(debug_) cout<<"Begin booking TH1 summary histos"<<endl;
    this->getBinning(aliDetId.subdetId(), XprimeResidual, nbins, xmin, xmax);
   
    sumContainer.sumXResiduals_ = tfd.make<TH1F>(Form("h_Xprime_%s_%d", aliTypeName, i),
						 sumTitle + xTitHists.ResXprimeHisto->GetXaxis()->GetTitle(),
						 nbins, xmin, xmax);

    for(int ihist=0;ihist<2;ihist++){
      sumContainer.sumXResiduals2_[ihist] = tfd.make<TH1F>(Form("h_Xprime_%i_%s_%d",ihist, aliTypeName, i),
							   sumTitle + xTitHists.ResXprimeHisto->GetXaxis()->GetTitle(),
							   nbins, xmin, xmax);
    
      sumContainer.sumResXprimeOverTheta2_[ihist] = tfd.make<TH1F>(Form("h_resXprimeOverTheta_%i_%s_%d",ihist, aliTypeName, i),
								   resxprimeoverthetaTitle,
								   100, -1000, 1000);
    }

    sumContainer.sumdttime_ = tfd.make<TH1F>(Form("h_dttime_%s_%d", aliTypeName, i),
					     dttimeTitle,
					     200, -100, 100);

    sumContainer.sumndt_ = tfd.make<TH1F>(Form("h_ndt_%s_%d", aliTypeName, i),
					     ndtTitle,
					     100, 0, 100);

    sumContainer.sumcharge_ = tfd.make<TH1F>(Form("h_charge_%s_%d", aliTypeName, i),
					     chargeTitle,
					     100, 0, 1000);

    sumContainer.sumnstrips_ = tfd.make<TH1F>(Form("h_nstrips_%s_%d", aliTypeName, i),
					      nstripsTitle,
					      10, 0, 10);
    
    sumContainer.sumdttimeerr_ = tfd.make<TH1F>(Form("h_dttimeerr_%s_%d", aliTypeName, i),
						dttimeerrTitle,
						100, 0, 100);
    
    
    sumContainer.sumhcaltime_ = tfd.make<TH1F>(Form("h_hcaltime_%s_%d", aliTypeName, i),
					       hcaltimeTitle,
					       200, -100, 100);
    
    if(debug_) cout<<"Booked summary histos A"<<endl;

    sumContainer.sumTanTrackAngle_ = tfd.make<TH1F>(Form("h_tanTrackAngle_%s_%d", aliTypeName, i),
						    taTitle,
						    100, -5, 5);

    sumContainer.sumdeltaTan_ = tfd.make<TH1F>(Form("h_deltatan_%s_%d", aliTypeName, i),
					       dtTitle,
					       100, -5, 5);

    sumContainer.sumTanLorentzAngle_ = tfd.make<TH1F>(Form("h_tanLorentzAngle_%s_%d", aliTypeName, i),
						      laTitle,
						      100, -0.5, 0.5);

  
    if(debug_) cout<<"Booked summary histos B"<<endl;

    sumContainer.sumResXprimeOverTheta_ = tfd.make<TH1F>(Form("h_resXprimeOverTheta_%s_%d", aliTypeName, i),
							 resxprimeoverthetaTitle,
							 100, -1000, 1000);

    sumContainer.sumResXprimeOverTheta_thick_ = tfd.make<TH1F>(Form("h_resXprimeOverThetaThick_%s_%d", aliTypeName, i),
							       resxprimeoverthetaTitle,
							       100, -1000, 1000);

    sumContainer.sumResXprimeOverTheta_thin_  = tfd.make<TH1F>(Form("h_resXprimeOverThetaThin_%s_%d", aliTypeName, i),
							       resxprimeoverthetaTitle,
							       100, -1000, 1000);
    
    sumContainer.sumuOrientation_ = tfd.make<TH1F>(Form("h_uOrientation_%s_%d", aliTypeName, i),
						   uoTitle,
						   3, -1.5, 1.5);

    sumContainer.sumvOrientation_ = tfd.make<TH1F>(Form("h_vOrientation_%s_%d", aliTypeName, i),
						   voTitle,
						   3, -1.5, 1.5);

//     sumContainer.sumvminusdusign_ = tfd.make<TH1F>(Form("h_vminusdusign_%s_%d", aliTypeName, i),
// 						   vdTitle,
// 						   5, -2.5, 2.5);
  
    sumContainer.sumwOrientation_ = tfd.make<TH1F>(Form("h_wOrientation_%s_%d", aliTypeName, i),
						   woTitle,
						   3, -1.5, 1.5);

    sumContainer.sumlocaltheta_ = tfd.make<TH1F>(Form("h_localtheta_%s_%d", aliTypeName, i),
						   ltTitle,
						   100, -5, 5);

    if(debug_) cout<<"Booked summary histos C"<<endl;

    sumContainer.sumlocalphi_ = tfd.make<TH1F>(Form("h_localphi_%s_%d", aliTypeName, i),
						   lpTitle,
						   100, -1, 1);

    /*
    for(int ihist=0;ihist<12;ihist++){
      sumContainer.sumdu_z_[ihist] = tfd.make<TH1F>(Form("h_duz_%i_%s_%d", ihist, aliTypeName, i),
						    Form("delta(u) for Module %i for %s %d in %s;", ihist, aliTypeName, i, typeName),
						    50, -500, 500);

      sumContainer.sumdw_z_[ihist] = tfd.make<TH1F>(Form("h_dwz_%i_%s_%d", ihist, aliTypeName, i),
						    Form("delta(w) for Module %i for %s %d in %s;", ihist, aliTypeName, i, typeName),
						    50, -1000, 1000);

    }
    */

    if(debug_) cout<<"Booked summary histos D"<<endl;
    /*
    sumContainer.sumHitX_ = tfd.make<TH1F>(Form("h_hitx_%s_%d", aliTypeName, i),
					   xTitle,
					   1000, -100, 100);

    sumContainer.sumHitY_ = tfd.make<TH1F>(Form("h_hity_%s_%d", aliTypeName, i),
					   yTitle,
					   1000, -100, 100);
    */
    sumContainer.sumHitZ_ = tfd.make<TH1F>(Form("h_hitz_%s_%d", aliTypeName, i),
					   zTitle,
					   100, -100, 100);

    sumContainer.sumHitEta_ = tfd.make<TH1F>(Form("h_hiteta_%s_%d", aliTypeName, i),
					     etaTitle,
					     100, -5, 5);
  

    if(debug_) cout<<"Booked all summary histos"<<endl;
    //cout<<"Make sumXResiduals "<<Form("h_Xprime_%s_%d", aliTypeName, i)<<" title "<<sumTitle<<endl;
    /*
    this->getBinning(aliDetId.subdetId(), NormXprimeResidual, nbins, xmin, xmax);
    sumContainer.sumNormXResiduals_ = tfd.make<TH1F>(Form("h_NormXprime_%s_%d",aliTypeName,i), 
						     sumTitle + xTitHists.NormResXprimeHisto->GetXaxis()->GetTitle(),
						     nbins, xmin, xmax);
    if (bookResidY) {
      this->getBinning(aliDetId.subdetId(), YprimeResidual, nbins, xmin, xmax);
      sumContainer.sumYResiduals_ = tfd.make<TH1F>(Form("h_Yprime_%s_%d",aliTypeName,i), 
						   sumTitle + xTitHists.ResYprimeHisto->GetXaxis()->GetTitle(),
						   nbins, xmin, xmax);
    
      this->getBinning(aliDetId.subdetId(), NormYprimeResidual, nbins, xmin, xmax);
      sumContainer.sumNormYResiduals_ = tfd.make<TH1F>(Form("h_NormYprime_%s_%d",aliTypeName,i), 
						       sumTitle + xTitHists.NormResYprimeHisto->GetXaxis()->GetTitle(),
						       nbins, xmin, xmax);
    }
    */
  }
  // If we are at the lowest level, we already sum up and fill the summary.

  // special case I: For DetUnits and Detwith  only one subcomponent start filling summary histos
  bool adet=(subtype == align::AlignableDet);
  bool adetu=(subtype  == align::AlignableDetUnit); 

  if(debug_) cout<<"subcompSize "<<subcompSize<<" alignableDet "<<adet<<" alignableDetUnit "<<adetu<<endl;
  if( (  subtype == align::AlignableDet && subcompSize == 1) || subtype  == align::AlignableDetUnit ) {
    if(debug_) cout<<"Sum histo type 1"<<endl;
   for(uint k = 0; k < aliSize; ++k) {
      if(debug_) cout<<"det "<<k<<endl;
      DetId detid = ali.components()[k]->id();
      ModuleHistos &histStruct = this->getHistStructFromMap(detid);
      if(debug_) cout<<"summarizeBinInContainer"<<endl;
      this->summarizeBinInContainer(k+1, detid.subdetId() ,sumContainer, histStruct );
      if(bookTH2_){
	if(debug_) cout<<"Add ResXprimeVsThetaHisto to sumResXprimeVsTheta"<<endl;
	sumContainer.sumResXprimeVsTheta_->Add(histStruct.ResXprimeVsThetaHisto);
	sumContainer.sumduvstrkangle_->Add(histStruct.duvstrkangleHisto);
	//sumContainer.sumduvslorangle_->Add(histStruct.duvslorangleHisto);
	//sumContainer.sumuz_->Add(histStruct.uzHisto);
	//sumContainer.sumvz_->Add(histStruct.vzHisto);
	//sumContainer.sumwz_->Add(histStruct.wzHisto);
	sumContainer.sumduvsdttime_->Add(histStruct.duvsdttimeHisto);
	sumContainer.sumdwvsdttime_->Add(histStruct.dwvsdttimeHisto);
	sumContainer.sumchargevsdttime_->Add(histStruct.chargevsdttimeHisto);
      }
      if(bookTH1_){
	sumContainer.sumXResiduals_->Add(histStruct.ResXprimeHisto);
	//vDeleteObjects_.push_back(histStruct.ResXprimeHisto);
	for(int ihist=0;ihist<2;ihist++){
	  sumContainer.sumXResiduals2_[ihist]->Add(histStruct.ResXprimeHisto2[ihist]);
	  sumContainer.sumResXprimeOverTheta2_[ihist]->Add(histStruct.ResXprimeOverThetaHisto2[ihist]);
	}
	sumContainer.sumTanTrackAngle_->Add(histStruct.tanTrackAngleHisto);
	sumContainer.sumTanLorentzAngle_->Add(histStruct.tanLorentzAngleHisto);
	sumContainer.sumdeltaTan_->Add(histStruct.deltaTanHisto);
	//sumContainer.sumNormXResiduals_->Add(histStruct.NormResXprimeHisto);
	sumContainer.sumResXprimeOverTheta_->Add(histStruct.ResXprimeOverThetaHisto);
	sumContainer.sumResXprimeOverTheta_thick_->Add(histStruct.ResXprimeOverThetaHisto_thick);
	sumContainer.sumResXprimeOverTheta_thin_ ->Add(histStruct.ResXprimeOverThetaHisto_thin);
	sumContainer.sumuOrientation_->Add(histStruct.uOrientationHisto);
	sumContainer.sumvOrientation_->Add(histStruct.vOrientationHisto);
	sumContainer.sumwOrientation_->Add(histStruct.wOrientationHisto);
	sumContainer.sumlocaltheta_->Add(histStruct.localthetaHisto);
	sumContainer.sumlocalphi_->Add(histStruct.localphiHisto);
	//sumContainer.sumvminusdusign_->Add(histStruct.vminusdusignHisto);
	sumContainer.sumHitEta_->Add(histStruct.HitEtaHisto);
	//sumContainer.sumHitX_->Add(histStruct.HitXHisto);
	//sumContainer.sumHitY_->Add(histStruct.HitYHisto);
	/*      
 	for(int ihist=0;ihist<12;ihist++){
	  sumContainer.sumdu_z_[ihist]->Add(histStruct.du_zHisto[ihist]);
	  sumContainer.sumdw_z_[ihist]->Add(histStruct.dw_zHisto[ihist]);
	}
	*/
	sumContainer.sumHitZ_->Add(histStruct.HitZHisto);
	sumContainer.sumdttime_->Add(histStruct.dttimeHisto);
	sumContainer.sumdttimeerr_->Add(histStruct.dttimeerrHisto);
	sumContainer.sumndt_->Add(histStruct.ndtHisto);
	sumContainer.sumcharge_->Add(histStruct.chargeHisto);
	sumContainer.sumnstrips_->Add(histStruct.nstripsHisto);
	sumContainer.sumhcaltime_->Add(histStruct.hcaltimeHisto);
	
	/*
	if( this->isPixel(detid.subdetId()) || stripYResiduals_ ) {
	  sumContainer.sumYResiduals_->Add(histStruct.ResYprimeHisto);
	  sumContainer.sumNormYResiduals_->Add(histStruct.NormResYprimeHisto);
	}
	*/
      }
    }
  } else if( subtype == align::AlignableDet && subcompSize > 1) { // fixed: was aliSize before
    if(debug_) cout<<"Sum histo type 2"<<endl;
    // special case II: Fill summary histos for dets with two detunits 
    for(uint k = 0; k < aliSize; ++k) {
      for(uint j = 0; j < subcompSize; ++j) { // assumes all have same size (as binning does)
	if(debug_) cout<<"det "<<k<<","<<j<<endl;
	DetId detid = ali.components()[k]->components()[j]->id();
	ModuleHistos &histStruct = this->getHistStructFromMap(detid);	
	this->summarizeBinInContainer(2*k+j+1, detid.subdetId() ,sumContainer, histStruct );
	if(bookTH2_){
	  if(debug_) cout<<"Add ResXprimeVsThetaHisto to sumResXprimeVsTheta"<<endl;
	  sumContainer.sumResXprimeVsTheta_->Add(histStruct.ResXprimeVsThetaHisto);
	  sumContainer.sumduvstrkangle_->Add(histStruct.duvstrkangleHisto);
	  //sumContainer.sumduvslorangle_->Add(histStruct.duvslorangleHisto);
	  //sumContainer.sumuz_->Add(histStruct.uzHisto);
	  //sumContainer.sumvz_->Add(histStruct.vzHisto);
	  //sumContainer.sumwz_->Add(histStruct.wzHisto);
	  sumContainer.sumduvsdttime_->Add(histStruct.duvsdttimeHisto);
	  sumContainer.sumdwvsdttime_->Add(histStruct.dwvsdttimeHisto);
	  sumContainer.sumchargevsdttime_->Add(histStruct.chargevsdttimeHisto);
	}
	if(bookTH1_){
	  if(debug_)cout<<"Begin summing containers"<<endl;
	  sumContainer.sumXResiduals_->Add( histStruct.ResXprimeHisto);
	  for(int ihist=0;ihist<2;ihist++){
	    sumContainer.sumXResiduals2_[ihist]->Add(histStruct.ResXprimeHisto2[ihist]);
	    sumContainer.sumResXprimeOverTheta2_[ihist]->Add(histStruct.ResXprimeOverThetaHisto2[ihist]);
	  }
	  sumContainer.sumTanTrackAngle_->Add(histStruct.tanTrackAngleHisto);
	  sumContainer.sumTanLorentzAngle_->Add(histStruct.tanLorentzAngleHisto);
	  sumContainer.sumdeltaTan_->Add(histStruct.deltaTanHisto);
	  //sumContainer.sumNormXResiduals_->Add(histStruct.NormResXprimeHisto);
	  sumContainer.sumResXprimeOverTheta_->Add(histStruct.ResXprimeOverThetaHisto);
	  sumContainer.sumResXprimeOverTheta_thick_->Add(histStruct.ResXprimeOverThetaHisto_thick);
	  sumContainer.sumResXprimeOverTheta_thin_ ->Add(histStruct.ResXprimeOverThetaHisto_thin);
	  sumContainer.sumuOrientation_->Add(histStruct.uOrientationHisto);
	  sumContainer.sumvOrientation_->Add(histStruct.vOrientationHisto);
	  sumContainer.sumwOrientation_->Add(histStruct.wOrientationHisto);
	  sumContainer.sumlocaltheta_->Add(histStruct.localthetaHisto);
	  sumContainer.sumlocalphi_->Add(histStruct.localphiHisto);
	  //sumContainer.sumvminusdusign_->Add(histStruct.vminusdusignHisto);	
	  sumContainer.sumHitEta_->Add(histStruct.HitEtaHisto);
	  //sumContainer.sumHitX_->Add(histStruct.HitXHisto);
	  //sumContainer.sumHitY_->Add(histStruct.HitYHisto);
	  sumContainer.sumdttime_->Add(histStruct.dttimeHisto);
	  sumContainer.sumdttimeerr_->Add(histStruct.dttimeerrHisto);
	  sumContainer.sumndt_->Add(histStruct.ndtHisto);
	  sumContainer.sumhcaltime_->Add(histStruct.hcaltimeHisto);
	  sumContainer.sumcharge_->Add(histStruct.chargeHisto);
	  sumContainer.sumnstrips_->Add(histStruct.nstripsHisto);
	  /*
	  for(int ihist=0;ihist<12;ihist++){
	    sumContainer.sumdu_z_[ihist]->Add(histStruct.du_zHisto[ihist]);
	    sumContainer.sumdw_z_[ihist]->Add(histStruct.dw_zHisto[ihist]);
	  }
	  */
	  sumContainer.sumHitZ_->Add(histStruct.HitZHisto);
	
	  /*
	  if( this->isPixel(detid.subdetId()) || stripYResiduals_ ) {
	    sumContainer.sumYResiduals_->Add( histStruct.ResYprimeHisto);
	    sumContainer.sumNormYResiduals_->Add( histStruct.NormResYprimeHisto);
	  }
	  */
	  if(debug_)cout<<"Finished summing containers"<<endl;
	}
      }
    }
  }else if(debug_) cout<<"skipping"<<endl;
  if(debug_) cout<<"returning sumContainer"<<endl;

  return sumContainer;
}


float 
TrackerOfflineValidation::Fwhm (const TH1* hist) const
{
  float max = hist->GetMaximum();
  int left = -1, right = -1;
  for(unsigned int i = 1, iEnd = hist->GetNbinsX(); i <= iEnd; ++i) {
    if(hist->GetBinContent(i) < max/2. && hist->GetBinContent(i+1) > max/2. && left == -1) {
      if(max/2. - hist->GetBinContent(i) < hist->GetBinContent(i+1) - max/2.) {
	left = i;
	++i;
      } else {
	left = i+1;
	++i;
      }
    }
    if(left != -1 && right == -1) {
      if(hist->GetBinContent(i) > max/2. && hist->GetBinContent(i+1) < max/2.) {
	if( hist->GetBinContent(i) - max/2. < max/2. - hist->GetBinContent(i+1)) {
	  right = i;
	} else {
	  right = i+1;
	}
	
      }
    }
  }
  return hist->GetXaxis()->GetBinCenter(right) - hist->GetXaxis()->GetBinCenter(left);
}

////////////////////////////////////////////////////////////////////////////////////
void 
TrackerOfflineValidation::fillTree(TTree &tree,
				   const std::map<int, TrackerOfflineValidation::ModuleHistos> &moduleHist_,
				   TkOffTreeVariables &treeMem, const TrackerGeometry &tkgeom)
{
 
  for(std::map<int, TrackerOfflineValidation::ModuleHistos>::const_iterator it = moduleHist_.begin(), 
	itEnd= moduleHist_.end(); it != itEnd;++it ) { 
    treeMem.clear(); // make empty/default
    //variables concerning the tracker components/hierarchy levels
    DetId detId_ = it->first;
    treeMem.moduleId = detId_;
    treeMem.subDetId = detId_.subdetId();
    treeMem.isDoubleSide =0;

    if(treeMem.subDetId == PixelSubdetector::PixelBarrel){
      PXBDetId pxbId(detId_); 
      treeMem.layer = pxbId.layer(); 
      treeMem.rod = pxbId.ladder();
  
    } else if(treeMem.subDetId == PixelSubdetector::PixelEndcap){
      PXFDetId pxfId(detId_); 
      treeMem.layer = pxfId.disk(); 
      treeMem.side = pxfId.side();
      treeMem.blade = pxfId.blade(); 
      treeMem.panel = pxfId.panel();

    } else if(treeMem.subDetId == StripSubdetector::TIB){
      TIBDetId tibId(detId_); 
      treeMem.layer = tibId.layer(); 
      treeMem.side = tibId.string()[0];
      treeMem.rod = tibId.string()[2]; 
      treeMem.outerInner = tibId.string()[1]; 
      treeMem.isStereo = tibId.stereo();
      treeMem.isDoubleSide = tibId.isDoubleSide();
    } else if(treeMem.subDetId == StripSubdetector::TID){
      TIDDetId tidId(detId_); 
      treeMem.layer = tidId.wheel(); 
      treeMem.side = tidId.side();
      treeMem.ring = tidId.ring(); 
      treeMem.outerInner = tidId.module()[0]; 
      treeMem.isStereo = tidId.stereo();
      treeMem.isDoubleSide = tidId.isDoubleSide();
    } else if(treeMem.subDetId == StripSubdetector::TOB){
      TOBDetId tobId(detId_); 
      treeMem.layer = tobId.layer(); 
      treeMem.side = tobId.rod()[0];
      treeMem.rod = tobId.rod()[1]; 
      treeMem.isStereo = tobId.stereo();
      treeMem.isDoubleSide = tobId.isDoubleSide();
    } else if(treeMem.subDetId == StripSubdetector::TEC) {
      TECDetId tecId(detId_); 
      treeMem.layer = tecId.wheel(); 
      treeMem.side  = tecId.side();
      treeMem.ring  = tecId.ring(); 
      treeMem.petal = tecId.petal()[1]; 
      treeMem.outerInner = tecId.petal()[0];
      treeMem.isStereo = tecId.stereo();
      treeMem.isDoubleSide = tecId.isDoubleSide(); 
    }
    
    //variables concerning the tracker geometry
    
    const Surface::PositionType &gPModule = tkgeom.idToDet(detId_)->position();
    treeMem.posPhi = gPModule.phi();
    treeMem.posEta = gPModule.eta();
    treeMem.posR   = gPModule.perp();
    treeMem.posX   = gPModule.x();
    treeMem.posY   = gPModule.y();
    treeMem.posZ   = gPModule.z();
 
    const Surface& surface =  tkgeom.idToDet(detId_)->surface();
    
    //global Orientation of local coordinate system of dets/detUnits   
    LocalPoint  lUDirection(1.,0.,0.), lVDirection(0.,1.,0.), lWDirection(0.,0.,1.);
    GlobalPoint gUDirection = surface.toGlobal(lUDirection),
                gVDirection = surface.toGlobal(lVDirection),
		gWDirection = surface.toGlobal(lWDirection);
    double dR(999.), dPhi(999.), dZ(999.);
    if(treeMem.subDetId==PixelSubdetector::PixelBarrel || treeMem.subDetId==StripSubdetector::TIB || treeMem.subDetId==StripSubdetector::TOB){
      dR = gWDirection.perp() - gPModule.perp();
      dPhi = deltaPhi(gUDirection.phi(),gPModule.phi());
      dZ = gVDirection.z() - gPModule.z();
      if(dZ>=0.)treeMem.rOrZDirection = 1; else treeMem.rOrZDirection = -1;
    }else if(treeMem.subDetId==PixelSubdetector::PixelEndcap){
      dR = gUDirection.perp() - gPModule.perp();
      dPhi = deltaPhi(gVDirection.phi(),gPModule.phi());
      dZ = gWDirection.z() - gPModule.z();
      if(dR>=0.)treeMem.rOrZDirection = 1; else treeMem.rOrZDirection = -1;
    }else if(treeMem.subDetId==StripSubdetector::TID || treeMem.subDetId==StripSubdetector::TEC){
      dR = gVDirection.perp() - gPModule.perp();
      dPhi = deltaPhi(gUDirection.phi(),gPModule.phi());
      dZ = gWDirection.z() - gPModule.z();
      if(dR>=0.)treeMem.rOrZDirection = 1; else treeMem.rOrZDirection = -1;
    }
    if(dR>=0.)treeMem.rDirection = 1; else treeMem.rDirection = -1;
    if(dPhi>=0.)treeMem.phiDirection = 1; else treeMem.phiDirection = -1;
    if(dZ>=0.)treeMem.zDirection = 1; else treeMem.zDirection = -1;
    
    
    //mean and RMS values (extracted from histograms(Xprime on module level)
    treeMem.entries = static_cast<UInt_t>(it->second.ResXprimeHisto->GetEntries());
    treeMem.meanX = it->second.ResXprimeHisto->GetMean();
    treeMem.rmsX  = it->second.ResXprimeHisto->GetRMS();
    //treeMem.sigmaX = Fwhm(it->second.ResXprimeHisto)/2.355;
    if (useFit_) {
      
      //call fit function which returns mean and sigma from the fit
      //for absolute residuals
      std::pair<float,float> fitResult1 = this->fitResiduals(it->second.ResXprimeHisto);
      treeMem.fitMeanX = fitResult1.first;
      treeMem.fitSigmaX = fitResult1.second;
      //for normalized residuals
      //std::pair<float,float> fitResult2 = this->fitResiduals(it->second.NormResXprimeHisto);
      //treeMem.fitMeanNormX = fitResult2.first;
      //treeMem.fitSigmaNormX = fitResult2.second;
    }
    //get median for absolute residuals
    treeMem.medianX   = this->getMedian(it->second.ResXprimeHisto);


    int numberOfBins=it->second.ResXprimeHisto->GetNbinsX();
    treeMem.numberOfUnderflows = it->second.ResXprimeHisto->GetBinContent(0);
    treeMem.numberOfOverflows = it->second.ResXprimeHisto->GetBinContent(numberOfBins+1);
    treeMem.numberOfOutliers =  it->second.ResXprimeHisto->GetBinContent(0)+it->second.ResXprimeHisto->GetBinContent(numberOfBins+1);
    //mean and RMS values (extracted from histograms(normalized Xprime on module level)
    //treeMem.meanNormX = it->second.NormResXprimeHisto->GetMean();
    //treeMem.rmsNormX = it->second.NormResXprimeHisto->GetRMS();

    double stats[20];
    //it->second.NormResXprimeHisto->GetStats(stats);
    // GF  treeMem.chi2PerDofX = stats[3]/(stats[0]-1);
    if (stats[0]) treeMem.chi2PerDofX = stats[3]/stats[0];
    
    //treeMem.sigmaNormX = Fwhm(it->second.NormResXprimeHisto)/2.355;
    treeMem.histNameX = it->second.ResXprimeHisto->GetName();
    //treeMem.histNameNormX = it->second.NormResXprimeHisto->GetName();
    

    // fill tree variables in local coordinates if set in cfg
    if(lCoorHistOn_) {
      //treeMem.meanLocalX = it->second.ResHisto->GetMean();
      //treeMem.rmsLocalX = it->second.ResHisto->GetRMS();
      //treeMem.meanNormLocalX = it->second.NormResHisto->GetMean();
      //treeMem.rmsNormLocalX = it->second.NormResHisto->GetRMS();

      //treeMem.histNameLocalX = it->second.ResHisto->GetName();
      //treeMem.histNameNormLocalX = it->second.NormResHisto->GetName();
    }

    // mean and RMS values in local y (extracted from histograms(normalized Yprime on module level)
    // might exist in pixel only
    /*   
    if (it->second.ResYprimeHisto) {//(stripYResiduals_){
      TH1 *h = it->second.ResYprimeHisto;
      treeMem.meanY = h->GetMean();
      treeMem.rmsY  = h->GetRMS();
      
      if (useFit_) { // fit function which returns mean and sigma from the fit
	std::pair<float,float> fitMeanSigma = this->fitResiduals(h);
	treeMem.fitMeanY  = fitMeanSigma.first;
	treeMem.fitSigmaY = fitMeanSigma.second;
      }
      //get median for absolute residuals
      treeMem.medianY   = this->getMedian(h);

      treeMem.histNameY = h->GetName();
    }
    if (it->second.NormResYprimeHisto) {
      //TH1 *h = it->second.NormResYprimeHisto;
      //treeMem.meanNormY = h->GetMean();
      //treeMem.rmsNormY  = h->GetRMS();
      h->GetStats(stats); // stats buffer defined above
      if (stats[0]) treeMem.chi2PerDofY = stats[3]/stats[0];

      if (useFit_) { // fit function which returns mean and sigma from the fit
	std::pair<float,float> fitMeanSigma = this->fitResiduals(h);
	treeMem.fitMeanNormY  = fitMeanSigma.first;
	treeMem.fitSigmaNormY = fitMeanSigma.second;
      }
      treeMem.histNameNormY = h->GetName();
    }
    */
    tree.Fill();
  }
}

std::pair<float,float> 
TrackerOfflineValidation::fitResiduals(TH1 *hist) const
{
  std::pair<float,float> fitResult(9999., 9999.);
  if (!hist || hist->GetEntries() < 20) return fitResult;

  float mean  = hist->GetMean();
  float sigma = hist->GetRMS();

  try { // for < CMSSW_2_2_0 since ROOT warnings from fit are converted to exceptions
    // Remove the try/catch for more recent CMSSW!
    // first fit: two RMS around mean
    TF1 func("tmp", "gaus", mean - 2.*sigma, mean + 2.*sigma); 
    if (0 == hist->Fit(&func,"QNR")) { // N: do not blow up file by storing fit!
      mean  = func.GetParameter(1);
      sigma = func.GetParameter(2);
      // second fit: three sigma of first fit around mean of first fit
      func.SetRange(mean - 3.*sigma, mean + 3.*sigma);
      // I: integral gives more correct results if binning is too wide
      // L: Likelihood can treat empty bins correctly (if hist not weighted...)
      if (0 == hist->Fit(&func, "Q0LR")) {
	if (hist->GetFunction(func.GetName())) { // Take care that it is later on drawn:
	  hist->GetFunction(func.GetName())->ResetBit(TF1::kNotDraw);
	}
	fitResult.first = func.GetParameter(1);
	fitResult.second = func.GetParameter(2);
      }
    }
  } catch (cms::Exception const & e) {
    edm::LogWarning("Alignment") << "@SUB=TrackerOfflineValidation::fitResiduals"
				 << "Caught this exception during ROOT fit: "
				 << e.what();
  }
  
  return fitResult;
}
float 
TrackerOfflineValidation::getMedian(const TH1 *histo) const
{

  float median = 999;
  int nbins = histo->GetNbinsX();

 
  //extract median from histogram
  double *x = new double[nbins];
  double *y = new double[nbins];
  for (int j = 0; j < nbins; j++) {
    x[j] = histo->GetBinCenter(j+1);
    y[j] = histo->GetBinContent(j+1);
  }
  median = TMath::Median(nbins, x, y);
  

  delete[] x; x = 0;
  delete [] y; y = 0;  

  return median;

}
//define this as a plug-in
DEFINE_FWK_MODULE(TrackerOfflineValidation);
