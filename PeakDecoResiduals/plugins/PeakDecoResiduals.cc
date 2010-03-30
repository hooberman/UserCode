// -*- C++ -*-
//
// Package:    PeakDecoResiduals
// Class:      PeakDecoResiduals
// 
/**\class PeakDecoResiduals PeakDecoResiduals.cc Alignment/Validator/src/PeakDecoResiduals.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Erik Butz
//         Created:  Tue Dec 11 14:03:05 CET 2007
// $Id: PeakDecoResiduals.cc,v 1.1 2010/03/17 15:46:48 benhoob Exp $
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
#include "Alignment/PeakDecoResiduals/interface/TrackerValidationVariables.h"
#include "Alignment/PeakDecoResiduals/interface/TkOffTreeVariables.h"
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/UtilAlgos/interface/TFileDirectory.h"
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

class PeakDecoResiduals : public edm::EDAnalyzer {
public:
  explicit PeakDecoResiduals(const edm::ParameterSet&);
  ~PeakDecoResiduals();
  
private:

  TH2F* hlumivsrun;      
  TH1F* hdttime[7];            
  TH1F* hdttimeerr[7];         
  TH1F* hndt[7];               
  TH2F* hduvsdtantheta[7];
  TH2F* hresxvsdtantheta[7];
  TH2F* hduvstrktheta[7];      
  TH2F* hduvsdtantheta_vw[7][2][2];  
  TH2F* hduvstrktheta_vw[7][2][2];   
  TH2F* hduvsdtantheta_vp_dt[7][8];  
  TH2F* hduvsdtantheta_vm_dt[7][8];
  TH2F* hchargevsdttime[7];    
  TH2F* hduvsdttime[7];        
  TH2F* hdwvsdttime[7];        
  TH1F* hhitx[7];              
  TH1F* hhity[7];
  TH1F* hhitz[7];
  TH1F* hcharge[7];            
  TH1F* hnstrips[7];           
  TH1F* huOrientation[7];      
  TH1F* hvOrientation[7];      
  TH1F* hwOrientation[7];      
  TH1F* hlocaltheta[7];        
  TH1F* hlocalphi[7];          
  TH1F* hdu[7];                
  TH1F* hdw[7];        
  TH1F* hresxoverdtanth[7];        
  TH1F* htantrktheta[7][3];       
  TH1F* htanlorangle[7][3];       
  TH1F* hdtantheta[7][3];         

  // ------------- private member function -------------
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  void BookHists(TFileDirectory &tfd);
  void makeTree(TFileDirectory &tfd);
   
  int  getSubdetInteger(const DetId& detid); 
  bool isBarrel(uint32_t subDetId);
  bool isEndCap(uint32_t subDetId);
  bool isPixel(uint32_t subDetId);
  bool isDetOrDetUnit(align::StructureType type);

  // From MillePedeAlignmentMonitor: Get Index for Arbitary vector<class> by name
  template <class OBJECT_TYPE>  
  int GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name);

  // ---------- member data ---------------------------


  const edm::ParameterSet parset_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  const TrackerGeometry *bareTkGeomPtr_; // ugly hack to book hists only once, but check 

  // parameters from cfg to steer
  bool debug_;
  bool runOnCosmics_;
  bool createTree_;

  //variables for tree
  TTree* outTree;
  Int_t   subdet_;
  Float_t du_;
  Float_t dw_;
  Float_t dtanth_;
  Float_t tantrk_;
  Float_t tanla_;
  Float_t dttime_;
  Float_t dttimeerr_;
  Int_t   ndt_;
  Float_t charge_;
  Int_t   nstrips_;
  Int_t   u_;
  Int_t   v_;
  Int_t   w_;
  Float_t x_;
  Float_t y_;
  Float_t z_;
  Float_t r_;
  Int_t   ring_;
  Int_t   wheel_; 
  Int_t   petal_;
  Int_t   side_;
  Int_t   stereo_;
  Int_t   ds_;
  Int_t   zplus_;
  Int_t   back_;
  Int_t   front_;
  Int_t   layer_;
  Int_t   rod_;
 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
template <class OBJECT_TYPE>  
int PeakDecoResiduals::GetIndex(const std::vector<OBJECT_TYPE*> &vec, const TString &name)
{
  int result = 0;
  for (typename std::vector<OBJECT_TYPE*>::const_iterator iter = vec.begin(), iterEnd = vec.end();
       iter != iterEnd; ++iter, ++result) {
    if (*iter && (*iter)->GetName() == name) return result;
  }
  edm::LogError("Alignment") << "@SUB=PeakDecoResiduals::GetIndex" << " could not find " << name;
  return -1;
}


// constructors and destructor
PeakDecoResiduals::PeakDecoResiduals(const edm::ParameterSet& iConfig)
  : parset_(iConfig), bareTkGeomPtr_(0), 
    debug_(parset_.getParameter<bool>("debug")),
    runOnCosmics_(parset_.getParameter<bool>("runOnCosmics")),
    createTree_(parset_.getParameter<bool>("createTree"))
  
{

  //now do what ever initialization is needed
  edm::Service<TFileService> fs;    
  TFileDirectory filedir = fs->mkdir("PeakDecoResiduals");  
  this->BookHists(filedir);
  if(createTree_){
    
    outTree = fs->make<TTree>("t","Tree");

    outTree->Branch("subdet",     &subdet_,      "subdet/I");
    outTree->Branch("du",         &du_,          "du/F");
    outTree->Branch("dw",         &dw_,          "dw/F");
    outTree->Branch("dtanth",     &dtanth_,      "dtanth/F");
    outTree->Branch("tantrk",     &tantrk_,      "tantrk/F");
    outTree->Branch("tanla",      &tanla_,       "tanla/F");
    outTree->Branch("dttime",     &dttime_,      "dttime/F");
    outTree->Branch("dttimeerr",  &dttimeerr_,   "dttimeerr/F");
    outTree->Branch("ndt",        &ndt_,         "ndt/I");
    outTree->Branch("charge",     &charge_,      "charge/F");
    outTree->Branch("nstrips",    &nstrips_,     "nstrips/I");
    outTree->Branch("u",          &u_,           "u/I");
    outTree->Branch("v",          &v_,           "v/I");
    outTree->Branch("w",          &w_,           "w/I");
    outTree->Branch("x",          &x_,           "x/F");
    outTree->Branch("y",          &y_,           "y/F");
    outTree->Branch("z",          &z_,           "z/F");
    outTree->Branch("r",          &r_,           "r/F");
    outTree->Branch("ring",       &ring_,        "ring/I");
    outTree->Branch("wheel",      &wheel_,       "wheel/I");
    outTree->Branch("petal",      &petal_,       "petal/I");
    outTree->Branch("side",       &side_,        "side/I");
    outTree->Branch("stereo",     &stereo_,      "stereo/I");  
    outTree->Branch("ds",         &ds_,          "ds/I");
    outTree->Branch("zplus",      &zplus_,       "zplus/I");
    outTree->Branch("back",       &back_,        "back/I");
    outTree->Branch("front",      &front_,       "front/I");
    outTree->Branch("layer",      &layer_,       "layer/I");
    outTree->Branch("rod",        &rod_,         "rod/I");
    
  }
}

PeakDecoResiduals::~PeakDecoResiduals()
{
  if(debug_) cout<<"PeakDecoResiduals::~PeakDecoResiduals"<<endl;
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //for( std::vector<TH1*>::const_iterator it = vDeleteObjects_.begin(), itEnd = vDeleteObjects_.end();it != itEnd; ++it) delete *it;
  //for( std::vector<TH2*>::const_iterator itTH2 = vDeleteObjectsTH2_.begin(), itEndTH2 = vDeleteObjectsTH2_.end();itTH2 != itEndTH2; ++itTH2) delete *itTH2;
}





void PeakDecoResiduals::BookHists(TFileDirectory &tfd){

  if(debug_) cout<<"PeakDecoResiduals::bookHists"<<endl;
  
  const char* det[7]  = {"pixba","pixec","TIB","TID","TOB","TECthick","TECthin"};
  const char* v[3]    = {"all","vp","vm"};
  const char* vdir[2] = {"vp","vm"};
  const char* wdir[2] = {"wp","wm"};

  hlumivsrun      = tfd.make<TH2F>("lumivsrun","lumivsrun",1000,123500,124500,1000,0,1000);

  for(int idet=0;idet<7;idet++){

    hdttime[idet]             = tfd.make<TH1F>(Form("dttime_%s",det[idet]),            
					 Form("dttime (%s)",det[idet]),100,-50,50);
    hdttimeerr[idet]          = tfd.make<TH1F>(Form("dttimeerr_%s",det[idet]),         
					 Form("dttimeerr (%s)",det[idet]),50,0,50);
    hndt[idet]                = tfd.make<TH1F>(Form("ndt_%s",det[idet]),               
					 Form("ndt (%s)",det[idet]),100,0,100);
    hduvsdtantheta[idet]      = tfd.make<TH2F>(Form("duvsdtantheta_%s",det[idet]),     
					 Form("duvsdtantheta (%s)",det[idet]),100,-2,2,200,-500,500);

    hresxvsdtantheta[idet]      = tfd.make<TH2F>(Form("resxvsdtantheta_%s",det[idet]),     
						 Form("resxvsdtantheta (%s)",det[idet]),100,-2,2,200,-500,500);


    hduvstrktheta[idet]       = tfd.make<TH2F>(Form("duvstrktheta_%s",det[idet]),      
					 Form("duvstrktheta (%s)",det[idet]),100,-2,2,200,-500,500);
    hchargevsdttime[idet]     = tfd.make<TH2F>(Form("chargevsdttime_%s",det[idet]),    
					 Form("chargevsdttime (%s)",det[idet]), 50, -25, 25, 100, 0, 1000);
    hduvsdttime[idet]         = tfd.make<TH2F>(Form("duvsdttime_%s",det[idet]),        
					 Form("duvsdttime (%s)",det[idet]), 50, -25, 25, 200, -500, 500);
    hdwvsdttime[idet]         = tfd.make<TH2F>(Form("dwvsdttime_%s",det[idet]),        
					 Form("dwvsdttime (%s)",det[idet]), 50, -25, 25, 200, -1000, 1000);
    hhitx[idet]               = tfd.make<TH1F>(Form("hitx_%s",det[idet]),              
					 Form("hitx (%s)",det[idet]),100,-5000,5000);
    hhity[idet]               = tfd.make<TH1F>(Form("hity_%s",det[idet]),              
					 Form("hity (%s)",det[idet]),100,-5000,5000);
    hhitz[idet]               = tfd.make<TH1F>(Form("hitz_%s",det[idet]),   
					 Form("hitz (%s)",det[idet]),100,-5000,5000);
    hcharge[idet]             = tfd.make<TH1F>(Form("charge_%s",det[idet]),   
					 Form("charge (%s)",det[idet]),100,0,1000);
    hnstrips[idet]            = tfd.make<TH1F>(Form("nstrips_%s",det[idet]),   
					 Form("nstrips (%s)",det[idet]),25,0,25);
    huOrientation[idet]       = tfd.make<TH1F>(Form("uOrientation_%s",det[idet]),   
					 Form("uOrientation (%s)",det[idet]),3,1.5,1.5);
    hvOrientation[idet]       = tfd.make<TH1F>(Form("vOrientation_%s",det[idet]),   
					 Form("vOrientation (%s)",det[idet]),3,1.5,1.5);
    hwOrientation[idet]       = tfd.make<TH1F>(Form("wOrientation_%s",det[idet]),   
					 Form("wOrientation (%s)",det[idet]),3,1.5,1.5);
    hlocaltheta[idet]         = tfd.make<TH1F>(Form("localtheta_%s",det[idet]),   
					 Form("localtheta (%s)",det[idet]),50,-2,2);
    hlocalphi[idet]           = tfd.make<TH1F>(Form("localphi_%s",det[idet]),   
					 Form("localphi (%s)",det[idet]),50,-1,1);
    hdu[idet]                 = tfd.make<TH1F>(Form("du_%s",det[idet]),   
					 Form("du (%s)",det[idet]),1000,-500,500);
    hdw[idet]                 = tfd.make<TH1F>(Form("dw_%s",det[idet]),   
					 Form("dw (%s)",det[idet]),1000,-1000,1000);

    hresxoverdtanth[idet]                 = tfd.make<TH1F>(Form("resxoverdtanth_%s",det[idet]),   
							   Form("resxoverdtanth (%s)",det[idet]),1000,-1000,1000);
    
    for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){
	hduvsdtantheta_vw[idet][i][j]   = tfd.make<TH2F>(Form("duvsdtantheta_%s_%s_%s",det[idet],vdir[i],wdir[j]),  
							 Form("duvsdtantheta V+ (%s) %s %s",det[idet],vdir[i],wdir[j]),100,-2,2,200,-500,500);
	hduvstrktheta_vw[idet][i][j]    = tfd.make<TH2F>(Form("duvstrktheta_%s_%s_%s",det[idet],vdir[i],wdir[j]),   
							 Form("duvstrktheta V+ (%s) %s %s",det[idet],vdir[i],wdir[j]),100,-2,2,200,-500,500);
      }
    }


    for(int idt=0;idt<8;idt++){
      hduvsdtantheta_vp_dt[idet][idt]   = tfd.make<TH2F>(Form("duvsdtantheta_vp_dt_%s_%i",det[idet],idt),  
						      Form("duvsdtantheta V+ (%s) (DT slice %i)",det[idet],idt),100,-2,2,200,-500,500);
      hduvsdtantheta_vm_dt[idet][idt]   = tfd.make<TH2F>(Form("duvsdtantheta_vm_dt_%s_%i",det[idet],idt),  
						      Form("duvsdtantheta V- (%s) (DT slice %i)",det[idet],idt),100,-2,2,200,-500,500);  
    }
    
    for(int iv=0;iv<3;iv++){
      htantrktheta[idet][iv]        = tfd.make<TH1F>(Form("tantrktheta_%s_%s",det[idet],v[iv]),   
						     Form("tantrktheta (%s) (%s)",det[idet],v[iv]),50,-2,2);
      htanlorangle[idet][iv]        = tfd.make<TH1F>(Form("tanlorangle_%s_%s",det[idet],v[iv]),   
						     Form("tanlorangle (%s) (%s)",det[idet],v[iv]),50,-2,2);
      hdtantheta[idet][iv]          = tfd.make<TH1F>(Form("dtantheta_%s_%s",det[idet],v[iv]),   
						     Form("dtantheta (%s) (%s)",det[idet],v[iv]),50,-2,2);
    }
  }

if(debug_)cout<<"Exit bookHists"<<endl;
}
				  
int getDTSlice(float t){
  
  int islice = -1;
  if( t > -20. && t < -12. ) islice = 0;
  if( t > -12. && t <  -8. ) islice = 1;
  if( t > -8.  && t <  -4. ) islice = 2;
  if( t > -4.  && t <   0. ) islice = 3;
  if( t >  0.  && t <   4. ) islice = 4;
  if( t >  4.  && t <   8. ) islice = 5;
  if( t >  8.  && t <  12. ) islice = 6;
  if( t > 12.  && t <  20. ) islice = 7;

  return islice;
}


bool PeakDecoResiduals::isBarrel(uint32_t subDetId)
{
  return (subDetId == StripSubdetector::TIB ||
	  subDetId == StripSubdetector::TOB ||
	  subDetId == PixelSubdetector::PixelBarrel );

}

bool PeakDecoResiduals::isEndCap(uint32_t subDetId)
{
  return ( subDetId == StripSubdetector::TID ||
	   subDetId == StripSubdetector::TEC ||
	   subDetId == PixelSubdetector::PixelEndcap);
}

bool PeakDecoResiduals::isPixel(uint32_t subDetId)
{
  return (subDetId == PixelSubdetector::PixelBarrel || subDetId == PixelSubdetector::PixelEndcap);
}


bool PeakDecoResiduals::isDetOrDetUnit(align::StructureType type)
{
  return ( type == align::AlignableDet || type == align::AlignableDetUnit);
}


int PeakDecoResiduals::getSubdetInteger(const DetId& detid){

  // get an integer depending on tracker subdetector
  // if no object exist, the reference is automatically created by the map
  // throw exception if non-tracker id is passed
  uint subdetid = detid.subdetId();

  if     (subdetid  == PixelSubdetector::PixelBarrel)   return 0;
  else if(subdetid  == PixelSubdetector::PixelEndcap)   return 1;
  else if(subdetid  == StripSubdetector::TIB)           return 2;
  else if(subdetid  == StripSubdetector::TID)           return 3;
  else if(subdetid  == StripSubdetector::TOB)           return 4;
  else if(subdetid  == StripSubdetector::TEC){
    TECDetId tecID(detid.rawId());
    int ring  = tecID.ring();
    if     (ring == 5 || ring ==6 || ring ==7)          return 5; //TEC thick
    else if(ring >0 && ring<5)                          return 6; //TEC thin
    else {cout<<"ERROR RINGNUM "<<ring<<endl;}
  }
  else {
    throw cms::Exception("Geometry Error")
      << "[PeakDecoResiduals] Error, tried to get reference for non-tracker subdet " << subdetid 
      << " from detector " << detid.det();
    return -9999;
  }
  
  return -9999;
}


// ------------ method called to for each event  ------------
void
PeakDecoResiduals::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  if(debug_) cout<<"PeakDecoResiduals::analyze"<<endl;
  int evt_run                       = iEvent.id().run()        ;
  int evt_lumiBlock                 = iEvent.luminosityBlock() ;  
  //int evt_event                     = iEvent.id().event()      ;

  //cout<<"PeakDecoResiduals Run: "<<evt_run<<" LumiBlock: "<<evt_lumiBlock<<endl;

  if(debug_)dout<<endl;
  
  hlumivsrun->Fill(evt_run,evt_lumiBlock);
  //using namespace edm;
  TrackerValidationVariables avalidator_(iSetup,parset_);
  edm::Service<TFileService> fs;
  
  std::vector<TrackerValidationVariables::AVHitStruct> v_hitstruct;
  //avalidator_.fillHitQuantities(iEvent,v_hitstruct);
  avalidator_.fillHitQuantities(iEvent,iSetup,v_hitstruct,runOnCosmics_);
  if(debug_)dout<<endl;
  
  // hit quantities: residuals, normalized residuals
  for (std::vector<TrackerValidationVariables::AVHitStruct>::const_iterator it = v_hitstruct.begin(),
  	 itEnd = v_hitstruct.end(); it != itEnd; ++it) {

    DetId detid(it->rawDetId);

    uint subdetid = detid.subdetId();
  
    if(subdetid  == PixelSubdetector::PixelBarrel){
      ring_    = -9999;
      wheel_   = -9999;
      petal_   = -9999;
      side_    = -9999;
      stereo_  = -9999;
      ds_      = -9999;
      zplus_   = -9999;
      back_    = -9999;
      front_   = -9999;
      layer_   = -9999;
      rod_     = -9999;
    }
    
    else if(subdetid  == PixelSubdetector::PixelEndcap){
      ring_    = -9999;
      wheel_   = -9999;
      petal_   = -9999;
      side_    = -9999;
      stereo_  = -9999;
      ds_      = -9999;
      zplus_   = -9999;
      back_    = -9999;
      front_   = -9999;
      layer_   = -9999;
      rod_     = -9999;
    }
    
    else if(subdetid  == StripSubdetector::TIB){
      TIBDetId tibID(detid.rawId());
      ring_    = -9999;
      wheel_   = -9999;
      petal_   = -9999;
      side_    = -9999;
      stereo_  = tibID.isStereo()      ? 1 : 0;
      ds_      = tibID.isDoubleSide()  ? 1 : 0;
      zplus_   = tibID.isZPlusSide()   ? 1 : 0;
      back_    = -9999;
      front_   = -9999;
      layer_   = tibID.layerNumber();
      rod_     = tibID.stringNumber();
    }
    
    else if(subdetid  == StripSubdetector::TOB){
      TOBDetId tobID(detid.rawId());
      ring_    = -9999;
      wheel_   = -9999;
      petal_   = -9999;
      side_    = -9999;
      stereo_  = tobID.isStereo()      ? 1 : 0;
      ds_      = tobID.isDoubleSide()  ? 1 : 0;
      zplus_   = tobID.isZPlusSide()   ? 1 : 0;
      back_    = -9999;
      front_   = -9999;
      layer_   = tobID.layerNumber();
      rod_     = tobID.rodNumber();
    }
    
    else if(subdetid  == StripSubdetector::TID){
      TIDDetId tidID(detid.rawId());
      ring_    = tidID.ringNumber(); 
      wheel_   = tidID.wheel();
      petal_   = -9999;
      side_    = tidID.side();
      stereo_  = tidID.isStereo()      ? 1 : 0;
      ds_      = tidID.isDoubleSide()  ? 1 : 0;
      zplus_   = tidID.isZPlusSide()   ? 1 : 0;
      back_    = tidID.isBackRing()    ? 1 : 0;
      front_   = tidID.isFrontRing()   ? 1 : 0;
      layer_   = -9999;
      rod_     = -9999;
    }
    
    else if(subdetid  == StripSubdetector::TEC){
      TECDetId tecID(detid.rawId());
      ring_    = tecID.ringNumber(); 
      wheel_   = tecID.wheelNumber();
      petal_   = tecID.petalNumber();
      side_    = tecID.side();
      stereo_  = tecID.isStereo()      ? 1 : 0;
      ds_      = tecID.isDoubleSide()  ? 1 : 0;
      zplus_   = tecID.isZPlusSide()   ? 1 : 0;
      back_    = tecID.isBackPetal()   ? 1 : 0;
      front_   = tecID.isFrontPetal()  ? 1 : 0;
      layer_   = -9999;
      rod_     = -9999;
    }

    else {
      throw cms::Exception("Geometry Error")
	<< "[PeakDecoResiduals] Error, tried to get reference for non-tracker subdet " << subdetid 
	<< " from detector " << detid.det();
    }
  
    int subdetint = this->getSubdetInteger(detid);

    if(it->nvalidmu==1){
      if(it->dttimeerr != -999)   hdttimeerr[subdetint]->Fill(it->dttimeerr);
      if(it->ndt       != -999)   hndt[subdetint]      ->Fill(it->ndt);
    }

    //for cosmics running, require single muon and cut on nDT hits and DT time error
    if(runOnCosmics_){
      if(it->nvalidmu != 1) continue;
      if(it->ndt < 25)      continue;
      if(it->dttimeerr >10) continue;
    }
   
    //if(it->nstrips == 1 || it->nstrips == 2)    continue;


    if(it->resXprime != -999. && it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.) {
      float dtantheta = it->tanTrackAngle-it->tanLorentzAngle;

      hduvsdtantheta[subdetint]   -> Fill(dtantheta,         10000*it->resXprime);
      hduvstrktheta[subdetint]    -> Fill(it->tanTrackAngle, 10000*it->resXprime);
      hresxvsdtantheta[subdetint] -> Fill(dtantheta,         10000*it->resX);

      int i = 0;
      int j = 0;
      if(it->vOrientation < 0) i = 1;
      if(it->wOrientation < 0) j = 1;
      
      hduvsdtantheta_vw[subdetint][i][j] -> Fill(dtantheta,         10000*it->resXprime);
      hduvstrktheta_vw[subdetint][i][j]  -> Fill(it->tanTrackAngle, 10000*it->resXprime);
    }

    if(runOnCosmics_){
      if(it->charge != -999. && it->dttime != -999.){
	hchargevsdttime[subdetint]->Fill(it->dttime,it->charge);
      }
      
      if(it->resXprime != -999. && it->dttime != -999. && it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.){
	float dtantheta = it->tanTrackAngle-it->tanLorentzAngle;
	if(dtantheta != 0) hdwvsdttime[subdetint]->Fill(it->dttime,  10000*it->resXprime/dtantheta);
	if(dtantheta > 0){
	  hduvsdttime[subdetint]->Fill(it->dttime,                   10000*it->resXprime);    
	}
	if(dtantheta < 0){
	  hduvsdttime[subdetint]->Fill(it->dttime,                 - 10000*it->resXprime);    
	}
      }
     
      if(getDTSlice(it->dttime) > -1 && it->resXprime != -999. && it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.){
	float dtantheta = it->tanTrackAngle-it->tanLorentzAngle;
	if(it->vOrientation < 0)	hduvsdtantheta_vm_dt[subdetint][getDTSlice(it->dttime)] -> Fill(dtantheta, 10000*it->resXprime);
	if(it->vOrientation > 0)	hduvsdtantheta_vp_dt[subdetint][getDTSlice(it->dttime)] -> Fill(dtantheta, 10000*it->resXprime);
      }
    }
    
    if(it->hitx != -999.) hhitx[subdetint]->Fill(it->hitx);
    if(it->hity != -999.) hhity[subdetint]->Fill(it->hity);
    if(it->hitz != -999.) hhitz[subdetint]->Fill(it->hitz);

    if(it->dttime != -999)  hdttime[subdetint]->  Fill(it->dttime);
    if(it->charge != -999.) hcharge[subdetint]->  Fill(it->charge);
    if(it->nstrips != -999) hnstrips[subdetint]-> Fill(it->nstrips);

    if(it->uOrientation != -999.)   huOrientation[subdetint]->Fill(it->uOrientation);
    if(it->vOrientation != -999.)   hvOrientation[subdetint]->Fill(it->vOrientation);
    if(it->wOrientation != -999.)   hwOrientation[subdetint]->Fill(it->wOrientation);
    if(it->localtheta != -999.)     hlocaltheta[subdetint]  ->Fill(tan(it->localtheta));
    if(it->localphi != -999.)       hlocalphi[subdetint]    ->Fill(cos(it->localphi));

    if(it->resXprime != -999. && it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.) {
      float dtantheta = it->tanTrackAngle-it->tanLorentzAngle;
      if(dtantheta != 0){
	hdw[subdetint]              ->Fill( 10000*it->resXprime/dtantheta);
	hresxoverdtanth[subdetint]  ->Fill( 10000*it->resX/dtantheta);
      }
      if(dtantheta > 0)  hdu[subdetint]->Fill( 10000*it->resXprime);
      if(dtantheta < 0)  hdu[subdetint]->Fill(-10000*it->resXprime);

    }
    int iv=1;
    if(it->vOrientation < 0) iv=2;
    
    if(it->tanTrackAngle != -999.)	                               htantrktheta[subdetint][iv]->Fill(it->tanTrackAngle);
    if(it->tanLorentzAngle != -999.)	                               htanlorangle[subdetint][iv]->Fill(it->tanLorentzAngle);
    if(it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.)     hdtantheta[subdetint][iv]  ->Fill(it->tanTrackAngle-it->tanLorentzAngle);
    
    if(it->tanTrackAngle != -999.)	                               htantrktheta[subdetint][0]->Fill(it->tanTrackAngle);
    if(it->tanLorentzAngle != -999.)	                               htanlorangle[subdetint][0]->Fill(it->tanLorentzAngle);
    if(it->tanTrackAngle != -999. && it->tanLorentzAngle != -999.)     hdtantheta[subdetint][0]  ->Fill(it->tanTrackAngle-it->tanLorentzAngle);
   
           
    if(createTree_){
      subdet_      = subdetint;
      du_          = 10000*it->resXprime;
      tantrk_      = it->tanTrackAngle;
      tanla_       = it->tanLorentzAngle;
      dtanth_      = it->tanTrackAngle-it->tanLorentzAngle;
      dw_          = du_/dtanth_;
      dttime_      = it->dttime;
      dttimeerr_   = it->dttimeerr;
      ndt_         = it->ndt;
      charge_      = it->charge;
      nstrips_     = it->nstrips;
      u_           = (int)it->uOrientation;
      v_           = (int)it->vOrientation;
      w_           = (int)it->wOrientation;
      x_           = it->hitx;
      y_           = it->hity;
      z_           = it->hitz;
      r_           = sqrt(pow(it->hitx,2) + pow(it->hity,2));
      outTree->Fill();
    }
  }
  if(debug_) cout<<"End PeakDecoResiduals::analyze"<<endl;

}



// ------------ method called once each job just after ending the event loop  ------------
void 
PeakDecoResiduals::endJob()
{

  if(debug_) cout<<"PeakDecoResiduals::endJob"<<endl;

  //if(createTree_){
  //outTree->Write();
  //}
   
}

//define this as a plug-in
DEFINE_FWK_MODULE(PeakDecoResiduals);
