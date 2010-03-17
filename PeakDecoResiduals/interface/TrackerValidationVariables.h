
#ifndef TrackerTrackerValidationVariables_h
#define TrackerTrackerValidationVariables_h
// system include files
#include <memory>
#include <vector>
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"

class MagneticField;
class TrackerGeometry;
class SiStripLorentzAngle;

class TrackerValidationVariables  {
 public:  
  struct AVHitStruct{
    AVHitStruct() : resX(-999.), resY(-999.), resErrX(-999.), resErrY(-999.), resXprime(-999.), resXprimeErr(-999.), 
	 resYprime(-999.), resYprimeErr(-999.), phi(-999.), eta(-999.), tanTrackAngle(-999.), tanLorentzAngle(-999.), 
	 uOrientation(-999.), vOrientation(-999.), wOrientation(-999.), dusign(-999.),
	 x(-999.), y(-999.), z(-999.), localtheta(-999.), localphi(-999.), dttime(-999.), dttimeerr(-999.), 
	 ndt(-999), ectype(-999), hcaltime(-999.), hcaltimeerr(-999.), nvalidmu(-999), charge(-999.), nstrips(-999), p(-999.),
	 rawDetId(0), overlapres(std::make_pair(0,-999.)) {}
    float resX;
    float resY;
    float resErrX;
    float resErrY;
    float resXprime;
    float resXprimeErr;
    float resYprime;
    float resYprimeErr;
    float phi;
    float eta;
    float tanTrackAngle;
    float tanLorentzAngle;
    float uOrientation;
    float vOrientation;
    float wOrientation;
    float dusign;
    float x;
    float y;
    float z;
    float localtheta;
    float localphi;
    float dttime;
    float dttimeerr;
    int ndt;
    int ectype;
    float hcaltime;
    float hcaltimeerr;
    int nvalidmu;
    float charge;
    int nstrips;
    float p;
    uint32_t rawDetId;
    std::pair<uint,float> overlapres;
  };
  struct AVTrackStruct{
  AVTrackStruct() : pt(0.), ptError(0.), px(0.), py(0.), pz(0.), eta(0.), phi(0.), kappa(0.),
	 chi2(0.), normchi2(0), d0(-999.), dz(-999.), charge(-999) {};
    float pt;
    float ptError;
    float px;
    float py;
    float pz;
    float eta;
    float phi;
    float kappa;
    float chi2;
    float normchi2;
    float d0;
    float dz;
    int charge;
  };
  struct timeStruct{
    timeStruct() :  nmu(-999), dttime(-999.), dttimeerr(-999.), ndt(-999), 
	 ecaltime(-999.), ecaltimeerr(-999.), ecalenergy(-999.),
	 hcaltime(-999.), hcaltimeerr(-999.), hcalenergy(-999.), nvalidmu(-999), p(-999.) {};

    int nmu;
    float dttime;
    float dttimeerr;
    int ndt;
    float ecaltime;
    float ecaltimeerr;
    float ecalenergy;
    float hcaltime;
    float hcaltimeerr;
    float hcalenergy;
    int nvalidmu;
    float p;
  };

  TrackerValidationVariables();
  TrackerValidationVariables(const edm::EventSetup&, const edm::ParameterSet&);
  ~TrackerValidationVariables();
  //void fillHitQuantities(const edm::Event&, std::vector<AVHitStruct> & v_avhitout );
  void fillHitQuantities(const edm::Event&, const edm::EventSetup&, std::vector<AVHitStruct> & v_avhitout, bool runOnCosmics_ );
  int getHisto(float z);
  void setTime(const edm::Event&, timeStruct &ts);
  //float getHCALTime(const edm::Event&);
  void fillTrackQuantities(const edm::Event&, std::vector<AVTrackStruct> & v_avtrackout );

 private:  
  const edm::ParameterSet conf_;
  edm::ESHandle<TrackerGeometry> tkGeom_;
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<SiStripLorentzAngle> SiStripLorentzAngle_;
  //edm::ESHandle<SiStripDetCabling> SiStripDetCabling_;

};
#endif
