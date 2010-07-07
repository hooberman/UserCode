// -*- C++ -*-
//
// Package:    PFRecHitAnalyzer
// Class:      PFRecHitAnalyzer
// 
/**\class PFRecHitAnalyzer PFRecHitAnalyzer.cc PFstudies/PFRecHitAnalyzer/src/PFRecHitAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu Jun 24 02:52:45 PDT 2010
// $Id: PFRecHitAnalyzer.cc,v 1.2 2010/06/30 16:11:51 benhoob Exp $
//
//

// include files
#include "PFstudies/PFRecHitAnalyzer/interface/PFRecHitAnalyzer.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

using namespace edm;
using namespace std;

//
// constants, enums and typedefs
//
typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;

//
// static data member definitions
//

//
// constructors and destructor
//
PFRecHitAnalyzer::PFRecHitAnalyzer(const edm::ParameterSet& iConfig)
{
  // input tags
  pf_ecal_tag_         	= iConfig.getParameter<edm::InputTag>("pf_ecal_tag"	);
  pf_hcal_tag_	        = iConfig.getParameter<edm::InputTag>("pf_hcal_tag"	);
  pf_hfem_tag_	        = iConfig.getParameter<edm::InputTag>("pf_hfem_tag"	);
  pf_hfhad_tag_	        = iConfig.getParameter<edm::InputTag>("pf_hfhad_tag"	);
  calomet_tag_	        = iConfig.getParameter<edm::InputTag>("calomet_tag"	);
  genmet_tag_	        = iConfig.getParameter<edm::InputTag>("genmet_tag"	);
  inputTagPFClustersECAL_  = iConfig.getParameter<edm::InputTag>("PFClustersECAL");
  inputTagPFClustersHCAL_  = iConfig.getParameter<edm::InputTag>("PFClustersHCAL");
  inputTagPFClustersHFEM_  = iConfig.getParameter<edm::InputTag>("PFClustersHFEM");
  inputTagPFClustersHFHAD_ = iConfig.getParameter<edm::InputTag>("PFClustersHFHAD");

  eb_threshold_	        = iConfig.getParameter<double>("eb_threshold"		);
  ee_threshold_	        = iConfig.getParameter<double>("ee_threshold"		);
  hb_threshold_	        = iConfig.getParameter<double>("hb_threshold"		);
  he_threshold_	        = iConfig.getParameter<double>("he_threshold"		);
  hfe_threshold_	= iConfig.getParameter<double>("hfe_threshold"		);
  hfh_threshold_	= iConfig.getParameter<double>("hfh_threshold"		);

  eb_cluster_threshold_	        = iConfig.getParameter<double>("eb_cluster_threshold"		);
  ee_cluster_threshold_	        = iConfig.getParameter<double>("ee_cluster_threshold"		);
  hb_cluster_threshold_	        = iConfig.getParameter<double>("hb_cluster_threshold"		);
  he_cluster_threshold_	        = iConfig.getParameter<double>("he_cluster_threshold"		);
  hfe_cluster_threshold_	= iConfig.getParameter<double>("hfe_cluster_threshold"		);
  hfh_cluster_threshold_	= iConfig.getParameter<double>("hfh_cluster_threshold"		);

  useClusters_         	= iConfig.getParameter<bool>("useClustersForMET");

  // register products
  produces<float>	("caloebsumet" ).setBranchAlias("calo_eb_sumet"		);
  produces<float>	("caloeesumet" ).setBranchAlias("calo_ee_sumet"		);
  produces<float>	("calohbsumet" ).setBranchAlias("calo_hb_sumet"		);
  produces<float>	("calohesumet" ).setBranchAlias("calo_he_sumet"		);
  produces<float>	("calohfesumet").setBranchAlias("calo_hfe_sumet"	);
  produces<float>	("calohfhsumet").setBranchAlias("calo_hfh_sumet"	);

  produces<float>	("pfebsumet").setBranchAlias("pf_eb_sumet"	);
  produces<float>	("pfeesumet").setBranchAlias("pf_ee_sumet"	);
  produces<float>	("pfhbsumet").setBranchAlias("pf_hb_sumet"	);
  produces<float>	("pfhesumet").setBranchAlias("pf_he_sumet"	);
  produces<float>	("pfhfesumet").setBranchAlias("pf_hfe_sumet"	);
  produces<float>	("pfhfhsumet").setBranchAlias("pf_hfh_sumet"	);

  produces<float>	("pfebmet").setBranchAlias("pf_eb_met"	);
  produces<float>	("pfeemet").setBranchAlias("pf_ee_met"	);
  produces<float>	("pfhbmet").setBranchAlias("pf_hb_met"	);
  produces<float>	("pfhemet").setBranchAlias("pf_he_met"	);
  produces<float>	("pfhfemet").setBranchAlias("pf_hfe_met"	);
  produces<float>	("pfhfhmet").setBranchAlias("pf_hfh_met"	);

  produces<float>	("pfclusebsumet").setBranchAlias("pfclus_eb_sumet"	);
  produces<float>	("pfcluseesumet").setBranchAlias("pfclus_ee_sumet"	);
  produces<float>	("pfclushbsumet").setBranchAlias("pfclus_hb_sumet"	);
  produces<float>	("pfclushesumet").setBranchAlias("pfclus_he_sumet"	);
  produces<float>	("pfclushfesumet").setBranchAlias("pfclus_hfe_sumet"	);
  produces<float>	("pfclushfhsumet").setBranchAlias("pfclus_hfh_sumet"	);

  produces<float>	("pfclusebmet").setBranchAlias("pfclus_eb_met"	);
  produces<float>	("pfcluseemet").setBranchAlias("pfclus_ee_met"	);
  produces<float>	("pfclushbmet").setBranchAlias("pfclus_hb_met"	);
  produces<float>	("pfclushemet").setBranchAlias("pfclus_he_met"	);
  produces<float>	("pfclushfemet").setBranchAlias("pfclus_hfe_met"	);
  produces<float>	("pfclushfhmet").setBranchAlias("pfclus_hfh_met"	);

  produces<float>	("calomet").setBranchAlias("calomet"		);
  produces<float>	("calosumet").setBranchAlias("calosumet"		);

  produces<float>	("pfmet").setBranchAlias("pfmet"		);
  produces<float>	("pfsumet").setBranchAlias("pfsumet"		);

  produces<float>	("pfclusmet").setBranchAlias("pfclusmet"		);
  produces<float>	("pfclussumet").setBranchAlias("pfclussumet"		);


  produces<float>    ("genmet").setBranchAlias("genmet"		);
  produces<float>    ("gensumet").setBranchAlias("gensumet"		);
  
  produces<std::vector<float> > ("pfrechitet"	).setBranchAlias("pf_rechit_et"		);
  produces<std::vector<float> > ("pfrechite"	).setBranchAlias("pf_rechit_e"		);
  produces<std::vector<float> > ("pfrechiteta"	).setBranchAlias("pf_rechit_eta"	);
  produces<std::vector<float> > ("pfrechitphi"	).setBranchAlias("pf_rechit_phi"	);
  produces<std::vector<float> > ("pfrechitdetid").setBranchAlias("pf_rechit_detid"      );

  produces<std::vector<float> > ("pfebrechitet"	).setBranchAlias("pf_ebrechit_et"	);
  produces<std::vector<float> > ("pfebrechite"	).setBranchAlias("pf_ebrechit_e"	);
  produces<std::vector<float> > ("pfebrechiteta").setBranchAlias("pf_ebrechit_eta"	);
  produces<std::vector<float> > ("pfebrechitphi").setBranchAlias("pf_ebrechit_phi"	);

  produces<std::vector<float> > ("pfeerechitet"	).setBranchAlias("pf_eerechit_et"	);
  produces<std::vector<float> > ("pfeerechite"	).setBranchAlias("pf_eerechit_e"	);
  produces<std::vector<float> > ("pfeerechiteta").setBranchAlias("pf_eerechit_eta"	);
  produces<std::vector<float> > ("pfeerechitphi").setBranchAlias("pf_eerechit_phi"	);

  produces<std::vector<float> > ("pfhbrechitet"	).setBranchAlias("pf_hbrechit_et"	);
  produces<std::vector<float> > ("pfhbrechite"	).setBranchAlias("pf_hbrechit_e"	);
  produces<std::vector<float> > ("pfhbrechiteta").setBranchAlias("pf_hbrechit_eta"	);
  produces<std::vector<float> > ("pfhbrechitphi").setBranchAlias("pf_hbrechit_phi"	);

  produces<std::vector<float> > ("pfherechitet"	).setBranchAlias("pf_herechit_et"	);
  produces<std::vector<float> > ("pfherechite"	).setBranchAlias("pf_herechit_e"	);
  produces<std::vector<float> > ("pfherechiteta").setBranchAlias("pf_herechit_eta"	);
  produces<std::vector<float> > ("pfherechitphi").setBranchAlias("pf_herechit_phi"	);

  produces<std::vector<float> > ("pfhferechitet" ).setBranchAlias("pf_hferechit_et"	);
  produces<std::vector<float> > ("pfhferechite"  ).setBranchAlias("pf_hferechit_e"	);
  produces<std::vector<float> > ("pfhferechiteta").setBranchAlias("pf_hferechit_eta"	);
  produces<std::vector<float> > ("pfhferechitphi").setBranchAlias("pf_hferechit_phi"	);

  produces<std::vector<float> > ("pfhfhrechitet" ).setBranchAlias("pf_hfhrechit_et"	);
  produces<std::vector<float> > ("pfhfhrechite" ).setBranchAlias("pf_hfhrechit_e"	);
  produces<std::vector<float> > ("pfhfhrechiteta").setBranchAlias("pf_hfhrechit_eta"	);
  produces<std::vector<float> > ("pfhfhrechitphi").setBranchAlias("pf_hfhrechit_phi"	);

  produces<std::vector<float> > ("pfclusteret"	).setBranchAlias("pf_cluster_et"	);
  produces<std::vector<float> > ("pfclustere"	).setBranchAlias("pf_cluster_e"		);
  produces<std::vector<float> > ("pfclustereta"	).setBranchAlias("pf_cluster_eta"	);
  produces<std::vector<float> > ("pfclusterphi"	).setBranchAlias("pf_cluster_phi"	);
  produces<std::vector<float> > ("pfclusterdetid").setBranchAlias("pf_cluster_detid"    );

  produces<std::vector<float> > ("pfebclusteret").setBranchAlias("pf_ebcluster_et"	);
  produces<std::vector<float> > ("pfebclustere"	).setBranchAlias("pf_ebcluster_e"	);
  produces<std::vector<float> > ("pfebclustereta").setBranchAlias("pf_ebcluster_eta"	);
  produces<std::vector<float> > ("pfebclusterphi").setBranchAlias("pf_ebcluster_phi"	);

  produces<std::vector<float> > ("pfeeclusteret").setBranchAlias("pf_eecluster_et"	);
  produces<std::vector<float> > ("pfeeclustere"	).setBranchAlias("pf_eecluster_e"	);
  produces<std::vector<float> > ("pfeeclustereta").setBranchAlias("pf_eecluster_eta"	);
  produces<std::vector<float> > ("pfeeclusterphi").setBranchAlias("pf_eecluster_phi"	);

  produces<std::vector<float> > ("pfhbclusteret").setBranchAlias("pf_hbcluster_et"	);
  produces<std::vector<float> > ("pfhbclustere"	).setBranchAlias("pf_hbcluster_e"	);
  produces<std::vector<float> > ("pfhbclustereta").setBranchAlias("pf_hbcluster_eta"	);
  produces<std::vector<float> > ("pfhbclusterphi").setBranchAlias("pf_hbcluster_phi"	);

  produces<std::vector<float> > ("pfheclusteret").setBranchAlias("pf_hecluster_et"	);
  produces<std::vector<float> > ("pfheclustere"	).setBranchAlias("pf_hecluster_e"	);
  produces<std::vector<float> > ("pfheclustereta").setBranchAlias("pf_hecluster_eta"	);
  produces<std::vector<float> > ("pfheclusterphi").setBranchAlias("pf_hecluster_phi"	);

  produces<std::vector<float> > ("pfhfeclusteret" ).setBranchAlias("pf_hfecluster_et"	);
  produces<std::vector<float> > ("pfhfeclustere"  ).setBranchAlias("pf_hfecluster_e"	);
  produces<std::vector<float> > ("pfhfeclustereta").setBranchAlias("pf_hfecluster_eta"	);
  produces<std::vector<float> > ("pfhfeclusterphi").setBranchAlias("pf_hfecluster_phi"	);

  produces<std::vector<float> > ("pfhfhclusteret" ).setBranchAlias("pf_hfhcluster_et"	);
  produces<std::vector<float> > ("pfhfhclustere" ).setBranchAlias("pf_hfhcluster_e"	);
  produces<std::vector<float> > ("pfhfhclustereta").setBranchAlias("pf_hfhcluster_eta"	);
  produces<std::vector<float> > ("pfhfhclusterphi").setBranchAlias("pf_hfhcluster_phi"	);
}


PFRecHitAnalyzer::~PFRecHitAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PFRecHitAnalyzer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<float> calo_eb_sumet		(new float);
  std::auto_ptr<float> calo_ee_sumet		(new float);
  std::auto_ptr<float> calo_hb_sumet		(new float);
  std::auto_ptr<float> calo_he_sumet		(new float);
  std::auto_ptr<float> calo_hfe_sumet	(new float);
  std::auto_ptr<float> calo_hfh_sumet	(new float);
  std::auto_ptr<float> pf_eb_sumet		(new float);
  std::auto_ptr<float> pf_ee_sumet		(new float);
  std::auto_ptr<float> pf_hb_sumet		(new float);
  std::auto_ptr<float> pf_he_sumet		(new float);
  std::auto_ptr<float> pf_hfe_sumet		(new float);
  std::auto_ptr<float> pf_hfh_sumet		(new float);
  std::auto_ptr<float> calomet		(new float);
  std::auto_ptr<float> calosumet		(new float);
  std::auto_ptr<float> pfmet			(new float);
  std::auto_ptr<float> pfsumet		(new float);
  std::auto_ptr<float> genmet		(new float);
  std::auto_ptr<float> gensumet		(new float);
  std::auto_ptr<float> pf_eb_met		(new float);
  std::auto_ptr<float> pf_ee_met		(new float);
  std::auto_ptr<float> pf_hb_met		(new float);
  std::auto_ptr<float> pf_he_met		(new float);
  std::auto_ptr<float> pf_hfe_met		(new float);
  std::auto_ptr<float> pf_hfh_met		(new float);
  std::auto_ptr<float> pfclus_eb_met		(new float);
  std::auto_ptr<float> pfclus_ee_met		(new float);
  std::auto_ptr<float> pfclus_hb_met		(new float);
  std::auto_ptr<float> pfclus_he_met		(new float);
  std::auto_ptr<float> pfclus_hfe_met		(new float);
  std::auto_ptr<float> pfclus_hfh_met		(new float);
  std::auto_ptr<float> pfclus_eb_sumet		(new float);
  std::auto_ptr<float> pfclus_ee_sumet		(new float);
  std::auto_ptr<float> pfclus_hb_sumet		(new float);
  std::auto_ptr<float> pfclus_he_sumet		(new float);
  std::auto_ptr<float> pfclus_hfe_sumet		(new float);
  std::auto_ptr<float> pfclus_hfh_sumet		(new float);
  std::auto_ptr<float> pfclusmet		(new float);
  std::auto_ptr<float> pfclussumet		(new float);

  //rechit quantities
  std::auto_ptr<std::vector<float> > pf_rechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_rechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_rechit_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_rechit_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<int>   > pf_rechit_detid		(new std::vector<int>  );
  std::auto_ptr<std::vector<float> > pf_ebrechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_ebrechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_ebrechit_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_ebrechit_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eerechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eerechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eerechit_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eerechit_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbrechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbrechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbrechit_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbrechit_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_herechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_herechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_herechit_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_herechit_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hferechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hferechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hferechit_eta	(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hferechit_phi	(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhrechit_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhrechit_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhrechit_eta	(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhrechit_phi	(new std::vector<float>);

  //pfcluster quantities
  std::auto_ptr<std::vector<float> > pf_cluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_cluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_cluster_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_cluster_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<int>   > pf_cluster_detid		(new std::vector<int>  );
  std::auto_ptr<std::vector<float> > pf_ebcluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_ebcluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_ebcluster_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_ebcluster_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eecluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eecluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eecluster_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_eecluster_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbcluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbcluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbcluster_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hbcluster_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hecluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hecluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hecluster_eta		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hecluster_phi		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfecluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfecluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfecluster_eta	(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfecluster_phi	(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhcluster_et		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhcluster_e		(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhcluster_eta	(new std::vector<float>);
  std::auto_ptr<std::vector<float> > pf_hfhcluster_phi	(new std::vector<float>);

  // get collections
  edm::Handle<reco::PFRecHitCollection> pf_ecal_h;
  iEvent.getByLabel(pf_ecal_tag_, pf_ecal_h);

  edm::Handle<reco::PFRecHitCollection> pf_hcal_h;
  iEvent.getByLabel(pf_hcal_tag_, pf_hcal_h);

  edm::Handle<reco::PFRecHitCollection> pf_hfem_h;
  iEvent.getByLabel(pf_hfem_tag_, pf_hfem_h);

  edm::Handle<reco::PFRecHitCollection> pf_hfhad_h;
  iEvent.getByLabel(pf_hfhad_tag_, pf_hfhad_h);

  edm::Handle<reco::CaloMETCollection> calomet_h;
  iEvent.getByLabel(calomet_tag_, calomet_h);
  reco::CaloMET metobj = calomet_h->front();

  edm::Handle<reco::GenMETCollection> genmet_h;
  iEvent.getByLabel(genmet_tag_, genmet_h);
  reco::GenMET genmetobj = genmet_h->front();

  bool found;

  //ECAL PFClusters
  edm::Handle< reco::PFClusterCollection > clustersECAL;
  found = iEvent.getByLabel(inputTagPFClustersECAL_, clustersECAL);      

  if(!found )
    LogError("PFBlockProducer")<<" cannot get ECAL clusters: "
                               <<inputTagPFClustersECAL_<<endl;

  //HCAL PFClusters
  edm::Handle< reco::PFClusterCollection > clustersHCAL;
  found = iEvent.getByLabel(inputTagPFClustersHCAL_, clustersHCAL);
  
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HCAL clusters: "
                               <<inputTagPFClustersHCAL_<<endl;

  //HFEM PFClusters
  edm::Handle< reco::PFClusterCollection > clustersHFEM;
  found = iEvent.getByLabel(inputTagPFClustersHFEM_, clustersHFEM);      
  
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HFEM clusters: "
                               <<inputTagPFClustersHFEM_<<endl;
  
  //HFHAD PFClusters
  edm::Handle< reco::PFClusterCollection > clustersHFHAD;
  found = iEvent.getByLabel(inputTagPFClustersHFHAD_, clustersHFHAD);      
  
  if(!found )
    LogError("PFBlockProducer")<<" cannot get HFHAD clusters: "
                               <<inputTagPFClustersHFHAD_<<endl;

  //http://cmslxr.fnal.gov/lxr/source/RecoParticleFlow/PFProducer/plugins/PFBlockProducer.cc
  //http://cmslxr.fnal.gov/lxr/source/RecoParticleFlow/PFProducer/src/PFBlockAlgo.cc

  // store genmet quantities
  *genmet   = genmetobj.pt();
  *gensumet = genmetobj.sumEt();

  // store calomet based quantities
  *calomet		= metobj.pt();
  *calosumet		= metobj.sumEt();
  *calo_eb_sumet	= metobj.emEtInEB();
  *calo_ee_sumet	= metobj.emEtInEE();
  *calo_hb_sumet	= metobj.hadEtInHB();
  *calo_he_sumet	= metobj.hadEtInHE();
  *calo_hfe_sumet	= metobj.emEtInHF();
  *calo_hfh_sumet	= metobj.hadEtInHF();


  // calculate MET and sumET using pf clusters
  float cl_met_x = 0.;
  float cl_met_y = 0.;

  // get EM energy in HF using pf clusters
  float cl_hfem_et   = 0.;
  float cl_hfemet_x  = 0.;
  float cl_hfemet_y  = 0.;

  // get HAD energy in HF using pf clusters
  float cl_hfhad_et  = 0.;
  float cl_hfhmet_x  = 0.;
  float cl_hfhmet_y  = 0.;

  // get energy HCAL barrel using pf clusters
  float cl_hb_et	= 0.;
  float cl_hbmet_x	= 0.;
  float cl_hbmet_y	= 0.;

  // get energy ECAL endcap using pf clusters
  float cl_he_et	= 0.;
  float cl_hemet_x	= 0.;
  float cl_hemet_y	= 0.;

  // get energy ECAL barrel using pf clusters
  float cl_eb_et	= 0.;
  float cl_ebmet_x	= 0.;
  float cl_ebmet_y	= 0.;

  // get energy ECAL endcap using pf clusters
  float cl_ee_et	= 0.;
  float cl_eemet_x	= 0.;
  float cl_eemet_y	= 0.;


  //---------------ECAL Barrel PFClusters---------------
  int i = 0;

  for (reco::PFClusterCollection::const_iterator it = clustersECAL->begin(); it != clustersECAL->end(); it++){

    if ( it->layer() != PFLayer::ECAL_BARREL ) continue;
    if ( it->energy() < eb_cluster_threshold_) continue;

    const REPPoint& cluster_pos = it->positionREP();	  
    double et = it->energy() / cosh( cluster_pos.eta() ); 
    
    //reco::PFCluster cl = clustersECAL->at(i);
    //cout << it << endl;
    //reco::PFCluster cl = it;
    printCluster( it );
    
    pf_ebcluster_et	->push_back(et			);
    pf_ebcluster_e	->push_back(it->energy()	);
    pf_ebcluster_eta	->push_back(cluster_pos.eta()	);
    pf_ebcluster_phi	->push_back(cluster_pos.phi()	);

    pf_cluster_et	->push_back(et			);
    pf_cluster_e	->push_back(it->energy()	);
    pf_cluster_eta	->push_back(cluster_pos.eta()	);
    pf_cluster_phi	->push_back(cluster_pos.phi()	);
    
    cl_met_x    -= et * cos(cluster_pos.phi());
    cl_met_y    -= et * sin(cluster_pos.phi());
    cl_ebmet_x  -= et * cos(cluster_pos.phi());
    cl_ebmet_y  -= et * sin(cluster_pos.phi());
    cl_eb_et    += et;      

    ++i;
  }

  //---------------ECAL Endcap PFClusters---------------

  for (reco::PFClusterCollection::const_iterator it = clustersECAL->begin(); it != clustersECAL->end(); it++){

    if ( it->layer() != PFLayer::ECAL_ENDCAP ) continue;
    if ( it->energy() < ee_cluster_threshold_) continue;

    const REPPoint cluster_pos = it->positionREP();	  
    double et = it->energy() / cosh( cluster_pos.eta() ); 
    
    //printCluster( it );

    pf_eecluster_et	->push_back(et			);
    pf_eecluster_e	->push_back(it->energy()	);
    pf_eecluster_eta	->push_back(cluster_pos.eta()	);
    pf_eecluster_phi	->push_back(cluster_pos.phi()	);

    pf_cluster_et	->push_back(et			);
    pf_cluster_e	->push_back(it->energy()	);
    pf_cluster_eta	->push_back(cluster_pos.eta()	);
    pf_cluster_phi	->push_back(cluster_pos.phi()	);
    
    cl_met_x    -= et * cos(cluster_pos.phi());
    cl_met_y    -= et * sin(cluster_pos.phi());
    cl_eemet_x  -= et * cos(cluster_pos.phi());
    cl_eemet_y  -= et * sin(cluster_pos.phi());
    cl_ee_et    += et;      
  }

  //---------------HCAL Barrel PFClusters---------------

  for (reco::PFClusterCollection::const_iterator it = clustersHCAL->begin(); it != clustersHCAL->end(); it++){

    if ( it->layer() != PFLayer::HCAL_BARREL1 ) continue;
    if ( it->energy() < hb_cluster_threshold_)  continue;

    const REPPoint cluster_pos = it->positionREP();	  
    double et = it->energy() / cosh( cluster_pos.eta() ); 
    
    //printCluster( it );

    pf_hbcluster_et	->push_back(et			);
    pf_hbcluster_e	->push_back(it->energy()	);
    pf_hbcluster_eta	->push_back(cluster_pos.eta()	);
    pf_hbcluster_phi	->push_back(cluster_pos.phi()	);

    pf_cluster_et	->push_back(et			);
    pf_cluster_e	->push_back(it->energy()	);
    pf_cluster_eta	->push_back(cluster_pos.eta()	);
    pf_cluster_phi	->push_back(cluster_pos.phi()	);
    
    cl_met_x    -= et * cos(cluster_pos.phi());
    cl_met_y    -= et * sin(cluster_pos.phi());
    cl_hbmet_x  -= et * cos(cluster_pos.phi());
    cl_hbmet_y  -= et * sin(cluster_pos.phi());
    cl_hb_et    += et;      
  }

  //---------------HCAL Endcap PFClusters---------------

  for (reco::PFClusterCollection::const_iterator it = clustersHCAL->begin(); it != clustersHCAL->end(); it++){

    if ( it->layer() != PFLayer::HCAL_ENDCAP ) continue;
    if ( it->energy() < he_cluster_threshold_) continue;

    const REPPoint cluster_pos = it->positionREP();	  
    double et = it->energy() / cosh( cluster_pos.eta() ); 

    //printCluster( it );

    pf_hecluster_et	->push_back(et			);
    pf_hecluster_e	->push_back(it->energy()	);
    pf_hecluster_eta	->push_back(cluster_pos.eta()	);
    pf_hecluster_phi	->push_back(cluster_pos.phi()	);

    pf_cluster_et	->push_back(et			);
    pf_cluster_e	->push_back(it->energy()	);
    pf_cluster_eta	->push_back(cluster_pos.eta()	);
    pf_cluster_phi	->push_back(cluster_pos.phi()	);
    
    cl_met_x    -= et * cos(cluster_pos.phi());
    cl_met_y    -= et * sin(cluster_pos.phi());
    cl_hemet_x  -= et * cos(cluster_pos.phi());
    cl_hemet_y  -= et * sin(cluster_pos.phi());
    cl_he_et    += et;      
  }

  //---------------HF HAD PFClusters---------------

  for (reco::PFClusterCollection::const_iterator it = clustersHFHAD->begin(); it != clustersHFHAD->end(); it++){

    if ( it->energy() < hfh_cluster_threshold_) continue;
 
    const REPPoint cluster_pos = it->positionREP();	  
    double et = it->energy() / cosh( cluster_pos.eta() ); 
  
    //printCluster( it );
  
    pf_hfhcluster_et	->push_back(et			);
    pf_hfhcluster_e	->push_back(it->energy()	);
    pf_hfhcluster_eta	->push_back(cluster_pos.eta()	);
    pf_hfhcluster_phi	->push_back(cluster_pos.phi()	);

    pf_cluster_et	->push_back(et			);
    pf_cluster_e	->push_back(it->energy()	);
    pf_cluster_eta	->push_back(cluster_pos.eta()	);
    pf_cluster_phi	->push_back(cluster_pos.phi()	);
    
    cl_met_x     -= et * cos(cluster_pos.phi());
    cl_met_y     -= et * sin(cluster_pos.phi());
    cl_hfhmet_x  -= et * cos(cluster_pos.phi());
    cl_hfhmet_y  -= et * sin(cluster_pos.phi());
    cl_hfhad_et  += et;      
  }

  //---------------HF EM PFClusters---------------

  for (reco::PFClusterCollection::const_iterator it = clustersHFEM->begin(); it != clustersHFEM->end(); it++){

    if ( it->energy() < hfe_cluster_threshold_) continue;

    const REPPoint cluster_pos = it->positionREP();	  
    double et = it->energy() / cosh( cluster_pos.eta() ); 
    
    //printCluster( it );

    pf_hfecluster_et	->push_back(et			);
    pf_hfecluster_e	->push_back(it->energy()	);
    pf_hfecluster_eta	->push_back(cluster_pos.eta()	);
    pf_hfecluster_phi	->push_back(cluster_pos.phi()	);

    pf_cluster_et	->push_back(et			);
    pf_cluster_e	->push_back(it->energy()	);
    pf_cluster_eta	->push_back(cluster_pos.eta()	);
    pf_cluster_phi	->push_back(cluster_pos.phi()	);
    
    cl_met_x     -= et * cos(cluster_pos.phi());
    cl_met_y     -= et * sin(cluster_pos.phi());
    cl_hfemet_x  -= et * cos(cluster_pos.phi());
    cl_hfemet_y  -= et * sin(cluster_pos.phi());
    cl_hfem_et   += et;      
  }

  *pfclus_hfe_sumet = cl_hfem_et;
  *pfclus_hfe_met   = sqrt(cl_hfemet_x * cl_hfemet_x + cl_hfemet_y * cl_hfemet_y);

  *pfclus_hfh_sumet = cl_hfhad_et;
  *pfclus_hfh_met   = sqrt(cl_hfhmet_x * cl_hfhmet_x + cl_hfhmet_y * cl_hfhmet_y);

  *pfclus_hb_sumet  = cl_hb_et;
  *pfclus_hb_met    = sqrt(cl_hbmet_x * cl_hbmet_x + cl_hbmet_y * cl_hbmet_y);

  *pfclus_he_sumet  = cl_he_et;
  *pfclus_he_met    = sqrt(cl_hemet_x * cl_hemet_x + cl_hemet_y * cl_hemet_y);

  *pfclus_eb_sumet  = cl_eb_et;
  *pfclus_eb_met    = sqrt(cl_ebmet_x * cl_ebmet_x + cl_ebmet_y * cl_ebmet_y);

  *pfclus_ee_sumet  = cl_ee_et;
  *pfclus_ee_met    = sqrt(cl_eemet_x * cl_eemet_x + cl_eemet_y * cl_eemet_y);

  *pfclussumet      = cl_eb_et + cl_ee_et + cl_hb_et + cl_he_et + cl_hfem_et + cl_hfhad_et;
  *pfclusmet        = sqrt(cl_met_x * cl_met_x + cl_met_y * cl_met_y);

  cout << " EB met "   << sqrt(cl_ebmet_x * cl_ebmet_x + cl_ebmet_y * cl_ebmet_y) 
       << " EB sumet " << cl_eb_et << endl;


  // calculate MET and sumET using pf rechits
  float met_x = 0.;
  float met_y = 0.;

  // get EM energy in HF using pf rechits
  float hfem_et   = 0.;
  float hfemet_x  = 0.;
  float hfemet_y  = 0.;

  // get HAD energy in HF using pf rechits
  float hfhad_et  = 0.;
  float hfhmet_x  = 0.;
  float hfhmet_y  = 0.;

  // get energy HCAL barrel using pf rechits
  float hb_et	= 0.;
  float hbmet_x	= 0.;
  float hbmet_y	= 0.;

  // get energy ECAL endcap using pf rechits
  float he_et	= 0.;
  float hemet_x	= 0.;
  float hemet_y	= 0.;

  // get energy ECAL barrel using pf rechits
  float eb_et	= 0.;
  float ebmet_x	= 0.;
  float ebmet_y	= 0.;

  // get energy ECAL endcap using pf rechits
  float ee_et	= 0.;
  float eemet_x	= 0.;
  float eemet_y	= 0.;

  for (reco::PFRecHitCollection::const_iterator it = pf_hfem_h->begin(); 
       it != pf_hfem_h->end(); it++)
    {
      if (it->energy() < hfe_threshold_)
        continue;

      pf_rechit_detid->push_back(it->detId());
	  
      float et = sqrt(it->pt2());
      hfem_et += et;

      const REPPoint rechit_pos = it->positionREP();	  

      pf_hferechit_et	->push_back(et			);
      pf_hferechit_e	->push_back(it->energy()	);
      pf_hferechit_eta	->push_back(rechit_pos.eta()	);
      pf_hferechit_phi	->push_back(rechit_pos.phi()	);
      pf_rechit_et	->push_back(et			);
      pf_rechit_e	->push_back(it->energy()	);
      pf_rechit_eta	->push_back(rechit_pos.eta()	);
      pf_rechit_phi	->push_back(rechit_pos.phi()	);

      met_x    -= et * cos(rechit_pos.phi());
      met_y    -= et * sin(rechit_pos.phi());
      hfemet_x -= et * cos(rechit_pos.phi());
      hfemet_y -= et * sin(rechit_pos.phi());
    }

  *pf_hfe_sumet = hfem_et;
  *pf_hfe_met   = sqrt(hfemet_x * hfemet_x + hfemet_y * hfemet_y);



  for (reco::PFRecHitCollection::const_iterator it = pf_hfhad_h->begin(); 
       it != pf_hfhad_h->end(); it++)
    {
      if (it->energy() < hfh_threshold_)
        continue;

      pf_rechit_detid->push_back(it->detId());

      float et = sqrt(it->pt2());
      hfhad_et += et;

      const REPPoint& rechit_pos = it->positionREP();

      pf_hfhrechit_et	->push_back(et			);
      pf_hfhrechit_e	->push_back(it->energy()        );
      pf_hfhrechit_eta	->push_back(rechit_pos.eta()	);
      pf_hfhrechit_phi	->push_back(rechit_pos.phi()	);
      pf_rechit_et	->push_back(et			);
      pf_rechit_e	->push_back(it->energy()	);
      pf_rechit_eta	->push_back(rechit_pos.eta()	);
      pf_rechit_phi	->push_back(rechit_pos.phi()	);

      met_x    -= et * cos(rechit_pos.phi());
      met_y    -= et * sin(rechit_pos.phi());
      hfhmet_x -= et * cos(rechit_pos.phi());
      hfhmet_y -= et * sin(rechit_pos.phi());
    }

  *pf_hfh_sumet = hfhad_et;
  *pf_hfh_met   = sqrt(hfhmet_x * hfhmet_x + hfhmet_y * hfhmet_y);
     

       
  for (reco::PFRecHitCollection::const_iterator it = pf_hcal_h->begin(); it != pf_hcal_h->end(); it++)
    {
      if (it->layer() == PFLayer::HCAL_ENDCAP){
        
        if (it->energy() < he_threshold_)
          continue;
      }


      else if (it->layer() == PFLayer::HCAL_BARREL1){
        
        if (it->energy() < he_threshold_)
          continue;
      }
      
      pf_rechit_detid->push_back(it->detId());

      float et = sqrt(it->pt2());
      
      reco::PFRecHit rhit = reco::PFRecHit(*it);
      rhit.calculatePositionREP();
      const REPPoint rechit_pos = it->positionREP();	  
      
      if (it->layer() == PFLayer::HCAL_ENDCAP)
        {
          he_et += et;
          pf_herechit_et	->push_back(et			);
          pf_herechit_e 	->push_back(it->energy()	);
          pf_herechit_eta	->push_back(rechit_pos.eta()	);
          pf_herechit_phi	->push_back(rechit_pos.phi()	);
          pf_rechit_et	        ->push_back(et			);
          pf_rechit_e	        ->push_back(it->energy()	);
          pf_rechit_eta	        ->push_back(rechit_pos.eta()	);
          pf_rechit_phi	        ->push_back(rechit_pos.phi()	);		  
          
          met_x   -= et * cos(rechit_pos.phi());
          met_y   -= et * sin(rechit_pos.phi());
          hemet_x -= et * cos(rechit_pos.phi());
          hemet_y -= et * sin(rechit_pos.phi());
        }
      
      else if (it->layer() == PFLayer::HCAL_BARREL1)
        {
          hb_et += et;
          pf_hbrechit_et	->push_back(et			);
          pf_hbrechit_e 	->push_back(it->energy()	);
          pf_hbrechit_eta	->push_back(rechit_pos.eta()	);
          pf_hbrechit_phi	->push_back(rechit_pos.phi()	);
          pf_rechit_et	->push_back(et			);
          pf_rechit_e 	->push_back(it->energy()	);
          pf_rechit_eta	->push_back(rechit_pos.eta()	);
          pf_rechit_phi	->push_back(rechit_pos.phi()	);		  
          
          met_x   -= et * cos(rechit_pos.phi());
          met_y   -= et * sin(rechit_pos.phi());
          hbmet_x -= et * cos(rechit_pos.phi());
          hbmet_y -= et * sin(rechit_pos.phi());
        }
    }


  *pf_hb_sumet = hb_et;
  *pf_hb_met   = sqrt(hbmet_x * hbmet_x + hbmet_y * hbmet_y);
  *pf_he_sumet = he_et;
  *pf_he_met   = sqrt(hemet_x * hemet_x + hemet_y * hemet_y);



  for (reco::PFRecHitCollection::const_iterator it = pf_ecal_h->begin(); it != pf_ecal_h->end(); it++)
    {
      float et = sqrt(it->pt2());

      if (it->layer() == PFLayer::ECAL_BARREL){
        
        if (it->energy() < eb_threshold_)
          continue;
      }


      else if (it->layer() == PFLayer::ECAL_ENDCAP){
        
        if (it->energy() < ee_threshold_)
          continue;
      }
      
      pf_rechit_detid->push_back(it->detId());

      reco::PFRecHit rhit = reco::PFRecHit(*it);
      rhit.calculatePositionREP();
      const REPPoint rechit_pos = it->positionREP();	  

      if (it->layer() == PFLayer::ECAL_BARREL)
        {
          eb_et += et;

          pf_ebrechit_et	->push_back(et			);
          pf_ebrechit_e	        ->push_back(it->energy()			);
          pf_ebrechit_eta	->push_back(rechit_pos.eta()	);
          pf_ebrechit_phi	->push_back(rechit_pos.phi()	);
          pf_rechit_et	->push_back(et			);
          pf_rechit_e	->push_back(it->energy()			);
          pf_rechit_eta	->push_back(rechit_pos.eta()	);
          pf_rechit_phi	->push_back(rechit_pos.phi()	);

          met_x   -= et * cos(rechit_pos.phi());
          met_y   -= et * sin(rechit_pos.phi());
          ebmet_x -= et * cos(rechit_pos.phi());
          ebmet_y -= et * sin(rechit_pos.phi());
        }
      else if (it->layer() == PFLayer::ECAL_ENDCAP)
        {
          ee_et += et;

          pf_eerechit_et	->push_back(et			);
          pf_eerechit_e	        ->push_back(it->energy()	);
          pf_eerechit_eta	->push_back(rechit_pos.eta()	);
          pf_eerechit_phi	->push_back(rechit_pos.phi()	);
          pf_rechit_et	->push_back(et			);
          pf_rechit_e	->push_back(it->energy()	);
          pf_rechit_eta	->push_back(rechit_pos.eta()	);
          pf_rechit_phi	->push_back(rechit_pos.phi()	);

          met_x   -= et * cos(rechit_pos.phi());
          met_y   -= et * sin(rechit_pos.phi());
          eemet_x -= et * cos(rechit_pos.phi());
          eemet_y -= et * sin(rechit_pos.phi());
        }
    }

  *pf_eb_sumet = eb_et;
  *pf_eb_met   = sqrt(ebmet_x * ebmet_x + ebmet_y * ebmet_y);
  *pf_ee_sumet = ee_et;
  *pf_ee_met   = sqrt(eemet_x * eemet_x + eemet_y * eemet_y);

  // combine to get met, sumet using pf rechits
  *pfsumet = eb_et + ee_et + hb_et + he_et + hfem_et + hfhad_et;
  *pfmet   = sqrt(met_x * met_x + met_y * met_y);

  // put branches in tree
  iEvent.put(calo_eb_sumet		, "caloebsumet"		);
  iEvent.put(calo_ee_sumet		, "caloeesumet"		);
  iEvent.put(calo_hb_sumet		, "calohbsumet"		);
  iEvent.put(calo_he_sumet		, "calohesumet"		);
  iEvent.put(calo_hfe_sumet		, "calohfesumet"	);
  iEvent.put(calo_hfh_sumet		, "calohfhsumet"	);
  iEvent.put(pf_eb_sumet		, "pfebsumet"		);
  iEvent.put(pf_ee_sumet		, "pfeesumet"		);
  iEvent.put(pf_hb_sumet		, "pfhbsumet"		);
  iEvent.put(pf_he_sumet		, "pfhesumet"		);
  iEvent.put(pf_hfe_sumet		, "pfhfesumet"		);
  iEvent.put(pf_hfh_sumet		, "pfhfhsumet"		);
  iEvent.put(pfclus_eb_sumet		, "pfclusebsumet"	);
  iEvent.put(pfclus_ee_sumet		, "pfcluseesumet"	);
  iEvent.put(pfclus_hb_sumet		, "pfclushbsumet"	);
  iEvent.put(pfclus_he_sumet		, "pfclushesumet"	);
  iEvent.put(pfclus_hfe_sumet		, "pfclushfesumet"	);
  iEvent.put(pfclus_hfh_sumet		, "pfclushfhsumet"	);
  iEvent.put(calomet			, "calomet"		);
  iEvent.put(calosumet		        , "calosumet"		);
  iEvent.put(pfmet			, "pfmet"		);
  iEvent.put(pfsumet			, "pfsumet"		);
  iEvent.put(pfclusmet			, "pfclusmet"		);
  iEvent.put(pfclussumet		, "pfclussumet"		);
  iEvent.put(genmet                     , "genmet"              );
  iEvent.put(gensumet                   , "gensumet"            );          
  iEvent.put(pf_eb_met	             	, "pfebmet"		);
  iEvent.put(pf_ee_met		        , "pfeemet"		);
  iEvent.put(pf_hb_met		        , "pfhbmet"		);
  iEvent.put(pf_he_met		        , "pfhemet"		);
  iEvent.put(pf_hfe_met		        , "pfhfemet"		);
  iEvent.put(pf_hfh_met		        , "pfhfhmet"		);
  iEvent.put(pfclus_eb_met	      	, "pfclusebmet"		);
  iEvent.put(pfclus_ee_met		, "pfcluseemet"		);
  iEvent.put(pfclus_hb_met		, "pfclushbmet"		);
  iEvent.put(pfclus_he_met		, "pfclushemet"		);
  iEvent.put(pfclus_hfe_met		, "pfclushfemet"	);
  iEvent.put(pfclus_hfh_met		, "pfclushfhmet"	);
  iEvent.put(pf_rechit_et		, "pfrechitet"		);
  iEvent.put(pf_rechit_e		, "pfrechite"		);
  iEvent.put(pf_rechit_eta		, "pfrechiteta"		);
  iEvent.put(pf_rechit_phi		, "pfrechitphi"		);
  iEvent.put(pf_ebrechit_et		, "pfebrechitet"	);
  iEvent.put(pf_ebrechit_e		, "pfebrechite"  	);
  iEvent.put(pf_ebrechit_eta		, "pfebrechiteta"	);
  iEvent.put(pf_ebrechit_phi		, "pfebrechitphi"	);
  iEvent.put(pf_eerechit_et		, "pfeerechitet"	);
  iEvent.put(pf_eerechit_e		, "pfeerechite"	        );
  iEvent.put(pf_eerechit_eta		, "pfeerechiteta"	);
  iEvent.put(pf_eerechit_phi		, "pfeerechitphi"	);
  iEvent.put(pf_hbrechit_et		, "pfhbrechitet"	);
  iEvent.put(pf_hbrechit_e		, "pfhbrechite"      	);
  iEvent.put(pf_hbrechit_eta		, "pfhbrechiteta"	);
  iEvent.put(pf_hbrechit_phi		, "pfhbrechitphi"	);
  iEvent.put(pf_herechit_et		, "pfherechitet"	);
  iEvent.put(pf_herechit_e		, "pfherechite"     	);
  iEvent.put(pf_herechit_eta		, "pfherechiteta"	);
  iEvent.put(pf_herechit_phi		, "pfherechitphi"	);
  iEvent.put(pf_hferechit_et		, "pfhferechitet"	);
  iEvent.put(pf_hferechit_e		, "pfhferechite"	);
  iEvent.put(pf_hferechit_eta	        , "pfhferechiteta"	);
  iEvent.put(pf_hferechit_phi	        , "pfhferechitphi"	);
  iEvent.put(pf_hfhrechit_et		, "pfhfhrechitet"	);
  iEvent.put(pf_hfhrechit_e		, "pfhfhrechite"	);
  iEvent.put(pf_hfhrechit_eta	        , "pfhfhrechiteta"	);
  iEvent.put(pf_hfhrechit_phi	        , "pfhfhrechitphi"	);
  iEvent.put(pf_cluster_et		, "pfclusteret"		);
  iEvent.put(pf_cluster_e		, "pfclustere"		);
  iEvent.put(pf_cluster_eta		, "pfclustereta"	);
  iEvent.put(pf_cluster_phi		, "pfclusterphi"	);
  iEvent.put(pf_ebcluster_et		, "pfebclusteret"	);
  iEvent.put(pf_ebcluster_e		, "pfebclustere"  	);
  iEvent.put(pf_ebcluster_eta		, "pfebclustereta"	);
  iEvent.put(pf_ebcluster_phi		, "pfebclusterphi"	);
  iEvent.put(pf_eecluster_et		, "pfeeclusteret"	);
  iEvent.put(pf_eecluster_e		, "pfeeclustere"	);
  iEvent.put(pf_eecluster_eta		, "pfeeclustereta"	);
  iEvent.put(pf_eecluster_phi		, "pfeeclusterphi"	);
  iEvent.put(pf_hbcluster_et		, "pfhbclusteret"	);
  iEvent.put(pf_hbcluster_e		, "pfhbclustere"      	);
  iEvent.put(pf_hbcluster_eta		, "pfhbclustereta"	);
  iEvent.put(pf_hbcluster_phi		, "pfhbclusterphi"	);
  iEvent.put(pf_hecluster_et		, "pfheclusteret"	);
  iEvent.put(pf_hecluster_e		, "pfheclustere"     	);
  iEvent.put(pf_hecluster_eta		, "pfheclustereta"	);
  iEvent.put(pf_hecluster_phi		, "pfheclusterphi"	);
  iEvent.put(pf_hfecluster_et		, "pfhfeclusteret"	);
  iEvent.put(pf_hfecluster_e		, "pfhfeclustere"	);
  iEvent.put(pf_hfecluster_eta	        , "pfhfeclustereta"	);
  iEvent.put(pf_hfecluster_phi	        , "pfhfeclusterphi"	);
  iEvent.put(pf_hfhcluster_et		, "pfhfhclusteret"	);
  iEvent.put(pf_hfhcluster_e		, "pfhfhclustere"	);
  iEvent.put(pf_hfhcluster_eta	        , "pfhfhclustereta"	);
  iEvent.put(pf_hfhcluster_phi	        , "pfhfhclusterphi"	);

}


void PFRecHitAnalyzer::printCluster( const reco::PFClusterCollection::const_iterator& cluster ) {
  
   const math::XYZPoint&  pos = cluster.position();
   const PFCluster::REPPoint&  posrep = cluster.positionREP();
   const std::vector< reco::PFRecHitFraction >& fracs = 
     cluster.recHitFractions();
   
   cout<<"PFcluster "
      <<", layer: "<<cluster.layer()
      <<"\tE = "<<cluster.energy()
      <<"\tXYZ: "
      <<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<" | "
      <<"\tREP: "
      <<posrep.Rho()<<","<<posrep.Eta()<<","<<posrep.Phi()<<" | "
      <<fracs.size()<<" rechits";
   
   for(unsigned i=0; i<fracs.size(); i++) {
     // PFRecHit is not available, print the detID
     if( !fracs[i].recHitRef().isAvailable() )
       cout<<cluster.printHitAndFraction(i)<<", ";
     else 
       cout<<fracs[i]<<", ";
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
PFRecHitAnalyzer::beginJob()
{
  cout << "PFRecHitAnalyzer::beginJob()" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFRecHitAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFRecHitAnalyzer);
