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
// $Id: PFRecHitAnalyzer.h,v 1.1 2010/06/30 14:58:37 benhoob Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"


//
// class declaration
//

class PFRecHitAnalyzer : public edm::EDProducer {
public:
  explicit PFRecHitAnalyzer(const edm::ParameterSet&);
  ~PFRecHitAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  //void printCluster( const reco::PFCluster& cluster );
  void printCluster( const reco::PFClusterCollection::const_iterator& cluster );
  
  // ----------member data ---------------------------
  edm::InputTag pf_ecal_tag_;
  edm::InputTag pf_hcal_tag_;
  edm::InputTag pf_hfem_tag_;
  edm::InputTag pf_hfhad_tag_;
  edm::InputTag calomet_tag_;
  edm::InputTag genmet_tag_;
  edm::InputTag inputTagPFClustersECAL_;
  edm::InputTag inputTagPFClustersHCAL_;
  edm::InputTag inputTagPFClustersHFEM_;
  edm::InputTag inputTagPFClustersHFHAD_;   
  
  double eb_threshold_;
  double ee_threshold_;
  double hb_threshold_;
  double he_threshold_;
  double hfe_threshold_;
  double hfh_threshold_;

  double eb_cluster_threshold_;
  double ee_cluster_threshold_;
  double hb_cluster_threshold_;
  double he_cluster_threshold_;
  double hfe_cluster_threshold_;
  double hfh_cluster_threshold_;

  bool useClusters_;
};
