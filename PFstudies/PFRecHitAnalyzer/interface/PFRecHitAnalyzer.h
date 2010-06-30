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
// $Id$
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
      
      // ----------member data ---------------------------
     edm::InputTag pf_ecal_tag_;
     edm::InputTag pf_hcal_tag_;
     edm::InputTag pf_hfem_tag_;
     edm::InputTag pf_hfhad_tag_;
     edm::InputTag calomet_tag_;
     edm::InputTag genmet_tag_;
     double eb_threshold_;
     double ee_threshold_;
     double hb_threshold_;
     double he_threshold_;
     double hfe_threshold_;
     double hfh_threshold_;
};
