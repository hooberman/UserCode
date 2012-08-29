#ifndef StopTree_H
#define StopTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"

#include "Math/LorentzVector.h"

#include <cmath>
#include "assert.h"

using namespace std;
using namespace ROOT::Math;

//
// Ntuple structure:
//
// Ntuple content:

class StopTree {

    public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        /// variables
	int           ndavtx_;
	float         ndavtxweight_;
	float         weight_;
	float         mutrigweight_;
	float         rhovor_;
	float         mgcor_;
	int           smu_;
	int           smu30_;
	int           trgmu30_;
	int           dil_;
	int           ngoodlep_;
	int           leptype_;
	int           id1_;
	int           id2_;
	int           npfjets30_;
	float 	      lep1chi2ndf_;	
	float 	      lep1dpt_;	
	float         t1metphicorr_;
	float         t1metphicorrmt_;
	float         t1metphicorrphi_;
	int           nbtagsssv_;
	int           nbtagsssvcorr_;
	float         pfcandpt10_;
	float         pfcandiso10_;
	int           nleps_;         
	LorentzVector lep1_;
	LorentzVector lep2_;
	LorentzVector pfjet1_;
	LorentzVector pfjet2_;
	LorentzVector pfjet3_;
	LorentzVector pfjet4_;
	LorentzVector pfjet5_;
	LorentzVector pfjet6_;

    public:
        /// this is the main element
        TTree *tree_;
        TFile *f_;

        /// hold the names of variables to facilitate things (filled during Init)
        vector<string> variables_;
	
        /// default constructor  
 StopTree() :  lep1Ptr_(&lep1_), lep2Ptr_(&lep2_), jet1Ptr_(&pfjet1_), jet2Ptr_(&pfjet2_), jet3Ptr_(&pfjet3_), jet4Ptr_(&pfjet4_), jet5Ptr_(&pfjet5_), jet6Ptr_(&pfjet6_) {}
        /// default destructor
        ~StopTree(){ 
	  cout << "~StopTree()" << endl;
	  if (f_) f_->Close();  
	  cout << "~StopTree() done" << endl;
	  
        };

        /// initialize varibles and fill list of available variables
        void InitVariables();

        /// load a StopTree
        void LoadTree(const char* file){
            f_ = TFile::Open(file);
            assert(f_);
            tree_ = dynamic_cast<TTree*>(f_->Get("t"));
            assert(tree_);
        }

        /// create a StopTree
        void CreateTree(){
            tree_ = new TTree("t","stop babytuple");
            f_ = 0;
            InitVariables();
            //book the branches
	    tree_->Branch("ndavtx", 		&ndavtx_, 		"ndavtx/I");	      	      
	    tree_->Branch("ndavtxweight", 	&ndavtxweight_, 	"ndavtxweight/F");    
	    tree_->Branch("weight", 		&weight_, 		"weight/F");	      	      
	    tree_->Branch("mutrigweight", 	&mutrigweight_, 	"mutrigweight/F");    
	    tree_->Branch("rhovor", 		&rhovor_, 		"rhovor/F");	      	      
	    tree_->Branch("mgcor", 		&mgcor_, 		"mgcor/F");	      	      
	    tree_->Branch("smu", 		&smu_, 			"smu/I");	      	      
	    tree_->Branch("smu30", 		&smu30_, 		"smu30/I");	      	      
	    tree_->Branch("trgmu30", 		&trgmu30_, 		"trgmu30/I");	      	      
	    tree_->Branch("dil", 		&dil_, 			"dil/I");	      	      
	    tree_->Branch("ngoodlep", 		&ngoodlep_, 		"ngoodlep/I");            
	    tree_->Branch("leptype", 		&leptype_, 		"leptype/I");	      	      
	    tree_->Branch("id1", 		&id1_, 			"id1/I");	      	      
	    tree_->Branch("id2", 		&id2_, 			"id2/I");	      	      
	    tree_->Branch("npfjets30", 		&npfjets30_, 		"npfjets30/I");          
	    tree_->Branch("lep1chi2ndf", 	&lep1chi2ndf_, 		"lep1chi2ndf/F");	      	      
	    tree_->Branch("lep1dpt", 		&lep1dpt_, 		"lep1dpt/F");	      	      
	    tree_->Branch("t1metphicorr", 	&t1metphicorr_, 	"t1metphicorr/F");	      	      
	    tree_->Branch("t1metphicorrmt", 	&t1metphicorrmt_, 	"t1metphicorrmt/F");          
	    tree_->Branch("t1metphicorrphi", 	&t1metphicorrphi_, 	"t1metphicorrphi/F");        
	    tree_->Branch("nbtagsssv", 		&nbtagsssv_, 		"nbtagsssv/I");          
	    tree_->Branch("nbtagsssvcorr", 	&nbtagsssvcorr_, 	"nbtagsssvcorr/I");  
	    tree_->Branch("pfcandpt10", 	&pfcandpt10_, 		"pfcandpt10/F");        
	    tree_->Branch("pfcandiso10", 	&pfcandiso10_, 		"pfcandiso10/F");      
	    tree_->Branch("nleps", 		&nleps_, 		"nleps/I");     
            tree_->Branch("lep1",   "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep1Ptr_);
            tree_->Branch("lep2",   "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &lep2Ptr_);
            tree_->Branch("pfjet1", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet1Ptr_);
            tree_->Branch("pfjet2", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet2Ptr_);
            tree_->Branch("pfjet3", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet3Ptr_);
            tree_->Branch("pfjet4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet4Ptr_);
            tree_->Branch("pfjet5", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet5Ptr_);
            tree_->Branch("pfjet6", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &jet6Ptr_);
             
        }

        // initialze a StopTree
        void InitTree(){
            assert(tree_);
            // don't forget to set pointers to zero before you set address
            // or you will fully appreciate that "ROOT sucks" :)
            InitVariables();
            //Set branch address
            Int_t currentState = gErrorIgnoreLevel;
            // gErrorIgnoreLevel = kError;
            gErrorIgnoreLevel = kBreak;
	    tree_->SetBranchAddress("ndavtx", 		&ndavtx_);	      	      
	    tree_->SetBranchAddress("ndavtxweight", 	&ndavtxweight_);    
	    tree_->SetBranchAddress("weight", 		&weight_);	      	      
	    tree_->SetBranchAddress("mutrigweight", 	&mutrigweight_);    
	    tree_->SetBranchAddress("rhovor", 		&rhovor_);	      	      
	    tree_->SetBranchAddress("mgcor", 		&mgcor_);	      	      
	    tree_->SetBranchAddress("smu", 		&smu_);	      	      
	    tree_->SetBranchAddress("smu30", 		&smu30_);	      	      
	    tree_->SetBranchAddress("trgmu30", 		&trgmu30_);	      	      
	    tree_->SetBranchAddress("dil", 		&dil_);	      	      
	    tree_->SetBranchAddress("ngoodlep", 	&ngoodlep_);            
	    tree_->SetBranchAddress("leptype", 		&leptype_);	      	      
	    tree_->SetBranchAddress("id1", 		&id1_);	      	      
	    tree_->SetBranchAddress("id2", 		&id2_);	      	      
	    tree_->SetBranchAddress("npfjets30", 	&npfjets30_);          
	    tree_->SetBranchAddress("lep1chi2ndf", 	&lep1chi2ndf_);	      	      
	    tree_->SetBranchAddress("lep1dpt", 		&lep1dpt_);	      	      
	    tree_->SetBranchAddress("t1metphicorr", 	&t1metphicorr_);	      	      
	    tree_->SetBranchAddress("t1metphicorrmt", 	&t1metphicorrmt_);          
	    tree_->SetBranchAddress("t1metphicorrphi", 	&t1metphicorrphi_);        
	    tree_->SetBranchAddress("nbtagsssv", 	&nbtagsssv_);          
	    tree_->SetBranchAddress("nbtagsssvcorr", 	&nbtagsssvcorr_);  
	    tree_->SetBranchAddress("pfcandpt10", 	&pfcandpt10_);        
	    tree_->SetBranchAddress("pfcandiso10", 	&pfcandiso10_);      
	    tree_->SetBranchAddress("nleps", 		&nleps_);     
            tree_->SetBranchAddress("lep1",   		&lep1Ptr_);
            tree_->SetBranchAddress("lep2",   		&lep2Ptr_);
            tree_->SetBranchAddress("pfjet1", 		&jet1Ptr_);
            tree_->SetBranchAddress("pfjet2", 		&jet2Ptr_);
            tree_->SetBranchAddress("pfjet3", 		&jet3Ptr_);
            tree_->SetBranchAddress("pfjet4", 		&jet4Ptr_);
            tree_->SetBranchAddress("pfjet5", 		&jet5Ptr_);
            tree_->SetBranchAddress("pfjet6", 		&jet6Ptr_);

            gErrorIgnoreLevel = currentState;
        }

        /// get a built in type variable by name
        double Get(string value);
        /// compare two StopTrees for a given event on a given level of precision; 
        /// returns the variables that failed the comparison 
        vector<string> Compare(StopTree* value, double prec=0.005);

    private:

        LorentzVector *lep1Ptr_;
        LorentzVector *lep2Ptr_;
        LorentzVector *jet1Ptr_;
        LorentzVector *jet2Ptr_;
        LorentzVector *jet3Ptr_;
        LorentzVector *jet4Ptr_;
        LorentzVector *jet5Ptr_;
        LorentzVector *jet6Ptr_;

}; 

inline void 
StopTree::InitVariables(){
    // create list of available variables
    if(variables_.empty()){
        //make sure that this is only done once
	variables_.push_back(string("ndavtx"		));
	variables_.push_back(string("ndavtxweight"	));
	variables_.push_back(string("weight"		));
	variables_.push_back(string("mutrigweight"	));
	variables_.push_back(string("rhovor"		));
	variables_.push_back(string("mgcor"		));
	variables_.push_back(string("smu"		));
	variables_.push_back(string("smu30"		));
	variables_.push_back(string("trgmu30"		));
	variables_.push_back(string("dil"		));
	variables_.push_back(string("ngoodlep"		));
	variables_.push_back(string("leptype"		));
	variables_.push_back(string("id1"		));
	variables_.push_back(string("id2"		));
	variables_.push_back(string("npfjets30"		));
	variables_.push_back(string("lep1chi2ndf"	));
	variables_.push_back(string("lep1dpt"		));
	variables_.push_back(string("t1metphicorr"	));
	variables_.push_back(string("t1metphicorrmt"	));
	variables_.push_back(string("t1metphicorrphi"	));
	variables_.push_back(string("nbtagsssv"		));
	variables_.push_back(string("nbtagsssvcorr"	));
	variables_.push_back(string("pfcandpt10"	));
	variables_.push_back(string("pfcandiso10"	));
	variables_.push_back(string("nleps"		));         
	variables_.push_back(string("lep1"		));
	variables_.push_back(string("lep2"		));
	variables_.push_back(string("pfjet1"		));
	variables_.push_back(string("pfjet2"		));
	variables_.push_back(string("pfjet3"		));
	variables_.push_back(string("pfjet4"		));
	variables_.push_back(string("pfjet5"		));
	variables_.push_back(string("pfjet6"		));
    }

    // inizialize variables
    ndavtx_		= 9;
    ndavtxweight_	= -999.;
    weight_		= -999.;
    mutrigweight_	= -999.;
    rhovor_		= -999.;
    mgcor_		= -999.;
    smu_		= 999;
    smu30_		= 999;
    trgmu30_		= 999;
    dil_		= 999;
    ngoodlep_		= 999;
    leptype_		= 999;
    id1_		= 999;
    id2_		= 999;
    npfjets30_		= 999;
    lep1chi2ndf_	= -999.;
    lep1dpt_		= -999.;
    t1metphicorr_	= -999.;
    t1metphicorrmt_	= -999.;
    t1metphicorrphi_	= -999.;
    nbtagsssv_		= 999;
    nbtagsssvcorr_	= 999;
    pfcandpt10_		= -999.;
    pfcandiso10_	= -999.;
    nleps_		= 999;         
    lep1_		= LorentzVector();
    lep2_		= LorentzVector();
    pfjet1_		= LorentzVector();
    pfjet2_		= LorentzVector();
    pfjet3_		= LorentzVector();
    pfjet4_		= LorentzVector();
    pfjet5_		= LorentzVector();
    pfjet6_		= LorentzVector();
    
}

    inline double
StopTree::Get(string value)
{

  if(value=="ndavtx" 		) { return this->ndavtx_;		}	      	      
  if(value=="ndavtxweight" 	) { return this->ndavtxweight_;		}    
  if(value=="weight" 		) { return this->weight_;		}	      	      
  if(value=="mutrigweight" 	) { return this->mutrigweight_;		}    
  if(value=="rhovor" 		) { return this->rhovor_;		}	      	      
  if(value=="mgcor" 		) { return this->mgcor_;		}	      	      
  if(value=="smu" 		) { return this->smu_;			}	      	      
  if(value=="smu30" 		) { return this->smu30_;		}	      	      
  if(value=="trgmu30" 		) { return this->trgmu30_;		}	      	      
  if(value=="dil" 		) { return this->dil_;			}	      	      
  if(value=="ngoodlep" 		) { return this->ngoodlep_;		}            
  if(value=="leptype" 		) { return this->leptype_;		}	      	      
  if(value=="id1" 		) { return this->id1_;			}	      	      
  if(value=="id2" 		) { return this->id2_;			}	      	      
  if(value=="npfjets30" 	) { return this->npfjets30_;		}          
  if(value=="lep1chi2ndf" 	) { return this->lep1chi2ndf_;		}	      	      
  if(value=="lep1dpt" 		) { return this->lep1dpt_;		}	      	      
  if(value=="t1metphicorr" 	) { return this->t1metphicorr_;		}	      	      
  if(value=="t1metphicorrmt" 	) { return this->t1metphicorrmt_;	}          
  if(value=="t1metphicorrphi" 	) { return this->t1metphicorrphi_;	}        
  if(value=="nbtagsssv" 	) { return this->nbtagsssv_;		}          
  if(value=="nbtagsssvcorr" 	) { return this->nbtagsssvcorr_;	}  
  if(value=="pfcandpt10" 	) { return this->pfcandpt10_;		}        
  if(value=="pfcandiso10" 	) { return this->pfcandiso10_;		}      
  if(value=="nleps" 		) { return this->nleps_;		}     
  /* if(value=="lep1"   		) { return this->lep1_;		} */
  /* if(value=="lep2"   		) { return this->lep2_;		} */
  /* if(value=="pfjet1" 		) { return this->pfjet1_;	} */
  /* if(value=="pfjet2" 		) { return this->pfjet2_;	} */
  /* if(value=="pfjet3" 		) { return this->pfjet3_;	} */
  /* if(value=="pfjet4" 		) { return this->pfjet4_;	} */
  /* if(value=="pfjet5" 		) { return this->pfjet5_;	} */
  /* if(value=="pfjet6" 		) { return this->pfjet6_;	} */

  return -9999.; 

}

inline vector<string> 
StopTree::Compare(StopTree* value, double prec){
    vector<string> fails;
    // this should always fit with ultimate precision
    /* if( this->event_ != value->event_ ){ fails.push_back( "event" ); } */
    /* if( this->run_   != value->run_   ){ fails.push_back( "run"   ); } */
    /* if( this->lumi_  != value->lumi_  ){ fails.push_back( "lumi"  ); } */

    // check within  (relative) precision
    for(vector<string>::const_iterator var=variables_.begin(); var!=variables_.end(); ++var){
        if( fabs(Get(*var)-value->Get(*var))/Get(*var)>prec ) fails.push_back(*var);
    }
    return fails;
}


#endif
