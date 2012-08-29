#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

#include <iostream>
#include "Math/LorentzVector.h"
 
#include <cmath>

using namespace std;

class StopTree;

class StopTreeLooper {

    public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        StopTreeLooper();
        ~StopTreeLooper();

        void setOutFileName(string filename); 
        void loop(TChain *chain, TString name);

	bool passSingleMuonSelection(const StopTree *sTree);

	float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );

    private:

	string m_outfilename_;



};

#endif
