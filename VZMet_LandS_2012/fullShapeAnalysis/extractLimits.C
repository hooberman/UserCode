#include "../../test/fitRvsCLs.C"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>

using namespace std;

char* version             = "V00-00-05";

bool fileInList(string thisfilename);

void extractLimits( bool print = false ){

  //------------------------------------------
  // create exclusion histogram
  //------------------------------------------

  TH2F* hexcl    = new TH2F( "hexcl"    , "hexcl"    , 41 , -5.0 , 405.0 , 41 , -5.0 , 405.0 );
  TH2F* hexp     = new TH2F( "hexp"     , "hexp"     , 41 , -5.0 , 405.0 , 41 , -5.0 , 405.0 );

  ofstream* doScript_failed = new ofstream();
  doScript_failed->open(Form("cards/%s/doLimits_failed.sh",version));

  //------------------------------------------
  // loop over CMSSM points
  //------------------------------------------

  vector<int>   mgvec;
  vector<int>   mlvec;
  vector<float> obsvec;
  vector<float> expvec;

  for( int mgbin = 1 ; mgbin <= hexcl->GetXaxis()->GetNbins() ; mgbin++ ){
    for( int mlbin = 1 ; mlbin <= hexcl->GetYaxis()->GetNbins() ; mlbin++ ){

      //------------------------------------------
      // restrict range
      //------------------------------------------

      if( mgbin > 38 ) continue;

      int mg = hexcl->GetXaxis()->GetBinCenter(mgbin);
      int ml = hexcl->GetXaxis()->GetBinCenter(mlbin);

      //if( mgbin == 15 && mlbin == 4  ) continue;
      //if( mgbin == 17 && mlbin == 6  ) continue;
      //if( mgbin == 39 && mlbin == 19 ) continue;

      hexcl->SetBinContent(mgbin,mlbin,0);
      hexp->SetBinContent(mgbin,mlbin,0);

      //------------------------------------------
      // open file, if available
      //------------------------------------------

      char* filename = Form("cards/%s/SMS_%i_%i_m2lnQ.root",version,mgbin,mlbin);

      bool found = fileInList( filename );
      if( !found ) continue;

      //------------------------------------------
      // check if point is excluded
      //------------------------------------------

      limitResult mylimit = run(filename,"plot");

      if( mylimit.obs < 1.e-10 ){
	*doScript_failed << Form("../../../../test/lands.exe -d SMS_%i_%i.txt -M Hybrid --freq  --nToysForCLsb 1500 --nToysForCLb 500  --scanRs 1 -vR [0.01,10.0,x1.1] -n SMS_%i_%i",mgbin,mlbin,mgbin,mlbin) << endl;
      }
      
      else{

	cout << "----------------------------------------" << endl;
	cout << " mgbin mlbin : " << mgbin << " " << mlbin << endl;
	cout << " mg    ml    : " << mg    << " " << ml    << endl;
	cout << " observed    : " << mylimit.obs           << endl;
	cout << " expected    : " << mylimit.exp           << endl;
	cout << "----------------------------------------" << endl;

	hexcl->SetBinContent(mgbin,mlbin,mylimit.obs);
	hexp->SetBinContent(mgbin,mlbin,mylimit.exp);

	mgvec.push_back(mg);
	mlvec.push_back(ml);
	obsvec.push_back(mylimit.obs);
	expvec.push_back(mylimit.exp);
      }
      
    }
  }

  for( int i = 0 ; i < mgvec.size() ; ++i ){
    cout << mgvec.at(i) << " " << mlvec.at(i) << " " << Form("%.2f",expvec.at(i)) << endl;
  }

  doScript_failed->close();
  
  //------------------------------------------
  // draw exclusion histogram
  //------------------------------------------

  TCanvas *can = new TCanvas("can","can",1000,800);
  can->cd();
  gPad->SetRightMargin(0.2);
  hexcl->GetXaxis()->SetLabelSize(0.03);
  hexcl->GetYaxis()->SetLabelSize(0.03);
  hexcl->GetXaxis()->SetTitle("gluino mass (GeV)");
  hexcl->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
  hexcl->Draw("colz");

  TFile* outfile = TFile::Open(Form("cards/%s/observed_limit.root",version),"RECREATE");
  hexcl->Write();
  hexp->Write();
  outfile->Close();

  if( print ){
    can->Print(Form("cards/%s/SMS.eps",version));
    can->Print(Form("cards/%s/SMS.png",version));
    gROOT->ProcessLine(Form(".! ps2pdf cards/%s/SMS.eps cards/%s/SMS.pdf",version,version));
  }

}

//------------------------------------------
// check if this file appears in file list
//------------------------------------------

bool fileInList(string thisfilename){

  ifstream* ifile = new ifstream();
  ifile->open(Form("cards/%s/file_list_CLs.txt",version));

  string filename;

  bool found = false;

  while( ifile->good() ){
    *ifile >> filename;
    if( filename == thisfilename ){
      found = true;
      break;
    }
  }

  ifile->close();
  delete ifile;

  return found;
}
