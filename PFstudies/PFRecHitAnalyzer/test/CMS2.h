// -*- C++ -*-
#ifndef CMS2_H
#define CMS2_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

#define PARANOIA

using namespace std; 
class CMS2 {
private: 
protected: 
	unsigned int index;
	float calo_eb_sumet_;
	TBranch *calo_eb_sumet_branch;
	bool calo_eb_sumet_isLoaded;
	float calo_ee_sumet_;
	TBranch *calo_ee_sumet_branch;
	bool calo_ee_sumet_isLoaded;
	float calo_hb_sumet_;
	TBranch *calo_hb_sumet_branch;
	bool calo_hb_sumet_isLoaded;
	float calo_he_sumet_;
	TBranch *calo_he_sumet_branch;
	bool calo_he_sumet_isLoaded;
	float calo_hfe_sumet_;
	TBranch *calo_hfe_sumet_branch;
	bool calo_hfe_sumet_isLoaded;
	float calo_hfh_sumet_;
	TBranch *calo_hfh_sumet_branch;
	bool calo_hfh_sumet_isLoaded;
	float calomet_;
	TBranch *calomet_branch;
	bool calomet_isLoaded;
	float calosumet_;
	TBranch *calosumet_branch;
	bool calosumet_isLoaded;
	float genmet_;
	TBranch *genmet_branch;
	bool genmet_isLoaded;
	float gensumet_;
	TBranch *gensumet_branch;
	bool gensumet_isLoaded;
	float pfclus_eb_met_;
	TBranch *pfclus_eb_met_branch;
	bool pfclus_eb_met_isLoaded;
	float pfclus_eb_sumet_;
	TBranch *pfclus_eb_sumet_branch;
	bool pfclus_eb_sumet_isLoaded;
	float pfclus_ee_met_;
	TBranch *pfclus_ee_met_branch;
	bool pfclus_ee_met_isLoaded;
	float pfclus_ee_sumet_;
	TBranch *pfclus_ee_sumet_branch;
	bool pfclus_ee_sumet_isLoaded;
	float pfclus_hb_met_;
	TBranch *pfclus_hb_met_branch;
	bool pfclus_hb_met_isLoaded;
	float pfclus_hb_sumet_;
	TBranch *pfclus_hb_sumet_branch;
	bool pfclus_hb_sumet_isLoaded;
	float pfclus_he_met_;
	TBranch *pfclus_he_met_branch;
	bool pfclus_he_met_isLoaded;
	float pfclus_he_sumet_;
	TBranch *pfclus_he_sumet_branch;
	bool pfclus_he_sumet_isLoaded;
	float pfclus_hfe_met_;
	TBranch *pfclus_hfe_met_branch;
	bool pfclus_hfe_met_isLoaded;
	float pfclus_hfe_sumet_;
	TBranch *pfclus_hfe_sumet_branch;
	bool pfclus_hfe_sumet_isLoaded;
	float pfclus_hfh_met_;
	TBranch *pfclus_hfh_met_branch;
	bool pfclus_hfh_met_isLoaded;
	float pfclus_hfh_sumet_;
	TBranch *pfclus_hfh_sumet_branch;
	bool pfclus_hfh_sumet_isLoaded;
	float pfclusmet_;
	TBranch *pfclusmet_branch;
	bool pfclusmet_isLoaded;
	float pfclussumet_;
	TBranch *pfclussumet_branch;
	bool pfclussumet_isLoaded;
	float pf_eb_met_;
	TBranch *pf_eb_met_branch;
	bool pf_eb_met_isLoaded;
	float pf_eb_sumet_;
	TBranch *pf_eb_sumet_branch;
	bool pf_eb_sumet_isLoaded;
	float pf_ee_met_;
	TBranch *pf_ee_met_branch;
	bool pf_ee_met_isLoaded;
	float pf_ee_sumet_;
	TBranch *pf_ee_sumet_branch;
	bool pf_ee_sumet_isLoaded;
	float pf_hb_met_;
	TBranch *pf_hb_met_branch;
	bool pf_hb_met_isLoaded;
	float pf_hb_sumet_;
	TBranch *pf_hb_sumet_branch;
	bool pf_hb_sumet_isLoaded;
	float pf_he_met_;
	TBranch *pf_he_met_branch;
	bool pf_he_met_isLoaded;
	float pf_he_sumet_;
	TBranch *pf_he_sumet_branch;
	bool pf_he_sumet_isLoaded;
	float pf_hfe_met_;
	TBranch *pf_hfe_met_branch;
	bool pf_hfe_met_isLoaded;
	float pf_hfe_sumet_;
	TBranch *pf_hfe_sumet_branch;
	bool pf_hfe_sumet_isLoaded;
	float pf_hfh_met_;
	TBranch *pf_hfh_met_branch;
	bool pf_hfh_met_isLoaded;
	float pf_hfh_sumet_;
	TBranch *pf_hfh_sumet_branch;
	bool pf_hfh_sumet_isLoaded;
	float pfmet_;
	TBranch *pfmet_branch;
	bool pfmet_isLoaded;
	float pfsumet_;
	TBranch *pfsumet_branch;
	bool pfsumet_isLoaded;
	vector<float> pf_cluster_detid_;
	TBranch *pf_cluster_detid_branch;
	bool pf_cluster_detid_isLoaded;
	vector<float> pf_cluster_e_;
	TBranch *pf_cluster_e_branch;
	bool pf_cluster_e_isLoaded;
	vector<float> pf_cluster_et_;
	TBranch *pf_cluster_et_branch;
	bool pf_cluster_et_isLoaded;
	vector<float> pf_cluster_eta_;
	TBranch *pf_cluster_eta_branch;
	bool pf_cluster_eta_isLoaded;
	vector<float> pf_cluster_phi_;
	TBranch *pf_cluster_phi_branch;
	bool pf_cluster_phi_isLoaded;
	vector<float> pf_ebcluster_e_;
	TBranch *pf_ebcluster_e_branch;
	bool pf_ebcluster_e_isLoaded;
	vector<float> pf_ebcluster_et_;
	TBranch *pf_ebcluster_et_branch;
	bool pf_ebcluster_et_isLoaded;
	vector<float> pf_ebcluster_eta_;
	TBranch *pf_ebcluster_eta_branch;
	bool pf_ebcluster_eta_isLoaded;
	vector<float> pf_ebcluster_phi_;
	TBranch *pf_ebcluster_phi_branch;
	bool pf_ebcluster_phi_isLoaded;
	vector<float> pf_ebrechit_e_;
	TBranch *pf_ebrechit_e_branch;
	bool pf_ebrechit_e_isLoaded;
	vector<float> pf_ebrechit_et_;
	TBranch *pf_ebrechit_et_branch;
	bool pf_ebrechit_et_isLoaded;
	vector<float> pf_ebrechit_eta_;
	TBranch *pf_ebrechit_eta_branch;
	bool pf_ebrechit_eta_isLoaded;
	vector<float> pf_ebrechit_phi_;
	TBranch *pf_ebrechit_phi_branch;
	bool pf_ebrechit_phi_isLoaded;
	vector<float> pf_eecluster_e_;
	TBranch *pf_eecluster_e_branch;
	bool pf_eecluster_e_isLoaded;
	vector<float> pf_eecluster_et_;
	TBranch *pf_eecluster_et_branch;
	bool pf_eecluster_et_isLoaded;
	vector<float> pf_eecluster_eta_;
	TBranch *pf_eecluster_eta_branch;
	bool pf_eecluster_eta_isLoaded;
	vector<float> pf_eecluster_phi_;
	TBranch *pf_eecluster_phi_branch;
	bool pf_eecluster_phi_isLoaded;
	vector<float> pf_eerechit_e_;
	TBranch *pf_eerechit_e_branch;
	bool pf_eerechit_e_isLoaded;
	vector<float> pf_eerechit_et_;
	TBranch *pf_eerechit_et_branch;
	bool pf_eerechit_et_isLoaded;
	vector<float> pf_eerechit_eta_;
	TBranch *pf_eerechit_eta_branch;
	bool pf_eerechit_eta_isLoaded;
	vector<float> pf_eerechit_phi_;
	TBranch *pf_eerechit_phi_branch;
	bool pf_eerechit_phi_isLoaded;
	vector<float> pf_hbcluster_e_;
	TBranch *pf_hbcluster_e_branch;
	bool pf_hbcluster_e_isLoaded;
	vector<float> pf_hbcluster_et_;
	TBranch *pf_hbcluster_et_branch;
	bool pf_hbcluster_et_isLoaded;
	vector<float> pf_hbcluster_eta_;
	TBranch *pf_hbcluster_eta_branch;
	bool pf_hbcluster_eta_isLoaded;
	vector<float> pf_hbcluster_phi_;
	TBranch *pf_hbcluster_phi_branch;
	bool pf_hbcluster_phi_isLoaded;
	vector<float> pf_hbrechit_e_;
	TBranch *pf_hbrechit_e_branch;
	bool pf_hbrechit_e_isLoaded;
	vector<float> pf_hbrechit_et_;
	TBranch *pf_hbrechit_et_branch;
	bool pf_hbrechit_et_isLoaded;
	vector<float> pf_hbrechit_eta_;
	TBranch *pf_hbrechit_eta_branch;
	bool pf_hbrechit_eta_isLoaded;
	vector<float> pf_hbrechit_phi_;
	TBranch *pf_hbrechit_phi_branch;
	bool pf_hbrechit_phi_isLoaded;
	vector<float> pf_hecluster_e_;
	TBranch *pf_hecluster_e_branch;
	bool pf_hecluster_e_isLoaded;
	vector<float> pf_hecluster_et_;
	TBranch *pf_hecluster_et_branch;
	bool pf_hecluster_et_isLoaded;
	vector<float> pf_hecluster_eta_;
	TBranch *pf_hecluster_eta_branch;
	bool pf_hecluster_eta_isLoaded;
	vector<float> pf_hecluster_phi_;
	TBranch *pf_hecluster_phi_branch;
	bool pf_hecluster_phi_isLoaded;
	vector<float> pf_herechit_e_;
	TBranch *pf_herechit_e_branch;
	bool pf_herechit_e_isLoaded;
	vector<float> pf_herechit_et_;
	TBranch *pf_herechit_et_branch;
	bool pf_herechit_et_isLoaded;
	vector<float> pf_herechit_eta_;
	TBranch *pf_herechit_eta_branch;
	bool pf_herechit_eta_isLoaded;
	vector<float> pf_herechit_phi_;
	TBranch *pf_herechit_phi_branch;
	bool pf_herechit_phi_isLoaded;
	vector<float> pf_hfecluster_e_;
	TBranch *pf_hfecluster_e_branch;
	bool pf_hfecluster_e_isLoaded;
	vector<float> pf_hfecluster_et_;
	TBranch *pf_hfecluster_et_branch;
	bool pf_hfecluster_et_isLoaded;
	vector<float> pf_hfecluster_eta_;
	TBranch *pf_hfecluster_eta_branch;
	bool pf_hfecluster_eta_isLoaded;
	vector<float> pf_hfecluster_phi_;
	TBranch *pf_hfecluster_phi_branch;
	bool pf_hfecluster_phi_isLoaded;
	vector<float> pf_hferechit_e_;
	TBranch *pf_hferechit_e_branch;
	bool pf_hferechit_e_isLoaded;
	vector<float> pf_hferechit_et_;
	TBranch *pf_hferechit_et_branch;
	bool pf_hferechit_et_isLoaded;
	vector<float> pf_hferechit_eta_;
	TBranch *pf_hferechit_eta_branch;
	bool pf_hferechit_eta_isLoaded;
	vector<float> pf_hferechit_phi_;
	TBranch *pf_hferechit_phi_branch;
	bool pf_hferechit_phi_isLoaded;
	vector<float> pf_hfhcluster_e_;
	TBranch *pf_hfhcluster_e_branch;
	bool pf_hfhcluster_e_isLoaded;
	vector<float> pf_hfhcluster_et_;
	TBranch *pf_hfhcluster_et_branch;
	bool pf_hfhcluster_et_isLoaded;
	vector<float> pf_hfhcluster_eta_;
	TBranch *pf_hfhcluster_eta_branch;
	bool pf_hfhcluster_eta_isLoaded;
	vector<float> pf_hfhcluster_phi_;
	TBranch *pf_hfhcluster_phi_branch;
	bool pf_hfhcluster_phi_isLoaded;
	vector<float> pf_hfhrechit_e_;
	TBranch *pf_hfhrechit_e_branch;
	bool pf_hfhrechit_e_isLoaded;
	vector<float> pf_hfhrechit_et_;
	TBranch *pf_hfhrechit_et_branch;
	bool pf_hfhrechit_et_isLoaded;
	vector<float> pf_hfhrechit_eta_;
	TBranch *pf_hfhrechit_eta_branch;
	bool pf_hfhrechit_eta_isLoaded;
	vector<float> pf_hfhrechit_phi_;
	TBranch *pf_hfhrechit_phi_branch;
	bool pf_hfhrechit_phi_isLoaded;
	vector<float> pf_rechit_detid_;
	TBranch *pf_rechit_detid_branch;
	bool pf_rechit_detid_isLoaded;
	vector<float> pf_rechit_e_;
	TBranch *pf_rechit_e_branch;
	bool pf_rechit_e_isLoaded;
	vector<float> pf_rechit_et_;
	TBranch *pf_rechit_et_branch;
	bool pf_rechit_et_isLoaded;
	vector<float> pf_rechit_eta_;
	TBranch *pf_rechit_eta_branch;
	bool pf_rechit_eta_isLoaded;
	vector<float> pf_rechit_phi_;
	TBranch *pf_rechit_phi_branch;
	bool pf_rechit_phi_isLoaded;
public: 
void Init(TTree *tree) {
  tree->SetMakeClass(1);
	calo_eb_sumet_branch = 0;
	if (tree->GetAlias("calo_eb_sumet") != 0) {
		calo_eb_sumet_branch = tree->GetBranch(tree->GetAlias("calo_eb_sumet"));
		calo_eb_sumet_branch->SetAddress(&calo_eb_sumet_);
	}
	if(calo_eb_sumet_branch == 0 ) {
	cout << "Branch calo_eb_sumet does not exist." << endl;
	}
	calo_ee_sumet_branch = 0;
	if (tree->GetAlias("calo_ee_sumet") != 0) {
		calo_ee_sumet_branch = tree->GetBranch(tree->GetAlias("calo_ee_sumet"));
		calo_ee_sumet_branch->SetAddress(&calo_ee_sumet_);
	}
	if(calo_ee_sumet_branch == 0 ) {
	cout << "Branch calo_ee_sumet does not exist." << endl;
	}
	calo_hb_sumet_branch = 0;
	if (tree->GetAlias("calo_hb_sumet") != 0) {
		calo_hb_sumet_branch = tree->GetBranch(tree->GetAlias("calo_hb_sumet"));
		calo_hb_sumet_branch->SetAddress(&calo_hb_sumet_);
	}
	if(calo_hb_sumet_branch == 0 ) {
	cout << "Branch calo_hb_sumet does not exist." << endl;
	}
	calo_he_sumet_branch = 0;
	if (tree->GetAlias("calo_he_sumet") != 0) {
		calo_he_sumet_branch = tree->GetBranch(tree->GetAlias("calo_he_sumet"));
		calo_he_sumet_branch->SetAddress(&calo_he_sumet_);
	}
	if(calo_he_sumet_branch == 0 ) {
	cout << "Branch calo_he_sumet does not exist." << endl;
	}
	calo_hfe_sumet_branch = 0;
	if (tree->GetAlias("calo_hfe_sumet") != 0) {
		calo_hfe_sumet_branch = tree->GetBranch(tree->GetAlias("calo_hfe_sumet"));
		calo_hfe_sumet_branch->SetAddress(&calo_hfe_sumet_);
	}
	if(calo_hfe_sumet_branch == 0 ) {
	cout << "Branch calo_hfe_sumet does not exist." << endl;
	}
	calo_hfh_sumet_branch = 0;
	if (tree->GetAlias("calo_hfh_sumet") != 0) {
		calo_hfh_sumet_branch = tree->GetBranch(tree->GetAlias("calo_hfh_sumet"));
		calo_hfh_sumet_branch->SetAddress(&calo_hfh_sumet_);
	}
	if(calo_hfh_sumet_branch == 0 ) {
	cout << "Branch calo_hfh_sumet does not exist." << endl;
	}
	calomet_branch = 0;
	if (tree->GetAlias("calomet") != 0) {
		calomet_branch = tree->GetBranch(tree->GetAlias("calomet"));
		calomet_branch->SetAddress(&calomet_);
	}
	if(calomet_branch == 0 ) {
	cout << "Branch calomet does not exist." << endl;
	}
	calosumet_branch = 0;
	if (tree->GetAlias("calosumet") != 0) {
		calosumet_branch = tree->GetBranch(tree->GetAlias("calosumet"));
		calosumet_branch->SetAddress(&calosumet_);
	}
	if(calosumet_branch == 0 ) {
	cout << "Branch calosumet does not exist." << endl;
	}
	genmet_branch = 0;
	if (tree->GetAlias("genmet") != 0) {
		genmet_branch = tree->GetBranch(tree->GetAlias("genmet"));
		genmet_branch->SetAddress(&genmet_);
	}
	if(genmet_branch == 0 ) {
	cout << "Branch genmet does not exist." << endl;
	}
	gensumet_branch = 0;
	if (tree->GetAlias("gensumet") != 0) {
		gensumet_branch = tree->GetBranch(tree->GetAlias("gensumet"));
		gensumet_branch->SetAddress(&gensumet_);
	}
	if(gensumet_branch == 0 ) {
	cout << "Branch gensumet does not exist." << endl;
	}
	pfclus_eb_met_branch = 0;
	if (tree->GetAlias("pfclus_eb_met") != 0) {
		pfclus_eb_met_branch = tree->GetBranch(tree->GetAlias("pfclus_eb_met"));
		pfclus_eb_met_branch->SetAddress(&pfclus_eb_met_);
	}
	if(pfclus_eb_met_branch == 0 ) {
	cout << "Branch pfclus_eb_met does not exist." << endl;
	}
	pfclus_eb_sumet_branch = 0;
	if (tree->GetAlias("pfclus_eb_sumet") != 0) {
		pfclus_eb_sumet_branch = tree->GetBranch(tree->GetAlias("pfclus_eb_sumet"));
		pfclus_eb_sumet_branch->SetAddress(&pfclus_eb_sumet_);
	}
	if(pfclus_eb_sumet_branch == 0 ) {
	cout << "Branch pfclus_eb_sumet does not exist." << endl;
	}
	pfclus_ee_met_branch = 0;
	if (tree->GetAlias("pfclus_ee_met") != 0) {
		pfclus_ee_met_branch = tree->GetBranch(tree->GetAlias("pfclus_ee_met"));
		pfclus_ee_met_branch->SetAddress(&pfclus_ee_met_);
	}
	if(pfclus_ee_met_branch == 0 ) {
	cout << "Branch pfclus_ee_met does not exist." << endl;
	}
	pfclus_ee_sumet_branch = 0;
	if (tree->GetAlias("pfclus_ee_sumet") != 0) {
		pfclus_ee_sumet_branch = tree->GetBranch(tree->GetAlias("pfclus_ee_sumet"));
		pfclus_ee_sumet_branch->SetAddress(&pfclus_ee_sumet_);
	}
	if(pfclus_ee_sumet_branch == 0 ) {
	cout << "Branch pfclus_ee_sumet does not exist." << endl;
	}
	pfclus_hb_met_branch = 0;
	if (tree->GetAlias("pfclus_hb_met") != 0) {
		pfclus_hb_met_branch = tree->GetBranch(tree->GetAlias("pfclus_hb_met"));
		pfclus_hb_met_branch->SetAddress(&pfclus_hb_met_);
	}
	if(pfclus_hb_met_branch == 0 ) {
	cout << "Branch pfclus_hb_met does not exist." << endl;
	}
	pfclus_hb_sumet_branch = 0;
	if (tree->GetAlias("pfclus_hb_sumet") != 0) {
		pfclus_hb_sumet_branch = tree->GetBranch(tree->GetAlias("pfclus_hb_sumet"));
		pfclus_hb_sumet_branch->SetAddress(&pfclus_hb_sumet_);
	}
	if(pfclus_hb_sumet_branch == 0 ) {
	cout << "Branch pfclus_hb_sumet does not exist." << endl;
	}
	pfclus_he_met_branch = 0;
	if (tree->GetAlias("pfclus_he_met") != 0) {
		pfclus_he_met_branch = tree->GetBranch(tree->GetAlias("pfclus_he_met"));
		pfclus_he_met_branch->SetAddress(&pfclus_he_met_);
	}
	if(pfclus_he_met_branch == 0 ) {
	cout << "Branch pfclus_he_met does not exist." << endl;
	}
	pfclus_he_sumet_branch = 0;
	if (tree->GetAlias("pfclus_he_sumet") != 0) {
		pfclus_he_sumet_branch = tree->GetBranch(tree->GetAlias("pfclus_he_sumet"));
		pfclus_he_sumet_branch->SetAddress(&pfclus_he_sumet_);
	}
	if(pfclus_he_sumet_branch == 0 ) {
	cout << "Branch pfclus_he_sumet does not exist." << endl;
	}
	pfclus_hfe_met_branch = 0;
	if (tree->GetAlias("pfclus_hfe_met") != 0) {
		pfclus_hfe_met_branch = tree->GetBranch(tree->GetAlias("pfclus_hfe_met"));
		pfclus_hfe_met_branch->SetAddress(&pfclus_hfe_met_);
	}
	if(pfclus_hfe_met_branch == 0 ) {
	cout << "Branch pfclus_hfe_met does not exist." << endl;
	}
	pfclus_hfe_sumet_branch = 0;
	if (tree->GetAlias("pfclus_hfe_sumet") != 0) {
		pfclus_hfe_sumet_branch = tree->GetBranch(tree->GetAlias("pfclus_hfe_sumet"));
		pfclus_hfe_sumet_branch->SetAddress(&pfclus_hfe_sumet_);
	}
	if(pfclus_hfe_sumet_branch == 0 ) {
	cout << "Branch pfclus_hfe_sumet does not exist." << endl;
	}
	pfclus_hfh_met_branch = 0;
	if (tree->GetAlias("pfclus_hfh_met") != 0) {
		pfclus_hfh_met_branch = tree->GetBranch(tree->GetAlias("pfclus_hfh_met"));
		pfclus_hfh_met_branch->SetAddress(&pfclus_hfh_met_);
	}
	if(pfclus_hfh_met_branch == 0 ) {
	cout << "Branch pfclus_hfh_met does not exist." << endl;
	}
	pfclus_hfh_sumet_branch = 0;
	if (tree->GetAlias("pfclus_hfh_sumet") != 0) {
		pfclus_hfh_sumet_branch = tree->GetBranch(tree->GetAlias("pfclus_hfh_sumet"));
		pfclus_hfh_sumet_branch->SetAddress(&pfclus_hfh_sumet_);
	}
	if(pfclus_hfh_sumet_branch == 0 ) {
	cout << "Branch pfclus_hfh_sumet does not exist." << endl;
	}
	pfclusmet_branch = 0;
	if (tree->GetAlias("pfclusmet") != 0) {
		pfclusmet_branch = tree->GetBranch(tree->GetAlias("pfclusmet"));
		pfclusmet_branch->SetAddress(&pfclusmet_);
	}
	if(pfclusmet_branch == 0 ) {
	cout << "Branch pfclusmet does not exist." << endl;
	}
	pfclussumet_branch = 0;
	if (tree->GetAlias("pfclussumet") != 0) {
		pfclussumet_branch = tree->GetBranch(tree->GetAlias("pfclussumet"));
		pfclussumet_branch->SetAddress(&pfclussumet_);
	}
	if(pfclussumet_branch == 0 ) {
	cout << "Branch pfclussumet does not exist." << endl;
	}
	pf_eb_met_branch = 0;
	if (tree->GetAlias("pf_eb_met") != 0) {
		pf_eb_met_branch = tree->GetBranch(tree->GetAlias("pf_eb_met"));
		pf_eb_met_branch->SetAddress(&pf_eb_met_);
	}
	if(pf_eb_met_branch == 0 ) {
	cout << "Branch pf_eb_met does not exist." << endl;
	}
	pf_eb_sumet_branch = 0;
	if (tree->GetAlias("pf_eb_sumet") != 0) {
		pf_eb_sumet_branch = tree->GetBranch(tree->GetAlias("pf_eb_sumet"));
		pf_eb_sumet_branch->SetAddress(&pf_eb_sumet_);
	}
	if(pf_eb_sumet_branch == 0 ) {
	cout << "Branch pf_eb_sumet does not exist." << endl;
	}
	pf_ee_met_branch = 0;
	if (tree->GetAlias("pf_ee_met") != 0) {
		pf_ee_met_branch = tree->GetBranch(tree->GetAlias("pf_ee_met"));
		pf_ee_met_branch->SetAddress(&pf_ee_met_);
	}
	if(pf_ee_met_branch == 0 ) {
	cout << "Branch pf_ee_met does not exist." << endl;
	}
	pf_ee_sumet_branch = 0;
	if (tree->GetAlias("pf_ee_sumet") != 0) {
		pf_ee_sumet_branch = tree->GetBranch(tree->GetAlias("pf_ee_sumet"));
		pf_ee_sumet_branch->SetAddress(&pf_ee_sumet_);
	}
	if(pf_ee_sumet_branch == 0 ) {
	cout << "Branch pf_ee_sumet does not exist." << endl;
	}
	pf_hb_met_branch = 0;
	if (tree->GetAlias("pf_hb_met") != 0) {
		pf_hb_met_branch = tree->GetBranch(tree->GetAlias("pf_hb_met"));
		pf_hb_met_branch->SetAddress(&pf_hb_met_);
	}
	if(pf_hb_met_branch == 0 ) {
	cout << "Branch pf_hb_met does not exist." << endl;
	}
	pf_hb_sumet_branch = 0;
	if (tree->GetAlias("pf_hb_sumet") != 0) {
		pf_hb_sumet_branch = tree->GetBranch(tree->GetAlias("pf_hb_sumet"));
		pf_hb_sumet_branch->SetAddress(&pf_hb_sumet_);
	}
	if(pf_hb_sumet_branch == 0 ) {
	cout << "Branch pf_hb_sumet does not exist." << endl;
	}
	pf_he_met_branch = 0;
	if (tree->GetAlias("pf_he_met") != 0) {
		pf_he_met_branch = tree->GetBranch(tree->GetAlias("pf_he_met"));
		pf_he_met_branch->SetAddress(&pf_he_met_);
	}
	if(pf_he_met_branch == 0 ) {
	cout << "Branch pf_he_met does not exist." << endl;
	}
	pf_he_sumet_branch = 0;
	if (tree->GetAlias("pf_he_sumet") != 0) {
		pf_he_sumet_branch = tree->GetBranch(tree->GetAlias("pf_he_sumet"));
		pf_he_sumet_branch->SetAddress(&pf_he_sumet_);
	}
	if(pf_he_sumet_branch == 0 ) {
	cout << "Branch pf_he_sumet does not exist." << endl;
	}
	pf_hfe_met_branch = 0;
	if (tree->GetAlias("pf_hfe_met") != 0) {
		pf_hfe_met_branch = tree->GetBranch(tree->GetAlias("pf_hfe_met"));
		pf_hfe_met_branch->SetAddress(&pf_hfe_met_);
	}
	if(pf_hfe_met_branch == 0 ) {
	cout << "Branch pf_hfe_met does not exist." << endl;
	}
	pf_hfe_sumet_branch = 0;
	if (tree->GetAlias("pf_hfe_sumet") != 0) {
		pf_hfe_sumet_branch = tree->GetBranch(tree->GetAlias("pf_hfe_sumet"));
		pf_hfe_sumet_branch->SetAddress(&pf_hfe_sumet_);
	}
	if(pf_hfe_sumet_branch == 0 ) {
	cout << "Branch pf_hfe_sumet does not exist." << endl;
	}
	pf_hfh_met_branch = 0;
	if (tree->GetAlias("pf_hfh_met") != 0) {
		pf_hfh_met_branch = tree->GetBranch(tree->GetAlias("pf_hfh_met"));
		pf_hfh_met_branch->SetAddress(&pf_hfh_met_);
	}
	if(pf_hfh_met_branch == 0 ) {
	cout << "Branch pf_hfh_met does not exist." << endl;
	}
	pf_hfh_sumet_branch = 0;
	if (tree->GetAlias("pf_hfh_sumet") != 0) {
		pf_hfh_sumet_branch = tree->GetBranch(tree->GetAlias("pf_hfh_sumet"));
		pf_hfh_sumet_branch->SetAddress(&pf_hfh_sumet_);
	}
	if(pf_hfh_sumet_branch == 0 ) {
	cout << "Branch pf_hfh_sumet does not exist." << endl;
	}
	pfmet_branch = 0;
	if (tree->GetAlias("pfmet") != 0) {
		pfmet_branch = tree->GetBranch(tree->GetAlias("pfmet"));
		pfmet_branch->SetAddress(&pfmet_);
	}
	if(pfmet_branch == 0 ) {
	cout << "Branch pfmet does not exist." << endl;
	}
	pfsumet_branch = 0;
	if (tree->GetAlias("pfsumet") != 0) {
		pfsumet_branch = tree->GetBranch(tree->GetAlias("pfsumet"));
		pfsumet_branch->SetAddress(&pfsumet_);
	}
	if(pfsumet_branch == 0 ) {
	cout << "Branch pfsumet does not exist." << endl;
	}
	pf_cluster_detid_branch = 0;
	if (tree->GetAlias("pf_cluster_detid") != 0) {
		pf_cluster_detid_branch = tree->GetBranch(tree->GetAlias("pf_cluster_detid"));
		pf_cluster_detid_branch->SetAddress(&pf_cluster_detid_);
	}
	if(pf_cluster_detid_branch == 0 ) {
	cout << "Branch pf_cluster_detid does not exist." << endl;
	}
	pf_cluster_e_branch = 0;
	if (tree->GetAlias("pf_cluster_e") != 0) {
		pf_cluster_e_branch = tree->GetBranch(tree->GetAlias("pf_cluster_e"));
		pf_cluster_e_branch->SetAddress(&pf_cluster_e_);
	}
	if(pf_cluster_e_branch == 0 ) {
	cout << "Branch pf_cluster_e does not exist." << endl;
	}
	pf_cluster_et_branch = 0;
	if (tree->GetAlias("pf_cluster_et") != 0) {
		pf_cluster_et_branch = tree->GetBranch(tree->GetAlias("pf_cluster_et"));
		pf_cluster_et_branch->SetAddress(&pf_cluster_et_);
	}
	if(pf_cluster_et_branch == 0 ) {
	cout << "Branch pf_cluster_et does not exist." << endl;
	}
	pf_cluster_eta_branch = 0;
	if (tree->GetAlias("pf_cluster_eta") != 0) {
		pf_cluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_cluster_eta"));
		pf_cluster_eta_branch->SetAddress(&pf_cluster_eta_);
	}
	if(pf_cluster_eta_branch == 0 ) {
	cout << "Branch pf_cluster_eta does not exist." << endl;
	}
	pf_cluster_phi_branch = 0;
	if (tree->GetAlias("pf_cluster_phi") != 0) {
		pf_cluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_cluster_phi"));
		pf_cluster_phi_branch->SetAddress(&pf_cluster_phi_);
	}
	if(pf_cluster_phi_branch == 0 ) {
	cout << "Branch pf_cluster_phi does not exist." << endl;
	}
	pf_ebcluster_e_branch = 0;
	if (tree->GetAlias("pf_ebcluster_e") != 0) {
		pf_ebcluster_e_branch = tree->GetBranch(tree->GetAlias("pf_ebcluster_e"));
		pf_ebcluster_e_branch->SetAddress(&pf_ebcluster_e_);
	}
	if(pf_ebcluster_e_branch == 0 ) {
	cout << "Branch pf_ebcluster_e does not exist." << endl;
	}
	pf_ebcluster_et_branch = 0;
	if (tree->GetAlias("pf_ebcluster_et") != 0) {
		pf_ebcluster_et_branch = tree->GetBranch(tree->GetAlias("pf_ebcluster_et"));
		pf_ebcluster_et_branch->SetAddress(&pf_ebcluster_et_);
	}
	if(pf_ebcluster_et_branch == 0 ) {
	cout << "Branch pf_ebcluster_et does not exist." << endl;
	}
	pf_ebcluster_eta_branch = 0;
	if (tree->GetAlias("pf_ebcluster_eta") != 0) {
		pf_ebcluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_ebcluster_eta"));
		pf_ebcluster_eta_branch->SetAddress(&pf_ebcluster_eta_);
	}
	if(pf_ebcluster_eta_branch == 0 ) {
	cout << "Branch pf_ebcluster_eta does not exist." << endl;
	}
	pf_ebcluster_phi_branch = 0;
	if (tree->GetAlias("pf_ebcluster_phi") != 0) {
		pf_ebcluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_ebcluster_phi"));
		pf_ebcluster_phi_branch->SetAddress(&pf_ebcluster_phi_);
	}
	if(pf_ebcluster_phi_branch == 0 ) {
	cout << "Branch pf_ebcluster_phi does not exist." << endl;
	}
	pf_ebrechit_e_branch = 0;
	if (tree->GetAlias("pf_ebrechit_e") != 0) {
		pf_ebrechit_e_branch = tree->GetBranch(tree->GetAlias("pf_ebrechit_e"));
		pf_ebrechit_e_branch->SetAddress(&pf_ebrechit_e_);
	}
	if(pf_ebrechit_e_branch == 0 ) {
	cout << "Branch pf_ebrechit_e does not exist." << endl;
	}
	pf_ebrechit_et_branch = 0;
	if (tree->GetAlias("pf_ebrechit_et") != 0) {
		pf_ebrechit_et_branch = tree->GetBranch(tree->GetAlias("pf_ebrechit_et"));
		pf_ebrechit_et_branch->SetAddress(&pf_ebrechit_et_);
	}
	if(pf_ebrechit_et_branch == 0 ) {
	cout << "Branch pf_ebrechit_et does not exist." << endl;
	}
	pf_ebrechit_eta_branch = 0;
	if (tree->GetAlias("pf_ebrechit_eta") != 0) {
		pf_ebrechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_ebrechit_eta"));
		pf_ebrechit_eta_branch->SetAddress(&pf_ebrechit_eta_);
	}
	if(pf_ebrechit_eta_branch == 0 ) {
	cout << "Branch pf_ebrechit_eta does not exist." << endl;
	}
	pf_ebrechit_phi_branch = 0;
	if (tree->GetAlias("pf_ebrechit_phi") != 0) {
		pf_ebrechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_ebrechit_phi"));
		pf_ebrechit_phi_branch->SetAddress(&pf_ebrechit_phi_);
	}
	if(pf_ebrechit_phi_branch == 0 ) {
	cout << "Branch pf_ebrechit_phi does not exist." << endl;
	}
	pf_eecluster_e_branch = 0;
	if (tree->GetAlias("pf_eecluster_e") != 0) {
		pf_eecluster_e_branch = tree->GetBranch(tree->GetAlias("pf_eecluster_e"));
		pf_eecluster_e_branch->SetAddress(&pf_eecluster_e_);
	}
	if(pf_eecluster_e_branch == 0 ) {
	cout << "Branch pf_eecluster_e does not exist." << endl;
	}
	pf_eecluster_et_branch = 0;
	if (tree->GetAlias("pf_eecluster_et") != 0) {
		pf_eecluster_et_branch = tree->GetBranch(tree->GetAlias("pf_eecluster_et"));
		pf_eecluster_et_branch->SetAddress(&pf_eecluster_et_);
	}
	if(pf_eecluster_et_branch == 0 ) {
	cout << "Branch pf_eecluster_et does not exist." << endl;
	}
	pf_eecluster_eta_branch = 0;
	if (tree->GetAlias("pf_eecluster_eta") != 0) {
		pf_eecluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_eecluster_eta"));
		pf_eecluster_eta_branch->SetAddress(&pf_eecluster_eta_);
	}
	if(pf_eecluster_eta_branch == 0 ) {
	cout << "Branch pf_eecluster_eta does not exist." << endl;
	}
	pf_eecluster_phi_branch = 0;
	if (tree->GetAlias("pf_eecluster_phi") != 0) {
		pf_eecluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_eecluster_phi"));
		pf_eecluster_phi_branch->SetAddress(&pf_eecluster_phi_);
	}
	if(pf_eecluster_phi_branch == 0 ) {
	cout << "Branch pf_eecluster_phi does not exist." << endl;
	}
	pf_eerechit_e_branch = 0;
	if (tree->GetAlias("pf_eerechit_e") != 0) {
		pf_eerechit_e_branch = tree->GetBranch(tree->GetAlias("pf_eerechit_e"));
		pf_eerechit_e_branch->SetAddress(&pf_eerechit_e_);
	}
	if(pf_eerechit_e_branch == 0 ) {
	cout << "Branch pf_eerechit_e does not exist." << endl;
	}
	pf_eerechit_et_branch = 0;
	if (tree->GetAlias("pf_eerechit_et") != 0) {
		pf_eerechit_et_branch = tree->GetBranch(tree->GetAlias("pf_eerechit_et"));
		pf_eerechit_et_branch->SetAddress(&pf_eerechit_et_);
	}
	if(pf_eerechit_et_branch == 0 ) {
	cout << "Branch pf_eerechit_et does not exist." << endl;
	}
	pf_eerechit_eta_branch = 0;
	if (tree->GetAlias("pf_eerechit_eta") != 0) {
		pf_eerechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_eerechit_eta"));
		pf_eerechit_eta_branch->SetAddress(&pf_eerechit_eta_);
	}
	if(pf_eerechit_eta_branch == 0 ) {
	cout << "Branch pf_eerechit_eta does not exist." << endl;
	}
	pf_eerechit_phi_branch = 0;
	if (tree->GetAlias("pf_eerechit_phi") != 0) {
		pf_eerechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_eerechit_phi"));
		pf_eerechit_phi_branch->SetAddress(&pf_eerechit_phi_);
	}
	if(pf_eerechit_phi_branch == 0 ) {
	cout << "Branch pf_eerechit_phi does not exist." << endl;
	}
	pf_hbcluster_e_branch = 0;
	if (tree->GetAlias("pf_hbcluster_e") != 0) {
		pf_hbcluster_e_branch = tree->GetBranch(tree->GetAlias("pf_hbcluster_e"));
		pf_hbcluster_e_branch->SetAddress(&pf_hbcluster_e_);
	}
	if(pf_hbcluster_e_branch == 0 ) {
	cout << "Branch pf_hbcluster_e does not exist." << endl;
	}
	pf_hbcluster_et_branch = 0;
	if (tree->GetAlias("pf_hbcluster_et") != 0) {
		pf_hbcluster_et_branch = tree->GetBranch(tree->GetAlias("pf_hbcluster_et"));
		pf_hbcluster_et_branch->SetAddress(&pf_hbcluster_et_);
	}
	if(pf_hbcluster_et_branch == 0 ) {
	cout << "Branch pf_hbcluster_et does not exist." << endl;
	}
	pf_hbcluster_eta_branch = 0;
	if (tree->GetAlias("pf_hbcluster_eta") != 0) {
		pf_hbcluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_hbcluster_eta"));
		pf_hbcluster_eta_branch->SetAddress(&pf_hbcluster_eta_);
	}
	if(pf_hbcluster_eta_branch == 0 ) {
	cout << "Branch pf_hbcluster_eta does not exist." << endl;
	}
	pf_hbcluster_phi_branch = 0;
	if (tree->GetAlias("pf_hbcluster_phi") != 0) {
		pf_hbcluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_hbcluster_phi"));
		pf_hbcluster_phi_branch->SetAddress(&pf_hbcluster_phi_);
	}
	if(pf_hbcluster_phi_branch == 0 ) {
	cout << "Branch pf_hbcluster_phi does not exist." << endl;
	}
	pf_hbrechit_e_branch = 0;
	if (tree->GetAlias("pf_hbrechit_e") != 0) {
		pf_hbrechit_e_branch = tree->GetBranch(tree->GetAlias("pf_hbrechit_e"));
		pf_hbrechit_e_branch->SetAddress(&pf_hbrechit_e_);
	}
	if(pf_hbrechit_e_branch == 0 ) {
	cout << "Branch pf_hbrechit_e does not exist." << endl;
	}
	pf_hbrechit_et_branch = 0;
	if (tree->GetAlias("pf_hbrechit_et") != 0) {
		pf_hbrechit_et_branch = tree->GetBranch(tree->GetAlias("pf_hbrechit_et"));
		pf_hbrechit_et_branch->SetAddress(&pf_hbrechit_et_);
	}
	if(pf_hbrechit_et_branch == 0 ) {
	cout << "Branch pf_hbrechit_et does not exist." << endl;
	}
	pf_hbrechit_eta_branch = 0;
	if (tree->GetAlias("pf_hbrechit_eta") != 0) {
		pf_hbrechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_hbrechit_eta"));
		pf_hbrechit_eta_branch->SetAddress(&pf_hbrechit_eta_);
	}
	if(pf_hbrechit_eta_branch == 0 ) {
	cout << "Branch pf_hbrechit_eta does not exist." << endl;
	}
	pf_hbrechit_phi_branch = 0;
	if (tree->GetAlias("pf_hbrechit_phi") != 0) {
		pf_hbrechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_hbrechit_phi"));
		pf_hbrechit_phi_branch->SetAddress(&pf_hbrechit_phi_);
	}
	if(pf_hbrechit_phi_branch == 0 ) {
	cout << "Branch pf_hbrechit_phi does not exist." << endl;
	}
	pf_hecluster_e_branch = 0;
	if (tree->GetAlias("pf_hecluster_e") != 0) {
		pf_hecluster_e_branch = tree->GetBranch(tree->GetAlias("pf_hecluster_e"));
		pf_hecluster_e_branch->SetAddress(&pf_hecluster_e_);
	}
	if(pf_hecluster_e_branch == 0 ) {
	cout << "Branch pf_hecluster_e does not exist." << endl;
	}
	pf_hecluster_et_branch = 0;
	if (tree->GetAlias("pf_hecluster_et") != 0) {
		pf_hecluster_et_branch = tree->GetBranch(tree->GetAlias("pf_hecluster_et"));
		pf_hecluster_et_branch->SetAddress(&pf_hecluster_et_);
	}
	if(pf_hecluster_et_branch == 0 ) {
	cout << "Branch pf_hecluster_et does not exist." << endl;
	}
	pf_hecluster_eta_branch = 0;
	if (tree->GetAlias("pf_hecluster_eta") != 0) {
		pf_hecluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_hecluster_eta"));
		pf_hecluster_eta_branch->SetAddress(&pf_hecluster_eta_);
	}
	if(pf_hecluster_eta_branch == 0 ) {
	cout << "Branch pf_hecluster_eta does not exist." << endl;
	}
	pf_hecluster_phi_branch = 0;
	if (tree->GetAlias("pf_hecluster_phi") != 0) {
		pf_hecluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_hecluster_phi"));
		pf_hecluster_phi_branch->SetAddress(&pf_hecluster_phi_);
	}
	if(pf_hecluster_phi_branch == 0 ) {
	cout << "Branch pf_hecluster_phi does not exist." << endl;
	}
	pf_herechit_e_branch = 0;
	if (tree->GetAlias("pf_herechit_e") != 0) {
		pf_herechit_e_branch = tree->GetBranch(tree->GetAlias("pf_herechit_e"));
		pf_herechit_e_branch->SetAddress(&pf_herechit_e_);
	}
	if(pf_herechit_e_branch == 0 ) {
	cout << "Branch pf_herechit_e does not exist." << endl;
	}
	pf_herechit_et_branch = 0;
	if (tree->GetAlias("pf_herechit_et") != 0) {
		pf_herechit_et_branch = tree->GetBranch(tree->GetAlias("pf_herechit_et"));
		pf_herechit_et_branch->SetAddress(&pf_herechit_et_);
	}
	if(pf_herechit_et_branch == 0 ) {
	cout << "Branch pf_herechit_et does not exist." << endl;
	}
	pf_herechit_eta_branch = 0;
	if (tree->GetAlias("pf_herechit_eta") != 0) {
		pf_herechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_herechit_eta"));
		pf_herechit_eta_branch->SetAddress(&pf_herechit_eta_);
	}
	if(pf_herechit_eta_branch == 0 ) {
	cout << "Branch pf_herechit_eta does not exist." << endl;
	}
	pf_herechit_phi_branch = 0;
	if (tree->GetAlias("pf_herechit_phi") != 0) {
		pf_herechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_herechit_phi"));
		pf_herechit_phi_branch->SetAddress(&pf_herechit_phi_);
	}
	if(pf_herechit_phi_branch == 0 ) {
	cout << "Branch pf_herechit_phi does not exist." << endl;
	}
	pf_hfecluster_e_branch = 0;
	if (tree->GetAlias("pf_hfecluster_e") != 0) {
		pf_hfecluster_e_branch = tree->GetBranch(tree->GetAlias("pf_hfecluster_e"));
		pf_hfecluster_e_branch->SetAddress(&pf_hfecluster_e_);
	}
	if(pf_hfecluster_e_branch == 0 ) {
	cout << "Branch pf_hfecluster_e does not exist." << endl;
	}
	pf_hfecluster_et_branch = 0;
	if (tree->GetAlias("pf_hfecluster_et") != 0) {
		pf_hfecluster_et_branch = tree->GetBranch(tree->GetAlias("pf_hfecluster_et"));
		pf_hfecluster_et_branch->SetAddress(&pf_hfecluster_et_);
	}
	if(pf_hfecluster_et_branch == 0 ) {
	cout << "Branch pf_hfecluster_et does not exist." << endl;
	}
	pf_hfecluster_eta_branch = 0;
	if (tree->GetAlias("pf_hfecluster_eta") != 0) {
		pf_hfecluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_hfecluster_eta"));
		pf_hfecluster_eta_branch->SetAddress(&pf_hfecluster_eta_);
	}
	if(pf_hfecluster_eta_branch == 0 ) {
	cout << "Branch pf_hfecluster_eta does not exist." << endl;
	}
	pf_hfecluster_phi_branch = 0;
	if (tree->GetAlias("pf_hfecluster_phi") != 0) {
		pf_hfecluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_hfecluster_phi"));
		pf_hfecluster_phi_branch->SetAddress(&pf_hfecluster_phi_);
	}
	if(pf_hfecluster_phi_branch == 0 ) {
	cout << "Branch pf_hfecluster_phi does not exist." << endl;
	}
	pf_hferechit_e_branch = 0;
	if (tree->GetAlias("pf_hferechit_e") != 0) {
		pf_hferechit_e_branch = tree->GetBranch(tree->GetAlias("pf_hferechit_e"));
		pf_hferechit_e_branch->SetAddress(&pf_hferechit_e_);
	}
	if(pf_hferechit_e_branch == 0 ) {
	cout << "Branch pf_hferechit_e does not exist." << endl;
	}
	pf_hferechit_et_branch = 0;
	if (tree->GetAlias("pf_hferechit_et") != 0) {
		pf_hferechit_et_branch = tree->GetBranch(tree->GetAlias("pf_hferechit_et"));
		pf_hferechit_et_branch->SetAddress(&pf_hferechit_et_);
	}
	if(pf_hferechit_et_branch == 0 ) {
	cout << "Branch pf_hferechit_et does not exist." << endl;
	}
	pf_hferechit_eta_branch = 0;
	if (tree->GetAlias("pf_hferechit_eta") != 0) {
		pf_hferechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_hferechit_eta"));
		pf_hferechit_eta_branch->SetAddress(&pf_hferechit_eta_);
	}
	if(pf_hferechit_eta_branch == 0 ) {
	cout << "Branch pf_hferechit_eta does not exist." << endl;
	}
	pf_hferechit_phi_branch = 0;
	if (tree->GetAlias("pf_hferechit_phi") != 0) {
		pf_hferechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_hferechit_phi"));
		pf_hferechit_phi_branch->SetAddress(&pf_hferechit_phi_);
	}
	if(pf_hferechit_phi_branch == 0 ) {
	cout << "Branch pf_hferechit_phi does not exist." << endl;
	}
	pf_hfhcluster_e_branch = 0;
	if (tree->GetAlias("pf_hfhcluster_e") != 0) {
		pf_hfhcluster_e_branch = tree->GetBranch(tree->GetAlias("pf_hfhcluster_e"));
		pf_hfhcluster_e_branch->SetAddress(&pf_hfhcluster_e_);
	}
	if(pf_hfhcluster_e_branch == 0 ) {
	cout << "Branch pf_hfhcluster_e does not exist." << endl;
	}
	pf_hfhcluster_et_branch = 0;
	if (tree->GetAlias("pf_hfhcluster_et") != 0) {
		pf_hfhcluster_et_branch = tree->GetBranch(tree->GetAlias("pf_hfhcluster_et"));
		pf_hfhcluster_et_branch->SetAddress(&pf_hfhcluster_et_);
	}
	if(pf_hfhcluster_et_branch == 0 ) {
	cout << "Branch pf_hfhcluster_et does not exist." << endl;
	}
	pf_hfhcluster_eta_branch = 0;
	if (tree->GetAlias("pf_hfhcluster_eta") != 0) {
		pf_hfhcluster_eta_branch = tree->GetBranch(tree->GetAlias("pf_hfhcluster_eta"));
		pf_hfhcluster_eta_branch->SetAddress(&pf_hfhcluster_eta_);
	}
	if(pf_hfhcluster_eta_branch == 0 ) {
	cout << "Branch pf_hfhcluster_eta does not exist." << endl;
	}
	pf_hfhcluster_phi_branch = 0;
	if (tree->GetAlias("pf_hfhcluster_phi") != 0) {
		pf_hfhcluster_phi_branch = tree->GetBranch(tree->GetAlias("pf_hfhcluster_phi"));
		pf_hfhcluster_phi_branch->SetAddress(&pf_hfhcluster_phi_);
	}
	if(pf_hfhcluster_phi_branch == 0 ) {
	cout << "Branch pf_hfhcluster_phi does not exist." << endl;
	}
	pf_hfhrechit_e_branch = 0;
	if (tree->GetAlias("pf_hfhrechit_e") != 0) {
		pf_hfhrechit_e_branch = tree->GetBranch(tree->GetAlias("pf_hfhrechit_e"));
		pf_hfhrechit_e_branch->SetAddress(&pf_hfhrechit_e_);
	}
	if(pf_hfhrechit_e_branch == 0 ) {
	cout << "Branch pf_hfhrechit_e does not exist." << endl;
	}
	pf_hfhrechit_et_branch = 0;
	if (tree->GetAlias("pf_hfhrechit_et") != 0) {
		pf_hfhrechit_et_branch = tree->GetBranch(tree->GetAlias("pf_hfhrechit_et"));
		pf_hfhrechit_et_branch->SetAddress(&pf_hfhrechit_et_);
	}
	if(pf_hfhrechit_et_branch == 0 ) {
	cout << "Branch pf_hfhrechit_et does not exist." << endl;
	}
	pf_hfhrechit_eta_branch = 0;
	if (tree->GetAlias("pf_hfhrechit_eta") != 0) {
		pf_hfhrechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_hfhrechit_eta"));
		pf_hfhrechit_eta_branch->SetAddress(&pf_hfhrechit_eta_);
	}
	if(pf_hfhrechit_eta_branch == 0 ) {
	cout << "Branch pf_hfhrechit_eta does not exist." << endl;
	}
	pf_hfhrechit_phi_branch = 0;
	if (tree->GetAlias("pf_hfhrechit_phi") != 0) {
		pf_hfhrechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_hfhrechit_phi"));
		pf_hfhrechit_phi_branch->SetAddress(&pf_hfhrechit_phi_);
	}
	if(pf_hfhrechit_phi_branch == 0 ) {
	cout << "Branch pf_hfhrechit_phi does not exist." << endl;
	}
	pf_rechit_detid_branch = 0;
	if (tree->GetAlias("pf_rechit_detid") != 0) {
		pf_rechit_detid_branch = tree->GetBranch(tree->GetAlias("pf_rechit_detid"));
		pf_rechit_detid_branch->SetAddress(&pf_rechit_detid_);
	}
	if(pf_rechit_detid_branch == 0 ) {
	cout << "Branch pf_rechit_detid does not exist." << endl;
	}
	pf_rechit_e_branch = 0;
	if (tree->GetAlias("pf_rechit_e") != 0) {
		pf_rechit_e_branch = tree->GetBranch(tree->GetAlias("pf_rechit_e"));
		pf_rechit_e_branch->SetAddress(&pf_rechit_e_);
	}
	if(pf_rechit_e_branch == 0 ) {
	cout << "Branch pf_rechit_e does not exist." << endl;
	}
	pf_rechit_et_branch = 0;
	if (tree->GetAlias("pf_rechit_et") != 0) {
		pf_rechit_et_branch = tree->GetBranch(tree->GetAlias("pf_rechit_et"));
		pf_rechit_et_branch->SetAddress(&pf_rechit_et_);
	}
	if(pf_rechit_et_branch == 0 ) {
	cout << "Branch pf_rechit_et does not exist." << endl;
	}
	pf_rechit_eta_branch = 0;
	if (tree->GetAlias("pf_rechit_eta") != 0) {
		pf_rechit_eta_branch = tree->GetBranch(tree->GetAlias("pf_rechit_eta"));
		pf_rechit_eta_branch->SetAddress(&pf_rechit_eta_);
	}
	if(pf_rechit_eta_branch == 0 ) {
	cout << "Branch pf_rechit_eta does not exist." << endl;
	}
	pf_rechit_phi_branch = 0;
	if (tree->GetAlias("pf_rechit_phi") != 0) {
		pf_rechit_phi_branch = tree->GetBranch(tree->GetAlias("pf_rechit_phi"));
		pf_rechit_phi_branch->SetAddress(&pf_rechit_phi_);
	}
	if(pf_rechit_phi_branch == 0 ) {
	cout << "Branch pf_rechit_phi does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		calo_eb_sumet_isLoaded = false;
		calo_ee_sumet_isLoaded = false;
		calo_hb_sumet_isLoaded = false;
		calo_he_sumet_isLoaded = false;
		calo_hfe_sumet_isLoaded = false;
		calo_hfh_sumet_isLoaded = false;
		calomet_isLoaded = false;
		calosumet_isLoaded = false;
		genmet_isLoaded = false;
		gensumet_isLoaded = false;
		pfclus_eb_met_isLoaded = false;
		pfclus_eb_sumet_isLoaded = false;
		pfclus_ee_met_isLoaded = false;
		pfclus_ee_sumet_isLoaded = false;
		pfclus_hb_met_isLoaded = false;
		pfclus_hb_sumet_isLoaded = false;
		pfclus_he_met_isLoaded = false;
		pfclus_he_sumet_isLoaded = false;
		pfclus_hfe_met_isLoaded = false;
		pfclus_hfe_sumet_isLoaded = false;
		pfclus_hfh_met_isLoaded = false;
		pfclus_hfh_sumet_isLoaded = false;
		pfclusmet_isLoaded = false;
		pfclussumet_isLoaded = false;
		pf_eb_met_isLoaded = false;
		pf_eb_sumet_isLoaded = false;
		pf_ee_met_isLoaded = false;
		pf_ee_sumet_isLoaded = false;
		pf_hb_met_isLoaded = false;
		pf_hb_sumet_isLoaded = false;
		pf_he_met_isLoaded = false;
		pf_he_sumet_isLoaded = false;
		pf_hfe_met_isLoaded = false;
		pf_hfe_sumet_isLoaded = false;
		pf_hfh_met_isLoaded = false;
		pf_hfh_sumet_isLoaded = false;
		pfmet_isLoaded = false;
		pfsumet_isLoaded = false;
		pf_cluster_detid_isLoaded = false;
		pf_cluster_e_isLoaded = false;
		pf_cluster_et_isLoaded = false;
		pf_cluster_eta_isLoaded = false;
		pf_cluster_phi_isLoaded = false;
		pf_ebcluster_e_isLoaded = false;
		pf_ebcluster_et_isLoaded = false;
		pf_ebcluster_eta_isLoaded = false;
		pf_ebcluster_phi_isLoaded = false;
		pf_ebrechit_e_isLoaded = false;
		pf_ebrechit_et_isLoaded = false;
		pf_ebrechit_eta_isLoaded = false;
		pf_ebrechit_phi_isLoaded = false;
		pf_eecluster_e_isLoaded = false;
		pf_eecluster_et_isLoaded = false;
		pf_eecluster_eta_isLoaded = false;
		pf_eecluster_phi_isLoaded = false;
		pf_eerechit_e_isLoaded = false;
		pf_eerechit_et_isLoaded = false;
		pf_eerechit_eta_isLoaded = false;
		pf_eerechit_phi_isLoaded = false;
		pf_hbcluster_e_isLoaded = false;
		pf_hbcluster_et_isLoaded = false;
		pf_hbcluster_eta_isLoaded = false;
		pf_hbcluster_phi_isLoaded = false;
		pf_hbrechit_e_isLoaded = false;
		pf_hbrechit_et_isLoaded = false;
		pf_hbrechit_eta_isLoaded = false;
		pf_hbrechit_phi_isLoaded = false;
		pf_hecluster_e_isLoaded = false;
		pf_hecluster_et_isLoaded = false;
		pf_hecluster_eta_isLoaded = false;
		pf_hecluster_phi_isLoaded = false;
		pf_herechit_e_isLoaded = false;
		pf_herechit_et_isLoaded = false;
		pf_herechit_eta_isLoaded = false;
		pf_herechit_phi_isLoaded = false;
		pf_hfecluster_e_isLoaded = false;
		pf_hfecluster_et_isLoaded = false;
		pf_hfecluster_eta_isLoaded = false;
		pf_hfecluster_phi_isLoaded = false;
		pf_hferechit_e_isLoaded = false;
		pf_hferechit_et_isLoaded = false;
		pf_hferechit_eta_isLoaded = false;
		pf_hferechit_phi_isLoaded = false;
		pf_hfhcluster_e_isLoaded = false;
		pf_hfhcluster_et_isLoaded = false;
		pf_hfhcluster_eta_isLoaded = false;
		pf_hfhcluster_phi_isLoaded = false;
		pf_hfhrechit_e_isLoaded = false;
		pf_hfhrechit_et_isLoaded = false;
		pf_hfhrechit_eta_isLoaded = false;
		pf_hfhrechit_phi_isLoaded = false;
		pf_rechit_detid_isLoaded = false;
		pf_rechit_e_isLoaded = false;
		pf_rechit_et_isLoaded = false;
		pf_rechit_eta_isLoaded = false;
		pf_rechit_phi_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (calo_eb_sumet_branch != 0) calo_eb_sumet();
	if (calo_ee_sumet_branch != 0) calo_ee_sumet();
	if (calo_hb_sumet_branch != 0) calo_hb_sumet();
	if (calo_he_sumet_branch != 0) calo_he_sumet();
	if (calo_hfe_sumet_branch != 0) calo_hfe_sumet();
	if (calo_hfh_sumet_branch != 0) calo_hfh_sumet();
	if (calomet_branch != 0) calomet();
	if (calosumet_branch != 0) calosumet();
	if (genmet_branch != 0) genmet();
	if (gensumet_branch != 0) gensumet();
	if (pfclus_eb_met_branch != 0) pfclus_eb_met();
	if (pfclus_eb_sumet_branch != 0) pfclus_eb_sumet();
	if (pfclus_ee_met_branch != 0) pfclus_ee_met();
	if (pfclus_ee_sumet_branch != 0) pfclus_ee_sumet();
	if (pfclus_hb_met_branch != 0) pfclus_hb_met();
	if (pfclus_hb_sumet_branch != 0) pfclus_hb_sumet();
	if (pfclus_he_met_branch != 0) pfclus_he_met();
	if (pfclus_he_sumet_branch != 0) pfclus_he_sumet();
	if (pfclus_hfe_met_branch != 0) pfclus_hfe_met();
	if (pfclus_hfe_sumet_branch != 0) pfclus_hfe_sumet();
	if (pfclus_hfh_met_branch != 0) pfclus_hfh_met();
	if (pfclus_hfh_sumet_branch != 0) pfclus_hfh_sumet();
	if (pfclusmet_branch != 0) pfclusmet();
	if (pfclussumet_branch != 0) pfclussumet();
	if (pf_eb_met_branch != 0) pf_eb_met();
	if (pf_eb_sumet_branch != 0) pf_eb_sumet();
	if (pf_ee_met_branch != 0) pf_ee_met();
	if (pf_ee_sumet_branch != 0) pf_ee_sumet();
	if (pf_hb_met_branch != 0) pf_hb_met();
	if (pf_hb_sumet_branch != 0) pf_hb_sumet();
	if (pf_he_met_branch != 0) pf_he_met();
	if (pf_he_sumet_branch != 0) pf_he_sumet();
	if (pf_hfe_met_branch != 0) pf_hfe_met();
	if (pf_hfe_sumet_branch != 0) pf_hfe_sumet();
	if (pf_hfh_met_branch != 0) pf_hfh_met();
	if (pf_hfh_sumet_branch != 0) pf_hfh_sumet();
	if (pfmet_branch != 0) pfmet();
	if (pfsumet_branch != 0) pfsumet();
	if (pf_cluster_detid_branch != 0) pf_cluster_detid();
	if (pf_cluster_e_branch != 0) pf_cluster_e();
	if (pf_cluster_et_branch != 0) pf_cluster_et();
	if (pf_cluster_eta_branch != 0) pf_cluster_eta();
	if (pf_cluster_phi_branch != 0) pf_cluster_phi();
	if (pf_ebcluster_e_branch != 0) pf_ebcluster_e();
	if (pf_ebcluster_et_branch != 0) pf_ebcluster_et();
	if (pf_ebcluster_eta_branch != 0) pf_ebcluster_eta();
	if (pf_ebcluster_phi_branch != 0) pf_ebcluster_phi();
	if (pf_ebrechit_e_branch != 0) pf_ebrechit_e();
	if (pf_ebrechit_et_branch != 0) pf_ebrechit_et();
	if (pf_ebrechit_eta_branch != 0) pf_ebrechit_eta();
	if (pf_ebrechit_phi_branch != 0) pf_ebrechit_phi();
	if (pf_eecluster_e_branch != 0) pf_eecluster_e();
	if (pf_eecluster_et_branch != 0) pf_eecluster_et();
	if (pf_eecluster_eta_branch != 0) pf_eecluster_eta();
	if (pf_eecluster_phi_branch != 0) pf_eecluster_phi();
	if (pf_eerechit_e_branch != 0) pf_eerechit_e();
	if (pf_eerechit_et_branch != 0) pf_eerechit_et();
	if (pf_eerechit_eta_branch != 0) pf_eerechit_eta();
	if (pf_eerechit_phi_branch != 0) pf_eerechit_phi();
	if (pf_hbcluster_e_branch != 0) pf_hbcluster_e();
	if (pf_hbcluster_et_branch != 0) pf_hbcluster_et();
	if (pf_hbcluster_eta_branch != 0) pf_hbcluster_eta();
	if (pf_hbcluster_phi_branch != 0) pf_hbcluster_phi();
	if (pf_hbrechit_e_branch != 0) pf_hbrechit_e();
	if (pf_hbrechit_et_branch != 0) pf_hbrechit_et();
	if (pf_hbrechit_eta_branch != 0) pf_hbrechit_eta();
	if (pf_hbrechit_phi_branch != 0) pf_hbrechit_phi();
	if (pf_hecluster_e_branch != 0) pf_hecluster_e();
	if (pf_hecluster_et_branch != 0) pf_hecluster_et();
	if (pf_hecluster_eta_branch != 0) pf_hecluster_eta();
	if (pf_hecluster_phi_branch != 0) pf_hecluster_phi();
	if (pf_herechit_e_branch != 0) pf_herechit_e();
	if (pf_herechit_et_branch != 0) pf_herechit_et();
	if (pf_herechit_eta_branch != 0) pf_herechit_eta();
	if (pf_herechit_phi_branch != 0) pf_herechit_phi();
	if (pf_hfecluster_e_branch != 0) pf_hfecluster_e();
	if (pf_hfecluster_et_branch != 0) pf_hfecluster_et();
	if (pf_hfecluster_eta_branch != 0) pf_hfecluster_eta();
	if (pf_hfecluster_phi_branch != 0) pf_hfecluster_phi();
	if (pf_hferechit_e_branch != 0) pf_hferechit_e();
	if (pf_hferechit_et_branch != 0) pf_hferechit_et();
	if (pf_hferechit_eta_branch != 0) pf_hferechit_eta();
	if (pf_hferechit_phi_branch != 0) pf_hferechit_phi();
	if (pf_hfhcluster_e_branch != 0) pf_hfhcluster_e();
	if (pf_hfhcluster_et_branch != 0) pf_hfhcluster_et();
	if (pf_hfhcluster_eta_branch != 0) pf_hfhcluster_eta();
	if (pf_hfhcluster_phi_branch != 0) pf_hfhcluster_phi();
	if (pf_hfhrechit_e_branch != 0) pf_hfhrechit_e();
	if (pf_hfhrechit_et_branch != 0) pf_hfhrechit_et();
	if (pf_hfhrechit_eta_branch != 0) pf_hfhrechit_eta();
	if (pf_hfhrechit_phi_branch != 0) pf_hfhrechit_phi();
	if (pf_rechit_detid_branch != 0) pf_rechit_detid();
	if (pf_rechit_e_branch != 0) pf_rechit_e();
	if (pf_rechit_et_branch != 0) pf_rechit_et();
	if (pf_rechit_eta_branch != 0) pf_rechit_eta();
	if (pf_rechit_phi_branch != 0) pf_rechit_phi();
}

	float &calo_eb_sumet()
	{
		if (not calo_eb_sumet_isLoaded) {
			if (calo_eb_sumet_branch != 0) {
				calo_eb_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_eb_sumet_)) {
					printf("branch calo_eb_sumet_branch contains a bad float: %f\n", calo_eb_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_eb_sumet_branch does not exist!\n");
				exit(1);
			}
			calo_eb_sumet_isLoaded = true;
		}
		return calo_eb_sumet_;
	}
	float &calo_ee_sumet()
	{
		if (not calo_ee_sumet_isLoaded) {
			if (calo_ee_sumet_branch != 0) {
				calo_ee_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_ee_sumet_)) {
					printf("branch calo_ee_sumet_branch contains a bad float: %f\n", calo_ee_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_ee_sumet_branch does not exist!\n");
				exit(1);
			}
			calo_ee_sumet_isLoaded = true;
		}
		return calo_ee_sumet_;
	}
	float &calo_hb_sumet()
	{
		if (not calo_hb_sumet_isLoaded) {
			if (calo_hb_sumet_branch != 0) {
				calo_hb_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_hb_sumet_)) {
					printf("branch calo_hb_sumet_branch contains a bad float: %f\n", calo_hb_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_hb_sumet_branch does not exist!\n");
				exit(1);
			}
			calo_hb_sumet_isLoaded = true;
		}
		return calo_hb_sumet_;
	}
	float &calo_he_sumet()
	{
		if (not calo_he_sumet_isLoaded) {
			if (calo_he_sumet_branch != 0) {
				calo_he_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_he_sumet_)) {
					printf("branch calo_he_sumet_branch contains a bad float: %f\n", calo_he_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_he_sumet_branch does not exist!\n");
				exit(1);
			}
			calo_he_sumet_isLoaded = true;
		}
		return calo_he_sumet_;
	}
	float &calo_hfe_sumet()
	{
		if (not calo_hfe_sumet_isLoaded) {
			if (calo_hfe_sumet_branch != 0) {
				calo_hfe_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_hfe_sumet_)) {
					printf("branch calo_hfe_sumet_branch contains a bad float: %f\n", calo_hfe_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_hfe_sumet_branch does not exist!\n");
				exit(1);
			}
			calo_hfe_sumet_isLoaded = true;
		}
		return calo_hfe_sumet_;
	}
	float &calo_hfh_sumet()
	{
		if (not calo_hfh_sumet_isLoaded) {
			if (calo_hfh_sumet_branch != 0) {
				calo_hfh_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calo_hfh_sumet_)) {
					printf("branch calo_hfh_sumet_branch contains a bad float: %f\n", calo_hfh_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calo_hfh_sumet_branch does not exist!\n");
				exit(1);
			}
			calo_hfh_sumet_isLoaded = true;
		}
		return calo_hfh_sumet_;
	}
	float &calomet()
	{
		if (not calomet_isLoaded) {
			if (calomet_branch != 0) {
				calomet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calomet_)) {
					printf("branch calomet_branch contains a bad float: %f\n", calomet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calomet_branch does not exist!\n");
				exit(1);
			}
			calomet_isLoaded = true;
		}
		return calomet_;
	}
	float &calosumet()
	{
		if (not calosumet_isLoaded) {
			if (calosumet_branch != 0) {
				calosumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(calosumet_)) {
					printf("branch calosumet_branch contains a bad float: %f\n", calosumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch calosumet_branch does not exist!\n");
				exit(1);
			}
			calosumet_isLoaded = true;
		}
		return calosumet_;
	}
	float &genmet()
	{
		if (not genmet_isLoaded) {
			if (genmet_branch != 0) {
				genmet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(genmet_)) {
					printf("branch genmet_branch contains a bad float: %f\n", genmet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genmet_branch does not exist!\n");
				exit(1);
			}
			genmet_isLoaded = true;
		}
		return genmet_;
	}
	float &gensumet()
	{
		if (not gensumet_isLoaded) {
			if (gensumet_branch != 0) {
				gensumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(gensumet_)) {
					printf("branch gensumet_branch contains a bad float: %f\n", gensumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gensumet_branch does not exist!\n");
				exit(1);
			}
			gensumet_isLoaded = true;
		}
		return gensumet_;
	}
	float &pfclus_eb_met()
	{
		if (not pfclus_eb_met_isLoaded) {
			if (pfclus_eb_met_branch != 0) {
				pfclus_eb_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_eb_met_)) {
					printf("branch pfclus_eb_met_branch contains a bad float: %f\n", pfclus_eb_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_eb_met_branch does not exist!\n");
				exit(1);
			}
			pfclus_eb_met_isLoaded = true;
		}
		return pfclus_eb_met_;
	}
	float &pfclus_eb_sumet()
	{
		if (not pfclus_eb_sumet_isLoaded) {
			if (pfclus_eb_sumet_branch != 0) {
				pfclus_eb_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_eb_sumet_)) {
					printf("branch pfclus_eb_sumet_branch contains a bad float: %f\n", pfclus_eb_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_eb_sumet_branch does not exist!\n");
				exit(1);
			}
			pfclus_eb_sumet_isLoaded = true;
		}
		return pfclus_eb_sumet_;
	}
	float &pfclus_ee_met()
	{
		if (not pfclus_ee_met_isLoaded) {
			if (pfclus_ee_met_branch != 0) {
				pfclus_ee_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_ee_met_)) {
					printf("branch pfclus_ee_met_branch contains a bad float: %f\n", pfclus_ee_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_ee_met_branch does not exist!\n");
				exit(1);
			}
			pfclus_ee_met_isLoaded = true;
		}
		return pfclus_ee_met_;
	}
	float &pfclus_ee_sumet()
	{
		if (not pfclus_ee_sumet_isLoaded) {
			if (pfclus_ee_sumet_branch != 0) {
				pfclus_ee_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_ee_sumet_)) {
					printf("branch pfclus_ee_sumet_branch contains a bad float: %f\n", pfclus_ee_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_ee_sumet_branch does not exist!\n");
				exit(1);
			}
			pfclus_ee_sumet_isLoaded = true;
		}
		return pfclus_ee_sumet_;
	}
	float &pfclus_hb_met()
	{
		if (not pfclus_hb_met_isLoaded) {
			if (pfclus_hb_met_branch != 0) {
				pfclus_hb_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_hb_met_)) {
					printf("branch pfclus_hb_met_branch contains a bad float: %f\n", pfclus_hb_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_hb_met_branch does not exist!\n");
				exit(1);
			}
			pfclus_hb_met_isLoaded = true;
		}
		return pfclus_hb_met_;
	}
	float &pfclus_hb_sumet()
	{
		if (not pfclus_hb_sumet_isLoaded) {
			if (pfclus_hb_sumet_branch != 0) {
				pfclus_hb_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_hb_sumet_)) {
					printf("branch pfclus_hb_sumet_branch contains a bad float: %f\n", pfclus_hb_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_hb_sumet_branch does not exist!\n");
				exit(1);
			}
			pfclus_hb_sumet_isLoaded = true;
		}
		return pfclus_hb_sumet_;
	}
	float &pfclus_he_met()
	{
		if (not pfclus_he_met_isLoaded) {
			if (pfclus_he_met_branch != 0) {
				pfclus_he_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_he_met_)) {
					printf("branch pfclus_he_met_branch contains a bad float: %f\n", pfclus_he_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_he_met_branch does not exist!\n");
				exit(1);
			}
			pfclus_he_met_isLoaded = true;
		}
		return pfclus_he_met_;
	}
	float &pfclus_he_sumet()
	{
		if (not pfclus_he_sumet_isLoaded) {
			if (pfclus_he_sumet_branch != 0) {
				pfclus_he_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_he_sumet_)) {
					printf("branch pfclus_he_sumet_branch contains a bad float: %f\n", pfclus_he_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_he_sumet_branch does not exist!\n");
				exit(1);
			}
			pfclus_he_sumet_isLoaded = true;
		}
		return pfclus_he_sumet_;
	}
	float &pfclus_hfe_met()
	{
		if (not pfclus_hfe_met_isLoaded) {
			if (pfclus_hfe_met_branch != 0) {
				pfclus_hfe_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_hfe_met_)) {
					printf("branch pfclus_hfe_met_branch contains a bad float: %f\n", pfclus_hfe_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_hfe_met_branch does not exist!\n");
				exit(1);
			}
			pfclus_hfe_met_isLoaded = true;
		}
		return pfclus_hfe_met_;
	}
	float &pfclus_hfe_sumet()
	{
		if (not pfclus_hfe_sumet_isLoaded) {
			if (pfclus_hfe_sumet_branch != 0) {
				pfclus_hfe_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_hfe_sumet_)) {
					printf("branch pfclus_hfe_sumet_branch contains a bad float: %f\n", pfclus_hfe_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_hfe_sumet_branch does not exist!\n");
				exit(1);
			}
			pfclus_hfe_sumet_isLoaded = true;
		}
		return pfclus_hfe_sumet_;
	}
	float &pfclus_hfh_met()
	{
		if (not pfclus_hfh_met_isLoaded) {
			if (pfclus_hfh_met_branch != 0) {
				pfclus_hfh_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_hfh_met_)) {
					printf("branch pfclus_hfh_met_branch contains a bad float: %f\n", pfclus_hfh_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_hfh_met_branch does not exist!\n");
				exit(1);
			}
			pfclus_hfh_met_isLoaded = true;
		}
		return pfclus_hfh_met_;
	}
	float &pfclus_hfh_sumet()
	{
		if (not pfclus_hfh_sumet_isLoaded) {
			if (pfclus_hfh_sumet_branch != 0) {
				pfclus_hfh_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclus_hfh_sumet_)) {
					printf("branch pfclus_hfh_sumet_branch contains a bad float: %f\n", pfclus_hfh_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclus_hfh_sumet_branch does not exist!\n");
				exit(1);
			}
			pfclus_hfh_sumet_isLoaded = true;
		}
		return pfclus_hfh_sumet_;
	}
	float &pfclusmet()
	{
		if (not pfclusmet_isLoaded) {
			if (pfclusmet_branch != 0) {
				pfclusmet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclusmet_)) {
					printf("branch pfclusmet_branch contains a bad float: %f\n", pfclusmet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclusmet_branch does not exist!\n");
				exit(1);
			}
			pfclusmet_isLoaded = true;
		}
		return pfclusmet_;
	}
	float &pfclussumet()
	{
		if (not pfclussumet_isLoaded) {
			if (pfclussumet_branch != 0) {
				pfclussumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfclussumet_)) {
					printf("branch pfclussumet_branch contains a bad float: %f\n", pfclussumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfclussumet_branch does not exist!\n");
				exit(1);
			}
			pfclussumet_isLoaded = true;
		}
		return pfclussumet_;
	}
	float &pf_eb_met()
	{
		if (not pf_eb_met_isLoaded) {
			if (pf_eb_met_branch != 0) {
				pf_eb_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_eb_met_)) {
					printf("branch pf_eb_met_branch contains a bad float: %f\n", pf_eb_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eb_met_branch does not exist!\n");
				exit(1);
			}
			pf_eb_met_isLoaded = true;
		}
		return pf_eb_met_;
	}
	float &pf_eb_sumet()
	{
		if (not pf_eb_sumet_isLoaded) {
			if (pf_eb_sumet_branch != 0) {
				pf_eb_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_eb_sumet_)) {
					printf("branch pf_eb_sumet_branch contains a bad float: %f\n", pf_eb_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eb_sumet_branch does not exist!\n");
				exit(1);
			}
			pf_eb_sumet_isLoaded = true;
		}
		return pf_eb_sumet_;
	}
	float &pf_ee_met()
	{
		if (not pf_ee_met_isLoaded) {
			if (pf_ee_met_branch != 0) {
				pf_ee_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_ee_met_)) {
					printf("branch pf_ee_met_branch contains a bad float: %f\n", pf_ee_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ee_met_branch does not exist!\n");
				exit(1);
			}
			pf_ee_met_isLoaded = true;
		}
		return pf_ee_met_;
	}
	float &pf_ee_sumet()
	{
		if (not pf_ee_sumet_isLoaded) {
			if (pf_ee_sumet_branch != 0) {
				pf_ee_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_ee_sumet_)) {
					printf("branch pf_ee_sumet_branch contains a bad float: %f\n", pf_ee_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ee_sumet_branch does not exist!\n");
				exit(1);
			}
			pf_ee_sumet_isLoaded = true;
		}
		return pf_ee_sumet_;
	}
	float &pf_hb_met()
	{
		if (not pf_hb_met_isLoaded) {
			if (pf_hb_met_branch != 0) {
				pf_hb_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_hb_met_)) {
					printf("branch pf_hb_met_branch contains a bad float: %f\n", pf_hb_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hb_met_branch does not exist!\n");
				exit(1);
			}
			pf_hb_met_isLoaded = true;
		}
		return pf_hb_met_;
	}
	float &pf_hb_sumet()
	{
		if (not pf_hb_sumet_isLoaded) {
			if (pf_hb_sumet_branch != 0) {
				pf_hb_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_hb_sumet_)) {
					printf("branch pf_hb_sumet_branch contains a bad float: %f\n", pf_hb_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hb_sumet_branch does not exist!\n");
				exit(1);
			}
			pf_hb_sumet_isLoaded = true;
		}
		return pf_hb_sumet_;
	}
	float &pf_he_met()
	{
		if (not pf_he_met_isLoaded) {
			if (pf_he_met_branch != 0) {
				pf_he_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_he_met_)) {
					printf("branch pf_he_met_branch contains a bad float: %f\n", pf_he_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_he_met_branch does not exist!\n");
				exit(1);
			}
			pf_he_met_isLoaded = true;
		}
		return pf_he_met_;
	}
	float &pf_he_sumet()
	{
		if (not pf_he_sumet_isLoaded) {
			if (pf_he_sumet_branch != 0) {
				pf_he_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_he_sumet_)) {
					printf("branch pf_he_sumet_branch contains a bad float: %f\n", pf_he_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_he_sumet_branch does not exist!\n");
				exit(1);
			}
			pf_he_sumet_isLoaded = true;
		}
		return pf_he_sumet_;
	}
	float &pf_hfe_met()
	{
		if (not pf_hfe_met_isLoaded) {
			if (pf_hfe_met_branch != 0) {
				pf_hfe_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_hfe_met_)) {
					printf("branch pf_hfe_met_branch contains a bad float: %f\n", pf_hfe_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfe_met_branch does not exist!\n");
				exit(1);
			}
			pf_hfe_met_isLoaded = true;
		}
		return pf_hfe_met_;
	}
	float &pf_hfe_sumet()
	{
		if (not pf_hfe_sumet_isLoaded) {
			if (pf_hfe_sumet_branch != 0) {
				pf_hfe_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_hfe_sumet_)) {
					printf("branch pf_hfe_sumet_branch contains a bad float: %f\n", pf_hfe_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfe_sumet_branch does not exist!\n");
				exit(1);
			}
			pf_hfe_sumet_isLoaded = true;
		}
		return pf_hfe_sumet_;
	}
	float &pf_hfh_met()
	{
		if (not pf_hfh_met_isLoaded) {
			if (pf_hfh_met_branch != 0) {
				pf_hfh_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_hfh_met_)) {
					printf("branch pf_hfh_met_branch contains a bad float: %f\n", pf_hfh_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfh_met_branch does not exist!\n");
				exit(1);
			}
			pf_hfh_met_isLoaded = true;
		}
		return pf_hfh_met_;
	}
	float &pf_hfh_sumet()
	{
		if (not pf_hfh_sumet_isLoaded) {
			if (pf_hfh_sumet_branch != 0) {
				pf_hfh_sumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pf_hfh_sumet_)) {
					printf("branch pf_hfh_sumet_branch contains a bad float: %f\n", pf_hfh_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfh_sumet_branch does not exist!\n");
				exit(1);
			}
			pf_hfh_sumet_isLoaded = true;
		}
		return pf_hfh_sumet_;
	}
	float &pfmet()
	{
		if (not pfmet_isLoaded) {
			if (pfmet_branch != 0) {
				pfmet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfmet_)) {
					printf("branch pfmet_branch contains a bad float: %f\n", pfmet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfmet_branch does not exist!\n");
				exit(1);
			}
			pfmet_isLoaded = true;
		}
		return pfmet_;
	}
	float &pfsumet()
	{
		if (not pfsumet_isLoaded) {
			if (pfsumet_branch != 0) {
				pfsumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(pfsumet_)) {
					printf("branch pfsumet_branch contains a bad float: %f\n", pfsumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pfsumet_branch does not exist!\n");
				exit(1);
			}
			pfsumet_isLoaded = true;
		}
		return pfsumet_;
	}
	vector<float> &pf_cluster_detid()
	{
		if (not pf_cluster_detid_isLoaded) {
			if (pf_cluster_detid_branch != 0) {
				pf_cluster_detid_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_cluster_detid_.begin(); i != pf_cluster_detid_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_cluster_detid_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_cluster_detid_branch does not exist!\n");
				exit(1);
			}
			pf_cluster_detid_isLoaded = true;
		}
		return pf_cluster_detid_;
	}
	vector<float> &pf_cluster_e()
	{
		if (not pf_cluster_e_isLoaded) {
			if (pf_cluster_e_branch != 0) {
				pf_cluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_cluster_e_.begin(); i != pf_cluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_cluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_cluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_cluster_e_isLoaded = true;
		}
		return pf_cluster_e_;
	}
	vector<float> &pf_cluster_et()
	{
		if (not pf_cluster_et_isLoaded) {
			if (pf_cluster_et_branch != 0) {
				pf_cluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_cluster_et_.begin(); i != pf_cluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_cluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_cluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_cluster_et_isLoaded = true;
		}
		return pf_cluster_et_;
	}
	vector<float> &pf_cluster_eta()
	{
		if (not pf_cluster_eta_isLoaded) {
			if (pf_cluster_eta_branch != 0) {
				pf_cluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_cluster_eta_.begin(); i != pf_cluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_cluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_cluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_cluster_eta_isLoaded = true;
		}
		return pf_cluster_eta_;
	}
	vector<float> &pf_cluster_phi()
	{
		if (not pf_cluster_phi_isLoaded) {
			if (pf_cluster_phi_branch != 0) {
				pf_cluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_cluster_phi_.begin(); i != pf_cluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_cluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_cluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_cluster_phi_isLoaded = true;
		}
		return pf_cluster_phi_;
	}
	vector<float> &pf_ebcluster_e()
	{
		if (not pf_ebcluster_e_isLoaded) {
			if (pf_ebcluster_e_branch != 0) {
				pf_ebcluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebcluster_e_.begin(); i != pf_ebcluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebcluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebcluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_ebcluster_e_isLoaded = true;
		}
		return pf_ebcluster_e_;
	}
	vector<float> &pf_ebcluster_et()
	{
		if (not pf_ebcluster_et_isLoaded) {
			if (pf_ebcluster_et_branch != 0) {
				pf_ebcluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebcluster_et_.begin(); i != pf_ebcluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebcluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebcluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_ebcluster_et_isLoaded = true;
		}
		return pf_ebcluster_et_;
	}
	vector<float> &pf_ebcluster_eta()
	{
		if (not pf_ebcluster_eta_isLoaded) {
			if (pf_ebcluster_eta_branch != 0) {
				pf_ebcluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebcluster_eta_.begin(); i != pf_ebcluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebcluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebcluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_ebcluster_eta_isLoaded = true;
		}
		return pf_ebcluster_eta_;
	}
	vector<float> &pf_ebcluster_phi()
	{
		if (not pf_ebcluster_phi_isLoaded) {
			if (pf_ebcluster_phi_branch != 0) {
				pf_ebcluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebcluster_phi_.begin(); i != pf_ebcluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebcluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebcluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_ebcluster_phi_isLoaded = true;
		}
		return pf_ebcluster_phi_;
	}
	vector<float> &pf_ebrechit_e()
	{
		if (not pf_ebrechit_e_isLoaded) {
			if (pf_ebrechit_e_branch != 0) {
				pf_ebrechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebrechit_e_.begin(); i != pf_ebrechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebrechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebrechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_ebrechit_e_isLoaded = true;
		}
		return pf_ebrechit_e_;
	}
	vector<float> &pf_ebrechit_et()
	{
		if (not pf_ebrechit_et_isLoaded) {
			if (pf_ebrechit_et_branch != 0) {
				pf_ebrechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebrechit_et_.begin(); i != pf_ebrechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebrechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebrechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_ebrechit_et_isLoaded = true;
		}
		return pf_ebrechit_et_;
	}
	vector<float> &pf_ebrechit_eta()
	{
		if (not pf_ebrechit_eta_isLoaded) {
			if (pf_ebrechit_eta_branch != 0) {
				pf_ebrechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebrechit_eta_.begin(); i != pf_ebrechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebrechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebrechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_ebrechit_eta_isLoaded = true;
		}
		return pf_ebrechit_eta_;
	}
	vector<float> &pf_ebrechit_phi()
	{
		if (not pf_ebrechit_phi_isLoaded) {
			if (pf_ebrechit_phi_branch != 0) {
				pf_ebrechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_ebrechit_phi_.begin(); i != pf_ebrechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_ebrechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_ebrechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_ebrechit_phi_isLoaded = true;
		}
		return pf_ebrechit_phi_;
	}
	vector<float> &pf_eecluster_e()
	{
		if (not pf_eecluster_e_isLoaded) {
			if (pf_eecluster_e_branch != 0) {
				pf_eecluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eecluster_e_.begin(); i != pf_eecluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eecluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eecluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_eecluster_e_isLoaded = true;
		}
		return pf_eecluster_e_;
	}
	vector<float> &pf_eecluster_et()
	{
		if (not pf_eecluster_et_isLoaded) {
			if (pf_eecluster_et_branch != 0) {
				pf_eecluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eecluster_et_.begin(); i != pf_eecluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eecluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eecluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_eecluster_et_isLoaded = true;
		}
		return pf_eecluster_et_;
	}
	vector<float> &pf_eecluster_eta()
	{
		if (not pf_eecluster_eta_isLoaded) {
			if (pf_eecluster_eta_branch != 0) {
				pf_eecluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eecluster_eta_.begin(); i != pf_eecluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eecluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eecluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_eecluster_eta_isLoaded = true;
		}
		return pf_eecluster_eta_;
	}
	vector<float> &pf_eecluster_phi()
	{
		if (not pf_eecluster_phi_isLoaded) {
			if (pf_eecluster_phi_branch != 0) {
				pf_eecluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eecluster_phi_.begin(); i != pf_eecluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eecluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eecluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_eecluster_phi_isLoaded = true;
		}
		return pf_eecluster_phi_;
	}
	vector<float> &pf_eerechit_e()
	{
		if (not pf_eerechit_e_isLoaded) {
			if (pf_eerechit_e_branch != 0) {
				pf_eerechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eerechit_e_.begin(); i != pf_eerechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eerechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eerechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_eerechit_e_isLoaded = true;
		}
		return pf_eerechit_e_;
	}
	vector<float> &pf_eerechit_et()
	{
		if (not pf_eerechit_et_isLoaded) {
			if (pf_eerechit_et_branch != 0) {
				pf_eerechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eerechit_et_.begin(); i != pf_eerechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eerechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eerechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_eerechit_et_isLoaded = true;
		}
		return pf_eerechit_et_;
	}
	vector<float> &pf_eerechit_eta()
	{
		if (not pf_eerechit_eta_isLoaded) {
			if (pf_eerechit_eta_branch != 0) {
				pf_eerechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eerechit_eta_.begin(); i != pf_eerechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eerechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eerechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_eerechit_eta_isLoaded = true;
		}
		return pf_eerechit_eta_;
	}
	vector<float> &pf_eerechit_phi()
	{
		if (not pf_eerechit_phi_isLoaded) {
			if (pf_eerechit_phi_branch != 0) {
				pf_eerechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_eerechit_phi_.begin(); i != pf_eerechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_eerechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_eerechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_eerechit_phi_isLoaded = true;
		}
		return pf_eerechit_phi_;
	}
	vector<float> &pf_hbcluster_e()
	{
		if (not pf_hbcluster_e_isLoaded) {
			if (pf_hbcluster_e_branch != 0) {
				pf_hbcluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbcluster_e_.begin(); i != pf_hbcluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbcluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbcluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_hbcluster_e_isLoaded = true;
		}
		return pf_hbcluster_e_;
	}
	vector<float> &pf_hbcluster_et()
	{
		if (not pf_hbcluster_et_isLoaded) {
			if (pf_hbcluster_et_branch != 0) {
				pf_hbcluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbcluster_et_.begin(); i != pf_hbcluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbcluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbcluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_hbcluster_et_isLoaded = true;
		}
		return pf_hbcluster_et_;
	}
	vector<float> &pf_hbcluster_eta()
	{
		if (not pf_hbcluster_eta_isLoaded) {
			if (pf_hbcluster_eta_branch != 0) {
				pf_hbcluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbcluster_eta_.begin(); i != pf_hbcluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbcluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbcluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hbcluster_eta_isLoaded = true;
		}
		return pf_hbcluster_eta_;
	}
	vector<float> &pf_hbcluster_phi()
	{
		if (not pf_hbcluster_phi_isLoaded) {
			if (pf_hbcluster_phi_branch != 0) {
				pf_hbcluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbcluster_phi_.begin(); i != pf_hbcluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbcluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbcluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hbcluster_phi_isLoaded = true;
		}
		return pf_hbcluster_phi_;
	}
	vector<float> &pf_hbrechit_e()
	{
		if (not pf_hbrechit_e_isLoaded) {
			if (pf_hbrechit_e_branch != 0) {
				pf_hbrechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbrechit_e_.begin(); i != pf_hbrechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbrechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbrechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_hbrechit_e_isLoaded = true;
		}
		return pf_hbrechit_e_;
	}
	vector<float> &pf_hbrechit_et()
	{
		if (not pf_hbrechit_et_isLoaded) {
			if (pf_hbrechit_et_branch != 0) {
				pf_hbrechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbrechit_et_.begin(); i != pf_hbrechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbrechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbrechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_hbrechit_et_isLoaded = true;
		}
		return pf_hbrechit_et_;
	}
	vector<float> &pf_hbrechit_eta()
	{
		if (not pf_hbrechit_eta_isLoaded) {
			if (pf_hbrechit_eta_branch != 0) {
				pf_hbrechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbrechit_eta_.begin(); i != pf_hbrechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbrechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbrechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hbrechit_eta_isLoaded = true;
		}
		return pf_hbrechit_eta_;
	}
	vector<float> &pf_hbrechit_phi()
	{
		if (not pf_hbrechit_phi_isLoaded) {
			if (pf_hbrechit_phi_branch != 0) {
				pf_hbrechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hbrechit_phi_.begin(); i != pf_hbrechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hbrechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hbrechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hbrechit_phi_isLoaded = true;
		}
		return pf_hbrechit_phi_;
	}
	vector<float> &pf_hecluster_e()
	{
		if (not pf_hecluster_e_isLoaded) {
			if (pf_hecluster_e_branch != 0) {
				pf_hecluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hecluster_e_.begin(); i != pf_hecluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hecluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hecluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_hecluster_e_isLoaded = true;
		}
		return pf_hecluster_e_;
	}
	vector<float> &pf_hecluster_et()
	{
		if (not pf_hecluster_et_isLoaded) {
			if (pf_hecluster_et_branch != 0) {
				pf_hecluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hecluster_et_.begin(); i != pf_hecluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hecluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hecluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_hecluster_et_isLoaded = true;
		}
		return pf_hecluster_et_;
	}
	vector<float> &pf_hecluster_eta()
	{
		if (not pf_hecluster_eta_isLoaded) {
			if (pf_hecluster_eta_branch != 0) {
				pf_hecluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hecluster_eta_.begin(); i != pf_hecluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hecluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hecluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hecluster_eta_isLoaded = true;
		}
		return pf_hecluster_eta_;
	}
	vector<float> &pf_hecluster_phi()
	{
		if (not pf_hecluster_phi_isLoaded) {
			if (pf_hecluster_phi_branch != 0) {
				pf_hecluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hecluster_phi_.begin(); i != pf_hecluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hecluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hecluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hecluster_phi_isLoaded = true;
		}
		return pf_hecluster_phi_;
	}
	vector<float> &pf_herechit_e()
	{
		if (not pf_herechit_e_isLoaded) {
			if (pf_herechit_e_branch != 0) {
				pf_herechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_herechit_e_.begin(); i != pf_herechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_herechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_herechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_herechit_e_isLoaded = true;
		}
		return pf_herechit_e_;
	}
	vector<float> &pf_herechit_et()
	{
		if (not pf_herechit_et_isLoaded) {
			if (pf_herechit_et_branch != 0) {
				pf_herechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_herechit_et_.begin(); i != pf_herechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_herechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_herechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_herechit_et_isLoaded = true;
		}
		return pf_herechit_et_;
	}
	vector<float> &pf_herechit_eta()
	{
		if (not pf_herechit_eta_isLoaded) {
			if (pf_herechit_eta_branch != 0) {
				pf_herechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_herechit_eta_.begin(); i != pf_herechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_herechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_herechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_herechit_eta_isLoaded = true;
		}
		return pf_herechit_eta_;
	}
	vector<float> &pf_herechit_phi()
	{
		if (not pf_herechit_phi_isLoaded) {
			if (pf_herechit_phi_branch != 0) {
				pf_herechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_herechit_phi_.begin(); i != pf_herechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_herechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_herechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_herechit_phi_isLoaded = true;
		}
		return pf_herechit_phi_;
	}
	vector<float> &pf_hfecluster_e()
	{
		if (not pf_hfecluster_e_isLoaded) {
			if (pf_hfecluster_e_branch != 0) {
				pf_hfecluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfecluster_e_.begin(); i != pf_hfecluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfecluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfecluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_hfecluster_e_isLoaded = true;
		}
		return pf_hfecluster_e_;
	}
	vector<float> &pf_hfecluster_et()
	{
		if (not pf_hfecluster_et_isLoaded) {
			if (pf_hfecluster_et_branch != 0) {
				pf_hfecluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfecluster_et_.begin(); i != pf_hfecluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfecluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfecluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_hfecluster_et_isLoaded = true;
		}
		return pf_hfecluster_et_;
	}
	vector<float> &pf_hfecluster_eta()
	{
		if (not pf_hfecluster_eta_isLoaded) {
			if (pf_hfecluster_eta_branch != 0) {
				pf_hfecluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfecluster_eta_.begin(); i != pf_hfecluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfecluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfecluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hfecluster_eta_isLoaded = true;
		}
		return pf_hfecluster_eta_;
	}
	vector<float> &pf_hfecluster_phi()
	{
		if (not pf_hfecluster_phi_isLoaded) {
			if (pf_hfecluster_phi_branch != 0) {
				pf_hfecluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfecluster_phi_.begin(); i != pf_hfecluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfecluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfecluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hfecluster_phi_isLoaded = true;
		}
		return pf_hfecluster_phi_;
	}
	vector<float> &pf_hferechit_e()
	{
		if (not pf_hferechit_e_isLoaded) {
			if (pf_hferechit_e_branch != 0) {
				pf_hferechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hferechit_e_.begin(); i != pf_hferechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hferechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hferechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_hferechit_e_isLoaded = true;
		}
		return pf_hferechit_e_;
	}
	vector<float> &pf_hferechit_et()
	{
		if (not pf_hferechit_et_isLoaded) {
			if (pf_hferechit_et_branch != 0) {
				pf_hferechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hferechit_et_.begin(); i != pf_hferechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hferechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hferechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_hferechit_et_isLoaded = true;
		}
		return pf_hferechit_et_;
	}
	vector<float> &pf_hferechit_eta()
	{
		if (not pf_hferechit_eta_isLoaded) {
			if (pf_hferechit_eta_branch != 0) {
				pf_hferechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hferechit_eta_.begin(); i != pf_hferechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hferechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hferechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hferechit_eta_isLoaded = true;
		}
		return pf_hferechit_eta_;
	}
	vector<float> &pf_hferechit_phi()
	{
		if (not pf_hferechit_phi_isLoaded) {
			if (pf_hferechit_phi_branch != 0) {
				pf_hferechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hferechit_phi_.begin(); i != pf_hferechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hferechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hferechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hferechit_phi_isLoaded = true;
		}
		return pf_hferechit_phi_;
	}
	vector<float> &pf_hfhcluster_e()
	{
		if (not pf_hfhcluster_e_isLoaded) {
			if (pf_hfhcluster_e_branch != 0) {
				pf_hfhcluster_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhcluster_e_.begin(); i != pf_hfhcluster_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhcluster_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhcluster_e_branch does not exist!\n");
				exit(1);
			}
			pf_hfhcluster_e_isLoaded = true;
		}
		return pf_hfhcluster_e_;
	}
	vector<float> &pf_hfhcluster_et()
	{
		if (not pf_hfhcluster_et_isLoaded) {
			if (pf_hfhcluster_et_branch != 0) {
				pf_hfhcluster_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhcluster_et_.begin(); i != pf_hfhcluster_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhcluster_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhcluster_et_branch does not exist!\n");
				exit(1);
			}
			pf_hfhcluster_et_isLoaded = true;
		}
		return pf_hfhcluster_et_;
	}
	vector<float> &pf_hfhcluster_eta()
	{
		if (not pf_hfhcluster_eta_isLoaded) {
			if (pf_hfhcluster_eta_branch != 0) {
				pf_hfhcluster_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhcluster_eta_.begin(); i != pf_hfhcluster_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhcluster_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhcluster_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hfhcluster_eta_isLoaded = true;
		}
		return pf_hfhcluster_eta_;
	}
	vector<float> &pf_hfhcluster_phi()
	{
		if (not pf_hfhcluster_phi_isLoaded) {
			if (pf_hfhcluster_phi_branch != 0) {
				pf_hfhcluster_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhcluster_phi_.begin(); i != pf_hfhcluster_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhcluster_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhcluster_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hfhcluster_phi_isLoaded = true;
		}
		return pf_hfhcluster_phi_;
	}
	vector<float> &pf_hfhrechit_e()
	{
		if (not pf_hfhrechit_e_isLoaded) {
			if (pf_hfhrechit_e_branch != 0) {
				pf_hfhrechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhrechit_e_.begin(); i != pf_hfhrechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhrechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhrechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_hfhrechit_e_isLoaded = true;
		}
		return pf_hfhrechit_e_;
	}
	vector<float> &pf_hfhrechit_et()
	{
		if (not pf_hfhrechit_et_isLoaded) {
			if (pf_hfhrechit_et_branch != 0) {
				pf_hfhrechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhrechit_et_.begin(); i != pf_hfhrechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhrechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhrechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_hfhrechit_et_isLoaded = true;
		}
		return pf_hfhrechit_et_;
	}
	vector<float> &pf_hfhrechit_eta()
	{
		if (not pf_hfhrechit_eta_isLoaded) {
			if (pf_hfhrechit_eta_branch != 0) {
				pf_hfhrechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhrechit_eta_.begin(); i != pf_hfhrechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhrechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhrechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_hfhrechit_eta_isLoaded = true;
		}
		return pf_hfhrechit_eta_;
	}
	vector<float> &pf_hfhrechit_phi()
	{
		if (not pf_hfhrechit_phi_isLoaded) {
			if (pf_hfhrechit_phi_branch != 0) {
				pf_hfhrechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_hfhrechit_phi_.begin(); i != pf_hfhrechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_hfhrechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_hfhrechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_hfhrechit_phi_isLoaded = true;
		}
		return pf_hfhrechit_phi_;
	}
	vector<float> &pf_rechit_detid()
	{
		if (not pf_rechit_detid_isLoaded) {
			if (pf_rechit_detid_branch != 0) {
				pf_rechit_detid_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_rechit_detid_.begin(); i != pf_rechit_detid_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_rechit_detid_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_rechit_detid_branch does not exist!\n");
				exit(1);
			}
			pf_rechit_detid_isLoaded = true;
		}
		return pf_rechit_detid_;
	}
	vector<float> &pf_rechit_e()
	{
		if (not pf_rechit_e_isLoaded) {
			if (pf_rechit_e_branch != 0) {
				pf_rechit_e_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_rechit_e_.begin(); i != pf_rechit_e_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_rechit_e_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_rechit_e_branch does not exist!\n");
				exit(1);
			}
			pf_rechit_e_isLoaded = true;
		}
		return pf_rechit_e_;
	}
	vector<float> &pf_rechit_et()
	{
		if (not pf_rechit_et_isLoaded) {
			if (pf_rechit_et_branch != 0) {
				pf_rechit_et_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_rechit_et_.begin(); i != pf_rechit_et_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_rechit_et_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_rechit_et_branch does not exist!\n");
				exit(1);
			}
			pf_rechit_et_isLoaded = true;
		}
		return pf_rechit_et_;
	}
	vector<float> &pf_rechit_eta()
	{
		if (not pf_rechit_eta_isLoaded) {
			if (pf_rechit_eta_branch != 0) {
				pf_rechit_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_rechit_eta_.begin(); i != pf_rechit_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_rechit_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_rechit_eta_branch does not exist!\n");
				exit(1);
			}
			pf_rechit_eta_isLoaded = true;
		}
		return pf_rechit_eta_;
	}
	vector<float> &pf_rechit_phi()
	{
		if (not pf_rechit_phi_isLoaded) {
			if (pf_rechit_phi_branch != 0) {
				pf_rechit_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pf_rechit_phi_.begin(); i != pf_rechit_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pf_rechit_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pf_rechit_phi_branch does not exist!\n");
				exit(1);
			}
			pf_rechit_phi_isLoaded = true;
		}
		return pf_rechit_phi_;
	}
};

#ifndef __CINT__
extern CMS2 cms2;
#endif

namespace tas {
	float &calo_eb_sumet();
	float &calo_ee_sumet();
	float &calo_hb_sumet();
	float &calo_he_sumet();
	float &calo_hfe_sumet();
	float &calo_hfh_sumet();
	float &calomet();
	float &calosumet();
	float &genmet();
	float &gensumet();
	float &pfclus_eb_met();
	float &pfclus_eb_sumet();
	float &pfclus_ee_met();
	float &pfclus_ee_sumet();
	float &pfclus_hb_met();
	float &pfclus_hb_sumet();
	float &pfclus_he_met();
	float &pfclus_he_sumet();
	float &pfclus_hfe_met();
	float &pfclus_hfe_sumet();
	float &pfclus_hfh_met();
	float &pfclus_hfh_sumet();
	float &pfclusmet();
	float &pfclussumet();
	float &pf_eb_met();
	float &pf_eb_sumet();
	float &pf_ee_met();
	float &pf_ee_sumet();
	float &pf_hb_met();
	float &pf_hb_sumet();
	float &pf_he_met();
	float &pf_he_sumet();
	float &pf_hfe_met();
	float &pf_hfe_sumet();
	float &pf_hfh_met();
	float &pf_hfh_sumet();
	float &pfmet();
	float &pfsumet();
	vector<float> &pf_cluster_detid();
	vector<float> &pf_cluster_e();
	vector<float> &pf_cluster_et();
	vector<float> &pf_cluster_eta();
	vector<float> &pf_cluster_phi();
	vector<float> &pf_ebcluster_e();
	vector<float> &pf_ebcluster_et();
	vector<float> &pf_ebcluster_eta();
	vector<float> &pf_ebcluster_phi();
	vector<float> &pf_ebrechit_e();
	vector<float> &pf_ebrechit_et();
	vector<float> &pf_ebrechit_eta();
	vector<float> &pf_ebrechit_phi();
	vector<float> &pf_eecluster_e();
	vector<float> &pf_eecluster_et();
	vector<float> &pf_eecluster_eta();
	vector<float> &pf_eecluster_phi();
	vector<float> &pf_eerechit_e();
	vector<float> &pf_eerechit_et();
	vector<float> &pf_eerechit_eta();
	vector<float> &pf_eerechit_phi();
	vector<float> &pf_hbcluster_e();
	vector<float> &pf_hbcluster_et();
	vector<float> &pf_hbcluster_eta();
	vector<float> &pf_hbcluster_phi();
	vector<float> &pf_hbrechit_e();
	vector<float> &pf_hbrechit_et();
	vector<float> &pf_hbrechit_eta();
	vector<float> &pf_hbrechit_phi();
	vector<float> &pf_hecluster_e();
	vector<float> &pf_hecluster_et();
	vector<float> &pf_hecluster_eta();
	vector<float> &pf_hecluster_phi();
	vector<float> &pf_herechit_e();
	vector<float> &pf_herechit_et();
	vector<float> &pf_herechit_eta();
	vector<float> &pf_herechit_phi();
	vector<float> &pf_hfecluster_e();
	vector<float> &pf_hfecluster_et();
	vector<float> &pf_hfecluster_eta();
	vector<float> &pf_hfecluster_phi();
	vector<float> &pf_hferechit_e();
	vector<float> &pf_hferechit_et();
	vector<float> &pf_hferechit_eta();
	vector<float> &pf_hferechit_phi();
	vector<float> &pf_hfhcluster_e();
	vector<float> &pf_hfhcluster_et();
	vector<float> &pf_hfhcluster_eta();
	vector<float> &pf_hfhcluster_phi();
	vector<float> &pf_hfhrechit_e();
	vector<float> &pf_hfhrechit_et();
	vector<float> &pf_hfhrechit_eta();
	vector<float> &pf_hfhrechit_phi();
	vector<float> &pf_rechit_detid();
	vector<float> &pf_rechit_e();
	vector<float> &pf_rechit_et();
	vector<float> &pf_rechit_eta();
	vector<float> &pf_rechit_phi();
}
#endif
