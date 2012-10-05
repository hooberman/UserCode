#include "CMS2.h"
CMS2 cms2;
namespace tas {
	float &calo_eb_sumet() { return cms2.calo_eb_sumet(); }
	float &calo_ee_sumet() { return cms2.calo_ee_sumet(); }
	float &calo_hb_sumet() { return cms2.calo_hb_sumet(); }
	float &calo_he_sumet() { return cms2.calo_he_sumet(); }
	float &calo_hfe_sumet() { return cms2.calo_hfe_sumet(); }
	float &calo_hfh_sumet() { return cms2.calo_hfh_sumet(); }
	float &calomet() { return cms2.calomet(); }
	float &calosumet() { return cms2.calosumet(); }
	float &genmet() { return cms2.genmet(); }
	float &gensumet() { return cms2.gensumet(); }
	float &pfclus_eb_met() { return cms2.pfclus_eb_met(); }
	float &pfclus_eb_sumet() { return cms2.pfclus_eb_sumet(); }
	float &pfclus_ee_met() { return cms2.pfclus_ee_met(); }
	float &pfclus_ee_sumet() { return cms2.pfclus_ee_sumet(); }
	float &pfclus_hb_met() { return cms2.pfclus_hb_met(); }
	float &pfclus_hb_sumet() { return cms2.pfclus_hb_sumet(); }
	float &pfclus_he_met() { return cms2.pfclus_he_met(); }
	float &pfclus_he_sumet() { return cms2.pfclus_he_sumet(); }
	float &pfclus_hfe_met() { return cms2.pfclus_hfe_met(); }
	float &pfclus_hfe_sumet() { return cms2.pfclus_hfe_sumet(); }
	float &pfclus_hfh_met() { return cms2.pfclus_hfh_met(); }
	float &pfclus_hfh_sumet() { return cms2.pfclus_hfh_sumet(); }
	float &pfclusmet() { return cms2.pfclusmet(); }
	float &pfclussumet() { return cms2.pfclussumet(); }
	float &pf_eb_met() { return cms2.pf_eb_met(); }
	float &pf_eb_sumet() { return cms2.pf_eb_sumet(); }
	float &pf_ee_met() { return cms2.pf_ee_met(); }
	float &pf_ee_sumet() { return cms2.pf_ee_sumet(); }
	float &pf_hb_met() { return cms2.pf_hb_met(); }
	float &pf_hb_sumet() { return cms2.pf_hb_sumet(); }
	float &pf_he_met() { return cms2.pf_he_met(); }
	float &pf_he_sumet() { return cms2.pf_he_sumet(); }
	float &pf_hfe_met() { return cms2.pf_hfe_met(); }
	float &pf_hfe_sumet() { return cms2.pf_hfe_sumet(); }
	float &pf_hfh_met() { return cms2.pf_hfh_met(); }
	float &pf_hfh_sumet() { return cms2.pf_hfh_sumet(); }
	float &pfmet() { return cms2.pfmet(); }
	float &pfsumet() { return cms2.pfsumet(); }
	vector<float> &pf_cluster_detid() { return cms2.pf_cluster_detid(); }
	vector<float> &pf_cluster_e() { return cms2.pf_cluster_e(); }
	vector<float> &pf_cluster_et() { return cms2.pf_cluster_et(); }
	vector<float> &pf_cluster_eta() { return cms2.pf_cluster_eta(); }
	vector<float> &pf_cluster_phi() { return cms2.pf_cluster_phi(); }
	vector<float> &pf_ebcluster_e() { return cms2.pf_ebcluster_e(); }
	vector<float> &pf_ebcluster_et() { return cms2.pf_ebcluster_et(); }
	vector<float> &pf_ebcluster_eta() { return cms2.pf_ebcluster_eta(); }
	vector<float> &pf_ebcluster_phi() { return cms2.pf_ebcluster_phi(); }
	vector<float> &pf_ebrechit_e() { return cms2.pf_ebrechit_e(); }
	vector<float> &pf_ebrechit_et() { return cms2.pf_ebrechit_et(); }
	vector<float> &pf_ebrechit_eta() { return cms2.pf_ebrechit_eta(); }
	vector<float> &pf_ebrechit_phi() { return cms2.pf_ebrechit_phi(); }
	vector<float> &pf_eecluster_e() { return cms2.pf_eecluster_e(); }
	vector<float> &pf_eecluster_et() { return cms2.pf_eecluster_et(); }
	vector<float> &pf_eecluster_eta() { return cms2.pf_eecluster_eta(); }
	vector<float> &pf_eecluster_phi() { return cms2.pf_eecluster_phi(); }
	vector<float> &pf_eerechit_e() { return cms2.pf_eerechit_e(); }
	vector<float> &pf_eerechit_et() { return cms2.pf_eerechit_et(); }
	vector<float> &pf_eerechit_eta() { return cms2.pf_eerechit_eta(); }
	vector<float> &pf_eerechit_phi() { return cms2.pf_eerechit_phi(); }
	vector<float> &pf_hbcluster_e() { return cms2.pf_hbcluster_e(); }
	vector<float> &pf_hbcluster_et() { return cms2.pf_hbcluster_et(); }
	vector<float> &pf_hbcluster_eta() { return cms2.pf_hbcluster_eta(); }
	vector<float> &pf_hbcluster_phi() { return cms2.pf_hbcluster_phi(); }
	vector<float> &pf_hbrechit_e() { return cms2.pf_hbrechit_e(); }
	vector<float> &pf_hbrechit_et() { return cms2.pf_hbrechit_et(); }
	vector<float> &pf_hbrechit_eta() { return cms2.pf_hbrechit_eta(); }
	vector<float> &pf_hbrechit_phi() { return cms2.pf_hbrechit_phi(); }
	vector<float> &pf_hecluster_e() { return cms2.pf_hecluster_e(); }
	vector<float> &pf_hecluster_et() { return cms2.pf_hecluster_et(); }
	vector<float> &pf_hecluster_eta() { return cms2.pf_hecluster_eta(); }
	vector<float> &pf_hecluster_phi() { return cms2.pf_hecluster_phi(); }
	vector<float> &pf_herechit_e() { return cms2.pf_herechit_e(); }
	vector<float> &pf_herechit_et() { return cms2.pf_herechit_et(); }
	vector<float> &pf_herechit_eta() { return cms2.pf_herechit_eta(); }
	vector<float> &pf_herechit_phi() { return cms2.pf_herechit_phi(); }
	vector<float> &pf_hfecluster_e() { return cms2.pf_hfecluster_e(); }
	vector<float> &pf_hfecluster_et() { return cms2.pf_hfecluster_et(); }
	vector<float> &pf_hfecluster_eta() { return cms2.pf_hfecluster_eta(); }
	vector<float> &pf_hfecluster_phi() { return cms2.pf_hfecluster_phi(); }
	vector<float> &pf_hferechit_e() { return cms2.pf_hferechit_e(); }
	vector<float> &pf_hferechit_et() { return cms2.pf_hferechit_et(); }
	vector<float> &pf_hferechit_eta() { return cms2.pf_hferechit_eta(); }
	vector<float> &pf_hferechit_phi() { return cms2.pf_hferechit_phi(); }
	vector<float> &pf_hfhcluster_e() { return cms2.pf_hfhcluster_e(); }
	vector<float> &pf_hfhcluster_et() { return cms2.pf_hfhcluster_et(); }
	vector<float> &pf_hfhcluster_eta() { return cms2.pf_hfhcluster_eta(); }
	vector<float> &pf_hfhcluster_phi() { return cms2.pf_hfhcluster_phi(); }
	vector<float> &pf_hfhrechit_e() { return cms2.pf_hfhrechit_e(); }
	vector<float> &pf_hfhrechit_et() { return cms2.pf_hfhrechit_et(); }
	vector<float> &pf_hfhrechit_eta() { return cms2.pf_hfhrechit_eta(); }
	vector<float> &pf_hfhrechit_phi() { return cms2.pf_hfhrechit_phi(); }
	vector<float> &pf_rechit_detid() { return cms2.pf_rechit_detid(); }
	vector<float> &pf_rechit_e() { return cms2.pf_rechit_e(); }
	vector<float> &pf_rechit_et() { return cms2.pf_rechit_et(); }
	vector<float> &pf_rechit_eta() { return cms2.pf_rechit_eta(); }
	vector<float> &pf_rechit_phi() { return cms2.pf_rechit_phi(); }
}