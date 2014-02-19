#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "style.hpp"
#include "timer.hpp"
#include "utils.hpp"

TH1D GetMCSum(std::vector<TH1D>& histos){
  if(histos.size()<2){
    return TH1D();
  }else{
    TH1D histo(histos.at(1));
    for(unsigned h(2); h<histos.size(); ++h){
      histo=histo+histos.at(h);
    }
    return histo;
  }
}

void GetHighLowRatio(TH1D& h, const double low, const double high,
                     double& ratio, double& uncert){
  const int bin_low(h.FindBin(low-0.000001)), bin_high(h.FindBin(high));
  const int num_bins(h.GetNbinsX());
  double uncert_high(0.0), uncert_low(0.0);
  const double integ_low(h.IntegralAndError(0, bin_low, uncert_low));
  const double integ_high(h.IntegralAndError(bin_high, num_bins+1, uncert_high));
  ratio=integ_high/integ_low;
  const double high_sq(integ_high*integ_high);
  const double low_sq(integ_low*integ_low);
  const double u_high_sq(uncert_high*uncert_high);
  const double u_low_sq(uncert_low*uncert_low);
  uncert=std::sqrt((high_sq*u_low_sq+low_sq*u_high_sq))/low_sq;
}

int main(){
  SetStyle();
  TH1::SetDefaultSumw2();

  TChain data("data","data");
  TChain ttbar("ttbar","ttbar");
  TChain qcd("qcd","qcd");
  TChain single_t_or_boson("single_t_or_boson","single_t_or_boson");
  TChain diboson("diboson","diboson");

  data.Add("reduced_trees/MET_*2012*1.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_FullLept*1.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_SemiLept*1.root/reduced_tree");
  qcd.Add("reduced_trees/TTJets_Hadronic*1.root/reduced_tree");
  qcd.Add("reduced_trees/BJets*1.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*channel*1.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*JetsToLNu_Tune*1.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/ZJetsToNuNu*1.root/reduced_tree");
  diboson.Add("reduced_trees/WH*1.root/reduced_tree");
  diboson.Add("reduced_trees/WW*1.root/reduced_tree");
  diboson.Add("reduced_trees/WZ*1.root/reduced_tree");
  diboson.Add("reduced_trees/ZZ*1.root/reduced_tree");
  std::vector<TChain*> chains(0);
  chains.push_back(&data);
  chains.push_back(&qcd);
  chains.push_back(&diboson);
  chains.push_back(&single_t_or_boson);
  chains.push_back(&ttbar);

  std::vector<std::string> names(0);

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    passes2CSVTCut(false), passesMETSig30Cut(false), passesMETCleaningCut(false),
    passesTriggerCut(false), passesNumJetsCut(false), passesMinDeltaPhiCut(false),
    passesLeptonVetoCut(false), passesIsoTrackVetoCut(false), passesDRCut(false),
    passesBTaggingCut(false), passesHiggsAvgMassCut(false),
    passesHiggsMassDiffCut(false), passesQCDTriggerCut(false);

  float pu_true_num_interactions(0.0);
  unsigned short num_primary_vertices(0);

  float highest_jet_pt(0.0), second_highest_jet_pt(0.0), third_highest_jet_pt(0.0),
    fourth_highest_jet_pt(0.0), fifth_highest_jet_pt(0.0);
  float highest_csv(0.0), second_highest_csv(0.0),
    third_highest_csv(0.0), fourth_highest_csv(0.0), fifth_highest_csv(0.0);
  float met_sig(0.0), met(0.0);
  unsigned short num_jets(0), num_b_tagged_jets(0);
  float min_delta_phi(0.0);
  unsigned short num_electrons(0), num_muons(0), num_taus(0), num_iso_tracks(0),
    num_leptons(0);
  float min_delta_R(0.0), max_delta_R(0.0);
  float average_higgs_mass(0.0), higgs_mass_difference(0.0);
  float ht_jets(0.0), ht_jets_met(0.0), ht_jets_leps(0.0), ht_jets_met_leps(0.0);
  float full_weight(0.0);
  float top_pt(0.0);

  std::vector<TH1D> h_pu_true_num_interactions(0);
  std::vector<TH1D> h_num_primary_vertices(0);
  std::vector<TH1D> h_highest_jet_pt(0);
  std::vector<TH1D> h_second_highest_jet_pt(0);
  std::vector<TH1D> h_third_highest_jet_pt(0);
  std::vector<TH1D> h_fourth_highest_jet_pt(0);
  std::vector<TH1D> h_fifth_highest_jet_pt(0);
  std::vector<TH1D> h_highest_csv(0);
  std::vector<TH1D> h_second_highest_csv(0);
  std::vector<TH1D> h_third_highest_csv(0);
  std::vector<TH1D> h_fourth_highest_csv(0);
  std::vector<TH1D> h_fifth_highest_csv(0);
  std::vector<TH1D> h_met_sig(0);
  std::vector<TH1D> h_met(0);
  std::vector<TH1D> h_num_jets(0);
  std::vector<TH1D> h_num_b_tagged_jets(0);
  std::vector<TH1D> h_min_delta_phi(0);
  std::vector<TH1D> h_num_electrons(0);
  std::vector<TH1D> h_num_muons(0);
  std::vector<TH1D> h_num_taus(0);
  std::vector<TH1D> h_num_iso_tracks(0);
  std::vector<TH1D> h_num_leptons(0);
  std::vector<TH1D> h_min_delta_R(0);
  std::vector<TH1D> h_max_delta_R(0);
  std::vector<TH1D> h_average_higgs_mass(0);
  std::vector<TH1D> h_higgs_mass_difference(0);
  std::vector<TH1D> h_ht_jets(0);
  std::vector<TH1D> h_ht_jets_met(0);
  std::vector<TH1D> h_ht_jets_leps(0);
  std::vector<TH1D> h_ht_jets_met_leps(0);
  std::vector<TH1D> h_top_pt(0);
  std::vector<TH1D> h_qcd_control_met_sig(0);
  std::vector<TH1D> h_qcd_control_met(0);
  std::vector<TH1D> h_qcd_control_highest_csv(0);
  std::vector<TH1D> h_qcd_control_second_highest_csv(0);
  std::vector<TH1D> h_qcd_control_third_highest_csv(0);
  std::vector<TH1D> h_qcd_control_fourth_highest_csv(0);
  std::vector<TH1D> h_qcd_control_fifth_highest_csv(0);
  std::vector<TH1D> h_1l_met_sig(0);
  std::vector<TH1D> h_1l_met(0);
  std::vector<TH1D> h_1l_highest_csv(0);
  std::vector<TH1D> h_1l_second_highest_csv(0);
  std::vector<TH1D> h_1l_third_highest_csv(0);
  std::vector<TH1D> h_1l_fourth_highest_csv(0);
  std::vector<TH1D> h_1l_fifth_highest_csv(0);

  for(unsigned chain_num(0); chain_num<chains.size(); ++chain_num){
    std::cout << chain_num << "/" << chains.size() << std::endl;
    TChain& chain(*chains.at(chain_num));
    chain.SetBranchStatus("*",0);
    setup(chain, "passesJSONCut", passesJSONCut);
    setup(chain, "passesPVCut", passesPVCut);
    setup(chain, "passesJet2PtCut", passesJet2PtCut);
    setup(chain, "passes2CSVTCut", passes2CSVTCut);
    setup(chain, "passesMETSig30Cut", passesMETSig30Cut);
    setup(chain, "passesMETCleaningCut", passesMETCleaningCut);
    setup(chain, "passesTriggerCut", passesTriggerCut);
    setup(chain, "passesNumJetsCut", passesNumJetsCut);
    setup(chain, "passesMinDeltaPhiCut", passesMinDeltaPhiCut);
    setup(chain, "passesLeptonVetoCut", passesLeptonVetoCut);
    setup(chain, "passesIsoTrackVetoCut", passesIsoTrackVetoCut);
    setup(chain, "passesDRCut", passesDRCut);
    setup(chain, "passesBTaggingCut", passesBTaggingCut);
    setup(chain, "passesHiggsAvgMassCut", passesHiggsAvgMassCut);
    setup(chain, "passesHiggsMassDiffCut", passesHiggsMassDiffCut);
    setup(chain, "passesQCDTriggerCut", passesQCDTriggerCut);
    setup(chain, "pu_true_num_interactions", pu_true_num_interactions);
    setup(chain, "num_primary_vertices", num_primary_vertices);
    setup(chain, "highest_jet_pt", highest_jet_pt);
    setup(chain, "second_highest_jet_pt", second_highest_jet_pt);
    setup(chain, "third_highest_jet_pt", third_highest_jet_pt);
    setup(chain, "fourth_highest_jet_pt", fourth_highest_jet_pt);
    setup(chain, "fifth_highest_jet_pt", fifth_highest_jet_pt);
    setup(chain, "highest_csv", highest_csv);
    setup(chain, "second_highest_csv", second_highest_csv);
    setup(chain, "third_highest_csv", third_highest_csv);
    setup(chain, "fourth_highest_csv", fourth_highest_csv);
    setup(chain, "fifth_highest_csv", fifth_highest_csv);
    setup(chain, "met_sig", met_sig);
    setup(chain, "met", met);
    setup(chain, "num_jets", num_jets);
    setup(chain, "num_b_tagged_jets", num_b_tagged_jets);
    setup(chain, "min_delta_phi", min_delta_phi);
    setup(chain, "num_electrons", num_electrons);
    setup(chain, "num_muons", num_muons);
    setup(chain, "num_taus", num_taus);
    setup(chain, "num_iso_tracks", num_iso_tracks);
    setup(chain, "num_leptons", num_leptons);
    setup(chain, "min_delta_R", min_delta_R);
    setup(chain, "max_delta_R", max_delta_R);
    setup(chain, "average_higgs_mass", average_higgs_mass);
    setup(chain, "higgs_mass_difference", higgs_mass_difference);
    setup(chain, "ht_jets", ht_jets);
    setup(chain, "ht_jets_met", ht_jets_met);
    setup(chain, "ht_jets_leps", ht_jets_leps);
    setup(chain, "ht_jets_met_leps", ht_jets_met_leps);
    setup(chain, "top_pt", top_pt);
    setup(chain, "full_weight", full_weight);

    const std::string name(chain.GetName());
    names.push_back(name);
    h_pu_true_num_interactions.push_back(TH1D(("h_pu_true_num_interactions"+name).c_str(), "PU True Num Interactions (baseline);Interactions;Events/1", 61, -0.5, 60.5));
    h_num_primary_vertices.push_back(TH1D(("h_num_primary_vertices"+name).c_str(), "Num Primary Vertices (baseline);Vertices;Events/1", 61, -0.5, 60.5));
    h_highest_jet_pt.push_back(TH1D(("h_highest_jet_pt"+name).c_str(), "Highest Jet Pt (baseline);Highest Jet Pt [GeV];Events/10 GeV", 50, 0.0, 500.0));
    h_second_highest_jet_pt.push_back(TH1D(("h_second_highest_jet_pt"+name).c_str(), "Second Highest Jet Pt (baseline);Second Highest Jet Pt [GeV];Events/10 GeV", 50, 0.0, 500.0));
    h_third_highest_jet_pt.push_back(TH1D(("h_third_highest_jet_pt"+name).c_str(), "Third Highest Jet Pt (baseline);Third Highest Jet Pt[GeV];Events/10 GeV", 50, 0.0, 500.0));
    h_fourth_highest_jet_pt.push_back(TH1D(("h_fourth_highest_jet_pt"+name).c_str(), "Fourth Highest Jet Pt (baseline);Fourth Highest Jet Pt [GeV];Events/10 GeV", 50, 0.0, 500.0));
    h_fifth_highest_jet_pt.push_back(TH1D(("h_fifth_highest_jet_pt"+name).c_str(), "Fifth Highest Jet Pt (baseline);Fifth Highest Jet Pt [GeV];Events/10 GeV", 50, 0.0, 500.0));
    h_highest_csv.push_back(TH1D(("h_highest_csv"+name).c_str(), "Highest CSV (baseline);Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_second_highest_csv.push_back(TH1D(("h_second_highest_csv"+name).c_str(), "Second Highest CSV (baseline);Second Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_third_highest_csv.push_back(TH1D(("h_third_highest_csv"+name).c_str(), "Third Highest CSV (baseline);Third Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_fourth_highest_csv.push_back(TH1D(("h_fourth_highest_csv"+name).c_str(), "Fourth Highest CSV (baseline);Fourth Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_fifth_highest_csv.push_back(TH1D(("h_fifth_highest_csv"+name).c_str(), "Fifth Highest CSV (baseline);Fifth Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_met_sig.push_back(TH1D(("h_met_sig"+name).c_str(), "S_{MET} (baseline);S_{MET};Events/10", 40, 0.0, 400.0));
    h_met.push_back(TH1D(("h_met"+name).c_str(), "MET (baseline);MET [GeV];Events/10 GeV", 40, 0.0, 400.0));
    h_num_jets.push_back(TH1D(("h_num_jets"+name).c_str(), "Num Jets (baseline);Num Jets;Events/1", 16, -0.5, 15.5));
    h_num_b_tagged_jets.push_back(TH1D(("h_num_b_tagged_jets"+name).c_str(), "Num b-Tagged Jets (baseline);Num b-Tagged Jets;Events/1", 16, -0.5, 15.5));
    h_min_delta_phi.push_back(TH1D(("h_min_delta_phi"+name).c_str(), "Min Delta Phi (baseline);Min Delta Phi;Events/0.1", 32, 0.0, 3.2));
    h_num_electrons.push_back(TH1D(("h_num_electrons"+name).c_str(), "Num Electrons (baseline);Num Electrons;Events/1", 16, -0.5, 15.5));
    h_num_muons.push_back(TH1D(("h_num_muons"+name).c_str(), "Num Muons (baseline);Num Muons;Events/1", 16, -0.5, 15.5));
    h_num_taus.push_back(TH1D(("h_num_taus"+name).c_str(), "Num Taus (baseline);Num Taus;Events/1", 16, -0.5, 15.5));
    h_num_iso_tracks.push_back(TH1D(("h_num_iso_tracks"+name).c_str(), "Num Iso Tracks (baseline);Num Iso Tracks;Events/1", 16, -0.5, 15.5));
    h_num_leptons.push_back(TH1D(("h_num_leptons"+name).c_str(), "Num Leptons (baseline);Num Leptons;Events/1", 16, -0.5, 15.5));
    h_min_delta_R.push_back(TH1D(("h_min_delta_R"+name).c_str(), "Min Delta R (baseline);Min Delta R;Events/0.1", 50, 0.0, 5.0));
    h_max_delta_R.push_back(TH1D(("h_max_delta_R"+name).c_str(), "Max Delta R (baseline);Max Delta R;Events/0.1", 50, 0.0, 5.0));
    h_average_higgs_mass.push_back(TH1D(("h_average_higgs_mass"+name).c_str(), "<m_{bb}> (baseline);<m_{bb}> [GeV];Events/5 GeV", 50, 0.0, 250.0));
    h_higgs_mass_difference.push_back(TH1D(("h_higgs_mass_difference"+name).c_str(), "#Delta m_{bb} (baseline);#Delta m_{bb} [GeV];Events/5 GeV", 50, 0.0, 250.0));
    h_ht_jets.push_back(TH1D(("h_ht_jets"+name).c_str(), "H_{T} (jets) (baseline);H_{T} [GeV];Events/10 GeV", 100, 0.0, 1000.0));
    h_ht_jets_met.push_back(TH1D(("h_ht_jets_met"+name).c_str(), "H_{T} (jets+MET) (baseline);H_{T} [GeV];Events/10 GeV", 100, 0.0, 1000.0));
    h_ht_jets_leps.push_back(TH1D(("h_ht_jets_leps"+name).c_str(), "H_{T} (jets+leps) (baseline);H_{T} [GeV];Events/10 GeV", 100, 0.0, 1000.0));
    h_ht_jets_met_leps.push_back(TH1D(("h_ht_jets_met_leps"+name).c_str(), "H_{T} (jets+leps+MET) (baseline);H_{T} [GeV];Events/10 GeV", 100, 0.0, 1000.0));
    h_top_pt.push_back(TH1D(("h_top_pt"+name).c_str(), "Top Pt (baseline);Top Pt [GeV];Events/10 GeV", 50, 0.0, 500.0));
    h_qcd_control_highest_csv.push_back(TH1D(("h_qcd_control_highest_csv"+name).c_str(), "Highest CSV (QCD control);Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_qcd_control_second_highest_csv.push_back(TH1D(("h_qcd_control_second_highest_csv"+name).c_str(), "Second Highest CSV (QCD control);Second Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_qcd_control_third_highest_csv.push_back(TH1D(("h_qcd_control_third_highest_csv"+name).c_str(), "Third Highest CSV (QCD control);Third Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_qcd_control_fourth_highest_csv.push_back(TH1D(("h_qcd_control_fourth_highest_csv"+name).c_str(), "Fourth Highest CSV (QCD control);Fourth Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_qcd_control_fifth_highest_csv.push_back(TH1D(("h_qcd_control_fifth_highest_csv"+name).c_str(), "Fifth Highest CSV (QCD control);Fifth Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_qcd_control_met_sig.push_back(TH1D(("h_qcd_control_met_sig"+name).c_str(), "S_{MET} (QCD control);S_{MET};Events/10", 40, 0.0, 400.0));
    h_qcd_control_met.push_back(TH1D(("h_qcd_control_met"+name).c_str(), "MET (QCD control);MET [GeV];Events/10 GeV", 40, 0.0, 400.0));    
    h_1l_highest_csv.push_back(TH1D(("h_1l_highest_csv"+name).c_str(), "Highest CSV (1l control);Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_1l_second_highest_csv.push_back(TH1D(("h_1l_second_highest_csv"+name).c_str(), "Second Highest CSV (1l control);Second Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_1l_third_highest_csv.push_back(TH1D(("h_1l_third_highest_csv"+name).c_str(), "Third Highest CSV (1l control);Third Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_1l_fourth_highest_csv.push_back(TH1D(("h_1l_fourth_highest_csv"+name).c_str(), "Fourth Highest CSV (1l control);Fourth Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_1l_fifth_highest_csv.push_back(TH1D(("h_1l_fifth_highest_csv"+name).c_str(), "Fifth Highest CSV (1l control);Fifth Highest CSV;Events/0.02", 50, 0.0, 1.0));
    h_1l_met_sig.push_back(TH1D(("h_1l_met_sig"+name).c_str(), "S_{MET} (1l control);S_{MET};Events/10", 40, 0.0, 400.0));
    h_1l_met.push_back(TH1D(("h_1l_met"+name).c_str(), "MET (1l control);MET [GeV];Events/10 GeV", 40, 0.0, 400.0));
    
    const int num_events(chain.GetEntries());
    Timer timer(num_events);
    timer.Start();
    for(int event(0); event<num_events; ++event){
      if(event%(1u<<16u)==0){
        timer.PrintRemainingTime();
      }
      chain.GetEntry(event);

      if(passesJSONCut && passesPVCut && second_highest_jet_pt>70.0
         /*&& passes2CSVTCut*/ && !passesMETSig30Cut && passesMETCleaningCut
         /*&& passesTriggerCut*/ && passesNumJetsCut && !passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut
	 && num_b_tagged_jets<3 && met>160.0 && passesQCDTriggerCut){
	h_qcd_control_highest_csv.at(chain_num).Fill(highest_csv, full_weight);
	h_qcd_control_second_highest_csv.at(chain_num).Fill(second_highest_csv, full_weight);
	h_qcd_control_third_highest_csv.at(chain_num).Fill(third_highest_csv, full_weight);
	h_qcd_control_fourth_highest_csv.at(chain_num).Fill(fourth_highest_csv, full_weight);
	h_qcd_control_fifth_highest_csv.at(chain_num).Fill(fifth_highest_csv, full_weight);
      }
      if(passesJSONCut && passesPVCut && second_highest_jet_pt>70.0
         /*&& passes2CSVTCut && !passesMETSig30Cut*/ && passesMETCleaningCut
         /*&& passesTriggerCut*/ && passesNumJetsCut && !passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut
	 && num_b_tagged_jets<3 && met>160.0 && passesQCDTriggerCut){
	h_qcd_control_met_sig.at(chain_num).Fill(met_sig, full_weight);
      }
      if(passesJSONCut && passesPVCut && second_highest_jet_pt>70.0
         /*&& passes2CSVTCut*/ && !passesMETSig30Cut && passesMETCleaningCut
         /*&& passesTriggerCut*/ && passesNumJetsCut && !passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut
	 && num_b_tagged_jets<3 /*&& met>160.0*/ && passesQCDTriggerCut){
	h_qcd_control_met.at(chain_num).Fill(met, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut /*&& passesMETSig30Cut*/ && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         /*&& passesLeptonVetoCut && passesIsoTrackVetoCut*/ && passesDRCut
	 && num_leptons==1 && num_taus==0){
	h_1l_met_sig.at(chain_num).Fill(met_sig, full_weight);
	h_1l_met.at(chain_num).Fill(met, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         /*&& passesLeptonVetoCut && passesIsoTrackVetoCut*/ && passesDRCut
	 && num_leptons==1 && num_taus==0){
	h_1l_highest_csv.at(chain_num).Fill(highest_csv, full_weight);
	h_1l_second_highest_csv.at(chain_num).Fill(second_highest_csv, full_weight);
	h_1l_third_highest_csv.at(chain_num).Fill(third_highest_csv, full_weight);
	h_1l_fourth_highest_csv.at(chain_num).Fill(fourth_highest_csv, full_weight);
	h_1l_fifth_highest_csv.at(chain_num).Fill(fifth_highest_csv, full_weight);
      }

      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
        h_pu_true_num_interactions.at(chain_num).Fill(pu_true_num_interactions, full_weight);
        h_num_primary_vertices.at(chain_num).Fill(num_primary_vertices, full_weight);
        h_third_highest_jet_pt.at(chain_num).Fill(third_highest_jet_pt, full_weight);
        h_fourth_highest_jet_pt.at(chain_num).Fill(fourth_highest_jet_pt, full_weight);
        h_fifth_highest_jet_pt.at(chain_num).Fill(fifth_highest_jet_pt, full_weight);
        h_third_highest_csv.at(chain_num).Fill(third_highest_csv, full_weight);
        h_fourth_highest_csv.at(chain_num).Fill(fourth_highest_csv, full_weight);
        h_fifth_highest_csv.at(chain_num).Fill(fifth_highest_csv, full_weight);
        h_min_delta_R.at(chain_num).Fill(min_delta_R, full_weight);
        h_ht_jets.at(chain_num).Fill(ht_jets, full_weight);
        h_ht_jets_met.at(chain_num).Fill(ht_jets_met, full_weight);
        h_ht_jets_leps.at(chain_num).Fill(ht_jets_leps, full_weight);
        h_ht_jets_met_leps.at(chain_num).Fill(ht_jets_met_leps, full_weight);
        h_top_pt.at(chain_num).Fill(top_pt, full_weight);
        h_average_higgs_mass.at(chain_num).Fill(average_higgs_mass, full_weight);
        h_higgs_mass_difference.at(chain_num).Fill(higgs_mass_difference, full_weight);
      }
      if(passesJSONCut && passesPVCut /*&& passesJet2PtCut*/
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
        h_highest_jet_pt.at(chain_num).Fill(highest_jet_pt, full_weight);
        h_second_highest_jet_pt.at(chain_num).Fill(second_highest_jet_pt, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         /*&& passes2CSVTCut*/ && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
        h_highest_csv.at(chain_num).Fill(highest_csv, full_weight);
        h_second_highest_csv.at(chain_num).Fill(second_highest_csv, full_weight);
        h_num_b_tagged_jets.at(chain_num).Fill(num_b_tagged_jets, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && /*passesMETSig30Cut &&*/ passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
        h_met_sig.at(chain_num).Fill(met_sig, full_weight);
        h_met.at(chain_num).Fill(met, full_weight);
      } 
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && /*passesNumJetsCut &&*/ passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut /*&& passesDRCut*/){
        h_num_jets.at(chain_num).Fill(num_jets, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut /*&& passesMinDeltaPhiCut*/
         && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
        h_min_delta_phi.at(chain_num).Fill(min_delta_phi, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         /*&& passesLeptonVetoCut && passesIsoTrackVetoCut*/ && passesDRCut){
        h_num_electrons.at(chain_num).Fill(num_electrons, full_weight);
        h_num_muons.at(chain_num).Fill(num_muons, full_weight);
        h_num_taus.at(chain_num).Fill(num_taus, full_weight);
        h_num_iso_tracks.at(chain_num).Fill(num_iso_tracks, full_weight);
        h_num_leptons.at(chain_num).Fill(num_leptons, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
         && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
         && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
         && passesLeptonVetoCut && passesIsoTrackVetoCut /*&& passesDRCut*/){
        h_max_delta_R.at(chain_num).Fill(max_delta_R, full_weight);
      }
      timer.Iterate();
    }
  }

  MakeRatioPlot(h_pu_true_num_interactions, names, "nm1/pu_true_num_interactions.pdf");
  MakeRatioPlot(h_num_primary_vertices, names, "nm1/num_primary_vertices.pdf");
  MakeRatioPlot(h_highest_jet_pt, names, "nm1/highest_jet_pt.pdf");
  MakeRatioPlot(h_second_highest_jet_pt, names, "nm1/second_highest_jet_pt.pdf");
  MakeRatioPlot(h_third_highest_jet_pt, names, "nm1/third_highest_jet_pt.pdf");
  MakeRatioPlot(h_fourth_highest_jet_pt, names, "nm1/fourth_highest_jet_pt.pdf");
  MakeRatioPlot(h_fifth_highest_jet_pt, names, "nm1/fifth_highest_jet_pt.pdf");
  MakeRatioPlot(h_highest_csv, names, "nm1/highest_csv.pdf");
  MakeRatioPlot(h_second_highest_csv, names, "nm1/second_highest_csv.pdf");
  MakeRatioPlot(h_third_highest_csv, names, "nm1/third_highest_csv.pdf");
  MakeRatioPlot(h_fourth_highest_csv, names, "nm1/fourth_highest_csv.pdf");
  MakeRatioPlot(h_fifth_highest_csv, names, "nm1/fifth_highest_csv.pdf");
  MakeRatioPlot(h_met_sig, names, "nm1/met_sig.pdf");
  MakeRatioPlot(h_met, names, "nm1/met.pdf");
  MakeRatioPlot(h_num_jets, names, "nm1/num_jets.pdf");
  MakeRatioPlot(h_num_b_tagged_jets, names, "nm1/num_b_tagged_jets.pdf");
  MakeRatioPlot(h_min_delta_phi, names, "nm1/min_delta_phi.pdf");
  MakeRatioPlot(h_num_electrons, names, "nm1/num_electrons.pdf");
  MakeRatioPlot(h_num_muons, names, "nm1/num_muons.pdf");
  MakeRatioPlot(h_num_taus, names, "nm1/num_taus.pdf");
  MakeRatioPlot(h_num_iso_tracks, names, "nm1/num_iso_tracks.pdf");
  MakeRatioPlot(h_num_leptons, names, "nm1/num_leptons.pdf");
  MakeRatioPlot(h_min_delta_R, names, "nm1/min_delta_R.pdf");
  MakeRatioPlot(h_max_delta_R, names, "nm1/max_delta_R.pdf");
  MakeRatioPlot(h_average_higgs_mass, names, "nm1/average_higgs_mass.pdf");
  MakeRatioPlot(h_higgs_mass_difference, names, "nm1/higgs_mass_difference.pdf");
  MakeRatioPlot(h_ht_jets, names, "nm1/ht_jets.pdf");
  MakeRatioPlot(h_ht_jets_met, names, "nm1/ht_jets_met.pdf");
  MakeRatioPlot(h_ht_jets_leps, names, "nm1/ht_jets_leps.pdf");
  MakeRatioPlot(h_ht_jets_met_leps, names, "nm1/ht_jets_met_leps.pdf");
  MakeRatioPlot(h_top_pt, names, "nm1/top_pt.pdf");
  MakeRatioPlot(h_qcd_control_highest_csv, names, "nm1/qcd_control_highest_csv.pdf");
  MakeRatioPlot(h_qcd_control_second_highest_csv, names, "nm1/qcd_control_second_highest_csv.pdf");
  MakeRatioPlot(h_qcd_control_third_highest_csv, names, "nm1/qcd_control_third_highest_csv.pdf");
  MakeRatioPlot(h_qcd_control_fourth_highest_csv, names, "nm1/qcd_control_fourth_highest_csv.pdf");
  MakeRatioPlot(h_qcd_control_fifth_highest_csv, names, "nm1/qcd_control_fifth_highest_csv.pdf");
  MakeRatioPlot(h_qcd_control_met_sig, names, "nm1/qcd_control_met_sig.pdf");
  MakeRatioPlot(h_qcd_control_met, names, "nm1/qcd_control_met.pdf");
  MakeRatioPlot(h_1l_highest_csv, names, "nm1/1l_highest_csv.pdf");
  MakeRatioPlot(h_1l_second_highest_csv, names, "nm1/1l_second_highest_csv.pdf");
  MakeRatioPlot(h_1l_third_highest_csv, names, "nm1/1l_third_highest_csv.pdf");
  MakeRatioPlot(h_1l_fourth_highest_csv, names, "nm1/1l_fourth_highest_csv.pdf");
  MakeRatioPlot(h_1l_fifth_highest_csv, names, "nm1/1l_fifth_highest_csv.pdf");
  MakeRatioPlot(h_1l_met_sig, names, "nm1/1l_met_sig.pdf");
  MakeRatioPlot(h_1l_met, names, "nm1/1l_met.pdf");

  //TH1D tot_min_delta_phi(GetMCSum(h_min_delta_phi));
  TH1D tot_min_delta_phi(h_min_delta_phi.at(0));
  const int the_bin(tot_min_delta_phi.FindBin(0.3));
  const int num_bins(tot_min_delta_phi.GetNbinsX());
  const double low_integral(tot_min_delta_phi.Integral(0, the_bin-1));
  const double high_integral(tot_min_delta_phi.Integral(the_bin, num_bins+1));
  std::cout << "Low integral: " << low_integral << std::endl;
  std::cout << "High integral: " << high_integral << std::endl;
}
