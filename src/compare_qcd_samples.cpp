//Makes comparison plots between high-stat QCD MC with only 2 gen. b quarks and low-stat QCD MC allowing more b quarks

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "style.hpp"
#include "timer.hpp"
#include "utils.hpp"

void compare(std::vector<TH1D>& h, const std::string& out_name){
  TCanvas canvas;
  h.at(0).SetLineColor(1);
  h.at(1).SetLineColor(2);
  h.at(0).SetStats(0);
  h.at(1).SetStats(0);
  h.at(0).Draw("hist");
  h.at(1).Draw("histsame");
  double the_max(1.1*std::max(get_maximum(h.at(0)), get_maximum(h.at(1))));
  h.at(0).SetMaximum(the_max);
  h.at(1).SetMaximum(the_max);
  TLegend legend(0.5, 0.7, 0.9, 0.85);
  legend.AddEntry(&h.at(0), "BJets", "lpe");
  legend.AddEntry(&h.at(1), "HT-binned", "lpe");
  legend.Draw("same");
  canvas.SetLogy(0);
  canvas.Print((out_name+".pdf").c_str());
  canvas.SetLogy(1);
  canvas.Print((out_name+"_log.pdf").c_str());
  normalize(h.at(0));
  normalize(h.at(1));
  the_max=1.1*std::max(get_maximum(h.at(0)), get_maximum(h.at(1)));
  h.at(0).SetMaximum(the_max);
  h.at(1).SetMaximum(the_max);
  h.at(0).Draw("hist");
  h.at(1).Draw("histsame");
  legend.Draw("same");
  canvas.SetLogy(0);
  canvas.Print((out_name+"_an.pdf").c_str());
  canvas.SetLogy(1);
  canvas.Print((out_name+"_an_log.pdf").c_str());
}

int main(){
  SetStyle();
  TChain qcd_2b("qcd_2b", "qcd_2b");
  TChain qcd_4b("qcd_4b", "qcd_4b");

  qcd_2b.Add("reduced_trees/BJets*1.root/reduced_tree");
  qcd_4b.Add("reduced_trees/QCD_HT*1.root/reduced_tree");

  std::cout << "TESTING:" << std::endl;
  double count(0.0), uncert(0.0);
  get_count_and_uncertainty(qcd_4b, "full_weight*(passesJSONCut)", count, uncert);
  get_count_and_uncertainty(qcd_4b, "full_weight*(passesJSONCut)", count, uncert);
  std::cout << count << " " << uncert << std::endl;

  std::vector<TChain*> chains(0);
  chains.push_back(&qcd_2b);
  chains.push_back(&qcd_4b);

  std::vector<std::string> names(0);

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    passes2CSVTCut(false), passesMETSig30Cut(false), passesMETCleaningCut(false),
    passesTriggerCut(false), passesNumJetsCut(false), passesMinDeltaPhiCut(false),
    passesLeptonVetoCut(false), passesIsoTrackVetoCut(false), passesDRCut(false),
    passesBTaggingCut(false), passesHiggsAvgMassCut(false),
    passesHiggsMassDiffCut(false), passesQCDTriggerCut(false);

  bool higgs_mass_signal_region(false), higgs_mass_sideband(false);

  bool passesBaselineSelection(false);

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
  short chargino_mass(0);

  std::vector<TH1D> h_num_b_tagged_jets(0);
  std::vector<TH1D> h_num_jets(0);
  std::vector<TH1D> h_met(0);
  std::vector<TH1D> h_met_sig(0);
  std::vector<TH1D> h_average_higgs_mass(0);
  std::vector<TH1D> h_higgs_mass_difference(0);
  std::vector<TH1D> h_max_delta_R(0);
  std::vector<TH1D> h_min_delta_phi(0);

  for(unsigned chain_num(0); chain_num<chains.size(); ++chain_num){
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
    setup(chain, "higgs_mass_signal_region", higgs_mass_signal_region);
    setup(chain, "higgs_mass_sideband", higgs_mass_sideband);
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
    setup(chain, "passesBaselineSelection", passesBaselineSelection);
    setup(chain, "chargino_mass", chargino_mass);

    const std::string name(chain.GetName());
    names.push_back(name);
    h_num_b_tagged_jets.push_back(TH1D(("h_num_b_tagged_jets"+name).c_str(), "Num. b-Tagged Jets;Num. b-Tagged Jets;Events/1", 16, -0.5, 15.5));
    h_num_jets.push_back(TH1D(("h_num_jets"+name).c_str(), "Num. Jets;Num. Tagged Jets;Events/1", 16, -0.5, 15.5));
    h_met.push_back(TH1D(("h_met"+name).c_str(), "MET;MET [GeV];Events/25 GeV", 20, 0., 500.0));
    h_met_sig.push_back(TH1D(("h_met_sig"+name).c_str(), "S_{MET};S_{MET};Events/10", 25, 0.0, 250.0));
    h_average_higgs_mass.push_back(TH1D(("h_average_higgs_mass"+name).c_str(), "<m_{bb}>;<m_{bb}> [GeV];Events/10 GeV", 20, 0.0, 400.0));
    h_higgs_mass_difference.push_back(TH1D(("h_higgs_mass_difference"+name).c_str(), "#Delta m_{bb};#Delta m_{bb} [GeV];Events/10 GeV", 20, 0.0, 400.0));
    h_max_delta_R.push_back(TH1D(("h_max_delta_R"+name).c_str(), "#Delta R_{max};#Delta R_{max};Events/0.5", 20, 0.0, 10.0));
    h_min_delta_phi.push_back(TH1D(("h_min_delta_phi"+name).c_str(), "#Delta #phi_{min};#Delta #phi_{min};Events/0.2", 16, 0.0, 1.6));
    const int num_events(chain.GetEntries());
    Timer timer(num_events);
    timer.Start();
    for(int event(0); event<num_events; ++event){
      if(event%(1u<<16u)==0){
        timer.PrintRemainingTime();
      }
      chain.GetEntry(event);
      timer.Iterate();

      if(chain_num==1){
	higgs_mass_signal_region=(average_higgs_mass>100.0 && average_higgs_mass<140.0 && higgs_mass_difference<20.0);
	higgs_mass_sideband=(average_higgs_mass<90.0 || average_higgs_mass>150.0 || higgs_mass_difference>30.0);
	num_leptons=num_electrons+num_muons+num_taus;
	top_pt=-999999.999999;
      }

      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut /*&& passesNumJetsCut && passesMinDeltaPhiCut*/ && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
	h_num_jets.at(chain_num).Fill(num_jets, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut /*&& passes2CSVTCut*/ && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut /*&& passesMinDeltaPhiCut*/ && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
	h_num_b_tagged_jets.at(chain_num).Fill(num_b_tagged_jets, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut /*&& passesMETSig30Cut*/ && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut /*&& passesMinDeltaPhiCut*/ && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
	h_met.at(chain_num).Fill(met, full_weight);
	h_met_sig.at(chain_num).Fill(met_sig, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut /*&& passesMinDeltaPhiCut*/ && passesLeptonVetoCut && passesIsoTrackVetoCut /*&& passesDRCut*/){
	h_max_delta_R.at(chain_num).Fill(max_delta_R, full_weight);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut /*&& passesMinDeltaPhiCut*/ && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
	h_min_delta_phi.at(chain_num).Fill(min_delta_phi, full_weight);
	h_average_higgs_mass.at(chain_num).Fill(average_higgs_mass, full_weight);
	h_higgs_mass_difference.at(chain_num).Fill(higgs_mass_difference, full_weight);
      }
    }
  }
  compare(h_num_b_tagged_jets, "crap/num_b_tagged_jets");
  compare(h_num_jets, "crap/num_jets");
  compare(h_met, "crap/met");
  compare(h_met_sig, "crap/met_sig");
  compare(h_average_higgs_mass, "crap/average_higgs_mass");
  compare(h_higgs_mass_difference, "crap/higgs_mass_difference");
  compare(h_max_delta_R, "crap/max_delta_R");
  compare(h_min_delta_phi, "crap/min_delta_phi");
}
