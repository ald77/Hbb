#include <string>
#include "TFile.h"
#include "TTree.h"
#include "timer.hpp"

int main(int argc, char *argv[]){
  if(argc<2) return 1;
  const unsigned short sel(atoi(argv[1]));
  std::string in_filename("");
  double lumi_weight_in(1.0);
  switch(sel){
  case 0:
    in_filename="reduced_trees/QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1945_v71_SyncSkim";
    lumi_weight_in=19399.*10360000.0/50129518.0;
    break;
  case 1:
    in_filename="reduced_trees/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1946_v71_SyncSkim";
    lumi_weight_in=19399.*276000.0/27062078.0;
    break;
  case 2:
    in_filename="reduced_trees/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1947_v71_SyncSkim";
    lumi_weight_in=19399.*8426.0/30599292.0;
    break;
  case 3:
    in_filename="reduced_trees/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1948_v71_SyncSkim";
    lumi_weight_in=19399.*204.0/13843863.0;
    break;
  default:
    return 2;
  }
  std::string out_filename(in_filename+"_mdpscaled.root");
  in_filename+=".root";

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    passes2CSVTCut(false), passesMETSig30Cut(false), passesMETCleaningCut(false),
    passesTriggerCut(false), passesNumJetsCut(false), passesMinDeltaPhiCut(false),
    passesLeptonVetoCut(false), passesIsoTrackVetoCut(false), passesDRCut(false),
    passesBTaggingCut(false), passesHiggsAvgMassCut(false),
    passesHiggsMassDiffCut(false), passesQCDTriggerCut(false), passesSingleMuTrigger(false);

  bool higgs_mass_signal_region(false), higgs_mass_sideband(false);

  bool passes4bSignalRegionCut(false), passes4bSidebandRegionCut(false);
  bool passes3bSignalRegionCut(false), passes3bSidebandRegionCut(false);
  bool passes2bSignalRegionCut(false), passes2bSidebandRegionCut(false);

  bool passesBaselineSelection(false), passesTriggerPlateauCuts(false);

  bool controlSampleOneLepton(false), controlSampleQCD(false);

  float pu_true_num_interactions(0.0);
  unsigned short num_primary_vertices(0);

  float highest_jet_pt(0.0), second_highest_jet_pt(0.0), third_highest_jet_pt(0.0),
    fourth_highest_jet_pt(0.0), fifth_highest_jet_pt(0.0);
  float highest_csv(0.0), second_highest_csv(0.0),
    third_highest_csv(0.0), fourth_highest_csv(0.0), fifth_highest_csv(0.0);
  float met_sig(0.0), met(0.0);
  unsigned short num_jets(0), num_b_tagged_jets(0);
  float min_delta_phi(0.0);
  unsigned short num_iso_tracks(0);
  unsigned short num_electrons(0), num_muons(0), num_taus(0), num_leptons(0);
  unsigned short num_noiso_electrons(0), num_noiso_muons(0), num_noiso_taus(0), num_noiso_leptons(0);
  unsigned short num_loose_electrons(0), num_loose_muons(0), num_loose_taus(0), num_loose_leptons(0);
  unsigned short num_medium_electrons(0), num_medium_muons(0), num_medium_taus(0), num_medium_leptons(0);
  unsigned short num_tight_electrons(0), num_tight_muons(0), num_tight_taus(0), num_tight_leptons(0);
  float min_delta_r(0.0), max_delta_r(0.0);
  float average_higgs_mass(0.0), higgs_mass_difference(0.0);
  float ht_jets(0.0), ht_jets_met(0.0), ht_jets_leps(0.0), ht_jets_met_leps(0.0);
  float full_weight(0.0), lumi_weight(0.0), top_pt_weight(0.0), top_pt_weight_official(0.0), pu_weight(0.0), trigger_weight(0.0);
  bool has_gluon_splitting(false);
  float top_pt(0.0);
  short lsp_mass(0), chargino_mass(0);

  UInt_t run(0), event(0), lumiblock(0);

  TFile file_old(in_filename.c_str(),"read");
  TFile file_new(out_filename.c_str(),"recreate");
  file_new.cd();
  TTree* reduced_tree_old(static_cast<TTree*>(file_old.Get("reduced_tree")));
  file_new.cd();

  reduced_tree_old->Branch("passesJSONCut", &passesJSONCut);
  reduced_tree_old->Branch("passesPVCut",&passesPVCut);
  reduced_tree_old->Branch("passesJet2PtCut",&passesJet2PtCut);
  reduced_tree_old->Branch("passes2CSVTCut",&passes2CSVTCut);
  reduced_tree_old->Branch("passesMETSig30Cut",&passesMETSig30Cut);
  reduced_tree_old->Branch("passesMETCleaningCut",&passesMETCleaningCut);
  reduced_tree_old->Branch("passesTriggerCut",&passesTriggerCut);
  reduced_tree_old->Branch("passesQCDTriggerCut",&passesQCDTriggerCut);
  reduced_tree_old->Branch("passesSingleMuTrigger",&passesSingleMuTrigger);
  reduced_tree_old->Branch("passesNumJetsCut",&passesNumJetsCut);
  reduced_tree_old->Branch("passesMinDeltaPhiCut",&passesMinDeltaPhiCut);
  reduced_tree_old->Branch("passesLeptonVetoCut",&passesLeptonVetoCut);
  reduced_tree_old->Branch("passesIsoTrackVetoCut",&passesIsoTrackVetoCut);
  reduced_tree_old->Branch("passesDRCut",&passesDRCut);
  reduced_tree_old->Branch("passesBTaggingCut",&passesBTaggingCut);
  reduced_tree_old->Branch("passesHiggsAvgMassCut",&passesHiggsAvgMassCut);
  reduced_tree_old->Branch("passesHiggsMassDiffCut",&passesHiggsMassDiffCut);

  reduced_tree_old->Branch("higgs_mass_signal_region",&higgs_mass_signal_region);
  reduced_tree_old->Branch("higgs_mass_sideband",&higgs_mass_sideband);
  reduced_tree_old->Branch("passes4bSignalRegionCut",&passes4bSignalRegionCut);
  reduced_tree_old->Branch("passes4bSidebandRegionCut",&passes4bSidebandRegionCut);
  reduced_tree_old->Branch("passes3bSignalRegionCut",&passes3bSignalRegionCut);
  reduced_tree_old->Branch("passes3bSidebandRegionCut",&passes3bSidebandRegionCut);
  reduced_tree_old->Branch("passes2bSignalRegionCut",&passes2bSignalRegionCut);
  reduced_tree_old->Branch("passes2bSidebandRegionCut",&passes2bSidebandRegionCut);

  reduced_tree_old->Branch("passesBaselineSelection",&passesBaselineSelection);
  reduced_tree_old->Branch("passesTriggerPlateauCuts",&passesTriggerPlateauCuts);

  reduced_tree_old->Branch("controlSampleOneLepton",&controlSampleOneLepton);
  reduced_tree_old->Branch("controlSampleQCD",&controlSampleQCD);

  reduced_tree_old->Branch("highest_jet_pt", &highest_jet_pt);
  reduced_tree_old->Branch("second_highest_jet_pt", &second_highest_jet_pt);
  reduced_tree_old->Branch("third_highest_jet_pt", &third_highest_jet_pt);
  reduced_tree_old->Branch("fourth_highest_jet_pt", &fourth_highest_jet_pt);
  reduced_tree_old->Branch("fifth_highest_jet_pt", &fifth_highest_jet_pt);

  reduced_tree_old->Branch("highest_csv", &highest_csv);
  reduced_tree_old->Branch("second_highest_csv", &second_highest_csv);
  reduced_tree_old->Branch("third_highest_csv", &third_highest_csv);
  reduced_tree_old->Branch("fourth_highest_csv", &fourth_highest_csv);
  reduced_tree_old->Branch("fifth_highest_csv", &fifth_highest_csv);

  reduced_tree_old->Branch("pu_true_num_interactions", &pu_true_num_interactions);
  reduced_tree_old->Branch("num_primary_vertices", &num_primary_vertices);

  reduced_tree_old->Branch("met_sig", &met_sig);
  reduced_tree_old->Branch("met", &met);

  reduced_tree_old->Branch("num_jets", &num_jets);
  reduced_tree_old->Branch("num_b_tagged_jets", &num_b_tagged_jets);

  reduced_tree_old->Branch("min_delta_phi", &min_delta_phi);

  reduced_tree_old->Branch("num_electrons", &num_electrons);
  reduced_tree_old->Branch("num_muons", &num_muons);
  reduced_tree_old->Branch("num_taus", &num_taus);
  reduced_tree_old->Branch("num_leptons", &num_leptons);

  reduced_tree_old->Branch("num_noiso_electrons", &num_noiso_electrons);
  reduced_tree_old->Branch("num_noiso_muons", &num_noiso_muons);
  reduced_tree_old->Branch("num_noiso_taus", &num_noiso_taus);
  reduced_tree_old->Branch("num_noiso_leptons", &num_noiso_leptons);

  reduced_tree_old->Branch("num_loose_electrons", &num_loose_electrons);
  reduced_tree_old->Branch("num_loose_muons", &num_loose_muons);
  reduced_tree_old->Branch("num_loose_taus", &num_loose_taus);
  reduced_tree_old->Branch("num_loose_leptons", &num_loose_leptons);

  reduced_tree_old->Branch("num_medium_electrons", &num_medium_electrons);
  reduced_tree_old->Branch("num_medium_muons", &num_medium_muons);
  reduced_tree_old->Branch("num_medium_taus", &num_medium_taus);
  reduced_tree_old->Branch("num_medium_leptons", &num_medium_leptons);

  reduced_tree_old->Branch("num_tight_electrons", &num_tight_electrons);
  reduced_tree_old->Branch("num_tight_muons", &num_tight_muons);
  reduced_tree_old->Branch("num_tight_taus", &num_tight_taus);
  reduced_tree_old->Branch("num_tight_leptons", &num_tight_leptons);

  reduced_tree_old->Branch("num_iso_tracks", &num_iso_tracks);

  reduced_tree_old->Branch("min_delta_R", &min_delta_r);
  reduced_tree_old->Branch("max_delta_R", &max_delta_r);

  reduced_tree_old->Branch("average_higgs_mass", &average_higgs_mass);
  reduced_tree_old->Branch("higgs_mass_difference", &higgs_mass_difference);

  reduced_tree_old->Branch("ht_jets", &ht_jets);
  reduced_tree_old->Branch("ht_jets_met", &ht_jets_met);
  reduced_tree_old->Branch("ht_jets_leps", &ht_jets_leps);
  reduced_tree_old->Branch("ht_jets_met_leps", &ht_jets_met_leps);

  reduced_tree_old->Branch("has_gluon_splitting", &has_gluon_splitting);
  reduced_tree_old->Branch("top_pt",&top_pt);

  reduced_tree_old->Branch("full_weight", &full_weight);
  reduced_tree_old->Branch("lumi_weight", &lumi_weight);
  reduced_tree_old->Branch("pu_weight", &pu_weight);
  reduced_tree_old->Branch("top_pt_weight", &top_pt_weight);
  reduced_tree_old->Branch("top_pt_weight_official", &top_pt_weight_official);
  reduced_tree_old->Branch("trigger_weight", &trigger_weight);
 
  reduced_tree_old->Branch("lsp_mass", &lsp_mass);
  reduced_tree_old->Branch("chargino_mass", &chargino_mass);

  reduced_tree_old->Branch("run", &run);
  reduced_tree_old->Branch("event", &event);
  reduced_tree_old->Branch("lumiblock", &lumiblock);

  file_new.cd();
  TTree* reduced_tree_new(reduced_tree_old->CloneTree(0));
  file_new.cd();
  reduced_tree_new->CopyAddresses(reduced_tree_old);
  file_new.cd();

  const unsigned num_entries(reduced_tree_old->GetEntries());
  Timer timer(num_entries);
  timer.Start();
  for(unsigned entry(0); entry<num_entries; ++entry){
    if(!(entry%(1<<14))) timer.PrintRemainingTime();
    reduced_tree_old->GetEntry(entry);
    //lumi_weight=lumi_weight_in;
    lumi_weight*=3.88;
    full_weight=pu_weight*lumi_weight*top_pt_weight_official*trigger_weight;
    reduced_tree_new->Fill();
    timer.Iterate();
  }
  file_new.cd();
  reduced_tree_new->Write();
  file_new.Close();
  file_old.Close();
}
