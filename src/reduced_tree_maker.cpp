#include "reduced_tree_maker.hpp"
#include <vector>
#include <string>
#include <set>
#include "timer.hpp"
#include "event_handler.hpp"

ReducedTreeMaker::ReducedTreeMaker(const std::string& in_file_name,
                                   const bool is_list,
                                   const double weight_in):
  EventHandler(in_file_name, is_list, weight_in, false){
}

void ReducedTreeMaker::MakeReducedTree(const std::string& out_file_name){
  TFile file(out_file_name.c_str(), "recreate");
  std::set<EventNumber> eventList;
  file.cd();

  const bool isRealData(sampleName.find("Run2012")!=std::string::npos);
  const bool isttbar(sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos);
  std::vector<float> dataDist(pu::RunsThrough203002, pu::RunsThrough203002+60);
  std::vector<float> MCDist(pu::Summer2012_S10, pu::Summer2012_S10+60);//QQQ this needs to change later for general pileup scenario
  reweight::LumiReWeighting lumiWeights(MCDist, dataDist);

  TTree reduced_tree("reduced_tree","reduced_tree");
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

  reduced_tree.Branch("passesJSONCut", &passesJSONCut);
  reduced_tree.Branch("passesPVCut",&passesPVCut);
  reduced_tree.Branch("passesJet2PtCut",&passesJet2PtCut);
  reduced_tree.Branch("passes2CSVTCut",&passes2CSVTCut);
  reduced_tree.Branch("passesMETSig30Cut",&passesMETSig30Cut);
  reduced_tree.Branch("passesMETCleaningCut",&passesMETCleaningCut);
  reduced_tree.Branch("passesTriggerCut",&passesTriggerCut);
  reduced_tree.Branch("passesQCDTriggerCut",&passesQCDTriggerCut);
  reduced_tree.Branch("passesSingleMuTrigger",&passesSingleMuTrigger);
  reduced_tree.Branch("passesNumJetsCut",&passesNumJetsCut);
  reduced_tree.Branch("passesMinDeltaPhiCut",&passesMinDeltaPhiCut);
  reduced_tree.Branch("passesLeptonVetoCut",&passesLeptonVetoCut);
  reduced_tree.Branch("passesIsoTrackVetoCut",&passesIsoTrackVetoCut);
  reduced_tree.Branch("passesDRCut",&passesDRCut);
  reduced_tree.Branch("passesBTaggingCut",&passesBTaggingCut);
  reduced_tree.Branch("passesHiggsAvgMassCut",&passesHiggsAvgMassCut);
  reduced_tree.Branch("passesHiggsMassDiffCut",&passesHiggsMassDiffCut);

  reduced_tree.Branch("higgs_mass_signal_region",&higgs_mass_signal_region);
  reduced_tree.Branch("higgs_mass_sideband",&higgs_mass_sideband);
  reduced_tree.Branch("passes4bSignalRegionCut",&passes4bSignalRegionCut);
  reduced_tree.Branch("passes4bSidebandRegionCut",&passes4bSidebandRegionCut);
  reduced_tree.Branch("passes3bSignalRegionCut",&passes3bSignalRegionCut);
  reduced_tree.Branch("passes3bSidebandRegionCut",&passes3bSidebandRegionCut);
  reduced_tree.Branch("passes2bSignalRegionCut",&passes2bSignalRegionCut);
  reduced_tree.Branch("passes2bSidebandRegionCut",&passes2bSidebandRegionCut);

  reduced_tree.Branch("passesBaselineSelection",&passesBaselineSelection);
  reduced_tree.Branch("passesTriggerPlateauCuts",&passesTriggerPlateauCuts);

  reduced_tree.Branch("controlSampleOneLepton",&controlSampleOneLepton);
  reduced_tree.Branch("controlSampleQCD",&controlSampleQCD);

  reduced_tree.Branch("highest_jet_pt", &highest_jet_pt);
  reduced_tree.Branch("second_highest_jet_pt", &second_highest_jet_pt);
  reduced_tree.Branch("third_highest_jet_pt", &third_highest_jet_pt);
  reduced_tree.Branch("fourth_highest_jet_pt", &fourth_highest_jet_pt);
  reduced_tree.Branch("fifth_highest_jet_pt", &fifth_highest_jet_pt);

  reduced_tree.Branch("highest_csv", &highest_csv);
  reduced_tree.Branch("second_highest_csv", &second_highest_csv);
  reduced_tree.Branch("third_highest_csv", &third_highest_csv);
  reduced_tree.Branch("fourth_highest_csv", &fourth_highest_csv);
  reduced_tree.Branch("fifth_highest_csv", &fifth_highest_csv);

  reduced_tree.Branch("pu_true_num_interactions", &pu_true_num_interactions);
  reduced_tree.Branch("num_primary_vertices", &num_primary_vertices);

  reduced_tree.Branch("met_sig", &met_sig);
  reduced_tree.Branch("met", &met);

  reduced_tree.Branch("num_jets", &num_jets);
  reduced_tree.Branch("num_b_tagged_jets", &num_b_tagged_jets);

  reduced_tree.Branch("min_delta_phi", &min_delta_phi);

  reduced_tree.Branch("num_electrons", &num_electrons);
  reduced_tree.Branch("num_muons", &num_muons);
  reduced_tree.Branch("num_taus", &num_taus);
  reduced_tree.Branch("num_leptons", &num_leptons);

  reduced_tree.Branch("num_noiso_electrons", &num_noiso_electrons);
  reduced_tree.Branch("num_noiso_muons", &num_noiso_muons);
  reduced_tree.Branch("num_noiso_taus", &num_noiso_taus);
  reduced_tree.Branch("num_noiso_leptons", &num_noiso_leptons);

  reduced_tree.Branch("num_loose_electrons", &num_loose_electrons);
  reduced_tree.Branch("num_loose_muons", &num_loose_muons);
  reduced_tree.Branch("num_loose_taus", &num_loose_taus);
  reduced_tree.Branch("num_loose_leptons", &num_loose_leptons);

  reduced_tree.Branch("num_medium_electrons", &num_medium_electrons);
  reduced_tree.Branch("num_medium_muons", &num_medium_muons);
  reduced_tree.Branch("num_medium_taus", &num_medium_taus);
  reduced_tree.Branch("num_medium_leptons", &num_medium_leptons);

  reduced_tree.Branch("num_tight_electrons", &num_tight_electrons);
  reduced_tree.Branch("num_tight_muons", &num_tight_muons);
  reduced_tree.Branch("num_tight_taus", &num_tight_taus);
  reduced_tree.Branch("num_tight_leptons", &num_tight_leptons);

  reduced_tree.Branch("num_iso_tracks", &num_iso_tracks);

  reduced_tree.Branch("min_delta_R", &min_delta_r);
  reduced_tree.Branch("max_delta_R", &max_delta_r);

  reduced_tree.Branch("average_higgs_mass", &average_higgs_mass);
  reduced_tree.Branch("higgs_mass_difference", &higgs_mass_difference);

  reduced_tree.Branch("ht_jets", &ht_jets);
  reduced_tree.Branch("ht_jets_met", &ht_jets_met);
  reduced_tree.Branch("ht_jets_leps", &ht_jets_leps);
  reduced_tree.Branch("ht_jets_met_leps", &ht_jets_met_leps);

  reduced_tree.Branch("has_gluon_splitting", &has_gluon_splitting);
  reduced_tree.Branch("top_pt",&top_pt);

  reduced_tree.Branch("full_weight", &full_weight);
  reduced_tree.Branch("lumi_weight", &lumi_weight);
  reduced_tree.Branch("pu_weight", &pu_weight);
  reduced_tree.Branch("top_pt_weight", &top_pt_weight);
  reduced_tree.Branch("top_pt_weight_official", &top_pt_weight_official);
  reduced_tree.Branch("trigger_weight", &trigger_weight);
 
  reduced_tree.Branch("lsp_mass", &lsp_mass);
  reduced_tree.Branch("chargino_mass", &chargino_mass);

  reduced_tree.Branch("run", &run);
  reduced_tree.Branch("event", &event);
  reduced_tree.Branch("lumiblock", &lumiblock);

  Timer timer(GetTotalEntries());
  timer.Start();
  for(int i(0); i<GetTotalEntries(); ++i){
    if(i%1000==0 && i!=0){
      timer.PrintRemainingTime();
    }
    timer.Iterate();
    GetEntry(i);

    std::pair<std::set<EventNumber>::iterator, bool> returnVal(eventList.insert(EventNumber(run, event, lumiblock)));
    if(!returnVal.second) continue;
    
    // Saving our cuts for the reduced tree
    passesTriggerPlateauCuts=PassesTriggerPlateauCuts();
    passesBaselineSelection=PassesBaselineSelection();
    passesJSONCut=PassesJSONCut();
    passesPVCut=PassesPVCut();
    passesJet2PtCut=PassesJet2PtCut();
    passes2CSVTCut=Passes2CSVTCut();
    passesMETSig30Cut=PassesMETSig30Cut();
    passesMETCleaningCut=PassesMETCleaningCut();
    passesTriggerCut=PassesTriggerCut();
    passesNumJetsCut=PassesNumJetsCut();
    passesMinDeltaPhiCut=PassesMinDeltaPhiCut();
    passesLeptonVetoCut=PassesLeptonVetoCut();
    passesIsoTrackVetoCut=PassesIsoTrackVetoCut();
    passesDRCut=PassesDRCut();
    passesBTaggingCut=PassesBTaggingCut();
    passesHiggsAvgMassCut=PassesHiggsAvgMassCut();
    passesHiggsMassDiffCut=PassesHiggsMassDiffCut();

    passesQCDTriggerCut=PassesQCDTriggerCut();
    passesSingleMuTrigger=PassesSpecificTrigger("HLT_IsoMu24_eta2p1");

    passes4bSignalRegionCut=PassesRegionACut();
    passes4bSidebandRegionCut=PassesRegionBCut();
    passes3bSignalRegionCut=PassesRegionC3bCut();
    passes3bSidebandRegionCut=PassesRegionD3bCut();
    passes2bSignalRegionCut=PassesRegionC2bCut();
    passes2bSidebandRegionCut=PassesRegionD2bCut();

    highest_jet_pt=GetHighestJetPt(1);
    second_highest_jet_pt=GetHighestJetPt(2);
    third_highest_jet_pt=GetHighestJetPt(3);
    fourth_highest_jet_pt=GetHighestJetPt(4);
    fifth_highest_jet_pt=GetHighestJetPt(5);

    highest_csv=GetHighestJetCSV(1);
    second_highest_csv=GetHighestJetCSV(2);
    third_highest_csv=GetHighestJetCSV(3);
    fourth_highest_csv=GetHighestJetCSV(4);
    fifth_highest_csv=GetHighestJetCSV(5);

    pu_true_num_interactions=GetNumInteractions();
    num_primary_vertices=GetNumVertices();

    met_sig=pfmets_fullSignif;
    met=pfTypeImets_et->at(0);

    num_jets=GetNumGoodJets();
    num_b_tagged_jets=GetNumBTaggedJets();

    min_delta_phi=GetMinDeltaPhiMET();

    num_electrons=GetNumElectrons(0);
    num_muons=GetNumMuons(0);
    num_taus=GetNumTaus(0);
    num_leptons=num_electrons+num_muons+num_taus;
    num_noiso_electrons=GetNumElectrons(0,false);
    num_noiso_muons=GetNumMuons(0,false);
    num_noiso_taus=GetNumTaus(0,false);
    num_noiso_leptons=num_noiso_electrons+num_noiso_muons+num_noiso_taus;
    num_loose_electrons=GetNumElectrons(1);
    num_loose_muons=GetNumMuons(1);
    num_loose_taus=GetNumTaus(1);
    num_loose_leptons=num_loose_electrons+num_loose_muons+num_loose_taus;
    num_medium_electrons=GetNumElectrons(2);
    num_medium_muons=GetNumMuons(2);
    num_medium_taus=GetNumTaus(2);
    num_medium_leptons=num_medium_electrons+num_medium_muons+num_medium_taus;
    num_tight_electrons=GetNumElectrons(3);
    num_tight_muons=GetNumMuons(3);
    num_tight_taus=GetNumTaus(3);
    num_tight_leptons=num_tight_electrons+num_tight_muons+num_tight_taus;
    num_iso_tracks=NewGetNumIsoTracks();

    min_delta_r=GetMinDR();
    max_delta_r=GetMaxDR();

    controlSampleOneLepton=PassesSingleLeptonControlCut();
    controlSampleQCD=PassesQCDControlCut();

    const std::pair<double, double> higgs_masses(GetHiggsMasses());
    average_higgs_mass=0.5*(higgs_masses.first+higgs_masses.second);
    higgs_mass_difference=fabs(higgs_masses.first-higgs_masses.second);

    higgs_mass_signal_region=passesHiggsAvgMassCut&&passesHiggsMassDiffCut;
    higgs_mass_sideband=(average_higgs_mass<90||average_higgs_mass>150)||higgs_mass_difference>30;

    ht_jets=GetHT(false, false);
    ht_jets_met=GetHT(true, false);
    ht_jets_leps=GetHT(false, true);
    ht_jets_met_leps=GetHT(true, true);

    double this_scale_factor(scaleFactor);
    if(sampleName.find("SMS-TChiHH")!=std::string::npos){
      this_scale_factor=look_up_scale_factor();
    }

    pu_weight=isRealData?1.0:GetPUWeight(lumiWeights);
    lumi_weight=this_scale_factor;
    top_pt_weight=isttbar?GetTopPtWeight():1.0;
    top_pt_weight_official=isttbar?GetTopPtWeightOfficial():1.0;
    trigger_weight=GetSbinWeight();
    full_weight=pu_weight*lumi_weight*top_pt_weight_official*trigger_weight;

    top_pt=GetTopPt();
    has_gluon_splitting=HasGluonSplitting();

    lsp_mass=GetLSPMass();
    chargino_mass=GetCharginoMass();

    reduced_tree.Fill(); 
  }

  reduced_tree.Write();
  file.Close();
}
