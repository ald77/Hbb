#include "reduced_tree_maker.hpp"
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
  TTree reduced_tree("reduced_tree","reduced_tree");
  bool passesPVCut(false), passesJet2PtCut(false), passes2CSVTCut(false),
    passesMETSig30Cut(false), passesMETCleaningCut(false), passesTriggerCut(false),
    passesNumJetsCut(false), passesMinDeltaPhiCut(false), passesLeptonVetoCut(false),
    passesIsoTrackVetoCut(false), passesDRCut(false), passesHiggsAvgMassCut(false),
    passesHiggsMassDiffCut(false);

  float highest_jet_pt(0.0), second_highest_jet_pt(0.0), third_highest_jet_pt(0.0),
    fourth_highest_jet_pt(0.0), fifth_highest_jet_pt(0.0);
  float highest_csv(0.0), second_highest_csv(0.0),
    third_highest_csv(0.0), fourth_highest_csv(0.0), fifth_highest_csv(0.0);
  float met_sig(0.0);
  unsigned short num_jets(0), num_b_tagged_jets(0);
  float min_delta_phi(0.0);
  unsigned short num_electrons(0), num_muons(0), num_taus(0), num_iso_tracks(0);
  float min_delta_r(0.0), max_delta_r(0.0);
  float average_higgs_mass(0.0), higgs_mass_difference(0.0);
  
  reduced_tree.Branch("passesPVCut",&passesPVCut,"passesPVCut/O");
  reduced_tree.Branch("passesJet2PtCut",&passesJet2PtCut,"passesJet2PtCut/O");
  reduced_tree.Branch("passes2CSVTCut",&passes2CSVTCut,"passes2CSVTCut/O");
  reduced_tree.Branch("passesMETSig30Cut",&passesMETSig30Cut,"passesMETSig30Cut/O");
  reduced_tree.Branch("passesMETCleaningCut",&passesMETCleaningCut,"passesMETCleaningCut/O");
  reduced_tree.Branch("passesTriggerCut",&passesTriggerCut,"passesTriggerCut/O");
  reduced_tree.Branch("passesNumJetsCut",&passesNumJetsCut,"passesNumJetsCut/O");
  reduced_tree.Branch("passesMinDeltaPhiCut",&passesMinDeltaPhiCut,"passesMinDeltaPhiCut/O");
  reduced_tree.Branch("passesLeptonVetoCut",&passesLeptonVetoCut,"passesLeptonVetoCut/O");
  reduced_tree.Branch("passesIsoTrackVetoCut",&passesIsoTrackVetoCut,"passesIsoTrackVetoCut/O");
  reduced_tree.Branch("passesDRCut",&passesDRCut,"passesDRCut/O");
  reduced_tree.Branch("passesHiggsAvgMassCut",&passesHiggsAvgMassCut,"passesHiggsAvgMassCut/O");
  reduced_tree.Branch("passesHiggsMassDiffCut",&passesHiggsMassDiffCut,"passesHiggsMassDiffCut/O");

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

  reduced_tree.Branch("met_sig", &met_sig);

  reduced_tree.Branch("num_jets", &num_jets);
  reduced_tree.Branch("num_b_tagged_jets", &num_b_tagged_jets);

  reduced_tree.Branch("min_delta_phi", &min_delta_phi);

  reduced_tree.Branch("num_electrons", &num_electrons);
  reduced_tree.Branch("num_muons", &num_muons);
  reduced_tree.Branch("num_taus", &num_taus);
  reduced_tree.Branch("num_iso_tracks", &num_iso_tracks);

  reduced_tree.Branch("min_delta_R", &min_delta_r);
  reduced_tree.Branch("max_delta_R", &max_delta_r);

  reduced_tree.Branch("average_higgs_mass", &average_higgs_mass);
  reduced_tree.Branch("higgs_mass_difference", &higgs_mass_difference);
  
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
    passesPVCut = PassesPVCut() ? true : false;
    passesJet2PtCut = PassesJet2PtCut() ? true : false;
    passes2CSVTCut = Passes2CSVTCut() ? true : false;
    passesMETSig30Cut = PassesMETSig30Cut() ? true : false;
    passesMETCleaningCut = PassesMETCleaningCut() ? true : false;
    passesTriggerCut = PassesTriggerCut() ? true : false;
    passesNumJetsCut = PassesNumJetsCut() ? true : false;
    passesMinDeltaPhiCut = PassesMinDeltaPhiCut() ? true : false;
    passesLeptonVetoCut = PassesLeptonVetoCut() ? true : false;
    passesIsoTrackVetoCut = PassesIsoTrackVetoCut() ? true : false;
    passesDRCut = PassesDRCut() ? true : false;
    passesHiggsAvgMassCut= PassesHiggsAvgMassCut() ? true : false;
    passesHiggsMassDiffCut= PassesHiggsMassDiffCut() ? true : false;

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

    met_sig=pfmets_fullSignif;

    num_jets=GetNumGoodJets();
    num_b_tagged_jets=GetNumBTaggedJets();

    min_delta_phi=GetMinDeltaPhiMET();

    num_electrons=GetNumVetoElectrons();
    num_muons=GetNumVetoMuons();
    num_taus=GetNumVetoTaus();
    num_iso_tracks=NewGetNumIsoTracks();

    min_delta_r=GetMinDR();
    max_delta_r=GetMaxDR();

    const std::pair<double, double> higgs_masses(GetHiggsMasses());
    average_higgs_mass=0.5*(higgs_masses.first+higgs_masses.second);
    higgs_mass_difference=fabs(higgs_masses.first-higgs_masses.second);

    reduced_tree.Fill(); 
  }

  reduced_tree.Write();
  file.Close();
}
