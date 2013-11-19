#include "reduced_tree_maker.hpp"
#include <string>
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
  TTree reducedTree("reducedTree","reducedTree");
  bool passesPVCut(false), passesJet2PtCut(false), passes2CSVTCut(false), passesMETSig30Cut(false),
    passesMETCleaningCut(false), passesTriggerCut(false), passesNumJetsCut(false),
    passesMinDeltaPhiCut(false), passesLeptonVetoCut(false), passesIsoTrackVetoCut(false), passesDRCut(false);
  reducedTree.Branch("passesPVCut",&passesPVCut,"passesPVCut/O");
  reducedTree.Branch("passesJet2PtCut",&passesJet2PtCut,"passesJet2PtCut/O");
  reducedTree.Branch("passes2CSVTCut",&passes2CSVTCut,"passes2CSVTCut/O");
  reducedTree.Branch("passesMETSig30Cut",&passesMETSig30Cut,"passesMETSig30Cut/O");
  reducedTree.Branch("passesMETCleaningCut",&passesMETCleaningCut,"passesMETCleaningCut/O");
  reducedTree.Branch("passesTriggerCut",&passesTriggerCut,"passesTriggerCut/O");
  reducedTree.Branch("passesNumJetsCut",&passesNumJetsCut,"passesNumJetsCut/O");
  reducedTree.Branch("passesMinDeltaPhiCut",&passesMinDeltaPhiCut,"passesMinDeltaPhiCut/O");
  reducedTree.Branch("passesLeptonVetoCut",&passesLeptonVetoCut,"passesLeptonVetoCut/O");
  reducedTree.Branch("passesIsoTrackVetoCut",&passesIsoTrackVetoCut,"passesIsoTrackVetoCut/O");
  reducedTree.Branch("passesDRCut",&passesDRCut,"passesDRCut/O");
  
  // Now we're ready to go
  
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
    reducedTree.Fill(); 
    // As it stands, we can't fill any more branches beyond this point .
    // We might want a separate function for reduced tree making, outside of MakePlots.
  }
}
