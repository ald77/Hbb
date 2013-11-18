#include "cfa_skimmer.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "event_handler.hpp"
#include "event_number.hpp"
#include "pu_constants.hpp"
#include "lumi_reweighting_stand_alone.hpp"
#include "timer.hpp"
#include "in_json_2012.hpp"

CfASkimmer::CfASkimmer(const std::string& in_file_name,
                       const bool is_list,
                       const double weight_in):
  EventHandler(in_file_name, is_list, weight_in, false){
}

void CfASkimmer::Skim(const std::string& out_file_name,
                      const int chargino_mass,
                      const int LSP_mass){
  GetTotalEntries();

  TFile skimFile(out_file_name.c_str(), "RECREATE");
  TDirectory *cfA_dir(skimFile.mkdir("configurableAnalysis", "configurableAnalysis"));
  cfA_dir->cd();
  TTree *skimTreeA(chainA.CloneTree(0));
  chainA.CopyAddresses(skimTreeA);
  TTree *skimTreeB(chainB.CloneTree(0));
  chainB.CopyAddresses(skimTreeB);
  skimFile.cd();

  std::set<EventNumber> eventList;
  unsigned int startCount(0), PVCount(0), METCleaningCount(0), JSONCount(0), TriggerCount(0), numJetsCount(0), CSVTCount(0), jet2PtCount(0), minDeltaPhiCount(0), leptonVetoCount(0), isoTrackVetoCount(0), bTagCount(0), higgsCount(0), DRCount(0), METSig30Count(0), METSig50Count(0), METSig100Count(0), METSig150Count(0);
  double startCountWeighted(0), PVCountWeighted(0), METCleaningCountWeighted(0), JSONCountWeighted(0), TriggerCountWeighted(0), numJetsCountWeighted(0), CSVTCountWeighted(0), jet2PtCountWeighted(0), minDeltaPhiCountWeighted(0), leptonVetoCountWeighted(0), isoTrackVetoCountWeighted(0), bTagCountWeighted(0), higgsCountWeighted(0), DRCountWeighted(0), METSig30CountWeighted(0), METSig50CountWeighted(0), METSig100CountWeighted(0), METSig150CountWeighted(0);

  std::vector<float> dataDist(pu::RunsThrough203002, pu::RunsThrough203002+60);
  std::vector<float> MCDist(pu::Summer2012_S10, pu::Summer2012_S10+60);//QQQ this needs to change later for general pileup scenario
  reweight::LumiReWeighting lumiWeights(MCDist, dataDist);

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
    const bool isttbar(sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos);
    const double localWeight(GetPUWeight(lumiWeights)*scaleFactor*(isttbar?GetTopPtWeight():1.0)*GetSbinWeight());
    
    ++startCount;
    startCountWeighted+=localWeight;

    if(Passes2CSVTCut()
       && PassesPVCut()
       && PassesJet2PtCut()
       && PassesMETSig30Cut()
       && PassesTChiMassCut(chargino_mass, LSP_mass)){
      skimTreeA->Fill();
      skimTreeB->Fill();
    }

    if(!PassesPVCut()) continue;
    ++PVCount;
    PVCountWeighted+=localWeight;
    if(!PassesJet2PtCut()) continue;
    ++jet2PtCount;
    jet2PtCountWeighted+=localWeight;
    if(!Passes2CSVTCut()) continue;
    ++CSVTCount;
    CSVTCountWeighted+=localWeight;
    if(!PassesMETSig30Cut()) continue;
    ++METSig30Count;
    METSig30CountWeighted+=localWeight;
    if(!PassesMETCleaningCut()) continue;
    ++METCleaningCount;
    METCleaningCountWeighted+=localWeight;
    if(!PassesJSONCut()) continue;
    ++JSONCount;
    JSONCountWeighted+=localWeight;
    if(!PassesTriggerCut()) continue;
    ++TriggerCount;
    TriggerCountWeighted+=localWeight;
    if(!PassesNumJetsCut()) continue;
    ++numJetsCount;
    numJetsCountWeighted+=localWeight;
    if(!PassesMinDeltaPhiCut()) continue;
    ++minDeltaPhiCount;
    minDeltaPhiCountWeighted+=localWeight;
    if(!PassesLeptonVetoCut()) continue;
    ++leptonVetoCount;
    leptonVetoCountWeighted+=localWeight;
    if(!PassesIsoTrackVetoCut()) continue;
    ++isoTrackVetoCount;
    isoTrackVetoCountWeighted+=localWeight;
    if(!PassesBTaggingCut()) continue;
    ++bTagCount;
    bTagCountWeighted+=localWeight;
    if(!PassesHiggsMassCut()) continue;
    ++higgsCount;
    higgsCountWeighted+=localWeight;
    if(!PassesDRCut()) continue;
    ++DRCount;
    DRCountWeighted+=localWeight;
    if(!PassesMETSig50Cut()) continue;
    ++METSig50Count;
    METSig50CountWeighted+=localWeight;
    if(!PassesMETSig100Cut()) continue;
    ++METSig100Count;
    METSig100CountWeighted+=localWeight;
    if(!PassesMETSig150Cut()) continue;
    ++METSig150Count;
    METSig150CountWeighted+=localWeight;
  }

  TTree cutFlow("cutFlow", "cutFlow");
  cutFlow.Branch("startCount", &startCount);
  cutFlow.Branch("PVCount", &PVCount);
  cutFlow.Branch("jet2PtCount", &jet2PtCount);
  cutFlow.Branch("CSVTCount", &CSVTCount);
  cutFlow.Branch("METSig30Count", &METSig30Count);
  cutFlow.Branch("METCleaningCount", &METCleaningCount);
  cutFlow.Branch("TriggerCount", &TriggerCount);
  cutFlow.Branch("numJetsCount", &numJetsCount);
  cutFlow.Branch("minDeltaPhiCount", &minDeltaPhiCount);
  cutFlow.Branch("leptonVetoCount", &leptonVetoCount);
  cutFlow.Branch("isoTrackVetoCount", &isoTrackVetoCount);
  cutFlow.Branch("bTagCount", &bTagCount);
  cutFlow.Branch("higgsCount", &higgsCount);
  cutFlow.Branch("DRCount", &DRCount);
  cutFlow.Branch("METSig50Count", &METSig50Count);
  cutFlow.Branch("METSig100Count", &METSig100Count);
  cutFlow.Branch("METSig150Count", &METSig150Count);
  cutFlow.Branch("startCountWeighted", &startCountWeighted);
  cutFlow.Branch("PVCountWeighted", &PVCountWeighted);
  cutFlow.Branch("jet2PtCountWeighted", &jet2PtCountWeighted);
  cutFlow.Branch("CSVTCountWeighted", &CSVTCountWeighted);
  cutFlow.Branch("METSig30CountWeighted", &METSig30CountWeighted);
  cutFlow.Branch("METCleaningCountWeighted", &METCleaningCountWeighted);
  cutFlow.Branch("TriggerCountWeighted", &TriggerCountWeighted);
  cutFlow.Branch("numJetsCountWeighted", &numJetsCountWeighted);
  cutFlow.Branch("minDeltaPhiCountWeighted", &minDeltaPhiCountWeighted);
  cutFlow.Branch("leptonVetoCountWeighted", &leptonVetoCountWeighted);
  cutFlow.Branch("isoTrackVetoCountWeighted", &isoTrackVetoCountWeighted);
  cutFlow.Branch("bTagCountWeighted", &bTagCountWeighted);
  cutFlow.Branch("higgsCountWeighted", &higgsCountWeighted);
  cutFlow.Branch("DRCountWeighted", &DRCountWeighted);
  cutFlow.Branch("METSig50CountWeighted", &METSig50CountWeighted);
  cutFlow.Branch("METSig100CountWeighted", &METSig100CountWeighted);
  cutFlow.Branch("METSig150CountWeighted", &METSig150CountWeighted);
  cutFlow.Fill();

  std::cout << "startCount " << startCount << std::endl;
  std::cout << "PVCount " << PVCount << std::endl;
  std::cout << "jet2PtCount " << jet2PtCount << std::endl;
  std::cout << "CSVTCount " << CSVTCount << std::endl;
  std::cout << "METSig30Count " << METSig30Count << std::endl;
  std::cout << "METCleaningCount " << METCleaningCount << std::endl;
  std::cout << "TriggerCount " << TriggerCount << std::endl;
  std::cout << "numJetsCount " << numJetsCount << std::endl;
  std::cout << "minDeltaPhiCount " << minDeltaPhiCount << std::endl;
  std::cout << "leptonVetoCount " << leptonVetoCount << std::endl;
  std::cout << "isoTrackVetoCount " << isoTrackVetoCount << std::endl;
  std::cout << "bTagCount " << bTagCount << std::endl;
  std::cout << "higgsCount " << higgsCount << std::endl;
  std::cout << "DRCount " << DRCount << std::endl;
  std::cout << "METSig50Count " << METSig50Count << std::endl;
  std::cout << "METSig100Count " << METSig100Count << std::endl;
  std::cout << "METSig150Count " << METSig150Count << std::endl;
  std::cout << "startCountWeighted " << startCountWeighted << std::endl;
  std::cout << "PVCountWeighted " << PVCountWeighted << std::endl;
  std::cout << "jet2PtCountWeighted " << jet2PtCountWeighted << std::endl;
  std::cout << "CSVTCountWeighted " << CSVTCountWeighted << std::endl;
  std::cout << "METSig30CountWeighted " << METSig30CountWeighted << std::endl;
  std::cout << "METCleaningCountWeighted " << METCleaningCountWeighted << std::endl;
  std::cout << "TriggerCountWeighted " << TriggerCountWeighted << std::endl;
  std::cout << "numJetsCountWeighted " << numJetsCountWeighted << std::endl;
  std::cout << "minDeltaPhiCountWeighted " << minDeltaPhiCountWeighted << std::endl;
  std::cout << "leptonVetoCountWeighted " << leptonVetoCountWeighted << std::endl;
  std::cout << "isoTrackVetoCountWeighted " << isoTrackVetoCountWeighted << std::endl;
  std::cout << "bTagCountWeighted " << bTagCountWeighted << std::endl;
  std::cout << "higgsCountWeighted " << higgsCountWeighted << std::endl;
  std::cout << "DRCountWeighted " << DRCountWeighted << std::endl;
  std::cout << "METSig50CountWeighted " << METSig50CountWeighted << std::endl;
  std::cout << "METSig100CountWeighted " << METSig100CountWeighted << std::endl;
  std::cout << "METSig150CountWeighted " << METSig150CountWeighted << std::endl;

  TTree timeInfo("timeInfo", "timeInfo");
  time_t my_timer_t;
  struct tm *timer_s(0);
  time(&my_timer_t);
  timer_s=localtime(&my_timer_t);
  timeInfo.Branch("year", &timer_s->tm_year);
  timeInfo.Branch("month", &timer_s->tm_mon);
  timeInfo.Branch("day", &timer_s->tm_mday);
  timeInfo.Branch("hour", &timer_s->tm_hour);
  timeInfo.Branch("minute", &timer_s->tm_min);
  timeInfo.Branch("second", &timer_s->tm_sec);
  timeInfo.Fill();

  cfA_dir->cd();
  skimTreeA->Write();
  skimTreeB->Write();
  skimFile.cd();
  cutFlow.Write();
  timeInfo.Write();
  skimFile.Close();
}
