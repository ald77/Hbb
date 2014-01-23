#include "cfa_plotter.hpp"
#include <climits>
#include <set>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "timer.hpp"
#include "event_handler.hpp"

CfAPlotter::CfAPlotter(const std::string& in_file_name,
                       const bool is_list,
                       const double weight_in):
  EventHandler(in_file_name, is_list, weight_in, true){
}

void CfAPlotter::MakePlots(const std::string& out_file_name){
  TH1::SetDefaultSumw2(true);
  std::set<EventNumber> eventList;
  unsigned int startCount(0), PVCount(0), METCleaningCount(0), TriggerCount(0), numJetsCount(0), CSVTCount(0), jet2PtCount(0), minDeltaPhiCount(0), leptonVetoCount(0), isoTrackVetoCount(0), bTagCount(0), higgsCount(0), DRCount(0), METSig30Count(0), METSig50Count(0), METSig100Count(0), METSig150Count(0);
  double startCountWeighted(0), PVCountWeighted(0), METCleaningCountWeighted(0), TriggerCountWeighted(0), numJetsCountWeighted(0), CSVTCountWeighted(0), jet2PtCountWeighted(0), minDeltaPhiCountWeighted(0), leptonVetoCountWeighted(0), isoTrackVetoCountWeighted(0), bTagCountWeighted(0), higgsCountWeighted(0), DRCountWeighted(0), METSig30CountWeighted(0), METSig50CountWeighted(0), METSig100CountWeighted(0), METSig150CountWeighted(0);

  unsigned int ACount(0), BCount(0), C3bCount(0), D3bCount(0), C2bCount(0), D2bCount(0);
  double ACountWeighted(0.0), BCountWeighted(0.0), C3bCountWeighted(0.0), D3bCountWeighted(0.0), C2bCountWeighted(0.0), D2bCountWeighted(0.0);

  unsigned int ADRInvCount(0), BDRInvCount(0), C3bDRInvCount(0), D3bDRInvCount(0), C2bDRInvCount(0), D2bDRInvCount(0);
  double ADRInvCountWeighted(0.0), BDRInvCountWeighted(0.0), C3bDRInvCountWeighted(0.0), D3bDRInvCountWeighted(0.0), C2bDRInvCountWeighted(0.0), D2bDRInvCountWeighted(0.0);

  unsigned int ASLCount(0), BSLCount(0), C3bSLCount(0), D3bSLCount(0), C2bSLCount(0), D2bSLCount(0);
  double ASLCountWeighted(0.0), BSLCountWeighted(0.0), C3bSLCountWeighted(0.0), D3bSLCountWeighted(0.0), C2bSLCountWeighted(0.0), D2bSLCountWeighted(0.0);

  unsigned int ASLsbin1Count(0), BSLsbin1Count(0), C3bSLsbin1Count(0), D3bSLsbin1Count(0), C2bSLsbin1Count(0), D2bSLsbin1Count(0);
  double ASLsbin1CountWeighted(0.0), BSLsbin1CountWeighted(0.0), C3bSLsbin1CountWeighted(0.0), D3bSLsbin1CountWeighted(0.0), C2bSLsbin1CountWeighted(0.0), D2bSLsbin1CountWeighted(0.0);

  unsigned int ASLsbin2Count(0), BSLsbin2Count(0), C3bSLsbin2Count(0), D3bSLsbin2Count(0), C2bSLsbin2Count(0), D2bSLsbin2Count(0);
  double ASLsbin2CountWeighted(0.0), BSLsbin2CountWeighted(0.0), C3bSLsbin2CountWeighted(0.0), D3bSLsbin2CountWeighted(0.0), C2bSLsbin2CountWeighted(0.0), D2bSLsbin2CountWeighted(0.0);

  unsigned int ASLsbin3Count(0), BSLsbin3Count(0), C3bSLsbin3Count(0), D3bSLsbin3Count(0), C2bSLsbin3Count(0), D2bSLsbin3Count(0);
  double ASLsbin3CountWeighted(0.0), BSLsbin3CountWeighted(0.0), C3bSLsbin3CountWeighted(0.0), D3bSLsbin3CountWeighted(0.0), C2bSLsbin3CountWeighted(0.0), D2bSLsbin3CountWeighted(0.0);

  unsigned int ASLsbin4Count(0), BSLsbin4Count(0), C3bSLsbin4Count(0), D3bSLsbin4Count(0), C2bSLsbin4Count(0), D2bSLsbin4Count(0);
  double ASLsbin4CountWeighted(0.0), BSLsbin4CountWeighted(0.0), C3bSLsbin4CountWeighted(0.0), D3bSLsbin4CountWeighted(0.0), C2bSLsbin4CountWeighted(0.0), D2bSLsbin4CountWeighted(0.0);

  unsigned int ADRInvsbin1Count(0), BDRInvsbin1Count(0), C3bDRInvsbin1Count(0), D3bDRInvsbin1Count(0), C2bDRInvsbin1Count(0), D2bDRInvsbin1Count(0);
  double ADRInvsbin1CountWeighted(0.0), BDRInvsbin1CountWeighted(0.0), C3bDRInvsbin1CountWeighted(0.0), D3bDRInvsbin1CountWeighted(0.0), C2bDRInvsbin1CountWeighted(0.0), D2bDRInvsbin1CountWeighted(0.0);

  unsigned int ADRInvsbin2Count(0), BDRInvsbin2Count(0), C3bDRInvsbin2Count(0), D3bDRInvsbin2Count(0), C2bDRInvsbin2Count(0), D2bDRInvsbin2Count(0);
  double ADRInvsbin2CountWeighted(0.0), BDRInvsbin2CountWeighted(0.0), C3bDRInvsbin2CountWeighted(0.0), D3bDRInvsbin2CountWeighted(0.0), C2bDRInvsbin2CountWeighted(0.0), D2bDRInvsbin2CountWeighted(0.0);

  unsigned int ADRInvsbin3Count(0), BDRInvsbin3Count(0), C3bDRInvsbin3Count(0), D3bDRInvsbin3Count(0), C2bDRInvsbin3Count(0), D2bDRInvsbin3Count(0);
  double ADRInvsbin3CountWeighted(0.0), BDRInvsbin3CountWeighted(0.0), C3bDRInvsbin3CountWeighted(0.0), D3bDRInvsbin3CountWeighted(0.0), C2bDRInvsbin3CountWeighted(0.0), D2bDRInvsbin3CountWeighted(0.0);

  unsigned int ADRInvsbin4Count(0), BDRInvsbin4Count(0), C3bDRInvsbin4Count(0), D3bDRInvsbin4Count(0), C2bDRInvsbin4Count(0), D2bDRInvsbin4Count(0);
  double ADRInvsbin4CountWeighted(0.0), BDRInvsbin4CountWeighted(0.0), C3bDRInvsbin4CountWeighted(0.0), D3bDRInvsbin4CountWeighted(0.0), C2bDRInvsbin4CountWeighted(0.0), D2bDRInvsbin4CountWeighted(0.0);

  unsigned int Asbin1Count(0), Bsbin1Count(0), C3bsbin1Count(0), D3bsbin1Count(0), C2bsbin1Count(0), D2bsbin1Count(0);
  double Asbin1CountWeighted(0.0), Bsbin1CountWeighted(0.0), C3bsbin1CountWeighted(0.0), D3bsbin1CountWeighted(0.0), C2bsbin1CountWeighted(0.0), D2bsbin1CountWeighted(0.0);

  unsigned int Asbin2Count(0), Bsbin2Count(0), C3bsbin2Count(0), D3bsbin2Count(0), C2bsbin2Count(0), D2bsbin2Count(0);
  double Asbin2CountWeighted(0.0), Bsbin2CountWeighted(0.0), C3bsbin2CountWeighted(0.0), D3bsbin2CountWeighted(0.0), C2bsbin2CountWeighted(0.0), D2bsbin2CountWeighted(0.0);

  unsigned int Asbin3Count(0), Bsbin3Count(0), C3bsbin3Count(0), D3bsbin3Count(0), C2bsbin3Count(0), D2bsbin3Count(0);
  double Asbin3CountWeighted(0.0), Bsbin3CountWeighted(0.0), C3bsbin3CountWeighted(0.0), D3bsbin3CountWeighted(0.0), C2bsbin3CountWeighted(0.0), D2bsbin3CountWeighted(0.0);

  unsigned int Asbin4Count(0), Bsbin4Count(0), C3bsbin4Count(0), D3bsbin4Count(0), C2bsbin4Count(0), D2bsbin4Count(0);
  double Asbin4CountWeighted(0.0), Bsbin4CountWeighted(0.0), C3bsbin4CountWeighted(0.0), D3bsbin4CountWeighted(0.0), C2bsbin4CountWeighted(0.0), D2bsbin4CountWeighted(0.0);

  std::vector<float> dataDist(pu::RunsThrough203002, pu::RunsThrough203002+60);
  std::vector<float> MCDist(pu::Summer2012_S10, pu::Summer2012_S10+60);//QQQ this needs to change later for general pileup scenario
  reweight::LumiReWeighting lumiWeights(MCDist, dataDist);

  const double metSigABCDBinEdges[6]={30.0,50.0,100.0,150.0,200.0,300.0};
  TH1D xx_metSig_SigOverSB_2b("xx_metSig_SigOverSB_4b","Mass Signal Window/Mass SB;MET Significance;Mass Signal Window/Mass SB", 5, metSigABCDBinEdges);
  TH1D xx_metSig_SigOverSB_3b("xx_metSig_SigOverSB_3b","Mass Signal Window/Mass SB;MET Significance;Mass Signal Window/Mass SB", 5, metSigABCDBinEdges);
  TH1D xx_metSig_SigOverSB_4b("xx_metSig_SigOverSB_2b","Mass Signal Window/Mass SB;MET Significance;Mass Signal Window/Mass SB", 5, metSigABCDBinEdges);

  TFile file(out_file_name.c_str(), "recreate");

  TCanvas xx_metSig_SigOverSB("xx_metSig_SigOverSB","xx_metSig_SigOverSB");

  TH1D xx_nm1_thirdBTag("xx_nm1_thirdBTag", "N-1 Third B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D xx_nm1_fourthBTag("xx_nm1_fourthBTag", "N-1 Fourth B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D xx_nm1_avgHiggsMass("xx_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D xx_nm1_higgsMassDiff("xx_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D xx_nm1_maxDeltaR("xx_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D xx_nm1_metSig("xx_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D xx_nm1_met("xx_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25 GeV", 40, 0.0, 1000.0);

  TH2D xx_MET_vs_minDeltaPhi("xx_MET_vs_minDeltaPhi","MET vs. Min. #Delta#phi;Min. #Delta#Phi;MET [GeV]",15,0.0,3.0,15,0.0,450.0);
  TH2D xx_METSig_vs_minDeltaPhi("xx_METSig_vs_minDeltaPhi","MET Significance vs. Min. #Delta#phi;Min. #Delta#Phi;MET Significance",15,0.0,3.0,15,0.0,300.0);

  TH2D xx_METOverSqrtHT_vs_METSig("xx_METOverSqrtHT_vs_METSig", "MET/#sqrt{HT} vs. MET Sig.;MET Significance; MET/#sqrt{HT} [#sqrt{GeV}]", 40, 0.0, 200.0, 24, 0.0, 12.0);
  TH1D xx_METOverSqrtHT_Over_METSig("xx_METOverSqrtHT_Over_METSig","(MET/#sqrt{HT})/MET Sig;(MET/#sqrt{HT})/MET Sig [#sqrt{GeV}];Events/19.4 fb^{-1}/0.01",30,0.0,0.3);
  TH1D xx_MET("xx_MET","MET;MET [GeV];Events/19.4 fb^{-1}/10 GeV",40,0.0,400.0);
  TH1D xx_ABCD("xx_ABCD", "ABCD;;Events/19.4 fb^{-1}",6,-0.5,5.5);
  xx_ABCD.GetXaxis()->SetBinLabel(1,"A");
  xx_ABCD.GetXaxis()->SetBinLabel(2,"B");
  xx_ABCD.GetXaxis()->SetBinLabel(3,"C(3b)");
  xx_ABCD.GetXaxis()->SetBinLabel(4,"D(3b)");
  xx_ABCD.GetXaxis()->SetBinLabel(5,"C(2b)");
  xx_ABCD.GetXaxis()->SetBinLabel(6,"D(2b)");

  TH1D xx_metSig_A("xx_metSig_A","MET Significance (Signal Mass Window, 4b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_B("xx_metSig_B","MET Significance (Mass Sideband, 4b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_C3b("xx_metSig_C3b","MET Significance (Signal Mass Window, 3b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_D3b("xx_metSig_D3b","MET Significance (Mass Sideband, 3b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_C2b("xx_metSig_C2b","MET Significance (Signal Mass Window, 2b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_D2b("xx_metSig_D2b","MET Significance (Mass Sideband, 2b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH2D xx_METToMETSigRatio_vs_numLowPtPfCands("xx_METToMETSigRatio_vs_numLowPtPfCands","MET/MET Sig. vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,80,0.0,8.0);
  TH2D xx_MET_vs_numLowPtPfCands("xx_MET_vs_numLowPtPfCands","MET vs. Low PT PF Cands.;Low PT PF Cands.;MET [GeV]",16,-0.5,15.5,40,0.0,600.0);
  TH2D xx_METSig_vs_numLowPtPfCands("xx_METSig_vs_numLowPtPfCands","MET Sig vs. Low PT PF Cands.;Low PT PF Cands.;MET Sig.",16,-0.5,15.5,40,0.0,400.0);
  TH2D xx_METToMETSigRatio_vs_fracMETFromLowPt("xx_METToMETSigRatio_vs_fracMETFromLowPt","MET/MET Sig. vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET/MET Sig. [GeV]",40,0.0,1.0,80,0.0,8.0);
  TH2D xx_MET_vs_fracMETFromLowPt("xx_MET_vs_fracMETFromLowPt","MET vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET [GeV]",40,0.0,1.0,40,0.0,600.0);
  TH2D xx_METSig_vs_fracMETFromLowPt("xx_METSig_vs_fracMETFromLowPt","MET Sig vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET Sig.",40,0.0,1.0,40,0.0,400.0);

  //Two b-jet control region
  TH1D twoB_nm1_thirdBTag("twoB_nm1_thirdBTag", "N-1 Third B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D twoB_nm1_fourthBTag("twoB_nm1_fourthBTag", "N-1 Fourth B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D twoB_nm1_avgHiggsMass("twoB_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D twoB_nm1_higgsMassDiff("twoB_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D twoB_nm1_maxDeltaR("twoB_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D twoB_nm1_metSig("twoB_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D twoB_nm1_met("twoB_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25 GeV", 40, 0.0, 1000.0);
  TH2D twoB_METOverSqrtHT_vs_METSig("twoB_METOverSqrtHT_vs_METSig", "MET/#sqrt{HT} vs. MET Sig.;MET Significance; MET/#sqrt{HT} [#sqrt{GeV}]", 40, 0.0, 200.0, 24, 0.0, 12.0);
  TH1D twoB_METOverSqrtHT_Over_METSig("twoB_METOverSqrtHT_Over_METSig","(MET/#sqrt{HT})/MET Sig;(MET/#sqrt{HT})/MET Sig [#sqrt{GeV}];Events/19.4 fb^{-1}/0.01",30,0.0,0.3);
  TH1D twoB_MET("twoB_MET","MET;MET [GeV];Events/19.4 fb^{-1}/10 GeV",40,0.0,400.0);

  TH2D twoB_MET_vs_minDeltaPhi("twoB_MET_vs_minDeltaPhi","MET vs. Min. #Delta#phi;Min. #Delta#Phi;MET [GeV]",15,0.0,3.0,15,0.0,450.0);
  TH2D twoB_METSig_vs_minDeltaPhi("twoB_METSig_vs_minDeltaPhi","MET Significance vs. Min. #Delta#phi;Min. #Delta#Phi;MET Significance",15,0.0,3.0,15,0.0,300.0);
  TH2D twoB_METToMETSigRatio_vs_numLowPtPfCands("twoB_METToMETSigRatio_vs_numLowPtPfCands","MET/MET Sig. vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,80,0.0,8.0);
  TH2D twoB_MET_vs_numLowPtPfCands("twoB_MET_vs_numLowPtPfCands","MET vs. Low PT PF Cands.;Low PT PF Cands.;MET [GeV]",16,-0.5,15.5,40,0.0,600.0);
  TH2D twoB_METSig_vs_numLowPtPfCands("twoB_METSig_vs_numLowPtPfCands","MET Sig vs. Low PT PF Cands.;Low PT PF Cands.;MET Sig.",16,-0.5,15.5,40,0.0,400.0);
  TH2D twoB_METToMETSigRatio_vs_fracMETFromLowPt("twoB_METToMETSigRatio_vs_fracMETFromLowPt","MET/MET Sig. vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET/MET Sig. [GeV]",40,0.0,1.0,80,0.0,8.0);
  TH2D twoB_MET_vs_fracMETFromLowPt("twoB_MET_vs_fracMETFromLowPt","MET vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET [GeV]",40,0.0,1.0,40,0.0,600.0);
  TH2D twoB_METSig_vs_fracMETFromLowPt("twoB_METSig_vs_fracMETFromLowPt","MET Sig vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET Sig.",40,0.0,1.0,40,0.0,400.0);

  //Light Higgs control region
  TH1D lightHiggs_nm1_thirdBTag("lightHiggs_nm1_thirdBTag", "N-1 Third B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D lightHiggs_nm1_fourthBTag("lightHiggs_nm1_fourthBTag", "N-1 Fourth B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D lightHiggs_nm1_avgHiggsMass("lightHiggs_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D lightHiggs_nm1_higgsMassDiff("lightHiggs_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D lightHiggs_nm1_maxDeltaR("lightHiggs_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D lightHiggs_nm1_metSig("lightHiggs_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D lightHiggs_nm1_met("lightHiggs_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH2D lightHiggs_METOverSqrtHT_vs_METSig("lightHiggs_METOverSqrtHT_vs_METSig", "MET/#sqrt{HT} vs. MET Sig.;MET Significance; MET/#sqrt{HT} [#sqrt{GeV}]", 40, 0.0, 200.0, 24, 0.0, 12.0);
  TH1D lightHiggs_METOverSqrtHT_Over_METSig("lightHiggs_METOverSqrtHT_Over_METSig","(MET/#sqrt{HT})/MET Sig;(MET/#sqrt{HT})/MET Sig [#sqrt{GeV}];Events/19.4 fb^{-1}/0.01",30,0.0,0.3);
  TH1D lightHiggs_MET("lightHiggs_MET","MET;MET [GeV];Events/19.4 fb^{-1}/25 GeV",40,0.0,400.0);

  TH2D lightHiggs_MET_vs_minDeltaPhi("lightHiggs_MET_vs_minDeltaPhi","MET vs. Min. #Delta#phi;Min. #Delta#Phi;MET [GeV]",15,0.0,3.0,15,0.0,450.0);
  TH2D lightHiggs_METSig_vs_minDeltaPhi("lightHiggs_METSig_vs_minDeltaPhi","MET Significance vs. Min. #Delta#phi;Min. #Delta#Phi;MET Significance",15,0.0,3.0,15,0.0,300.0);
  TH2D lightHiggs_METToMETSigRatio_vs_numLowPtPfCands("lightHiggs_METToMETSigRatio_vs_numLowPtPfCands","MET/MET Sig. vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,80,0.0,8.0);
  TH2D lightHiggs_MET_vs_numLowPtPfCands("lightHiggs_MET_vs_numLowPtPfCands","MET vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,40,0.0,600.0);
  TH2D lightHiggs_METSig_vs_numLowPtPfCands("lightHiggs_METSig_vs_numLowPtPfCands","MET Sig vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,40,0.0,400.0);
  TH2D lightHiggs_METToMETSigRatio_vs_fracMETFromLowPt("lightHiggs_METToMETSigRatio_vs_fracMETFromLowPt","MET/MET Sig. vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET/MET Sig. [GeV]",40,0.0,1.0,80,0.0,8.0);
  TH2D lightHiggs_MET_vs_fracMETFromLowPt("lightHiggs_MET_vs_fracMETFromLowPt","MET vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET [GeV]",40,0.0,1.0,40,0.0,600.0);
  TH2D lightHiggs_METSig_vs_fracMETFromLowPt("lightHiggs_METSig_vs_fracMETFromLowPt","MET Sig vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET Sig.",40,0.0,1.0,40,0.0,400.0);

  TH1D sorb_MET("sorb_MET","S/#sqrt{B};MET [GeV];S/#sqrt{B}",1024,0.0,350.0);
  TH1D sorb_METSig("sorb_METSig","S/#sqrt{B};MET Significance;S/#sqrt{B}",1024,0.0,150.0);
  TH1D sorb_METSig2012("sorb_METSig2012","S/#sqrt{B};MET Significance 2012;S/#sqrt{B}",1024,0.0,150.0);
  TH1D sorb_METOverSqrtHT("sorb_METOverSqrtHT","S/#sqrt{B};MET/#sqrt{HT} [#sqrt{GeV}];S/#sqrt{B}",1024,0.0,15.0);

  TH1D npv_before("npv_before", "Num. Primary Vertices Before Pileup Reweighting;Vertices;Events/19.4 fb^{-1}", 61, -0.5, 60.5);
  TH1D npv_after("npv_after", "Num. Primary Vertices After Pileup Reweighting;Vertices;Events/19.4 fb^{-1}", 61, -0.5, 60.5);

  TH1D topPt_before("topPt_before", "Top p_{T} Before Top p_{T} Reweighting;top p_{T} [GeV];Events/5 GeV/19.4 fb^{-1}", 50, 0.0, 500.0);
  TH1D topPt_after("topPt_after", "Top p_{T} After Top p_{T} Reweighting;top p_{T} [GeV];Events/5 GeV/19.4 fb^{-1}", 50, 0.0, 500.0);

  TH1D xx_origin1("xx_origin1", "Origin of Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D xx_origin2("xx_origin2", "Origin of Second Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D xx_origin3("xx_origin3", "Origin of Third Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D xx_origin4("xx_origin4", "Origin of Fourth Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  CfAPlots::FixBinLabels(xx_origin1);
  CfAPlots::FixBinLabels(xx_origin2);
  CfAPlots::FixBinLabels(xx_origin3);
  CfAPlots::FixBinLabels(xx_origin4);

  TH1D xx_beta1("xx_beta1", "beta for Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D xx_beta2("xx_beta2", "beta for Second Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D xx_beta3("xx_beta3", "beta for Third Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D xx_beta4("xx_beta4", "beta for Fourth Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);

  TH1D yy_origin1("yy_origin1", "Origin of Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D yy_origin2("yy_origin2", "Origin of Second Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D yy_origin3("yy_origin3", "Origin of Third Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D yy_origin4("yy_origin4", "Origin of Fourth Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  CfAPlots::FixBinLabels(yy_origin1);
  CfAPlots::FixBinLabels(yy_origin2);
  CfAPlots::FixBinLabels(yy_origin3);
  CfAPlots::FixBinLabels(yy_origin4);

  TH1D yy_beta1("yy_beta1", "beta for Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D yy_beta2("yy_beta2", "beta for Second Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D yy_beta3("yy_beta3", "beta for Third Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D yy_beta4("yy_beta4", "beta for Fourth Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);

  TH1D xx_Sbins_4bSig("xx_Sbins_4bSig", "MET Significance Counts (4b, mass signal window);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_4bSB("yy_Sbins_4bSB", "MET Significance Counts (4b, mass sideband);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_3bSig("yy_Sbins_3bSig", "MET Significance Counts (3b, mass signal window);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_3bSB("yy_Sbins_3bSB", "MET Significance Counts (3b, mass sideband);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_2bSig("yy_Sbins_2bSig", "MET Significance Counts (2b, mass signal window);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_2bSB("yy_Sbins_2bSB", "MET Significance Counts (2b, mass sideband);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  CfAPlots::FixSbinLabels(xx_Sbins_4bSig);
  CfAPlots::FixSbinLabels(yy_Sbins_4bSB);
  CfAPlots::FixSbinLabels(yy_Sbins_3bSig);
  CfAPlots::FixSbinLabels(yy_Sbins_3bSB);
  CfAPlots::FixSbinLabels(yy_Sbins_2bSig);
  CfAPlots::FixSbinLabels(yy_Sbins_2bSB);

  TH1D sl_Sbins_4bSig("sl_Sbins_4bSig", "MET Significance Counts (4b, mass signal window, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_4bSB("sl_Sbins_4bSB", "MET Significance Counts (4b, mass sideband, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_3bSig("sl_Sbins_3bSig", "MET Significance Counts (3b, mass signal window, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_3bSB("sl_Sbins_3bSB", "MET Significance Counts (3b, mass sideband, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_2bSig("sl_Sbins_2bSig", "MET Significance Counts (2b, mass signal window, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_2bSB("sl_Sbins_2bSB", "MET Significance Counts (2b, mass sideband, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  CfAPlots::FixSbinLabels(sl_Sbins_4bSig);
  CfAPlots::FixSbinLabels(sl_Sbins_4bSB);
  CfAPlots::FixSbinLabels(sl_Sbins_3bSig);
  CfAPlots::FixSbinLabels(sl_Sbins_3bSB);
  CfAPlots::FixSbinLabels(sl_Sbins_2bSig);
  CfAPlots::FixSbinLabels(sl_Sbins_2bSB);

  TH1D sbin1_4b_sb_nm1_thirdCSV("sbin1_4b_sb_nm1_thirdCSV", "N-1 Third CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D sbin1_4b_sb_nm1_fourthCSV("sbin1_4b_sb_nm1_fourthCSV", "N-1 Fourth CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D sbin1_4b_sb_nm1_avgHiggsMass("sbin1_4b_sb_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D sbin1_4b_sb_nm1_higgsMassDiff("sbin1_4b_sb_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D sbin1_4b_sb_nm1_maxDeltaR("sbin1_4b_sb_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D sbin1_4b_sb_nm1_metSig("sbin1_4b_sb_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D sbin1_4b_sb_nm1_met("sbin1_4b_sb_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH1D sbin1_4b_sb_nm1_numJets("sbin1_4b_sb_nm1_numJets", "N-1 Number of Jets;Number of Jets;Events/19.4 fb^{-1}", 16, -0.5, 15.5);
  TH1D sbin1_4b_sb_nm1_minDeltaPhi("sbin1_4b_sb_nm1_minDeltaPhi", "N-1 min. #Delta#phi;min. #Delta#phi;Events/19.4 fb^{-1}", 30, 0.0, 4.0*atan(1.0));
  
  TH1D single_lepton_nm1_thirdCSV("single_lepton_nm1_thirdCSV", "N-1 Third CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D single_lepton_nm1_fourthCSV("single_lepton_nm1_fourthCSV", "N-1 Fourth CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D single_lepton_nm1_avgHiggsMass("single_lepton_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D single_lepton_nm1_higgsMassDiff("single_lepton_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D single_lepton_nm1_maxDeltaR("single_lepton_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D single_lepton_nm1_metSig("single_lepton_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D single_lepton_nm1_met("single_lepton_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH1D single_lepton_nm1_numJets("single_lepton_nm1_numJets", "N-1 Number of Jets;Number of Jets;Events/19.4 fb^{-1}", 16, -0.5, 15.5);
  TH1D single_lepton_nm1_minDeltaPhi("single_lepton_nm1_minDeltaPhi", "N-1 min. #Delta#phi;min. #Delta#phi;Events/19.4 fb^{-1}", 30, 0.0, 4.0*atan(1.0));

  TH1D sig_sb_nm1_thirdCSV("sig_sb_nm1_thirdCSV", "N-1 Third CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D sig_sb_nm1_fourthCSV("sig_sb_nm1_fourthCSV", "N-1 Fourth CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D sig_sb_nm1_avgHiggsMass("sig_sb_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D sig_sb_nm1_higgsMassDiff("sig_sb_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D sig_sb_nm1_maxDeltaR("sig_sb_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D sig_sb_nm1_metSig("sig_sb_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D sig_sb_nm1_met("sig_sb_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH1D sig_sb_nm1_numJets("sig_sb_nm1_numJets", "N-1 Number of Jets;Number of Jets;Events/19.4 fb^{-1}", 16, -0.5, 15.5);
  TH1D sig_sb_nm1_minDeltaPhi("sig_sb_nm1_minDeltaPhi", "N-1 min. #Delta#phi;min. #Delta#phi;Events/19.4 fb^{-1}", 30, 0.0, 4.0*atan(1.0));

  TH1D minimax_mbb("minimax_mbb", "Minimax m_{bb};Minimax m_{bb};Events/19.4 fb^{-1}/10", 50, 0.0, 500.0);
  TH1D maximin_mbb("maximin_mbb", "Maximin m_{bb};Maximin m_{bb};Events/19.4 fb^{-1}/10", 50, 0.0, 500.0);

  // jet flavor comparisons
  TH1D jet_parton_Id_CSV1("jet_parton_Id_CSV1",";jets_AK5PF_partonId;Jets",5,0.5,5.5);
  TH1D jet_parton_Id_CSV2("jet_parton_Id_CSV2",";jets_AK5PF_partonId;Jets",5,0.5,5.5);
  TH1D jet_parton_Id_CSV3("jet_parton_Id_CSV3",";jets_AK5PF_partonId;Jets",5,0.5,5.5);
  TH1D jet_parton_Id_CSV4("jet_parton_Id_CSV4",";jets_AK5PF_partonId;Jets",5,0.5,5.5);
  const int n_Ids(5);
  char partonFlavour_names[n_Ids][32] = {"X","lf","c","b","g"};
  for (int i=1;i<=n_Ids;i++) {
    jet_parton_Id_CSV1.GetXaxis()->SetBinLabel(i,partonFlavour_names[i-1]);
    jet_parton_Id_CSV2.GetXaxis()->SetBinLabel(i,partonFlavour_names[i-1]);
    jet_parton_Id_CSV3.GetXaxis()->SetBinLabel(i,partonFlavour_names[i-1]);
    jet_parton_Id_CSV4.GetXaxis()->SetBinLabel(i,partonFlavour_names[i-1]);
  }
  

  // Define reduced tree and all variables to save
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
  for(int i(0); i<GetTotalEntries() && i<50000; ++i){
    if(i%1000==0 && i!=0){
      timer.PrintRemainingTime();
    }
    timer.Iterate();
    GetEntry(i);

    const bool isRealData(sampleName.find("Run2012")!=std::string::npos);
    const bool isttbar(sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos);
    const double localWeight((isRealData?1.0:GetPUWeight(lumiWeights))*scaleFactor*(true && isttbar?GetTopPtWeight():1.0)*GetSbinWeight()*(HasGluonSplitting()?1.0:1.0));
    const double nonpileupweight(scaleFactor*(true && isttbar?GetTopPtWeight():1.0)*(HasGluonSplitting()?1.0:1.0)*GetSbinWeight());
    const double notopweight((isRealData?1.0:GetPUWeight(lumiWeights))*scaleFactor*GetSbinWeight());

    if(!PassesJSONCut()) continue;
    if(!PassesTChiMassCut(300,1)
       && sampleName.find("SMS-TChiZH")!=std::string::npos){
      continue;
    }

    std::pair<std::set<EventNumber>::iterator, bool> returnVal(eventList.insert(EventNumber(run, event, lumiblock)));
    if(!returnVal.second) continue;
    ++startCount;
    startCountWeighted+=localWeight;

    if(PassesPVCut() && PassesJet2PtCut() && Passes2CSVTCut() && PassesMETSig30Cut() && PassesMETCleaningCut() && PassesTriggerCut() && PassesNumJetsCut() && PassesMinDeltaPhiCut() && PassesLeptonVetoCut() && PassesIsoTrackVetoCut() && PassesDRCut()){
      double npv(-1.0);
      for(unsigned int bc(0); bc<PU_bunchCrossing->size(); ++bc){
        if(PU_bunchCrossing->at(bc)==0){
          npv = PU_TrueNumInteractions->at(bc);
        }
      }
      npv_before.Fill(npv, nonpileupweight);
      npv_after.Fill(npv, localWeight);
      double topPt(-1.0);
      for ( unsigned int mc=0; mc< mc_doc_id->size(); mc++ ) {
        if ( mc_doc_id->at(mc)==6 ) { topPt = mc_doc_pt->at(mc); break; }
      }
      topPt_before.Fill(topPt, notopweight);
      topPt_after.Fill(topPt, localWeight);

      minimax_mbb.Fill(GetMinimaxMbb(), localWeight);
      maximin_mbb.Fill(GetMaximinMbb(), localWeight);
    }

    if(!betaUpToDate) GetBeta();
    if(PassesRegionACut()){
      const std::vector<std::pair<int,int> > bo(GetBOrigins());
      xx_origin1.Fill(CfAPlots::GetType(bo.at(0)),localWeight);
      xx_origin2.Fill(CfAPlots::GetType(bo.at(1)),localWeight);
      xx_origin3.Fill(CfAPlots::GetType(bo.at(2)),localWeight);
      xx_origin4.Fill(CfAPlots::GetType(bo.at(3)),localWeight);
      std::vector<std::pair<double, double> > CSV_beta(0);
      for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
        if(isGoodJet(jet, true, 20.0, 2.4, false)){
          CSV_beta.push_back(std::pair<double,double>(beta.at(jet),
                                                      jets_AK5PF_btag_secVertexCombined->at(jet)));
        }
      }
      std::sort(CSV_beta.begin(), CSV_beta.end(), std::greater<std::pair<double, double> >());
      xx_beta1.Fill(CSV_beta.at(0).first,localWeight);
      xx_beta2.Fill(CSV_beta.at(1).first,localWeight);
      xx_beta3.Fill(CSV_beta.at(2).first,localWeight);
      xx_beta4.Fill(CSV_beta.at(3).first,localWeight);
    }

    if(PassesRegionBCut()){
      const std::vector<std::pair<int,int> > bo(GetBOrigins());
      yy_origin1.Fill(CfAPlots::GetType(bo.at(0)),localWeight);
      yy_origin2.Fill(CfAPlots::GetType(bo.at(1)),localWeight);
      yy_origin3.Fill(CfAPlots::GetType(bo.at(2)),localWeight);
      yy_origin4.Fill(CfAPlots::GetType(bo.at(3)),localWeight);
      std::vector<std::pair<double, double> > CSV_beta(0);
      for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
        if(isGoodJet(jet, true, 20.0, 2.4, false)){
          CSV_beta.push_back(std::pair<double,double>(beta.at(jet),
                                                      jets_AK5PF_btag_secVertexCombined->at(jet)));
        }
      }
      std::sort(CSV_beta.begin(), CSV_beta.end(), std::greater<std::pair<double, double> >());
      yy_beta1.Fill(CSV_beta.at(0).first,localWeight);
      yy_beta2.Fill(CSV_beta.at(1).first,localWeight);
      yy_beta3.Fill(CSV_beta.at(2).first,localWeight);
      yy_beta4.Fill(CSV_beta.at(3).first,localWeight);
    }

    double SbinToFill(0.0);
    if(pfmets_fullSignif>30.0 && pfmets_fullSignif<50.0) SbinToFill=1.0;
    if(pfmets_fullSignif>50.0 && pfmets_fullSignif<100.0) SbinToFill=2.0;
    if(pfmets_fullSignif>100.0 && pfmets_fullSignif<150.0) SbinToFill=3.0;
    if(pfmets_fullSignif>150.0) SbinToFill=4.0;

    if(PassesRegionACut()) xx_Sbins_4bSig.Fill(SbinToFill, localWeight);
    if(PassesRegionBCut()) yy_Sbins_4bSB.Fill(SbinToFill, localWeight);
    if(PassesRegionC3bCut()) yy_Sbins_3bSig.Fill(SbinToFill, localWeight);
    if(PassesRegionD3bCut()) yy_Sbins_3bSB.Fill(SbinToFill, localWeight);
    if(PassesRegionC2bCut()) yy_Sbins_2bSig.Fill(SbinToFill, localWeight);
    if(PassesRegionD2bCut()) yy_Sbins_2bSB.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionACut()) sl_Sbins_4bSig.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionBCut()) sl_Sbins_4bSB.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionC3bCut()) sl_Sbins_3bSig.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionD3bCut()) sl_Sbins_3bSB.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionC2bCut()) sl_Sbins_2bSig.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionD2bCut()) sl_Sbins_2bSB.Fill(SbinToFill, localWeight);

    const uint_least32_t failure_code(GetCutFailCode());
    uint_least32_t masked_fc(failure_code & ~kMETSig100 & ~kMETSig150);
    if((masked_fc & ~kNumJets) == kMETSig50) sbin1_4b_sb_nm1_numJets.Fill(GetNumGoodJets(),localWeight);
    if((masked_fc & ~kMinDeltaPhi) == kMETSig50) sbin1_4b_sb_nm1_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),localWeight);
    if((masked_fc & ~k3rdBTag & ~k4thBTag) == kMETSig50) sbin1_4b_sb_nm1_thirdCSV.Fill(GetHighestCSV(3),localWeight);
    if((masked_fc & ~k4thBTag) == kMETSig50) sbin1_4b_sb_nm1_fourthCSV.Fill(GetHighestCSV(4),localWeight);
    const std::pair<double,double> higgsMasses(GetHiggsMasses());
    if((masked_fc & ~kHiggsMassDiff) == kMETSig50) sbin1_4b_sb_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
    if((masked_fc & ~kHiggsAvgMass) == kMETSig50) sbin1_4b_sb_nm1_higgsMassDiff.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
    if((masked_fc & ~kDeltaR) == kMETSig50) sbin1_4b_sb_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
    if(masked_fc == kMETSig50) sbin1_4b_sb_nm1_met.Fill(pfmets_et->at(0),localWeight);
    if(masked_fc == kMETSig50) sbin1_4b_sb_nm1_metSig.Fill(pfmets_fullSignif,localWeight);

    masked_fc=(failure_code & ~kMETSig50 & ~kMETSig100 & ~kMETSig150);
    if(PassesSingleLeptonCut()){
      if((masked_fc & ~kNumJets) == kGood) single_lepton_nm1_numJets.Fill(GetNumGoodJets(),localWeight);
      if((masked_fc & ~kMinDeltaPhi) == kGood) single_lepton_nm1_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),localWeight);
      if((masked_fc & ~k3rdBTag & ~k4thBTag) == kGood) single_lepton_nm1_thirdCSV.Fill(GetHighestCSV(3),localWeight);
      if((masked_fc & ~k4thBTag) == kGood) single_lepton_nm1_fourthCSV.Fill(GetHighestCSV(4),localWeight);
      if((masked_fc & ~kHiggsMassDiff) == kGood) single_lepton_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
      if((masked_fc & ~kHiggsAvgMass) == kGood) single_lepton_nm1_higgsMassDiff.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
      if((masked_fc & ~kDeltaR) == kGood) single_lepton_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
      if((masked_fc & ~kMETSig30) == kGood) single_lepton_nm1_met.Fill(pfmets_et->at(0),localWeight);
      if((masked_fc & ~kMETSig30) == kGood) single_lepton_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
    }

    masked_fc = failure_code & ~(kMETSig150 | kMETSig100 | kMETSig50 | kHiggsAvgMass | kHiggsMassDiff | k4thBTag | k3rdBTag);
    if(!(masked_fc & ~k3rdBTag)) sig_sb_nm1_thirdCSV.Fill(GetHighestCSV(3),localWeight);
    if(!(masked_fc & ~k4thBTag)) sig_sb_nm1_fourthCSV.Fill(GetHighestCSV(4),localWeight);
    if(!(masked_fc & ~kHiggsAvgMass)) sig_sb_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
    if(!(masked_fc & ~kHiggsMassDiff)) sig_sb_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
    if(!(masked_fc & ~kDeltaR)) sig_sb_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
    if(!masked_fc){
      sig_sb_nm1_met.Fill(pfmets_et->at(0),localWeight);
      sig_sb_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
    }
    if(!(masked_fc & ~kNumJets)) sig_sb_nm1_numJets.Fill(GetNumGoodJets(),localWeight);
    if(!(masked_fc & ~kMinDeltaPhi)) sig_sb_nm1_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),localWeight);

    if(PassesPVCut() && PassesMETCleaningCut() && PassesTriggerCut() && PassesNumJetsCut() && Passes2CSVTCut() && PassesJet2PtCut() && PassesMinDeltaPhiCut() && PassesLeptonVetoCut() && PassesIsoTrackVetoCut()){

      if(PassesBTaggingCut() && PassesHiggsMassCut() && PassesDRCut()){
        //MET S/sqrt(B) plots go here!
        sorb_MET.Fill(pfmets_et->at(0), localWeight);
        sorb_METSig.Fill(pfmets_fullSignif, localWeight);
        if(GetcfAVersion()>=69){
          sorb_METSig.Fill(pfmets_fullSignif, localWeight);
        }
        sorb_METOverSqrtHT.Fill(pfmets_et->at(0)/sqrt(GetHT()),localWeight);
      }

      if(GetNumCSVMJets()<=2 && GetNumCSVLJets()<=3){
        //2b control region plots go here!
        if(PassesMETSig30Cut() && PassesHiggsMassCut() && PassesDRCut()){
          twoB_nm1_thirdBTag.Fill(GetHighestCSV(3),localWeight);
          twoB_nm1_fourthBTag.Fill(GetHighestCSV(4),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesDRCut()){
          twoB_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsAvgMassCut() && PassesDRCut()){
          twoB_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsMassCut()){
          twoB_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
        }
        if(PassesHiggsMassCut() && PassesDRCut()){
          twoB_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
          twoB_nm1_met.Fill(pfmets_et->at(0),localWeight);
          twoB_METOverSqrtHT_vs_METSig.Fill(pfmets_fullSignif, pfmets_et->at(0)/sqrt(GetHT()), localWeight);
          twoB_METOverSqrtHT_Over_METSig.Fill((pfmets_et->at(0)/sqrt(GetHT()))/pfmets_fullSignif, localWeight);
          twoB_MET.Fill(pfmets_et->at(0),localWeight);
          twoB_MET_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),pfmets_et->at(0),localWeight);
          twoB_METSig_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),pfmets_fullSignif,localWeight);
          twoB_METToMETSigRatio_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
          twoB_MET_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0),localWeight);
          twoB_METSig_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_fullSignif,localWeight);
          twoB_METToMETSigRatio_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
          twoB_MET_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0),localWeight);
          twoB_METSig_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_fullSignif,localWeight);
        }
      }

      if(0.5*(higgsMasses.first+higgsMasses.second)<75.0){
        //Higgs mass control region plots go here!
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesDRCut()){
          lightHiggs_nm1_thirdBTag.Fill(GetHighestCSV(3),localWeight);
          lightHiggs_nm1_fourthBTag.Fill(GetHighestCSV(4),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesDRCut() && PassesBTaggingCut()){
          lightHiggs_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
        }
        if(PassesMETSig30Cut() && PassesDRCut() && PassesBTaggingCut()){
          lightHiggs_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesBTaggingCut()){
          lightHiggs_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
        }
        if(PassesHiggsMassDiffCut() && PassesDRCut() && PassesBTaggingCut()){
          lightHiggs_METOverSqrtHT_vs_METSig.Fill(pfmets_fullSignif, pfmets_et->at(0)/sqrt(GetHT()), localWeight);
          lightHiggs_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
          lightHiggs_nm1_met.Fill(pfmets_et->at(0),localWeight);
          lightHiggs_METOverSqrtHT_Over_METSig.Fill((pfmets_et->at(0)/sqrt(GetHT()))/pfmets_fullSignif, localWeight);
          lightHiggs_MET.Fill(pfmets_et->at(0),localWeight);
          lightHiggs_MET_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),pfmets_et->at(0),localWeight);
          lightHiggs_METSig_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),pfmets_fullSignif,localWeight);
          lightHiggs_METToMETSigRatio_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
          lightHiggs_MET_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0),localWeight);
          lightHiggs_METSig_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_fullSignif,localWeight);
          lightHiggs_METToMETSigRatio_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
          lightHiggs_MET_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0),localWeight);
          lightHiggs_METSig_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_fullSignif,localWeight);
        }
      }

      //Signal region plots go here!
      if(PassesMETSig30Cut() && PassesHiggsMassCut() && PassesDRCut()){
        unsigned int CSV1(sortedBJetCache[0].GetIndex()), CSV2(sortedBJetCache[1].GetIndex()),
                     CSV3(sortedBJetCache[2].GetIndex()), CSV4(sortedBJetCache[3].GetIndex());
        jet_parton_Id_CSV1.Fill(GetPartonIdBin(jets_AK5PF_parton_Id->at(CSV1)),localWeight);
        jet_parton_Id_CSV2.Fill(GetPartonIdBin(jets_AK5PF_parton_Id->at(CSV2)),localWeight);
        jet_parton_Id_CSV3.Fill(GetPartonIdBin(jets_AK5PF_parton_Id->at(CSV3)),localWeight);
        jet_parton_Id_CSV4.Fill(GetPartonIdBin(jets_AK5PF_parton_Id->at(CSV4)),localWeight);
        xx_nm1_thirdBTag.Fill(GetHighestCSV(3),localWeight);
        if(GetNumCSVMJets()>=3){
          xx_nm1_fourthBTag.Fill(GetHighestCSV(4),localWeight);
        }
      }
      if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesBTaggingCut() && PassesDRCut()){
        xx_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
      }
      if(PassesMETSig30Cut() && PassesHiggsAvgMassCut() && PassesBTaggingCut() && PassesDRCut()){
        xx_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
      }
      if(PassesMETSig30Cut() && PassesHiggsMassCut() && PassesBTaggingCut()){
        xx_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
      }
      if(PassesHiggsMassCut() && PassesBTaggingCut() && PassesDRCut()){
        xx_METOverSqrtHT_vs_METSig.Fill(pfmets_fullSignif, pfmets_et->at(0)/sqrt(GetHT()), localWeight);
        xx_METOverSqrtHT_Over_METSig.Fill((pfmets_et->at(0)/sqrt(GetHT()))/pfmets_fullSignif, localWeight);
        xx_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
        xx_nm1_met.Fill(pfmets_et->at(0),localWeight);
        xx_MET.Fill(pfmets_et->at(0),localWeight);
        xx_MET_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),pfmets_et->at(0),localWeight);
        xx_METSig_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(UINT_MAX),pfmets_fullSignif,localWeight);
      }
    }

    xx_METToMETSigRatio_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
    xx_MET_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0),localWeight);
    xx_METSig_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_fullSignif,localWeight);
    xx_METToMETSigRatio_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
    xx_MET_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0),localWeight);
    xx_METSig_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_fullSignif,localWeight);

    unsigned int this_sbin(0);
    if(pfmets_fullSignif>30.0 && pfmets_fullSignif<50.0){
      this_sbin=1;
    }else if(pfmets_fullSignif>50.0 && pfmets_fullSignif<100.0){
      this_sbin=2;
    }else if(pfmets_fullSignif>100.0 && pfmets_fullSignif<150.0){
      this_sbin=3;
    }else if(pfmets_fullSignif>150.0){
      this_sbin=4;
    }
    if(PassesRegionACut()){
      switch(this_sbin){
      case 1: ++Asbin1Count;
        Asbin1CountWeighted+=localWeight;
        break;
      case 2: ++Asbin2Count;
        Asbin2CountWeighted+=localWeight;
        break;
      case 3: ++Asbin3Count;
        Asbin3CountWeighted+=localWeight;
        break;
      case 4: ++Asbin4Count;
        Asbin4CountWeighted+=localWeight;
      }
      ++ACount;
      ACountWeighted+=localWeight;
      xx_ABCD.Fill(0.0,localWeight);
      xx_metSig_A.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionBCut()){
      switch(this_sbin){
      case 1: ++Bsbin1Count;
        Bsbin1CountWeighted+=localWeight;
        break;
      case 2: ++Bsbin2Count;
        Bsbin2CountWeighted+=localWeight;
        break;
      case 3: ++Bsbin3Count;
        Bsbin3CountWeighted+=localWeight;
        break;
      case 4: ++Bsbin4Count;
        Bsbin4CountWeighted+=localWeight;
      }
      ++BCount;
      BCountWeighted+=localWeight;
      xx_ABCD.Fill(1.0,localWeight);
      xx_metSig_B.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionC3bCut()){
      switch(this_sbin){
      case 1: ++C3bsbin1Count;
        C3bsbin1CountWeighted+=localWeight;
        break;
      case 2: ++C3bsbin2Count;
        C3bsbin2CountWeighted+=localWeight;
        break;
      case 3: ++C3bsbin3Count;
        C3bsbin3CountWeighted+=localWeight;
        break;
      case 4: ++C3bsbin4Count;
        C3bsbin4CountWeighted+=localWeight;
      }
      ++C3bCount;
      C3bCountWeighted+=localWeight;
      xx_ABCD.Fill(2.0,localWeight);
      xx_metSig_C3b.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionD3bCut()){
      switch(this_sbin){
      case 1: ++D3bsbin1Count;
        D3bsbin1CountWeighted+=localWeight;
        break;
      case 2: ++D3bsbin2Count;
        D3bsbin2CountWeighted+=localWeight;
        break;
      case 3: ++D3bsbin3Count;
        D3bsbin3CountWeighted+=localWeight;
        break;
      case 4: ++D3bsbin4Count;
        D3bsbin4CountWeighted+=localWeight;
      }
      ++D3bCount;
      D3bCountWeighted+=localWeight;
      xx_ABCD.Fill(3.0,localWeight);
      xx_metSig_D3b.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionC2bCut()){
      switch(this_sbin){
      case 1: ++C2bsbin1Count;
        C2bsbin1CountWeighted+=localWeight;
        break;
      case 2: ++C2bsbin2Count;
        C2bsbin2CountWeighted+=localWeight;
        break;
      case 3: ++C2bsbin3Count;
        C2bsbin3CountWeighted+=localWeight;
        break;
      case 4: ++C2bsbin4Count;
        C2bsbin4CountWeighted+=localWeight;
      }
      ++C2bCount;
      C2bCountWeighted+=localWeight;
      xx_ABCD.Fill(4.0,localWeight);
      xx_metSig_C2b.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionD2bCut()){
      switch(this_sbin){
      case 1: ++D2bsbin1Count;
        D2bsbin1CountWeighted+=localWeight;
        break;
      case 2: ++D2bsbin2Count;
        D2bsbin2CountWeighted+=localWeight;
        break;
      case 3: ++D2bsbin3Count;
        D2bsbin3CountWeighted+=localWeight;
        break;
      case 4: ++D2bsbin4Count;
        D2bsbin4CountWeighted+=localWeight;
      }
      ++D2bCount;
      D2bCountWeighted+=localWeight;
      xx_ABCD.Fill(5.0,localWeight);
      xx_metSig_D2b.Fill(pfmets_fullSignif,localWeight);
    }

    if(PassesSingleLeptonRegionACut()){
      ++ASLCount;
      ASLCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++ASLsbin1Count;
        ASLsbin1CountWeighted+=localWeight;
        break;
      case 2: ++ASLsbin2Count;
        ASLsbin2CountWeighted+=localWeight;
        break;
      case 3: ++ASLsbin3Count;
        ASLsbin3CountWeighted+=localWeight;
        break;
      case 4: ++ASLsbin4Count;
        ASLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionBCut()){
      ++BSLCount;
      BSLCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++BSLsbin1Count;
        BSLsbin1CountWeighted+=localWeight;
        break;
      case 2: ++BSLsbin2Count;
        BSLsbin2CountWeighted+=localWeight;
        break;
      case 3: ++BSLsbin3Count;
        BSLsbin3CountWeighted+=localWeight;
        break;
      case 4: ++BSLsbin4Count;
        BSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionC3bCut()){
      ++C3bSLCount;
      C3bSLCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++C3bSLsbin1Count;
        C3bSLsbin1CountWeighted+=localWeight;
        break;
      case 2: ++C3bSLsbin2Count;
        C3bSLsbin2CountWeighted+=localWeight;
        break;
      case 3: ++C3bSLsbin3Count;
        C3bSLsbin3CountWeighted+=localWeight;
        break;
      case 4: ++C3bSLsbin4Count;
        C3bSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionD3bCut()){
      ++D3bSLCount;
      D3bSLCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++D3bSLsbin1Count;
        D3bSLsbin1CountWeighted+=localWeight;
        break;
      case 2: ++D3bSLsbin2Count;
        D3bSLsbin2CountWeighted+=localWeight;
        break;
      case 3: ++D3bSLsbin3Count;
        D3bSLsbin3CountWeighted+=localWeight;
        break;
      case 4: ++D3bSLsbin4Count;
        D3bSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionC2bCut()){
      ++C2bSLCount;
      C2bSLCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++C2bSLsbin1Count;
        C2bSLsbin1CountWeighted+=localWeight;
        break;
      case 2: ++C2bSLsbin2Count;
        C2bSLsbin2CountWeighted+=localWeight;
        break;
      case 3: ++C2bSLsbin3Count;
        C2bSLsbin3CountWeighted+=localWeight;
        break;
      case 4: ++C2bSLsbin4Count;
        C2bSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionD2bCut()){
      ++D2bSLCount;
      D2bSLCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++D2bSLsbin1Count;
        D2bSLsbin1CountWeighted+=localWeight;
        break;
      case 2: ++D2bSLsbin2Count;
        D2bSLsbin2CountWeighted+=localWeight;
        break;
      case 3: ++D2bSLsbin3Count;
        D2bSLsbin3CountWeighted+=localWeight;
        break;
      case 4: ++D2bSLsbin4Count;
        D2bSLsbin4CountWeighted+=localWeight;
      }
    }

    if(PassesInvertedDRRegionACut()){
      ++ADRInvCount;
      ADRInvCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++ADRInvsbin1Count;
        ADRInvsbin1CountWeighted+=localWeight;
        break;
      case 2: ++ADRInvsbin2Count;
        ADRInvsbin2CountWeighted+=localWeight;
        break;
      case 3: ++ADRInvsbin3Count;
        ADRInvsbin3CountWeighted+=localWeight;
        break;
      case 4: ++ADRInvsbin4Count;
        ADRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionBCut()){
      ++BDRInvCount;
      BDRInvCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++BDRInvsbin1Count;
        BDRInvsbin1CountWeighted+=localWeight;
        break;
      case 2: ++BDRInvsbin2Count;
        BDRInvsbin2CountWeighted+=localWeight;
        break;
      case 3: ++BDRInvsbin3Count;
        BDRInvsbin3CountWeighted+=localWeight;
        break;
      case 4: ++BDRInvsbin4Count;
        BDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionC3bCut()){
      ++C3bDRInvCount;
      C3bDRInvCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++C3bDRInvsbin1Count;
        C3bDRInvsbin1CountWeighted+=localWeight;
        break;
      case 2: ++C3bDRInvsbin2Count;
        C3bDRInvsbin2CountWeighted+=localWeight;
        break;
      case 3: ++C3bDRInvsbin3Count;
        C3bDRInvsbin3CountWeighted+=localWeight;
        break;
      case 4: ++C3bDRInvsbin4Count;
        C3bDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionD3bCut()){
      ++D3bDRInvCount;
      D3bDRInvCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++D3bDRInvsbin1Count;
        D3bDRInvsbin1CountWeighted+=localWeight;
        break;
      case 2: ++D3bDRInvsbin2Count;
        D3bDRInvsbin2CountWeighted+=localWeight;
        break;
      case 3: ++D3bDRInvsbin3Count;
        D3bDRInvsbin3CountWeighted+=localWeight;
        break;
      case 4: ++D3bDRInvsbin4Count;
        D3bDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionC2bCut()){
      ++C2bDRInvCount;
      C2bDRInvCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++C2bDRInvsbin1Count;
        C2bDRInvsbin1CountWeighted+=localWeight;
        break;
      case 2: ++C2bDRInvsbin2Count;
        C2bDRInvsbin2CountWeighted+=localWeight;
        break;
      case 3: ++C2bDRInvsbin3Count;
        C2bDRInvsbin3CountWeighted+=localWeight;
        break;
      case 4: ++C2bDRInvsbin4Count;
        C2bDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionD2bCut()){
      ++D2bDRInvCount;
      D2bDRInvCountWeighted+=localWeight;
      switch(this_sbin){
      case 1: ++D2bDRInvsbin1Count;
        D2bDRInvsbin1CountWeighted+=localWeight;
        break;
      case 2: ++D2bDRInvsbin2Count;
        D2bDRInvsbin2CountWeighted+=localWeight;
        break;
      case 3: ++D2bDRInvsbin3Count;
        D2bDRInvsbin3CountWeighted+=localWeight;
        break;
      case 4: ++D2bDRInvsbin4Count;
        D2bDRInvsbin4CountWeighted+=localWeight;
      }
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
  for(int i(1); i<=xx_metSig_A.GetNbinsX(); ++i){
    xx_metSig_A.SetBinContent(i,xx_metSig_A.GetBinContent(i)*10.0/xx_metSig_A.GetBinWidth(i));
    xx_metSig_B.SetBinContent(i,xx_metSig_B.GetBinContent(i)*10.0/xx_metSig_B.GetBinWidth(i));
    xx_metSig_C3b.SetBinContent(i,xx_metSig_C3b.GetBinContent(i)*10.0/xx_metSig_C3b.GetBinWidth(i));
    xx_metSig_D3b.SetBinContent(i,xx_metSig_D3b.GetBinContent(i)*10.0/xx_metSig_D3b.GetBinWidth(i));
    xx_metSig_C2b.SetBinContent(i,xx_metSig_C2b.GetBinContent(i)*10.0/xx_metSig_C2b.GetBinWidth(i));
    xx_metSig_D2b.SetBinContent(i,xx_metSig_D2b.GetBinContent(i)*10.0/xx_metSig_D2b.GetBinWidth(i));
  }

  xx_metSig_SigOverSB_4b.Divide(&xx_metSig_A,&xx_metSig_B);
  xx_metSig_SigOverSB_3b.Divide(&xx_metSig_C3b,&xx_metSig_D3b);
  xx_metSig_SigOverSB_2b.Divide(&xx_metSig_C2b,&xx_metSig_D2b);
  xx_metSig_SigOverSB_4b.SetLineColor(1);
  xx_metSig_SigOverSB_3b.SetLineColor(2);
  xx_metSig_SigOverSB_2b.SetLineColor(3);
  xx_metSig_SigOverSB_2b.Draw("e1p");
  xx_metSig_SigOverSB_3b.Draw("e1psame");
  xx_metSig_SigOverSB_4b.Draw("e1psame");
  xx_metSig_SigOverSB.Write();

  double CtotCountWeighted(C3bCountWeighted+C2bCountWeighted);
  double DtotCountWeighted(D3bCountWeighted+D2bCountWeighted);

  double kappa3b((ACountWeighted/BCountWeighted)/(C3bCountWeighted/D3bCountWeighted));
  double kappa2b((ACountWeighted/BCountWeighted)/(C2bCountWeighted/D2bCountWeighted));
  double kappatot((ACountWeighted/BCountWeighted)/(CtotCountWeighted/DtotCountWeighted));

  TTree ABCD("ABCD","ABCD");
  ABCD.Branch("A", &ACount);
  ABCD.Branch("B", &BCount);
  ABCD.Branch("C3b", &C3bCount);
  ABCD.Branch("D3b", &D3bCount);
  ABCD.Branch("C2b", &C2bCount);
  ABCD.Branch("D2b", &D2bCount);
  ABCD.Branch("AWeighted", &ACountWeighted);
  ABCD.Branch("BWeighted", &BCountWeighted);
  ABCD.Branch("C3bWeighted", &C3bCountWeighted);
  ABCD.Branch("D3bWeighted", &D3bCountWeighted);
  ABCD.Branch("C2bWeighted", &C2bCountWeighted);
  ABCD.Branch("D2bWeighted", &D2bCountWeighted);
  ABCD.Branch("CtotWeighted", &CtotCountWeighted);
  ABCD.Branch("DtotWeighted", &DtotCountWeighted);
  ABCD.Branch("A_DRInv", &ADRInvCount);
  ABCD.Branch("B_DRInv", &BDRInvCount);
  ABCD.Branch("C3b_DRInv", &C3bDRInvCount);
  ABCD.Branch("D3b_DRInv", &D3bDRInvCount);
  ABCD.Branch("C2b_DRInv", &C2bDRInvCount);
  ABCD.Branch("D2b_DRInv", &D2bDRInvCount);
  ABCD.Branch("AWeighted_DRInv", &ADRInvCountWeighted);
  ABCD.Branch("BWeighted_DRInv", &BDRInvCountWeighted);
  ABCD.Branch("C3bWeighted_DRInv", &C3bDRInvCountWeighted);
  ABCD.Branch("D3bWeighted_DRInv", &D3bDRInvCountWeighted);
  ABCD.Branch("C2bWeighted_DRInv", &C2bDRInvCountWeighted);
  ABCD.Branch("D2bWeighted_DRInv", &D2bDRInvCountWeighted);
  ABCD.Branch("A_SL", &ASLCount);
  ABCD.Branch("B_SL", &BSLCount);
  ABCD.Branch("C3b_SL", &C3bSLCount);
  ABCD.Branch("D3b_SL", &D3bSLCount);
  ABCD.Branch("C2b_SL", &C2bSLCount);
  ABCD.Branch("D2b_SL", &D2bSLCount);
  ABCD.Branch("AWeighted_SL", &ASLCountWeighted);
  ABCD.Branch("BWeighted_SL", &BSLCountWeighted);
  ABCD.Branch("C3bWeighted_SL", &C3bSLCountWeighted);
  ABCD.Branch("D3bWeighted_SL", &D3bSLCountWeighted);
  ABCD.Branch("C2bWeighted_SL", &C2bSLCountWeighted);
  ABCD.Branch("D2bWeighted_SL", &D2bSLCountWeighted);
  ABCD.Branch("A_SLsbin1", &ASLsbin1Count);
  ABCD.Branch("B_SLsbin1", &BSLsbin1Count);
  ABCD.Branch("C3b_SLsbin1", &C3bSLsbin1Count);
  ABCD.Branch("D3b_SLsbin1", &D3bSLsbin1Count);
  ABCD.Branch("C2b_SLsbin1", &C2bSLsbin1Count);
  ABCD.Branch("D2b_SLsbin1", &D2bSLsbin1Count);
  ABCD.Branch("AWeighted_SLsbin1", &ASLsbin1CountWeighted);
  ABCD.Branch("BWeighted_SLsbin1", &BSLsbin1CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin1", &C3bSLsbin1CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin1", &D3bSLsbin1CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin1", &C2bSLsbin1CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin1", &D2bSLsbin1CountWeighted);
  ABCD.Branch("A_SLsbin2", &ASLsbin2Count);
  ABCD.Branch("B_SLsbin2", &BSLsbin2Count);
  ABCD.Branch("C3b_SLsbin2", &C3bSLsbin2Count);
  ABCD.Branch("D3b_SLsbin2", &D3bSLsbin2Count);
  ABCD.Branch("C2b_SLsbin2", &C2bSLsbin2Count);
  ABCD.Branch("D2b_SLsbin2", &D2bSLsbin2Count);
  ABCD.Branch("AWeighted_SLsbin2", &ASLsbin2CountWeighted);
  ABCD.Branch("BWeighted_SLsbin2", &BSLsbin2CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin2", &C3bSLsbin2CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin2", &D3bSLsbin2CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin2", &C2bSLsbin2CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin2", &D2bSLsbin2CountWeighted);
  ABCD.Branch("A_SLsbin3", &ASLsbin3Count);
  ABCD.Branch("B_SLsbin3", &BSLsbin3Count);
  ABCD.Branch("C3b_SLsbin3", &C3bSLsbin3Count);
  ABCD.Branch("D3b_SLsbin3", &D3bSLsbin3Count);
  ABCD.Branch("C2b_SLsbin3", &C2bSLsbin3Count);
  ABCD.Branch("D2b_SLsbin3", &D2bSLsbin3Count);
  ABCD.Branch("AWeighted_SLsbin3", &ASLsbin3CountWeighted);
  ABCD.Branch("BWeighted_SLsbin3", &BSLsbin3CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin3", &C3bSLsbin3CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin3", &D3bSLsbin3CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin3", &C2bSLsbin3CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin3", &D2bSLsbin3CountWeighted);
  ABCD.Branch("A_SLsbin4", &ASLsbin4Count);
  ABCD.Branch("B_SLsbin4", &BSLsbin4Count);
  ABCD.Branch("C3b_SLsbin4", &C3bSLsbin4Count);
  ABCD.Branch("D3b_SLsbin4", &D3bSLsbin4Count);
  ABCD.Branch("C2b_SLsbin4", &C2bSLsbin4Count);
  ABCD.Branch("D2b_SLsbin4", &D2bSLsbin4Count);
  ABCD.Branch("AWeighted_SLsbin4", &ASLsbin4CountWeighted);
  ABCD.Branch("BWeighted_SLsbin4", &BSLsbin4CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin4", &C3bSLsbin4CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin4", &D3bSLsbin4CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin4", &C2bSLsbin4CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin4", &D2bSLsbin4CountWeighted);
  ABCD.Branch("A_DRInvsbin1", &ADRInvsbin1Count);
  ABCD.Branch("B_DRInvsbin1", &BDRInvsbin1Count);
  ABCD.Branch("C3b_DRInvsbin1", &C3bDRInvsbin1Count);
  ABCD.Branch("D3b_DRInvsbin1", &D3bDRInvsbin1Count);
  ABCD.Branch("C2b_DRInvsbin1", &C2bDRInvsbin1Count);
  ABCD.Branch("D2b_DRInvsbin1", &D2bDRInvsbin1Count);
  ABCD.Branch("AWeighted_DRInvsbin1", &ADRInvsbin1CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin1", &BDRInvsbin1CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin1", &C3bDRInvsbin1CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin1", &D3bDRInvsbin1CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin1", &C2bDRInvsbin1CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin1", &D2bDRInvsbin1CountWeighted);
  ABCD.Branch("A_DRInvsbin2", &ADRInvsbin2Count);
  ABCD.Branch("B_DRInvsbin2", &BDRInvsbin2Count);
  ABCD.Branch("C3b_DRInvsbin2", &C3bDRInvsbin2Count);
  ABCD.Branch("D3b_DRInvsbin2", &D3bDRInvsbin2Count);
  ABCD.Branch("C2b_DRInvsbin2", &C2bDRInvsbin2Count);
  ABCD.Branch("D2b_DRInvsbin2", &D2bDRInvsbin2Count);
  ABCD.Branch("AWeighted_DRInvsbin2", &ADRInvsbin2CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin2", &BDRInvsbin2CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin2", &C3bDRInvsbin2CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin2", &D3bDRInvsbin2CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin2", &C2bDRInvsbin2CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin2", &D2bDRInvsbin2CountWeighted);
  ABCD.Branch("A_DRInvsbin3", &ADRInvsbin3Count);
  ABCD.Branch("B_DRInvsbin3", &BDRInvsbin3Count);
  ABCD.Branch("C3b_DRInvsbin3", &C3bDRInvsbin3Count);
  ABCD.Branch("D3b_DRInvsbin3", &D3bDRInvsbin3Count);
  ABCD.Branch("C2b_DRInvsbin3", &C2bDRInvsbin3Count);
  ABCD.Branch("D2b_DRInvsbin3", &D2bDRInvsbin3Count);
  ABCD.Branch("AWeighted_DRInvsbin3", &ADRInvsbin3CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin3", &BDRInvsbin3CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin3", &C3bDRInvsbin3CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin3", &D3bDRInvsbin3CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin3", &C2bDRInvsbin3CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin3", &D2bDRInvsbin3CountWeighted);
  ABCD.Branch("A_DRInvsbin4", &ADRInvsbin4Count);
  ABCD.Branch("B_DRInvsbin4", &BDRInvsbin4Count);
  ABCD.Branch("C3b_DRInvsbin4", &C3bDRInvsbin4Count);
  ABCD.Branch("D3b_DRInvsbin4", &D3bDRInvsbin4Count);
  ABCD.Branch("C2b_DRInvsbin4", &C2bDRInvsbin4Count);
  ABCD.Branch("D2b_DRInvsbin4", &D2bDRInvsbin4Count);
  ABCD.Branch("AWeighted_DRInvsbin4", &ADRInvsbin4CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin4", &BDRInvsbin4CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin4", &C3bDRInvsbin4CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin4", &D3bDRInvsbin4CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin4", &C2bDRInvsbin4CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin4", &D2bDRInvsbin4CountWeighted);
  ABCD.Branch("A_sbin1", &Asbin1Count);
  ABCD.Branch("B_sbin1", &Bsbin1Count);
  ABCD.Branch("C3b_sbin1", &C3bsbin1Count);
  ABCD.Branch("D3b_sbin1", &D3bsbin1Count);
  ABCD.Branch("C2b_sbin1", &C2bsbin1Count);
  ABCD.Branch("D2b_sbin1", &D2bsbin1Count);
  ABCD.Branch("AWeighted_sbin1", &Asbin1CountWeighted);
  ABCD.Branch("BWeighted_sbin1", &Bsbin1CountWeighted);
  ABCD.Branch("C3bWeighted_sbin1", &C3bsbin1CountWeighted);
  ABCD.Branch("D3bWeighted_sbin1", &D3bsbin1CountWeighted);
  ABCD.Branch("C2bWeighted_sbin1", &C2bsbin1CountWeighted);
  ABCD.Branch("D2bWeighted_sbin1", &D2bsbin1CountWeighted);
  ABCD.Branch("A_sbin2", &Asbin2Count);
  ABCD.Branch("B_sbin2", &Bsbin2Count);
  ABCD.Branch("C3b_sbin2", &C3bsbin2Count);
  ABCD.Branch("D3b_sbin2", &D3bsbin2Count);
  ABCD.Branch("C2b_sbin2", &C2bsbin2Count);
  ABCD.Branch("D2b_sbin2", &D2bsbin2Count);
  ABCD.Branch("AWeighted_sbin2", &Asbin2CountWeighted);
  ABCD.Branch("BWeighted_sbin2", &Bsbin2CountWeighted);
  ABCD.Branch("C3bWeighted_sbin2", &C3bsbin2CountWeighted);
  ABCD.Branch("D3bWeighted_sbin2", &D3bsbin2CountWeighted);
  ABCD.Branch("C2bWeighted_sbin2", &C2bsbin2CountWeighted);
  ABCD.Branch("D2bWeighted_sbin2", &D2bsbin2CountWeighted);
  ABCD.Branch("A_sbin3", &Asbin3Count);
  ABCD.Branch("B_sbin3", &Bsbin3Count);
  ABCD.Branch("C3b_sbin3", &C3bsbin3Count);
  ABCD.Branch("D3b_sbin3", &D3bsbin3Count);
  ABCD.Branch("C2b_sbin3", &C2bsbin3Count);
  ABCD.Branch("D2b_sbin3", &D2bsbin3Count);
  ABCD.Branch("AWeighted_sbin3", &Asbin3CountWeighted);
  ABCD.Branch("BWeighted_sbin3", &Bsbin3CountWeighted);
  ABCD.Branch("C3bWeighted_sbin3", &C3bsbin3CountWeighted);
  ABCD.Branch("D3bWeighted_sbin3", &D3bsbin3CountWeighted);
  ABCD.Branch("C2bWeighted_sbin3", &C2bsbin3CountWeighted);
  ABCD.Branch("D2bWeighted_sbin3", &D2bsbin3CountWeighted);
  ABCD.Branch("A_sbin4", &Asbin4Count);
  ABCD.Branch("B_sbin4", &Bsbin4Count);
  ABCD.Branch("C3b_sbin4", &C3bsbin4Count);
  ABCD.Branch("D3b_sbin4", &D3bsbin4Count);
  ABCD.Branch("C2b_sbin4", &C2bsbin4Count);
  ABCD.Branch("D2b_sbin4", &D2bsbin4Count);
  ABCD.Branch("AWeighted_sbin4", &Asbin4CountWeighted);
  ABCD.Branch("BWeighted_sbin4", &Bsbin4CountWeighted);
  ABCD.Branch("C3bWeighted_sbin4", &C3bsbin4CountWeighted);
  ABCD.Branch("D3bWeighted_sbin4", &D3bsbin4CountWeighted);
  ABCD.Branch("C2bWeighted_sbin4", &C2bsbin4CountWeighted);
  ABCD.Branch("D2bWeighted_sbin4", &D2bsbin4CountWeighted);
  ABCD.Branch("kappa3b", &kappa3b);
  ABCD.Branch("kappa2b", &kappa2b);
  ABCD.Branch("kappatot", &kappatot);
  ABCD.Fill();

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

  file.Write();
  file.Close();
}

namespace CfAPlots{
  void FixSbinLabels(TH1D &h){
    h.GetXaxis()->SetBinLabel(1,"1");
    h.GetXaxis()->SetBinLabel(2,"2");
    h.GetXaxis()->SetBinLabel(3,"3");
    h.GetXaxis()->SetBinLabel(4,"4");
  }

  void FixBinLabels(TH1D &h){
    h.GetXaxis()->SetBinLabel(1, "H->b");
    h.GetXaxis()->SetBinLabel(2, "t->b");
    h.GetXaxis()->SetBinLabel(3, "W->b");
    h.GetXaxis()->SetBinLabel(4, "g->b");
    h.GetXaxis()->SetBinLabel(5, "Z->b");
    h.GetXaxis()->SetBinLabel(6, "other->b");
    h.GetXaxis()->SetBinLabel(7, "W->c");
    h.GetXaxis()->SetBinLabel(8, "g->c");
    h.GetXaxis()->SetBinLabel(9, "Z->c");
    h.GetXaxis()->SetBinLabel(10, "other->c");
    h.GetXaxis()->SetBinLabel(11, "W->uds");
    h.GetXaxis()->SetBinLabel(12, "g->uds");
    h.GetXaxis()->SetBinLabel(13, "Z->uds");
    h.GetXaxis()->SetBinLabel(14, "other->uds");
    h.GetXaxis()->SetBinLabel(15, "gluon");
    h.GetXaxis()->SetBinLabel(16, "electron");
    h.GetXaxis()->SetBinLabel(17, "muon");
    h.GetXaxis()->SetBinLabel(18, "tau");
    h.GetXaxis()->SetBinLabel(19, "other");
  }

  int GetType(const std::pair<int,int> &b_origin){
    if(b_origin.first==5){
      if(b_origin.second==25){
        return 1;
      }else if(b_origin.second==6){
        return 2;
      }else if(b_origin.second==24){
        return 3;
      }else if(b_origin.second==21){
        return 4;
      }else if(b_origin.second==23){
        return 5;
      }else{
        return 6;
      }
    }else if(b_origin.first==4){
      if(b_origin.second==24){
        return 7;
      }else if(b_origin.second==21){
        return 8;
      }else if(b_origin.second==23){
        return 9;
      }else{
        return 10;
      }
    }else if(b_origin.first>=1 && b_origin.first<=3){
      if(b_origin.second==24){
        return 11;
      }else if(b_origin.second==21){
        return 12;
      }else if(b_origin.second==23){
        return 13;
      }else{
        return 14;
      }
    }else if(b_origin.first==21){
      return 15;
    }else if(b_origin.first==11){
      return 16;
    }else if(b_origin.first==13){
      return 17;
    }else if(b_origin.first==15){
      return 18;
    }else{
      return 19;
    }
  }
}
