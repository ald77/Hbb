#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH2D.h"
#include "TBox.h"
#include "TLine.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "style.hpp"
#include "timer.hpp"
#include "utils.hpp"
#include "plotter.hpp"

void compare_histos(TH1D h1, TH1D h2, TH1D h3, const std::string out_name,
		    const bool normalize=false){
  TCanvas c;
  h1.SetLineColor(2);
  h2.SetLineColor(3);
  h3.SetLineColor(4);
  h1.SetFillColor(0);
  h2.SetFillColor(0);
  h3.SetFillColor(0);
  h1.SetStats(0);
  h2.SetStats(0);
  h3.SetStats(0);
  h1.Draw("hist");
  h2.Draw("histsame");
  h3.Draw("histsame");
  if(normalize){
    h1.Scale(1.0/h1.Integral("width"));
    h2.Scale(1.0/h2.Integral("width"));
    h3.Scale(1.0/h3.Integral("width"));
    h1.GetYaxis()->SetTitle("Normalized Frequency");
    h2.GetYaxis()->SetTitle("Normalized Frequency");
    h3.GetYaxis()->SetTitle("Normalized Frequency");
  }
  TLegend l(0.7, 0.5, 0.95, 0.85);
  l.AddEntry(&h1, "2 b-tags", "l");
  l.AddEntry(&h2, "3 b-tags", "l");
  l.AddEntry(&h3, "4 b-tags", "l");
  l.Draw("same");
  c.Print(out_name.c_str());
}

void sig_sb_scat(TH2D& hist, const std::string& out, const bool print_num=false){
  hist.SetStats(0);
  TCanvas c;
  c.cd();
  c.GetPad(0)->SetTopMargin(0.95);

  const int bx1(1), bx2(18), bx3(21), bx4(28), bx5(31), bx6(-1);
  const int by1(1), by2(4), by3(6), by4(7), by5(-1);
  const double inner(hist.Integral(bx3, bx4, by1, by2));
  const double outer1(hist.Integral(bx1, bx2, by1, by3));
  const double outer2(hist.Integral(bx5, bx6, by1, by3));
  const double outer3(hist.Integral(bx1, bx6, by4, by5));
  const double outer(outer1+outer2+outer3);

  const double scale(4096.0/hist.Integral());
  std::ostringstream oss_opt("");
  oss_opt << "scat=" << scale;
  hist.Draw(oss_opt.str().c_str());

  TLine v1(90.0, 0.0, 90.0, 30.0);
  TLine v2(100.0, 0.0, 100.0, 20.0);
  TLine v3(140.0, 0.0, 140.0, 20.0);
  TLine v4(150.0, 0.0, 150.0, 30.0);
  TLine h1(100.0, 20.0, 140.0, 20.0);
  TLine h2(90.0, 30.0, 150.0, 30.0);
  TBox b1(100.0, 0.0, 140.0, 20.0);
  TBox b2(0.0, 0.0, 90.0, 30.0);
  TBox b3(150.0, 0.0, 250.0, 30.0);
  TBox b4(0.0, 30.0, 250.0, 120.0);
  b1.SetFillColor(3);
  b2.SetFillColor(2);
  b3.SetFillColor(2);
  b4.SetFillColor(2);
  b1.SetFillStyle(3003);
  b2.SetFillStyle(3003);
  b3.SetFillStyle(3003);
  b4.SetFillStyle(3003);
  v1.Draw("same");
  v2.Draw("same");
  v3.Draw("same");
  v4.Draw("same");
  h1.Draw("same");
  h2.Draw("same");
  b1.Draw("same");
  b1.Draw("same");
  b2.Draw("same");
  b3.Draw("same");
  b4.Draw("same");

  std::ostringstream oss_inner("");
  std::ostringstream oss_outer("");
  oss_inner << std::fixed << std::setprecision(0) << inner;
  oss_outer << std::fixed << std::setprecision(0) << outer;

  TPaveText text_inner(100.0, 0.0, 140.0, 20.0);
  text_inner.SetX1(100.0);
  text_inner.SetX2(140.0);
  text_inner.SetY1(0.0);
  text_inner.SetY2(20.0);
  text_inner.SetBorderSize(0);
  text_inner.SetShadowColor(0);
  text_inner.SetLineColor(0);
  text_inner.SetLineWidth(0);
  TPaveText text_outer(150.0, 30.0, 190.0, 50.0);
  text_outer.SetX1(150.0);
  text_outer.SetX2(190.0);
  text_outer.SetY1(30.0);
  text_outer.SetY2(50.0);
  text_outer.SetBorderSize(0);
  text_outer.SetShadowColor(0);
  text_outer.SetLineColor(0);
  text_outer.SetLineWidth(0);
  text_inner.AddText(oss_inner.str().c_str());
  text_outer.AddText(oss_outer.str().c_str());
  text_inner.SetFillStyle(0);
  text_inner.SetTextColor(kRed);
  text_outer.SetFillStyle(0);
  text_outer.SetTextColor(kGreen);
  hist.Draw((oss_opt.str()+"same").c_str());
  if(print_num){
    text_inner.Draw("same");
    text_outer.Draw("same");
  }
  c.Print(out.c_str());
}

void scat_plot(TH2D& hist, const std::string& out, const bool print_num=false){
  hist.SetStats(0);
  TCanvas c;
  c.cd();
  c.GetPad(0)->SetTopMargin(0.95);

  const int bx1(1), bx2(18), bx3(21), bx4(28), bx5(31), bx6(-1);
  const int by1(1), by2(4), by3(6), by4(7), by5(-1);
  const double inner(hist.Integral(bx3, bx4, by1, by2));
  const double outer1(hist.Integral(bx1, bx2, by1, by3));
  const double outer2(hist.Integral(bx5, bx6, by1, by3));
  const double outer3(hist.Integral(bx1, bx6, by4, by5));
  const double outer(outer1+outer2+outer3);

  const double scale(4096.0/hist.Integral());
  std::ostringstream oss_opt("");
  oss_opt << "scat=" << scale;
  hist.Draw(oss_opt.str().c_str());

  TLine v1(90.0, 0.0, 90.0, 30.0);
  TLine v2(100.0, 0.0, 100.0, 20.0);
  TLine v3(140.0, 0.0, 140.0, 20.0);
  TLine v4(150.0, 0.0, 150.0, 30.0);
  TLine h1(100.0, 20.0, 140.0, 20.0);
  TLine h2(90.0, 30.0, 150.0, 30.0);
  TBox b1(100.0, 0.0, 140.0, 20.0);
  TBox b2(0.0, 0.0, 90.0, 30.0);
  TBox b3(150.0, 0.0, 250.0, 30.0);
  TBox b4(0.0, 30.0, 250.0, 120.0);
  b1.SetFillColor(3);
  b2.SetFillColor(2);
  b3.SetFillColor(2);
  b4.SetFillColor(2);
  b1.SetFillStyle(3003);
  b2.SetFillStyle(3003);
  b3.SetFillStyle(3003);
  b4.SetFillStyle(3003);
  /*v1.Draw("same");
  v2.Draw("same");
  v3.Draw("same");
  v4.Draw("same");
  h1.Draw("same");
  h2.Draw("same");
  b1.Draw("same");
  b1.Draw("same");
  b2.Draw("same");
  b3.Draw("same");
  b4.Draw("same");*/

  std::ostringstream oss_inner("");
  std::ostringstream oss_outer("");
  oss_inner << std::fixed << std::setprecision(0) << inner;
  oss_outer << std::fixed << std::setprecision(0) << outer;

  TPaveText text_inner(100.0, 0.0, 140.0, 20.0);
  text_inner.SetX1(100.0);
  text_inner.SetX2(140.0);
  text_inner.SetY1(0.0);
  text_inner.SetY2(20.0);
  text_inner.SetBorderSize(0);
  text_inner.SetShadowColor(0);
  text_inner.SetLineColor(0);
  text_inner.SetLineWidth(0);
  TPaveText text_outer(150.0, 30.0, 190.0, 50.0);
  text_outer.SetX1(150.0);
  text_outer.SetX2(190.0);
  text_outer.SetY1(30.0);
  text_outer.SetY2(50.0);
  text_outer.SetBorderSize(0);
  text_outer.SetShadowColor(0);
  text_outer.SetLineColor(0);
  text_outer.SetLineWidth(0);
  text_inner.AddText(oss_inner.str().c_str());
  text_outer.AddText(oss_outer.str().c_str());
  text_inner.SetFillStyle(0);
  text_inner.SetTextColor(kRed);
  text_outer.SetFillStyle(0);
  text_outer.SetTextColor(kGreen);
  hist.Draw((oss_opt.str()+"same").c_str());
  if(print_num){
    text_inner.Draw("same");
    text_outer.Draw("same");
  }
  c.Print(out.c_str());
}

unsigned count_points(const TGraph& graph, const double xmin, const double xmax,
		      const double ymin, const double ymax){
  unsigned count(0);
  const int n_points(graph.GetN());
  for(int point(0); point<n_points; ++point){
    double x(0.0), y(0.0);
    graph.GetPoint(point, x, y);
    if(x>=xmin && x<xmax && y>=ymin && y<ymax) ++count;
  }
  return count;
}

void sig_sb_scat_data(TGraph& graph, const std::string& out, const bool print_num=false){
  graph.SetMarkerStyle(20);
  TH2D hist("h",";<m_{bb}> [GeV];#Delta m_{bb} [GeV]", 50, 0.0, 250.0, 24, 0.0, 120.0);
  hist.SetStats(0);
  TCanvas c;
  c.cd();
  c.GetPad(0)->SetTopMargin(0.95);

  const double dbl_max(std::numeric_limits<double>::max());
  const double inner(count_points(graph, 100.0, 140.0, 0.0, 20.0));
  const double outer1(count_points(graph, 0.0, 90.0, 0.0, 30.0));
  const double outer2(count_points(graph, 150.0, dbl_max, 0.0, 30.0));
  const double outer3(count_points(graph, 0.0, dbl_max, 30.0, dbl_max));
  const double outer(outer1+outer2+outer3);

  hist.Draw();

  TLine v1(90.0, 0.0, 90.0, 30.0);
  TLine v2(100.0, 0.0, 100.0, 20.0);
  TLine v3(140.0, 0.0, 140.0, 20.0);
  TLine v4(150.0, 0.0, 150.0, 30.0);
  TLine h1(100.0, 20.0, 140.0, 20.0);
  TLine h2(90.0, 30.0, 150.0, 30.0);
  TBox b1(100.0, 0.0, 140.0, 20.0);
  TBox b2(0.0, 0.0, 90.0, 30.0);
  TBox b3(150.0, 0.0, 250.0, 30.0);
  TBox b4(0.0, 30.0, 250.0, 120.0);
  b1.SetFillColor(3);
  b2.SetFillColor(2);
  b3.SetFillColor(2);
  b4.SetFillColor(2);
  b1.SetFillStyle(3003);
  b2.SetFillStyle(3003);
  b3.SetFillStyle(3003);
  b4.SetFillStyle(3003);
  v1.Draw("same");
  v2.Draw("same");
  v3.Draw("same");
  v4.Draw("same");
  h1.Draw("same");
  h2.Draw("same");
  b1.Draw("same");
  b1.Draw("same");
  b2.Draw("same");
  b3.Draw("same");
  b4.Draw("same");

  std::ostringstream oss_inner("");
  std::ostringstream oss_outer("");
  oss_inner << std::fixed << std::setprecision(0) << inner;
  oss_outer << std::fixed << std::setprecision(0) << outer;

  TPaveText text_inner(100.0, 0.0, 140.0, 20.0);
  text_inner.SetX1(100.0);
  text_inner.SetX2(140.0);
  text_inner.SetY1(0.0);
  text_inner.SetY2(20.0);
  text_inner.SetBorderSize(0);
  text_inner.SetShadowColor(0);
  text_inner.SetLineColor(0);
  text_inner.SetLineWidth(0);
  TPaveText text_outer(150.0, 30.0, 190.0, 50.0);
  text_outer.SetX1(150.0);
  text_outer.SetX2(190.0);
  text_outer.SetY1(30.0);
  text_outer.SetY2(50.0);
  text_outer.SetBorderSize(0);
  text_outer.SetShadowColor(0);
  text_outer.SetLineColor(0);
  text_outer.SetLineWidth(0);
  text_inner.AddText(oss_inner.str().c_str());
  text_outer.AddText(oss_outer.str().c_str());
  text_inner.SetFillStyle(0);
  text_inner.SetTextColor(kRed);
  text_outer.SetFillStyle(0);
  text_outer.SetTextColor(kGreen);
  graph.Draw("p");
  if(print_num){
    text_inner.Draw("same");
    text_outer.Draw("same");
  }
  c.Print(out.c_str());
}

void scat_plot_data(TGraph& graph, const std::string& out, const bool print_num=false){
  graph.SetMarkerStyle(20);
  TH2D hist("h",";<m_{bb}> [GeV];#Delta R", 35, 0.0, 250.0, 35, 0.0, 5.0);
  hist.SetStats(0);
  TCanvas c;
  c.cd();
  c.GetPad(0)->SetTopMargin(0.95);

  const double dbl_max(std::numeric_limits<double>::max());
  const double inner(count_points(graph, 100.0, 140.0, 0.0, 20.0));
  const double outer1(count_points(graph, 0.0, 90.0, 0.0, 30.0));
  const double outer2(count_points(graph, 150.0, dbl_max, 0.0, 30.0));
  const double outer3(count_points(graph, 0.0, dbl_max, 30.0, dbl_max));
  const double outer(outer1+outer2+outer3);

  hist.Draw();

  TLine v1(90.0, 0.0, 90.0, 30.0);
  TLine v2(100.0, 0.0, 100.0, 20.0);
  TLine v3(140.0, 0.0, 140.0, 20.0);
  TLine v4(150.0, 0.0, 150.0, 30.0);
  TLine h1(100.0, 20.0, 140.0, 20.0);
  TLine h2(90.0, 30.0, 150.0, 30.0);
  TBox b1(100.0, 0.0, 140.0, 20.0);
  TBox b2(0.0, 0.0, 90.0, 30.0);
  TBox b3(150.0, 0.0, 250.0, 30.0);
  TBox b4(0.0, 30.0, 250.0, 120.0);
  b1.SetFillColor(3);
  b2.SetFillColor(2);
  b3.SetFillColor(2);
  b4.SetFillColor(2);
  b1.SetFillStyle(3003);
  b2.SetFillStyle(3003);
  b3.SetFillStyle(3003);
  b4.SetFillStyle(3003);
  /*v1.Draw("same");
  v2.Draw("same");
  v3.Draw("same");
  v4.Draw("same");
  h1.Draw("same");
  h2.Draw("same");
  b1.Draw("same");
  b1.Draw("same");
  b2.Draw("same");
  b3.Draw("same");
  b4.Draw("same");*/

  std::ostringstream oss_inner("");
  std::ostringstream oss_outer("");
  oss_inner << std::fixed << std::setprecision(0) << inner;
  oss_outer << std::fixed << std::setprecision(0) << outer;

  TPaveText text_inner(100.0, 0.0, 140.0, 20.0);
  text_inner.SetX1(100.0);
  text_inner.SetX2(140.0);
  text_inner.SetY1(0.0);
  text_inner.SetY2(20.0);
  text_inner.SetBorderSize(0);
  text_inner.SetShadowColor(0);
  text_inner.SetLineColor(0);
  text_inner.SetLineWidth(0);
  TPaveText text_outer(150.0, 30.0, 190.0, 50.0);
  text_outer.SetX1(150.0);
  text_outer.SetX2(190.0);
  text_outer.SetY1(30.0);
  text_outer.SetY2(50.0);
  text_outer.SetBorderSize(0);
  text_outer.SetShadowColor(0);
  text_outer.SetLineColor(0);
  text_outer.SetLineWidth(0);
  text_inner.AddText(oss_inner.str().c_str());
  text_outer.AddText(oss_outer.str().c_str());
  text_inner.SetFillStyle(0);
  text_inner.SetTextColor(kRed);
  text_outer.SetFillStyle(0);
  text_outer.SetTextColor(kGreen);
  graph.Draw("p");
  if(false && print_num){
    text_inner.Draw("same");
    text_outer.Draw("same");
  }
  c.Print(out.c_str());
}

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
  TChain signal250("signal_m=250_GeV", "signal_m=250_GeV");
  TChain signal350("signal_m=350_GeV", "signal_m=350_GeV");

  data.Add("reduced_trees/MET_*2012*1_SyncSkim.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_FullLept*1_SyncSkim.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_SemiLept*1_SyncSkim.root/reduced_tree");
  qcd.Add("reduced_trees/TTJets_Ha2dronic*1_SyncSkim.root/reduced_tree");
  qcd.Add("reduced_trees/BJets*1_SyncSkim.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*channel*1_SyncSkim.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*JetsToLNu_Tune*1_SyncSkim.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/ZJetsToNuNu*1_SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/WH*1_SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/WW*1_SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/WZ*1_SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/ZZ*1_SyncSkim.root/reduced_tree");
  signal250.Add("reduced_trees/SMS-TChiHH_2b2b_2J_mChargino-250_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71_SyncSkim.root/reduced_tree");
  signal350.Add("reduced_trees/SMS-TChiHH_2b2b_2J_mChargino-350_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71_SyncSkim.root/reduced_tree");
  std::vector<TChain*> chains(0);
  chains.push_back(&data);
  chains.push_back(&qcd);
  chains.push_back(&diboson);
  chains.push_back(&single_t_or_boson);
  chains.push_back(&ttbar);
  chains.push_back(&signal250);
  chains.push_back(&signal350);

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
  unsigned short num_loose_electrons(0), num_loose_muons(0), num_loose_taus(0),
    num_loose_leptons(0);
  float min_delta_R(0.0), max_delta_R(0.0);
  float average_higgs_mass(0.0), higgs_mass_difference(0.0);
  float ht_jets(0.0), ht_jets_met(0.0), ht_jets_leps(0.0), ht_jets_met_leps(0.0);
  float full_weight(0.0);
  float top_pt(0.0);
  short chargino_mass(0);

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
  std::vector<TH1D> h_num_iso_tracks(0);
  std::vector<TH1D> h_num_electrons(0);
  std::vector<TH1D> h_num_muons(0);
  std::vector<TH1D> h_num_taus(0);
  std::vector<TH1D> h_num_leptons(0);
  std::vector<TH1D> h_num_loose_leptons(0);
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
  std::vector<TH1D> h_2b_avg(0);
  std::vector<TH1D> h_2b_diff(0);
  std::vector<TH1D> h_2b_sbin(0);
  std::vector<TH1D> h_3b_avg(0);
  std::vector<TH1D> h_3b_diff(0);
  std::vector<TH1D> h_3b_sbin(0);
  std::vector<TH1D> h_4b_avg(0);
  std::vector<TH1D> h_4b_diff(0);
  std::vector<TH1D> h_4b_sbin(0);
  std::vector<TH2D> h_diff_vs_avg(0);
  std::vector<TH2D> h_2b_diff_vs_avg(0);
  std::vector<TH2D> h_3b_diff_vs_avg(0);
  std::vector<TH2D> h_4b_diff_vs_avg(0);
  std::vector<TH1D> h_faildr_num_b_tagged_jets(0);
  std::vector<TH1D> h_faildr_met_sig(0);
  std::vector<TH1D> h_faildr_met(0);
  std::vector<TH1D> h_faildr_average_higgs_mass(0);
  std::vector<TH1D> h_faildr_higgs_mass_difference(0);
  std::vector<TH1D> h_passdr_num_b_tagged_jets(0);
  std::vector<TH1D> h_passdr_met_sig(0);
  std::vector<TH1D> h_passdr_met(0);
  std::vector<TH1D> h_passdr_average_higgs_mass(0);
  std::vector<TH1D> h_passdr_higgs_mass_difference(0);
  std::vector<TH1D> h_nodr_num_b_tagged_jets(0);
  std::vector<TH1D> h_nodr_met_sig(0);
  std::vector<TH1D> h_nodr_met(0);
  std::vector<TH1D> h_nodr_average_higgs_mass(0);
  std::vector<TH1D> h_nodr_higgs_mass_difference(0);
  std::vector<TH1D> h_maxdr_2b_sig(0);
  std::vector<TH1D> h_maxdr_3b_sig(0);
  std::vector<TH1D> h_maxdr_4b_sig(0);
  std::vector<TH1D> h_maxdr_2b_sb(0);
  std::vector<TH1D> h_maxdr_3b_sb(0);
  std::vector<TH1D> h_maxdr_4b_sb(0);
  std::vector<TH1D> h_mbb_2b_sb(0);
  std::vector<TH1D> h_mbb_3b_sb(0);
  std::vector<TH1D> h_mbb_4b_sb(0);
  std::vector<TH1D> h_mbb_2b_sig(0);
  std::vector<TH1D> h_mbb_3b_sig(0);
  std::vector<TH1D> h_mbb_4b_sig(0);
  std::vector<TH1D> h_mbb_nodr_2b_sb(0);
  std::vector<TH1D> h_mbb_nodr_3b_sb(0);
  std::vector<TH1D> h_mbb_nodr_4b_sb(0);
  std::vector<TH1D> h_mbb_nodr_2b_sig(0);
  std::vector<TH1D> h_mbb_nodr_3b_sig(0);
  std::vector<TH1D> h_mbb_nodr_4b_sig(0);
  std::vector<TH1D> h_mbb_1l_2b_sb(0);
  std::vector<TH1D> h_mbb_1l_3b_sb(0);
  std::vector<TH1D> h_mbb_1l_4b_sb(0);
  std::vector<TH1D> h_mbb_1l_2b_sig(0);
  std::vector<TH1D> h_mbb_1l_3b_sig(0);
  std::vector<TH1D> h_mbb_1l_4b_sig(0);
  std::vector<TH1D> h_mbb_1l_nodr_2b_sb(0);
  std::vector<TH1D> h_mbb_1l_nodr_3b_sb(0);
  std::vector<TH1D> h_mbb_1l_nodr_4b_sb(0);
  std::vector<TH1D> h_mbb_1l_nodr_2b_sig(0);
  std::vector<TH1D> h_mbb_1l_nodr_3b_sig(0);
  std::vector<TH1D> h_mbb_1l_nodr_4b_sig(0);

  std::vector<TH2D> h_dr_vs_avg(0);
  TGraph g_dr_vs_avg(0);
  TGraph g_dr_vs_avg_2b_sig(0);
  TGraph g_dr_vs_avg_3b_sig(0);
  TGraph g_dr_vs_avg_4b_sig(0);
  TGraph g_dr_vs_avg_2b_sb(0);
  TGraph g_dr_vs_avg_3b_sb(0);
  TGraph g_dr_vs_avg_4b_sb(0);

  TGraph g_faildr_diff_vs_avg_data(0);
  TGraph g_passdr_diff_vs_avg_data(0);
  TGraph g_nodr_diff_vs_avg_data(0);

  TGraph g_diff_vs_avg_data(0);
  TGraph g_2b_diff_vs_avg_data(0);
  TGraph g_3b_diff_vs_avg_data(0);
  TGraph g_4b_diff_vs_avg_data(0);

  for(unsigned chain_num(0); chain_num<chains.size(); ++chain_num){
    const double stupid_factor(chain_num==0?19307.0/5208.0:1.0);
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
    setup(chain, "num_iso_tracks", num_iso_tracks);
    setup(chain, "num_electrons", num_electrons);
    setup(chain, "num_muons", num_muons);
    setup(chain, "num_taus", num_taus);
    setup(chain, "num_leptons", num_leptons);
    setup(chain, "num_loose_electrons", num_loose_electrons);
    setup(chain, "num_loose_muons", num_loose_muons);
    setup(chain, "num_loose_taus", num_loose_taus);
    setup(chain, "num_loose_leptons", num_loose_leptons);
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
    h_num_leptons.push_back(TH1D(("h_num_leptons"+name).c_str(), "Num Veto Leptons (baseline);Num Leptons;Events/1", 16, -0.5, 15.5));
    h_num_loose_leptons.push_back(TH1D(("h_num_loose_leptons"+name).c_str(), "Num Loose Leptons (baseline);Num Leptons;Events/1", 16, -0.5, 15.5));
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

    h_2b_avg.push_back(TH1D(("h_2b_avg"+name).c_str(), "<m_{bb}> (baseline);<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_2b_diff.push_back(TH1D(("h_2b_diff"+name).c_str(), "#Delta m_{bb} (baseline);#Delta m_{bb} [GeV];Events/5 GeV", 24, 0.0, 120.0));
    h_2b_sbin.push_back(TH1D(("h_2b_sbin"+name).c_str(), "S_{MET} bin Counts (2b full selection);S_{MET} bin;Events", 4, 0.5, 4.5));
    h_3b_avg.push_back(TH1D(("h_3b_avg"+name).c_str(), "<m_{bb}> (baseline);<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_3b_diff.push_back(TH1D(("h_3b_diff"+name).c_str(), "#Delta m_{bb} (baseline);#Delta m_{bb} [GeV];Events/5 GeV", 24, 0.0, 120.0));
    h_3b_sbin.push_back(TH1D(("h_3b_sbin"+name).c_str(), "S_{MET} bin Counts (3b full selection);S_{MET} bin;Events", 4, 0.5, 4.5));
    h_4b_avg.push_back(TH1D(("h_4b_avg"+name).c_str(), "<m_{bb}> (baseline);<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_4b_diff.push_back(TH1D(("h_4b_diff"+name).c_str(), "#Delta m_{bb} (baseline);#Delta m_{bb} [GeV];Events/5 GeV", 24, 0.0, 120.0));
    h_4b_sbin.push_back(TH1D(("h_4b_sbin"+name).c_str(), "S_{MET} bin Counts (4b full selection);S_{MET} bin;Events", 4, 0.5, 4.5));

    h_diff_vs_avg.push_back(TH2D(("h_diff_vs_avg"+name).c_str(), ";<m_{bb}> [GeV];#Delta m_{bb} [GeV]", 50, 0.0, 250.0, 24, 0.0, 120.0));
    h_2b_diff_vs_avg.push_back(TH2D(("h_2b_diff_vs_avg"+name).c_str(), ";<m_{bb}> [GeV];#Delta m_{bb} [GeV]", 50, 0.0, 250.0, 24, 0.0, 120.0));
    h_3b_diff_vs_avg.push_back(TH2D(("h_3b_diff_vs_avg"+name).c_str(), ";<m_{bb}> [GeV];#Delta m_{bb} [GeV]", 50, 0.0, 250.0, 24, 0.0, 120.0));
    h_4b_diff_vs_avg.push_back(TH2D(("h_4b_diff_vs_avg"+name).c_str(), ";<m_{bb}> [GeV];#Delta m_{bb} [GeV]", 50, 0.0, 250.0, 24, 0.0, 120.0));

    h_faildr_num_b_tagged_jets.push_back(TH1D(("h_faildr_num_b_tagged_jets"+name).c_str(), "Num b-Tagged Jets (full sel.);Num b-Tagged jets;Events/1", 16, -0.5, 15.5));
    h_faildr_met_sig.push_back(TH1D(("h_faildr_met_sig"+name).c_str(), "S_{MET} (full sel.);S_{MET};Events/10", 40, 0.0, 400.0));
    h_faildr_met.push_back(TH1D(("h_faildr_met"+name).c_str(), "MET (full sel.);MET [GeV];Events/10 GeV", 40, 0.0, 400.0));
    h_faildr_average_higgs_mass.push_back(TH1D(("h_faildr_average_higgs_mass"+name).c_str(), "<m_{bb}> (full sel.);<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_faildr_higgs_mass_difference.push_back(TH1D(("h_faildr_higgs_mass_difference"+name).c_str(), "#Delta m_{bb} (full sel.);#Delta m_{bb} [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_passdr_num_b_tagged_jets.push_back(TH1D(("h_passdr_num_b_tagged_jets"+name).c_str(), "Num b-Tagged Jets (full sel.);Num b-Tagged jets;Events/1", 16, -0.5, 15.5));
    h_passdr_met_sig.push_back(TH1D(("h_passdr_met_sig"+name).c_str(), "S_{MET} (full sel.);S_{MET};Events/10", 40, 0.0, 400.0));
    h_passdr_met.push_back(TH1D(("h_passdr_met"+name).c_str(), "MET (full sel.);MET [GeV];Events/10 GeV", 40, 0.0, 400.0));
    h_passdr_average_higgs_mass.push_back(TH1D(("h_passdr_average_higgs_mass"+name).c_str(), "<m_{bb}> (full sel.);<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_passdr_higgs_mass_difference.push_back(TH1D(("h_passdr_higgs_mass_difference"+name).c_str(), "#Delta m_{bb} (full sel.);#Delta m_{bb} [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_nodr_num_b_tagged_jets.push_back(TH1D(("h_nodr_num_b_tagged_jets"+name).c_str(), "Num b-Tagged Jets (full sel.);Num b-Tagged jets;Events/1", 16, -0.5, 15.5));
    h_nodr_met_sig.push_back(TH1D(("h_nodr_met_sig"+name).c_str(), "S_{MET} (full sel.);S_{MET};Events/10", 40, 0.0, 400.0));
    h_nodr_met.push_back(TH1D(("h_nodr_met"+name).c_str(), "MET (full sel.);MET [GeV];Events/10 GeV", 40, 0.0, 400.0));
    h_nodr_average_higgs_mass.push_back(TH1D(("h_nodr_average_higgs_mass"+name).c_str(), "<m_{bb}> (full sel.);<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_nodr_higgs_mass_difference.push_back(TH1D(("h_nodr_higgs_mass_difference"+name).c_str(), "#Delta m_{bb} (full sel.);#Delta m_{bb} [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_maxdr_2b_sig.push_back(TH1D(("h_maxdr_2b_sig"+name).c_str(), "Max #Delta R (2b SIG);Max #Delta R;Events/0.2", 25, 0.0, 5.0));
    h_maxdr_3b_sig.push_back(TH1D(("h_maxdr_3b_sig"+name).c_str(), "Max #Delta R (3b SIG);Max #Delta R;Events/0.2", 25, 0.0, 5.0));
    h_maxdr_4b_sig.push_back(TH1D(("h_maxdr_4b_sig"+name).c_str(), "Max #Delta R (>3b SIG);Max #Delta R;Events/0.2", 25, 0.0, 5.0));
    h_maxdr_2b_sb.push_back(TH1D(("h_maxdr_2b_sb"+name).c_str(), "Max #Delta R (2b SB);Max #Delta R;Events/0.2", 25, 0.0, 5.0));
    h_maxdr_3b_sb.push_back(TH1D(("h_maxdr_3b_sb"+name).c_str(), "Max #Delta R (3b SB);Max #Delta R;Events/0.2", 25, 0.0, 5.0));
    h_maxdr_4b_sb.push_back(TH1D(("h_maxdr_4b_sb"+name).c_str(), "Max #Delta R (>3b SB);Max #Delta R;Events/0.2", 25, 0.0, 5.0));

    h_dr_vs_avg.push_back(TH2D(("h_dr_vs_avg"+name).c_str(), ";<m_{bb}> [GeV];Max #Delta R", 35, 0.0, 250.0, 35, 0.0, 5.0));

    h_mbb_2b_sb.push_back(TH1D(("h_mbb_2b_sb"+name).c_str(), "2b SB;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_3b_sb.push_back(TH1D(("h_mbb_3b_sb"+name).c_str(), "3b SB;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_4b_sb.push_back(TH1D(("h_mbb_4b_sb"+name).c_str(), "4b SB;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_2b_sig.push_back(TH1D(("h_mbb_2b_sig"+name).c_str(), "2b SIG;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_3b_sig.push_back(TH1D(("h_mbb_3b_sig"+name).c_str(), "3b SIG;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_4b_sig.push_back(TH1D(("h_mbb_4b_sig"+name).c_str(), "4b SIG;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_nodr_2b_sb.push_back(TH1D(("h_mbb_nodr_2b_sb"+name).c_str(), "2b SB no #Delta R cut;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_nodr_3b_sb.push_back(TH1D(("h_mbb_nodr_3b_sb"+name).c_str(), "3b SB no #Delta R cut;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_nodr_4b_sb.push_back(TH1D(("h_mbb_nodr_4b_sb"+name).c_str(), "4b SB no #Delta R cut;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_nodr_2b_sig.push_back(TH1D(("h_mbb_nodr_2b_sig"+name).c_str(), "2b SIG no #Delta R cut;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_nodr_3b_sig.push_back(TH1D(("h_mbb_nodr_3b_sig"+name).c_str(), "3b SIG no #Delta R cut;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_nodr_4b_sig.push_back(TH1D(("h_mbb_nodr_4b_sig"+name).c_str(), "4b SIG no #Delta R cut;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_2b_sb.push_back(TH1D(("h_mbb_1l_2b_sb"+name).c_str(), "2b SB 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_3b_sb.push_back(TH1D(("h_mbb_1l_3b_sb"+name).c_str(), "3b SB 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_4b_sb.push_back(TH1D(("h_mbb_1l_4b_sb"+name).c_str(), "4b SB 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_2b_sig.push_back(TH1D(("h_mbb_1l_2b_sig"+name).c_str(), "2b SIG 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_3b_sig.push_back(TH1D(("h_mbb_1l_3b_sig"+name).c_str(), "3b SIG 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_4b_sig.push_back(TH1D(("h_mbb_1l_4b_sig"+name).c_str(), "4b SIG 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_nodr_2b_sb.push_back(TH1D(("h_mbb_1l_nodr_2b_sb"+name).c_str(), "2b SB no #Delta R cut, 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_nodr_3b_sb.push_back(TH1D(("h_mbb_1l_nodr_3b_sb"+name).c_str(), "3b SB no #Delta R cut, 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_nodr_4b_sb.push_back(TH1D(("h_mbb_1l_nodr_4b_sb"+name).c_str(), "4b SB no #Delta R cut, 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_nodr_2b_sig.push_back(TH1D(("h_mbb_1l_nodr_2b_sig"+name).c_str(), "2b SIG no #Delta R cut, 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_nodr_3b_sig.push_back(TH1D(("h_mbb_1l_nodr_3b_sig"+name).c_str(), "3b SIG no #Delta R cut, 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));
    h_mbb_1l_nodr_4b_sig.push_back(TH1D(("h_mbb_1l_nodr_4b_sig"+name).c_str(), "4b SIG no #Delta R cut, 1l;<m_{bb}> [GeV];Events/10 GeV", 25, 0.0, 250.0));

    const int num_events(chain.GetEntries());
    Timer timer(num_events);
    timer.Start();
    for(int event(0); event<num_events; ++event){
      if(event%(1u<<16u)==0){
        timer.PrintRemainingTime();
      }
      chain.GetEntry(event);
      if(chain_num==chains.size()-1) full_weight*=0.004187617/1.33586290758103132e-05;
      if(chain_num==chains.size()-2) full_weight*=0.009816715/3.10530595015734434e-05;

      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
	if(num_b_tagged_jets==2){
	  if(higgs_mass_difference>30.0){
	    h_mbb_2b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_2b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets==3){
	  if(higgs_mass_difference>30.0){
	    h_mbb_3b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_3b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets>=4){
	  if(higgs_mass_difference>30.0){
	    h_mbb_4b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_4b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut){
	if(num_b_tagged_jets==2){
	  if(higgs_mass_difference>30.0){
	    h_mbb_nodr_2b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_nodr_2b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets==3){
	  if(higgs_mass_difference>30.0){
	    h_mbb_nodr_3b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_nodr_3b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets>=4){
	  if(higgs_mass_difference>30.0){
	    h_mbb_nodr_4b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_nodr_4b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && (num_loose_electrons+num_loose_muons==1 && num_taus==0) && passesDRCut){
	if(num_b_tagged_jets==2){
	  if(higgs_mass_difference>30.0){
	    h_mbb_1l_2b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_1l_2b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets==3){
	  if(higgs_mass_difference>30.0){
	    h_mbb_1l_3b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_1l_3b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets>=4){
	  if(higgs_mass_difference>30.0){
	    h_mbb_1l_4b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_1l_4b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && (num_loose_electrons+num_loose_muons==1 && num_taus==0)){
	if(num_b_tagged_jets==2){
	  if(higgs_mass_difference>30.0){
	    h_mbb_1l_nodr_2b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_1l_nodr_2b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets==3){
	  if(higgs_mass_difference>30.0){
	    h_mbb_1l_nodr_3b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_1l_nodr_3b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}else if(num_b_tagged_jets>=4){
	  if(higgs_mass_difference>30.0){
	    h_mbb_1l_nodr_4b_sb.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }else if(higgs_mass_difference<20.0){
	    h_mbb_1l_nodr_4b_sig.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	}
      }

      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut){
	if(num_b_tagged_jets==2){
	  if(higgs_mass_signal_region){
	    h_maxdr_2b_sig.at(chain_num).Fill(max_delta_R, full_weight);
	  }else if(higgs_mass_sideband){
	    h_maxdr_2b_sb.at(chain_num).Fill(max_delta_R, full_weight);
	  }
	}else if(num_b_tagged_jets==3){
	  if(higgs_mass_signal_region){
	    h_maxdr_3b_sig.at(chain_num).Fill(max_delta_R, full_weight);
	  }else if(higgs_mass_sideband){
	    h_maxdr_3b_sb.at(chain_num).Fill(max_delta_R, full_weight);
	  }
	}else if(num_b_tagged_jets>=4){
	  if(higgs_mass_signal_region){
	    h_maxdr_4b_sig.at(chain_num).Fill(max_delta_R, full_weight);
	  }else if(higgs_mass_sideband){
	    h_maxdr_4b_sb.at(chain_num).Fill(max_delta_R, full_weight);
	  }
	}

	h_dr_vs_avg.at(chain_num).Fill(average_higgs_mass, max_delta_R, full_weight);
	if(chain_num==0){
	  add_point(g_dr_vs_avg, average_higgs_mass, max_delta_R);
	  if(num_b_tagged_jets==2){
	    if(higgs_mass_signal_region){
	      add_point(g_dr_vs_avg_2b_sig, average_higgs_mass, max_delta_R);
	    }else if(higgs_mass_sideband){
	      add_point(g_dr_vs_avg_2b_sb, average_higgs_mass, max_delta_R);
	    }
	  }else if(num_b_tagged_jets==3){
	    if(higgs_mass_signal_region){
	      add_point(g_dr_vs_avg_3b_sig, average_higgs_mass, max_delta_R);
	    }else if(higgs_mass_sideband){
	      add_point(g_dr_vs_avg_3b_sb, average_higgs_mass, max_delta_R);
	    }
	  }else if(num_b_tagged_jets>=4){
	    if(higgs_mass_signal_region){
	      add_point(g_dr_vs_avg_4b_sig, average_higgs_mass, max_delta_R);
	    }else if(higgs_mass_sideband){
	      add_point(g_dr_vs_avg_4b_sb, average_higgs_mass, max_delta_R);
	    }
	  }
	}
      }

      if(passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut && passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesDRCut){
	unsigned short sbin(0);
	if(met_sig<50.0){
	  sbin=1;
	}else if(met_sig<100.0){
	  sbin=2;
	}else if(met_sig<150.0){
	  sbin=3;
	}else{
	  sbin=4;
	}
	if(num_b_tagged_jets==2){
	  h_2b_avg.at(chain_num).Fill(average_higgs_mass, full_weight);
	  h_2b_diff.at(chain_num).Fill(higgs_mass_difference, full_weight);
	  h_2b_sbin.at(chain_num).Fill(sbin, full_weight);
	  h_2b_diff_vs_avg.at(chain_num).Fill(average_higgs_mass, higgs_mass_difference, full_weight);
	  if(chain_num==0){
	    g_2b_diff_vs_avg_data.SetPoint(g_2b_diff_vs_avg_data.GetN(),
					   average_higgs_mass,
					   higgs_mass_difference);
	  }
	}else if(num_b_tagged_jets==3){
	  h_3b_avg.at(chain_num).Fill(average_higgs_mass, full_weight);
	  h_3b_diff.at(chain_num).Fill(higgs_mass_difference, full_weight);
	  h_3b_sbin.at(chain_num).Fill(sbin, full_weight);
	  h_3b_diff_vs_avg.at(chain_num).Fill(average_higgs_mass, higgs_mass_difference, full_weight);
	  if(chain_num==0){
	    g_3b_diff_vs_avg_data.SetPoint(g_3b_diff_vs_avg_data.GetN(),
					   average_higgs_mass,
					   higgs_mass_difference);
	  }
	}else if(num_b_tagged_jets>=4){
	  h_4b_avg.at(chain_num).Fill(average_higgs_mass, full_weight);
	  h_4b_diff.at(chain_num).Fill(higgs_mass_difference, full_weight);
	  h_4b_sbin.at(chain_num).Fill(sbin, full_weight);
	  h_4b_diff_vs_avg.at(chain_num).Fill(average_higgs_mass, higgs_mass_difference, full_weight);
	  if(chain_num==0){
	    g_4b_diff_vs_avg_data.SetPoint(g_4b_diff_vs_avg_data.GetN(),
					   average_higgs_mass,
					   higgs_mass_difference);
	  }
	}
	h_diff_vs_avg.at(chain_num).Fill(average_higgs_mass, higgs_mass_difference, full_weight);
	if(chain_num==0){
	  g_diff_vs_avg_data.SetPoint(g_diff_vs_avg_data.GetN(),
				      average_higgs_mass,
				      higgs_mass_difference);
	}
      }
    
      if(passesJSONCut && passesPVCut && second_highest_jet_pt>70.0
	 /*&& passes2CSVTCut*/ && !passesMETSig30Cut && passesMETCleaningCut
	 /*&& passesTriggerCut*/ && passesNumJetsCut && !passesMinDeltaPhiCut
	 && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut
	 && num_b_tagged_jets<3 && met>160.0 && passesQCDTriggerCut){
	h_qcd_control_highest_csv.at(chain_num).Fill(highest_csv, full_weight*stupid_factor);
	h_qcd_control_second_highest_csv.at(chain_num).Fill(second_highest_csv, full_weight*stupid_factor);
	h_qcd_control_third_highest_csv.at(chain_num).Fill(third_highest_csv, full_weight*stupid_factor);
	h_qcd_control_fourth_highest_csv.at(chain_num).Fill(fourth_highest_csv, full_weight*stupid_factor);
	h_qcd_control_fifth_highest_csv.at(chain_num).Fill(fifth_highest_csv, full_weight*stupid_factor);
      }
      if(passesJSONCut && passesPVCut && second_highest_jet_pt>70.0
	 /*&& passes2CSVTCut && !passesMETSig30Cut*/ && passesMETCleaningCut
	 /*&& passesTriggerCut*/ && passesNumJetsCut && !passesMinDeltaPhiCut
	 && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut
	 && num_b_tagged_jets<3 && met>160.0 && passesQCDTriggerCut){
	h_qcd_control_met_sig.at(chain_num).Fill(met_sig, full_weight*stupid_factor);
      }
      if(passesJSONCut && passesPVCut && second_highest_jet_pt>70.0
	 /*&& passes2CSVTCut*/ && !passesMETSig30Cut && passesMETCleaningCut
	 /*&& passesTriggerCut*/ && passesNumJetsCut && !passesMinDeltaPhiCut
	 && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut
	 && num_b_tagged_jets<3 /*&& met>160.0*/ && passesQCDTriggerCut){
	h_qcd_control_met.at(chain_num).Fill(met, full_weight*stupid_factor);
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
	 && passes2CSVTCut /*&& passesMETSig30Cut*/ && passesMETCleaningCut
	 && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
	 /*&& passesLeptonVetoCut && passesIsoTrackVetoCut*/ && passesDRCut
	 && num_leptons==1 && num_taus==0){
	h_1l_met_sig.at(chain_num).Fill(met_sig, full_weight*stupid_factor);
	h_1l_met.at(chain_num).Fill(met, full_weight*stupid_factor);
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
	if(num_taus==0){
	  h_num_leptons.at(chain_num).Fill(num_electrons+num_muons, full_weight);
	  h_num_loose_leptons.at(chain_num).Fill(num_loose_electrons+num_loose_muons, full_weight);
	}
      }
      if(passesJSONCut && passesPVCut && passesJet2PtCut
	 && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut
	 && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
	 && passesLeptonVetoCut && passesIsoTrackVetoCut /*&& passesDRCut*/){
	h_max_delta_R.at(chain_num).Fill(max_delta_R, full_weight);
      }

      if(passesJSONCut && passesPVCut && passesJet2PtCut
	 && passes2CSVTCut && passesMETCleaningCut
	 && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut
	 && passesLeptonVetoCut && passesIsoTrackVetoCut){
	if(!passesDRCut){
	  if(passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesMETSig30Cut){
	    h_faildr_num_b_tagged_jets.at(chain_num).Fill(num_b_tagged_jets, full_weight);
	  }
	  if(passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesBTaggingCut){
	    h_faildr_met_sig.at(chain_num).Fill(met_sig, full_weight);
	    h_faildr_met.at(chain_num).Fill(met, full_weight);
	  }
	  if(passesBTaggingCut && passesHiggsMassDiffCut && passesMETSig30Cut){
	    h_faildr_average_higgs_mass.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	  if(passesHiggsAvgMassCut && passesBTaggingCut && passesMETSig30Cut){
	    h_faildr_higgs_mass_difference.at(chain_num).Fill(higgs_mass_difference, full_weight);
	  }
	  if(passesMETSig30Cut && passesBTaggingCut && chain_num==0){
	    g_faildr_diff_vs_avg_data.SetPoint(g_faildr_diff_vs_avg_data.GetN(),
					       average_higgs_mass,
					       higgs_mass_difference);

	  }
	}else{
	  if(passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesMETSig30Cut){
	    h_passdr_num_b_tagged_jets.at(chain_num).Fill(num_b_tagged_jets, full_weight);
	  }
	  if(passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesBTaggingCut){
	    h_passdr_met_sig.at(chain_num).Fill(met_sig, full_weight);
	    h_passdr_met.at(chain_num).Fill(met, full_weight);
	  }
	  if(passesBTaggingCut && passesHiggsMassDiffCut && passesMETSig30Cut){
	    h_passdr_average_higgs_mass.at(chain_num).Fill(average_higgs_mass, full_weight);
	  }
	  if(passesHiggsAvgMassCut && passesBTaggingCut && passesMETSig30Cut){
	    h_passdr_higgs_mass_difference.at(chain_num).Fill(higgs_mass_difference, full_weight);
	  }
	  if(passesMETSig30Cut && num_b_tagged_jets>=4 && chain_num==0){
	    g_passdr_diff_vs_avg_data.SetPoint(g_passdr_diff_vs_avg_data.GetN(),
					       average_higgs_mass,
					       higgs_mass_difference);

	  }
	}
	if(passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesMETSig30Cut){
	  h_nodr_num_b_tagged_jets.at(chain_num).Fill(num_b_tagged_jets, full_weight);
	}
	if(passesHiggsAvgMassCut && passesHiggsMassDiffCut && passesBTaggingCut){
	  h_nodr_met_sig.at(chain_num).Fill(met_sig, full_weight);
	  h_nodr_met.at(chain_num).Fill(met, full_weight);
	}
	if(passesBTaggingCut && passesHiggsMassDiffCut && passesMETSig30Cut){
	  h_nodr_average_higgs_mass.at(chain_num).Fill(average_higgs_mass, full_weight);
	}
	if(passesHiggsAvgMassCut && passesBTaggingCut && passesMETSig30Cut){
	  h_nodr_higgs_mass_difference.at(chain_num).Fill(higgs_mass_difference, full_weight);
	}
	  if(passesMETSig30Cut && passesBTaggingCut && chain_num==0){
	    g_nodr_diff_vs_avg_data.SetPoint(g_nodr_diff_vs_avg_data.GetN(),
					       average_higgs_mass,
					       higgs_mass_difference);

	  }
      }
      timer.Iterate();
    }
  }

  std::vector<std::string> sub_names(names.begin()+1, names.end());

  plotter plot;
  plot.set_mc_names(sub_names);
  plot_data_mc(plot, h_mbb_2b_sb, "nm1/mbb_2b_sb.pdf");
  plot_data_mc(plot, h_mbb_3b_sb, "nm1/mbb_3b_sb.pdf");
  plot_data_mc(plot, h_mbb_4b_sb, "nm1/mbb_4b_sb.pdf");
  plot_data_mc(plot, h_mbb_2b_sig, "nm1/mbb_2b_sig.pdf");
  plot_data_mc(plot, h_mbb_3b_sig, "nm1/mbb_3b_sig.pdf");
  plot_data_mc(plot, h_mbb_4b_sig, "nm1/mbb_4b_sig.pdf");
  plot_data_mc(plot, h_mbb_nodr_2b_sb, "nm1/mbb_nodr_2b_sb.pdf");
  plot_data_mc(plot, h_mbb_nodr_3b_sb, "nm1/mbb_nodr_3b_sb.pdf");
  plot_data_mc(plot, h_mbb_nodr_4b_sb, "nm1/mbb_nodr_4b_sb.pdf");
  plot_data_mc(plot, h_mbb_nodr_2b_sig, "nm1/mbb_nodr_2b_sig.pdf");
  plot_data_mc(plot, h_mbb_nodr_3b_sig, "nm1/mbb_nodr_3b_sig.pdf");
  plot_data_mc(plot, h_mbb_nodr_4b_sig, "nm1/mbb_nodr_4b_sig.pdf");
  plot_data_mc(plot, h_mbb_1l_2b_sb, "nm1/mbb_1l_2b_sb.pdf");
  plot_data_mc(plot, h_mbb_1l_3b_sb, "nm1/mbb_1l_3b_sb.pdf");
  plot_data_mc(plot, h_mbb_1l_4b_sb, "nm1/mbb_1l_4b_sb.pdf");
  plot_data_mc(plot, h_mbb_1l_2b_sig, "nm1/mbb_1l_2b_sig.pdf");
  plot_data_mc(plot, h_mbb_1l_3b_sig, "nm1/mbb_1l_3b_sig.pdf");
  plot_data_mc(plot, h_mbb_1l_4b_sig, "nm1/mbb_1l_4b_sig.pdf");
  plot_data_mc(plot, h_mbb_1l_nodr_2b_sb, "nm1/mbb_1l_nodr_2b_sb.pdf");
  plot_data_mc(plot, h_mbb_1l_nodr_3b_sb, "nm1/mbb_1l_nodr_3b_sb.pdf");
  plot_data_mc(plot, h_mbb_1l_nodr_4b_sb, "nm1/mbb_1l_nodr_4b_sb.pdf");
  plot_data_mc(plot, h_mbb_1l_nodr_2b_sig, "nm1/mbb_1l_nodr_2b_sig.pdf");
  plot_data_mc(plot, h_mbb_1l_nodr_3b_sig, "nm1/mbb_1l_nodr_3b_sig.pdf");
  plot_data_mc(plot, h_mbb_1l_nodr_4b_sig, "nm1/mbb_1l_nodr_4b_sig.pdf");
  plot_data_mc(plot, h_pu_true_num_interactions, "nm1/pu_true_num_interactions.pdf");
  plot_data_mc(plot, h_num_primary_vertices, "nm1/num_primary_vertices.pdf");
  plot_data_mc(plot, h_highest_jet_pt, "nm1/highest_jet_pt.pdf");
  plot_data_mc(plot, h_second_highest_jet_pt, "nm1/second_highest_jet_pt.pdf");
  plot_data_mc(plot, h_third_highest_jet_pt, "nm1/third_highest_jet_pt.pdf");
  plot_data_mc(plot, h_fourth_highest_jet_pt, "nm1/fourth_highest_jet_pt.pdf");
  plot_data_mc(plot, h_fifth_highest_jet_pt, "nm1/fifth_highest_jet_pt.pdf");
  plot_data_mc(plot, h_highest_csv, "nm1/highest_csv.pdf");
  plot_data_mc(plot, h_second_highest_csv, "nm1/second_highest_csv.pdf");
  plot_data_mc(plot, h_third_highest_csv, "nm1/third_highest_csv.pdf");
  plot_data_mc(plot, h_fourth_highest_csv, "nm1/fourth_highest_csv.pdf");
  plot_data_mc(plot, h_fifth_highest_csv, "nm1/fifth_highest_csv.pdf");
  plot_data_mc(plot, h_met_sig, "nm1/met_sig.pdf");
  plot_data_mc(plot, h_met, "nm1/met.pdf");
  plot_data_mc(plot, h_num_jets, "nm1/num_jets.pdf");
  plot_data_mc(plot, h_num_b_tagged_jets, "nm1/num_b_tagged_jets.pdf");
  plot_data_mc(plot, h_min_delta_phi, "nm1/min_delta_phi.pdf");
  plot_data_mc(plot, h_num_electrons, "nm1/num_electrons.pdf");
  plot_data_mc(plot, h_num_muons, "nm1/num_muons.pdf");
  plot_data_mc(plot, h_num_taus, "nm1/num_taus.pdf");
  plot_data_mc(plot, h_num_iso_tracks, "nm1/num_iso_tracks.pdf");
  plot_data_mc(plot, h_num_leptons, "nm1/num_leptons.pdf");
  plot_data_mc(plot, h_num_loose_leptons, "nm1/num_loose_leptons.pdf");
  plot_data_mc(plot, h_min_delta_R, "nm1/min_delta_R.pdf");
  plot_data_mc(plot, h_max_delta_R, "nm1/max_delta_R.pdf");
  plot_data_mc(plot, h_average_higgs_mass, "nm1/average_higgs_mass.pdf");
  plot_data_mc(plot, h_higgs_mass_difference, "nm1/higgs_mass_difference.pdf");
  plot_data_mc(plot, h_ht_jets, "nm1/ht_jets.pdf");
  plot_data_mc(plot, h_ht_jets_met, "nm1/ht_jets_met.pdf");
  plot_data_mc(plot, h_ht_jets_leps, "nm1/ht_jets_leps.pdf");
  plot_data_mc(plot, h_ht_jets_met_leps, "nm1/ht_jets_met_leps.pdf");
  plot_data_mc(plot, h_top_pt, "nm1/top_pt.pdf");
  plot_data_mc(plot, h_qcd_control_highest_csv, "nm1/qcd_control_highest_csv.pdf");
  plot_data_mc(plot, h_qcd_control_second_highest_csv, "nm1/qcd_control_second_highest_csv.pdf");
  plot_data_mc(plot, h_qcd_control_third_highest_csv, "nm1/qcd_control_third_highest_csv.pdf");
  plot_data_mc(plot, h_qcd_control_fourth_highest_csv, "nm1/qcd_control_fourth_highest_csv.pdf");
  plot_data_mc(plot, h_qcd_control_fifth_highest_csv, "nm1/qcd_control_fifth_highest_csv.pdf");
  plot_data_mc(plot, h_qcd_control_met_sig, "nm1/qcd_control_met_sig.pdf");
  plot_data_mc(plot, h_qcd_control_met, "nm1/qcd_control_met.pdf");
  plot_data_mc(plot, h_1l_highest_csv, "nm1/1l_highest_csv.pdf");
  plot_data_mc(plot, h_1l_second_highest_csv, "nm1/1l_second_highest_csv.pdf");
  plot_data_mc(plot, h_1l_third_highest_csv, "nm1/1l_third_highest_csv.pdf");
  plot_data_mc(plot, h_1l_fourth_highest_csv, "nm1/1l_fourth_highest_csv.pdf");
  plot_data_mc(plot, h_1l_fifth_highest_csv, "nm1/1l_fifth_highest_csv.pdf");
  plot_data_mc(plot, h_1l_met_sig, "nm1/1l_met_sig.pdf");
  plot_data_mc(plot, h_1l_met, "nm1/1l_met.pdf");
  plot_data_mc_sig(plot, h_2b_sbin, "nm1/2b_sbin.pdf");
  plot_data_mc_sig(plot, h_3b_sbin, "nm1/3b_sbin.pdf");
  plot_data_mc_sig(plot, h_4b_sbin, "nm1/4b_sbin.pdf");
  plot_data_mc(plot, h_faildr_num_b_tagged_jets, "nm1/faildr_num_b_tagged_jets.pdf");
  plot_data_mc(plot, h_faildr_met_sig, "nm1/faildr_met_sig.pdf");
  plot_data_mc(plot, h_faildr_met, "nm1/faildr_met.pdf");
  plot_data_mc(plot, h_faildr_average_higgs_mass, "nm1/faildr_average_higgs_mass.pdf");
  plot_data_mc(plot, h_faildr_higgs_mass_difference, "nm1/faildr_higgs_mass_difference.pdf");
  plot_data_mc(plot, h_passdr_num_b_tagged_jets, "nm1/passdr_num_b_tagged_jets.pdf");
  plot_data_mc(plot, h_passdr_met_sig, "nm1/passdr_met_sig.pdf");
  plot_data_mc(plot, h_passdr_met, "nm1/passdr_met.pdf");
  plot_data_mc(plot, h_passdr_average_higgs_mass, "nm1/passdr_average_higgs_mass.pdf");
  plot_data_mc(plot, h_passdr_higgs_mass_difference, "nm1/passdr_higgs_mass_difference.pdf");
  plot_data_mc(plot, h_nodr_num_b_tagged_jets, "nm1/nodr_num_b_tagged_jets.pdf");
  plot_data_mc(plot, h_nodr_met_sig, "nm1/nodr_met_sig.pdf");
  plot_data_mc(plot, h_nodr_met, "nm1/nodr_met.pdf");
  plot_data_mc(plot, h_nodr_average_higgs_mass, "nm1/nodr_average_higgs_mass.pdf");
  plot_data_mc(plot, h_nodr_higgs_mass_difference, "nm1/nodr_higgs_mass_difference.pdf");
  plot_data_mc(plot, h_maxdr_2b_sig, "nm1/maxdr_2b_sig.pdf");
  plot_data_mc(plot, h_maxdr_3b_sig, "nm1/maxdr_3b_sig.pdf");
  plot_data_mc(plot, h_maxdr_4b_sig, "nm1/maxdr_4b_sig.pdf");
  plot_data_mc(plot, h_maxdr_2b_sb, "nm1/maxdr_2b_sb.pdf");
  plot_data_mc(plot, h_maxdr_3b_sb, "nm1/maxdr_3b_sb.pdf");
  plot_data_mc(plot, h_maxdr_4b_sb, "nm1/maxdr_4b_sb.pdf");
  
  for(unsigned chain_num(0); chain_num<chains.size(); ++chain_num){
    TCanvas cc;
    h_diff_vs_avg.at(chain_num).Draw("SCAT");
    cc.Print(("nm1/diff_vs_avg"+names.at(chain_num)+".pdf").c_str());
  }

  TH1D h_2b_avg_sm(std::accumulate(h_2b_avg.begin()+2, h_2b_avg.end()-2, h_2b_avg.at(1)));
  TH1D h_3b_avg_sm(std::accumulate(h_3b_avg.begin()+2, h_3b_avg.end()-2, h_3b_avg.at(1)));
  TH1D h_4b_avg_sm(std::accumulate(h_4b_avg.begin()+2, h_4b_avg.end()-2, h_4b_avg.at(1)));
  TH1D h_2b_diff_sm(std::accumulate(h_2b_diff.begin()+2, h_2b_diff.end()-2, h_2b_diff.at(1)));
  TH1D h_3b_diff_sm(std::accumulate(h_3b_diff.begin()+2, h_3b_diff.end()-2, h_3b_diff.at(1)));
  TH1D h_4b_diff_sm(std::accumulate(h_4b_diff.begin()+2, h_4b_diff.end()-2, h_4b_diff.at(1)));
  compare_histos(h_2b_avg_sm, h_3b_avg_sm, h_4b_avg_sm, "nm1/avg_compare.pdf");
  compare_histos(h_2b_avg_sm, h_3b_avg_sm, h_4b_avg_sm, "nm1/avg_compare_norm.pdf", true);
  compare_histos(h_2b_diff_sm, h_3b_diff_sm, h_4b_diff_sm, "nm1/diff_compare.pdf");
  compare_histos(h_2b_diff_sm, h_3b_diff_sm, h_4b_diff_sm, "nm1/diff_compare_norm.pdf", true);

  compare_histos(h_2b_avg.at(1), h_3b_avg.at(1), h_4b_avg.at(1), "nm1/avg_compare_qcd.pdf");
  compare_histos(h_2b_avg.at(1), h_3b_avg.at(1), h_4b_avg.at(1), "nm1/avg_compare_qcd_norm.pdf", true);
  compare_histos(h_2b_diff.at(1), h_3b_diff.at(1), h_4b_diff.at(1), "nm1/diff_compare_qcd.pdf");
  compare_histos(h_2b_diff.at(1), h_3b_diff.at(1), h_4b_diff.at(1), "nm1/diff_compare_qcd_norm.pdf", true);

  TH2D h_diff_vs_avg_sm(std::accumulate(h_diff_vs_avg.begin()+2, h_diff_vs_avg.end()-2, h_diff_vs_avg.at(1)));
  TH2D h_2b_diff_vs_avg_sm(std::accumulate(h_2b_diff_vs_avg.begin()+2, h_2b_diff_vs_avg.end()-2, h_2b_diff_vs_avg.at(1)));
  TH2D h_3b_diff_vs_avg_sm(std::accumulate(h_3b_diff_vs_avg.begin()+2, h_3b_diff_vs_avg.end()-2, h_3b_diff_vs_avg.at(1)));
  TH2D h_4b_diff_vs_avg_sm(std::accumulate(h_4b_diff_vs_avg.begin()+2, h_4b_diff_vs_avg.end()-2, h_4b_diff_vs_avg.at(1)));

  TH2D h_dr_vs_avg_sm(std::accumulate(h_dr_vs_avg.begin()+2, h_dr_vs_avg.end()-2, h_dr_vs_avg.at(1)));
  scat_plot(h_dr_vs_avg_sm, "nm1/h_dr_vs_avg_sm.pdf");

  sig_sb_scat(h_diff_vs_avg_sm, "nm1/diff_vs_avg_sm.pdf");
  sig_sb_scat(h_2b_diff_vs_avg_sm, "nm1/2b_diff_vs_avg_sm.pdf");
  sig_sb_scat(h_3b_diff_vs_avg_sm, "nm1/3b_diff_vs_avg_sm.pdf");
  sig_sb_scat(h_4b_diff_vs_avg_sm, "nm1/4b_diff_vs_avg_sm.pdf");
  sig_sb_scat(h_2b_diff_vs_avg_sm, "nm1/2b_diff_vs_avg_sig250.pdf");
  sig_sb_scat(h_3b_diff_vs_avg_sm, "nm1/3b_diff_vs_avg_sig250.pdf");
  sig_sb_scat(*(h_diff_vs_avg.end()-2), "nm1/diff_vs_avg_sig250.pdf");
  sig_sb_scat(*(h_2b_diff_vs_avg.end()-2), "nm1/2b_diff_vs_avg_sig250.pdf");
  sig_sb_scat(*(h_3b_diff_vs_avg.end()-2), "nm1/3b_diff_vs_avg_sig250.pdf");
  sig_sb_scat(*(h_4b_diff_vs_avg.end()-2), "nm1/4b_diff_vs_avg_sig250.pdf");
  sig_sb_scat(*(h_diff_vs_avg.end()-1), "nm1/diff_vs_avg_sig350.pdf");
  sig_sb_scat(*(h_2b_diff_vs_avg.end()-1), "nm1/2b_diff_vs_avg_sig350.pdf");
  sig_sb_scat(*(h_3b_diff_vs_avg.end()-1), "nm1/3b_diff_vs_avg_sig350.pdf");
  sig_sb_scat(*(h_4b_diff_vs_avg.end()-1), "nm1/4b_diff_vs_avg_sig350.pdf");
  sig_sb_scat(h_diff_vs_avg.at(1), "nm1/diff_vs_avg_qcd.pdf");
  sig_sb_scat(h_2b_diff_vs_avg.at(1), "nm1/2b_diff_vs_avg_qcd.pdf");
  sig_sb_scat(h_3b_diff_vs_avg.at(1), "nm1/3b_diff_vs_avg_qcd.pdf");
  sig_sb_scat(h_4b_diff_vs_avg.at(1), "nm1/4b_diff_vs_avg_qcd.pdf");
  sig_sb_scat(h_diff_vs_avg_sm, "nm1/diff_vs_avg_sm_count.pdf", true);
  sig_sb_scat(h_2b_diff_vs_avg_sm, "nm1/2b_diff_vs_avg_sm_count.pdf", true);
  sig_sb_scat(h_3b_diff_vs_avg_sm, "nm1/3b_diff_vs_avg_sm_count.pdf", true);
  sig_sb_scat(h_4b_diff_vs_avg_sm, "nm1/4b_diff_vs_avg_sm_count.pdf", true);
  sig_sb_scat(h_2b_diff_vs_avg_sm, "nm1/2b_diff_vs_avg_sig250_count.pdf", true);
  sig_sb_scat(h_3b_diff_vs_avg_sm, "nm1/3b_diff_vs_avg_sig250_count.pdf", true);
  sig_sb_scat(*(h_diff_vs_avg.end()-2), "nm1/diff_vs_avg_sig250_count.pdf", true);
  sig_sb_scat(*(h_2b_diff_vs_avg.end()-2), "nm1/2b_diff_vs_avg_sig250_count.pdf", true);
  sig_sb_scat(*(h_3b_diff_vs_avg.end()-2), "nm1/3b_diff_vs_avg_sig250_count.pdf", true);
  sig_sb_scat(*(h_4b_diff_vs_avg.end()-2), "nm1/4b_diff_vs_avg_sig250_count.pdf", true);
  sig_sb_scat(*(h_diff_vs_avg.end()-1), "nm1/diff_vs_avg_sig350_count.pdf", true);
  sig_sb_scat(*(h_2b_diff_vs_avg.end()-1), "nm1/2b_diff_vs_avg_sig350_count.pdf", true);
  sig_sb_scat(*(h_3b_diff_vs_avg.end()-1), "nm1/3b_diff_vs_avg_sig350_count.pdf", true);
  sig_sb_scat(*(h_4b_diff_vs_avg.end()-1), "nm1/4b_diff_vs_avg_sig350_count.pdf", true);
  sig_sb_scat(h_diff_vs_avg.at(1), "nm1/diff_vs_avg_qcd_count.pdf", true);
  sig_sb_scat(h_2b_diff_vs_avg.at(1), "nm1/2b_diff_vs_avg_qcd_count.pdf", true);
  sig_sb_scat(h_3b_diff_vs_avg.at(1), "nm1/3b_diff_vs_avg_qcd_count.pdf", true);
  sig_sb_scat(h_4b_diff_vs_avg.at(1), "nm1/4b_diff_vs_avg_qcd_count.pdf", true);
  sig_sb_scat_data(g_diff_vs_avg_data, "nm1/diff_vs_avg_data.pdf");
  sig_sb_scat_data(g_2b_diff_vs_avg_data, "nm1/2b_diff_vs_avg_data.pdf");
  sig_sb_scat_data(g_3b_diff_vs_avg_data, "nm1/3b_diff_vs_avg_data.pdf");
  sig_sb_scat_data(g_4b_diff_vs_avg_data, "nm1/4b_diff_vs_avg_data.pdf");
  sig_sb_scat_data(g_diff_vs_avg_data, "nm1/diff_vs_avg_data_count.pdf", true);
  sig_sb_scat_data(g_2b_diff_vs_avg_data, "nm1/2b_diff_vs_avg_data_count.pdf", true);
  sig_sb_scat_data(g_3b_diff_vs_avg_data, "nm1/3b_diff_vs_avg_data_count.pdf", true);
  sig_sb_scat_data(g_4b_diff_vs_avg_data, "nm1/4b_diff_vs_avg_data_count.pdf", true);

  sig_sb_scat_data(g_faildr_diff_vs_avg_data, "nm1/faildr_diff_vs_avg_data_count.pdf", true);
  sig_sb_scat_data(g_passdr_diff_vs_avg_data, "nm1/passdr_diff_vs_avg_data_count.pdf", true);
  sig_sb_scat_data(g_nodr_diff_vs_avg_data, "nm1/nodr_diff_vs_avg_data_count.pdf", true);

  scat_plot_data(g_dr_vs_avg, "nm1/dr_vs_avg.pdf");
  scat_plot_data(g_dr_vs_avg_2b_sig, "nm1/dr_vs_avg_2b_sig.pdf");
  scat_plot_data(g_dr_vs_avg_2b_sb, "nm1/dr_vs_avg_2b_sb.pdf");
  scat_plot_data(g_dr_vs_avg_3b_sig, "nm1/dr_vs_avg_3b_sig.pdf");
  scat_plot_data(g_dr_vs_avg_3b_sb, "nm1/dr_vs_avg_3b_sb.pdf");
  scat_plot_data(g_dr_vs_avg_4b_sig, "nm1/dr_vs_avg_4b_sig.pdf");
  scat_plot_data(g_dr_vs_avg_4b_sb, "nm1/dr_vs_avg_4b_sb.pdf");

  //TH1D tot_min_delta_phi(GetMCSum(h_min_delta_phi));
  TH1D tot_min_delta_phi(h_min_delta_phi.at(0));
  const int the_bin(tot_min_delta_phi.FindBin(0.3));
  const int num_bins(tot_min_delta_phi.GetNbinsX());
  const double low_integral(tot_min_delta_phi.Integral(0, the_bin-1));
  const double high_integral(tot_min_delta_phi.Integral(the_bin, num_bins+1));
  std::cout << "Low integral: " << low_integral << std::endl;
  std::cout << "High integral: " << high_integral << std::endl;

  std::cout << (h_4b_diff_vs_avg.end()-1)->Integral() << std::endl;
}
