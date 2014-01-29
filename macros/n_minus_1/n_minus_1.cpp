#include <fstream>
#include <vector>
#include <math.h> 
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "cut.hpp"
#include "my_style.hpp"

ofstream logEff;

TChain ttbar_ch("reduced_tree");
TChain qcd_ch("reduced_tree");
TChain znn_ch("reduced_tree");
TChain wt_ch("reduced_tree");
TChain other_ch("reduced_tree");
TChain data_obs_ch("reduced_tree");
TChain m350_ch("reduced_tree");

TTree* ttbar_;
TTree* qcd_;
TTree* znn_;
TTree* wt_;
TTree* other_;
TTree* data_obs_;
TTree* m350_;


void n_minus_1_plot(Cut* cut, TString title, int nbinsx, double xlow, double xup,
		    TString plotname="plot", TString plotdir="",
		    TString options="baselineCuts:plotSig:plotData:plotLog")
{

  // Parse options
  //bool printStat = options.Contains("printStat") && (!options.Contains("!printStat"));
  //bool printCard = options.Contains("printCard") && (!options.Contains("!printCard"));
  //bool plotUOflow = options.Contains("plotUOflow") && (!options.Contains("!plotUOflow")); // FIXME: not yet implemented
  bool baselineCuts = options.Contains("baselineCuts") && (!options.Contains("!baselineCuts"));
  bool plotSig = options.Contains("plotSig") && (!options.Contains("!plotSig"));
  bool plotData = options.Contains("plotData") && (!options.Contains("!plotData"));
  bool plotLog = options.Contains("plotLog") && (!options.Contains("!plotLog"));

  // Store baseline counts
  TH1D * httbar_b = new TH1D("ttbar_b" , title, nbinsx, xlow, xup);
  TH1D * hqcd_b = new TH1D("qcd_b" , title, nbinsx, xlow, xup);
  TH1D * hznn_b = new TH1D("znn_b" , title, nbinsx, xlow, xup);
  TH1D * hwt_b = new TH1D("wt_b" , title, nbinsx, xlow, xup);
  TH1D * hother_b = new TH1D("other_b" , title, nbinsx, xlow, xup);
  TH1D * hm350_b = new TH1D("m350_b" , title, nbinsx, xlow, xup);
  TCut baseline_;
  if (cut->GetVar().EqualTo("met_sig")) baseline_="full_weight";
  else baseline_ = "(passesMETSig30Cut)*full_weight";
  ttbar_->Project("ttbar_b",cut->GetVar(),baseline_);
  qcd_->Project("qcd_b",cut->GetVar(),baseline_);
  znn_->Project("znn_b",cut->GetVar(),baseline_);
  wt_->Project("wt_b",cut->GetVar(),baseline_);
  other_->Project("other_b",cut->GetVar(),baseline_);
  m350_->Project("m350_b",cut->GetVar(),baseline_);
  vector<double> baseline_counts;
  baseline_counts.push_back(httbar_b->Integral());
  baseline_counts.push_back(hqcd_b->Integral());
  baseline_counts.push_back(hznn_b->Integral());
  baseline_counts.push_back(hwt_b->Integral());
  baseline_counts.push_back(hother_b->Integral());
  baseline_counts.push_back(baseline_counts[0]+baseline_counts[1]+
			    baseline_counts[2]+baseline_counts[3]+
			    baseline_counts[4]);
  baseline_counts.push_back(hm350_b->Integral());
  httbar_b->Reset();
  hqcd_b->Reset();
  hznn_b->Reset();
  hwt_b->Reset();
  hother_b->Reset();
  hm350_b->Reset();
  //Do it agian, with the cut
  ttbar_->Project("ttbar_b",cut->GetVar(),Form("(%s)",cut->Get_Var_Cut().Data())*baseline_);
  qcd_->Project("qcd_b",cut->GetVar(),Form("(%s)",cut->Get_Var_Cut().Data())*baseline_);
  znn_->Project("znn_b",cut->GetVar(),Form("(%s)",cut->Get_Var_Cut().Data())*baseline_);
  wt_->Project("wt_b",cut->GetVar(),Form("(%s)",cut->Get_Var_Cut().Data())*baseline_);
  other_->Project("other_b",cut->GetVar(),Form("(%s)",cut->Get_Var_Cut().Data())*baseline_);
  m350_->Project("m350_b",cut->GetVar(),Form("(%s)",cut->Get_Var_Cut().Data())*baseline_);
  vector<double> baseline_counts_a;
  baseline_counts_a.push_back(httbar_b->Integral());
  baseline_counts_a.push_back(hqcd_b->Integral());
  baseline_counts_a.push_back(hznn_b->Integral());
  baseline_counts_a.push_back(hwt_b->Integral());
  baseline_counts_a.push_back(hother_b->Integral());
  baseline_counts_a.push_back(baseline_counts_a[0]+baseline_counts_a[1]+
			      baseline_counts_a[2]+baseline_counts_a[3]+
			      baseline_counts_a[4]);
  baseline_counts_a.push_back(hm350_b->Integral());

  // Setup style
  cout << "Setting tdr style."  << endl;
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  setTDRStyle(tdrStyle);
  tdrStyle->cd();

  // Book histograms
  TH1D * httbar = new TH1D("ttbar" , title, nbinsx, xlow, xup);
  TH1D * hqcd = new TH1D("qcd" , title, nbinsx, xlow, xup);
  TH1D * hznn = new TH1D("znn" , title, nbinsx, xlow, xup);
  TH1D * hwt = new TH1D("wt" , title, nbinsx, xlow, xup);
  TH1D * hother = new TH1D("other" , title, nbinsx, xlow, xup);
  TH1D * hmc_exp = new TH1D("mc_exp" , title, nbinsx, xlow, xup);
  TH1D * hdata_obs = new TH1D("data_obs" , title, nbinsx, xlow, xup);  
  TH1D * hm350 = new TH1D("m350" , title, nbinsx, xlow, xup);

  // Load Trees, then fill histograms
  TTree* ttbar;
  TTree* qcd;
  TTree* znn;
  TTree* wt;
  TTree* other;
  TTree* data_obs;
  TTree* m350;

  if(baselineCuts) {
    ttbar_->Project("ttbar",cut->GetVar(),baseline_);
    qcd_->Project("qcd",cut->GetVar(),baseline_);
    znn_->Project("znn",cut->GetVar(),baseline_);
    wt_->Project("wt",cut->GetVar(),baseline_);
    other_->Project("other",cut->GetVar(),baseline_);
    if (plotData) {
      data_obs_->Project("data_obs",cut->GetVar(),baseline_);
    }
    if (plotSig) {
      m350_->Project("m350",cut->GetVar(),baseline_);
    }
  }
  else {
    ttbar = (TTree*) ttbar_->CopyTree(cut->Get_N_Minus_1_Cuts());
    cout << "... DONE: ttbar copy tree. N=" << ttbar->GetEntries() << endl;
    ttbar->Project("ttbar",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    qcd = (TTree*) qcd_->CopyTree(cut->Get_N_Minus_1_Cuts());
    cout << "... DONE: qcd copy tree. N=" << qcd->GetEntries() << endl;
    qcd->Project("qcd",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    znn = (TTree*) znn_->CopyTree(cut->Get_N_Minus_1_Cuts());
    cout << "... DONE: znn copy tree. N=" << znn->GetEntries() << endl;
    znn->Project("znn",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    wt = (TTree*) wt_->CopyTree(cut->Get_N_Minus_1_Cuts());
    cout << "... DONE: wt copy tree. N=" << wt->GetEntries() << endl;
    wt->Project("wt",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    other = (TTree*) other_->CopyTree(cut->Get_N_Minus_1_Cuts());
    cout << "... DONE: other copy tree. N=" << other->GetEntries() << endl;
    other->Project("other",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    if (plotData) {
      data_obs = (TTree*) data_obs_->CopyTree(cut->Get_N_Minus_1_Cuts());
      cout << "... DONE: data_obs copy tree. N=" << data_obs->GetEntries() << endl;
      data_obs->Project("data_obs",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    }
    if (plotSig) {
      m350 = (TTree*) m350_->CopyTree(cut->Get_N_Minus_1_Cuts());
      cout << "... DONE: m350 copy tree. N=" << m350->GetEntries() << endl;
      m350->Project("m350",cut->GetVar(),cut->Get_N_Minus_1_Cuts()*"full_weight");
    }
  }

  // Add up MC histograms
  hmc_exp->Add(httbar);
  hmc_exp->Add(hqcd);
  hmc_exp->Add(hznn);
  hmc_exp->Add(hwt);
  hmc_exp->Add(hother);
  cout << "... DONE: add all backgrounds to mc_exp." << endl;

  // store counts for eff calculation
  vector<double> n_minus_1;
  n_minus_1.push_back(httbar->Integral());
  n_minus_1.push_back(hqcd->Integral());
  n_minus_1.push_back(hznn->Integral());
  n_minus_1.push_back(hwt->Integral());
  n_minus_1.push_back(hother->Integral());
  n_minus_1.push_back(hmc_exp->Integral());
  n_minus_1.push_back(hm350->Integral());

  vector<double> after_cut; // need a smarter way to save these?
  after_cut.push_back(6.30543868345375813);
  after_cut.push_back(0);
  after_cut.push_back(0.0519032645970582962);
  after_cut.push_back(0);
  after_cut.push_back(0);
  after_cut.push_back(6.35734194805081643);
  after_cut.push_back(5.18291499199057171);

  // save to efficiency log file
  logEff << cut->Get_Var_Cut() << endl;
  logEff << "Baseline Efficiency: " << endl;
  for (unsigned int sample(0);sample<after_cut.size();sample++) {
    logEff << baseline_counts_a[sample]/baseline_counts[sample];
    if (sample<after_cut.size()-1) logEff << " && ";
    else logEff << " \\" << endl;
  }
  logEff << "N minus 1 Efficiency: " << endl;
  for (unsigned int sample(0);sample<after_cut.size();sample++) {
    logEff << after_cut[sample]/n_minus_1[sample];
    if (sample<after_cut.size()-1) logEff << " && ";
    else logEff << " \\" << endl;
  }

  THStack * hs = new THStack("hs", "");
  hs->Add(httbar);
  hs->Add(hqcd);
  hs->Add(hznn);
  hs->Add(hwt);
  hs->Add(hother);
  // if (plotSig) hs->Add(hm350);

  // Setup canvas and pads
  TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
  TPad * pad1 = new TPad("pad1", "top pad" , 0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy(plotLog);

  // Setup histogram styles
  set_style(httbar, "ttbar");
  set_style(hqcd, "qcd");
  set_style(hznn, "znn");
  set_style(hwt, "wt");
  set_style(hother, "other");
  set_style(hdata_obs, "data_obs");
  set_style(hm350, "m350");

  // Setup auxiliary histograms (ratios, errors, etc)
  TH1D * staterr = (TH1D *) hmc_exp->Clone("staterr");
  staterr->Sumw2();
  //staterr->SetFillColor(kRed);
  staterr->SetFillColor(kGray+3);
  staterr->SetMarkerSize(0);
  staterr->SetFillStyle(3013);

  TH1D * ratio = (TH1D *) hdata_obs->Clone("ratio");
  ratio->Sumw2();
  ratio->SetMarkerSize(0.8);
  //ratio->SetMarkerSize(0.5);
  ratio->Divide(hdata_obs, hmc_exp, 1., 1., "");

  TH1D * ratiostaterr = (TH1D *) hmc_exp->Clone("ratiostaterr");
  ratiostaterr->Sumw2();
  ratiostaterr->SetStats(0);
  ratiostaterr->SetTitle(title);
  ratiostaterr->GetYaxis()->SetTitle("Data/MC");
  ratiostaterr->SetMaximum(2.2);
  ratiostaterr->SetMinimum(0);
  ratiostaterr->SetMarkerSize(0);
  //ratiostaterr->SetFillColor(kRed);
  ratiostaterr->SetFillColor(kGray+3);
  ratiostaterr->SetFillStyle(3013);
  ratiostaterr->GetXaxis()->SetLabelSize(0.12);
  ratiostaterr->GetXaxis()->SetTitleSize(0.14);
  ratiostaterr->GetXaxis()->SetTitleOffset(1.10);
  ratiostaterr->GetYaxis()->SetLabelSize(0.10);
  ratiostaterr->GetYaxis()->SetTitleSize(0.12);
  ratiostaterr->GetYaxis()->SetTitleOffset(0.6);
  ratiostaterr->GetYaxis()->SetNdivisions(505);

  TLine* ratiounity = new TLine(xlow,1,xup,1);
  ratiounity->SetLineStyle(2);

  for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
    ratiostaterr->SetBinContent(i, 1.0);
    if (hmc_exp->GetBinContent(i) > 1e-6) { //< not empty
      double binerror = hmc_exp->GetBinError(i) / hmc_exp->GetBinContent(i);
      ratiostaterr->SetBinError(i, binerror);
    } else {
      ratiostaterr->SetBinError(i, 999.);
    }
  }

  TH1D * ratiosysterr = (TH1D *) ratiostaterr->Clone("ratiosysterr");
  ratiosysterr->Sumw2();
  ratiosysterr->SetMarkerSize(0);
  ratiosysterr->SetFillColor(kYellow-4);
  //ratiosysterr->SetFillStyle(3002);
  ratiosysterr->SetFillStyle(1001);

  for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
    if (hmc_exp->GetBinContent(i) > 1e-6) { //< not empty
      double binerror2 = (pow(hmc_exp->GetBinError(i), 2) +
			  pow(httbar->GetBinContent(i), 2) +
			  pow(hqcd->GetBinContent(i), 2) +
			  pow(hznn->GetBinContent(i), 2) +
			  pow(hwt->GetBinContent(i), 2) +
			  pow(hother->GetBinContent(i), 2)
			  );

      double binerror = sqrt(binerror2);
      ratiosysterr->SetBinError(i, binerror / hmc_exp->GetBinContent(i));
    }
  }
  
  // Setup legends
  TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
  set_style(leg1);
  if (plotData) leg1->AddEntry(hdata_obs, "Data", "pel");
  if (plotSig) leg1->AddEntry(hm350, "m_{#chi} = 350 GeV", "l");
  leg1->AddEntry(httbar, "t#bar{t}", "f");
  leg1->AddEntry(hqcd, "qcd", "f");
  leg1->AddEntry(hznn, "znn", "f");

  TLegend * leg2 = new TLegend(0.72, 0.68, 0.94, 0.92);
  set_style(leg2);
  leg2->AddEntry(hwt, "wt", "f");
  leg2->AddEntry(hother, "other", "f");
  leg2->AddEntry(staterr, "MC uncert. (stat)", "f");

  TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
  set_style(ratioleg);
  ratioleg->SetTextSize(0.07);
  ratioleg->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");

  // Draw MC signal and backgrounds
  cout << "MakePlots(): Drawing..." << endl;
  pad1->cd();
  if (plotLog) pad1->SetLogy();

  double ymax = TMath::Max(hdata_obs->GetMaximum(), hs->GetMaximum());
  hs->SetMaximum(ymax * 1.7 + (ymax>1 ? sqrt(ymax) : 0.));
  if (plotLog) hs->SetMaximum(ymax * 200 + (ymax>1 ? sqrt(ymax) : 0.));
  hs->SetMinimum(0.01);
  hs->Draw("hist");
  hs->GetXaxis()->SetLabelSize(0);
  double binwidth = (xup - xlow) / nbinsx;
  TString ytitle = Form("Events / %.3f", binwidth);
  hs->GetYaxis()->SetTitle(ytitle);

  staterr->Draw("e2 same");
  if (plotSig) {
    hm350->SetLineColor(2);
    hm350->SetLineWidth(3);
    hm350->SetFillColor(0);
    hm350->Draw("hist same");
  }
  if (plotData) {
    hdata_obs->Draw("e1 same");
  }

  // Draw legends
  leg1->Draw();
  leg2->Draw();
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextFont(62);
  latex->SetTextSize(0.052);
  latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 8 TeV, L = 19.4 fb^{-1}");
  // something more?
  //latex->DrawLatex(0.19, 0.79, "More schtuff");

  // Draw ratio
  pad2->cd();
  pad2->SetGridy(0);
  ratiostaterr->Draw("e2");
  //ratiosysterr->Draw("e2 same");
  ratiostaterr->Draw("e2 same");
  ratiounity->Draw();
  ratio->Draw("e1 same");
  ratioleg->Draw();

  // Kolmogorov-Smirnov test and Chi2 test
  TPaveText * pave = new TPaveText(0.18, 0.86, 0.28, 0.96, "brNDC");
  if (plotData) {
    pave->SetLineColor(0);
    pave->SetFillColor(0);
    pave->SetShadowColor(0);
    pave->SetBorderSize(1);
    double nchisq = hdata_obs->Chi2Test(hmc_exp, "UWCHI2/NDF"); // MC uncert. (stat)
    //double kolprob = hdata_obs->KolmogorovTest(hmc_exp); // MC uncert. (stat)
    TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f", nchisq));
    //TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f", nchisq, kolprob));
    text->SetTextFont(62);
    text->SetTextSize(0.07);
    //text->SetTextSize(0.06);
    pave->Draw();
  }

  // Print
  cout << "MakePlots(): Printing..." << endl;
  pad1->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();
  pad2->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();
  c1->cd();

  gPad->Print(plotdir+plotname+".png");
  gPad->Print(plotdir+plotname+".pdf");

  TFile* outrootfile = TFile::Open(plotdir+plotname+".root", "RECREATE");
  httbar->Write();
  hqcd->Write();
  hznn->Write();
  hwt->Write();
  hother->Write();
  hm350->Write();
  hmc_exp->Write();
  hdata_obs->Write();
  outrootfile->Close();

  // Clean up
  delete staterr;
  delete ratio;
  delete ratiostaterr;
  delete ratiosysterr;
  delete leg1;
  delete leg2;
  delete ratioleg;
  delete latex;
  delete pave;
  delete hs;
  delete pad1;
  delete pad2;
  delete c1;

  delete httbar;
  delete hqcd;
  delete hznn;
  delete hwt;
  delete hother;
  delete hm350;
  delete httbar_b;
  delete hqcd_b;
  delete hznn_b;
  delete hwt_b;
  delete hother_b;
  delete hm350_b;
  delete hdata_obs;
  delete hmc_exp;

  cout << "MakePlots(): DONE!" << endl;

  return;
}

int main() {

  logEff.open("efficiencies.log", std::ofstream::out | std::ofstream::app);

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);
  TString plotdir = "plots/";

  if (gSystem->AccessPathName(plotdir))
    gSystem->mkdir(plotdir);

  

  // This is for calculating efficiency relative to the baseline
  TCut baseline("passesJSONCut&&passesPVCut&&passesJet2PtCut&&passes2CSVTCut&&passesMETCleaningCut&&passesTriggerCut");

  ttbar_ch.Add("../../reduced_trees/TTJets_FullLeptMGDecays*.root");
  ttbar_ch.Add("../../reduced_trees/TTJets_SemiLeptMGDecays*v71.root");
  qcd_ch.Add("../../reduced_trees/BJets*.root");
  qcd_ch.Add("../../reduced_trees/TTJets_Had*.root");
  znn_ch.Add("../../reduced_trees/ZJets*.root");
  wt_ch.Add("../../reduced_trees/W*Jets*187*.root");
  wt_ch.Add("../../reduced_trees/T_*.root");
  wt_ch.Add("../../reduced_trees/Tbar_*.root");
  other_ch.Add("../../reduced_trees/ZZ*.root");
  // other_ch.Add("../../reduced_trees/ZH*.root");
  other_ch.Add("../../reduced_trees/WZ*.root");
  other_ch.Add("../../reduced_trees/WH*.root");
  other_ch.Add("../../reduced_trees/WW*.root");
  data_obs_ch.Add("../../reduced_trees/MET_Run2012*.root");
  m350_ch.Add("../../reduced_trees/SMS-TChiHH_2b2b_2J_mChargino-350_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71_SyncSkim.root");


  // First do a quick baseline "skim" (not that quick) and load trees
  cout << "Applying baseline cuts: \n" << (string)baseline << endl;
  ttbar_ = (TTree*) ttbar_ch.CopyTree(baseline);
  cout << "... DONE: ttbar copy tree. N=" << ttbar_->GetEntries() << endl;
  qcd_ = (TTree*) qcd_ch.CopyTree(baseline);
  cout << "... DONE: qcd copy tree. N=" << qcd_->GetEntries() << endl;
  znn_ = (TTree*) znn_ch.CopyTree(baseline);
  cout << "... DONE: znn copy tree. N=" << znn_->GetEntries() << endl;
  wt_ = (TTree*) wt_ch.CopyTree(baseline);
  cout << "... DONE: wt copy tree. N=" << wt_->GetEntries() << endl;
  other_ = (TTree*) other_ch.CopyTree(baseline);
  cout << "... DONE: other copy tree. N=" << other_->GetEntries() << endl;
  
  data_obs_ = (TTree*) data_obs_ch.CopyTree(baseline);
  cout << "... DONE: data_obs copy tree. N=" << data_obs_->GetEntries() << endl;  
  m350_ = (TTree*) m350_ch.CopyTree(baseline);
  cout << "... DONE: m350 copy tree. N=" << m350_->GetEntries() << endl;

  
  // The N Minus 1 plots
  Cut* max_delta_R_cut = new Cut("max_delta_R");
  cout << "Plotting " << max_delta_R_cut->GetVar() << endl;
  cout << "The cut on this var is " << (string)max_delta_R_cut->Get_Var_Cut() << endl;
  cout << "Cuts: " << (string)max_delta_R_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(max_delta_R_cut, 
		 ";#DeltaR_{max};Events / 0.2", 30, 0., 6.,
		 "max_delta_R_n_minus_1", plotdir,
		 "plotSig");
  Cut* met_sig_cut = new Cut("met_sig");
  cout << "Plotting " << met_sig_cut->GetVar() << endl;
  cout << "Cuts: " << (string)met_sig_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(met_sig_cut, 
		 ";S_{MET};Events / 20", 10, 0., 200.,
		 "met_sig_n_minus_1", plotdir,
		 "plotSig");
  Cut* average_higgs_mass_cut = new Cut("average_higgs_mass");
  cout << "Plotting " << average_higgs_mass_cut->GetVar() << endl;
  cout << "Cuts: " << (string)average_higgs_mass_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(average_higgs_mass_cut, 
		 ";<m_{bb}> [GeV];Events / 10 GeV", 20, 0., 200.,
		 "average_higgs_mass_n_minus_1", plotdir,
		 "plotSig");
  Cut* min_delta_phi_cut = new Cut("min_delta_phi");
  cout << "Plotting " << min_delta_phi_cut->GetVar() << endl;
  cout << "Cuts: " << (string)min_delta_phi_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(min_delta_phi_cut, 
		 ";#Delta#phi_{min}(jet,MET);Events / 0.2", 15, 0., 3.,
		 "min_delta_phi_n_minus_1", plotdir,
		 "plotSig");
  Cut* num_iso_tracks_cut = new Cut("num_iso_tracks");
  cout << "Plotting " << num_iso_tracks_cut->GetVar() << endl;
  cout << "Cuts: " << (string)num_iso_tracks_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(num_iso_tracks_cut, 
		 ";n_{IsoTracks};Events", 5, 0., 5.,
		 "num_iso_tracks_n_minus_1", plotdir,
		 "plotSig");
  Cut* higgs_mass_difference_cut = new Cut("higgs_mass_difference");
  cout << "Plotting " << higgs_mass_difference_cut->GetVar() << endl;
  cout << "Cuts: " << (string)higgs_mass_difference_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(higgs_mass_difference_cut,
		 ";|#Deltam_{bb}| [GeV];Events / 10 GeV", 10, 0., 100.,
		 "higgs_mass_difference_n_minus_1", plotdir,
		 "plotSig");
  Cut* num_jets_cut = new Cut("num_jets");
  cout << "Plotting " << num_jets_cut->GetVar() << endl;
  cout << "Cuts: " << (string)num_jets_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(num_jets_cut, 
		 ";n_{jets} (p_{T} > 20 GeV);Events", 8, 0., 8.,
		 "num_jets_n_minus_1", plotdir,
		 "plotSig");
  Cut* third_highest_csv_cut = new Cut("third_highest_csv");
  cout << "Plotting " << third_highest_csv_cut->GetVar() << endl;
  cout << "Cuts: " << (string)third_highest_csv_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(third_highest_csv_cut, 
		 ";CSV_{3};Events / 0.1", 10, 0., 1.,
		 "third_highest_csv_n_minus_1", plotdir,
		 "plotSig");
  Cut* fourth_highest_csv_cut = new Cut("fourth_highest_csv");
  cout << "Plotting " << fourth_highest_csv_cut->GetVar() << endl;
  cout << "Cuts: " << (string)fourth_highest_csv_cut->Get_N_Minus_1_Cuts() << endl;
  n_minus_1_plot(fourth_highest_csv_cut,
		 ";CSV_{4};Events / 0.1", 10, 0., 1.,
		 "fourth_highest_csv_n_minus_1", plotdir,
		 "plotSig");

  // The Plots with just the baseline selection
  /*
  n_minus_1_plot(max_delta_R_cut, 
		 ";#DeltaR_{max};Events / 0.2", 30, 0., 6.,
		 "max_delta_R_baseline", plotdir,
		 "plotSig:baselineCuts");
  n_minus_1_plot(met_sig_cut, 
		 ";S_{MET};Events / 20", 10, 0., 200.,
		 "met_sig_baseline", plotdir,
		 "plotSig:baselineCuts:plotLog");
  n_minus_1_plot(average_higgs_mass_cut, 
		 ";<m_{bb}> [GeV];Events / 10 GeV", 20, 0., 200.,
		 "average_higgs_mass_baseline", plotdir,
		 "plotSig:baselineCuts");
  n_minus_1_plot(min_delta_phi_cut, 
		 ";#Delta#phi_{min}(jet,MET);Events / 0.2", 15, 0., 3.,
		 "min_delta_phi_baseline", plotdir,
		 "plotSig:baselineCuts");
  n_minus_1_plot(num_iso_tracks_cut, 
		 ";n_{IsoTracks};Events", 5, 0., 5.,
		 "num_iso_tracks_baseline", plotdir,
		 "plotSig:baselineCuts");
   n_minus_1_plot(higgs_mass_difference_cut,
		 ";|#Deltam_{bb}| [GeV];Events / 10 GeV", 10, 0., 100.,
		 "higgs_mass_difference_baseline", plotdir,
		 "plotSig:baselineCuts");
  n_minus_1_plot(num_jets_cut, 
		 ";n_{jets} (p_{T} > 20 GeV);Events", 8, 0., 8.,
		 "num_jets_baseline", plotdir,
		 "plotSig:baselineCuts");
  n_minus_1_plot(third_highest_csv_cut, 
		 ";CSV_{3};Events / 0.1", 10, 0., 1.,
		 "third_highest_csv_baseline", plotdir,
		 "plotSig:baselineCuts:plotLog");
  n_minus_1_plot(fourth_highest_csv_cut,
		 ";CSV_{4};Events / 0.1", 10, 0., 1.,
		 "fourth_highest_csv_baselinelt plo", plotdir,
		 "plotSig:baselineCuts:plotLog");*/


  logEff.close();

  return 0;
}
