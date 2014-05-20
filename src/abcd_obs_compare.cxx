#include <cmath>
#include <iostream>
#include <string>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "style.hpp"

void set(TH1D& h, unsigned short bin, double count);
void set(TH1D& h, unsigned short bin, double count, double err);
void set(TChain& chain, TH1D& h, const std::string& cut);
double get(TChain& chain, const std::string& cut);
void draw(TH1D count, TH1D abcd, TH1D sig, const std::string&);

int main(){
  SetStyle();
  TH1D on_2b_count("on_2b_count","S_{MET} Bin Counts (2b full sel, #Delta R<2.2);S_{MET} bin;Events",
		   4, 0.5, 4.5);
  TH1D on_2b_abcd("on_2b_abcd","S_{MET} Bin Counts (2b full sel, #Delta R<2.2);S_{MET} bin;Events",
		  4, 0.5, 4.5);
  TH1D on_2b_sig("on_2b_sig","S_{MET} Bin Counts (2b full sel, #Delta R<2.2);S_{MET} bin;Events",
		 4, 0.5, 4.5);
  TH1D on_3b_count("on_3b_count","S_{MET} Bin Counts (3b full sel, #Delta R<2.2);S_{MET} bin;Events",
		   4, 0.5, 4.5);
  TH1D on_3b_abcd("on_3b_abcd","S_{MET} Bin Counts (3b full sel, #Delta R<2.2);S_{MET} bin;Events",
		  4, 0.5, 4.5);
  TH1D on_3b_sig("on_3b_sig","S_{MET} Bin Counts (3b full sel, #Delta R<2.2);S_{MET} bin;Events",
		 4, 0.5, 4.5);
  TH1D on_4b_count("on_4b_count","S_{MET} Bin Counts (4b full sel, #Delta R<2.2);S_{MET} bin;Events",
		   4, 0.5, 4.5);
  TH1D on_4b_abcd("on_4b_abcd","S_{MET} Bin Counts (4b full sel, #Delta R<2.2);S_{MET} bin;Events",
		  4, 0.5, 4.5);
  TH1D on_4b_sig("on_4b_sig","S_{MET} Bin Counts (4b full sel, #Delta R<2.2);S_{MET} bin;Events",
		 4, 0.5, 4.5);

  set(on_2b_count, 1, 60.0);
  set(on_2b_count, 2, 105.0);
  set(on_2b_count, 3, 25.0);
  set(on_2b_count, 4, 15.0);
  set(on_3b_count, 1, 4.0);
  set(on_3b_count, 2, 15.0);
  set(on_3b_count, 3, 1.0);
  set(on_3b_count, 4, 0.0);
  set(on_4b_count, 1, 4.0);
  set(on_4b_count, 2, 7.0);
  set(on_4b_count, 3, 3.0);
  set(on_4b_count, 4, 0.0);
  set(on_2b_abcd, 1, 60.0);
  set(on_2b_abcd, 2, 105.0);
  set(on_2b_abcd, 3, 25.0);
  set(on_2b_abcd, 4, 15.0);
  set(on_3b_abcd, 1, 6.64, 1.27);
  set(on_3b_abcd, 2, 11.4, 1.71);
  set(on_3b_abcd, 3, 2.24, 0.710);
  set(on_3b_abcd, 4, 1.54, 0.731);
  set(on_4b_abcd, 1, 3.02, 0.73);
  set(on_4b_abcd, 2, 4.78, 0.94);
  set(on_4b_abcd, 3, 0.622, 0.308);
  set(on_4b_abcd, 4, 0.441, 0.336);

  TH1D off_2b_count("off_2b_count","S_{MET} Bin Counts (2b full sel, no #Delta R cut);S_{MET} bin;Events",
		    4, 0.5, 4.5);
  TH1D off_2b_abcd("off_2b_abcd","S_{MET} Bin Counts (2b full sel, no #Delta R cut);S_{MET} bin;Events",
		   4, 0.5, 4.5);
  TH1D off_2b_sig("off_2b_sig","S_{MET} Bin Counts (2b full sel, no #Delta R cut);S_{MET} bin;Events",
		  4, 0.5, 4.5);
  TH1D off_3b_count("off_3b_count","S_{MET} Bin Counts (3b full sel, no #Delta R cut);S_{MET} bin;Events",
		    4, 0.5, 4.5);
  TH1D off_3b_abcd("off_3b_abcd","S_{MET} Bin Counts (3b full sel, no #Delta R cut);S_{MET} bin;Events",
		   4, 0.5, 4.5);
  TH1D off_3b_sig("off_3b_sig","S_{MET} Bin Counts (3b full sel, no #Delta R cut);S_{MET} bin;Events",
		  4, 0.5, 4.5);
  TH1D off_4b_count("off_4b_count","S_{MET} Bin Counts (4b full sel, no #Delta R cut);S_{MET} bin;Events",
		    4, 0.5, 4.5);
  TH1D off_4b_abcd("off_4b_abcd","S_{MET} Bin Counts (4b full sel, no #Delta R cut);S_{MET} bin;Events",
		   4, 0.5, 4.5);
  TH1D off_4b_sig("off_4b_sig","S_{MET} Bin Counts (4b full sel, no #Delta R cut);S_{MET} bin;Events",
		  4, 0.5, 4.5);

  set(off_2b_count, 1, 325.0);
  set(off_2b_count, 2, 504.0);
  set(off_2b_count, 3, 111.0);
  set(off_2b_count, 4, 30.0);
  set(off_3b_count, 1, 44.0);
  set(off_3b_count, 2, 72.0);
  set(off_3b_count, 3, 13.0);
  set(off_3b_count, 4, 1.0);
  set(off_4b_count, 1, 22.0);
  set(off_4b_count, 2, 25.0);
  set(off_4b_count, 3, 7.0);
  set(off_4b_count, 4, 0.0);
  set(off_2b_abcd, 1, 325.0);
  set(off_2b_abcd, 2, 504.0);
  set(off_2b_abcd, 3, 111.0);
  set(off_2b_abcd, 4, 30.0);
  set(off_3b_abcd, 1, 35.3, 3.6);
  set(off_3b_abcd, 2, 56.6, 4.6);
  set(off_3b_abcd, 3, 9.94, 1.8);
  set(off_3b_abcd, 4, 2.92, 0.95);
  set(off_4b_abcd, 1, 15.3, 2.1);
  set(off_4b_abcd, 2, 17.3, 2.2);
  set(off_4b_abcd, 3, 2.78, 0.85);
  set(off_4b_abcd, 4, 1.95, 0.73);

  TChain chain("","");
  chain.Add("reduced_trees/SMS-TChiHH_2b2b_2J_mChargino-250_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71_SyncSkim.root/reduced_tree");
  set(chain, on_2b_sig, "num_b_tagged_jets==2 && passesDRCut");
  set(chain, on_3b_sig, "num_b_tagged_jets==3 && passesDRCut");
  set(chain, on_4b_sig, "num_b_tagged_jets>=4 && passesDRCut");
  set(chain, off_2b_sig, "num_b_tagged_jets==2 && !passesDRCut");
  set(chain, off_3b_sig, "num_b_tagged_jets==3 && !passesDRCut");
  set(chain, off_4b_sig, "num_b_tagged_jets>=4 && !passesDRCut");

  on_2b_abcd.SetFillColor(4);
  on_3b_abcd.SetFillColor(4);
  on_4b_abcd.SetFillColor(4);
  off_2b_abcd.SetFillColor(4);
  off_3b_abcd.SetFillColor(4);
  off_4b_abcd.SetFillColor(4);
  on_2b_sig.SetFillColor(2);
  on_3b_sig.SetFillColor(2);
  on_4b_sig.SetFillColor(2);
  off_2b_sig.SetFillColor(2);
  off_3b_sig.SetFillColor(2);
  off_4b_sig.SetFillColor(2);

  draw(on_2b_count, on_2b_abcd, on_2b_sig, "crap/on_2b.pdf");
  draw(on_3b_count, on_3b_abcd, on_3b_sig, "crap/on_3b.pdf");
  draw(on_4b_count, on_4b_abcd, on_4b_sig, "crap/on_4b.pdf");
  draw(off_2b_count, off_2b_abcd, off_2b_sig, "crap/off_2b.pdf");
  draw(off_3b_count, off_3b_abcd, off_3b_sig, "crap/off_3b.pdf");
  draw(off_4b_count, off_4b_abcd, off_4b_sig, "crap/off_4b.pdf");
}

void set(TH1D& h, unsigned short bin, double count){
  h.SetBinContent(bin, count);
  h.SetBinError(bin, sqrt(count));
}

void set(TH1D& h, unsigned short bin, double count, double err){
  h.SetBinContent(bin, count);
  h.SetBinError(bin, err);
}

void set(TChain& chain, TH1D& h, const std::string& cut){
  const std::string base("0.004187617/1.33586290758103132e-05*full_weight*(passesJSONCut && passesPVCut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passes2CSVTCut && passesJet2PtCut && passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut && passesMETSig30Cut && passesHiggsAvgMassCut && passesHiggsMassDiffCut && ");
  h.SetBinContent(1, get(chain, base+cut+" && met_sig>30.0 && met_sig<50.0)"));
  h.SetBinContent(2, get(chain, base+cut+" && met_sig>50.0 && met_sig<100.0)"));
  h.SetBinContent(3, get(chain, base+cut+" && met_sig>100.0 && met_sig<150.0)"));
  h.SetBinContent(4, get(chain, base+cut+" && met_sig>150.0)"));
}

double get(TChain& chain, const std::string& cut){
  TH1D hqqq("hqqq", "hqqq", 1, -9999., 9999.);
  chain.Project("hqqq", "passesPVCut", cut.c_str());
  return hqqq.Integral();
}

void draw(TH1D count, TH1D abcd, TH1D sig, const std::string& outname){
  for(int bin(0); bin<5; ++bin){
    sig.SetBinContent(bin,sig.GetBinContent(bin)+abcd.GetBinContent(bin));
  }
  count.SetStats(0);
  abcd.SetStats(0);
  sig.SetStats(0);
  abcd.SetFillColor(kCyan);
  sig.SetFillColor(kRed);
  count.SetMinimum(0.0);
  sig.SetMinimum(0.0);
  abcd.SetMinimum(0.0);
  double max(0.0);
  double this_max(count.GetBinContent(count.GetMaximumBin()));
  if(this_max>max) max=this_max;
  this_max=abcd.GetBinContent(abcd.GetMaximumBin());
  if(this_max>max) max=this_max;
  this_max=sig.GetBinContent(sig.GetMaximumBin());
  if(this_max>max) max=this_max;
  max*=1.1;

  count.SetMaximum(max);
  abcd.SetMaximum(max);
  sig.SetMaximum(max);

  TCanvas canvas;
  sig.Draw("hist");
  abcd.Draw("histsame");
  count.Draw("samee1p");

  TLegend legend(0.7,0.5,0.95,0.85);
  legend.AddEntry(&count, "CMS data", "lpe");
  legend.AddEntry(&abcd, "ABCD pred.", "f");
  legend.AddEntry(&sig, "Signal (m=250 GeV)", "f");
  legend.Draw("same");

  canvas.Print(outname.c_str());
}
