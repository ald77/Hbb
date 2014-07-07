#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "timer.hpp"
#include "style.hpp"

int main(){
  SetStyle();
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);
  TFile fast_file("reduced_trees/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V7C_FSIM-v2_AODSIM_UCSB1976_v71.root","read");
  TFile full_file("reduced_trees/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1850_v71.root","read");

  TChain* fast_chain(static_cast<TChain*>(fast_file.Get("reduced_tree")));
  TChain* full_chain(static_cast<TChain*>(full_file.Get("reduced_tree")));

  //const std::string selection("full_weight");
  //const std::string selection("full_weight*(passesTriggerCut && passesTriggerPlateauCuts)");
  const std::string selection("full_weight*(passesBaselineSelection==1)");

  TCanvas canvas("canvas","canvas");
  const int num_leaves(fast_chain->GetListOfLeaves()->GetSize());
  Timer timer(num_leaves);
  timer.Start();
  for(int leaf(0); leaf<num_leaves; ++leaf){
    timer.PrintRemainingTime();
    const std::string leaf_name(static_cast<TLeaf*>((fast_chain->GetListOfLeaves()->At(leaf)))->GetBranch()->GetName());
    if(leaf_name.find("trigger")==std::string::npos){
      fast_chain->Draw((leaf_name+">>fast_"+leaf_name).c_str(), selection.c_str());
      full_chain->Draw((leaf_name+">>full_"+leaf_name).c_str(), selection.c_str());
      TH1F* h_fast(static_cast<TH1F*>(gDirectory->Get(("fast_"+leaf_name).c_str())));
      TH1F* h_full(static_cast<TH1F*>(gDirectory->Get(("full_"+leaf_name).c_str())));
      h_fast->Scale(1.0/h_fast->Integral("width"));
      h_full->Scale(1.0/h_full->Integral("width"));
      h_fast->SetLineColor(1);
      h_full->SetLineColor(2);
      h_full->Draw("hist");
      h_fast->Draw("histsame");
      TLegend legend(0.7,0.7,0.95,0.85);
      legend.AddEntry(h_fast, "Fast Sim", "lpe");
      legend.AddEntry(h_full, "Full Sim", "lpe");
      legend.Draw("same");
      canvas.SetLogx(0);
      canvas.SetLogy(0);
      canvas.Print(("pdfs/baseline/"+leaf_name+"_lin.pdf").c_str());
      canvas.SetLogy(1);
      canvas.Print(("pdfs/baseline/"+leaf_name+"_log.pdf").c_str());
      canvas.SetLogx(1);
      canvas.Print(("pdfs/baseline/"+leaf_name+"_llg.pdf").c_str());
    }
    timer.Iterate();
  }
}
