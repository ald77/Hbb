#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "style.hpp"

int main(){
  SetStyle();
  TH1::SetDefaultSumw2();

  TChain data("data","data");
  TChain ttbar("ttbar","ttbar");
  TChain qcd("qcd","qcd");
  TChain single_t_or_boson("single_t_or_boson","single_t_or_boson");
  TChain diboson("diboson","diboson");

  data.Add("reduced_trees_test/MET_*2012*.root/reduced_tree");
  ttbar.Add("reduced_trees_test/TTJets_FullLept*.root/reduced_tree");
  ttbar.Add("reduced_trees_test/TTJets_SemiLept*.root/reduced_tree");
  qcd.Add("reduced_trees_test/TTJets_Hadronic*.root/reduced_tree");
  qcd.Add("reduced_trees_test/BJets*.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees_test/*channel*.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees_test/*JetsTo*.root/reduced_tree");
  diboson.Add("reduced_trees_test/WH*.root/reduced_tree");
  diboson.Add("reduced_trees_test/WW*.root/reduced_tree");
  diboson.Add("reduced_trees_test/WZ*.root/reduced_tree");
  diboson.Add("reduced_trees_test/ZZ*.root/reduced_tree");

  std::vector<TChain*> chains(0);
  chains.push_back(&data);
  chains.push_back(&qcd);
  chains.push_back(&diboson);
  chains.push_back(&single_t_or_boson);
  chains.push_back(&ttbar);

  std::vector<TH1D> h_before(0), h_after(0);
  std::vector<TH1D> h_top_pt_before(0), h_top_pt_after(0);

  for(unsigned histo(0); histo<chains.size(); ++histo){
    const std::string name(chains.at(histo)->GetName());
    h_before.push_back(TH1D(("h_"+name+"_before").c_str(), "Good Vertices Before Pileup Reweighting (baseline selection);Vertices;Events/1", 61, -0.5, 60.5));
    h_after.push_back(TH1D(("h_"+name+"_after").c_str(), "Good Vertices After Pileup Reweighting (baseline selection);Vertices;Events/1", 61, -0.5, 60.5));
    h_before.at(histo).SetLineColor(histo+1);
    h_before.at(histo).SetFillColor(histo+1);
    h_after.at(histo).SetLineColor(histo+1);
    h_after.at(histo).SetFillColor(histo+1);
  }

  THStack s_before("s_before","Good Vertices Before Pileup Reweighting (baseline selection);Vertices;Events/1");
  THStack s_after("s_after","Good Vertices After Pileup Reweighting (baseline selection);Vertices;Events/1");
  TH1D tot_before(h_before.at(0));
  TH1D tot_after(h_after.at(0));
  bool passesBaselineSelection(false), passesTriggerPlateauCuts(false);

  float pu_true_num_interactions(0.0);
  unsigned short num_primary_vertices(0);

  float full_weight(0.0), lumi_weight(0.0), top_pt_weight(0.0), pu_weight(0.0), trigger_weight(0.0);

  for(unsigned chain(0); chain<chains.size(); ++chain){
    chains.at(chain)->SetBranchStatus("*",0);
    chains.at(chain)->SetBranchStatus("passesTriggerPlateauCuts",1);
    chains.at(chain)->SetBranchStatus("passesBaselineSelection",1);
    chains.at(chain)->SetBranchStatus("pu_true_num_interactions",1);
    chains.at(chain)->SetBranchStatus("num_primary_vertices",1);
    chains.at(chain)->SetBranchStatus("full_weight",1);
    chains.at(chain)->SetBranchStatus("lumi_weight",1);
    chains.at(chain)->SetBranchStatus("pu_weight",1);
    chains.at(chain)->SetBranchStatus("top_pt_weight",1);
    chains.at(chain)->SetBranchStatus("trigger_weight",1);
    chains.at(chain)->SetBranchAddress("passesTriggerPlateauCuts",&passesTriggerPlateauCuts);
    chains.at(chain)->SetBranchAddress("passesBaselineSelection",&passesBaselineSelection);
    chains.at(chain)->SetBranchAddress("pu_true_num_interactions",&pu_true_num_interactions);
    chains.at(chain)->SetBranchAddress("num_primary_vertices",&num_primary_vertices);
    chains.at(chain)->SetBranchAddress("full_weight", &full_weight);
    chains.at(chain)->SetBranchAddress("lumi_weight", &lumi_weight);
    chains.at(chain)->SetBranchAddress("pu_weight", &pu_weight);
    chains.at(chain)->SetBranchAddress("top_pt_weight", &top_pt_weight);
    chains.at(chain)->SetBranchAddress("trigger_weight", &trigger_weight);
    int num_entries(chains.at(chain)->GetEntries());

    for(int entry(0); entry<num_entries; ++entry){
      chains.at(chain)->GetEntry(entry);

      if(passesBaselineSelection){
	const double no_pu_weight(lumi_weight*top_pt_weight*trigger_weight);
	h_before.at(chain).Fill(num_primary_vertices,no_pu_weight);
	h_after.at(chain).Fill(num_primary_vertices,full_weight);
      }
    }

    if(chain!=0){
      s_before.Add(&h_before.at(chain));
      s_after.Add(&h_after.at(chain));
      tot_before.Add(&h_before.at(chain));
      tot_after.Add(&h_after.at(chain));
    }
  }

  const double before_max(std::max(s_before.GetMaximum(), h_before.at(0).GetMaximum()));
  s_before.SetMaximum(before_max);
  h_before.at(0).SetMaximum(before_max);

  const double after_max(std::max(s_after.GetMaximum(), h_after.at(0).GetMaximum()));
  s_after.SetMaximum(after_max);
  h_after.at(0).SetMaximum(after_max);

  TH1D rat_before("rat_before",";Vertices;MC/Data",61,-0.5,60.5);
  std::cout << tot_before.GetNbinsX() << " " << h_before.at(0).GetNbinsX() << std::endl;
  rat_before.Divide(&tot_before, &h_before.at(0));
  std::cout << "NO" << std::endl;
  //rat_before.SetMinimum(0.0);
  //rat_before.SetMaximum(2.0);
  rat_before.SetStats(0);
  TH1D rat_after("rat_before",";Vertices;MC/Data",61,-0.5,60.5);
  rat_after.Divide(&tot_after, &h_after.at(0));
  //rat_after.SetMinimum(0.0);
  //rat_after.SetMaximum(2.0);
  rat_after.SetStats(0);

  const double chisq_before(h_before.at(0).Chi2Test(&tot_before, "UWOFUFP"));
  const double chisq_after(h_after.at(0).Chi2Test(&tot_after, "UWOFUFP"));
  TPaveText text_before(0.7,0.6,0.9,0.7);
  TPaveText text_after(0.7,0.6,0.9,0.7);
  std::ostringstream oss(std::string(""));
  oss << "p-value: " << chisq_before << std::endl;
  text_before.AddText(oss.str().c_str());
  oss.str()="";
  oss << "p-value: " << chisq_after << std::endl;
  text_after.AddText(oss.str().c_str());

  const double div(1.0/3.0);

  TLegend legend(0.7, 0.5, 0.95, 1.0-0.15/(1.0-div));

  for(unsigned histo(0); histo<h_before.size(); ++histo){
    if(histo==0){
      legend.AddEntry(&h_before.at(histo), chains.at(histo)->GetName(), "lpe");
    }else{
      legend.AddEntry(&h_before.at(histo), chains.at(histo)->GetName(), "f");
    }
    h_before.at(histo).SetTitleSize(0,"X");
    h_before.at(histo).SetTitleSize(0.08,"Y");
    h_before.at(histo).SetTitleOffset(0.0,"Y");
    h_after.at(histo).SetTitleSize(0,"X");
  }
  tot_before.SetTitleSize(0,"X");
  tot_before.SetTitleSize(0.08,"Y");
  tot_before.SetTitleOffset(0.0,"Y");
  tot_after.SetTitleSize(0,"X");
  rat_before.SetTitleOffset(0.0,"Y");
  rat_before.SetTitleSize(0.08,"Y");

  TCanvas canvas_before;
  canvas_before.Divide(2);
  canvas_before.cd(1);
  canvas_before.GetPad(1)->SetPad(0.0,div,1.0,1.0);
  canvas_before.GetPad(1)->SetMargin(0.15,0.05,0.0,0.15/(1.0-div));
  canvas_before.GetPad(1)->SetLogy(1);
  s_before.Draw("hist");
  h_before.at(0).Draw("samee1p");
  legend.Draw("same");
  //text_before.Draw("same");
  canvas_before.cd(2);
  canvas_before.GetPad(2)->SetPad(0.0,0.0,1.0,div);
  canvas_before.GetPad(2)->SetMargin(0.15,0.05,0.15/div,0.0);
  canvas_before.GetPad(2)->SetGridy(1);
  canvas_before.GetPad(2)->SetLogy(1);
  rat_before.Draw("e1p");
  canvas_before.cd(0);
  canvas_before.Print("before.pdf");

  TCanvas canvas_after;
  canvas_after.Divide(2);
  canvas_after.cd(1);
  canvas_after.GetPad(1)->SetPad(0.0,div,1.0,1.0);
  canvas_after.GetPad(1)->SetMargin(0.15,0.05,0.0,0.15/(1.0-div));
  canvas_after.GetPad(1)->SetLogy(1);
  s_after.Draw("hist");
  h_after.at(0).Draw("samee1p");
  legend.Draw("same");
  //text_after.Draw("same");
  canvas_after.cd(2);
  canvas_after.GetPad(2)->SetPad(0.0,0.0,1.0,div);
  canvas_after.GetPad(2)->SetMargin(0.15,0.05,0.15/div,0.0);
  canvas_after.GetPad(2)->SetGridy(1);
  canvas_after.GetPad(2)->SetLogy(1);
  rat_after.Draw("e1p");
  canvas_after.cd(0);
  canvas_after.Print("after.pdf");
}
