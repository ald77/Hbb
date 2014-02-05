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

void MakeRatioPlot(std::vector<TH1D>& histos, std::vector<std::string>& names, const std::string& out_name){
  if(histos.size()==0 || names.size()!=histos.size()) return;
  const double div(1.0/3.0);
  std::ostringstream oss("");
  oss << histos.at(0).GetTitle() << ";" << histos.at(0).GetXaxis()->GetTitle() << ";" << histos.at(0).GetYaxis()->GetTitle() << std::endl;
  THStack s((std::string("s_")+histos.at(0).GetName()).c_str(),
	    oss.str().c_str());
  std::vector<double> edges(histos.at(0).GetNbinsX()+1);
  for(unsigned bin(0); bin<edges.size(); ++bin){
    edges.at(bin)=histos.at(0).GetBinLowEdge(bin+1);
  }
  TH1D* tot(static_cast<TH1D*>(histos.at(0).Clone("abcdefgh")));

  for(int bin(0); bin<=tot->GetNbinsX()+1; ++bin){
    tot->SetBinContent(bin,0.0);
  }
  for(unsigned histo(0); histo<histos.size(); ++histo){
    histos.at(histo).SetTitleSize(0,"X");
    histos.at(histo).SetTitleSize(0.08,"Y");
    //histos.at(histo).SetTitleOffset(0.0,"Y");
  }
  for(unsigned histo(1); histo<histos.size(); ++histo){
    histos.at(histo).SetFillColor(histo+1);
    s.Add(&histos.at(histo));
    tot->Add(&histos.at(histo));
  }

  const double the_max(std::max(s.GetMaximum(), histos.at(0).GetMaximum()));
  s.SetMaximum(the_max);
  histos.at(0).SetMaximum(the_max);

  TH1D rat(histos.at(0));
  rat.SetTitle("");
  rat.GetYaxis()->SetTitle("MC/Data Ratio");
  rat.Divide(tot, &histos.at(0));
  rat.SetStats(0);

  rat.SetTitleSize(0.12,"X");
  rat.SetTitleSize(0.12,"Y");
  rat.SetLabelSize(0.1,"X");
  rat.SetLabelSize(0.1,"Y");
  rat.SetTitleOffset(0.32,"Y");

  TLegend l(0.7, 0.5, 0.95, 1.0-0.15/(1.0-div));
  for(unsigned histo(0); histo<histos.size(); ++histo){
    std::string marker_string(histo==0?"lpe":"f");
    l.AddEntry(&histos.at(histo), names.at(histo).c_str(), marker_string.c_str());
  }

  TCanvas c("c","c");
  c.Divide(2);
  c.cd(1);
  c.GetPad(1)->SetPad(0.0,div,1.0,1.0);
  c.GetPad(1)->SetMargin(0.15,0.05,0.0,0.15/(1.0-div));
  c.GetPad(1)->SetLogy(1);
  s.Draw("hist");
  histos.at(0).Draw("samee1p");
  l.Draw("same");
  c.cd(2);
  c.GetPad(2)->SetPad(0.0,0.0,1.0,div);
  c.GetPad(2)->SetMargin(0.15,0.05,0.15/div,0.0);
  c.GetPad(2)->SetGridy(1);
  c.GetPad(2)->SetLogy(1);
  rat.Draw("e1p");
  c.cd(0);
  c.Print(out_name.c_str());
}

int main(){
  SetStyle();
  TH1::SetDefaultSumw2();

  TChain data("data","data");
  TChain ttbar("ttbar","ttbar");
  TChain qcd("qcd","qcd");
  TChain single_t_or_boson("single_t_or_boson","single_t_or_boson");
  TChain diboson("diboson","diboson");

  data.Add("reduced_trees/MET_*2012*SyncSkim.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_FullLept*SyncSkim.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_SemiLept*SyncSkim.root/reduced_tree");
  qcd.Add("reduced_trees/TTJets_Hadronic*SyncSkim.root/reduced_tree");
  qcd.Add("reduced_trees/BJets*SyncSkim.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*channel*SyncSkim.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*JetsTo*SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/WH*SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/WW*SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/WZ*SyncSkim.root/reduced_tree");
  diboson.Add("reduced_trees/ZZ*SyncSkim.root/reduced_tree");
  std::vector<TChain*> chains(0);
  chains.push_back(&data);
  chains.push_back(&qcd);
  chains.push_back(&diboson);
  chains.push_back(&single_t_or_boson);
  chains.push_back(&ttbar);

  std::vector<TH1D> h_pv_nopileup(0), h_pv_notoppt(0), h_pv(0);
  std::vector<TH1D> h_btags_nopileup(0), h_btags_notoppt(0), h_btags(0);
  std::vector<TH1D> h_mbb_nopileup(0), h_mbb_notoppt(0), h_mbb(0);
  std::vector<TH1D> h_n_jet(0);
  unsigned num_bins(6);
  for(unsigned histo(0); histo<num_bins; ++histo){
    std::ostringstream oss("");
    oss << histo;
    if(histo==num_bins-1){
      oss << "+";
    }
    h_n_jet.push_back(TH1D(("h_"+oss.str()+"_jet").c_str(),
			   ("Good Vertices in "+oss.str()+" b-Tag Events (baseline, w/o num. jets cut));Vertices;Events/1").c_str(),
			   61, -0.5, 60.5));
    h_n_jet.at(histo).SetLineColor(histo+1);
  }

  std::vector<std::string> names(0);
  for(unsigned histo(0); histo<chains.size(); ++histo){
    const std::string name(chains.at(histo)->GetName());
    h_pv_nopileup.push_back(TH1D(("h_"+name+"_pv_nopileup").c_str(),"PVs (no pileup reweighting, baseline selection);Vertices;Events/1", 61, -0.5, 60.5));
    h_pv_notoppt.push_back(TH1D(("h_"+name+"_pv_notoppt").c_str(),"PVs (no top pt-based reweighting, baseline selection);Vertices;Events/1", 61, -0.5, 60.5));
    h_pv.push_back(TH1D(("h_"+name+"_pv").c_str(),"PVs (baseline selection);Vertices;Events/1", 61, -0.5, 60.5));
    h_btags_nopileup.push_back(TH1D(("h_"+name+"_btags_nopileup").c_str(),"b-tags (no pileup reweighting, baseline selection);b-tags;Events/1", 9, -0.5, 8.5));
    h_btags_notoppt.push_back(TH1D(("h_"+name+"_btags_notoppt").c_str(),"b-tags (no top pt-based reweighting, baseline selection);b-tags;Events/1", 9, -0.5, 8.5));
    h_btags.push_back(TH1D(("h_"+name+"_btags").c_str(),"b-tags (baseline selection);b-tags;Events/1", 9, -0.5, 8.5));
    h_mbb_nopileup.push_back(TH1D(("h_"+name+"_mbb_nopileup").c_str(),"<m_{bb}> (no pileup reweighting, baseline selection);m_{bb} [GeV];Events/5 GeV", 50, 0.0, 250.0));
    h_mbb_notoppt.push_back(TH1D(("h_"+name+"_mbb_notoppt").c_str(),"<m_{bb}> (no top pt-based reweighting, baseline selection);m_{bb} [GeV];Events/5 GeV", 50, 0.0, 250.0));
    h_mbb.push_back(TH1D(("h_"+name+"_mbb").c_str(),"<m_{bb}> (baseline selection);m_{bb} [GeV];Events/5 GeV", 50, 0.0, 250.0));
    names.push_back(name);
  }

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    passes2CSVTCut(false), passesMETSig30Cut(false), passesMETCleaningCut(false),
    passesTriggerCut(false), passesNumJetsCut(false), passesMinDeltaPhiCut(false),
    passesLeptonVetoCut(false), passesIsoTrackVetoCut(false), passesDRCut(false),
    passesBTaggingCut(false), passesHiggsAvgMassCut(false),
    passesHiggsMassDiffCut(false), passesQCDTriggerCut(false);

  bool passesTriggerPlateauCuts(false), passesBaselineSelection(false);

  float pu_true_num_interactions(0.0);
  unsigned short num_primary_vertices(0);

  float full_weight(0.0), lumi_weight(0.0), top_pt_weight(0.0), pu_weight(0.0), trigger_weight(0.0);
  float average_higgs_mass(0.0);

  unsigned short num_b_tagged_jets(0);

  for(unsigned chain(0); chain<chains.size(); ++chain){
    chains.at(chain)->SetBranchStatus("*",1);
    chains.at(chain)->SetBranchStatus("local_trigger_name",0);
    chains.at(chain)->SetBranchStatus("local_trigger_decision",0);
    chains.at(chain)->SetBranchAddress("passesTriggerPlateauCuts",&passesTriggerPlateauCuts);
    chains.at(chain)->SetBranchAddress("passesBaselineSelection",&passesBaselineSelection);
    chains.at(chain)->SetBranchAddress("passesJSONCut", &passesJSONCut);
    chains.at(chain)->SetBranchAddress("passesPVCut", &passesPVCut);
    chains.at(chain)->SetBranchAddress("passesJet2PtCut", &passesJet2PtCut);
    chains.at(chain)->SetBranchAddress("passes2CSVTCut", &passes2CSVTCut);
    chains.at(chain)->SetBranchAddress("passesMETSig30Cut", &passesMETSig30Cut);
    chains.at(chain)->SetBranchAddress("passesMETCleaningCut", &passesMETCleaningCut);
    chains.at(chain)->SetBranchAddress("passesTriggerCut", &passesTriggerCut);
    chains.at(chain)->SetBranchAddress("passesNumJetsCut", &passesNumJetsCut);
    chains.at(chain)->SetBranchAddress("passesMinDeltaPhiCut", &passesMinDeltaPhiCut);
    chains.at(chain)->SetBranchAddress("passesLeptonVetoCut", &passesLeptonVetoCut);
    chains.at(chain)->SetBranchAddress("passesIsoTrackVetoCut", &passesIsoTrackVetoCut);
    chains.at(chain)->SetBranchAddress("passesDRCut", &passesDRCut);
    chains.at(chain)->SetBranchAddress("passesBTaggingCut", &passesBTaggingCut);
    chains.at(chain)->SetBranchAddress("passesHiggsAvgMassCut", &passesHiggsAvgMassCut);
    chains.at(chain)->SetBranchAddress("passesHiggsMassDiffCut", &passesHiggsMassDiffCut);
    chains.at(chain)->SetBranchAddress("passesQCDTriggerCut", &passesQCDTriggerCut);

    chains.at(chain)->SetBranchAddress("pu_true_num_interactions",&pu_true_num_interactions);
    chains.at(chain)->SetBranchAddress("average_higgs_mass",&average_higgs_mass);
    chains.at(chain)->SetBranchAddress("num_primary_vertices",&num_primary_vertices);
    chains.at(chain)->SetBranchAddress("full_weight", &full_weight);
    chains.at(chain)->SetBranchAddress("lumi_weight", &lumi_weight);
    chains.at(chain)->SetBranchAddress("pu_weight", &pu_weight);
    chains.at(chain)->SetBranchAddress("top_pt_weight", &top_pt_weight);
    chains.at(chain)->SetBranchAddress("trigger_weight", &trigger_weight);
    chains.at(chain)->SetBranchAddress("num_b_tagged_jets", &num_b_tagged_jets);
    int num_entries(chains.at(chain)->GetEntries());

    for(int entry(0); entry<num_entries; ++entry){
      chains.at(chain)->GetEntry(entry);

      if(passesBaselineSelection){
	const double no_pu_weight(lumi_weight*top_pt_weight*trigger_weight);
	const double no_toppt_weight(lumi_weight*pu_weight*trigger_weight);
	h_pv_nopileup.at(chain).Fill(num_primary_vertices, no_pu_weight);
	h_pv_notoppt.at(chain).Fill(num_primary_vertices, no_toppt_weight);
	h_pv.at(chain).Fill(num_primary_vertices, full_weight);
	h_btags_nopileup.at(chain).Fill(num_b_tagged_jets, no_pu_weight);
	h_btags_notoppt.at(chain).Fill(num_b_tagged_jets, no_toppt_weight);
	h_btags.at(chain).Fill(num_b_tagged_jets, full_weight);
	h_mbb_nopileup.at(chain).Fill(average_higgs_mass, no_pu_weight);
	h_mbb_notoppt.at(chain).Fill(average_higgs_mass, no_toppt_weight);
	h_mbb.at(chain).Fill(average_higgs_mass, full_weight);
      }
      if(chain==0){
	if(passesJSONCut && passesPVCut && passesJet2PtCut && passesMETSig30Cut
	   && passesMETCleaningCut && passesTriggerCut && passesMinDeltaPhiCut
	   && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut){
	  unsigned num_b_tags(num_b_tagged_jets);
	  if(num_b_tags>=num_bins) num_b_tags=num_bins-1;
	  h_n_jet.at(num_b_tags).Fill(num_primary_vertices,full_weight);
	}
      }
    }
  }
  MakeRatioPlot(h_pv, names, "pv.pdf");
  MakeRatioPlot(h_pv_nopileup, names, "pv_nopileup.pdf");
  MakeRatioPlot(h_pv_notoppt, names, "pv_notoppt.pdf");
  MakeRatioPlot(h_btags, names, "btags.pdf");
  MakeRatioPlot(h_btags_nopileup, names, "btags_nopileup.pdf");
  MakeRatioPlot(h_btags_notoppt, names, "btags_notoppt.pdf");
  MakeRatioPlot(h_mbb, names, "mbb.pdf");
  MakeRatioPlot(h_mbb_nopileup, names, "mbb_nopileup.pdf");
  MakeRatioPlot(h_mbb_notoppt, names, "mbb_notoppt.pdf");
}
