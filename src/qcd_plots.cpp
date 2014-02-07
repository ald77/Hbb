#include <cmath>
#include <iostream>
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "style.hpp"
#include "timer.hpp"

void get_independence_model(const TH2D& in, TH2D& out){
  out=in;
  const unsigned x_bins(in.GetNbinsX()), y_bins(in.GetNbinsY());

  TH1D* x_proj(in.ProjectionX("_px",0,-1,"e"));
  TH1D* y_proj(in.ProjectionY("_py",0,-1,"e"));
  const double integral(in.Integral(0, -1, 0, -1));
  for(unsigned x_bin(0); x_bin<=x_bins+1; ++x_bin){
    for(unsigned y_bin(0); y_bin<=y_bins+1; ++y_bin){
      const double x_val(x_proj->GetBinContent(x_bin));
      const double y_val(y_proj->GetBinContent(y_bin));
      const double x_err(x_proj->GetBinError(x_bin));
      const double y_err(y_proj->GetBinError(y_bin));
      out.SetBinContent(x_bin, y_bin, x_val*y_val/integral);
      out.SetBinError(x_bin, y_bin,
		      sqrt(x_val*x_val*y_err*y_err+y_val*y_val*x_err*x_err)/integral);
    }
  }
}

void draw(TH2D& h, const std::string& name){
  TCanvas c;
  h.SetStats(0);
  h.Draw("colz");
  h.Draw("textesame");
  c.Print(name.c_str());
}

template<typename T>
void setup(TChain& chain, const std::string& name, T& variable){
  chain.SetBranchStatus(name.c_str(), 1);
  chain.SetBranchAddress(name.c_str(), &variable);
}

int main(){
  SetStyle();
  TH1::SetDefaultSumw2();

  TChain qcd("qcd","qcd");
  qcd.Add("reduced_trees/QCD*1.root/reduced_tree");

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    passesMETCleaningCut(false), passesNumJetsCut(false), passesLeptonVetoCut(false),
    passesIsoTrackVetoCut(false), passesDRCut(false), passesQCDTriggerCut(false);

  unsigned short num_b_tagged_jets(0);

  float met_sig(0.0);
  float average_higgs_mass(0.0);
  float full_weight(0.0);

  TH2D btags_metsig("btags_metsig","b-Tagged Jets vs. SMET Bin (baseline w/o MDP, SMET, b-Tagging);SMET Bin;b-Tags", 4, 0.5, 4.5, 9, -0.5, 8.5);
  TH2D btags_mbb("btags_mbb","b-Tagged Jets vs. m_{bb} (baseline w/o MDP, SMET, b-Tagging);<m_{bb}> [GeV];b-Tags", 6, 0.0, 300.0, 9, -0.5, 8.5);
  TH2D metsig_mbb("metsig_mbb","SMET Bin vs. m_{bb} (baseline w/o MDP, SMET, b-Tagging);<m_{bb}> [GeV];SMET Bin", 6, 0.0, 300.0, 4, 0.5, 4.5);

  qcd.SetBranchStatus("*",0);
  setup(qcd, "passesJSONCut", passesJSONCut);
  setup(qcd, "passesPVCut", passesPVCut);
  setup(qcd, "passesJet2PtCut", passesJet2PtCut);
  setup(qcd, "passesMETCleaningCut", passesMETCleaningCut);
  setup(qcd, "passesNumJetsCut", passesNumJetsCut);
  setup(qcd, "passesLeptonVetoCut", passesLeptonVetoCut);
  setup(qcd, "passesIsoTrackVetoCut", passesIsoTrackVetoCut);
  setup(qcd, "passesDRCut", passesDRCut);
  setup(qcd, "passesQCDTriggerCut", passesQCDTriggerCut);
  setup(qcd, "met_sig", met_sig);
  setup(qcd, "full_weight", full_weight);
  setup(qcd, "num_b_tagged_jets", num_b_tagged_jets);
  setup(qcd, "average_higgs_mass", average_higgs_mass);

  const int num_entries(qcd.GetEntries());
  Timer timer(num_entries);
  timer.Start();
  for(int entry(0); entry<num_entries; ++entry){
    if(entry%(1u << 16u)==0){
      timer.PrintRemainingTime();
    }
    qcd.GetEntry(entry);
    if(passesJSONCut && passesPVCut && passesJet2PtCut && passesMETCleaningCut
       && (true || passesQCDTriggerCut) && passesNumJetsCut && passesLeptonVetoCut
       && passesIsoTrackVetoCut && passesDRCut){
      int smet_bin(0);
      if(met_sig<0.0){
	smet_bin=-1;
      }else if(met_sig<30.0){
	smet_bin=0;
      }else if(met_sig<50.0){
	smet_bin=1;
      }else if(met_sig<100.0){
	smet_bin=2;
      }else if(met_sig<150.0){
	smet_bin=3;
      }else{
	smet_bin=4;
      }
      btags_metsig.Fill(smet_bin, num_b_tagged_jets, full_weight);
      btags_mbb.Fill(average_higgs_mass, num_b_tagged_jets, full_weight);
      metsig_mbb.Fill(average_higgs_mass, smet_bin, full_weight);
    }
    timer.Iterate();
  }
  TH2D *btags_metsig_indep(static_cast<TH2D*>(btags_metsig.Clone("btags_metsig_indep")));
  TH2D *btags_mbb_indep(static_cast<TH2D*>(btags_mbb.Clone("btags_mbb_indep")));
  TH2D *metsig_mbb_indep(static_cast<TH2D*>(metsig_mbb.Clone("metsig_mbb_indep")));
  get_independence_model(btags_metsig, *btags_metsig_indep);
  get_independence_model(btags_mbb, *btags_mbb_indep);
  get_independence_model(metsig_mbb, *metsig_mbb_indep);

  TCanvas c;
  draw(btags_metsig, "normal_btags_metsig.pdf");
  draw(*btags_metsig_indep, "indep_btags_metsig.pdf");
  draw(btags_mbb, "normal_btags_mbb.pdf");
  draw(*btags_mbb_indep, "indep_btags_mbb.pdf");
  draw(metsig_mbb, "normal_metsig_mbb.pdf");
  draw(*metsig_mbb_indep, "indep_metsig_mbb.pdf");
}
