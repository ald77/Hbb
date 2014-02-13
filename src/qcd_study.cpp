#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "style.hpp"
#include "timer.hpp"
#include "utils.hpp"

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

  data.Add("reduced_trees/MET_*2012*1.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_FullLept*1.root/reduced_tree");
  ttbar.Add("reduced_trees/TTJets_SemiLept*1.root/reduced_tree");
  qcd.Add("reduced_trees/TTJets_Hadronic*1.root/reduced_tree");
  qcd.Add("reduced_trees/BJets*1.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*channel*1.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/*JetsToLNu_Tune*1.root/reduced_tree");
  single_t_or_boson.Add("reduced_trees/ZJetsToNuNu*1.root/reduced_tree");
  diboson.Add("reduced_trees/WH*1.root/reduced_tree");
  diboson.Add("reduced_trees/WW*1.root/reduced_tree");
  diboson.Add("reduced_trees/WZ*1.root/reduced_tree");
  diboson.Add("reduced_trees/ZZ*1.root/reduced_tree");
  std::vector<TChain*> chains(0);
  chains.push_back(&data);
  chains.push_back(&qcd);
  chains.push_back(&diboson);
  chains.push_back(&single_t_or_boson);
  chains.push_back(&ttbar);

  std::vector<std::string> names(0);

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    /*passesMETSig30Cut(false),*/ passesMETCleaningCut(false), passesNumJetsCut(false),
    passesMinDeltaPhiCut(false), passesLeptonVetoCut(false),
    passesIsoTrackVetoCut(false), passesDRCut(false), passesHiggsAvgMassCut(false),
    passesHiggsMassDiffCut(false), passesQCDTriggerCut(false);

  float met(0.0), met_sig(0.0), min_delta_phi(0.0), average_higgs_mass(0.0),
    second_highest_jet_pt(0.0), higgs_mass_difference(0.0), full_weight(0.0);

  unsigned short num_b_tagged_jets(0);

  std::vector<TH1D> h_min_delta_phi(0), h_average_higgs_mass(0), h_higgs_mass_diff(0),
    h_metsig(0);
  std::vector<TH1D> h_min_delta_phi_sbin0(0), h_min_delta_phi_sbin1(0),
    h_min_delta_phi_sbin2(0), h_min_delta_phi_sbin3(0), h_min_delta_phi_sbin4(0);

  for(unsigned chain_num(0); chain_num<chains.size(); ++chain_num){
    std::cout << chain_num << "/" << chains.size() << std::endl;
    TChain& chain(*chains.at(chain_num));
    chain.SetBranchStatus("*",0);
    setup(chain, "passesJSONCut", passesJSONCut);
    setup(chain, "passesPVCut", passesPVCut);
    setup(chain, "passesJet2PtCut", passesJet2PtCut);
    setup(chain, "passesMETCleaningCut", passesMETCleaningCut);
    setup(chain, "passesNumJetsCut", passesNumJetsCut);
    setup(chain, "passesLeptonVetoCut", passesLeptonVetoCut);
    setup(chain, "passesIsoTrackVetoCut", passesIsoTrackVetoCut);
    setup(chain, "passesDRCut", passesDRCut);
    setup(chain, "passesHiggsAvgMassCut", passesHiggsAvgMassCut);
    setup(chain, "passesHiggsMassDiffCut", passesHiggsMassDiffCut);
    setup(chain, "passesQCDTriggerCut", passesQCDTriggerCut);
    setup(chain, "met", met);
    setup(chain, "second_highest_jet_pt", second_highest_jet_pt);
    setup(chain, "met_sig", met_sig);
    setup(chain, "min_delta_phi", min_delta_phi);
    setup(chain, "average_higgs_mass", average_higgs_mass);
    setup(chain, "higgs_mass_difference", higgs_mass_difference);
    setup(chain, "full_weight", full_weight);

    const std::string name(chain.GetName());
    names.push_back(name);
    
    h_min_delta_phi.push_back(TH1D(("h_min_delta_phi_"+name).c_str(),
				   "Min. Delta Phi (QCD Control);Min Delta Phi;Events/0.1",
				   32, 0.0, 3.2));
    h_average_higgs_mass.push_back(TH1D(("h_average_higgs_mass_"+name).c_str(),
					"<m_{bb}> (QCD Control);<m_{bb}> [GeV];Events/5 GeV",
					50, 0.0, 250.0));
    h_higgs_mass_diff.push_back(TH1D(("h_higgs_mass_diff_"+name).c_str(),
				     "#Delta m_{bb} (QCD Control);#Delta m_{bb} [GeV];Events/1 GeV",
				     50, 0.0, 50.0));
    h_metsig.push_back(TH1D(("h_metsig_"+name).c_str(),
			    "S_{MET} (QCD Control);S_{MET};Events/10", 40, 0.0, 400.0));
    h_min_delta_phi_sbin0.push_back(TH1D(("h_min_delta_phi_sbin0_"+name).c_str(),
					 "Min. Delta Phi (S-bin 0) (QCD Control);Min Delta Phi;Events/0.1",
					 32, 0.0, 3.2));
    h_min_delta_phi_sbin1.push_back(TH1D(("h_min_delta_phi_sbin1_"+name).c_str(),
					 "Min. Delta Phi (S-bin 1) (QCD Control);Min Delta Phi;Events/0.1",
					 32, 0.0, 3.2));
    h_min_delta_phi_sbin2.push_back(TH1D(("h_min_delta_phi_sbin2_"+name).c_str(),
					 "Min. Delta Phi (S-bin 2) (QCD Control);Min Delta Phi;Events/0.1",
					 32, 0.0, 3.2));
    h_min_delta_phi_sbin3.push_back(TH1D(("h_min_delta_phi_sbin3_"+name).c_str(),
					 "Min. Delta Phi (S-bin 3) (QCD Control);Min Delta Phi;Events/0.1",
					 32, 0.0, 3.2));
    h_min_delta_phi_sbin4.push_back(TH1D(("h_min_delta_phi_sbin4_"+name).c_str(),
					 "Min. Delta Phi (S-bin 4) (QCD Control);Min Delta Phi;Events/0.1",
					 32, 0.0, 3.2));
    
    const int num_events(chain.GetEntries());
    Timer timer(num_events);
    timer.Start();
    for(int event(0); event<num_events; ++event){
      if(event%(1u<<16u)==0){
	timer.PrintRemainingTime();
      }
      chain.GetEntry(event);
      const double stupid_factor(chain_num==0?19307.0/5208.0:1.0);
      full_weight*=stupid_factor;
      if(passesJSONCut && passesPVCut && passesJet2PtCut && passesMETCleaningCut
	 && passesNumJetsCut && passesLeptonVetoCut && passesIsoTrackVetoCut
	 && /*passesDRCut &&*/ passesQCDTriggerCut && num_b_tagged_jets<3 && met>160.0 && second_highest_jet_pt>70.0){
	if(/*passesHiggsMassDiffCut && !passesMETSig30Cut &&*/ !passesMinDeltaPhiCut){
	  h_average_higgs_mass.at(chain_num).Fill(average_higgs_mass, full_weight);
	}
	if(/*passesHiggsAvgMassCut && !passesMETSig30Cut && */!passesMinDeltaPhiCut){
	  h_higgs_mass_diff.at(chain_num).Fill(higgs_mass_difference, full_weight);
	}
	if(/*passesHiggsAvgMassCut && passesHiggsMassDiffCut
	     && */!passesMinDeltaPhiCut){
	  h_metsig.at(chain_num).Fill(met_sig, full_weight);
	}
	if(/*passesHiggsAvgMassCut && passesHiggsMassDiffCut && !passesMETSig30Cut*/ true){
	  h_min_delta_phi.at(chain_num).Fill(min_delta_phi, full_weight);
	}
	if(/*passesHiggsAvgMassCut && passesHiggsMassDiffCut*/ true){
	  if(met_sig<30.0){
	    h_min_delta_phi_sbin0.at(chain_num).Fill(min_delta_phi, full_weight);
	  }else if(met_sig<50.0){
	    h_min_delta_phi_sbin1.at(chain_num).Fill(min_delta_phi, full_weight);
	  }else if(met_sig<100.0){
	    h_min_delta_phi_sbin2.at(chain_num).Fill(min_delta_phi, full_weight);
	  }else if(met_sig<150.0){
	    h_min_delta_phi_sbin3.at(chain_num).Fill(min_delta_phi, full_weight);
	  }else{
	    h_min_delta_phi_sbin4.at(chain_num).Fill(min_delta_phi, full_weight);
	  }
	}
      }
      timer.Iterate();
    }
  }

  MakeRatioPlot(h_min_delta_phi, names, "qcd/h_min_delta_phi.pdf");
  MakeRatioPlot(h_average_higgs_mass, names, "qcd/h_average_higgs_mass.pdf");
  MakeRatioPlot(h_higgs_mass_diff, names, "qcd/h_higgs_mass_diff.pdf");
  MakeRatioPlot(h_metsig, names, "qcd/h_metsig.pdf");
  MakeRatioPlot(h_min_delta_phi_sbin0, names, "qcd/h_min_delta_phi_sbin0.pdf");
  MakeRatioPlot(h_min_delta_phi_sbin1, names, "qcd/h_min_delta_phi_sbin1.pdf");
  MakeRatioPlot(h_min_delta_phi_sbin2, names, "qcd/h_min_delta_phi_sbin2.pdf");
  MakeRatioPlot(h_min_delta_phi_sbin3, names, "qcd/h_min_delta_phi_sbin3.pdf");
  MakeRatioPlot(h_min_delta_phi_sbin4, names, "qcd/h_min_delta_phi_sbin4.pdf");
  double rat_tot(0.0), uncert_tot(0.0);
  double rat_sbin0(0.0), uncert_sbin0(0.0);
  double rat_sbin1(0.0), uncert_sbin1(0.0);
  double rat_sbin2(0.0), uncert_sbin2(0.0);
  double rat_sbin3(0.0), uncert_sbin3(0.0);
  double rat_sbin4(0.0), uncert_sbin4(0.0);
  const double low(0.3), high(0.3);
  GetHighLowRatio(h_min_delta_phi.at(0), low, high, rat_tot, uncert_tot);
  GetHighLowRatio(h_min_delta_phi_sbin0.at(0), low, high, rat_sbin0, uncert_sbin0);
  GetHighLowRatio(h_min_delta_phi_sbin1.at(0), low, high, rat_sbin1, uncert_sbin1);
  GetHighLowRatio(h_min_delta_phi_sbin2.at(0), low, high, rat_sbin2, uncert_sbin2);
  GetHighLowRatio(h_min_delta_phi_sbin3.at(0), low, high, rat_sbin3, uncert_sbin3);
  GetHighLowRatio(h_min_delta_phi_sbin4.at(0), low, high, rat_sbin4, uncert_sbin4);
  double mc_rat_tot(0.0), mc_uncert_tot(0.0);
  double mc_rat_sbin0(0.0), mc_uncert_sbin0(0.0);
  double mc_rat_sbin1(0.0), mc_uncert_sbin1(0.0);
  double mc_rat_sbin2(0.0), mc_uncert_sbin2(0.0);
  double mc_rat_sbin3(0.0), mc_uncert_sbin3(0.0);
  double mc_rat_sbin4(0.0), mc_uncert_sbin4(0.0);
  GetHighLowRatio(h_min_delta_phi.at(1), low, high, mc_rat_tot, mc_uncert_tot);
  GetHighLowRatio(h_min_delta_phi_sbin0.at(1), low, high,
		  mc_rat_sbin0, mc_uncert_sbin0);
  GetHighLowRatio(h_min_delta_phi_sbin1.at(1), low, high,
		  mc_rat_sbin1, mc_uncert_sbin1);
  GetHighLowRatio(h_min_delta_phi_sbin2.at(1), low, high,
		  mc_rat_sbin2, mc_uncert_sbin2);
  GetHighLowRatio(h_min_delta_phi_sbin3.at(1), low, high,
		  mc_rat_sbin3, mc_uncert_sbin3);
  GetHighLowRatio(h_min_delta_phi_sbin4.at(1), low, high,
		  mc_rat_sbin4, mc_uncert_sbin4);

  std::cout << rat_tot << " " << uncert_tot << std::endl;
  std::cout << rat_sbin0 << " " << uncert_sbin0 << std::endl;
  std::cout << rat_sbin1 << " " << uncert_sbin1 << std::endl;
  std::cout << rat_sbin2 << " " << uncert_sbin2 << std::endl;
  std::cout << rat_sbin3 << " " << uncert_sbin3 << std::endl;
  std::cout << rat_sbin4 << " " << uncert_sbin4 << std::endl;
  std::cout << mc_rat_tot << " " << mc_uncert_tot << std::endl;
  std::cout << mc_rat_sbin0 << " " << mc_uncert_sbin0 << std::endl;
  std::cout << mc_rat_sbin1 << " " << mc_uncert_sbin1 << std::endl;
  std::cout << mc_rat_sbin2 << " " << mc_uncert_sbin2 << std::endl;
  std::cout << mc_rat_sbin3 << " " << mc_uncert_sbin3 << std::endl;
  std::cout << mc_rat_sbin4 << " " << mc_uncert_sbin4 << std::endl;


}
