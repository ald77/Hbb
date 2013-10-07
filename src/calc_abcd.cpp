#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "weights.hpp"

const unsigned int max_draws(10000);
TRandom3 rng(0);

void GetFiles(std::vector<TFile*> &files, const std::vector<std::string> &names){
  for(std::vector<TFile*>::size_type file(0); file<files.size(); ++file){
    if(files.at(file)!=NULL && files.at(file)->IsOpen() && !files.at(file)->IsZombie()){
      files.at(file)->Close();
    }   
    delete files.at(file);
    files.at(file)=NULL;
  }
  files.clear();
  files.resize(names.size(), NULL);
  for(std::vector<std::string>::size_type name(0); name<names.size(); ++name){
    files.at(name)=new TFile(names.at(name).c_str(),"read");
  }
}

void GetTrees(std::vector<TTree*> &trees, const std::vector<TFile*> &files){
  for(std::vector<TTree*>::size_type tree(0); tree<trees.size(); ++tree){
    if(trees.at(tree)!=NULL){
      delete trees.at(tree);
      trees.at(tree)=NULL;
    }
  }
  trees.clear();
  trees.resize(files.size(), NULL);
  for(std::vector<TFile*>::size_type file(0); file<files.size(); ++file){
    if(files.at(file)!=NULL && files.at(file)->IsOpen() && !files.at(file)->IsZombie()){
      files.at(file)->GetObject("ABCD",trees.at(file));
    }
  }
}

void GetValues(std::vector<double> &A, std::vector<double> &B, std::vector<double> &C3,
	       std::vector<double> &D3, std::vector<double> &C2, std::vector<double> &D2,
	       const std::vector<TTree*> &trees){
  A.clear();
  B.clear();
  C3.clear();
  D3.clear();
  C2.clear();
  D2.clear();
  A.resize(trees.size(),0.0);
  B.resize(trees.size(),0.0);
  C3.resize(trees.size(),0.0);
  D3.resize(trees.size(),0.0);
  C2.resize(trees.size(),0.0);
  D2.resize(trees.size(),0.0);
  for(std::vector<TTree*>::size_type tree(0); tree<trees.size(); ++tree){
    if(trees.at(tree)!=NULL && trees.at(tree)->GetEntries()>0){
      trees.at(tree)->SetBranchAddress("AWeighted",&A.at(tree));
      trees.at(tree)->SetBranchAddress("BWeighted",&B.at(tree));
      trees.at(tree)->SetBranchAddress("C3bWeighted",&C3.at(tree));
      trees.at(tree)->SetBranchAddress("D3bWeighted",&D3.at(tree));
      trees.at(tree)->SetBranchAddress("C2bWeighted",&C2.at(tree));
      trees.at(tree)->SetBranchAddress("D2bWeighted",&D2.at(tree));
      
      trees.at(tree)->GetEntry(0);
    }
  }
}

void GetWeights(std::vector<double> &weights, const std::vector<std::string> &names,
		const WeightCalculator &weightCalc){
  weights.clear();
  weights.resize(names.size(), 0.0);
  for(std::vector<std::string>::size_type name(0); name<names.size(); ++name){
    weights.at(name)=weightCalc.GetWeight(names.at(name));
  }
}

void KillTrees(std::vector<TTree*> &trees){
  for(std::vector<TTree*>::size_type tree(0); tree<trees.size(); ++tree){
    if(trees.at(tree)!=NULL){
      delete trees.at(tree);
      trees.at(tree)=NULL;
    }
  }
}

void KillFiles(std::vector<TFile*> &files){
  for(std::vector<TFile*>::size_type file(0); file<files.size(); ++file){
    if(files.at(file)!=NULL && files.at(file)->IsOpen()){
      files.at(file)->Close();
    }
    delete files.at(file);
    files.at(file)=NULL;
  }
}

double GetRandom(const double kp1, const double weight){
  const double quantile(rng.Rndm());
  double left(0.0), right(DBL_MAX);
  double middle(left+0.5*(right-left));
  while(left<middle && middle<right){
    if(TMath::Gamma(kp1, middle)<quantile){
      left=middle;
    }else{
      right=middle;
    }
    middle=left+0.5*(right-left);
  }
  return middle*weight;
}

void GetUncertaintiesOld(double &up, double &down, const double center,
		      const std::vector<double> vals){
  if(vals.size()!=0){
    const double frac(TMath::Erf(1.0/sqrt(2.0)));
    const std::vector<double>::size_type delta(static_cast<int>(ceil(frac*(vals.size()-1))));
    double min_diff(DBL_MAX);
    double left(0.0), right(0.0);
    for(std::vector<double>::size_type val(0); val+delta<vals.size(); ++val){
      std::vector<double>::size_type val2(val+delta);
      const double diff(vals.at(val2)-vals.at(val));
      if(diff<min_diff){
	min_diff=diff;
	right=vals.at(val2);
	left=vals.at(val);
      }
    }
    if(center<left || center>right){
      std::cout << left << " " << center << " " << right << std::endl;
    }
    up=right-center;
    down=center-left;
    if(up<0.0) up=0.0;
    if(down<0.0) down=0.0;
  }else{
    up=0.0;
    down=0.0;
  }
}

void GetUncertainties(double &up, double &down, const double center,
		      const std::vector<double> vals){
  unsigned int best_pos(-1);
  double best_diff(DBL_MAX);
  for(unsigned int pos(0); pos<vals.size(); ++pos){
    const double diff(fabs(center-vals.at(pos)));
    if(diff<best_diff){
      best_diff=diff;
      best_pos=pos;
    }
  }
  const double frac(TMath::Erf(1.0/sqrt(2.0)));
  const std::vector<double>::size_type delta(static_cast<int>(0.5*ceil(frac*(vals.size()-1))));
  unsigned int upper(best_pos+delta), lower(best_pos-delta);
  if(lower>best_pos) lower=best_pos;
  if(upper>=vals.size()) upper=vals.size()-1;
  up=vals.at(upper)-vals.at(best_pos);
  down=vals.at(best_pos)-vals.at(lower);
}

void GetSummedVals(double &up, double &down, double &center,
		   const std::vector<double> &val,
		   const std::vector<double> &weights,
		   const std::vector<double>::size_type &low,
		   const std::vector<double>::size_type &high){
  center=0.0;
  for(std::vector<double>::size_type sample(low); sample<high; ++sample){
    center+=val.at(sample);
  }
  std::vector<double> draws(0);
  for(unsigned int draw(0); draw<max_draws; ++draw){
    double thisVal(0.0);
    for(std::vector<double>::size_type sample(low); sample<high; ++sample){
      thisVal+=GetRandom(val.at(sample)/weights.at(sample), weights.at(sample));
    }
    draws.push_back(thisVal);
  }
  std::sort(draws.begin(), draws.end());
  GetUncertainties(up, down, center, draws);
}

void GetSummedKappas(double &up, double &down, double &center, const std::vector<double> &A,
		      const std::vector<double> &B, const std::vector<double> &C,
		     const std::vector<double> &D, const std::vector<double> &weights,
		     const std::vector<double>::size_type &low,
		     const std::vector<double>::size_type &high){
  std::vector<double> draws(0);
  for(unsigned int draw(0); draw<max_draws; ++draw){
    double sumA(0.0), sumB(0.0), sumC(0.0), sumD(0.0);
    for(std::vector<double>::size_type sample(low); sample<A.size() && sample<high; ++sample){
      sumA+=GetRandom(A.at(sample)/weights.at(sample), weights.at(sample));
      sumB+=GetRandom(B.at(sample)/weights.at(sample), weights.at(sample));
      sumC+=GetRandom(C.at(sample)/weights.at(sample), weights.at(sample));
      sumD+=GetRandom(D.at(sample)/weights.at(sample), weights.at(sample));
    }
    if(sumB*sumC>0.0){
      draws.push_back(sumA*sumD/(sumB*sumC));
    }else{
      draws.push_back(DBL_MAX);
    }
  }

  double sumA(0.0), sumB(0.0), sumC(0.0), sumD(0.0);
  for(std::vector<double>::size_type sample(low); sample<high; ++sample){
    sumA+=A.at(sample);
    sumB+=B.at(sample);
    sumC+=C.at(sample);
    sumD+=D.at(sample);
  }
  if(sumB*sumC>0.0){
    center=sumA*sumD/(sumB*sumC);
  }else{
    center=DBL_MAX;
  }
  std::sort(draws.begin(), draws.end());
  GetUncertainties(up, down, center, draws);
}

void GetSummedAPred(double &up, double &down, double &center,
		    const std::vector<double> &B, const std::vector<double> &C,
		    const std::vector<double> &D, const std::vector<double> &weights,
		    const std::vector<double>::size_type &low,
		    const std::vector<double>::size_type &high){
  std::vector<double> draws(0);
  for(unsigned int draw(0); draw<max_draws; ++draw){
    double sumB(0.0), sumC(0.0), sumD(0.0);
    for(std::vector<double>::size_type sample(low); sample<B.size() && sample<high; ++sample){
      sumB+=GetRandom(B.at(sample)/weights.at(sample), weights.at(sample));
      sumC+=GetRandom(C.at(sample)/weights.at(sample), weights.at(sample));
      sumD+=GetRandom(D.at(sample)/weights.at(sample), weights.at(sample));
    }
    if(sumD>0.0){
      draws.push_back(sumB*sumC/sumD);
    }else{
      draws.push_back(DBL_MAX);
    }
  }

  double sumB(0.0), sumC(0.0), sumD(0.0);
  for(std::vector<double>::size_type sample(low); sample<high; ++sample){
    sumB+=B.at(sample);
    sumC+=C.at(sample);
    sumD+=D.at(sample);
  }
  if(sumD>0.0){
    center=sumB*sumC/sumD;
  }else{
    center=DBL_MAX;
  }
  std::sort(draws.begin(), draws.end());
  GetUncertainties(up, down, center, draws);
}

void PrintLine(std::string name, const double A, const double Aup, const double Adown, const double B, const double Bup, const double Bdown, const double C3, const double C3up, const double C3down, const double D3, const double D3up, const double D3down, const double C2, const double C2up, const double C2down, const double D2, const double D2up, const double D2down, const double kappa3, const double kappa3up, const double kappa3down, const double kappa2, const double kappa2up, const double kappa2down, const double pred23, const double predup23, const double preddown23, const double pred24, const double predup24, const double preddown24, const double pred34, const double predup34, const double preddown34){
  std::cout << name << " & $"
	    << std::setprecision(3) << A << "_{-"
	    << std::setprecision(2) << Adown << "}^{+" << Aup << "}$ & $"
	    << std::setprecision(3) << B << "_{-"
	    << std::setprecision(2) << Bdown << "}^{+" << Bup << "}$ & $"
	    << std::setprecision(3) << C3 << "_{-"
	    << std::setprecision(2) << C3down << "}^{+" << C3up << "}$ & $"
	    << std::setprecision(3) << D3 << "_{-"
	    << std::setprecision(2) << D3down << "}^{+" << D3up << "}$ & $"
	    << std::setprecision(3) << C2 << "_{-"
	    << std::setprecision(2) << C2down << "}^{+" << C2up << "}$ & $"
	    << std::setprecision(3) << D2 << "_{-"
	    << std::setprecision(2) << D2down << "}^{+" << D2up << "}$ & $"
	    << std::setprecision(3) << kappa3 << "_{-"
	    << std::setprecision(2) << kappa3down << "}^{+" << kappa3up << "}$ & $"
	    << std::setprecision(3) << kappa2 << "_{-"
	    << std::setprecision(2) << kappa2down << "}^{+" << kappa2up << "}$ & $"
	    << std::setprecision(3) << pred23 << "_{-"
	    << std::setprecision(2) << preddown23 << "}^{+" << predup23 << "}$ & $"
	    << std::setprecision(3) << pred24 << "_{-"
	    << std::setprecision(2) << preddown24 << "}^{+" << predup24 << "}$ & $"
	    << std::setprecision(3) << pred34 << "_{-"
	    << std::setprecision(2) << preddown34 << "}^{+" << predup34 << "}$"
	    << std::endl;
}

int main(){
  {TTree crap;}
  std::vector<std::string> names(0);
  names.push_back("raw_plots_and_values/MET_Run2012A-13Jul2012-v1_AOD_UCSB1852_v71_SyncSkim.root");//0
  names.push_back("raw_plots_and_values/MET_Run2012B-13Jul2012-v1_AOD_UCSB1853_v71_SyncSkim.root");//1
  names.push_back("raw_plots_and_values/MET_Run2012C-24Aug2012-v1_AOD_UCSB1854_v71_SyncSkim.root");//2
  names.push_back("raw_plots_and_values/MET_Run2012C-PromptReco-v2_AOD_UCSB1867_v71_SyncSkim.root");//3
  names.push_back("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1869_v71_SyncSkim.root");//4
  names.push_back("raw_plots_and_values/MET_Run2012D-PromptReco-v1_AOD_UCSB1870_v71_SyncSkim.root");//5
  names.push_back("raw_plots_and_values/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1903_v71_SyncSkim.root");//6
  names.push_back("raw_plots_and_values/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_UCSB1897_v71_SyncSkim.root");//7
  names.push_back("raw_plots_and_values/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1904_v71_SyncSkim.root");//8
  names.push_back("raw_plots_and_values/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1898_v71_SyncSkim.root");//9
  names.push_back("raw_plots_and_values/QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1905_v71_SyncSkim.root");//10
  names.push_back("raw_plots_and_values/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1899_v71_SyncSkim.root");//11
  names.push_back("raw_plots_and_values/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1900_v71_SyncSkim.root");//12
  names.push_back("raw_plots_and_values/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1901_v71_SyncSkim.root");//13
  names.push_back("raw_plots_and_values/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1902_v71_SyncSkim.root");//14
  names.push_back("raw_plots_and_values/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v2_AODSIM_UCSB1883_v71_SyncSkim.root");//15
  names.push_back("raw_plots_and_values/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_SyncSkim.root");//16
  names.push_back("raw_plots_and_values/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1880_v71_SyncSkim.root");//17
  names.push_back("raw_plots_and_values/TTH_HToBB_M-125_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1855_v71_SyncSkim.root");//18
  names.push_back("raw_plots_and_values/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1857_v71_SyncSkim.root");//19
  names.push_back("raw_plots_and_values/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1856_v71_SyncSkim.root");//20
  names.push_back("raw_plots_and_values/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1864_v71_SyncSkim.root");//21
  names.push_back("raw_plots_and_values/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1865_v71_SyncSkim.root");//22
  names.push_back("raw_plots_and_values/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1866_v71_SyncSkim.root");//23
  names.push_back("raw_plots_and_values/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1860_v71_SyncSkim.root");//24
  names.push_back("raw_plots_and_values/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1861_v71_SyncSkim.root");//25
  names.push_back("raw_plots_and_values/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1862_v71_SyncSkim.root");//26
  names.push_back("raw_plots_and_values/W2JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1877_v71_SyncSkim.root");//27
  names.push_back("raw_plots_and_values/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1878_v71_SyncSkim.root");//28
  names.push_back("raw_plots_and_values/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1879_v71_SyncSkim.root");//29
  names.push_back("raw_plots_and_values/ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1887_v71_SyncSkim.root");//30
  names.push_back("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1889_v71_SyncSkim.root");//31
  names.push_back("raw_plots_and_values/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1888_v71_SyncSkim.root");//32
  names.push_back("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1891_v71_SyncSkim.root");//33
  names.push_back("raw_plots_and_values/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1890_v71_SyncSkim.root");//34
  names.push_back("raw_plots_and_values/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1858_v71_SyncSkim.root");//35
  names.push_back("raw_plots_and_values/ZH_ZToBB_HToBB_M-125_8TeV-powheg-herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1868_v71_SyncSkim.root");//36
  names.push_back("raw_plots_and_values/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1874_v71_SyncSkim.root");//37
  names.push_back("raw_plots_and_values/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1876_v71_SyncSkim.root");//38
  names.push_back("raw_plots_and_values/WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1859_v71_SyncSkim.root");//39

  std::vector<TFile*> files(0);
  GetFiles(files, names);

  std::vector<TTree*> trees(0);
  GetTrees(trees, files);

  std::vector<double> A(0), B(0), C3(0), D3(0), C2(0), D2(0);
  GetValues(A, B, C3, D3, C2, D2, trees);

  WeightCalculator weightCalc(19399.0);
  std::vector<double> weights(0);
  GetWeights(weights, names, weightCalc);

  std::vector<std::string> tex(0);
  std::vector<unsigned int> upper(0), lower(0);
  tex.push_back("Data"); upper.push_back(6); lower.push_back(0);
  tex.push_back("QCD"); upper.push_back(15); lower.push_back(6);
  tex.push_back("$t\\overline{t}$ (2l)"); upper.push_back(16); lower.push_back(15);
  tex.push_back("$t\\overline{t}$ (1l)"); upper.push_back(17); lower.push_back(16);
  tex.push_back("$t\\overline{t}$ (0l)"); upper.push_back(18); lower.push_back(17);
  tex.push_back("t\\overline{t}H$"); upper.push_back(19); lower.push_back(18);
  tex.push_back("$t\\overline{t}V$"); upper.push_back(21); lower.push_back(19);
  tex.push_back("$t$"); upper.push_back(27); lower.push_back(21);
  tex.push_back("$V$"); upper.push_back(35); lower.push_back(27);
  tex.push_back("$VH$"); upper.push_back(37); lower.push_back(35);
  tex.push_back("$VV$"); upper.push_back(39); lower.push_back(37);
  tex.push_back("SM total"); upper.push_back(39); lower.push_back(15);
  tex.push_back("SM total (no QCD)"); upper.push_back(39); lower.push_back(15);

  std::vector<double> Aval(tex.size()), Aup(tex.size()), Adown(tex.size()), Bval(tex.size()), Bup(tex.size()), Bdown(tex.size()), C3val(tex.size()),
    C3up(tex.size()), C3down(tex.size()), D3val(tex.size()), D3up(tex.size()), D3down(tex.size()), C2val(tex.size()), C2up(tex.size()), C2down(tex.size()),
    D2val(tex.size()), D2up(tex.size()), D2down(tex.size()), kappa3val(tex.size()), kappa3up(tex.size()), kappa3down(tex.size()),
    kappa2val(tex.size()), kappa2up(tex.size()), kappa2down(tex.size());

  std::vector<double> pred23(tex.size()), predup23(tex.size()), preddown23(tex.size());
  std::vector<double> pred24(tex.size()), predup24(tex.size()), preddown24(tex.size());
  std::vector<double> pred34(tex.size()), predup34(tex.size()), preddown34(tex.size());

  for(unsigned int i(0); i<Aval.size(); ++i){
    GetSummedVals(Aup.at(i), Adown.at(i), Aval.at(i), A,
		  weights, lower.at(i), upper.at(i));
    GetSummedVals(Bup.at(i), Bdown.at(i), Bval.at(i), B,
		  weights, lower.at(i), upper.at(i));
    GetSummedVals(C3up.at(i), C3down.at(i), C3val.at(i), C3,
		  weights, lower.at(i), upper.at(i));
    GetSummedVals(D3up.at(i), D3down.at(i), D3val.at(i), D3,
		  weights, lower.at(i), upper.at(i));
    GetSummedVals(C2up.at(i), C2down.at(i), C2val.at(i), C2,
		  weights, lower.at(i), upper.at(i));
    GetSummedVals(D2up.at(i), D2down.at(i), D2val.at(i), D2,
		  weights, lower.at(i), upper.at(i));
    GetSummedKappas(kappa3up.at(i), kappa3down.at(i), kappa3val.at(i), A, B, C3, D3,
		    weights, lower.at(i), upper.at(i));
    GetSummedKappas(kappa2up.at(i), kappa2down.at(i), kappa2val.at(i), A, B, C2, D2,
		    weights, lower.at(i), upper.at(i));
    GetSummedAPred(predup23.at(i), preddown23.at(i), pred23.at(i), D3, C2, D2, weights,
		   lower.at(i), upper.at(i));
    GetSummedAPred(predup24.at(i), preddown24.at(i), pred24.at(i), B, C2, D2, weights,
		   lower.at(i), upper.at(i));
    GetSummedAPred(predup34.at(i), preddown34.at(i), pred34.at(i), B, C3, D3, weights,
		   lower.at(i), upper.at(i));
    PrintLine( tex.at(i), Aval.at(i), Aup.at(i), Adown.at(i), Bval.at(i), Bup.at(i),
	       Bdown.at(i), C3val.at(i), C3up.at(i), C3down.at(i), D3val.at(i),
	       D3up.at(i), D3down.at(i), C2val.at(i), C2up.at(i), C2down.at(i),
	       D2val.at(i), D2up.at(i), D2down.at(i), kappa3val.at(i), kappa3up.at(i),
	       kappa3down.at(i), kappa2val.at(i), kappa2up.at(i), kappa2down.at(i),
	       pred23.at(i), predup23.at(i), preddown23.at(i), pred24.at(i),
	       predup24.at(i), preddown24.at(i), pred34.at(i), predup34.at(i),
	       preddown34.at(i));
  }

  KillTrees(trees);
  KillFiles(files);
}
