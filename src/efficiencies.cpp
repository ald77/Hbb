#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include "TChain.h"
#include "timer.hpp"

int main(){
  std::vector<std::vector<std::string> > names(2);
  //names.at(0).push_back("reduced_trees/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71.root/reduced_tree");
  //names.at(0).push_back("reduced_trees/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext1-v1_AODSIM_UCSB1962_v71.root/reduced_tree");
  //names.at(0).push_back("reduced_trees/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext2-v1_AODSIM_UCSB1959_v71.root/reduced_tree");
  names.at(0).push_back("reduced_trees/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1850_v71.root/reduced_tree");
  //names.at(1).push_back("reduced_trees/SMS-TChiHH_2b2b_2J_mChargino-350_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1871_v71_SyncSkim.root/reduced_tree");
  names.at(1).push_back("reduced_trees/SMS-TChiHH_2b2b_2J_mChargino-130to500_mLSP-1_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19-v1_AODSIM_UCSB1992_v71.root/reduced_tree");
  std::vector<TChain*> chains(names.size(),NULL);
  for(unsigned int chain(0); chain<chains.size(); ++chain){
    chains.at(chain)=new TChain(names.at(chain).at(0).c_str(), names.at(chain).at(0).c_str());
    for(unsigned int name(0); name<names.at(chain).size(); ++name){
      chains.at(chain)->Add(names.at(chain).at(name).c_str());
      std::cout << "Adding to chain " << chain << ": " << names.at(chain).at(name) << std::endl;
    }
  }

  bool passesJSONCut(false), passesPVCut(false), passesJet2PtCut(false),
    passes2CSVTCut(false), passesMETSig30Cut(false), passesMETCleaningCut(false),
    passesTriggerCut(false), passesNumJetsCut(false), passesMinDeltaPhiCut(false),
    passesLeptonVetoCut(false), passesIsoTrackVetoCut(false), passesDRCut(false),
    passesHiggsAvgMassCut(false), passesHiggsMassDiffCut(false);

  bool passesBTaggingCut(false);//Crap, I forgot to put this in the reduced trees...

  float lumi_weight(0.0);

  float lsp_mass(0.0), chargino_mass(0.0);

  unsigned short num_b_tagged_jets(0);

  std::vector<bool*> cuts(0);
  cuts.push_back(&passesJSONCut);
  cuts.push_back(&passes2CSVTCut);
  cuts.push_back(&passesMETSig30Cut);
  cuts.push_back(&passesTriggerCut);
  cuts.push_back(&passesPVCut);
  cuts.push_back(&passesMETCleaningCut);
  cuts.push_back(&passesNumJetsCut);
  cuts.push_back(&passesJet2PtCut);
  cuts.push_back(&passesLeptonVetoCut);
  cuts.push_back(&passesIsoTrackVetoCut);
  cuts.push_back(&passesMinDeltaPhiCut);
  cuts.push_back(&passesBTaggingCut);
  cuts.push_back(&passesHiggsAvgMassCut);
  cuts.push_back(&passesHiggsMassDiffCut);
  cuts.push_back(&passesDRCut);

  std::vector<std::vector<unsigned int> > passes_initial(chains.size(),std::vector<unsigned int>(cuts.size(),0u));
  std::vector<std::vector<unsigned int> > passes_in_place(chains.size(),std::vector<unsigned int>(cuts.size(),0u));
  std::vector<std::vector<unsigned int> > passes_nm1(chains.size(),std::vector<unsigned int>(cuts.size(),0u));
  std::vector<std::vector<unsigned int> > totals_initial(chains.size(),std::vector<unsigned int>(cuts.size(),0u));
  std::vector<std::vector<unsigned int> > totals_in_place(chains.size(),std::vector<unsigned int>(cuts.size(),0u));
  std::vector<std::vector<unsigned int> > totals_nm1(chains.size(),std::vector<unsigned int>(cuts.size(),0u));
  std::vector<double> weights(chains.size(),1.0);

  for(unsigned int chain(0); chain<chains.size(); ++chain){
    chains.at(chain)->SetBranchStatus("*",0);
    chains.at(chain)->SetBranchStatus("passesJSONCut", 1);
    chains.at(chain)->SetBranchStatus("passesPVCut", 1);
    chains.at(chain)->SetBranchStatus("passesJet2PtCut", 1);
    chains.at(chain)->SetBranchStatus("passes2CSVTCut", 1);
    chains.at(chain)->SetBranchStatus("passesMETSig30Cut", 1);
    chains.at(chain)->SetBranchStatus("passesMETCleaningCut", 1);
    chains.at(chain)->SetBranchStatus("passesTriggerCut", 1);
    chains.at(chain)->SetBranchStatus("passesNumJetsCut", 1);
    chains.at(chain)->SetBranchStatus("passesMinDeltaPhiCut", 1);
    chains.at(chain)->SetBranchStatus("passesLeptonVetoCut", 1);
    chains.at(chain)->SetBranchStatus("passesIsoTrackVetoCut", 1);
    chains.at(chain)->SetBranchStatus("passesDRCut", 1);
    chains.at(chain)->SetBranchStatus("passesHiggsAvgMassCut", 1);
    chains.at(chain)->SetBranchStatus("passesHiggsMassDiffCut", 1);

    chains.at(chain)->SetBranchStatus("lumi_weight", 1);

    chains.at(chain)->SetBranchStatus("lsp_mass", 1);
    chains.at(chain)->SetBranchStatus("chargino_mass", 1);

    chains.at(chain)->SetBranchStatus("num_b_tagged_jets", 1);
    const int num_entries(chains.at(chain)->GetEntries());
    std::cout << "Num. entries: " << num_entries << std::endl;
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
    chains.at(chain)->SetBranchAddress("passesHiggsAvgMassCut", &passesHiggsAvgMassCut);
    chains.at(chain)->SetBranchAddress("passesHiggsMassDiffCut", &passesHiggsMassDiffCut);

    chains.at(chain)->SetBranchAddress("lumi_weight", &lumi_weight);

    chains.at(chain)->SetBranchAddress("lsp_mass", &lsp_mass);
    chains.at(chain)->SetBranchAddress("chargino_mass", &chargino_mass);

    chains.at(chain)->SetBranchAddress("num_b_tagged_jets", &num_b_tagged_jets);

    Timer timer(num_entries);
    timer.Start();
    for(int entry(0); entry<num_entries; ++entry){
      if(entry%(1<<18)==0 && entry!=0){
        timer.PrintRemainingTime();
      }
      chains.at(chain)->GetEntry(entry);
      weights.at(chain)=lumi_weight;
      timer.Iterate();
      if(chargino_mass!=-1 && chargino_mass!=250) continue;
      if(lsp_mass!=-1 && lsp_mass!=1) continue;
      for(unsigned int cut(0); cut<cuts.size(); ++cut){
        ++totals_initial.at(chain).at(cut);
        bool passes_before(true), passes_after(true);
        for(unsigned int cut_before(0); cut_before<cut && passes_before==true; ++cut_before){
          if(!(*(cuts.at(cut_before)))){
            passes_before=false;
          }
        }
        for(unsigned int cut_after(cut+1); cut_after<cuts.size() && passes_after==true; ++cut_after){
          if(!(*(cuts.at(cut_after)))){
            passes_after=false;
          }
        }
        
        if(passes_before){
          ++totals_in_place.at(chain).at(cut);
          if(passes_after){
            ++totals_nm1.at(chain).at(cut);
          }
        }
        if(*(cuts.at(cut))){
          ++passes_initial.at(chain).at(cut);
          if(passes_before){
            ++passes_in_place.at(chain).at(cut);
            if(passes_after){
              ++passes_nm1.at(chain).at(cut);
            }
          }
        }
      }
    }

    for(unsigned int cut(0); cut<cuts.size(); ++cut){
      std::cout << "chain=" << chain << ", cut=" << cut << "\n"
                << passes_initial.at(chain).at(cut) << "/" << totals_initial.at(chain).at(cut) << "\n"
                << passes_in_place.at(chain).at(cut) << "/" << totals_in_place.at(chain).at(cut) << "\n"
                << passes_nm1.at(chain).at(cut) << "/" << totals_nm1.at(chain).at(cut) << "\n\n";
    }
  }
  std::vector<std::vector<unsigned int> > sample_lists(2);

  for(unsigned int list(0); list<sample_lists.size(); ++list){
    std::vector<double> total_events_initial(cuts.size(),0.0);
    std::vector<double> total_events_in_place(cuts.size(),0.0);
    std::vector<double> total_events_nm1(cuts.size(),0.0);
    std::vector<double> eff_initial(cuts.size(),0.0);
    std::vector<double> eff_in_place(cuts.size(),0.0);
    std::vector<double> eff_nm1(cuts.size(),0.0);
    std::vector<double> uncert_initial(cuts.size(),0.0);
    std::vector<double> uncert_in_place(cuts.size(),0.0);
    std::vector<double> uncert_nm1(cuts.size(),0.0);
    for(unsigned int sample(0); sample<sample_lists.at(list).size(); ++sample){
      for(unsigned int cut(0); cut<cuts.size(); ++cut){
        const unsigned int index(sample_lists.at(list).at(sample));
        const double weight_squared(weights.at(index)*weights.at(index));
        total_events_initial.at(cut)+=totals_initial.at(index).at(cut)*weights.at(index);
        total_events_in_place.at(cut)+=totals_in_place.at(index).at(cut)*weights.at(index);
        total_events_nm1.at(cut)+=totals_nm1.at(index).at(cut)*weights.at(index);
        const double ratio_initial(passes_initial.at(index).at(cut)/(static_cast<double>(totals_initial.at(index).at(cut))));
        const double ratio_in_place(passes_in_place.at(index).at(cut)/(static_cast<double>(totals_in_place.at(index).at(cut))));
        const double ratio_nm1(passes_nm1.at(index).at(cut)/(static_cast<double>(totals_nm1.at(index).at(cut))));
        eff_initial.at(cut)+=passes_initial.at(index).at(cut)*weights.at(index);
        eff_in_place.at(cut)+=passes_in_place.at(index).at(cut)*weights.at(index);
        eff_nm1.at(cut)+=passes_nm1.at(index).at(cut)*weights.at(index);
        uncert_initial.at(cut)+=ratio_initial*(1.0-ratio_initial)*totals_initial.at(index).at(cut)*weight_squared;
        uncert_in_place.at(cut)+=ratio_in_place*(1.0-ratio_in_place)*totals_in_place.at(index).at(cut)*weight_squared;
        uncert_nm1.at(cut)+=ratio_nm1*(1.0-ratio_nm1)*totals_nm1.at(index).at(cut)*weight_squared;
      }
    }

    for(unsigned int cut(0); cut<cuts.size(); ++cut){
      eff_initial.at(cut)/=total_events_initial.at(cut);
      eff_in_place.at(cut)/=total_events_in_place.at(cut);
      eff_nm1.at(cut)/=total_events_nm1.at(cut);
      uncert_initial.at(cut)/=(total_events_initial.at(cut)*total_events_initial.at(cut));
      uncert_in_place.at(cut)/=(total_events_in_place.at(cut)*total_events_in_place.at(cut));
      uncert_nm1.at(cut)/=(total_events_nm1.at(cut)*total_events_nm1.at(cut));
      uncert_initial.at(cut)=sqrt(uncert_initial.at(cut));
      uncert_in_place.at(cut)=sqrt(uncert_in_place.at(cut));
      uncert_nm1.at(cut)=sqrt(uncert_nm1.at(cut));

      std::cout << "list=" << list << ", cut=" << cut << std::endl;
      std::cout << eff_initial.at(cut) << "+-" << uncert_initial.at(cut) << std::endl;
      std::cout << eff_in_place.at(cut) << "+-" << uncert_in_place.at(cut) << std::endl;
      std::cout << eff_nm1.at(cut) << "+-" << uncert_nm1.at(cut) << std::endl;
    }
  }
}
