#include "qcd_systematic.hpp"

#include <string>
#include "TChain.h"
#include "style.hpp"
#include "timer.hpp"
#include "utils.hpp"

int main(){
  SetStyle();
  const std::string cut_baseline("passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && !passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut");
  const std::string cut_2b("num_b_tagged_jets==2");
  const std::string cut_3b("num_b_tagged_jets==3");
  const std::string cut_4b("num_b_tagged_jets>=4");
  const std::string cut_sig("higgs_mass_signal_region");
  const std::string cut_sb("higgs_mass_sideband");
  const std::string cut_sbin1("met_sig>=30.0 && met_sig<50.0");
  const std::string cut_sbin2("met_sig>=50.0 && met_sig<100.0");
  const std::string cut_sbin3("met_sig>=100.0 && met_sig<150.0");
  const std::string cut_sbin4("met_sig>=150.0");
  
  TChain chain("chain","chain");

  /*chain.Add("reduced_trees/BJets*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/T*channel*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTH*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTJets_FullLept*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTJets_SemiLept*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTJets_Hadronic*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTW*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTZ*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/W*JetsToLNu_TuneZ2Star*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/WH_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/WW_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/WZ_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/ZH_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/ZZ_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/ZJetsToNuNu*SyncSkim.root/reduced_tree");*/

  chain.Add("reduced_trees/MET_Run2012*SyncSkim.root/reduced_tree");

  double counts[3][2][4], uncerts[3][2][4];

  Timer timer(24);
  timer.Start();
  for(unsigned short nb_bin(0); nb_bin<3; ++nb_bin){
    std::string cut_nb("");
    switch(nb_bin){
    case 0:
      cut_nb=cut_2b;
      break;
    case 1:
      cut_nb=cut_3b;
      break;
    case 2:
      cut_nb=cut_4b;
      break;
    default:
      break;
    }
    for(unsigned short mbb_bin(0); mbb_bin<2; ++mbb_bin){
      std::string cut_mbb("");
      switch(mbb_bin){
      case 0:
	cut_mbb=cut_sig;
	break;
      case 1:
	cut_mbb=cut_sb;
	break;
      default:
      break;
      }
      for(unsigned short smet_bin(0); smet_bin<4; ++smet_bin){
	std::string cut_smet("");
	switch(smet_bin){
	case 0:
	  cut_smet=cut_sbin1;
	  break;
	case 1:
	  cut_smet=cut_sbin2;
	  break;
	case 2:
	  cut_smet=cut_sbin3;
	  break;
	case 3:
	  cut_smet=cut_sbin4;
	  break;
	default:
	  break;
	}

	const std::string cut_full("full_weight*("+cut_baseline+"&&"+cut_nb+"&&"+cut_mbb+"&&"+cut_smet+")");

	get_count_and_uncertainty(chain, cut_full, counts[nb_bin][mbb_bin][smet_bin], uncerts[nb_bin][mbb_bin][smet_bin]);
	timer.Iterate();
	timer.PrintRemainingTime();
      }
    }
  }
}
