#include "cut.hpp"
#include <map>
#include <string>

using namespace std;

map<TString, TString> Cut::varToBoolTable;
//map<TString, double> Cut::varToValTable;

Cut::Cut(TString varIn, double minValIn, double maxValIn) :
  var(varIn),
  minVal(minValIn),
  maxVal(maxValIn)
{
  SetCuts();
}

TString Cut::GetVar() const {
  return var;
}

double Cut::GetMinVal() const {
  return minVal;
}

double Cut::GetMaxVal() const {
  return maxVal;
}

TString Cut::Get_Var_Cut() const{ // returns the cut on this variable
  for(map<TString, TString>::iterator it(varToBoolTable.begin());
      it!=varToBoolTable.end(); ++it){
    if(var.EqualTo(it->first)){
      return it->second.Data();
    }
  }
  return "";
}

TCut Cut::Get_N_Minus_1_Cuts() const{ // returns all other cuts

  TCut baseline("passesJSONCut&&passesPVCut&&passesJet2PtCut&&passes2CSVTCut&&passesMETCleaningCut&&passesTriggerCut");

  //TCut three_tag("third_highest_csv>0.679");
  //TCut four_tag("fourth_highest_csv>0.244");

  // TCut full(baseline+three_tag+four_tag+"passesMETSig30Cut&&passesNumJetsCut&&passesMinDeltaPhiCut&&passesLeptonVetoCut&&passesIsoTrackVetoCut&&passesDRCut&&passesHiggsAvgMassCut&passesHiggsMassDiffCut"); 

  TCut cut_n_minus_1(baseline);

  // Build the TCut, skipping the cut on the variable to plot
  for(map<TString, TString>::iterator it(varToBoolTable.begin());
      it!=varToBoolTable.end(); ++it){
    if(var.EqualTo(it->first)){
      continue;
    }
    if(var.Contains("csv")&&it->first.Contains("csv")) { // for b-tagging, take off 2 cuts
      continue;      
    }
    cut_n_minus_1+=it->second;
  }
  //cut_n_minus_1*=full_weight; // to apply weights here
  return cut_n_minus_1;

}


void Cut::SetCuts(){
  varToBoolTable["met_sig"]="passesMETSig30Cut";
  varToBoolTable["met"]="passesMETSig30Cut";
  varToBoolTable["num_jets"]="passesNumJetsCut";
  varToBoolTable["min_delta_phi"]="passesMinDeltaPhiCut";
  varToBoolTable["num_electrons"]="passesLeptonVetoCut";
  varToBoolTable["num_muons"]="passesLeptonVetoCut";
  varToBoolTable["num_taus"]="passesLeptonVetoCut";
  varToBoolTable["num_iso_tracks"]="passesIsoTrackVetoCut";
  varToBoolTable["max_delta_R"]="passesDRCut";
  varToBoolTable["third_highest_csv"]="third_highest_csv>0.679";
  varToBoolTable["fourth_highest_csv"]="fourth_highest_csv>0.244";
  varToBoolTable["average_higgs_mass"]="passesHiggsAvgMassCut";
  varToBoolTable["higgs_mass_difference"]="passesHiggsMassDiffCut";
}

