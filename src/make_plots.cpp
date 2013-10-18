#include <iostream>
#include <string>
#include <unistd.h>
#include "TH1.h"
#include "event_handler.hpp"
#include "weights.hpp"

int main(int argc, char *argv[]){
  TH1::SetDefaultSumw2(true);
  std::string inFilename("");
  bool iscfA(false);
  int c(0);
  while((c=getopt(argc, argv, "i:c"))!=-1){
    switch(c){
    case 'i':
      inFilename=optarg;
      break;
    case 'c':
      iscfA=true;
      break;
    }
  }
  
  std::string outFilename("");
  if(iscfA){
    outFilename="raw_plots_and_values/"+inFilename+".root";
    inFilename="/net/cms2/cms2r0/cfA/"+inFilename+"/cfA_"+inFilename+"*.root";
  }else{
    std::string baseName(inFilename);
    size_t pos(baseName.find(".root"));
    if(pos!=std::string::npos){
      baseName.erase(pos);
    }
    pos=baseName.rfind("/");
    if(pos!=std::string::npos){
      if(pos!=baseName.size()-1){
	baseName.erase(0,pos+1);
      }else{
	baseName.append("file_name_ended_with_slash");
      }
    }
    outFilename="raw_plots_and_values/"+baseName+".root";
    std::cout << inFilename << "\n" << baseName << "\n" << outFilename << "\n";
  }

  WeightCalculator w(19399);
  EventHandler eH(inFilename, false, w.GetWeight(inFilename), true); // turning on fastMode!
  //EventHandler eH(inFilename, true, w.GetWeight("/cms2r0/cfA/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71/cfA_TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1884_v71_f1_1_cRU.root"));
  eH.MakePlots(outFilename);
}
