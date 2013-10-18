#include "event_handler.hpp"
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <utility>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <climits>
#include <cfloat>
#include <ctime>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <functional>
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "pu_constants.hpp"
#include "lumi_reweighting_stand_alone.hpp"
#include "cfa.hpp"
#include "event_number.hpp"
#include "b_jet.hpp"
#include "timer.hpp"
#include "math.hpp"
#include "in_json_2012.hpp"

const double EventHandler::CSVTCut(0.898);
const double EventHandler::CSVMCut(0.679);
const double EventHandler::CSVLCut(0.244);
const std::vector<std::vector<int> > VRunLumiPrompt(MakeVRunLumi("Golden"));
const std::vector<std::vector<int> > VRunLumi24Aug(MakeVRunLumi("24Aug"));
const std::vector<std::vector<int> > VRunLumi13Jul(MakeVRunLumi("13Jul"));

EventHandler::EventHandler(const std::string &fileName, const bool isList, const double scaleFactorIn):
  cfA(fileName, isList),
  higgsBJetPairing(std::make_pair(TLorentzVector(0.0,0.0,0.0,0.0),TLorentzVector(0.0,0.0,0.0,0.0)),std::make_pair(TLorentzVector(0.0,0.0,0.0,0.0),TLorentzVector(0.0,0.0,0.0,0.0))),
  sortedBJetCache(0),
  higgsPairingUpToDate(false),
  bJetsUpToDate(false),
  betaUpToDate(false),
  scaleFactor(scaleFactorIn),
  beta(0){
}

void EventHandler::SetScaleFactor(const double crossSection, const double luminosity, int numEntries){
  const int maxEntriesA(chainA.GetEntries()), maxEntriesB(chainB.GetEntries());
  if(maxEntriesA!=maxEntriesB){
    fprintf(stderr,"Error: Chains have different numbers of entries.\n");
  }
  if(maxEntriesA==0 || maxEntriesB==0){
    fprintf(stderr, "Error: Empty chains.\n");
  }
  if(numEntries<0){
    numEntries=maxEntriesA;
  }
  if(numEntries>0 && luminosity>0.0 && crossSection>0.0){
    scaleFactor=luminosity*crossSection/static_cast<double>(numEntries);
  }else{
    scaleFactor=1.0;
  }
}

void EventHandler::GetEntry(const unsigned int entry){
  cfA::GetEntry(entry);
  higgsPairingUpToDate=false;
  bJetsUpToDate=false;
  betaUpToDate=false;
}

int EventHandler::GetcfAVersion() const{
  size_t pos(sampleName.rfind("_v"));
  if(pos!=std::string::npos && pos<(sampleName.size()-2)){
    std::istringstream iss(sampleName.substr(pos+2));
    int version(0);
    iss >> version;
    return version;
  }else{
    return 0;
  }
}

void EventHandler::GetBeta(const std::string which) const{
  betaUpToDate=true;

  //Clear out the vector before starting a new event!
  beta.clear();

  if (GetcfAVersion()<69){
    beta.resize(jets_AK5PF_pt->size(), 0.0);
  }else{
    int totjet = 0;
    int matches = 0;
    for (unsigned int ijet=0; ijet<jets_AK5PF_pt->size(); ++ijet) {
      const float pt = jets_AK5PF_pt->at(ijet);
      const float eta = fabs(jets_AK5PF_eta->at(ijet));
      
      int i = 0;
      totjet++;
      for (std::vector<std::vector<float> >::const_iterator itr = puJet_rejectionBeta->begin(); itr != puJet_rejectionBeta->end(); ++itr, ++i) {
        int j = 0;
        float mypt = 0;
        float myeta = 0;
        float mybeta = 0;
        float result = 0;
        float tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0;
        for ( std::vector<float>::const_iterator it = itr->begin(); it != itr->end(); ++it, ++j) {
	  
          if ( (j%6)==0 ) mypt = *it;  
          if ( (j%6)==1 ) myeta = fabs(*it); 
          if ( (j%6)==2 ) tmp1 = *it;  
          if ( (j%6)==3 ) tmp2 = *it;  
          if ( (j%6)==4 ) tmp3 = *it;  
          if ( (j%6)==5 ) tmp4 = *it;  
	  
          //if ( which == "beta" )                 result = tmp1; 
          //else if ( which == "betaStar" )        result = tmp2; 
          //else if ( which == "betaClassic" )     result = tmp3; 
          //else if ( which == "betaStarClassic" ) result = tmp4; 
          //else result = -5; //Don't assert..
	  
          if ( which.compare("beta")==0 )                 result = tmp1; 
          else if ( which.compare("betaStar")==0 )        result = tmp2; 
	  else if ( which.compare("betaClassic")==0 )     result = tmp3; 
	  else if ( which.compare("betaStarClassic")==0 ) result = tmp4; 
          else result = -5; //Don't assert..
	  
        }//vector of info of each jet
        if ( mypt == pt && myeta == eta ) {
          matches++;
          mybeta = result;
          beta.push_back(mybeta);
          break;
        }     
      }//vector of jets
    } //ijet
  }
}

void EventHandler::SetScaleFactor(const double scaleFactorIn){
  scaleFactor=scaleFactorIn;
}

bool EventHandler::PassesPVCut() const{
  if(beamSpot_x->size()<1 || pv_x->size()<1) return false;
  const double pv_rho(sqrt(pv_x->at(0)*pv_x->at(0) + pv_y->at(0)*pv_y->at(0)));
  if(pv_ndof->at(0)>4 && fabs(pv_z->at(0))<24. && pv_rho<2.0 && pv_isFake->at(0)==0) return true;
  return false;
}

bool EventHandler::PassesMETCleaningCut() const{
  for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
    if(isProblemJet(jet)) return false;
  }
  if(pfTypeImets_et->at(0)>2.0*pfmets_et->at(0)) return false;
  return cschalofilter_decision
    && (hbhefilter_decision || sampleName.find("TChihh")!=std::string::npos || sampleName.find("HbbHbb")!=std::string::npos)
    && hcallaserfilter_decision 
    && ecalTPfilter_decision 
    && trackingfailurefilter_decision 
    && eebadscfilter_decision 
    && (ecallaserfilter_decision || sampleName.find("_v66")!=std::string::npos)
    && greedymuonfilter_decision 
    && inconsistentPFmuonfilter_decision 
    && scrapingVeto_decision
    && PassesBadJetFilter()
    && GetPBNR()>=1;
}

bool EventHandler::PassesTriggerCut() const{
  for(unsigned int a=0; a<trigger_name->size(); ++a){
    if((trigger_name->at(a).find("HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v")!=std::string::npos
	|| trigger_name->at(a).find("HLT_PFMET150_v")!=std::string::npos
	|| trigger_name->at(a).find("HLT_DiCentralPFJet30_PFMHT80_v")!=std::string::npos)
       && trigger_prescalevalue->at(a)==1 && trigger_decision->at(a)==1){
      return true;
    }
  }
  return false;
}

bool EventHandler::PassesNumJetsCut() const{
  const int numGoodJets(GetNumGoodJets());
  if(numGoodJets==4 || numGoodJets==5){
    return true;
  }else{
    return false;
  }
}

bool EventHandler::PassesJet2PtCut() const{
  std::vector<double> pts(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i)){
      pts.push_back(jets_AK5PF_pt->at(i));
    }
  }
  std::sort(pts.begin(), pts.end(), std::greater<double>());
  if(pts.size()>1 && pts.at(1)>50.0){
    return true;
  }else{
    return false;
  }
}

bool EventHandler::PassesHiggsMassCut() const{
  return PassesHiggsAvgMassCut() && PassesHiggsMassDiffCut();
}

bool EventHandler::PassesHiggsAvgMassCut() const{
  const std::pair<double,double> higgsMasses(GetHiggsMasses());
  const double avg(0.5*(higgsMasses.first+higgsMasses.second));
  return avg>100.0 && avg<140.0;
}

bool EventHandler::PassesHiggsMassDiffCut() const{
  const std::pair<double,double> higgsMasses(GetHiggsMasses());
  return fabs(higgsMasses.first-higgsMasses.second)<20.0;
}

bool EventHandler::PassesDRCut() const{
  return GetMaxDR()<2.2;
}

bool EventHandler::PassesInvertedDRCut() const{
  return GetMaxDR()>2.2;
}

bool EventHandler::PassesMETSig50Cut() const{
  return pfmets_fullSignif>50.0;
}

bool EventHandler::PassesMETSig80Cut() const{
  return pfmets_fullSignif>80.0;
}

bool EventHandler::PassesMETSig100Cut() const{
  return pfmets_fullSignif>100.0;
}

bool EventHandler::PassesMETSig150Cut() const{
  return pfmets_fullSignif>150.0;
}

bool EventHandler::PassesJSONCut() const{
  if(sampleName.find("Run2012")!=std::string::npos){
    if(sampleName.find("PromptReco")!=std::string::npos
       &&!inJSON(VRunLumiPrompt, run, lumiblock)) return false;
    if(sampleName.find("24Aug")!=std::string::npos
       && !inJSON(VRunLumi24Aug, run, lumiblock)) return false;
    if(sampleName.find("13Jul")!=std::string::npos
       && !inJSON(VRunLumi13Jul, run, lumiblock)) return false;
    return true;
  }else{
    return true;
  }
}

uint_least32_t EventHandler::GetCutFailCode() const{
  if(!betaUpToDate) GetBeta();
  if(!bJetsUpToDate) GetSortedBJets();
  if(!higgsPairingUpToDate) GetHiggsBJetPairing();
  uint_least32_t fail_code(0);
  if(!PassesMETSig150Cut()) fail_code |= 0x00001;
  if(!PassesMETSig100Cut()) fail_code |= 0x00002;
  if(!PassesMETSig50Cut()) fail_code |= 0x00004;
  if(!PassesMETSig30Cut()) fail_code |= 0x00008;
  if(!PassesDRCut()) fail_code |= 0x00010;
  if(!PassesHiggsAvgMassCut()) fail_code |= 0x00020;
  if(!PassesHiggsMassDiffCut()) fail_code |= 0x00040;
  if(!PassesBTaggingCut()) fail_code |= 0x00080;
  if(!(GetNumCSVTJets()>=2 && GetNumCSVMJets()>=3)) fail_code |= 0x00100;
  if(!(GetNumCSVTJets()>=2)) fail_code |= 0x00200;
  if(!PassesIsoTrackVetoCut()) fail_code |= 0x00400;
  if(!PassesLeptonVetoCut()) fail_code |= 0x00800;
  if(!PassesMinDeltaPhiCut()) fail_code |= 0x01000;
  if(!PassesNumJetsCut()) fail_code |= 0x02000;
  if(!PassesTriggerCut()) fail_code |= 0x04000;
  if(!PassesJSONCut()) fail_code |= 0x08000;
  if(!PassesMETCleaningCut()) fail_code |= 0x10000;
  if(!PassesPVCut()) fail_code |= 0x20000;
  if(!PassesTChiZHMassCut()) fail_code |= 0x40000;
  if(!PassesTChiHHMassCut()) fail_code |= 0x80000;
  return fail_code;
}

bool EventHandler::PassesRegionACut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || !PassesBTaggingCut()) return false;
  return true;
}

bool EventHandler::PassesRegionBCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || !PassesBTaggingCut()) return false;
  return true;
}

bool EventHandler::PassesRegionC3bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || GetNumCSVMJets()<3 || GetNumCSVLJets()>=4) return false;
  return true;
}

bool EventHandler::PassesRegionD3bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || GetNumCSVMJets()<3 || GetNumCSVLJets()>=4) return false;
  return true;
}

bool EventHandler::PassesRegionC2bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || GetNumCSVMJets()>=3) return false;
  return true;
}

bool EventHandler::PassesRegionD2bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || GetNumCSVMJets()>=3) return false;
  return true;
}

bool EventHandler::PassesInvertedDRRegionACut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesInvertedDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || !PassesBTaggingCut()) return false;
  return true;
}

bool EventHandler::PassesInvertedDRRegionBCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesInvertedDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || !PassesBTaggingCut()) return false;
  return true;
}

bool EventHandler::PassesInvertedDRRegionC3bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesInvertedDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || GetNumCSVMJets()<3 || GetNumCSVLJets()>=4) return false;
  return true;
}

bool EventHandler::PassesInvertedDRRegionD3bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesInvertedDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || GetNumCSVMJets()<3 || GetNumCSVLJets()>=4) return false;
  return true;
}

bool EventHandler::PassesInvertedDRRegionC2bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesInvertedDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || GetNumCSVMJets()>=3) return false;
  return true;
}

bool EventHandler::PassesInvertedDRRegionD2bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesLeptonVetoCut()) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  if(!PassesInvertedDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || GetNumCSVMJets()>=3) return false;
  return true;
}

bool EventHandler::PassesSingleLeptonRegionACut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesSingleLeptonCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || !PassesBTaggingCut()) return false;
  return true;
}

bool EventHandler::PassesSingleLeptonRegionBCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesSingleLeptonCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || !PassesBTaggingCut()) return false;
  return true;
}

bool EventHandler::PassesSingleLeptonRegionC3bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesSingleLeptonCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || GetNumCSVMJets()<3 || GetNumCSVLJets()>=4) return false;
  return true;
}

bool EventHandler::PassesSingleLeptonRegionD3bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesSingleLeptonCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || GetNumCSVMJets()<3 || GetNumCSVLJets()>=4) return false;
  return true;
}

bool EventHandler::PassesSingleLeptonRegionC2bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesSingleLeptonCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  if(!PassesHiggsMassCut() || GetNumCSVMJets()>=3) return false;
  return true;
}

bool EventHandler::PassesSingleLeptonRegionD2bCut() const{
  if(!PassesPVCut()) return false;
  if(!PassesMETCleaningCut()) return false;
  if(!PassesTriggerCut()) return false;
  if(!PassesNumJetsCut()) return false;
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMinDeltaPhiCut()) return false;
  if(!PassesSingleLeptonCut()) return false;
  if(!PassesDRCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  const std::pair<double, double> higgs_masses(GetHiggsMasses());
  const double avg(0.5*(higgs_masses.first+higgs_masses.second));
  const double delta(fabs(higgs_masses.first-higgs_masses.second));
  if((delta<30.0 && avg>90.0 && avg<150.0) || GetNumCSVMJets()>=3) return false;
  return true;
}

bool EventHandler::PassesBadJetFilter() const{
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i,false,30.0,DBL_MAX) && !isGoodJet(i,true,30.0,DBL_MAX)) return false;
  }
  return true;
}

bool EventHandler::HasGluonSplitting() const{
  for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
    if(TMath::Nint(fabs(jets_AK5PF_parton_Id->at(jet)))==21
       && TMath::Nint(fabs(jets_AK5PF_partonFlavour->at(jet)))==5){
      return true;
    }
  }
  return false;
}

int EventHandler::GetPBNR() const{
  //RA2 particle-based noise rejection                                                                                                                                           

  bool nhBad=false;
  bool phBad=false;
  for (unsigned int it = 0; it<jets_AK5PF_pt->size(); it++) {
    //cfA version from Keith                                                                                                                                                     
    double NHF = jets_AK5PF_neutralHadE->at(it)/(jets_AK5PF_energy->at(it)*jets_AK5PF_corrFactorRaw->at(it));
    double PEF = jets_AK5PF_photonEnergy->at(it)/(jets_AK5PF_energy->at(it)*jets_AK5PF_corrFactorRaw->at(it));
    if (NHF > 0.9)  nhBad = true;
    if (PEF > 0.95) phBad = true;
  }

  if (nhBad && phBad) return -3;
  else if (phBad) return -2;
  else if (nhBad) return -1;
  return 1;
}

double EventHandler::GetHT(const bool useMET, const bool useLeps) const{
  double HT(0.0);
  if(useMET && pfmets_et->size()>0) HT+=pfmets_et->at(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i)) HT+=jets_AK5PF_pt->at(i);
  }
  if(useLeps){
    for(unsigned int i(0); i<pf_els_pt->size(); ++i){
      if(isVetoElectron(i)) HT+=pf_els_pt->at(i);
    }
    for(unsigned int i(0); i<pf_mus_pt->size(); ++i){
      if(isVetoMuon(i)) HT+=pf_mus_pt->at(i);
    }
    for(unsigned int i(0); i<taus_pt->size(); ++i){
      if(isVetoTau(i)) HT+=taus_pt->at(i);
    }
  }
  return HT;
}

unsigned int EventHandler::GetNumLowPtPfCands(const double ptThresh) const{
  unsigned int cands(0);
  for(unsigned int i(0); i<pfcand_pt->size(); ++i){
    if(pfcand_pt->at(i)<ptThresh) ++cands;
  }
  return cands;
}

double EventHandler::GetMaxDR() const{
  GetHiggsBJetPairing();
  if(sortedBJetCache.size()<4){
    return DBL_MAX;
  }else{
    const double dRa(higgsBJetPairing.first.first.DeltaR(higgsBJetPairing.first.second));
    const double dRb(higgsBJetPairing.second.first.DeltaR(higgsBJetPairing.second.second));
    return dRa>dRb?dRa:dRb;
  }
}

std::pair<double, double> EventHandler::GetHiggsMasses() const{
  GetHiggsBJetPairing();
  const double massA((higgsBJetPairing.first.first+higgsBJetPairing.first.second).M());
  const double massB((higgsBJetPairing.second.first+higgsBJetPairing.second.second).M());
  return (massA>massB)?(std::make_pair(massA,massB)):(std::make_pair(massB,massA));
}

double EventHandler::GetHiggsDeltaR() const{
  GetHiggsBJetPairing();
  return (higgsBJetPairing.first.first+higgsBJetPairing.first.second).DeltaR(higgsBJetPairing.second.first+higgsBJetPairing.second.second);
}

double EventHandler::GetMETOfLowPtPfCands(const double ptThresh) const{
  double px(0.0), py(0.0);
  for(unsigned int i(0); i<pfcand_pt->size(); ++i){
    if(pfcand_pt->at(i)<ptThresh){
      px+=pfcand_px->at(i);
      py+=pfcand_py->at(i);
    }
  }
  return sqrt(px*px+py*py);
}

void EventHandler::GetSortedBJets() const{
  if(!bJetsUpToDate){
    sortedBJetCache.clear();
    for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
      if(isGoodJet(i)){
	sortedBJetCache.push_back(BJet(TLorentzVector(jets_AK5PF_px->at(i),jets_AK5PF_py->at(i),jets_AK5PF_pz->at(i),jets_AK5PF_energy->at(i)),jets_AK5PF_btag_secVertexCombined->at(i))) \
	  ;
      }
    }
    std::sort(sortedBJetCache.begin(),sortedBJetCache.end(), std::greater<BJet>());
    bJetsUpToDate=true;
  }
}

void EventHandler::GetHiggsBJetPairing() const{
  if(!higgsPairingUpToDate){
    if(!bJetsUpToDate) GetSortedBJets();
    if(sortedBJetCache.size()<4){
      higgsBJetPairing=std::make_pair(std::make_pair(TLorentzVector(0.0,0.0,0.0,0.0),TLorentzVector(0.0,0.0,0.0,0.0)),std::make_pair(TLorentzVector(0.0,0.0,0.0,0.0),TLorentzVector(0.0,0.0,0.0,0.0)));
    }else{
      // Compute Higgs masses
      // Three pairings
      const double m1a((sortedBJetCache.at(0).GetLorentzVector()+sortedBJetCache.at(1).GetLorentzVector()).M());
      const double m1b((sortedBJetCache.at(2).GetLorentzVector()+sortedBJetCache.at(3).GetLorentzVector()).M());
      const double m2a((sortedBJetCache.at(0).GetLorentzVector()+sortedBJetCache.at(2).GetLorentzVector()).M());
      const double m2b((sortedBJetCache.at(1).GetLorentzVector()+sortedBJetCache.at(3).GetLorentzVector()).M());
      const double m3a((sortedBJetCache.at(0).GetLorentzVector()+sortedBJetCache.at(3).GetLorentzVector()).M());
      const double m3b((sortedBJetCache.at(1).GetLorentzVector()+sortedBJetCache.at(2).GetLorentzVector()).M());
      
      const double delta1(fabs(m1a-m1b)), delta2(fabs(m2a-m2b)), delta3(fabs(m3a-m3b));
      
      if(delta1<=delta2 && delta1<=delta3){
	higgsBJetPairing=std::make_pair(std::make_pair(sortedBJetCache.at(0).GetLorentzVector(),sortedBJetCache.at(1).GetLorentzVector()),std::make_pair(sortedBJetCache.at(2).GetLorentzVector(),sortedBJetCache.at(3).GetLorentzVector()));
      }else if(delta2<=delta1 && delta2<=delta3){
	higgsBJetPairing=std::make_pair(std::make_pair(sortedBJetCache.at(0).GetLorentzVector(),sortedBJetCache.at(2).GetLorentzVector()),std::make_pair(sortedBJetCache.at(1).GetLorentzVector(),sortedBJetCache.at(3).GetLorentzVector()));
      }else if(delta3<=delta1 && delta3<=delta2){
	higgsBJetPairing=std::make_pair(std::make_pair(sortedBJetCache.at(0).GetLorentzVector(),sortedBJetCache.at(3).GetLorentzVector()),std::make_pair(sortedBJetCache.at(1).GetLorentzVector(),sortedBJetCache.at(2).GetLorentzVector()));
      }
    }
  }
  higgsPairingUpToDate=true;
}

double EventHandler::GetHighestCSV(unsigned int pos) const{
  GetSortedBJets();
  --pos;
  if(pos>=sortedBJetCache.size()){
    return -DBL_MAX;
  }else{
    return sortedBJetCache.at(pos).GetBTag();
  }
}

bool EventHandler::PassesTChiZHMassCut() const{
  if(sampleName.find("TChiZH")==std::string::npos) return true;
  if(model_params->find("chargino300")==std::string::npos) return false;
  if(model_params->find("bino1_")==std::string::npos) return false;
  return true;
}

bool EventHandler::PassesTChiHHMassCut(int mChi, int mLSP) const{
  if (mChi<0||mLSP<0) return true;
  // parse the model_params string in the new signal samples for mass points
  char s_mChi[10], s_mLSP[10];
  sprintf(s_mChi, "chargino%d_", mChi);
  sprintf(s_mLSP, "bino%d_", mLSP);

  if(sampleName.find("SMS-TChiHH")==std::string::npos) return true;
  if(model_params->find(s_mChi)==std::string::npos) return false;
  if(model_params->find(s_mLSP)==std::string::npos) return false;
  return true;
}

int GetSimpleParticle(const double &id){
  const int iid(static_cast<int>(fabs(id)));
  switch(iid){
  case 1:
  case 2:
  case 3:
    return 1;
  case 4:
    return 4;
  case 5:
    return 5;
  case 6:
    return 6;
  case 11:
    return 11;
  case 13:
    return 13;
  case 15:
    return 15;
  case 12:
  case 14:
  case 16:
    return 12;
  case 21:
    return 21;
  case 22:
    return 22;
  case 23:
    return 23;
  case 24:
    return 24;
  case 25:
    return 25;
  default:
    return 0;
  }
}

std::vector<std::pair<int, int> > EventHandler::GetBOrigins() const{
  if(!bJetsUpToDate) GetSortedBJets();
  std::vector<std::pair<int, int> > x(0);
  for(unsigned int jet(0); jet<sortedBJetCache.size(); ++jet){
    double minDeltaR(DBL_MAX);
    unsigned int bestJet(-1);
    for(unsigned int mc(0); mc<mc_doc_pt->size(); ++mc){
      const double thisDeltaR(Math::GetDeltaR(sortedBJetCache.at(jet).GetLorentzVector().Phi(),
					      sortedBJetCache.at(jet).GetLorentzVector().Eta(),
					      mc_doc_phi->at(mc),
					      mc_doc_eta->at(mc)));
      if(thisDeltaR<minDeltaR){
	minDeltaR=thisDeltaR;
	bestJet=mc;
      }
    }
    if(bestJet<mc_doc_id->size()){
      const int id(GetSimpleParticle(mc_doc_id->at(bestJet)));
      const int mom(GetSimpleParticle(mc_doc_mother_id->at(bestJet)));
      x.push_back(std::pair<int,int>(id, mom));
    }else{
      x.push_back(std::pair<int, int>(0,0));
    }
  }
  return x;
}

void FixSbinLabels(TH1D &h){
  h.GetXaxis()->SetBinLabel(1,"1");
  h.GetXaxis()->SetBinLabel(2,"2");
  h.GetXaxis()->SetBinLabel(3,"3");
  h.GetXaxis()->SetBinLabel(4,"4");
}

void FixBinLabels(TH1D &h){
  h.GetXaxis()->SetBinLabel(1, "H->b");
  h.GetXaxis()->SetBinLabel(2, "t->b");
  h.GetXaxis()->SetBinLabel(3, "W->b");
  h.GetXaxis()->SetBinLabel(4, "g->b");
  h.GetXaxis()->SetBinLabel(5, "Z->b");
  h.GetXaxis()->SetBinLabel(6, "other->b");
  h.GetXaxis()->SetBinLabel(7, "W->c");
  h.GetXaxis()->SetBinLabel(8, "g->c");
  h.GetXaxis()->SetBinLabel(9, "Z->c");
  h.GetXaxis()->SetBinLabel(10, "other->c");
  h.GetXaxis()->SetBinLabel(11, "W->uds");
  h.GetXaxis()->SetBinLabel(12, "g->uds");
  h.GetXaxis()->SetBinLabel(13, "Z->uds");
  h.GetXaxis()->SetBinLabel(14, "other->uds");
  h.GetXaxis()->SetBinLabel(15, "gluon");
  h.GetXaxis()->SetBinLabel(16, "electron");
  h.GetXaxis()->SetBinLabel(17, "muon");
  h.GetXaxis()->SetBinLabel(18, "tau");
  h.GetXaxis()->SetBinLabel(19, "other");
}

int GetType(const std::pair<int,int> &bo){
  if(bo.first==5){
    if(bo.second==25){
      return 1;
    }else if(bo.second==6){
      return 2;
    }else if(bo.second==24){
      return 3;
    }else if(bo.second==21){
      return 4;
    }else if(bo.second==23){
      return 5;
    }else{
      return 6;
    }
  }else if(bo.first==4){
    if(bo.second==24){
      return 7;
    }else if(bo.second==21){
      return 8;
    }else if(bo.second==23){
      return 9;
    }else{
      return 10;
    }
  }else if(bo.first>=1 && bo.first<=3){
    if(bo.second==24){
      return 11;
    }else if(bo.second==21){
      return 12;
    }else if(bo.second==23){
      return 13;
    }else{
      return 14;
    }
  }else if(bo.first==21){
    return 15;
  }else if(bo.first==11){
    return 16;
  }else if(bo.first==13){
    return 17;
  }else if(bo.first==15){
    return 18;
  }else{
    return 19;
  }
}

void EventHandler::MakePlots(const std::string &outFileName){
  TH1::SetDefaultSumw2(true);
  std::set<EventNumber> eventList;
  unsigned int startCount(0), PVCount(0), METCleaningCount(0), TriggerCount(0), numJetsCount(0), CSVTCount(0), jet2PtCount(0), minDeltaPhiCount(0), leptonVetoCount(0), isoTrackVetoCount(0), bTagCount(0), higgsCount(0), DRCount(0), METSig30Count(0), METSig50Count(0), METSig100Count(0), METSig150Count(0);
  double startCountWeighted(0), PVCountWeighted(0), METCleaningCountWeighted(0), TriggerCountWeighted(0), numJetsCountWeighted(0), CSVTCountWeighted(0), jet2PtCountWeighted(0), minDeltaPhiCountWeighted(0), leptonVetoCountWeighted(0), isoTrackVetoCountWeighted(0), bTagCountWeighted(0), higgsCountWeighted(0), DRCountWeighted(0), METSig30CountWeighted(0), METSig50CountWeighted(0), METSig100CountWeighted(0), METSig150CountWeighted(0);

  unsigned int ACount(0), BCount(0), C3bCount(0), D3bCount(0), C2bCount(0), D2bCount(0);
  double ACountWeighted(0.0), BCountWeighted(0.0), C3bCountWeighted(0.0), D3bCountWeighted(0.0), C2bCountWeighted(0.0), D2bCountWeighted(0.0);
  double AUncert(0.0), BUncert(0.0), C3bUncert(0.0), D3bUncert(0.0), C2bUncert(0.0), D2bUncert(0.0);

  unsigned int ADRInvCount(0), BDRInvCount(0), C3bDRInvCount(0), D3bDRInvCount(0), C2bDRInvCount(0), D2bDRInvCount(0);
  double ADRInvCountWeighted(0.0), BDRInvCountWeighted(0.0), C3bDRInvCountWeighted(0.0), D3bDRInvCountWeighted(0.0), C2bDRInvCountWeighted(0.0), D2bDRInvCountWeighted(0.0);
  double ADRInvUncert(0.0), BDRInvUncert(0.0), C3bDRInvUncert(0.0), D3bDRInvUncert(0.0), C2bDRInvUncert(0.0), D2bDRInvUncert(0.0);

  unsigned int ASLCount(0), BSLCount(0), C3bSLCount(0), D3bSLCount(0), C2bSLCount(0), D2bSLCount(0);
  double ASLCountWeighted(0.0), BSLCountWeighted(0.0), C3bSLCountWeighted(0.0), D3bSLCountWeighted(0.0), C2bSLCountWeighted(0.0), D2bSLCountWeighted(0.0);
  double ASLUncert(0.0), BSLUncert(0.0), C3bSLUncert(0.0), D3bSLUncert(0.0), C2bSLUncert(0.0), D2bSLUncert(0.0);

  unsigned int ASLsbin1Count(0), BSLsbin1Count(0), C3bSLsbin1Count(0), D3bSLsbin1Count(0), C2bSLsbin1Count(0), D2bSLsbin1Count(0);
  double ASLsbin1CountWeighted(0.0), BSLsbin1CountWeighted(0.0), C3bSLsbin1CountWeighted(0.0), D3bSLsbin1CountWeighted(0.0), C2bSLsbin1CountWeighted(0.0), D2bSLsbin1CountWeighted(0.0);

  unsigned int ASLsbin2Count(0), BSLsbin2Count(0), C3bSLsbin2Count(0), D3bSLsbin2Count(0), C2bSLsbin2Count(0), D2bSLsbin2Count(0);
  double ASLsbin2CountWeighted(0.0), BSLsbin2CountWeighted(0.0), C3bSLsbin2CountWeighted(0.0), D3bSLsbin2CountWeighted(0.0), C2bSLsbin2CountWeighted(0.0), D2bSLsbin2CountWeighted(0.0);

  unsigned int ASLsbin3Count(0), BSLsbin3Count(0), C3bSLsbin3Count(0), D3bSLsbin3Count(0), C2bSLsbin3Count(0), D2bSLsbin3Count(0);
  double ASLsbin3CountWeighted(0.0), BSLsbin3CountWeighted(0.0), C3bSLsbin3CountWeighted(0.0), D3bSLsbin3CountWeighted(0.0), C2bSLsbin3CountWeighted(0.0), D2bSLsbin3CountWeighted(0.0);

  unsigned int ASLsbin4Count(0), BSLsbin4Count(0), C3bSLsbin4Count(0), D3bSLsbin4Count(0), C2bSLsbin4Count(0), D2bSLsbin4Count(0);
  double ASLsbin4CountWeighted(0.0), BSLsbin4CountWeighted(0.0), C3bSLsbin4CountWeighted(0.0), D3bSLsbin4CountWeighted(0.0), C2bSLsbin4CountWeighted(0.0), D2bSLsbin4CountWeighted(0.0);

  unsigned int ADRInvsbin1Count(0), BDRInvsbin1Count(0), C3bDRInvsbin1Count(0), D3bDRInvsbin1Count(0), C2bDRInvsbin1Count(0), D2bDRInvsbin1Count(0);
  double ADRInvsbin1CountWeighted(0.0), BDRInvsbin1CountWeighted(0.0), C3bDRInvsbin1CountWeighted(0.0), D3bDRInvsbin1CountWeighted(0.0), C2bDRInvsbin1CountWeighted(0.0), D2bDRInvsbin1CountWeighted(0.0);

  unsigned int ADRInvsbin2Count(0), BDRInvsbin2Count(0), C3bDRInvsbin2Count(0), D3bDRInvsbin2Count(0), C2bDRInvsbin2Count(0), D2bDRInvsbin2Count(0);
  double ADRInvsbin2CountWeighted(0.0), BDRInvsbin2CountWeighted(0.0), C3bDRInvsbin2CountWeighted(0.0), D3bDRInvsbin2CountWeighted(0.0), C2bDRInvsbin2CountWeighted(0.0), D2bDRInvsbin2CountWeighted(0.0);

  unsigned int ADRInvsbin3Count(0), BDRInvsbin3Count(0), C3bDRInvsbin3Count(0), D3bDRInvsbin3Count(0), C2bDRInvsbin3Count(0), D2bDRInvsbin3Count(0);
  double ADRInvsbin3CountWeighted(0.0), BDRInvsbin3CountWeighted(0.0), C3bDRInvsbin3CountWeighted(0.0), D3bDRInvsbin3CountWeighted(0.0), C2bDRInvsbin3CountWeighted(0.0), D2bDRInvsbin3CountWeighted(0.0);

  unsigned int ADRInvsbin4Count(0), BDRInvsbin4Count(0), C3bDRInvsbin4Count(0), D3bDRInvsbin4Count(0), C2bDRInvsbin4Count(0), D2bDRInvsbin4Count(0);
  double ADRInvsbin4CountWeighted(0.0), BDRInvsbin4CountWeighted(0.0), C3bDRInvsbin4CountWeighted(0.0), D3bDRInvsbin4CountWeighted(0.0), C2bDRInvsbin4CountWeighted(0.0), D2bDRInvsbin4CountWeighted(0.0);

  unsigned int Asbin1Count(0), Bsbin1Count(0), C3bsbin1Count(0), D3bsbin1Count(0), C2bsbin1Count(0), D2bsbin1Count(0);
  double Asbin1CountWeighted(0.0), Bsbin1CountWeighted(0.0), C3bsbin1CountWeighted(0.0), D3bsbin1CountWeighted(0.0), C2bsbin1CountWeighted(0.0), D2bsbin1CountWeighted(0.0);

  unsigned int Asbin2Count(0), Bsbin2Count(0), C3bsbin2Count(0), D3bsbin2Count(0), C2bsbin2Count(0), D2bsbin2Count(0);
  double Asbin2CountWeighted(0.0), Bsbin2CountWeighted(0.0), C3bsbin2CountWeighted(0.0), D3bsbin2CountWeighted(0.0), C2bsbin2CountWeighted(0.0), D2bsbin2CountWeighted(0.0);

  unsigned int Asbin3Count(0), Bsbin3Count(0), C3bsbin3Count(0), D3bsbin3Count(0), C2bsbin3Count(0), D2bsbin3Count(0);
  double Asbin3CountWeighted(0.0), Bsbin3CountWeighted(0.0), C3bsbin3CountWeighted(0.0), D3bsbin3CountWeighted(0.0), C2bsbin3CountWeighted(0.0), D2bsbin3CountWeighted(0.0);

  unsigned int Asbin4Count(0), Bsbin4Count(0), C3bsbin4Count(0), D3bsbin4Count(0), C2bsbin4Count(0), D2bsbin4Count(0);
  double Asbin4CountWeighted(0.0), Bsbin4CountWeighted(0.0), C3bsbin4CountWeighted(0.0), D3bsbin4CountWeighted(0.0), C2bsbin4CountWeighted(0.0), D2bsbin4CountWeighted(0.0);

  std::vector<float> dataDist(pu::RunsThrough203002, pu::RunsThrough203002+60);
  std::vector<float> MCDist(pu::Summer2012_S10, pu::Summer2012_S10+60);//QQQ this needs to change later for general pileup scenario
  reweight::LumiReWeighting lumiWeights(MCDist, dataDist);

  const double metSigABCDBinEdges[6]={30.0,50.0,100.0,150.0,200.0,300.0};
  TH1D xx_metSig_SigOverSB_2b("xx_metSig_SigOverSB_4b","Mass Signal Window/Mass SB;MET Significance;Mass Signal Window/Mass SB", 5, metSigABCDBinEdges);
  TH1D xx_metSig_SigOverSB_3b("xx_metSig_SigOverSB_3b","Mass Signal Window/Mass SB;MET Significance;Mass Signal Window/Mass SB", 5, metSigABCDBinEdges);
  TH1D xx_metSig_SigOverSB_4b("xx_metSig_SigOverSB_2b","Mass Signal Window/Mass SB;MET Significance;Mass Signal Window/Mass SB", 5, metSigABCDBinEdges);

  TFile file(outFileName.c_str(), "recreate");

  TCanvas xx_metSig_SigOverSB("xx_metSig_SigOverSB","xx_metSig_SigOverSB");

  TH1D xx_nm1_thirdBTag("xx_nm1_thirdBTag", "N-1 Third B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D xx_nm1_fourthBTag("xx_nm1_fourthBTag", "N-1 Fourth B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D xx_nm1_avgHiggsMass("xx_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D xx_nm1_higgsMassDiff("xx_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D xx_nm1_maxDeltaR("xx_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D xx_nm1_metSig("xx_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D xx_nm1_met("xx_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25 GeV", 40, 0.0, 1000.0);

  TH2D xx_MET_vs_minDeltaPhi("xx_MET_vs_minDeltaPhi","MET vs. Min. #Delta#phi;Min. #Delta#Phi;MET [GeV]",15,0.0,3.0,15,0.0,450.0);
  TH2D xx_METSig_vs_minDeltaPhi("xx_METSig_vs_minDeltaPhi","MET Significance vs. Min. #Delta#phi;Min. #Delta#Phi;MET Significance",15,0.0,3.0,15,0.0,300.0);

  TH2D xx_METOverSqrtHT_vs_METSig("xx_METOverSqrtHT_vs_METSig", "MET/#sqrt{HT} vs. MET Sig.;MET Significance; MET/#sqrt{HT} [#sqrt{GeV}]", 40, 0.0, 200.0, 24, 0.0, 12.0);
  TH1D xx_METOverSqrtHT_Over_METSig("xx_METOverSqrtHT_Over_METSig","(MET/#sqrt{HT})/MET Sig;(MET/#sqrt{HT})/MET Sig [#sqrt{GeV}];Events/19.4 fb^{-1}/0.01",30,0.0,0.3);
  TH1D xx_MET("xx_MET","MET;MET [GeV];Events/19.4 fb^{-1}/10 GeV",40,0.0,400.0);
  TH1D xx_ABCD("xx_ABCD", "ABCD;;Events/19.4 fb^{-1}",6,-0.5,5.5);
  xx_ABCD.GetXaxis()->SetBinLabel(1,"A");
  xx_ABCD.GetXaxis()->SetBinLabel(2,"B");
  xx_ABCD.GetXaxis()->SetBinLabel(3,"C(3b)");
  xx_ABCD.GetXaxis()->SetBinLabel(4,"D(3b)");
  xx_ABCD.GetXaxis()->SetBinLabel(5,"C(2b)");
  xx_ABCD.GetXaxis()->SetBinLabel(6,"D(2b)");

  TH1D xx_metSig_A("xx_metSig_A","MET Significance (Signal Mass Window, 4b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_B("xx_metSig_B","MET Significance (Mass Sideband, 4b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_C3b("xx_metSig_C3b","MET Significance (Signal Mass Window, 3b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_D3b("xx_metSig_D3b","MET Significance (Mass Sideband, 3b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_C2b("xx_metSig_C2b","MET Significance (Signal Mass Window, 2b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH1D xx_metSig_D2b("xx_metSig_D2b","MET Significance (Mass Sideband, 2b);MET Significance;Events/19.4 fb^{-1}/10",5,metSigABCDBinEdges);
  TH2D xx_METToMETSigRatio_vs_numLowPtPfCands("xx_METToMETSigRatio_vs_numLowPtPfCands","MET/MET Sig. vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,80,0.0,8.0);
  TH2D xx_MET_vs_numLowPtPfCands("xx_MET_vs_numLowPtPfCands","MET vs. Low PT PF Cands.;Low PT PF Cands.;MET [GeV]",16,-0.5,15.5,40,0.0,600.0);
  TH2D xx_METSig_vs_numLowPtPfCands("xx_METSig_vs_numLowPtPfCands","MET Sig vs. Low PT PF Cands.;Low PT PF Cands.;MET Sig.",16,-0.5,15.5,40,0.0,400.0);
  TH2D xx_METToMETSigRatio_vs_fracMETFromLowPt("xx_METToMETSigRatio_vs_fracMETFromLowPt","MET/MET Sig. vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET/MET Sig. [GeV]",40,0.0,1.0,80,0.0,8.0);
  TH2D xx_MET_vs_fracMETFromLowPt("xx_MET_vs_fracMETFromLowPt","MET vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET [GeV]",40,0.0,1.0,40,0.0,600.0);
  TH2D xx_METSig_vs_fracMETFromLowPt("xx_METSig_vs_fracMETFromLowPt","MET Sig vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET Sig.",40,0.0,1.0,40,0.0,400.0);

  //Two b-jet control region
  TH1D twoB_nm1_thirdBTag("twoB_nm1_thirdBTag", "N-1 Third B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D twoB_nm1_fourthBTag("twoB_nm1_fourthBTag", "N-1 Fourth B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D twoB_nm1_avgHiggsMass("twoB_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D twoB_nm1_higgsMassDiff("twoB_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D twoB_nm1_maxDeltaR("twoB_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D twoB_nm1_metSig("twoB_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D twoB_nm1_met("twoB_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25 GeV", 40, 0.0, 1000.0);
  TH2D twoB_METOverSqrtHT_vs_METSig("twoB_METOverSqrtHT_vs_METSig", "MET/#sqrt{HT} vs. MET Sig.;MET Significance; MET/#sqrt{HT} [#sqrt{GeV}]", 40, 0.0, 200.0, 24, 0.0, 12.0);
  TH1D twoB_METOverSqrtHT_Over_METSig("twoB_METOverSqrtHT_Over_METSig","(MET/#sqrt{HT})/MET Sig;(MET/#sqrt{HT})/MET Sig [#sqrt{GeV}];Events/19.4 fb^{-1}/0.01",30,0.0,0.3);
  TH1D twoB_MET("twoB_MET","MET;MET [GeV];Events/19.4 fb^{-1}/10 GeV",40,0.0,400.0);

  TH2D twoB_MET_vs_minDeltaPhi("twoB_MET_vs_minDeltaPhi","MET vs. Min. #Delta#phi;Min. #Delta#Phi;MET [GeV]",15,0.0,3.0,15,0.0,450.0);
  TH2D twoB_METSig_vs_minDeltaPhi("twoB_METSig_vs_minDeltaPhi","MET Significance vs. Min. #Delta#phi;Min. #Delta#Phi;MET Significance",15,0.0,3.0,15,0.0,300.0);
  TH2D twoB_METToMETSigRatio_vs_numLowPtPfCands("twoB_METToMETSigRatio_vs_numLowPtPfCands","MET/MET Sig. vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,80,0.0,8.0);
  TH2D twoB_MET_vs_numLowPtPfCands("twoB_MET_vs_numLowPtPfCands","MET vs. Low PT PF Cands.;Low PT PF Cands.;MET [GeV]",16,-0.5,15.5,40,0.0,600.0);
  TH2D twoB_METSig_vs_numLowPtPfCands("twoB_METSig_vs_numLowPtPfCands","MET Sig vs. Low PT PF Cands.;Low PT PF Cands.;MET Sig.",16,-0.5,15.5,40,0.0,400.0);
  TH2D twoB_METToMETSigRatio_vs_fracMETFromLowPt("twoB_METToMETSigRatio_vs_fracMETFromLowPt","MET/MET Sig. vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET/MET Sig. [GeV]",40,0.0,1.0,80,0.0,8.0);
  TH2D twoB_MET_vs_fracMETFromLowPt("twoB_MET_vs_fracMETFromLowPt","MET vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET [GeV]",40,0.0,1.0,40,0.0,600.0);
  TH2D twoB_METSig_vs_fracMETFromLowPt("twoB_METSig_vs_fracMETFromLowPt","MET Sig vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET Sig.",40,0.0,1.0,40,0.0,400.0);

  //Light Higgs control region
  TH1D lightHiggs_nm1_thirdBTag("lightHiggs_nm1_thirdBTag", "N-1 Third B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D lightHiggs_nm1_fourthBTag("lightHiggs_nm1_fourthBTag", "N-1 Fourth B-Tag;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D lightHiggs_nm1_avgHiggsMass("lightHiggs_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D lightHiggs_nm1_higgsMassDiff("lightHiggs_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D lightHiggs_nm1_maxDeltaR("lightHiggs_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D lightHiggs_nm1_metSig("lightHiggs_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D lightHiggs_nm1_met("lightHiggs_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH2D lightHiggs_METOverSqrtHT_vs_METSig("lightHiggs_METOverSqrtHT_vs_METSig", "MET/#sqrt{HT} vs. MET Sig.;MET Significance; MET/#sqrt{HT} [#sqrt{GeV}]", 40, 0.0, 200.0, 24, 0.0, 12.0);
  TH1D lightHiggs_METOverSqrtHT_Over_METSig("lightHiggs_METOverSqrtHT_Over_METSig","(MET/#sqrt{HT})/MET Sig;(MET/#sqrt{HT})/MET Sig [#sqrt{GeV}];Events/19.4 fb^{-1}/0.01",30,0.0,0.3);
  TH1D lightHiggs_MET("lightHiggs_MET","MET;MET [GeV];Events/19.4 fb^{-1}/25 GeV",40,0.0,400.0);

  TH2D lightHiggs_MET_vs_minDeltaPhi("lightHiggs_MET_vs_minDeltaPhi","MET vs. Min. #Delta#phi;Min. #Delta#Phi;MET [GeV]",15,0.0,3.0,15,0.0,450.0);
  TH2D lightHiggs_METSig_vs_minDeltaPhi("lightHiggs_METSig_vs_minDeltaPhi","MET Significance vs. Min. #Delta#phi;Min. #Delta#Phi;MET Significance",15,0.0,3.0,15,0.0,300.0);
  TH2D lightHiggs_METToMETSigRatio_vs_numLowPtPfCands("lightHiggs_METToMETSigRatio_vs_numLowPtPfCands","MET/MET Sig. vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,80,0.0,8.0);
  TH2D lightHiggs_MET_vs_numLowPtPfCands("lightHiggs_MET_vs_numLowPtPfCands","MET vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,40,0.0,600.0);
  TH2D lightHiggs_METSig_vs_numLowPtPfCands("lightHiggs_METSig_vs_numLowPtPfCands","MET Sig vs. Low PT PF Cands.;Low PT PF Cands.;MET/MET Sig. [GeV]",16,-0.5,15.5,40,0.0,400.0);
  TH2D lightHiggs_METToMETSigRatio_vs_fracMETFromLowPt("lightHiggs_METToMETSigRatio_vs_fracMETFromLowPt","MET/MET Sig. vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET/MET Sig. [GeV]",40,0.0,1.0,80,0.0,8.0);
  TH2D lightHiggs_MET_vs_fracMETFromLowPt("lightHiggs_MET_vs_fracMETFromLowPt","MET vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET [GeV]",40,0.0,1.0,40,0.0,600.0);
  TH2D lightHiggs_METSig_vs_fracMETFromLowPt("lightHiggs_METSig_vs_fracMETFromLowPt","MET Sig vs. Frac. of MET from Low PT PF Cands;|MET(PF<20)|/|MET|;MET Sig.",40,0.0,1.0,40,0.0,400.0);

  TH1D sorb_MET("sorb_MET","S/#sqrt{B};MET [GeV];S/#sqrt{B}",1024,0.0,350.0);
  TH1D sorb_METSig("sorb_METSig","S/#sqrt{B};MET Significance;S/#sqrt{B}",1024,0.0,150.0);
  TH1D sorb_METSig2012("sorb_METSig2012","S/#sqrt{B};MET Significance 2012;S/#sqrt{B}",1024,0.0,150.0);
  TH1D sorb_METOverSqrtHT("sorb_METOverSqrtHT","S/#sqrt{B};MET/#sqrt{HT} [#sqrt{GeV}];S/#sqrt{B}",1024,0.0,15.0);

  TH1D npv_before("npv_before", "Num. Primary Vertices Before Pileup Reweighting;Vertices;Events/19.4 fb^{-1}", 61, -0.5, 60.5);
  TH1D npv_after("npv_after", "Num. Primary Vertices After Pileup Reweighting;Vertices;Events/19.4 fb^{-1}", 61, -0.5, 60.5);

  TH1D topPt_before("topPt_before", "Top p_{T} Before Top p_{T} Reweighting;top p_{T} [GeV];Events/5 GeV/19.4 fb^{-1}", 50, 0.0, 500.0);
  TH1D topPt_after("topPt_after", "Top p_{T} After Top p_{T} Reweighting;top p_{T} [GeV];Events/5 GeV/19.4 fb^{-1}", 50, 0.0, 500.0);

  TH1D xx_origin1("xx_origin1", "Origin of Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D xx_origin2("xx_origin2", "Origin of Second Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D xx_origin3("xx_origin3", "Origin of Third Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D xx_origin4("xx_origin4", "Origin of Fourth Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  FixBinLabels(xx_origin1);
  FixBinLabels(xx_origin2);
  FixBinLabels(xx_origin3);
  FixBinLabels(xx_origin4);

  TH1D xx_beta1("xx_beta1", "beta for Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D xx_beta2("xx_beta2", "beta for Second Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D xx_beta3("xx_beta3", "beta for Third Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D xx_beta4("xx_beta4", "beta for Fourth Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);

  TH1D yy_origin1("yy_origin1", "Origin of Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D yy_origin2("yy_origin2", "Origin of Second Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D yy_origin3("yy_origin3", "Origin of Third Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  TH1D yy_origin4("yy_origin4", "Origin of Fourth Highest CSV Jet;Origin;Events/19.4 fb^{-1}", 19, 0.5, 19.5);
  FixBinLabels(yy_origin1);
  FixBinLabels(yy_origin2);
  FixBinLabels(yy_origin3);
  FixBinLabels(yy_origin4);

  TH1D yy_beta1("yy_beta1", "beta for Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D yy_beta2("yy_beta2", "beta for Second Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D yy_beta3("yy_beta3", "beta for Third Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);
  TH1D yy_beta4("yy_beta4", "beta for Fourth Highest CSV Jet;beta;Events/0.05/19.4 fb^{-1}", 20, 0.0, 1.0);

  TH1D xx_Sbins_4bSig("xx_Sbins_4bSig", "MET Significance Counts (4b, mass signal window);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_4bSB("yy_Sbins_4bSB", "MET Significance Counts (4b, mass sideband);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_3bSig("yy_Sbins_3bSig", "MET Significance Counts (3b, mass signal window);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_3bSB("yy_Sbins_3bSB", "MET Significance Counts (3b, mass sideband);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_2bSig("yy_Sbins_2bSig", "MET Significance Counts (2b, mass signal window);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D yy_Sbins_2bSB("yy_Sbins_2bSB", "MET Significance Counts (2b, mass sideband);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  FixSbinLabels(xx_Sbins_4bSig);
  FixSbinLabels(yy_Sbins_4bSB);
  FixSbinLabels(yy_Sbins_3bSig);
  FixSbinLabels(yy_Sbins_3bSB);
  FixSbinLabels(yy_Sbins_2bSig);
  FixSbinLabels(yy_Sbins_2bSB);

  TH1D sl_Sbins_4bSig("sl_Sbins_4bSig", "MET Significance Counts (4b, mass signal window, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_4bSB("sl_Sbins_4bSB", "MET Significance Counts (4b, mass sideband, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_3bSig("sl_Sbins_3bSig", "MET Significance Counts (3b, mass signal window, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_3bSB("sl_Sbins_3bSB", "MET Significance Counts (3b, mass sideband, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_2bSig("sl_Sbins_2bSig", "MET Significance Counts (2b, mass signal window, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  TH1D sl_Sbins_2bSB("sl_Sbins_2bSB", "MET Significance Counts (2b, mass sideband, 1l);S bin;Events/19.4 fb^{-1}", 4, 0.5, 4.5);
  FixSbinLabels(sl_Sbins_4bSig);
  FixSbinLabels(sl_Sbins_4bSB);
  FixSbinLabels(sl_Sbins_3bSig);
  FixSbinLabels(sl_Sbins_3bSB);
  FixSbinLabels(sl_Sbins_2bSig);
  FixSbinLabels(sl_Sbins_2bSB);

  TH1D sbin1_4b_sb_nm1_thirdCSV("sbin1_4b_sb_nm1_thirdCSV", "N-1 Third CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D sbin1_4b_sb_nm1_fourthCSV("sbin1_4b_sb_nm1_fourthCSV", "N-1 Fourth CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D sbin1_4b_sb_nm1_avgHiggsMass("sbin1_4b_sb_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D sbin1_4b_sb_nm1_higgsMassDiff("sbin1_4b_sb_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D sbin1_4b_sb_nm1_maxDeltaR("sbin1_4b_sb_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D sbin1_4b_sb_nm1_metSig("sbin1_4b_sb_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D sbin1_4b_sb_nm1_met("sbin1_4b_sb_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH1D sbin1_4b_sb_nm1_numJets("sbin1_4b_sb_nm1_numJets", "N-1 Number of Jets;Number of Jets;Events/19.4 fb^{-1}", 16, -0.5, 15.5);
  TH1D sbin1_4b_sb_nm1_minDeltaPhi("sbin1_4b_sb_nm1_minDeltaPhi", "N-1 min. #Delta#phi;min. #Delta#phi;Events/19.4 fb^{-1}", 30, 0.0, 4.0*atan(1.0));
  
  TH1D single_lepton_nm1_thirdCSV("single_lepton_nm1_thirdCSV", "N-1 Third CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D single_lepton_nm1_fourthCSV("single_lepton_nm1_fourthCSV", "N-1 Fourth CSV;CSV;Events/19.4 fb^{-1}/0.04", 25, 0.0, 1.0);
  TH1D single_lepton_nm1_avgHiggsMass("single_lepton_nm1_avgHiggsMass", "N-1 Avg. Higgs Mass;Avg. Higgs Mass [GeV];Events/19.4 fb^{-1}/10 GeV", 30, 0, 300.0);
  TH1D single_lepton_nm1_higgsMassDiff("single_lepton_nm1_higgsMassDiff", "N-1 Higgs Mass Diff.;Higgs Mass Diff [GeV];Events/19.4 fb^{-1}/5 GeV", 20, 0, 100.0);
  TH1D single_lepton_nm1_maxDeltaR("single_lepton_nm1_maxDeltaR", "N-1 Max #Delta R;Max #Delta R;Events/19.4 fb^{-1}/0.2", 25, 0.0, 5.0);
  TH1D single_lepton_nm1_metSig("single_lepton_nm1_metSig", "N-1 MET Significance;MET Significance;Events/19.4 fb^{-1}/10", 30, 0.0, 300.0);
  TH1D single_lepton_nm1_met("single_lepton_nm1_met", "N-1 MET;MET [GeV];Events/19.4 fb^{-1}/25", 40, 0.0, 1000.0);
  TH1D single_lepton_nm1_numJets("single_lepton_nm1_numJets", "N-1 Number of Jets;Number of Jets;Events/19.4 fb^{-1}", 16, -0.5, 15.5);
  TH1D single_lepton_nm1_minDeltaPhi("single_lepton_nm1_minDeltaPhi", "N-1 min. #Delta#phi;min. #Delta#phi;Events/19.4 fb^{-1}", 30, 0.0, 4.0*atan(1.0));
  
  Timer timer(GetTotalEntries());
  timer.Start();
  for(int i(0); i<GetTotalEntries(); ++i){
    if(i%1000==0 && i!=0){
      timer.PrintRemainingTime();
    }
    timer.Iterate();
    GetEntry(i);

    const bool isRealData(sampleName.find("Run2012")!=std::string::npos);
    const bool isttbar(sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos);
    const double localWeight((isRealData?1.0:GetPUWeight(lumiWeights))*scaleFactor*(true && isttbar?GetTopPtWeight():1.0)*GetSbinWeight()*(HasGluonSplitting()?1.0:1.0));
    const double nonpileupweight(scaleFactor*(true && isttbar?GetTopPtWeight():1.0)*(HasGluonSplitting()?1.0:1.0)*GetSbinWeight());
    const double notopweight((isRealData?1.0:GetPUWeight(lumiWeights))*scaleFactor*GetSbinWeight());

    if(!PassesJSONCut()) continue;
    if(!PassesTChiZHMassCut()) continue;

    std::pair<std::set<EventNumber>::iterator, bool> returnVal(eventList.insert(EventNumber(run, event, lumiblock)));
    if(!returnVal.second) continue;
    ++startCount;
    startCountWeighted+=localWeight;

    if(PassesPVCut() && PassesJet2PtCut() && Passes2CSVTCut() && PassesMETSig30Cut() && PassesMETCleaningCut() && PassesTriggerCut() && PassesNumJetsCut() && PassesMinDeltaPhiCut() && PassesLeptonVetoCut() && PassesIsoTrackVetoCut() && PassesDRCut()){
      double npv(-1.0);
      for(unsigned int bc(0); bc<PU_bunchCrossing->size(); ++bc){
	if(PU_bunchCrossing->at(bc)==0){
	  npv = PU_TrueNumInteractions->at(bc);
	}
      }
      npv_before.Fill(npv, nonpileupweight);
      npv_after.Fill(npv, localWeight);
      double topPt(-1.0);
      for ( unsigned int mc=0; mc< mc_doc_id->size(); mc++ ) {
	if ( mc_doc_id->at(mc)==6 ) { topPt = mc_doc_pt->at(mc); break; }
      }
      topPt_before.Fill(topPt, notopweight);
      topPt_after.Fill(topPt, localWeight);
    }

    if(!betaUpToDate) GetBeta();
    if(PassesRegionACut()){
      const std::vector<std::pair<int,int> > bo(GetBOrigins());
      xx_origin1.Fill(GetType(bo.at(0)),localWeight);
      xx_origin2.Fill(GetType(bo.at(1)),localWeight);
      xx_origin3.Fill(GetType(bo.at(2)),localWeight);
      xx_origin4.Fill(GetType(bo.at(3)),localWeight);
      std::vector<std::pair<double, double> > CSV_beta(0);
      for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
	if(isGoodJet(jet, true, 20.0, 2.4, false)){
	  CSV_beta.push_back(std::pair<double,double>(beta.at(jet),
						      jets_AK5PF_btag_secVertexCombined->at(jet)));
	}
      }
      std::sort(CSV_beta.begin(), CSV_beta.end(), std::greater<std::pair<double, double> >());
      xx_beta1.Fill(CSV_beta.at(0).first,localWeight);
      xx_beta2.Fill(CSV_beta.at(1).first,localWeight);
      xx_beta3.Fill(CSV_beta.at(2).first,localWeight);
      xx_beta4.Fill(CSV_beta.at(3).first,localWeight);
    }

    if(PassesRegionBCut()){
      const std::vector<std::pair<int,int> > bo(GetBOrigins());
      yy_origin1.Fill(GetType(bo.at(0)),localWeight);
      yy_origin2.Fill(GetType(bo.at(1)),localWeight);
      yy_origin3.Fill(GetType(bo.at(2)),localWeight);
      yy_origin4.Fill(GetType(bo.at(3)),localWeight);
      std::vector<std::pair<double, double> > CSV_beta(0);
      for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
	if(isGoodJet(jet, true, 20.0, 2.4, false)){
	  CSV_beta.push_back(std::pair<double,double>(beta.at(jet),
						      jets_AK5PF_btag_secVertexCombined->at(jet)));
	}
      }
      std::sort(CSV_beta.begin(), CSV_beta.end(), std::greater<std::pair<double, double> >());
      yy_beta1.Fill(CSV_beta.at(0).first,localWeight);
      yy_beta2.Fill(CSV_beta.at(1).first,localWeight);
      yy_beta3.Fill(CSV_beta.at(2).first,localWeight);
      yy_beta4.Fill(CSV_beta.at(3).first,localWeight);
    }

    double SbinToFill(0.0);
    if(pfmets_fullSignif>30.0 && pfmets_fullSignif<50.0) SbinToFill=1.0;
    if(pfmets_fullSignif>50.0 && pfmets_fullSignif<100.0) SbinToFill=2.0;
    if(pfmets_fullSignif>100.0 && pfmets_fullSignif<150.0) SbinToFill=3.0;
    if(pfmets_fullSignif>150.0) SbinToFill=4.0;

    if(PassesRegionACut()) xx_Sbins_4bSig.Fill(SbinToFill, localWeight);
    if(PassesRegionBCut()) yy_Sbins_4bSB.Fill(SbinToFill, localWeight);
    if(PassesRegionC3bCut()) yy_Sbins_3bSig.Fill(SbinToFill, localWeight);
    if(PassesRegionD3bCut()) yy_Sbins_3bSB.Fill(SbinToFill, localWeight);
    if(PassesRegionC2bCut()) yy_Sbins_2bSig.Fill(SbinToFill, localWeight);
    if(PassesRegionD2bCut()) yy_Sbins_2bSB.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionACut()) sl_Sbins_4bSig.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionBCut()) sl_Sbins_4bSB.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionC3bCut()) sl_Sbins_3bSig.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionD3bCut()) sl_Sbins_3bSB.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionC2bCut()) sl_Sbins_2bSig.Fill(SbinToFill, localWeight);
    if(PassesSingleLeptonRegionD2bCut()) sl_Sbins_2bSB.Fill(SbinToFill, localWeight);

    const uint_least32_t failure_code(GetCutFailCode());
    uint_least32_t masked_fc(failure_code & ~kMETSig100 & ~kMETSig150);
    if((masked_fc & ~kNumJets) == kMETSig50) sbin1_4b_sb_nm1_numJets.Fill(GetNumGoodJets(),localWeight);
    if((masked_fc & ~kMinDeltaPhi) == kMETSig50) sbin1_4b_sb_nm1_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),localWeight);
    if((masked_fc & ~k3rdBTag & ~k4thBTag) == kMETSig50) sbin1_4b_sb_nm1_thirdCSV.Fill(GetHighestCSV(3),localWeight);
    if((masked_fc & ~k4thBTag) == kMETSig50) sbin1_4b_sb_nm1_fourthCSV.Fill(GetHighestCSV(4),localWeight);
    const std::pair<double,double> higgsMasses(GetHiggsMasses());
    if((masked_fc & ~kHiggsMassDiff) == kMETSig50) sbin1_4b_sb_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
    if((masked_fc & ~kHiggsAvgMass) == kMETSig50) sbin1_4b_sb_nm1_higgsMassDiff.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
    if((masked_fc & ~kDeltaR) == kMETSig50) sbin1_4b_sb_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
    if(masked_fc == kMETSig50) sbin1_4b_sb_nm1_met.Fill(pfmets_et->at(0),localWeight);
    if(masked_fc == kMETSig50) sbin1_4b_sb_nm1_metSig.Fill(pfmets_fullSignif,localWeight);

    masked_fc=(failure_code & ~kMETSig50 & ~kMETSig100 & ~kMETSig150);
    if(PassesSingleLeptonCut()){
      if((masked_fc & ~kNumJets) == kGood) single_lepton_nm1_numJets.Fill(GetNumGoodJets(),localWeight);
      if((masked_fc & ~kMinDeltaPhi) == kGood) single_lepton_nm1_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),localWeight);
      if((masked_fc & ~k3rdBTag & ~k4thBTag) == kGood) single_lepton_nm1_thirdCSV.Fill(GetHighestCSV(3),localWeight);
      if((masked_fc & ~k4thBTag) == kGood) single_lepton_nm1_fourthCSV.Fill(GetHighestCSV(4),localWeight);
      if((masked_fc & ~kHiggsMassDiff) == kGood) single_lepton_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
      if((masked_fc & ~kHiggsAvgMass) == kGood) single_lepton_nm1_higgsMassDiff.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
      if((masked_fc & ~kDeltaR) == kGood) single_lepton_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
      if((masked_fc & ~kMETSig30) == kGood) single_lepton_nm1_met.Fill(pfmets_et->at(0),localWeight);
      if((masked_fc & ~kMETSig30) == kGood) single_lepton_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
    }

    if(PassesPVCut() && PassesMETCleaningCut() && PassesTriggerCut() && PassesNumJetsCut() && Passes2CSVTCut() && PassesJet2PtCut() && PassesMinDeltaPhiCut() && PassesLeptonVetoCut() && PassesIsoTrackVetoCut()){

      if(PassesBTaggingCut() && PassesHiggsMassCut() && PassesDRCut()){
	//MET S/sqrt(B) plots go here!
	sorb_MET.Fill(pfmets_et->at(0), localWeight);
	sorb_METSig.Fill(pfmets_fullSignif, localWeight);
	if(GetcfAVersion()>=69){
	  sorb_METSig.Fill(pfmets_fullSignif, localWeight);
	}
	sorb_METOverSqrtHT.Fill(pfmets_et->at(0)/sqrt(GetHT()),localWeight);
      }

      if(GetNumCSVMJets()<=2 && GetNumCSVLJets()<=3){
	//2b control region plots go here!
	if(PassesMETSig30Cut() && PassesHiggsMassCut() && PassesDRCut()){
	  twoB_nm1_thirdBTag.Fill(GetHighestCSV(3),localWeight);
	  twoB_nm1_fourthBTag.Fill(GetHighestCSV(4),localWeight);
	}
	if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesDRCut()){
	  twoB_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
	}
	if(PassesMETSig30Cut() && PassesHiggsAvgMassCut() && PassesDRCut()){
	  twoB_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
	}
	if(PassesMETSig30Cut() && PassesHiggsMassCut()){
	  twoB_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
	}
	if(PassesHiggsMassCut() && PassesDRCut()){
	  twoB_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
	  twoB_nm1_met.Fill(pfmets_et->at(0),localWeight);
	  twoB_METOverSqrtHT_vs_METSig.Fill(pfmets_fullSignif, pfmets_et->at(0)/sqrt(GetHT()), localWeight);
	  twoB_METOverSqrtHT_Over_METSig.Fill((pfmets_et->at(0)/sqrt(GetHT()))/pfmets_fullSignif, localWeight);
	  twoB_MET.Fill(pfmets_et->at(0),localWeight);
	  twoB_MET_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),pfmets_et->at(0),localWeight);
	  twoB_METSig_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),pfmets_fullSignif,localWeight);
	  twoB_METToMETSigRatio_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
	  twoB_MET_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0),localWeight);
	  twoB_METSig_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_fullSignif,localWeight);
	  twoB_METToMETSigRatio_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
	  twoB_MET_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0),localWeight);
	  twoB_METSig_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_fullSignif,localWeight);
	}
      }

      if(0.5*(higgsMasses.first+higgsMasses.second)<75.0){
	//Higgs mass control region plots go here!
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesDRCut()){
          lightHiggs_nm1_thirdBTag.Fill(GetHighestCSV(3),localWeight);
          lightHiggs_nm1_fourthBTag.Fill(GetHighestCSV(4),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesDRCut() && PassesBTaggingCut()){
          lightHiggs_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
        }
        if(PassesMETSig30Cut() && PassesDRCut() && PassesBTaggingCut()){
          lightHiggs_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
        }
        if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesBTaggingCut()){
          lightHiggs_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
        }
        if(PassesHiggsMassDiffCut() && PassesDRCut() && PassesBTaggingCut()){
	  lightHiggs_METOverSqrtHT_vs_METSig.Fill(pfmets_fullSignif, pfmets_et->at(0)/sqrt(GetHT()), localWeight);
          lightHiggs_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
          lightHiggs_nm1_met.Fill(pfmets_et->at(0),localWeight);
	  lightHiggs_METOverSqrtHT_Over_METSig.Fill((pfmets_et->at(0)/sqrt(GetHT()))/pfmets_fullSignif, localWeight);
	  lightHiggs_MET.Fill(pfmets_et->at(0),localWeight);
	  lightHiggs_MET_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),pfmets_et->at(0),localWeight);
	  lightHiggs_METSig_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),pfmets_fullSignif,localWeight);
	  lightHiggs_METToMETSigRatio_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
	  lightHiggs_MET_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0),localWeight);
	  lightHiggs_METSig_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_fullSignif,localWeight);
	  lightHiggs_METToMETSigRatio_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
	  lightHiggs_MET_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0),localWeight);
	  lightHiggs_METSig_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_fullSignif,localWeight);
        }
      }

      //Signal region plots go here!
      if(PassesMETSig30Cut() && PassesHiggsMassCut() && PassesDRCut()){
	xx_nm1_thirdBTag.Fill(GetHighestCSV(3),localWeight);
	if(GetNumCSVMJets()>=3){
	  xx_nm1_fourthBTag.Fill(GetHighestCSV(4),localWeight);
	}
      }
      if(PassesMETSig30Cut() && PassesHiggsMassDiffCut() && PassesBTaggingCut() && PassesDRCut()){
	xx_nm1_avgHiggsMass.Fill(0.5*(higgsMasses.first+higgsMasses.second),localWeight);
      }
      if(PassesMETSig30Cut() && PassesHiggsAvgMassCut() && PassesBTaggingCut() && PassesDRCut()){
	xx_nm1_higgsMassDiff.Fill(fabs(higgsMasses.first-higgsMasses.second),localWeight);
      }
      if(PassesMETSig30Cut() && PassesHiggsMassCut() && PassesBTaggingCut()){
	xx_nm1_maxDeltaR.Fill(GetMaxDR(),localWeight);
      }
      if(PassesHiggsMassCut() && PassesBTaggingCut() && PassesDRCut()){
	xx_METOverSqrtHT_vs_METSig.Fill(pfmets_fullSignif, pfmets_et->at(0)/sqrt(GetHT()), localWeight);
	xx_METOverSqrtHT_Over_METSig.Fill((pfmets_et->at(0)/sqrt(GetHT()))/pfmets_fullSignif, localWeight);
	xx_nm1_metSig.Fill(pfmets_fullSignif,localWeight);
	xx_nm1_met.Fill(pfmets_et->at(0),localWeight);
	xx_MET.Fill(pfmets_et->at(0),localWeight);
	xx_MET_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),pfmets_et->at(0),localWeight);
	xx_METSig_vs_minDeltaPhi.Fill(GetMinDeltaPhiMET(3),pfmets_fullSignif,localWeight);
      }
    }

    xx_METToMETSigRatio_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
    xx_MET_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_et->at(0),localWeight);
    xx_METSig_vs_numLowPtPfCands.Fill(static_cast<double>(GetNumLowPtPfCands(20.0)),pfmets_fullSignif,localWeight);
    xx_METToMETSigRatio_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0)/pfmets_fullSignif,localWeight);
    xx_MET_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_et->at(0),localWeight);
    xx_METSig_vs_fracMETFromLowPt.Fill(GetMETOfLowPtPfCands(20.0)/pfmets_et->at(0),pfmets_fullSignif,localWeight);

    unsigned int this_sbin(0);
    if(pfmets_fullSignif>30.0 && pfmets_fullSignif<50.0){
      this_sbin=1;
    }else if(pfmets_fullSignif>50.0 && pfmets_fullSignif<100.0){
      this_sbin=2;
    }else if(pfmets_fullSignif>100.0 && pfmets_fullSignif<150.0){
      this_sbin=3;
    }else if(pfmets_fullSignif>150.0){
      this_sbin=4;
    }
    if(PassesRegionACut()){
      switch(this_sbin){
      case 1: ++Asbin1Count;
	Asbin1CountWeighted+=localWeight;
	break;
      case 2: ++Asbin2Count;
	Asbin2CountWeighted+=localWeight;
	break;
      case 3: ++Asbin3Count;
	Asbin3CountWeighted+=localWeight;
	break;
      case 4: ++Asbin4Count;
	Asbin4CountWeighted+=localWeight;
      }
      ++ACount;
      ACountWeighted+=localWeight;
      AUncert+=localWeight*localWeight;
      xx_ABCD.Fill(0.0,localWeight);
      xx_metSig_A.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionBCut()){
      switch(this_sbin){
      case 1: ++Bsbin1Count;
	Bsbin1CountWeighted+=localWeight;
	break;
      case 2: ++Bsbin2Count;
	Bsbin2CountWeighted+=localWeight;
	break;
      case 3: ++Bsbin3Count;
	Bsbin3CountWeighted+=localWeight;
	break;
      case 4: ++Bsbin4Count;
	Bsbin4CountWeighted+=localWeight;
      }
      ++BCount;
      BCountWeighted+=localWeight;
      BUncert+=localWeight*localWeight;
      xx_ABCD.Fill(1.0,localWeight);
      xx_metSig_B.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionC3bCut()){
      switch(this_sbin){
      case 1: ++C3bsbin1Count;
	C3bsbin1CountWeighted+=localWeight;
	break;
      case 2: ++C3bsbin2Count;
	C3bsbin2CountWeighted+=localWeight;
	break;
      case 3: ++C3bsbin3Count;
	C3bsbin3CountWeighted+=localWeight;
	break;
      case 4: ++C3bsbin4Count;
	C3bsbin4CountWeighted+=localWeight;
      }
      ++C3bCount;
      C3bCountWeighted+=localWeight;
      C3bUncert+=localWeight*localWeight;
      xx_ABCD.Fill(2.0,localWeight);
      xx_metSig_C3b.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionD3bCut()){
      switch(this_sbin){
      case 1: ++D3bsbin1Count;
	D3bsbin1CountWeighted+=localWeight;
	break;
      case 2: ++D3bsbin2Count;
	D3bsbin2CountWeighted+=localWeight;
	break;
      case 3: ++D3bsbin3Count;
	D3bsbin3CountWeighted+=localWeight;
	break;
      case 4: ++D3bsbin4Count;
	D3bsbin4CountWeighted+=localWeight;
      }
      ++D3bCount;
      D3bCountWeighted+=localWeight;
      D3bUncert+=localWeight*localWeight;
      xx_ABCD.Fill(3.0,localWeight);
      xx_metSig_D3b.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionC2bCut()){
      switch(this_sbin){
      case 1: ++C2bsbin1Count;
	C2bsbin1CountWeighted+=localWeight;
	break;
      case 2: ++C2bsbin2Count;
	C2bsbin2CountWeighted+=localWeight;
	break;
      case 3: ++C2bsbin3Count;
	C2bsbin3CountWeighted+=localWeight;
	break;
      case 4: ++C2bsbin4Count;
	C2bsbin4CountWeighted+=localWeight;
      }
      ++C2bCount;
      C2bCountWeighted+=localWeight;
      C2bUncert+=localWeight*localWeight;
      xx_ABCD.Fill(4.0,localWeight);
      xx_metSig_C2b.Fill(pfmets_fullSignif,localWeight);
    }
    if(PassesRegionD2bCut()){
      switch(this_sbin){
      case 1: ++D2bsbin1Count;
	D2bsbin1CountWeighted+=localWeight;
	break;
      case 2: ++D2bsbin2Count;
	D2bsbin2CountWeighted+=localWeight;
	break;
      case 3: ++D2bsbin3Count;
	D2bsbin3CountWeighted+=localWeight;
	break;
      case 4: ++D2bsbin4Count;
	D2bsbin4CountWeighted+=localWeight;
      }
      ++D2bCount;
      D2bCountWeighted+=localWeight;
      D2bUncert+=localWeight*localWeight;
      xx_ABCD.Fill(5.0,localWeight);
      xx_metSig_D2b.Fill(pfmets_fullSignif,localWeight);
    }

    if(PassesSingleLeptonRegionACut()){
      ++ASLCount;
      ASLCountWeighted+=localWeight;
      ASLUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++ASLsbin1Count;
	ASLsbin1CountWeighted+=localWeight;
	break;
      case 2: ++ASLsbin2Count;
	ASLsbin2CountWeighted+=localWeight;
	break;
      case 3: ++ASLsbin3Count;
	ASLsbin3CountWeighted+=localWeight;
	break;
      case 4: ++ASLsbin4Count;
	ASLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionBCut()){
      ++BSLCount;
      BSLCountWeighted+=localWeight;
      BSLUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++BSLsbin1Count;
	BSLsbin1CountWeighted+=localWeight;
	break;
      case 2: ++BSLsbin2Count;
	BSLsbin2CountWeighted+=localWeight;
	break;
      case 3: ++BSLsbin3Count;
	BSLsbin3CountWeighted+=localWeight;
	break;
      case 4: ++BSLsbin4Count;
	BSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionC3bCut()){
      ++C3bSLCount;
      C3bSLCountWeighted+=localWeight;
      C3bSLUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++C3bSLsbin1Count;
	C3bSLsbin1CountWeighted+=localWeight;
	break;
      case 2: ++C3bSLsbin2Count;
	C3bSLsbin2CountWeighted+=localWeight;
	break;
      case 3: ++C3bSLsbin3Count;
	C3bSLsbin3CountWeighted+=localWeight;
	break;
      case 4: ++C3bSLsbin4Count;
	C3bSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionD3bCut()){
      ++D3bSLCount;
      D3bSLCountWeighted+=localWeight;
      D3bSLUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++D3bSLsbin1Count;
	D3bSLsbin1CountWeighted+=localWeight;
	break;
      case 2: ++D3bSLsbin2Count;
	D3bSLsbin2CountWeighted+=localWeight;
	break;
      case 3: ++D3bSLsbin3Count;
	D3bSLsbin3CountWeighted+=localWeight;
	break;
      case 4: ++D3bSLsbin4Count;
	D3bSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionC2bCut()){
      ++C2bSLCount;
      C2bSLCountWeighted+=localWeight;
      C2bSLUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++C2bSLsbin1Count;
	C2bSLsbin1CountWeighted+=localWeight;
	break;
      case 2: ++C2bSLsbin2Count;
	C2bSLsbin2CountWeighted+=localWeight;
	break;
      case 3: ++C2bSLsbin3Count;
	C2bSLsbin3CountWeighted+=localWeight;
	break;
      case 4: ++C2bSLsbin4Count;
	C2bSLsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesSingleLeptonRegionD2bCut()){
      ++D2bSLCount;
      D2bSLCountWeighted+=localWeight;
      D2bSLUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++D2bSLsbin1Count;
	D2bSLsbin1CountWeighted+=localWeight;
	break;
      case 2: ++D2bSLsbin2Count;
	D2bSLsbin2CountWeighted+=localWeight;
	break;
      case 3: ++D2bSLsbin3Count;
	D2bSLsbin3CountWeighted+=localWeight;
	break;
      case 4: ++D2bSLsbin4Count;
	D2bSLsbin4CountWeighted+=localWeight;
      }
    }

    if(PassesInvertedDRRegionACut()){
      ++ADRInvCount;
      ADRInvCountWeighted+=localWeight;
      ADRInvUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++ADRInvsbin1Count;
	ADRInvsbin1CountWeighted+=localWeight;
	break;
      case 2: ++ADRInvsbin2Count;
	ADRInvsbin2CountWeighted+=localWeight;
	break;
      case 3: ++ADRInvsbin3Count;
	ADRInvsbin3CountWeighted+=localWeight;
	break;
      case 4: ++ADRInvsbin4Count;
	ADRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionBCut()){
      ++BDRInvCount;
      BDRInvCountWeighted+=localWeight;
      BDRInvUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++BDRInvsbin1Count;
	BDRInvsbin1CountWeighted+=localWeight;
	break;
      case 2: ++BDRInvsbin2Count;
	BDRInvsbin2CountWeighted+=localWeight;
	break;
      case 3: ++BDRInvsbin3Count;
	BDRInvsbin3CountWeighted+=localWeight;
	break;
      case 4: ++BDRInvsbin4Count;
	BDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionC3bCut()){
      ++C3bDRInvCount;
      C3bDRInvCountWeighted+=localWeight;
      C3bDRInvUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++C3bDRInvsbin1Count;
	C3bDRInvsbin1CountWeighted+=localWeight;
	break;
      case 2: ++C3bDRInvsbin2Count;
	C3bDRInvsbin2CountWeighted+=localWeight;
	break;
      case 3: ++C3bDRInvsbin3Count;
	C3bDRInvsbin3CountWeighted+=localWeight;
	break;
      case 4: ++C3bDRInvsbin4Count;
	C3bDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionD3bCut()){
      ++D3bDRInvCount;
      D3bDRInvCountWeighted+=localWeight;
      D3bDRInvUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++D3bDRInvsbin1Count;
	D3bDRInvsbin1CountWeighted+=localWeight;
	break;
      case 2: ++D3bDRInvsbin2Count;
	D3bDRInvsbin2CountWeighted+=localWeight;
	break;
      case 3: ++D3bDRInvsbin3Count;
	D3bDRInvsbin3CountWeighted+=localWeight;
	break;
      case 4: ++D3bDRInvsbin4Count;
	D3bDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionC2bCut()){
      ++C2bDRInvCount;
      C2bDRInvCountWeighted+=localWeight;
      C2bDRInvUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++C2bDRInvsbin1Count;
	C2bDRInvsbin1CountWeighted+=localWeight;
	break;
      case 2: ++C2bDRInvsbin2Count;
	C2bDRInvsbin2CountWeighted+=localWeight;
	break;
      case 3: ++C2bDRInvsbin3Count;
	C2bDRInvsbin3CountWeighted+=localWeight;
	break;
      case 4: ++C2bDRInvsbin4Count;
	C2bDRInvsbin4CountWeighted+=localWeight;
      }
    }
    if(PassesInvertedDRRegionD2bCut()){
      ++D2bDRInvCount;
      D2bDRInvCountWeighted+=localWeight;
      D2bDRInvUncert+=localWeight*localWeight;
      switch(this_sbin){
      case 1: ++D2bDRInvsbin1Count;
	D2bDRInvsbin1CountWeighted+=localWeight;
	break;
      case 2: ++D2bDRInvsbin2Count;
	D2bDRInvsbin2CountWeighted+=localWeight;
	break;
      case 3: ++D2bDRInvsbin3Count;
	D2bDRInvsbin3CountWeighted+=localWeight;
	break;
      case 4: ++D2bDRInvsbin4Count;
	D2bDRInvsbin4CountWeighted+=localWeight;
      }
    }

    if(!PassesPVCut()) continue;
    ++PVCount;
    PVCountWeighted+=localWeight;
    if(!PassesJet2PtCut()) continue;
    ++jet2PtCount;
    jet2PtCountWeighted+=localWeight;
    if(!Passes2CSVTCut()) continue;
    ++CSVTCount;
    CSVTCountWeighted+=localWeight;
    if(!PassesMETSig30Cut()) continue;
    ++METSig30Count;
    METSig30CountWeighted+=localWeight;
    if(!PassesMETCleaningCut()) continue;
    ++METCleaningCount;
    METCleaningCountWeighted+=localWeight;
    if(!PassesTriggerCut()) continue;
    ++TriggerCount;
    TriggerCountWeighted+=localWeight;
    if(!PassesNumJetsCut()) continue;
    ++numJetsCount;
    numJetsCountWeighted+=localWeight;
    if(!PassesMinDeltaPhiCut()) continue;
    ++minDeltaPhiCount;
    minDeltaPhiCountWeighted+=localWeight;
    if(!PassesLeptonVetoCut()) continue;
    ++leptonVetoCount;
    leptonVetoCountWeighted+=localWeight;
    if(!PassesIsoTrackVetoCut()) continue;
    ++isoTrackVetoCount;
    isoTrackVetoCountWeighted+=localWeight;
    if(!PassesBTaggingCut()) continue;
    ++bTagCount;
    bTagCountWeighted+=localWeight;
    if(!PassesHiggsMassCut()) continue;
    ++higgsCount;
    higgsCountWeighted+=localWeight;
    if(!PassesDRCut()) continue;
    ++DRCount;
    DRCountWeighted+=localWeight;
    if(!PassesMETSig50Cut()) continue;
    ++METSig50Count;
    METSig50CountWeighted+=localWeight;
    if(!PassesMETSig100Cut()) continue;
    ++METSig100Count;
    METSig100CountWeighted+=localWeight;
    if(!PassesMETSig150Cut()) continue;
    ++METSig150Count;
    METSig150CountWeighted+=localWeight;
  }

  for(int i(1); i<=xx_metSig_A.GetNbinsX(); ++i){
    xx_metSig_A.SetBinContent(i,xx_metSig_A.GetBinContent(i)*10.0/xx_metSig_A.GetBinWidth(i));
    xx_metSig_B.SetBinContent(i,xx_metSig_B.GetBinContent(i)*10.0/xx_metSig_B.GetBinWidth(i));
    xx_metSig_C3b.SetBinContent(i,xx_metSig_C3b.GetBinContent(i)*10.0/xx_metSig_C3b.GetBinWidth(i));
    xx_metSig_D3b.SetBinContent(i,xx_metSig_D3b.GetBinContent(i)*10.0/xx_metSig_D3b.GetBinWidth(i));
    xx_metSig_C2b.SetBinContent(i,xx_metSig_C2b.GetBinContent(i)*10.0/xx_metSig_C2b.GetBinWidth(i));
    xx_metSig_D2b.SetBinContent(i,xx_metSig_D2b.GetBinContent(i)*10.0/xx_metSig_D2b.GetBinWidth(i));
  }

  xx_metSig_SigOverSB_4b.Divide(&xx_metSig_A,&xx_metSig_B);
  xx_metSig_SigOverSB_3b.Divide(&xx_metSig_C3b,&xx_metSig_D3b);
  xx_metSig_SigOverSB_2b.Divide(&xx_metSig_C2b,&xx_metSig_D2b);
  xx_metSig_SigOverSB_4b.SetLineColor(1);
  xx_metSig_SigOverSB_3b.SetLineColor(2);
  xx_metSig_SigOverSB_2b.SetLineColor(3);
  xx_metSig_SigOverSB_2b.Draw("e1p");
  xx_metSig_SigOverSB_3b.Draw("e1psame");
  xx_metSig_SigOverSB_4b.Draw("e1psame");
  xx_metSig_SigOverSB.Write();

  double CtotCountWeighted(C3bCountWeighted+C2bCountWeighted);
  double DtotCountWeighted(D3bCountWeighted+D2bCountWeighted);
  double CtotUncert(C3bUncert+C2bUncert);
  double DtotUncert(D3bUncert+D2bUncert);

  double Asq(ACountWeighted*ACountWeighted);
  double Bsq(BCountWeighted*BCountWeighted);
  double C3bsq(C3bCountWeighted*C3bCountWeighted);
  double D3bsq(D3bCountWeighted*D3bCountWeighted);
  double C2bsq(C2bCountWeighted*C2bCountWeighted);
  double D2bsq(D2bCountWeighted*D2bCountWeighted);
  double Ctotsq(CtotCountWeighted*CtotCountWeighted);
  double Dtotsq(DtotCountWeighted*DtotCountWeighted);

  double kappa3b((ACountWeighted/BCountWeighted)/(C3bCountWeighted/D3bCountWeighted));
  double kappa2b((ACountWeighted/BCountWeighted)/(C2bCountWeighted/D2bCountWeighted));
  double kappatot((ACountWeighted/BCountWeighted)/(CtotCountWeighted/DtotCountWeighted));

  double numeratorA3b(AUncert*Bsq*C3bsq*D3bsq);
  double numeratorB3b(Asq*BUncert*C3bsq*D3bsq);
  double numeratorC3b(Asq*Bsq*C3bUncert*D3bsq);
  double numeratorD3b(Asq*Bsq*C3bsq*D3bUncert);
  double denominator3b(Bsq*Bsq*C3bsq*C3bsq);
  double kappaUncert3b((numeratorA3b+numeratorB3b+numeratorC3b+numeratorD3b)/denominator3b);

  double numeratorA2b(AUncert*Bsq*C2bsq*D2bsq);
  double numeratorB2b(Asq*BUncert*C2bsq*D2bsq);
  double numeratorC2b(Asq*Bsq*C2bUncert*D2bsq);
  double numeratorD2b(Asq*Bsq*C2bsq*D2bUncert);
  double denominator2b(Bsq*Bsq*C2bsq*C2bsq);
  double kappaUncert2b((numeratorA2b+numeratorB2b+numeratorC2b+numeratorD2b)/denominator2b);

  double numeratorAtot(AUncert*Bsq*Ctotsq*Dtotsq);
  double numeratorBtot(Asq*BUncert*Ctotsq*Dtotsq);
  double numeratorCtot(Asq*Bsq*CtotUncert*Dtotsq);
  double numeratorDtot(Asq*Bsq*Ctotsq*DtotUncert);
  double denominatortot(Bsq*Bsq*Ctotsq*Ctotsq);
  double kappaUncerttot((numeratorAtot+numeratorBtot+numeratorCtot+numeratorDtot)/denominatortot);

  AUncert=sqrt(AUncert);
  BUncert=sqrt(BUncert);
  C3bUncert=sqrt(C3bUncert);
  D3bUncert=sqrt(D3bUncert);
  C2bUncert=sqrt(C2bUncert);
  D2bUncert=sqrt(D2bUncert);
  CtotUncert=sqrt(CtotUncert);
  DtotUncert=sqrt(DtotUncert);

  TTree ABCD("ABCD","ABCD");
  ABCD.Branch("A", &ACount);
  ABCD.Branch("B", &BCount);
  ABCD.Branch("C3b", &C3bCount);
  ABCD.Branch("D3b", &D3bCount);
  ABCD.Branch("C2b", &C2bCount);
  ABCD.Branch("D2b", &D2bCount);
  ABCD.Branch("AWeighted", &ACountWeighted);
  ABCD.Branch("BWeighted", &BCountWeighted);
  ABCD.Branch("C3bWeighted", &C3bCountWeighted);
  ABCD.Branch("D3bWeighted", &D3bCountWeighted);
  ABCD.Branch("C2bWeighted", &C2bCountWeighted);
  ABCD.Branch("D2bWeighted", &D2bCountWeighted);
  ABCD.Branch("CtotWeighted", &CtotCountWeighted);
  ABCD.Branch("DtotWeighted", &DtotCountWeighted);
  ABCD.Branch("AUncert", &AUncert);
  ABCD.Branch("BUncert", &BUncert);
  ABCD.Branch("C3bUncert", &C3bUncert);
  ABCD.Branch("D3bUncert", &D3bUncert);
  ABCD.Branch("C2bUncert", &C2bUncert);
  ABCD.Branch("D2bUncert", &D2bUncert);
  ABCD.Branch("CtotUncert", &CtotUncert);
  ABCD.Branch("DtotUncert", &DtotUncert);
  ABCD.Branch("A_DRInv", &ADRInvCount);
  ABCD.Branch("B_DRInv", &BDRInvCount);
  ABCD.Branch("C3b_DRInv", &C3bDRInvCount);
  ABCD.Branch("D3b_DRInv", &D3bDRInvCount);
  ABCD.Branch("C2b_DRInv", &C2bDRInvCount);
  ABCD.Branch("D2b_DRInv", &D2bDRInvCount);
  ABCD.Branch("AWeighted_DRInv", &ADRInvCountWeighted);
  ABCD.Branch("BWeighted_DRInv", &BDRInvCountWeighted);
  ABCD.Branch("C3bWeighted_DRInv", &C3bDRInvCountWeighted);
  ABCD.Branch("D3bWeighted_DRInv", &D3bDRInvCountWeighted);
  ABCD.Branch("C2bWeighted_DRInv", &C2bDRInvCountWeighted);
  ABCD.Branch("D2bWeighted_DRInv", &D2bDRInvCountWeighted);
  ABCD.Branch("AUncert_DRInv", &ADRInvUncert);
  ABCD.Branch("BDRInvUncert_DRInv", &BDRInvUncert);
  ABCD.Branch("C3bUncert_DRInv", &C3bDRInvUncert);
  ABCD.Branch("D3bUncert_DRInv", &D3bDRInvUncert);
  ABCD.Branch("C2bUncert_DRInv", &C2bDRInvUncert);
  ABCD.Branch("D2bUncert_DRInv", &D2bDRInvUncert);
  ABCD.Branch("A_SL", &ASLCount);
  ABCD.Branch("B_SL", &BSLCount);
  ABCD.Branch("C3b_SL", &C3bSLCount);
  ABCD.Branch("D3b_SL", &D3bSLCount);
  ABCD.Branch("C2b_SL", &C2bSLCount);
  ABCD.Branch("D2b_SL", &D2bSLCount);
  ABCD.Branch("AWeighted_SL", &ASLCountWeighted);
  ABCD.Branch("BWeighted_SL", &BSLCountWeighted);
  ABCD.Branch("C3bWeighted_SL", &C3bSLCountWeighted);
  ABCD.Branch("D3bWeighted_SL", &D3bSLCountWeighted);
  ABCD.Branch("C2bWeighted_SL", &C2bSLCountWeighted);
  ABCD.Branch("D2bWeighted_SL", &D2bSLCountWeighted);
  ABCD.Branch("AUncert_SL", &ASLUncert);
  ABCD.Branch("BSLUncert_SL", &BSLUncert);
  ABCD.Branch("C3bUncert_SL", &C3bSLUncert);
  ABCD.Branch("D3bUncert_SL", &D3bSLUncert);
  ABCD.Branch("C2bUncert_SL", &C2bSLUncert);
  ABCD.Branch("D2bUncert_SL", &D2bSLUncert);
  ABCD.Branch("A_SLsbin1", &ASLsbin1Count);
  ABCD.Branch("B_SLsbin1", &BSLsbin1Count);
  ABCD.Branch("C3b_SLsbin1", &C3bSLsbin1Count);
  ABCD.Branch("D3b_SLsbin1", &D3bSLsbin1Count);
  ABCD.Branch("C2b_SLsbin1", &C2bSLsbin1Count);
  ABCD.Branch("D2b_SLsbin1", &D2bSLsbin1Count);
  ABCD.Branch("AWeighted_SLsbin1", &ASLsbin1CountWeighted);
  ABCD.Branch("BWeighted_SLsbin1", &BSLsbin1CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin1", &C3bSLsbin1CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin1", &D3bSLsbin1CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin1", &C2bSLsbin1CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin1", &D2bSLsbin1CountWeighted);
  ABCD.Branch("A_SLsbin2", &ASLsbin2Count);
  ABCD.Branch("B_SLsbin2", &BSLsbin2Count);
  ABCD.Branch("C3b_SLsbin2", &C3bSLsbin2Count);
  ABCD.Branch("D3b_SLsbin2", &D3bSLsbin2Count);
  ABCD.Branch("C2b_SLsbin2", &C2bSLsbin2Count);
  ABCD.Branch("D2b_SLsbin2", &D2bSLsbin2Count);
  ABCD.Branch("AWeighted_SLsbin2", &ASLsbin2CountWeighted);
  ABCD.Branch("BWeighted_SLsbin2", &BSLsbin2CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin2", &C3bSLsbin2CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin2", &D3bSLsbin2CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin2", &C2bSLsbin2CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin2", &D2bSLsbin2CountWeighted);
  ABCD.Branch("A_SLsbin3", &ASLsbin3Count);
  ABCD.Branch("B_SLsbin3", &BSLsbin3Count);
  ABCD.Branch("C3b_SLsbin3", &C3bSLsbin3Count);
  ABCD.Branch("D3b_SLsbin3", &D3bSLsbin3Count);
  ABCD.Branch("C2b_SLsbin3", &C2bSLsbin3Count);
  ABCD.Branch("D2b_SLsbin3", &D2bSLsbin3Count);
  ABCD.Branch("AWeighted_SLsbin3", &ASLsbin3CountWeighted);
  ABCD.Branch("BWeighted_SLsbin3", &BSLsbin3CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin3", &C3bSLsbin3CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin3", &D3bSLsbin3CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin3", &C2bSLsbin3CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin3", &D2bSLsbin3CountWeighted);
  ABCD.Branch("A_SLsbin4", &ASLsbin4Count);
  ABCD.Branch("B_SLsbin4", &BSLsbin4Count);
  ABCD.Branch("C3b_SLsbin4", &C3bSLsbin4Count);
  ABCD.Branch("D3b_SLsbin4", &D3bSLsbin4Count);
  ABCD.Branch("C2b_SLsbin4", &C2bSLsbin4Count);
  ABCD.Branch("D2b_SLsbin4", &D2bSLsbin4Count);
  ABCD.Branch("AWeighted_SLsbin4", &ASLsbin4CountWeighted);
  ABCD.Branch("BWeighted_SLsbin4", &BSLsbin4CountWeighted);
  ABCD.Branch("C3bWeighted_SLsbin4", &C3bSLsbin4CountWeighted);
  ABCD.Branch("D3bWeighted_SLsbin4", &D3bSLsbin4CountWeighted);
  ABCD.Branch("C2bWeighted_SLsbin4", &C2bSLsbin4CountWeighted);
  ABCD.Branch("D2bWeighted_SLsbin4", &D2bSLsbin4CountWeighted);
  ABCD.Branch("A_DRInvsbin1", &ADRInvsbin1Count);
  ABCD.Branch("B_DRInvsbin1", &BDRInvsbin1Count);
  ABCD.Branch("C3b_DRInvsbin1", &C3bDRInvsbin1Count);
  ABCD.Branch("D3b_DRInvsbin1", &D3bDRInvsbin1Count);
  ABCD.Branch("C2b_DRInvsbin1", &C2bDRInvsbin1Count);
  ABCD.Branch("D2b_DRInvsbin1", &D2bDRInvsbin1Count);
  ABCD.Branch("AWeighted_DRInvsbin1", &ADRInvsbin1CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin1", &BDRInvsbin1CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin1", &C3bDRInvsbin1CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin1", &D3bDRInvsbin1CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin1", &C2bDRInvsbin1CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin1", &D2bDRInvsbin1CountWeighted);
  ABCD.Branch("A_DRInvsbin2", &ADRInvsbin2Count);
  ABCD.Branch("B_DRInvsbin2", &BDRInvsbin2Count);
  ABCD.Branch("C3b_DRInvsbin2", &C3bDRInvsbin2Count);
  ABCD.Branch("D3b_DRInvsbin2", &D3bDRInvsbin2Count);
  ABCD.Branch("C2b_DRInvsbin2", &C2bDRInvsbin2Count);
  ABCD.Branch("D2b_DRInvsbin2", &D2bDRInvsbin2Count);
  ABCD.Branch("AWeighted_DRInvsbin2", &ADRInvsbin2CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin2", &BDRInvsbin2CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin2", &C3bDRInvsbin2CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin2", &D3bDRInvsbin2CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin2", &C2bDRInvsbin2CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin2", &D2bDRInvsbin2CountWeighted);
  ABCD.Branch("A_DRInvsbin3", &ADRInvsbin3Count);
  ABCD.Branch("B_DRInvsbin3", &BDRInvsbin3Count);
  ABCD.Branch("C3b_DRInvsbin3", &C3bDRInvsbin3Count);
  ABCD.Branch("D3b_DRInvsbin3", &D3bDRInvsbin3Count);
  ABCD.Branch("C2b_DRInvsbin3", &C2bDRInvsbin3Count);
  ABCD.Branch("D2b_DRInvsbin3", &D2bDRInvsbin3Count);
  ABCD.Branch("AWeighted_DRInvsbin3", &ADRInvsbin3CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin3", &BDRInvsbin3CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin3", &C3bDRInvsbin3CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin3", &D3bDRInvsbin3CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin3", &C2bDRInvsbin3CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin3", &D2bDRInvsbin3CountWeighted);
  ABCD.Branch("A_DRInvsbin4", &ADRInvsbin4Count);
  ABCD.Branch("B_DRInvsbin4", &BDRInvsbin4Count);
  ABCD.Branch("C3b_DRInvsbin4", &C3bDRInvsbin4Count);
  ABCD.Branch("D3b_DRInvsbin4", &D3bDRInvsbin4Count);
  ABCD.Branch("C2b_DRInvsbin4", &C2bDRInvsbin4Count);
  ABCD.Branch("D2b_DRInvsbin4", &D2bDRInvsbin4Count);
  ABCD.Branch("AWeighted_DRInvsbin4", &ADRInvsbin4CountWeighted);
  ABCD.Branch("BWeighted_DRInvsbin4", &BDRInvsbin4CountWeighted);
  ABCD.Branch("C3bWeighted_DRInvsbin4", &C3bDRInvsbin4CountWeighted);
  ABCD.Branch("D3bWeighted_DRInvsbin4", &D3bDRInvsbin4CountWeighted);
  ABCD.Branch("C2bWeighted_DRInvsbin4", &C2bDRInvsbin4CountWeighted);
  ABCD.Branch("D2bWeighted_DRInvsbin4", &D2bDRInvsbin4CountWeighted);
  ABCD.Branch("A_sbin1", &Asbin1Count);
  ABCD.Branch("B_sbin1", &Bsbin1Count);
  ABCD.Branch("C3b_sbin1", &C3bsbin1Count);
  ABCD.Branch("D3b_sbin1", &D3bsbin1Count);
  ABCD.Branch("C2b_sbin1", &C2bsbin1Count);
  ABCD.Branch("D2b_sbin1", &D2bsbin1Count);
  ABCD.Branch("AWeighted_sbin1", &Asbin1CountWeighted);
  ABCD.Branch("BWeighted_sbin1", &Bsbin1CountWeighted);
  ABCD.Branch("C3bWeighted_sbin1", &C3bsbin1CountWeighted);
  ABCD.Branch("D3bWeighted_sbin1", &D3bsbin1CountWeighted);
  ABCD.Branch("C2bWeighted_sbin1", &C2bsbin1CountWeighted);
  ABCD.Branch("D2bWeighted_sbin1", &D2bsbin1CountWeighted);
  ABCD.Branch("A_sbin2", &Asbin2Count);
  ABCD.Branch("B_sbin2", &Bsbin2Count);
  ABCD.Branch("C3b_sbin2", &C3bsbin2Count);
  ABCD.Branch("D3b_sbin2", &D3bsbin2Count);
  ABCD.Branch("C2b_sbin2", &C2bsbin2Count);
  ABCD.Branch("D2b_sbin2", &D2bsbin2Count);
  ABCD.Branch("AWeighted_sbin2", &Asbin2CountWeighted);
  ABCD.Branch("BWeighted_sbin2", &Bsbin2CountWeighted);
  ABCD.Branch("C3bWeighted_sbin2", &C3bsbin2CountWeighted);
  ABCD.Branch("D3bWeighted_sbin2", &D3bsbin2CountWeighted);
  ABCD.Branch("C2bWeighted_sbin2", &C2bsbin2CountWeighted);
  ABCD.Branch("D2bWeighted_sbin2", &D2bsbin2CountWeighted);
  ABCD.Branch("A_sbin3", &Asbin3Count);
  ABCD.Branch("B_sbin3", &Bsbin3Count);
  ABCD.Branch("C3b_sbin3", &C3bsbin3Count);
  ABCD.Branch("D3b_sbin3", &D3bsbin3Count);
  ABCD.Branch("C2b_sbin3", &C2bsbin3Count);
  ABCD.Branch("D2b_sbin3", &D2bsbin3Count);
  ABCD.Branch("AWeighted_sbin3", &Asbin3CountWeighted);
  ABCD.Branch("BWeighted_sbin3", &Bsbin3CountWeighted);
  ABCD.Branch("C3bWeighted_sbin3", &C3bsbin3CountWeighted);
  ABCD.Branch("D3bWeighted_sbin3", &D3bsbin3CountWeighted);
  ABCD.Branch("C2bWeighted_sbin3", &C2bsbin3CountWeighted);
  ABCD.Branch("D2bWeighted_sbin3", &D2bsbin3CountWeighted);
  ABCD.Branch("A_sbin4", &Asbin4Count);
  ABCD.Branch("B_sbin4", &Bsbin4Count);
  ABCD.Branch("C3b_sbin4", &C3bsbin4Count);
  ABCD.Branch("D3b_sbin4", &D3bsbin4Count);
  ABCD.Branch("C2b_sbin4", &C2bsbin4Count);
  ABCD.Branch("D2b_sbin4", &D2bsbin4Count);
  ABCD.Branch("AWeighted_sbin4", &Asbin4CountWeighted);
  ABCD.Branch("BWeighted_sbin4", &Bsbin4CountWeighted);
  ABCD.Branch("C3bWeighted_sbin4", &C3bsbin4CountWeighted);
  ABCD.Branch("D3bWeighted_sbin4", &D3bsbin4CountWeighted);
  ABCD.Branch("C2bWeighted_sbin4", &C2bsbin4CountWeighted);
  ABCD.Branch("D2bWeighted_sbin4", &D2bsbin4CountWeighted);
  ABCD.Branch("kappa3b", &kappa3b);
  ABCD.Branch("kappa2b", &kappa2b);
  ABCD.Branch("kappatot", &kappatot);
  ABCD.Branch("kappaUncert3b", &kappaUncert3b);
  ABCD.Branch("kappaUncert2b", &kappaUncert2b);
  ABCD.Branch("kappaUncerttot", &kappaUncerttot);
  ABCD.Fill();

  TTree cutFlow("cutFlow", "cutFlow");
  cutFlow.Branch("startCount", &startCount);
  cutFlow.Branch("PVCount", &PVCount);
  cutFlow.Branch("jet2PtCount", &jet2PtCount);
  cutFlow.Branch("CSVTCount", &CSVTCount);
  cutFlow.Branch("METSig30Count", &METSig30Count);
  cutFlow.Branch("METCleaningCount", &METCleaningCount);
  cutFlow.Branch("TriggerCount", &TriggerCount);
  cutFlow.Branch("numJetsCounts", &numJetsCount);
  cutFlow.Branch("minDeltaPhiCount", &minDeltaPhiCount);
  cutFlow.Branch("leptonVetoCount", &leptonVetoCount);
  cutFlow.Branch("isoTrackVetoCount", &isoTrackVetoCount);
  cutFlow.Branch("bTagCount", &bTagCount);
  cutFlow.Branch("higgsCount", &higgsCount);
  cutFlow.Branch("DRCount", &DRCount);
  cutFlow.Branch("METSig50Count", &METSig50Count);
  cutFlow.Branch("METSig100Count", &METSig100Count);
  cutFlow.Branch("METSig150Count", &METSig150Count);
  cutFlow.Branch("startCountWeighted", &startCountWeighted);
  cutFlow.Branch("PVCountWeighted", &PVCountWeighted);
  cutFlow.Branch("jet2PtCountWeighted", &jet2PtCountWeighted);
  cutFlow.Branch("CSVTCountWeighted", &CSVTCountWeighted);
  cutFlow.Branch("METSig30CountWeighted", &METSig30CountWeighted);
  cutFlow.Branch("METCleaningCountWeighted", &METCleaningCountWeighted);
  cutFlow.Branch("TriggerCountWeighted", &TriggerCountWeighted);
  cutFlow.Branch("numJetsCountWeighteds", &numJetsCountWeighted);
  cutFlow.Branch("minDeltaPhiCountWeighted", &minDeltaPhiCountWeighted);
  cutFlow.Branch("leptonVetoCountWeighted", &leptonVetoCountWeighted);
  cutFlow.Branch("isoTrackVetoCountWeighted", &isoTrackVetoCountWeighted);
  cutFlow.Branch("bTagCountWeighted", &bTagCountWeighted);
  cutFlow.Branch("higgsCountWeighted", &higgsCountWeighted);
  cutFlow.Branch("DRCountWeighted", &DRCountWeighted);
  cutFlow.Branch("METSig50CountWeighted", &METSig50CountWeighted);
  cutFlow.Branch("METSig100CountWeighted", &METSig100CountWeighted);
  cutFlow.Branch("METSig150CountWeighted", &METSig150CountWeighted);
  cutFlow.Fill();

  std::cout << "startCount " << startCount << std::endl;
  std::cout << "PVCount " << PVCount << std::endl;
  std::cout << "METCleaningCount " << METCleaningCount << std::endl;
  std::cout << "TriggerCount " << TriggerCount << std::endl;
  std::cout << "numJetsCounts " << numJetsCount << std::endl;
  std::cout << "CSVTCount " << CSVTCount << std::endl;
  std::cout << "jet2PtCount " << jet2PtCount << std::endl;
  std::cout << "minDeltaPhiCount " << minDeltaPhiCount << std::endl;
  std::cout << "leptonVetoCount " << leptonVetoCount << std::endl;
  std::cout << "isoTrackVetoCount " << isoTrackVetoCount << std::endl;
  std::cout << "bTagCount " << bTagCount << std::endl;
  std::cout << "higgsCount " << higgsCount << std::endl;
  std::cout << "DRCount " << DRCount << std::endl;
  std::cout << "METSig50Count " << METSig50Count << std::endl;
  std::cout << "METSig100Count " << METSig100Count << std::endl;
  std::cout << "METSig150Count " << METSig150Count << std::endl;
  std::cout << "startCountWeighted " << startCountWeighted << std::endl;
  std::cout << "PVCountWeighted " << PVCountWeighted << std::endl;
  std::cout << "METCleaningCountWeighted " << METCleaningCountWeighted << std::endl;
  std::cout << "TriggerCountWeighted " << TriggerCountWeighted << std::endl;
  std::cout << "numJetsCountWeighteds " << numJetsCountWeighted << std::endl;
  std::cout << "CSVTCountWeighted " << CSVTCountWeighted << std::endl;
  std::cout << "jet2PtCountWeighted " << jet2PtCountWeighted << std::endl;
  std::cout << "minDeltaPhiCountWeighted " << minDeltaPhiCountWeighted << std::endl;
  std::cout << "leptonVetoCountWeighted " << leptonVetoCountWeighted << std::endl;
  std::cout << "isoTrackVetoCountWeighted " << isoTrackVetoCountWeighted << std::endl;
  std::cout << "bTagCountWeighted " << bTagCountWeighted << std::endl;
  std::cout << "higgsCountWeighted " << higgsCountWeighted << std::endl;
  std::cout << "DRCountWeighted " << DRCountWeighted << std::endl;
  std::cout << "METSig50CountWeighted " << METSig50CountWeighted << std::endl;
  std::cout << "METSig100CountWeighted " << METSig100CountWeighted << std::endl;
  std::cout << "METSig150CountWeighted " << METSig150CountWeighted << std::endl;

  file.Write();
  file.Close();
}

bool EventHandler::Passes2CSVTCut() const{
  return GetNumCSVTJets()>=2;
}

bool EventHandler::PassesMETSig30Cut() const{
  return pfmets_fullSignif>30.0;
}

bool EventHandler::PassesLeptonVetoCut() const{
  return GetNumVetoElectrons()==0 && GetNumVetoMuons()==0 && GetNumVetoTaus()==0;
}
bool EventHandler::PassesSingleLeptonCut() const{
  return GetNumVetoElectrons()+GetNumVetoMuons()==1 && GetNumVetoTaus()==0;
}

bool EventHandler::PassesIsoTrackVetoCut() const{
  return NewGetNumIsoTracks(10.0)==0;
}

bool EventHandler::PassesBTaggingCut() const{
  return GetNumCSVTJets()>=2 && GetNumCSVMJets()>=3 && GetNumCSVLJets()>=4;
}

int EventHandler::GetNumGoodJets() const{
  int numGoodJets(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i)) ++numGoodJets;
  }
  return numGoodJets;
}

int EventHandler::GetNumCSVTJets() const{
  int numPassing(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i) && jets_AK5PF_btag_secVertexCombined->at(i)>CSVTCut){
      ++numPassing;
    }
  }
  return numPassing;
}

int EventHandler::GetNumCSVMJets() const{
  int numPassing(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i) && jets_AK5PF_btag_secVertexCombined->at(i)>CSVMCut){
      ++numPassing;
    }
  }
  return numPassing;
}

int EventHandler::GetNumCSVLJets() const{
  int numPassing(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i) && jets_AK5PF_btag_secVertexCombined->at(i)>CSVLCut){
      ++numPassing;
    }
  }
  return numPassing;
}

double EventHandler::GetMinDeltaPhiMET(const unsigned int maxjets) const{
  std::vector<std::pair<double, double> > jets(0);
  for(unsigned int i(0); i<jets_AK5PF_phi->size(); ++i){
    if(isGoodJet(i, false, 20.0, 5.0, false)){
      jets.push_back(std::make_pair(jets_AK5PF_pt->at(i),jets_AK5PF_phi->at(i)));
    }
  }

  std::sort(jets.begin(), jets.end(), std::greater<std::pair<double, double> >());

  double mindp(DBL_MAX);
  for(unsigned int i(0); i<jets.size() && i<maxjets; ++i){
    const double thisdp(fabs((Math::GetAbsDeltaPhi(jets.at(i).second, pfTypeImets_phi->at(0)))));
    if(thisdp<mindp){
      mindp=thisdp;
    }
  }
  
  return mindp;
}

bool EventHandler::PassesMinDeltaPhiCut() const{
  if(pfmets_fullSignif<50.0){ 
    return GetMinDeltaPhiMET(UINT_MAX)>0.5;
  }else{
    return GetMinDeltaPhiMET(UINT_MAX)>0.3;
  }
}

bool EventHandler::isProblemJet(const unsigned int ijet) const{
  return jets_AK5PF_pt->at(ijet)>50.0
    && fabs(jets_AK5PF_eta->at(ijet))>0.9
    && fabs(jets_AK5PF_eta->at(ijet))<1.9
    && jets_AK5PF_chg_Mult->at(ijet)-jets_AK5PF_neutral_Mult->at(ijet)>=40;
}

bool EventHandler::isGoodJet(const unsigned int ijet, const bool jetid, const double ptThresh, const double etaThresh, const bool doBeta) const{
  if(jets_AK5PF_pt->at(ijet)<ptThresh || fabs(jets_AK5PF_eta->at(ijet))>etaThresh) return false;
  if( jetid && !jetPassLooseID(ijet) ) return false;
  if(GetcfAVersion()<69||sampleName.find("SMS-TChiHH")!=std::string::npos) return true;
  if(!betaUpToDate) GetBeta();
  if(beta.at(ijet)<0.2 && doBeta) return false;
  return true;
}

bool EventHandler::jetPassLooseID(const unsigned int ijet) const{
  //want the uncorrected energy
  const double jetenergy = jets_AK5PF_energy->at(ijet) * jets_AK5PF_corrFactorRaw->at(ijet);
  const int numConst = static_cast<int>(jets_AK5PF_mu_Mult->at(ijet)+jets_AK5PF_neutral_Mult->at(ijet)+jets_AK5PF_chg_Mult->at(ijet)); //stealing from Keith
  
  if (jetenergy>0.0) {
    if (jets_AK5PF_neutralHadE->at(ijet) /jetenergy <= 0.99
	&& jets_AK5PF_neutralEmE->at(ijet) / jetenergy <= 0.99
	&& numConst >= 2
	&& ( fabs(jets_AK5PF_eta->at(ijet))>=2.4
	     || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chgHadE->at(ijet)/jetenergy>0))
	&& ( fabs(jets_AK5PF_eta->at(ijet))>=2.4
	     || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chgEmE->at(ijet)/jetenergy<0.99))
	&& ( fabs(jets_AK5PF_eta->at(ijet))>=2.4
	     || (fabs(jets_AK5PF_eta->at(ijet))<2.4 && jets_AK5PF_chg_Mult->at(ijet)>0))){
      return true;
    }
  }
  return false;
}

bool EventHandler::isVetoElectron(const unsigned int k) const{
  //if(k>pf_else_pt->size()) return false;
  if (fabs(pf_els_scEta->at(k)) >= 2.5 ) return false;
  if (pf_els_pt->at(k) < 10) return false;
  if(pf_els_isEB->at(k)){
    if ( fabs(pf_els_dEtaIn->at(k)) > 0.007)  return false;
    if ( fabs(pf_els_dPhiIn->at(k)) > 0.8)  return false;
    if (pf_els_sigmaIEtaIEta->at(k) > 0.01) return false;
    if (pf_els_hadOverEm->at(k) > 0.15) return false;
  }else if(pf_els_isEE->at(k)){
    if ( fabs(pf_els_dEtaIn->at(k)) > 0.01)  return false;
    if ( fabs(pf_els_dPhiIn->at(k)) > 0.7)  return false;
    if (pf_els_sigmaIEtaIEta->at(k) > 0.03) return false;
  }else{
    fprintf(stderr, "Warning: Electron is not in barrel or endcap.\n");
    return false;
  }
  const double beamx(beamSpot_x->at(0)), beamy(beamSpot_y->at(0)); 
  const double d0(pf_els_d0dum->at(k)-beamx*sin(pf_els_tk_phi->at(k))+beamy*cos(pf_els_tk_phi->at(k)));
  if ( fabs(d0) >= 0.04 ) return false;
  if ( fabs(pf_els_vz->at(k) - pv_z->at(0) ) >= 0.2 ) return false;

  if(GetElectronRelIso(k)>=0.15) return false;
  return true;
}

double EventHandler::GetElectronRelIso(const unsigned int k) const{
  const double rho(rho_kt6PFJetsForIsolation2012);
  // get effective area from delR=0.3 2011 data table for neutral+gamma based on supercluster eta pf_els_scEta->at(k)
  double AE(0.10); 
  const double abseta(fabs(pf_els_scEta->at(k)));
  if      ( abseta < 1.0 ) AE = 0.13;
  else if ( abseta >=1.0 && abseta <1.479) AE = 0.14;
  else if ( abseta >=1.479&&abseta <2.0)   AE = 0.07;
  else if ( abseta >=2.0 && abseta <2.2) AE = 0.09;
  else if ( abseta >=2.2 && abseta <2.3) AE = 0.11;
  else if ( abseta >=2.3 && abseta <2.4) AE = 0.11;
  else if ( abseta >=2.4 ) AE = 0.14;

  const double eleIso(pf_els_PFphotonIsoR03->at(k) + pf_els_PFneutralHadronIsoR03->at(k) - rho*AE);
  return ( pf_els_PFchargedHadronIsoR03->at(k) + ( eleIso > 0 ? eleIso : 0.0 ) )/pf_els_pt->at(k);
}

bool EventHandler::isVetoMuon(const unsigned int k) const{
  if (fabs(pf_mus_eta->at(k)) >= 2.4 ) return false;
  if (pf_mus_pt->at(k) < 10) return false;
  if ( !pf_mus_id_GlobalMuonPromptTight->at(k)) return false;
  // GlobalMuonPromptTight includes: isGlobal, globalTrack()->normalizedChi2() < 10, numberOfValidMuonHits() > 0
  if ( pf_mus_numberOfMatchedStations->at(k) <= 1 ) return false;
  const double beamx (beamSpot_x->at(0)), beamy(beamSpot_y->at(0));   
  const double d0 = pf_mus_tk_d0dum->at(k)-beamx*sin(pf_mus_tk_phi->at(k))+beamy*cos(pf_mus_tk_phi->at(k));
  const double pf_mus_vz = pf_mus_tk_vz->at(k);
  const double pf_mus_dz_vtx = fabs(pf_mus_vz-pv_z->at(0));
  if (fabs(d0)>=0.2 || pf_mus_dz_vtx>=0.5) return false;
  if ( !pf_mus_tk_numvalPixelhits->at(k)) return false;
  if ( pf_mus_tk_LayersWithMeasurement->at(k) <= 5 ) return false;
  
  double isoNeutral(pf_mus_pfIsolationR04_sumNeutralHadronEt->at(k) + pf_mus_pfIsolationR04_sumPhotonEt->at(k) - 0.5*pf_mus_pfIsolationR04_sumPUPt->at(k));
  if(isoNeutral<0.0) isoNeutral=0.0;
  const double pf_mus_rel_iso((pf_mus_pfIsolationR04_sumChargedHadronPt->at(k) + isoNeutral) / pf_mus_pt->at(k));
  if (pf_mus_rel_iso > 0.2) return false;
  return true;
}

bool EventHandler::isVetoTau(const unsigned int k) const{
  if (taus_pt->at(k)<20) return false; // Updated 6/24
  if (fabs(taus_eta->at(k)) > 2.4) return false;
  if (taus_byLooseIsolationDeltaBetaCorr->at(k) <= 0) return false;
  return true;
}

int EventHandler::GetNumVetoElectrons() const{
  int count(0);
  for(unsigned int i(0); i<pf_els_pt->size(); ++i){
    if(isVetoElectron(i)) ++count;
  }
  return count;
}

int EventHandler::GetNumVetoMuons() const{
  int count(0);
  for(unsigned int i(0); i<pf_mus_pt->size(); ++i){
    if(isVetoMuon(i)) ++count;
  }
  return count;
}

int EventHandler::GetNumVetoTaus() const{
  int count(0);
  for(unsigned int i(0); i<taus_pt->size(); ++i){
    if(isVetoTau(i)) ++count;
  }
  return count;
}

bool EventHandler::isIsoTrack(const unsigned int itracks, const double ptThresh) const{
  if(!isQualityTrack(itracks)) return false;
  if (fabs(tracks_vz->at(itracks) - pv_z->at(0) ) >= 0.05) return false;
  if(tracks_pt->at(itracks)<ptThresh) return false;
  double isosum=0;
  for (unsigned int jtracks=0; jtracks<tracks_chi2->size(); jtracks++) {
    if (itracks==jtracks) continue;  //don't count yourself
    if (!isQualityTrack(jtracks)) continue;
    if(Math::GetDeltaR(tracks_phi->at(itracks),tracks_eta->at(itracks),tracks_phi->at(jtracks),tracks_eta->at(jtracks))>0.3) continue;
    //cut on dz of this track
    if ( fabs( tracks_vz->at(jtracks) - pv_z->at(0)) >= 0.05) continue;
    isosum += tracks_pt->at(jtracks);
  }
  if ( isosum / tracks_pt->at(itracks) > 0.05) return false;
  return true;
}

bool EventHandler::isQualityTrack(const unsigned int k) const{
  if (fabs(tracks_eta->at(k))>2.4 ) return false;
  if (!(tracks_highPurity->at(k)>0)) return false;
  return true;  
}

int EventHandler::GetNumIsoTracks(const double ptThresh) const{
  int count(0);
  for(unsigned int i(0); i<tracks_pt->size(); ++i){
    if(isIsoTrack(i,ptThresh)) ++count;
  }
  return count;
}

int EventHandler::NewGetNumIsoTracks(const double ptThresh) const{
  int nisotracks=0;
  if ( GetcfAVersion() < 71 ) return nisotracks;
  for ( unsigned int itrack = 0 ; itrack < isotk_pt->size() ; ++itrack) {
    if ( (isotk_pt->at(itrack) >= ptThresh) &&
         (isotk_iso->at(itrack) /isotk_pt->at(itrack) < 0.1 ) &&
         ( fabs(isotk_dzpv->at(itrack)) <0.1) && //important to apply this; was not applied at ntuple creation
         ( fabs(isotk_eta->at(itrack)) < 2.4 ) //this is more of a sanity check
	 ){
      ++nisotracks;
    }
  }
  return nisotracks;
}

double EventHandler::GetPUWeight(reweight::LumiReWeighting &lumiWeights) const{
  double npv(-1.0);
  for(unsigned int i(0); i<PU_bunchCrossing->size(); ++i){
    if(PU_bunchCrossing->at(i)==0){
      npv = PU_TrueNumInteractions->at(i);
    }
  }
  return lumiWeights.weight(npv);
}

double EventHandler::GetSbinWeight() const{
  if(sampleName.find("Run2012")!=std::string::npos){
    return 1.0;
  }else if(pfmets_fullSignif>30.0 && pfmets_fullSignif<50.0){
    return 0.804;
  }else if(pfmets_fullSignif>50.0 && pfmets_fullSignif<100.0){
    return 0.897;
  }else if(pfmets_fullSignif>100.0){
    return 0.944;
  }else{
    return 1.0;
  }
}

double EventHandler::GetTopPtWeight() const{
  double topPt(-1.0);

  //only for use with ttbar
  //(2 string comparisons for every event is not so good for performance, but oh well)
  if (sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos) {

    //code from Darren converted to cfA
    for ( unsigned int i=0; i< mc_doc_id->size(); i++ ) {
      // look for the *top* (not antitop) pT
      if ( mc_doc_id->at(i)==6 ) { topPt = mc_doc_pt->at(i); break; }
    }

    if (topPt<0) return 1;
    // note -- we should remove this, now that it is demonstrated that powheg and mc@nlo behave the same as MG
    if ( sampleName.find("madgraph")!=std::string::npos) { //only return weight for madgraph

      const  double p0 = 1.18246e+00;
      const  double p1 = 4.63312e+02;
      const  double p2 = 2.10061e-06;

      double x=topPt;
      if ( x>p1 ) x = p1; //use flat factor above 463 GeV
      double result = p0 + p2 * x * ( x - 2 * p1 );
      return double(result);
    }
  }

  return 1;
}

void EventHandler::Skim(const std::string &skimFileName,
			const int chargino_mass,
			const int LSP_mass){
  GetTotalEntries();

  TFile skimFile(skimFileName.c_str(), "RECREATE");
  TDirectory *cfA_dir(skimFile.mkdir("configurableAnalysis", "configurableAnalysis"));
  cfA_dir->cd();
  TTree *skimTreeA(chainA.CloneTree(0));
  chainA.CopyAddresses(skimTreeA);
  TTree *skimTreeB(chainB.CloneTree(0));
  chainB.CopyAddresses(skimTreeB);
  skimFile.cd();

  std::set<EventNumber> eventList;
  unsigned int startCount(0), PVCount(0), METCleaningCount(0), JSONCount(0), TriggerCount(0), numJetsCount(0), CSVTCount(0), jet2PtCount(0), minDeltaPhiCount(0), leptonVetoCount(0), isoTrackVetoCount(0), bTagCount(0), higgsCount(0), DRCount(0), METSig30Count(0), METSig50Count(0), METSig100Count(0), METSig150Count(0);
  double startCountWeighted(0), PVCountWeighted(0), METCleaningCountWeighted(0), JSONCountWeighted(0), TriggerCountWeighted(0), numJetsCountWeighted(0), CSVTCountWeighted(0), jet2PtCountWeighted(0), minDeltaPhiCountWeighted(0), leptonVetoCountWeighted(0), isoTrackVetoCountWeighted(0), bTagCountWeighted(0), higgsCountWeighted(0), DRCountWeighted(0), METSig30CountWeighted(0), METSig50CountWeighted(0), METSig100CountWeighted(0), METSig150CountWeighted(0);

  std::vector<float> dataDist(pu::RunsThrough203002, pu::RunsThrough203002+60);
  std::vector<float> MCDist(pu::Summer2012, pu::Summer2012+60);//QQQ this needs to change later for general pileup scenario
  reweight::LumiReWeighting lumiWeights(MCDist, dataDist);

  Timer timer(GetTotalEntries());
  timer.Start();
  for(int i(0); i<GetTotalEntries(); ++i){
    if(i%1000==0 && i!=0){
      timer.PrintRemainingTime();
    }
    timer.Iterate();

    GetEntry(i);
    std::pair<std::set<EventNumber>::iterator, bool> returnVal(eventList.insert(EventNumber(run, event, lumiblock)));
    if(!returnVal.second) continue;
    const bool isttbar(sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos);
    const double localWeight(GetPUWeight(lumiWeights)*scaleFactor*(isttbar?GetTopPtWeight():1.0)*GetSbinWeight());
    
    // Select mass points
    if(!PassesTChiZHMassCut()) continue;
    ++startCount;
    startCountWeighted+=localWeight;

    if(Passes2CSVTCut()
       && PassesPVCut()
       && PassesJet2PtCut()
       && PassesMETSig30Cut()
       && PassesTChiHHMassCut(chargino_mass, LSP_mass)){
      skimTreeA->Fill();
      skimTreeB->Fill();
    }

    if(!PassesPVCut()) continue;
    ++PVCount;
    PVCountWeighted+=localWeight;
    if(!PassesJet2PtCut()) continue;
    ++jet2PtCount;
    jet2PtCountWeighted+=localWeight;
    if(!Passes2CSVTCut()) continue;
    ++CSVTCount;
    CSVTCountWeighted+=localWeight;
    if(!PassesMETSig30Cut()) continue;
    ++METSig30Count;
    METSig30CountWeighted+=localWeight;
    if(!PassesMETCleaningCut()) continue;
    ++METCleaningCount;
    METCleaningCountWeighted+=localWeight;
    if(!PassesJSONCut()) continue;
    ++JSONCount;
    JSONCountWeighted+=localWeight;
    if(!PassesTriggerCut()) continue;
    ++TriggerCount;
    TriggerCountWeighted+=localWeight;
    if(!PassesNumJetsCut()) continue;
    ++numJetsCount;
    numJetsCountWeighted+=localWeight;
    if(!PassesMinDeltaPhiCut()) continue;
    ++minDeltaPhiCount;
    minDeltaPhiCountWeighted+=localWeight;
    if(!PassesLeptonVetoCut()) continue;
    ++leptonVetoCount;
    leptonVetoCountWeighted+=localWeight;
    if(!PassesIsoTrackVetoCut()) continue;
    ++isoTrackVetoCount;
    isoTrackVetoCountWeighted+=localWeight;
    if(!PassesBTaggingCut()) continue;
    ++bTagCount;
    bTagCountWeighted+=localWeight;
    if(!PassesHiggsMassCut()) continue;
    ++higgsCount;
    higgsCountWeighted+=localWeight;
    if(!PassesDRCut()) continue;
    ++DRCount;
    DRCountWeighted+=localWeight;
    if(!PassesMETSig50Cut()) continue;
    ++METSig50Count;
    METSig50CountWeighted+=localWeight;
    if(!PassesMETSig100Cut()) continue;
    ++METSig100Count;
    METSig100CountWeighted+=localWeight;
    if(!PassesMETSig150Cut()) continue;
    ++METSig150Count;
    METSig150CountWeighted+=localWeight;
  }

  TTree cutFlow("cutFlow", "cutFlow");
  cutFlow.Branch("startCount", &startCount);
  cutFlow.Branch("CSVTCount", &CSVTCount);
  cutFlow.Branch("PVCount", &PVCount);
  cutFlow.Branch("TriggerCount", &TriggerCount);
  cutFlow.Branch("numJetsCounts", &numJetsCount);
  cutFlow.Branch("jet2PtCount", &jet2PtCount);
  cutFlow.Branch("JSONCount", &JSONCount);
  cutFlow.Branch("minDeltaPhiCount", &minDeltaPhiCount);
  cutFlow.Branch("METSig30Count", &METSig30Count);
  cutFlow.Branch("METCleaningCount", &METCleaningCount);
  cutFlow.Branch("leptonVetoCount", &leptonVetoCount);
  cutFlow.Branch("isoTrackVetoCount", &isoTrackVetoCount);
  cutFlow.Branch("bTagCount", &bTagCount);
  cutFlow.Branch("higgsCount", &higgsCount);
  cutFlow.Branch("DRCount", &DRCount);
  cutFlow.Branch("METSig50Count", &METSig50Count);
  cutFlow.Branch("METSig100Count", &METSig100Count);
  cutFlow.Branch("METSig150Count", &METSig150Count);
  cutFlow.Branch("startCountWeighted", &startCountWeighted);
  cutFlow.Branch("PVCountWeighted", &PVCountWeighted);
  cutFlow.Branch("JSONCountWeighted", &JSONCountWeighted);
  cutFlow.Branch("METCleaningCountWeighted", &METCleaningCountWeighted);
  cutFlow.Branch("TriggerCountWeighted", &TriggerCountWeighted);
  cutFlow.Branch("numJetsCountWeighteds", &numJetsCountWeighted);
  cutFlow.Branch("CSVTCountWeighted", &CSVTCountWeighted);
  cutFlow.Branch("jet2PtCountWeighted", &jet2PtCountWeighted);
  cutFlow.Branch("minDeltaPhiCountWeighted", &minDeltaPhiCountWeighted);
  cutFlow.Branch("METSig30CountWeighted", &METSig30CountWeighted);
  cutFlow.Branch("leptonVetoCountWeighted", &leptonVetoCountWeighted);
  cutFlow.Branch("isoTrackVetoCountWeighted", &isoTrackVetoCountWeighted);
  cutFlow.Branch("bTagCountWeighted", &bTagCountWeighted);
  cutFlow.Branch("higgsCountWeighted", &higgsCountWeighted);
  cutFlow.Branch("DRCountWeighted", &DRCountWeighted);
  cutFlow.Branch("METSig50CountWeighted", &METSig50CountWeighted);
  cutFlow.Branch("METSig100CountWeighted", &METSig100CountWeighted);
  cutFlow.Branch("METSig150CountWeighted", &METSig150CountWeighted);
  cutFlow.Fill();

  TTree timeInfo("timeInfo", "timeInfo");
  time_t my_timer_t;
  struct tm *timer_s(0);
  time(&my_timer_t);
  timer_s=localtime(&my_timer_t);
  timeInfo.Branch("year", &timer_s->tm_year);
  timeInfo.Branch("month", &timer_s->tm_mon);
  timeInfo.Branch("day", &timer_s->tm_mday);
  timeInfo.Branch("hour", &timer_s->tm_hour);
  timeInfo.Branch("minute", &timer_s->tm_min);
  timeInfo.Branch("second", &timer_s->tm_sec);
  timeInfo.Fill();

  cfA_dir->cd();
  skimTreeA->Write();
  skimTreeB->Write();
  skimFile.cd();
  cutFlow.Write();
  timeInfo.Write();
  skimFile.Close();
}
