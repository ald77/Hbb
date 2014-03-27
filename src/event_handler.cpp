#include "event_handler.hpp"
#include <cassert>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <utility>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <limits>
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

EventHandler::EventHandler(const std::string &fileName, const bool isList, const double scaleFactorIn, const bool fastMode):
  cfA(fileName, isList),
  higgsBJetPairing(std::make_pair(TLorentzVector(0.0,0.0,0.0,0.0),TLorentzVector(0.0,0.0,0.0,0.0)),std::make_pair(TLorentzVector(0.0,0.0,0.0,0.0),TLorentzVector(0.0,0.0,0.0,0.0))),
  sortedBJetCache(0),
  higgsPairingUpToDate(false),
  bJetsUpToDate(false),
  betaUpToDate(false),
  scaleFactor(scaleFactorIn),
  beta(0){
  if (fastMode) { // turn off unnecessary branches
    chainA.SetBranchStatus("els_*",0);
    chainA.SetBranchStatus("triggerobject_*",0);
    chainA.SetBranchStatus("standalone_t*",0);
    chainA.SetBranchStatus("L1trigger_*",0);
    chainA.SetBranchStatus("passprescale*",0);
    chainA.SetBranchStatus("jets_AK5PFclean_*",0);
    chainA.SetBranchStatus("softjetUp_*",0);
    chainA.SetBranchStatus("pdfweights_*",0);
    chainA.SetBranchStatus("photon_*",0);
    chainB.SetBranchStatus("Ntcmets",0);
    chainB.SetBranchStatus("tcmets_*",0);
    chainB.SetBranchStatus("Nphotons",0);
    chainB.SetBranchStatus("photons_*",0);
    chainB.SetBranchStatus("Npf_photons",0);
    chainB.SetBranchStatus("pf_photons_*",0);
    chainB.SetBranchStatus("Nmus",0);
    chainB.SetBranchStatus("mus_*",0);
    chainB.SetBranchStatus("Nels",0);
    chainB.SetBranchStatus("els_*",0);
    chainB.SetBranchStatus("Nmets*",0);
    chainB.SetBranchStatus("mets*",0);
    chainB.SetBranchStatus("mets_AK5_et",1);
    chainB.SetBranchStatus("Njets_AK5PFclean",0);
    chainB.SetBranchStatus("jets_AK5PFclean_*",0);
    chainB.SetBranchStatus("Nmc*",0);
    chainB.SetBranchStatus("mc_*",0);
    chainB.SetBranchStatus("Nmc_doc*",1);
    chainB.SetBranchStatus("mc_doc*",1);
  }
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
  return isGoodVertex(0);
}

bool EventHandler::PassesMETCleaningCut() const{
  for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
    if(isProblemJet(jet)) return false;
  }
  //if(pfTypeImets_et->at(0)>2.0*pfmets_et->at(0)) return false;
  if(pfTypeImets_et->at(0)>2.0*mets_AK5_et->at(0)) return false; // updated 11/14 (JB-F)
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

bool EventHandler::PassesQCDTriggerCut() const{
  for(unsigned int a=0; a<trigger_name->size(); ++a){
    if((trigger_name->at(a).find("HLT_DiCentralPFJet50_PFMET80_v")!=std::string::npos
        || trigger_name->at(a).find("HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v")!=std::string::npos)
       && trigger_prescalevalue->at(a)==1 && trigger_decision->at(a)==1){
      return true;
    }
  }
  return false;
}

bool EventHandler::PassesSpecificTrigger(const std::string trigger) const{ // just check if a specific trigger fired
  for(unsigned int a=0; a<trigger_name->size(); ++a){
    if(trigger_name->at(a).find(trigger)!=std::string::npos
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
  uint_least32_t fail_code(kGood);
  if(!PassesMETSig150Cut()) fail_code |= kMETSig150;
  if(!PassesMETSig100Cut()) fail_code |= kMETSig100;
  if(!PassesMETSig50Cut()) fail_code |= kMETSig50;
  if(!PassesMETSig30Cut()) fail_code |= kMETSig30;
  if(!PassesDRCut()) fail_code |= kDeltaR;
  if(!PassesHiggsAvgMassCut()) fail_code |= kHiggsAvgMass;
  if(!PassesHiggsMassDiffCut()) fail_code |= kHiggsMassDiff;
  if(!PassesBTaggingCut()) fail_code |= k4thBTag;
  if(!(GetNumCSVTJets()>=2 && GetNumCSVMJets()>=3)) fail_code |= k3rdBTag;
  if(!(GetNumCSVTJets()>=2)) fail_code |= k2ndBTag;
  if(!PassesIsoTrackVetoCut()) fail_code |= kIsoTrackVeto;
  if(!PassesLeptonVetoCut()) fail_code |= kLeptonVeto;
  if(!PassesMinDeltaPhiCut()) fail_code |= kMinDeltaPhi;
  if(!PassesNumJetsCut()) fail_code |= kNumJets;
  if(!PassesTriggerCut()) fail_code |= kTrigger;
  if(!PassesJSONCut()) fail_code |= kJSON;
  if(!PassesMETCleaningCut()) fail_code |= kMETCleaning;
  if(!PassesPVCut()) fail_code |= kPV;
  if(!PassesTChiMassCut()) fail_code |= kTChiMassCut;
  return fail_code;
}

bool EventHandler::PassesTriggerPlateauCuts() const{
  if(!Passes2CSVTCut()) return false;
  if(!PassesJet2PtCut()) return false;
  if(!PassesMETSig30Cut()) return false;
  return true;
}

bool EventHandler::PassesBaselineSelection() const{
  if(!PassesJSONCut()) return false;
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
  return true;
}

bool EventHandler::PassesSingleLeptonControlCut() const {
  if(GetNumElectrons(1)+GetNumMuons(1)!=1) return false;
  if(GetNumTaus()>0) return false;
  if(!PassesIsoTrackVetoCut()) return false;
  return true;
}

bool EventHandler::PassesQCDControlCut() const {
  if(!PassesQCDTriggerCut()) return false;
  if(GetNumCSVMJets()>2) return false;
  if(GetMinDeltaPhiMET(UINT_MAX)>0.3) return false;
  return true;
}

bool EventHandler::PassesRegionACut() const{
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(!PassesJSONCut()) return false;
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
  if(useMET && pfTypeImets_et->size()>0) HT+=pfTypeImets_et->at(0);
  for(unsigned int i(0); i<jets_AK5PF_pt->size(); ++i){
    if(isGoodJet(i)) HT+=jets_AK5PF_pt->at(i);
  }
  if(useLeps){
    for(unsigned int i(0); i<pf_els_pt->size(); ++i){
      if(isElectron(i)) HT+=pf_els_pt->at(i);
    }
    for(unsigned int i(0); i<pf_mus_pt->size(); ++i){
      if(isMuon(i)) HT+=pf_mus_pt->at(i);
    }
    for(unsigned int i(0); i<taus_pt->size(); ++i){
      if(isTau(i)) HT+=taus_pt->at(i);
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

double EventHandler::GetMinDR() const{
  GetHiggsBJetPairing();
  if(sortedBJetCache.size()<4){
    return DBL_MAX;
  }else{
    const double dRa(higgsBJetPairing.first.first.DeltaR(higgsBJetPairing.first.second));
    const double dRb(higgsBJetPairing.second.first.DeltaR(higgsBJetPairing.second.second));
    return dRa<dRb?dRa:dRb;
  }
}

std::pair<double, double> EventHandler::GetHiggsMasses() const{
  GetHiggsBJetPairing();
  const double massA((higgsBJetPairing.first.first+higgsBJetPairing.first.second).M());
  const double massB((higgsBJetPairing.second.first+higgsBJetPairing.second.second).M());
  return (massA>massB)?(std::make_pair(massA,massB)):(std::make_pair(massB,massA));
}

double EventHandler::GetMinimaxMbb() const{
  GetHiggsBJetPairing();
  const double mass_1a((higgsBJetPairing.first.first
                        +higgsBJetPairing.first.second).M());
  const double mass_1b((higgsBJetPairing.second.first
                        +higgsBJetPairing.second.second).M());
  const double mass_2a((higgsBJetPairing.first.first
                        +higgsBJetPairing.second.first).M());
  const double mass_2b((higgsBJetPairing.first.second
                        +higgsBJetPairing.second.second).M());
  const double mass_3a((higgsBJetPairing.first.first
                        +higgsBJetPairing.second.second).M());
  const double mass_3b((higgsBJetPairing.first.second
                        +higgsBJetPairing.second.first).M());
  const double mass_1(std::max(mass_1a, mass_1b));
  const double mass_2(std::max(mass_2a, mass_2b));
  const double mass_3(std::max(mass_3a, mass_3b));
  return std::min(mass_1, std::min(mass_2, mass_3));
}

double EventHandler::GetMaximinMbb() const{
  GetHiggsBJetPairing();
  const double mass_1a((higgsBJetPairing.first.first
                        +higgsBJetPairing.first.second).M());
  const double mass_1b((higgsBJetPairing.second.first
                        +higgsBJetPairing.second.second).M());
  const double mass_2a((higgsBJetPairing.first.first
                        +higgsBJetPairing.second.first).M());
  const double mass_2b((higgsBJetPairing.first.second
                        +higgsBJetPairing.second.second).M());
  const double mass_3a((higgsBJetPairing.first.first
                        +higgsBJetPairing.second.second).M());
  const double mass_3b((higgsBJetPairing.first.second
                        +higgsBJetPairing.second.first).M());
  const double mass_1(std::min(mass_1a, mass_1b));
  const double mass_2(std::min(mass_2a, mass_2b));
  const double mass_3(std::min(mass_3a, mass_3b));
  return std::max(mass_1, std::max(mass_2, mass_3));
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
        sortedBJetCache.push_back(BJet(TLorentzVector(jets_AK5PF_px->at(i),jets_AK5PF_py->at(i),jets_AK5PF_pz->at(i),jets_AK5PF_energy->at(i)),jets_AK5PF_btag_secVertexCombined->at(i),i,static_cast<unsigned int>(jets_AK5PF_parton_Id->at(i)),static_cast<unsigned int>(jets_AK5PF_parton_motherId->at(i))));
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

int EventHandler::GetPartonIdBin(float partonId) const{
  if (fabs(partonId)>=1.&&fabs(partonId)<=3.) return 2;
  else if (fabs(partonId)==4.) return 3;
  else if (fabs(partonId)==5.) return 4;
  else if (fabs(partonId)==21.) return 5;
  else return 1;
}

bool EventHandler::PassesTChiMassCut(int mChi, int mLSP) const{
  if (mChi<0||mLSP<0) return true;
  // parse the model_params string in the new signal samples for mass points
  char s_mChi[64], s_mLSP[64];
  sprintf(s_mChi, "chargino%d_", mChi);
  sprintf(s_mLSP, "bino%d_", mLSP);

  if(sampleName.find("SMS-TChiHH")==std::string::npos
     && sampleName.find("SMS_TChiZH")){
    return true;
  }
  if(model_params->find(s_mChi)==std::string::npos) return false;
  if(model_params->find(s_mLSP)==std::string::npos) return false;
  return true;
}

int EventHandler::GetCharginoMass() const{
  return GetMass("chargino");
}

int EventHandler::GetLSPMass() const{
  return GetMass("bino");
}

int EventHandler::GetMass(const std::string& token) const{
  std::string::size_type pos(model_params->find(token));
  if(pos==std::string::npos){
    return -1;
  }else{
    pos+=token.size();
    const std::string::size_type end(model_params->find("_",pos));
    if(end==std::string::npos){
      return -1;
    }else{
      return atoi((model_params->substr(pos, end-pos+1)).c_str());
    }
  }
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

std::vector<std::pair<int, int> > EventHandler::GetBOrigins_new() const{
  if(!bJetsUpToDate) GetSortedBJets();
  std::vector<std::pair<int, int> > x(0);
  for(unsigned int jet(0); jet<sortedBJetCache.size(); ++jet){
    /*
      printf("i: %d, CSV: %3.2f, Id: %d, momId: %d\n", 
      sortedBJetCache[jet].GetIndex(), sortedBJetCache[jet].GetBTag(),
      sortedBJetCache[jet].GetPartonId(), sortedBJetCache[jet].GetMotherId());
    */
    const unsigned int id(sortedBJetCache.at(jet).GetPartonId());
    const unsigned int mom(sortedBJetCache.at(jet).GetMotherId());
    x.push_back(std::pair<int,int>(id, mom));
  }
  return x;
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

bool EventHandler::Passes2CSVTCut() const{
  return GetNumCSVTJets()>=2;
}

bool EventHandler::PassesMETSig30Cut() const{
  return pfmets_fullSignif>30.0;
}

bool EventHandler::PassesLeptonVetoCut() const{
  return GetNumElectrons()==0 && GetNumMuons()==0 && GetNumTaus()==0;
}
bool EventHandler::PassesSingleLeptonCut() const{
  return GetNumElectrons()+GetNumMuons()==1 && GetNumTaus()==0;
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

int EventHandler::GetNumBTaggedJets() const{
  std::vector<double> csvs(0);
  for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
    if(isGoodJet(jet)){
      csvs.push_back(jets_AK5PF_btag_secVertexCombined->at(jet));
    }
  }
  std::sort(csvs.begin(), csvs.end(), std::greater<double>());
  if(csvs.size()<=0 || csvs.at(0)<CSVTCut){
    return 0;
  }else if(csvs.size()<=1 || csvs.at(1)<CSVTCut){
    return 1;
  }else if(csvs.size()<=2 || csvs.at(2)<CSVMCut){
    return 2;
  }else if(csvs.size()<=3 || csvs.at(3)<CSVLCut){
    return 3;
  }else{
    unsigned int jet(3);
    for(; jet<csvs.size() && csvs.at(jet)>CSVLCut; ++jet){
    }
    return jet;
  }
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

bool EventHandler::isElectron(const unsigned int k,
			      const unsigned short level,
			      const bool use_iso) const{
  //N.B.: cut does not have the fabs(1/E-1/p) and conversion rejection cuts from the EGamma POG!!!
  const double dmax(std::numeric_limits<double>::max());
  double pt_cut(10.0); //Not actually part of the EGamma ID
  double eta_cut(0.007), phi_cut(0.8), sigmaietaieta_cut(0.01), h_over_e_cut(0.15), d0_cut(0.04), dz_cut(0.2), iso_cut(0.15);
  switch(level){
  case 1:
    pt_cut=20.0;
    if(pf_els_isEB->at(k)){
      eta_cut=0.007; phi_cut=0.15; sigmaietaieta_cut=0.01; h_over_e_cut=0.12;
      d0_cut=0.02; dz_cut=0.2; iso_cut=0.15;
    }else if(pf_els_isEE->at(k)){
      eta_cut=0.009; phi_cut=0.10; sigmaietaieta_cut=0.03; h_over_e_cut=0.1;
      d0_cut=0.02; dz_cut=0.2; iso_cut=(pf_els_pt->at(k)>20.0?0.15:0.10);
    }
    break;
  case 2:
    pt_cut=20.0;
    if(pf_els_isEB->at(k)){
      eta_cut=0.004; phi_cut=0.06; sigmaietaieta_cut=0.01; h_over_e_cut=0.12;
      d0_cut=0.02; dz_cut=0.1; iso_cut=0.15;
    }else if(pf_els_isEE->at(k)){
      eta_cut=0.007; phi_cut=0.03; sigmaietaieta_cut=0.03; h_over_e_cut=0.1;
      d0_cut=0.02; dz_cut=0.1; iso_cut=(pf_els_pt->at(k)>20.0?0.15:0.10);
    }
    break;
  case 3:
    pt_cut=20.0;
    if(pf_els_isEB->at(k)){
      eta_cut=0.004; phi_cut=0.03; sigmaietaieta_cut=0.01; h_over_e_cut=0.12;
      d0_cut=0.02; dz_cut=0.1; iso_cut=0.10;
    }else if(pf_els_isEE->at(k)){
      eta_cut=0.005; phi_cut=0.02; sigmaietaieta_cut=0.03; h_over_e_cut=0.1;
      d0_cut=0.02; dz_cut=0.1; iso_cut=(pf_els_pt->at(k)>20.0?0.10:0.07);
    }
    break;
  case 0: //intentionally falling through to default "veto" case
  default:
    pt_cut=10.0;
    if(pf_els_isEB->at(k)){
      eta_cut=0.007; phi_cut=0.8; sigmaietaieta_cut=0.01; h_over_e_cut=0.15;
      d0_cut=0.04; dz_cut=0.2; iso_cut=0.15;
    }else if(pf_els_isEE->at(k)){
      eta_cut=0.01; phi_cut=0.7; sigmaietaieta_cut=0.03; h_over_e_cut=dmax;
      d0_cut=0.04; dz_cut=0.2; iso_cut=0.15;
    }
    break;
  }
  //if(k>pf_else_pt->size()) return false;
  if (fabs(pf_els_scEta->at(k)) >= 2.5 ) return false;
  if (pf_els_pt->at(k) < pt_cut) return false;

  if ( fabs(pf_els_dEtaIn->at(k)) > eta_cut)  return false;
  if ( fabs(pf_els_dPhiIn->at(k)) > phi_cut)  return false;
  if (pf_els_sigmaIEtaIEta->at(k) > sigmaietaieta_cut) return false;
  if (pf_els_hadOverEm->at(k) > h_over_e_cut) return false;

  const double beamx(beamSpot_x->at(0)), beamy(beamSpot_y->at(0)); 
  const double d0(pf_els_d0dum->at(k)-beamx*sin(pf_els_tk_phi->at(k))+beamy*cos(pf_els_tk_phi->at(k)));
  if ( fabs(d0) >= d0_cut ) return false;
  if ( fabs(pf_els_vz->at(k) - pv_z->at(0) ) >= dz_cut ) return false;

  if(GetElectronRelIso(k)>=iso_cut && use_iso) return false;
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

bool EventHandler::isMuon(const unsigned int k,
			  const unsigned short level,
			  const bool use_iso) const{
  double pt_thresh(10.0);
  if(level>=1) pt_thresh=20.0;

  if (fabs(pf_mus_eta->at(k)) >= 2.4 ) return false;
  if (pf_mus_pt->at(k) < pt_thresh) return false;
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
  if (pf_mus_rel_iso > 0.2 && use_iso) return false;
  return true;
}

bool EventHandler::isTau(const unsigned int k,
			 const unsigned short level,
			 const bool require_iso) const{
  double pt_cut(20.0);
  if(level>0) pt_cut=30.0;

  if (taus_pt->at(k)<pt_cut) return false;
  if (fabs(taus_eta->at(k)) > 2.4) return false;
  if (require_iso && (taus_byLooseIsolationDeltaBetaCorr->at(k) <= 0)) return false;
  return true;
}

int EventHandler::GetNumElectrons(const unsigned short level, const bool use_iso) const{
  int count(0);
  for(unsigned int i(0); i<pf_els_pt->size(); ++i){
    if(isElectron(i,level, use_iso)) ++count;
  }
  return count;
}

int EventHandler::GetNumMuons(const unsigned short level, const bool use_iso) const{
  int count(0);
  for(unsigned int i(0); i<pf_mus_pt->size(); ++i){
    if(isMuon(i,level,use_iso)) ++count;
  }
  return count;
}

int EventHandler::GetNumTaus(const unsigned short level, const bool use_iso) const{
  int count(0);
  for(unsigned int i(0); i<taus_pt->size(); ++i){
    if(isTau(i,level,use_iso)) ++count;
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

double EventHandler::GetNumInteractions() const{
  double npv(-1.0);
  for(unsigned int i(0); i<PU_bunchCrossing->size(); ++i){
    if(PU_bunchCrossing->at(i)==0){
      npv = PU_TrueNumInteractions->at(i);
    }
  }
  return npv;
}

bool EventHandler::isGoodVertex(const unsigned int vertex) const{
  const double pv_rho(sqrt(pv_x->at(vertex)*pv_x->at(vertex) + pv_y->at(vertex)*pv_y->at(vertex)));
  return pv_ndof->at(vertex)>4 && fabs(pv_z->at(vertex))<24. && pv_rho<2.0 && pv_isFake->at(vertex)==0;
}

unsigned short EventHandler::GetNumVertices() const{
  unsigned short num_vertices(0);
  for(unsigned int vertex(0); vertex<pv_x->size(); ++vertex){
    if(isGoodVertex(vertex)) ++num_vertices;
  }
  return num_vertices;
  //return pv_x->size();
}

double EventHandler::GetPUWeight(reweight::LumiReWeighting &lumiWeights) const{
  return lumiWeights.weight(GetNumInteractions());
}

double EventHandler::GetSbinWeight() const{
  if(sampleName.find("Run2012")!=std::string::npos){
    return 1.0;
  }else if(pfmets_fullSignif<50.0){
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

double EventHandler::look_up_scale_factor() const{
  if(GetLSPMass()!=1) return scaleFactor;
  std::string::size_type last_decimal(model_params->rfind("."));
  std::string::size_type last_space(model_params->rfind(" ",last_decimal));
  std::string cross_section_str(model_params->substr(last_space+1));
  double cross_section(atof(cross_section_str.c_str()));
  double num_events(0);
  switch(GetCharginoMass()){
  case 130:
    num_events=47890.0;
    break;
  case 150:
    num_events=415152.0;
    break;
  case 175:
    num_events=399589.0;
    break;
  case 200:
    num_events=387968.0;
    break;
  case 225:
    num_events=379120.0;
    break;
  case 250:
    num_events=152428.0;
    break;
  case 275:
    num_events=150227.0;
    break;
  case 300:
    num_events=146559.0;
    break;
  case 325:
    num_events=142987.0;
    break;
  case 350:
    num_events=80305.0;
    break;
  case 375:
    num_events=78645.0;
    break;
  case 400:
    num_events=78026.0;
    break;
  case 425:
    num_events=77734.0;
    break;
  case 450:
    num_events=77444.0;
    break;
  case 475:
    num_events=77069.0;
    break;
  case 500:
    num_events=76193.0;
    break;
  default:
    return scaleFactor;
  }
  //Should probably not hard code the lumi but this requires changing a few other things...
  return 19.399*cross_section/num_events;
}

double EventHandler::GetTopPt() const {
  double topPt(-1);
  for(unsigned int i(0); i<mc_doc_id->size(); ++i){
    // look for the *top* (not antitop) pT
    if(mc_doc_id->at(i)==6) topPt = mc_doc_pt->at(i);
    if(topPt>=0) break; //check to see if we're done
  }
  return topPt;
}

double EventHandler::GetTopPtWeightOfficial() const{ // New 11/07
  double topweight(-1.0);

  //official recipe from
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  //(2 string comparisons for every event is not so good for performance, but oh well)
  if (sampleName.find("TTJets")!=std::string::npos || sampleName.find("TT_")!=std::string::npos) {
    double topPt(-1);
    double topbarPt(-1);
    for(unsigned int i(0); i<mc_doc_id->size(); ++i){
      // look for the *top* (not antitop) pT
      if(mc_doc_id->at(i)==6) topPt = mc_doc_pt->at(i);
      if(mc_doc_id->at(i)==-6) topbarPt = mc_doc_pt->at(i);
      if(topPt>=0 && topbarPt>=0) break; //check to see if we're done
    }

    //SF(x)=exp(a+bx)
    const double a(0.156); //combined 8 TeV values
    const double b(-0.00137);

    //important choice -- I've decided to use the 400 GeV value for values above 400 GeV
    //an alternative would be to just blindly use the formula
    if (topPt >400) topPt=400;
    if (topbarPt >400) topbarPt=400;

    const double SFt(exp(a + b*topPt));
    const double SFtbar(exp(a + b*topbarPt));

    topweight = sqrt( SFt * SFtbar );
  }

  if(topweight<0) topweight=1;
  return topweight;
}

double EventHandler::GetHighestJetPt(const unsigned int nth_highest) const{
  std::vector<double> pts(0);
  for(unsigned int jet(0); jet<jets_AK5PF_pt->size(); ++jet){
    if(isGoodJet(jet)){
      pts.push_back(jets_AK5PF_pt->at(jet));
    }
  }
  std::sort(pts.begin(), pts.end(), std::greater<double>());
  if(nth_highest<=pts.size()){
    return pts.at(nth_highest-1);
  }else{
    return 0.0;
  }
}

double EventHandler::GetHighestJetCSV(const unsigned int nth_highest) const{
  std::vector<double> csvs(0);
  for(unsigned int jet(0); jet<jets_AK5PF_btag_secVertexCombined->size(); ++jet){
    if(isGoodJet(jet)){
      csvs.push_back(jets_AK5PF_btag_secVertexCombined->at(jet));
    }
  }
  std::sort(csvs.begin(), csvs.end(), std::greater<double>());
  if(nth_highest<=csvs.size()){
    return csvs.at(nth_highest-1);
  }else{
    return 0.0;
  }
}
