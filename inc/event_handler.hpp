#ifndef H_EVENTHANDLER
#define H_EVENTHANDLER

#include <cfloat>
#include <vector>
#include <string>
#include <utility>
#include <limits>
#include <stdint.h>
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "pu_constants.hpp"
#include "lumi_reweighting_stand_alone.hpp"
#include "cfa.hpp"
#include "b_jet.hpp"
#include "event_number.hpp"

class EventHandler : public cfA{
public:
  EventHandler(const std::string &, const bool, const double=1.0, const bool=false);

  void SetScaleFactor(const double);
  void SetScaleFactor(const double, const double, const int);

  enum FailureType{
    kGood = 0u,
    kMETSig150 = 1u<<0,
    kMETSig100 = 1u<<1,
    kMETSig50 = 1u<<2,
    kMETSig30 = 1u<<3,
    kDeltaR = 1u<<4,
    kHiggsAvgMass = 1u<<5,
    kHiggsMassDiff = 1u<<6,
    k4thBTag = 1u<<7,
    k3rdBTag = 1u<<8,
    k2ndBTag = 1u<<9,
    kIsoTrackVeto = 1u<<10,
    kLeptonVeto = 1u<<11,
    kMinDeltaPhi = 1u<<12,
    kNumJets = 1u<<13,
    kTrigger = 1u<<14,
    kJSON = 1u<<15,
    kMETCleaning = 1u<<16,
    kPV = 1u<<17,
    kTChiMassCut = 1u<<18
  };

protected:
  mutable std::pair<std::pair<TLorentzVector, TLorentzVector>, std::pair<TLorentzVector, TLorentzVector> > higgsBJetPairing;//caching for efficiency
  mutable std::vector<BJet> sortedBJetCache;//caching for efficiency
  mutable bool higgsPairingUpToDate, bJetsUpToDate;//cached value correct
  mutable bool betaUpToDate;
  static const double CSVTCut, CSVMCut, CSVLCut;
  double scaleFactor;
  mutable std::vector<double> beta;

  int GetcfAVersion() const;

  void GetEntry(const unsigned int);

  void GetBeta(const std::string which="beta") const;

  double GetNumInteractions() const;
  unsigned short GetNumVertices() const;
  double GetPUWeight(reweight::LumiReWeighting &) const;

  bool isGoodVertex(const unsigned int) const;

  bool PassesBaselineSelection() const;
  bool PassesTriggerPlateauCuts() const;
  bool PassesPVCut() const;
  bool PassesMETCleaningCut() const;
  bool PassesTriggerCut() const;
  bool PassesQCDTriggerCut() const;
  bool PassesSpecificTrigger(const std::string) const;
  bool PassesNumJetsCut() const;
  bool Passes2CSVTCut() const;
  bool PassesJet2PtCut() const;
  bool PassesMinDeltaPhiCut() const;
  bool PassesMETSig30Cut() const;
  bool PassesLeptonVetoCut() const;
  bool PassesIsoTrackVetoCut() const;
  bool PassesBTaggingCut() const;
  bool PassesHiggsAvgMassCut() const;
  bool PassesHiggsMassDiffCut() const;
  bool PassesHiggsMassCut() const;
  bool PassesDRCut() const;
  bool PassesInvertedDRCut() const;
  bool PassesMETSig50Cut() const;
  bool PassesMETSig80Cut() const;
  bool PassesMETSig100Cut() const;
  bool PassesMETSig150Cut() const;

  bool PassesSingleLeptonCut() const;
  bool PassesSingleLeptonControlCut() const;
  bool PassesQCDControlCut() const;
  bool PassesJSONCut() const;

  uint_least32_t GetCutFailCode() const;

  bool PassesRegionACut() const;
  bool PassesRegionBCut() const;
  bool PassesRegionC3bCut() const;
  bool PassesRegionD3bCut() const;
  bool PassesRegionC2bCut() const;
  bool PassesRegionD2bCut() const;

  bool PassesInvertedDRRegionACut() const;
  bool PassesInvertedDRRegionBCut() const;
  bool PassesInvertedDRRegionC3bCut() const;
  bool PassesInvertedDRRegionD3bCut() const;
  bool PassesInvertedDRRegionC2bCut() const;
  bool PassesInvertedDRRegionD2bCut() const;

  bool PassesSingleLeptonRegionACut() const;
  bool PassesSingleLeptonRegionBCut() const;
  bool PassesSingleLeptonRegionC3bCut() const;
  bool PassesSingleLeptonRegionD3bCut() const;
  bool PassesSingleLeptonRegionC2bCut() const;
  bool PassesSingleLeptonRegionD2bCut() const;

  bool PassesBadJetFilter() const;

  bool PassesTChiMassCut(int=-1, int=1) const;

  int GetMass(const std::string&) const;
  int GetCharginoMass() const;
  int GetLSPMass() const;

  void GetHiggsBJetPairing() const;
  void GetSortedBJets() const;

  std::pair<double, double> GetHiggsMasses() const;
  double GetMinimaxMbb() const;
  double GetMaximinMbb() const;
  double GetHiggsDeltaR() const;

  unsigned int GetNumLowPtPfCands(const double=20.0) const;
  double GetMETOfLowPtPfCands(const double=20.0) const;

  int GetPBNR() const;
  double GetMinDeltaPhiMET(const unsigned int=std::numeric_limits<unsigned int>::max()) const;

  int GetNumGoodJets() const;
  int GetNumCSVTJets() const;
  int GetNumCSVMJets() const;
  int GetNumCSVLJets() const;
  int GetNumBTaggedJets() const;

  bool isGoodJet(const unsigned int, const bool=true, const double=20.0,
                 const double=2.4, const bool=true) const;
  bool isProblemJet(const unsigned int) const;
  bool jetPassLooseID(const unsigned int) const;

  bool isElectron(const unsigned int, const unsigned short=0, const bool=true) const;
  bool isMuon(const unsigned int, const unsigned short=0, const bool=true) const;
  bool isTau(const unsigned int, const unsigned short=0, const bool=true) const;

  int GetNumElectrons(const unsigned short=0, const bool=true) const;
  int GetNumMuons(const unsigned short=0, const bool=true) const;
  int GetNumTaus(const unsigned short=0, const bool=true) const;

  bool isIsoTrack(const unsigned int, const double=10.0) const;
  bool isQualityTrack(const unsigned int) const;

  int GetNumIsoTracks(const double=10.0) const;
  int NewGetNumIsoTracks(const double=10.0) const;

  double GetElectronRelIso(const unsigned int) const;

  double GetSbinWeight() const;
  double GetTopPtWeight() const;
  double GetTopPt() const;
  double GetTopPtWeightOfficial() const;

  double GetMinDR() const;
  double GetMaxDR() const;
  double GetHT(const bool=true, const bool=false) const;
  double GetHighestCSV(const unsigned int=1) const;
  
  int GetPartonIdBin(const float=0.) const;

  bool HasGluonSplitting() const;

  double look_up_scale_factor() const;

  std::vector<std::pair<int,int> > GetBOrigins() const;
  std::vector<std::pair<int,int> > GetBOrigins_new() const;

  double GetHighestJetPt(const unsigned int=1) const;
  double GetHighestJetCSV(const unsigned int=1) const;
};

#endif
