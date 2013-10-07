#ifndef H_EVENTHANDLER
#define H_EVENTHANDLER

#include <vector>
#include <string>
#include <utility>
#include <cfloat>
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
  EventHandler(const std::string &, const bool, const double=1.0);
  void Skim(const std::string &);
  void MakePlots(const std::string &);

  void SetScaleFactor(const double);
  void SetScaleFactor(const double, const double, const int);

private:
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

  double GetPUWeight(reweight::LumiReWeighting &) const;

  bool PassesPVCut() const;
  bool PassesMETCleaningCut() const;
  bool PassesTriggerCut() const;
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
  bool PassesMETSig50Cut() const;
  bool PassesMETSig100Cut() const;
  bool PassesMETSig150Cut() const;

  bool PassesRegionACut() const;
  bool PassesRegionBCut() const;
  bool PassesRegionC3bCut() const;
  bool PassesRegionD3bCut() const;
  bool PassesRegionC2bCut() const;
  bool PassesRegionD2bCut() const;

  bool PassesBadJetFilter() const;

  bool PassesTChiZHMassCut() const;

  void GetHiggsBJetPairing() const;
  void GetSortedBJets() const;

  std::pair<double, double> GetHiggsMasses() const;
  double GetHiggsDeltaR() const;

  unsigned int GetNumLowPtPfCands(const double=20.0) const;
  double GetMETOfLowPtPfCands(const double=20.0) const;

  int GetPBNR() const;
  double GetMinDeltaPhiMET(const unsigned int) const;

  int GetNumGoodJets() const;
  int GetNumCSVTJets() const;
  int GetNumCSVMJets() const;
  int GetNumCSVLJets() const;

  bool isGoodJet(const unsigned int, const bool=true, const double=20.0, const double=2.4) const;
  bool isProblemJet(const unsigned int) const;
  bool jetPassLooseID(const unsigned int) const;

  bool isVetoElectron(const unsigned int) const;
  bool isVetoMuon(const unsigned int) const;
  bool isVetoTau(const unsigned int) const;

  int GetNumVetoElectrons() const;
  int GetNumVetoMuons() const;
  int GetNumVetoTaus() const;

  bool isIsoTrack(const unsigned int, const double=10.0) const;
  bool isQualityTrack(const unsigned int) const;

  int GetNumIsoTracks(const double=10.0) const;
  int NewGetNumIsoTracks(const double=10.0) const;

  double GetElectronRelIso(const unsigned int) const;

  double GetSbinWeight() const;
  double GetTopPtWeight() const;

  double GetMaxDR() const;
  double GetHT(const bool=true, const bool=false) const;
  double GetHighestBTag(const unsigned int=1) const;

  bool HasGluonSplitting() const;

  std::vector<std::pair<int,int> > GetBOrigins() const;
};

#endif
