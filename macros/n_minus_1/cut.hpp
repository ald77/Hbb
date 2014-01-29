#ifndef __cut_h__
#define __cut_h__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>

#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "THStack.h"
#include "TCut.h"
#include "TString.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TString.h"

using namespace std;

class Cut {
 public:
  Cut(TString var, double minVal=-99999., double maxVal=99999.);
  ~Cut();

  TString GetVar() const;
  double GetMinVal() const;
  double GetMaxVal() const;
  TString Get_Var_Cut() const;
  TCut Get_N_Minus_1_Cuts() const; 
  //double Get_Cut_Val() const;

 private:
  TString var;
  double minVal;
  double maxVal;
  static map<TString, TString> varToBoolTable;
  //static map<TString, double> varToValTable;

  void SetCuts();
  void SetVals();
};


#endif
