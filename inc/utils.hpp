#ifndef H_UTILS
#define H_UTILS

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "THStack.h"
#include "TChain.h"

template<typename T>
void setup(TChain& chain, const std::string& name, T& variable){
  chain.SetBranchStatus(name.c_str(), 1);
  chain.SetBranchAddress(name.c_str(), &variable);
}

void MakeRatioPlot(std::vector<TH1D>& histos, std::vector<std::string>& names,
		   const std::string& out_name);

#endif
