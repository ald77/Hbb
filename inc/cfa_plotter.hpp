#ifndef H_CFA_PLOTTER
#define H_CFA_PLOTTER

#include <string>
#include "event_handler.hpp"

class CfAPlotter : public EventHandler{
public:
  CfAPlotter(const std::string& in_file_name,
             const bool is_list,
             const double weight_in=1.0);

  void MakePlots(const std::string& out_file_name);
};

namespace CfAPlots{
  void FixBinLabels_new(TH1D &h);
  void FixBinLabels(TH1D &h);
  void FixSbinLabels(TH1D &h);
  int GetType_new(const std::pair<int,int> &b_origin);
  int GetType(const std::pair<int,int> &b_origin);
}

#endif
