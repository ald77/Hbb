#ifndef H_PLOTTER
#define H_PLOTTER

#include <vector>
#include <string>
#include "TH1.h"

class plotter{
public:
  plotter();

  void plot_data_mc(TH1D& data,
		    std::vector<TH1D>& mc,
		    const std::string& output_name) const;
  void plot_data_mc_sig(TH1D& data,
			std::vector<TH1D>& mc,
			TH1D& signal,
			const std::string& output_name) const;

  void set_mc_names(const std::vector<std::string>& mc_names);

private:
  unsigned short canvas_width_, canvas_height_;
  float legend_left_, legend_right_, legend_bottom_, legend_top_;
  float vertical_divide_;
  std::vector<std::string> mc_names_;
};

void plot_data_mc(const plotter& plot,
		    std::vector<TH1D>& histos,
		    const std::string& output_name);

void plot_data_mc_sig(const plotter& plot,
			std::vector<TH1D>& histos,
			const std::string& output_name);

#endif
