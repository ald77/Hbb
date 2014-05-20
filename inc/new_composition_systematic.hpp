#ifndef H_COMPOSITION_SYSTEMATIC
#define H_COMPOSITION_SYSTEMATIC

#include <array>
#include <vector>
#include <string>
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

void make_plot(TGraphAsymmErrors& g1,
	       TGraphAsymmErrors& g2,
	       TGraphAsymmErrors& g3,
	       const std::string& name);

void set_content(TGraphAsymmErrors& h,
		 const std::vector<float> qcd[3][2][4],
		 const std::vector<float> ttbar[3][2][4],
		 const unsigned num_b,
		 const unsigned den_b,
		 const unsigned sbin);

void merge(std::vector<float> out[3][2][4],
	   const std::vector<float> a[3][2][4],
	   const std::vector<float> b[3][2][4],
	   const float x);

void fix_sizes(std::vector<float> out[3][2][4],
	       const std::vector<float> a[3][2][4],
	       const std::vector<float> b[3][2][4]);

void get_log_kappa(std::vector<float>& lk,
		     const std::vector<float> count[3][2][4],
		     const unsigned num_idx,
		     const unsigned den_idx,
		     const unsigned sbin);

void make_band_plot(std::vector<float>& lk231,
		    std::vector<float>& lk232,
		    std::vector<float>& lk233,
		    std::vector<float>& lk234,
		    std::vector<float>& lk241,
		    std::vector<float>& lk242,
		    std::vector<float>& lk243,
		    std::vector<float>& lk244,
		    const std::string& title);

void get_val_and_uncert(float vals[8],
			float ups[8],
			float downs[8],
			unsigned idx,
			std::vector<float>& vec);

void get_uncertainties(float& up,
		       float& down,
		       const float center,
		       std::vector<float>& vec);

void get_independence_model(std::vector<float> out[3][2][4],
			    const std::vector<float> in[3][2][4]);

void get_nb_sbin_model(std::vector<float> out[3][2][4],
			    const std::vector<float> in[3][2][4]);

void print_table(std::vector<float> x[3][2][4],
		 const std::string& file_name);

#endif
