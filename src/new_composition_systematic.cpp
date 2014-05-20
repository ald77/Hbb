#include "new_composition_systematic.hpp"

#include <random>
#include <array>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "TH1D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "style.hpp"
#include "utils.hpp"
#include "math.hpp"

int main(){
  const unsigned num_reps(100000);
  const unsigned num_bins(100);
  SetStyle();
  gStyle->SetPadTopMargin(0.05);  
  gROOT->ForceStyle();

  std::mt19937 rng;
  initialize_prng(rng);
  
  std::ifstream file_in("counts_norm_fix.log");
  float qcd_mean[3][2][4], ttbar_mean[3][2][4];
  for(unsigned k(0); k<4; ++k){
    for(unsigned i(0); i<3; ++i){
      for(unsigned j(0); j<2; ++j){
	file_in >> qcd_mean[i][j][k];
      }
    }
  }
  for(unsigned k(0); k<4; ++k){
    for(unsigned i(0); i<3; ++i){
      for(unsigned j(0); j<2; ++j){
	file_in >> ttbar_mean[i][j][k];
      }
    }
  }
  file_in.close();

  std::vector<float> qcd_central[3][2][4], ttbar_central[3][2][4];

  std::gamma_distribution<float> qcd_gd[3][2][4], ttbar_gd[3][2][4];
  for(unsigned i(0); i<3; ++i){
    for(unsigned j(0); j<2; ++j){
      for(unsigned k(0); k<4; ++k){
	const std::gamma_distribution<float> qcd_gd_temp(qcd_mean[i][j][k]+1.0, 1.0);
	qcd_gd[i][j][k].param(qcd_gd_temp.param());
	const std::gamma_distribution<float> ttbar_gd_temp(ttbar_mean[i][j][k]+1.0, 1.0);
	ttbar_gd[i][j][k].param(ttbar_gd_temp.param());
	qcd_central[i][j][k].push_back(qcd_mean[i][j][k]);
	ttbar_central[i][j][k].push_back(ttbar_mean[i][j][k]);
      }
    }
  }

  std::vector<float> qcd_dist[3][2][4], ttbar_dist[3][2][4];
  for(unsigned i(0); i<3; ++i){
    for(unsigned j(0); j<2; ++j){
      for(unsigned k(0); k<4; ++k){
	qcd_dist[i][j][k].resize(num_reps);
	ttbar_dist[i][j][k].resize(num_reps);
	for(unsigned rep(0); rep<num_reps; ++rep){
	  qcd_dist[i][j][k][rep]=qcd_gd[i][j][k](rng);
	  ttbar_dist[i][j][k][rep]=ttbar_gd[i][j][k](rng);
	}
      }
    }
  }

  TGraphAsymmErrors kappa230(num_bins);
  TGraphAsymmErrors kappa231(num_bins);
  TGraphAsymmErrors kappa232(num_bins);
  TGraphAsymmErrors kappa233(num_bins);
  TGraphAsymmErrors kappa234(num_bins);
  TGraphAsymmErrors kappa240(num_bins);
  TGraphAsymmErrors kappa241(num_bins);
  TGraphAsymmErrors kappa242(num_bins);
  TGraphAsymmErrors kappa243(num_bins);
  TGraphAsymmErrors kappa244(num_bins);
  TGraphAsymmErrors kappa340(num_bins);
  TGraphAsymmErrors kappa341(num_bins);
  TGraphAsymmErrors kappa342(num_bins);
  TGraphAsymmErrors kappa343(num_bins);
  TGraphAsymmErrors kappa344(num_bins);

  for(unsigned model(0); model<4; ++model){
    switch(model){
    case 0:
      
      break;
    case 1:
      break;
    case 2:
      break;
    case 3:
      break;
    default:
      break;
    }

    std::vector<float> ttbar_indep[3][2][4], qcd_indep[3][2][4];
    get_independence_model(ttbar_indep, ttbar_dist);
    get_nb_sbin_model(qcd_indep, qcd_dist);

    set_content(kappa230, qcd_dist, ttbar_indep, 2, 3, 0);
    set_content(kappa231, qcd_dist, ttbar_indep, 2, 3, 1);
    set_content(kappa232, qcd_dist, ttbar_indep, 2, 3, 2);
    set_content(kappa233, qcd_dist, ttbar_indep, 2, 3, 3);
    set_content(kappa234, qcd_dist, ttbar_indep, 2, 3, 4);
    set_content(kappa240, qcd_dist, ttbar_indep, 2, 4, 0);
    set_content(kappa241, qcd_dist, ttbar_indep, 2, 4, 1);
    set_content(kappa242, qcd_dist, ttbar_indep, 2, 4, 2);
    set_content(kappa243, qcd_dist, ttbar_indep, 2, 4, 3);
    set_content(kappa244, qcd_dist, ttbar_indep, 2, 4, 4);
    set_content(kappa340, qcd_dist, ttbar_indep, 3, 4, 0);
    set_content(kappa341, qcd_dist, ttbar_indep, 3, 4, 1);
    set_content(kappa342, qcd_dist, ttbar_indep, 3, 4, 2);
    set_content(kappa343, qcd_dist, ttbar_indep, 3, 4, 3);
    set_content(kappa344, qcd_dist, ttbar_indep, 3, 4, 4);
    make_plot(kappa230, kappa240, kappa340, "0_raw.pdf");
    make_plot(kappa231, kappa241, kappa341, "1_raw.pdf");
    make_plot(kappa232, kappa242, kappa342, "2_raw.pdf");
    make_plot(kappa233, kappa243, kappa343, "3_raw.pdf");
    make_plot(kappa234, kappa244, kappa344, "4_raw.pdf");
    set_content(kappa230, qcd_indep, ttbar_indep, 2, 3, 0);
    set_content(kappa231, qcd_indep, ttbar_indep, 2, 3, 1);
    set_content(kappa232, qcd_indep, ttbar_indep, 2, 3, 2);
    set_content(kappa233, qcd_indep, ttbar_indep, 2, 3, 3);
    set_content(kappa234, qcd_indep, ttbar_indep, 2, 3, 4);
    set_content(kappa240, qcd_indep, ttbar_indep, 2, 4, 0);
    set_content(kappa241, qcd_indep, ttbar_indep, 2, 4, 1);
    set_content(kappa242, qcd_indep, ttbar_indep, 2, 4, 2);
    set_content(kappa243, qcd_indep, ttbar_indep, 2, 4, 3);
    set_content(kappa244, qcd_indep, ttbar_indep, 2, 4, 4);
    set_content(kappa340, qcd_indep, ttbar_indep, 3, 4, 0);
    set_content(kappa341, qcd_indep, ttbar_indep, 3, 4, 1);
    set_content(kappa342, qcd_indep, ttbar_indep, 3, 4, 2);
    set_content(kappa343, qcd_indep, ttbar_indep, 3, 4, 3);
    set_content(kappa344, qcd_indep, ttbar_indep, 3, 4, 4);
    make_plot(kappa230, kappa240, kappa340, "0_fix.pdf");
    make_plot(kappa231, kappa241, kappa341, "1_fix.pdf");
    make_plot(kappa232, kappa242, kappa342, "2_fix.pdf");
    make_plot(kappa233, kappa243, kappa343, "3_fix.pdf");
    make_plot(kappa234, kappa244, kappa344, "4_fix.pdf");

    std::vector<float> lk231, lk232, lk233, lk234;
    std::vector<float> lk241, lk242, lk243, lk244;

    get_log_kappa(lk231, ttbar_dist, 1, 0, 0);
    get_log_kappa(lk232, ttbar_dist, 1, 0, 1);
    get_log_kappa(lk233, ttbar_dist, 1, 0, 2);
    get_log_kappa(lk234, ttbar_dist, 1, 0, 3);
    get_log_kappa(lk241, ttbar_dist, 2, 0, 0);
    get_log_kappa(lk242, ttbar_dist, 2, 0, 1);
    get_log_kappa(lk243, ttbar_dist, 2, 0, 2);
    get_log_kappa(lk244, ttbar_dist, 2, 0, 3);
    make_band_plot(lk231, lk232, lk233, lk234, lk241, lk242, lk243, lk244, "ttbar_band.pdf");
    get_log_kappa(lk231, qcd_dist, 1, 0, 0);
    get_log_kappa(lk232, qcd_dist, 1, 0, 1);
    get_log_kappa(lk233, qcd_dist, 1, 0, 2);
    get_log_kappa(lk234, qcd_dist, 1, 0, 3);
    get_log_kappa(lk241, qcd_dist, 2, 0, 0);
    get_log_kappa(lk242, qcd_dist, 2, 0, 1);
    get_log_kappa(lk243, qcd_dist, 2, 0, 2);
    get_log_kappa(lk244, qcd_dist, 2, 0, 3);
    make_band_plot(lk231, lk232, lk233, lk234, lk241, lk242, lk243, lk244, "qcd_band.pdf");

    print_table(ttbar_dist, "ttbar_dist_table.tex");
    print_table(qcd_dist, "qcd_dist_table.tex");
    print_table(ttbar_indep, "ttbar_indep_table.tex");
    print_table(qcd_indep, "qcd_indep_table.tex");
  }
}

void set_content(TGraphAsymmErrors& h,
		 const std::vector<float> qcd[3][2][4],
		 const std::vector<float> ttbar[3][2][4],
		 const unsigned num_b,
		 const unsigned den_b,
		 const unsigned sbin){
  std::vector<float> counts[3][2][4];
  fix_sizes(counts, qcd, ttbar);
  std::vector<float> log_kappa;
  const unsigned num_bins(h.GetN());
  const float width(1.0/num_bins);
  for(unsigned bin(1); bin<=num_bins; ++bin){
    const float x((bin-0.5)*width);
    merge(counts, qcd, ttbar, x);
    get_log_kappa(log_kappa, counts, num_b-2, den_b-2, sbin);

    const float y(Math::Mean(log_kappa.begin(), log_kappa.end()));
    float uncert_up(0.0), uncert_down(0.0);
    get_uncertainties(uncert_up, uncert_down, y, log_kappa);
    h.SetPoint(bin, x, y);
    h.SetPointError(bin, 0.5*width, 0.5*width, uncert_down, uncert_up);
  }
}

void fix_sizes(std::vector<float> out[3][2][4],
	       const std::vector<float> a[3][2][4],
	       const std::vector<float> b[3][2][4]){
  for(unsigned i(0); i<3; ++i){
    for(unsigned j(0); j<2; ++j){
      for(unsigned k(0); k<4; ++k){
	out[i][j][k].resize(std::min(a[i][j][k].size(), b[i][j][k].size()));
      }
    }
  }
}

void merge(std::vector<float> out[3][2][4],
	   const std::vector<float> a[3][2][4],
	   const std::vector<float> b[3][2][4],
	   const float x){
  const unsigned the_size(out[0][0][0].size());
  for(unsigned i(0); i<3; ++i){
    for(unsigned j(0); j<2; ++j){
      for(unsigned k(0); k<4; ++k){
	for(unsigned l(0); l<the_size; ++l){
	  out[i][j][k][l]=x*a[i][j][k][l]+(1.0-x)*b[i][j][k][l];
	}
      }
    }
  }
}

void get_log_kappa(std::vector<float>& lk,
		   const std::vector<float> count[3][2][4],
		   const unsigned num_idx,
		   const unsigned den_idx,
		   const unsigned sbin){
  lk.resize(count[0][0][0].size());
  float a(0.0), b(0.0), c(0.0), d(0.0);
  if(sbin<4){
    for(unsigned i(0); i<count[0][0][0].size(); ++i){
      a=count[num_idx][0][sbin].at(i);
      b=count[num_idx][1][sbin].at(i);
      c=count[den_idx][0][sbin].at(i);
      d=count[den_idx][1][sbin].at(i);
      lk.at(i)=log(a)+log(d)-log(b)-log(c);
    }
  }else{
    for(unsigned i(0); i<count[0][0][0].size(); ++i){
      a=count[num_idx][0][0].at(i)
	+count[num_idx][0][1].at(i)
	+count[num_idx][0][2].at(i)
	+count[num_idx][0][3].at(i);
      b=count[num_idx][1][0].at(i)
	+count[num_idx][1][1].at(i)
	+count[num_idx][1][2].at(i)
	+count[num_idx][1][3].at(i);
      c=count[den_idx][0][0].at(i)
	+count[den_idx][0][1].at(i)
	+count[den_idx][0][2].at(i)
	+count[den_idx][0][3].at(i);
      d=count[den_idx][1][0].at(i)
	+count[den_idx][1][1].at(i)
	+count[den_idx][1][2].at(i)
	+count[den_idx][1][3].at(i);
      lk.at(i)=log(a)+log(d)-log(b)-log(c);
    }
  }
}

void make_band_plot(std::vector<float>& ilk231,
		    std::vector<float>& ilk232,
		    std::vector<float>& ilk233,
		    std::vector<float>& ilk234,
		    std::vector<float>& ilk241,
		    std::vector<float>& ilk242,
		    std::vector<float>& ilk243,
		    std::vector<float>& ilk244,
		    const std::string& out_file){
  TH1D h("h", ";;ln(#kappa)", 8, -0.501, 7.501);
  h.GetXaxis()->SetBinLabel(1, "#kappa_{23}, bin 1");
  h.GetXaxis()->SetBinLabel(2, "#kappa_{23}, bin 2");
  h.GetXaxis()->SetBinLabel(3, "#kappa_{23}, bin 3");
  h.GetXaxis()->SetBinLabel(4, "#kappa_{23}, bin 4");
  h.GetXaxis()->SetBinLabel(5, "#kappa_{24}, bin 1");
  h.GetXaxis()->SetBinLabel(6, "#kappa_{24}, bin 2");
  h.GetXaxis()->SetBinLabel(7, "#kappa_{24}, bin 3");
  h.GetXaxis()->SetBinLabel(8, "#kappa_{24}, bin 4");
  TGraphAsymmErrors g(8);
  float vals[8], ups[8], downs[8], xs[8], zeros[8];
  for(unsigned i(0); i<8; ++i){
    xs[i]=i;
    zeros[i]=0.5;
  }
  get_val_and_uncert(vals, ups, downs, 0, ilk231);
  get_val_and_uncert(vals, ups, downs, 1, ilk232);
  get_val_and_uncert(vals, ups, downs, 2, ilk233);
  get_val_and_uncert(vals, ups, downs, 3, ilk234);
  get_val_and_uncert(vals, ups, downs, 4, ilk241);
  get_val_and_uncert(vals, ups, downs, 5, ilk242);
  get_val_and_uncert(vals, ups, downs, 6, ilk243);
  get_val_and_uncert(vals, ups, downs, 7, ilk244);
  TGraphAsymmErrors graph(8, xs, vals,
			  zeros, zeros,
			  downs, ups);
  const std::string title(";;ln(#kappa)");
  graph.SetTitle(title.c_str());
  const float max(get_maximum_with_error(graph));
  const float min(get_minimum_with_error(graph));
  h.SetMaximum(max+0.01);
  h.SetMinimum(min-0.01);
  h.SetStats(0);

  TCanvas canvas;
  h.Draw();
  graph.Draw("p");
  canvas.Print(out_file.c_str());
}

void get_val_and_uncert(float vals[8],
			float ups[8],
			float downs[8],
			unsigned idx,
			std::vector<float>& vec){
  vals[idx]=Math::Mean(vec.begin(), vec.end());
  get_uncertainties(ups[idx], downs[idx], vals[idx], vec);
}

void get_uncertainties(float& up,
		       float& down,
		       const float center,
		       std::vector<float>& vec){
  std::sort(vec.begin(), vec.end());
  const unsigned the_size(vec.size());
  const unsigned dpos(floor(0.682689492137086*the_size));
  float min_diff(std::numeric_limits<float>::max());
  unsigned best_low(0);
  float diff(0.0);
  for(unsigned low(0); low+dpos<the_size; ++low){
    diff=vec.at(low+dpos)-vec.at(low);
    if(diff<min_diff){
      best_low=low;
      min_diff=diff;
    }
  }
  up=vec.at(best_low+dpos)-center;
  down=center-vec.at(best_low);
  if(up<0.0) up=0.0;
  if(down<0.0) down=0.0;
}

void make_plot(TGraphAsymmErrors& g1,
	       TGraphAsymmErrors& g2,
	       TGraphAsymmErrors& g3,
	       const std::string& name){
  std::string title(";QCD Fraction;ln(#kappa)");
  TH1D h("crap", title.c_str(), 1, 0.0, 1.0);
  h.SetStats(0);
  h.SetLineColor(0);
  TCanvas canvas;
  g1.SetFillStyle(3003);
  g2.SetFillStyle(3004);
  g3.SetFillStyle(3005);
  g1.SetFillColor(2);
  g2.SetFillColor(3);
  g3.SetFillColor(4);
  g1.SetLineColor(2);
  g2.SetLineColor(3);
  g3.SetLineColor(4);
  double max(get_maximum_with_error(g1));
  max=std::max(max, get_maximum_with_error(g2));
  max=std::max(max, get_maximum_with_error(g3));
  double min(get_minimum_with_error(g1));
  min=std::min(min, get_minimum_with_error(g2));
  min=std::min(min, get_minimum_with_error(g3));
  h.SetMinimum(min-0.0001);
  h.SetMaximum(max+0.0001);
  h.Draw();
  g1.Draw("3same");
  g1.Draw("LXsame");
  g2.Draw("3same");
  g2.Draw("LXsame");
  g3.Draw("3same");
  g3.Draw("LXsame");
  TLegend legend(0.15,0.7,0.5,0.95);
  legend.AddEntry(&g1, "#kappa_{23}", "lf");
  legend.AddEntry(&g2, "#kappa_{24}", "lf");
  legend.AddEntry(&g3, "#kappa_{34}", "lf");
  legend.Draw("same");
  canvas.Print(name.c_str());
}

void get_independence_model(std::vector<float> out[3][2][4],
			    const std::vector<float> in[3][2][4]){
  for(unsigned index3(0); index3<4; ++index3){
    for(unsigned index1(0); index1<3; ++index1){
      for(unsigned index2(0); index2<2; ++index2){
	out[index1][index2][index3].resize(in[index1][index2][index3].size());
      }
    }
  }
  for(unsigned i(0); i<in[0][0][0].size(); ++i){
    float out1[3], out2[2], out3[4];
    std::vector<float> sum_vec(0);
    for(unsigned index1(0); index1<3; ++index1){
      out1[index1]=0.0;
      sum_vec.clear();
      for(unsigned index2(0); index2<2; ++index2){
	for(unsigned index3(0); index3<4; ++index3){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
      }
      out1[index1]=Math::Sum(sum_vec.begin(), sum_vec.end());
    }
    for(unsigned index2(0); index2<2; ++index2){
      out2[index2]=0.0;
      sum_vec.clear();
      for(unsigned index1(0); index1<3; ++index1){
	for(unsigned index3(0); index3<4; ++index3){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
      }
      out2[index2]=Math::Sum(sum_vec.begin(), sum_vec.end());
    }
    for(unsigned index3(0); index3<4; ++index3){
      out3[index3]=0.0;
      sum_vec.clear();
      for(unsigned index1(0); index1<3; ++index1){
	for(unsigned index2(0); index2<2; ++index2){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
      }
      out3[index3]=Math::Sum(sum_vec.begin(), sum_vec.end());
    }
    sum_vec.clear();
    for(unsigned index3(0); index3<4; ++index3){
      for(unsigned index1(0); index1<3; ++index1){
	for(unsigned index2(0); index2<2; ++index2){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
      }
    }
    float sum_sq_inv(Math::Sum(sum_vec.begin(), sum_vec.end()));
    sum_sq_inv=1.0/(sum_sq_inv*sum_sq_inv);
    for(unsigned index3(0); index3<4; ++index3){
      for(unsigned index1(0); index1<3; ++index1){
	for(unsigned index2(0); index2<2; ++index2){
	  out[index1][index2][index3].at(i)=out1[index1]*out2[index2]*out3[index3]*sum_sq_inv;
	}
      }
    }
  }
}

void get_mbb_sbin_model(std::vector<float> out[3][2][4],
			const std::vector<float> in[3][2][4]){
  for(unsigned index1(0); index1<3; ++index1){
    for(unsigned index2(0); index2<2; ++index2){
      for(unsigned index3(0); index3<4; ++index3){
	out[index1][index2][index3].resize(in[index1][index2][index3].size());
      }
    }
  }
  for(unsigned i(0); i<in[0][0][0].size(); ++i){
    for(unsigned index1(0); index1<3; ++index1){
      float out2[2], out3[4];
      std::vector<float> sum_vec(0);
      for(unsigned index2(0); index2<2; ++index2){
	out2[index2]=0.0;
	sum_vec.clear();
	for(unsigned index3(0); index3<4; ++index3){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
	out2[index2]=Math::Sum(sum_vec.begin(), sum_vec.end());
      }
      for(unsigned index3(0); index3<4; ++index3){
	out3[index3]=0.0;
	sum_vec.clear();
	for(unsigned index2(0); index2<2; ++index2){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
	out3[index3]=Math::Sum(sum_vec.begin(), sum_vec.end());
      }
      sum_vec.clear();
      for(unsigned index3(0); index3<4; ++index3){
	for(unsigned index2(0); index2<2; ++index2){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
      }
      const float sum(Math::Sum(sum_vec.begin(), sum_vec.end()));
      for(unsigned index3(0); index3<4; ++index3){
	for(unsigned index2(0); index2<2; ++index2){
	  out[index1][index2][index3].at(i)=out2[index1]*out3[index3]/sum;
	}
      }
    }
  }
}

void get_nb_sbin_model(std::vector<float> out[3][2][4],
		       const std::vector<float> in[3][2][4]){
  for(unsigned index3(0); index3<4; ++index3){
    for(unsigned index1(0); index1<3; ++index1){
      for(unsigned index2(0); index2<2; ++index2){
	out[index1][index2][index3].resize(in[index1][index2][index3].size());
      }
    }
  }
  for(unsigned i(0); i<in[0][0][0].size(); ++i){
    for(unsigned index2(0); index2<2; ++index2){
      float out1[3], out3[4];
      std::vector<float> sum_vec(0);
      for(unsigned index1(0); index1<3; ++index1){
	out1[index1]=0.0;
	sum_vec.clear();
	for(unsigned index3(0); index3<4; ++index3){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
	out1[index1]=Math::Sum(sum_vec.begin(), sum_vec.end());
      }
      for(unsigned index3(0); index3<4; ++index3){
	out3[index3]=0.0;
	sum_vec.clear();
	for(unsigned index1(0); index1<3; ++index1){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
	out3[index3]=Math::Sum(sum_vec.begin(), sum_vec.end());
      }
      sum_vec.clear();
      for(unsigned index3(0); index3<4; ++index3){
	for(unsigned index1(0); index1<3; ++index1){
	  sum_vec.push_back(in[index1][index2][index3].at(i));
	}
      }
      const float sum(Math::Sum(sum_vec.begin(), sum_vec.end()));
      for(unsigned index3(0); index3<4; ++index3){
	for(unsigned index1(0); index1<3; ++index1){
	  out[index1][index2][index3].at(i)=out1[index1]*out3[index3]/sum;
	}
      }
    }
  }
}

void print_table(std::vector<float> x[3][2][4],
		 const std::string& file_name){
  std::ofstream ofs(file_name);
  for(unsigned k(0); k<4; ++k){
    for(unsigned j(0); j<2; ++j){
      for(unsigned i(0); i<3; ++i){
	const float mid(Math::HalfSampleMode(x[i][j][k].begin(), x[i][j][k].end()));
	float down(0.0), up(0.0);
	get_uncertainties(up, down, mid, x[i][j][k]);
	ofs << '$' << fix_width(mid,4) << "^{+" << fix_width(up,4) << "}_{-" << fix_width(down,4) << "}$";
	if(!(j==1 && i==2)) ofs << " & ";
      }
    }
    if(k!=3) ofs << " \\\\";
    ofs << std::endl;
  }
  ofs.close();
}

void print_table(std::vector<float> x[3][2][4],
		 const float mid[3][2][4],
		 const std::string& file_name){
  std::ofstream ofs(file_name);
  for(unsigned k(0); k<4; ++k){
    for(unsigned j(0); j<2; ++j){
      for(unsigned i(0); i<3; ++i){
	float down(0.0), up(0.0);
	get_uncertainties(up, down, mid[i][j][k], x[i][j][k]);
	ofs << '$' << fix_width(mid[i][j][k],4) << "^{+" << fix_width(up,4) << "}_{-" << fix_width(down,4) << "}$";
	if(!(j==1 && i==2)) ofs << " & ";
      }
    }
    if(k!=3) ofs << " \\\\";
    ofs << std::endl;
  }
  ofs.close();
}

//All independent: 7
//nb-sbin independent: 12
//mbb-sbin independent: 15
//mbb-nb independent: 16
//No independent: 24
