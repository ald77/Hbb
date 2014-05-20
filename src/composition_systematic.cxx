#include "composition_systematic.hpp"

#include <string>
#include <algorithm>
#include <fstream>
#include "TH1D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "style.hpp"
#include "utils.hpp"
#include "math.hpp"

int main(){
  SetStyle();
  gStyle->SetPadTopMargin(0.05);  
  gROOT->ForceStyle();

  std::ifstream file_in("counts_fixed_other_bjets.txt");
  double qcd[3][2][5], other_background[3][2][5];
  for(unsigned k(0); k<4; ++k){
    for(unsigned i(0); i<3; ++i){
      for(unsigned j(0); j<2; ++j){
	file_in >> qcd[i][j][k];
      }
    }
  }
  for(unsigned k(0); k<4; ++k){
    for(unsigned i(0); i<3; ++i){
      for(unsigned j(0); j<2; ++j){
	file_in >> other_background[i][j][k];
      }
    }
  }
  file_in.close();
  get_total(qcd);
  get_total(other_background);
  
  const unsigned num_bins(1000);
  const double upper(8.0);
  TH1D k230("k230", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k231("k231", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k232("k232", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k233("k233", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k234("k234", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k240("k240", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k241("k241", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k242("k242", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k243("k243", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k244("k244", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k340("k340", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k341("k341", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k342("k342", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k343("k343", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);
  TH1D k344("k344", ";QCD Multiplier;log(#kappa)", num_bins, 0.0, upper);

  set_content(k230, qcd, other_background, 0, 0);
  set_content(k231, qcd, other_background, 0, 1);
  set_content(k232, qcd, other_background, 0, 2);
  set_content(k233, qcd, other_background, 0, 3);
  set_content(k234, qcd, other_background, 0, 4);
  set_content(k240, qcd, other_background, 1, 0);
  set_content(k241, qcd, other_background, 1, 1);
  set_content(k242, qcd, other_background, 1, 2);
  set_content(k243, qcd, other_background, 1, 3);
  set_content(k244, qcd, other_background, 1, 4);
  set_content(k340, qcd, other_background, 2, 0);
  set_content(k341, qcd, other_background, 2, 1);
  set_content(k342, qcd, other_background, 2, 2);
  set_content(k343, qcd, other_background, 2, 3);
  set_content(k344, qcd, other_background, 2, 4);
  k230.SetLineColor(1);
  k231.SetLineColor(1);
  k232.SetLineColor(1);
  k233.SetLineColor(1);
  k234.SetLineColor(1);
  k240.SetLineColor(2);
  k241.SetLineColor(2);
  k242.SetLineColor(2);
  k243.SetLineColor(2);
  k244.SetLineColor(2);
  k340.SetLineColor(3);
  k341.SetLineColor(3);
  k342.SetLineColor(3);
  k343.SetLineColor(3);
  k344.SetLineColor(3);

  TCanvas canvas;
  TLegend legend(0.5, 0.7, 0.95, 0.85);
  k230.Draw("c");
  k240.Draw("csame");
  k340.Draw("csame");
  legend.AddEntry(&k230, "#kappa_{23}", "l");
  legend.AddEntry(&k240, "#kappa_{24}", "l");
  legend.AddEntry(&k340, "#kappa_{34}", "l");
  plot(legend, k230, k240, k340, "compo/0.pdf");
  plot(legend, k231, k241, k341, "compo/1.pdf");
  plot(legend, k232, k242, k342, "compo/2.pdf");
  plot(legend, k233, k243, k343, "compo/3.pdf");
  plot(legend, k234, k244, k344, "compo/total.pdf");
}

void plot(TLegend& l, TH1D& h1, TH1D& h2, TH1D& h3, const std::string& name){
  const double ma1(get_maximum(h1));
  const double ma2(get_maximum(h2));
  const double ma3(get_maximum(h3));
  const double the_max(std::max(ma1,std::max(ma2,ma3)));
  h1.SetMaximum(the_max);
  h2.SetMaximum(the_max);
  h3.SetMaximum(the_max);
  const double mi1(get_minimum(h1));
  const double mi2(get_minimum(h2));
  const double mi3(get_minimum(h3));
  const double the_min(std::min(mi1,std::min(mi2,mi3)));
  h1.SetMinimum(the_min);
  h2.SetMinimum(the_min);
  h3.SetMinimum(the_min);
  h1.SetStats(0);
  h2.SetStats(0);
  h3.SetStats(0);
  TCanvas canvas;
  h1.Draw("c");
  h2.Draw("csame");
  h3.Draw("csame");
  l.Draw("same");
  canvas.Print(name.c_str());
}

void set_content(TH1D& histo,
		 const double qcd[3][2][5],
		 const double other_background[3][2][5],
		 const unsigned short kappa_type,
		 const unsigned short sbin){
  const int num_bins(histo.GetNbinsX());
  unsigned short num_b_index_1(0), num_b_index_2(0);
  switch(kappa_type){
  case 0:
    num_b_index_1=0;
    num_b_index_2=1;
    break;
  case 1:
    num_b_index_1=0;
    num_b_index_2=2;
    break;
  case 2:
    num_b_index_1=1;
    num_b_index_2=2;
    break;
  default:
    return;
  }
  for(int bin(1); bin<=num_bins; ++bin){
    const double x(histo.GetBinCenter(bin));
    const double a(other_background[num_b_index_2][0][sbin]
		   +x*qcd[num_b_index_2][0][sbin]);
    const double b(other_background[num_b_index_2][1][sbin]
		   +x*qcd[num_b_index_2][1][sbin]);
    const double c(other_background[num_b_index_1][0][sbin]
		   +x*qcd[num_b_index_1][0][sbin]);
    const double d(other_background[num_b_index_1][1][sbin]
		   +x*qcd[num_b_index_1][1][sbin]);
    std::vector<double> v(0);
    v.push_back(log(a));
    v.push_back(-log(b));
    v.push_back(-log(c));
    v.push_back(log(d));
    histo.SetBinContent(bin, Math::Sum(v.begin(), v.end()));
  }
}

void get_total(double x[3][2][5]){
  for(unsigned i(0); i<3; ++i){
    for(unsigned j(0); j<2; ++j){
      x[i][j][4]=x[i][j][0]+x[i][j][1]+x[i][j][2]+x[i][j][3];
    }
  }
}
