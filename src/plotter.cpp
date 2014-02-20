#include <algorithm>
#include <vector>
#include <string>
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "plotter.hpp"
#include "utils.hpp"
#include "math.hpp"

void plot_data_mc(const plotter& plot,
		  std::vector<TH1D>& histos,
		  const std::string& output_name){
  if(histos.size()>0){
    std::vector<TH1D> mc(histos.begin()+1, histos.end());
    plot.plot_data_mc(histos.at(0), mc, output_name);
  }
}

void plot_data_mc_sig(const plotter& plot,
		      std::vector<TH1D>& histos,
		      const std::string& output_name){
  if(histos.size()>1){
    TH1D data(histos.at(0));
    std::vector<TH1D> mc(histos.begin()+1, histos.end()-1);
    TH1D sig(histos.at(histos.size()-1));
    plot.plot_data_mc_sig(data, mc, sig, output_name);
  }
}

plotter::plotter():
  canvas_width_(600),
  canvas_height_(600),
  legend_left_(0.7),
  legend_right_(0.95),
  legend_bottom_(0.5),
  legend_top_(0.85),
  vertical_divide_(1.0/3.0),
  mc_names_(0){
}

void plotter::set_mc_names(const std::vector<std::string>& mc_names){
  mc_names_=mc_names;
}
		    
void plotter::plot_data_mc(TH1D& data,
			   std::vector<TH1D>& mc,
			   const std::string& output_name) const{
  assign_color(data, 1);
  assign_colors(mc);

  const std::string name(data.GetName());
  std::string full_title("");
  get_full_title(data, full_title);
  const unsigned num_bins(data.GetNbinsX());
  std::vector<double> bin_edges(num_bins+1);
  for(std::vector<double>::size_type bin(0); bin<bin_edges.size(); ++bin){
    bin_edges.at(bin)=data.GetBinLowEdge(bin+1);
  }

  THStack stack(("s_"+name).c_str(), full_title.c_str());
  TH1D total(("h_"+name).c_str(), full_title.c_str(), num_bins, &bin_edges.at(0));
  assign_color(total, 1);

  for(std::vector<TH1D>::size_type histo(0); histo<mc.size(); ++histo){
    stack.Add(&mc.at(histo));
    total.Add(&mc.at(histo));
  }

  const double the_max(std::max(data.GetMaximum(), total.GetMaximum()));
  stack.SetMaximum(the_max);
  total.SetMaximum(the_max);

  if(mc.size()>0){
    double the_min(mc.at(0).GetMinimum(0.0));
    for(std::vector<TH1D>::size_type histo(1); histo<mc.size(); ++histo){
      const double this_min(mc.at(histo).GetMinimum(0.0));
      if(this_min<the_min) the_min=this_min;
    }
    stack.SetMinimum(the_min);
    total.SetMinimum(the_min);
  }

  TH1D residuals(data);
  residuals.SetTitle("");
  residuals.GetYaxis()->SetTitle("Residual");
  residuals.SetStats(0);
  for(unsigned bin(0); bin<=num_bins+1; ++bin){
    const double xd(data.GetBinContent(bin));
    const double xm(total.GetBinContent(bin));
    const double sd(fabs(data.GetBinError(bin)));
    const double sm(fabs(total.GetBinError(bin)));
    const double denom(Math::add_in_quadrature(sd, sm));
    double content(0.0);
    if(denom!=0.0){
      content=(xd-xm)/denom;
    }
    residuals.SetBinContent(bin, content);
    residuals.SetBinError(bin, 0.0);
  }

  residuals.SetTitleSize(0.12,"X");
  residuals.SetTitleSize(0.12,"Y");
  residuals.SetLabelSize(0.1,"X");
  residuals.SetLabelSize(0.1,"Y");
  residuals.SetTitleOffset(0.32,"Y");
  residuals.SetMarkerStyle(20);

  TLegend legend(legend_left_, legend_bottom_, legend_right_,
		 1.0-(1.0-legend_top_)/(1.0-vertical_divide_));
  legend.AddEntry(&data, "CMS Data", "lpe");
  for(std::vector<TH1D>::size_type histo(0); histo<mc.size(); ++histo){
    legend.AddEntry(&mc.at(histo), mc_names_.at(histo).c_str(), "f");
  }

  TCanvas canvas("","", canvas_width_, canvas_height_);
  canvas.Divide(2);
  canvas.cd(1);
  canvas.GetPad(1)->SetPad(0.0, vertical_divide_, 1.0, 1.0);
  canvas.GetPad(1)->SetMargin(0.15, 0.05, 0.0, 0.15/(1.0-vertical_divide_));
  canvas.GetPad(1)->SetLogy(1);
  stack.Draw("hist");
  data.Draw("samee1");
  legend.Draw("same");
  canvas.cd(2);
  canvas.GetPad(2)->SetPad(0.0, 0.0, 1.0, vertical_divide_);
  canvas.GetPad(2)->SetMargin(0.15, 0.05, 0.15/vertical_divide_, 0.0);
  canvas.GetPad(2)->SetGridy(1);
  residuals.Draw("histp");
  canvas.cd(0);
  canvas.Print(output_name.c_str());
}

void plotter::plot_data_mc_sig(TH1D& data,
			       std::vector<TH1D>& mc,
			       TH1D& sig,
			       const std::string& output_name) const{
  assign_color(data, 1);
  assign_colors(mc);
  assign_color(sig, 2);

  const std::string name(data.GetName());
  std::string full_title("");
  get_full_title(data, full_title);
  const unsigned num_bins(data.GetNbinsX());
  std::vector<double> bin_edges(num_bins+1);
  for(std::vector<double>::size_type bin(0); bin<bin_edges.size(); ++bin){
    bin_edges.at(bin)=data.GetBinLowEdge(bin+1);
  }

  THStack stack(("s_"+name).c_str(), full_title.c_str());
  TH1D total(("h_"+name).c_str(), full_title.c_str(), num_bins, &bin_edges.at(0));
  assign_color(total, 1);

  for(std::vector<TH1D>::size_type histo(0); histo<mc.size(); ++histo){
    stack.Add(&mc.at(histo));
    total.Add(&mc.at(histo));
  }

  const double the_max(std::max(data.GetMaximum(), total.GetMaximum()));
  stack.SetMaximum(the_max);
  total.SetMaximum(the_max);

  if(mc.size()>0){
    double the_min(mc.at(0).GetMinimum(0.0));
    for(std::vector<TH1D>::size_type histo(1); histo<mc.size(); ++histo){
      const double this_min(mc.at(histo).GetMinimum(0.0));
      if(this_min<the_min) the_min=this_min;
    }
    stack.SetMinimum(the_min);
    total.SetMinimum(the_min);
  }

  TH1D residuals(data);
  residuals.SetTitle("");
  residuals.GetYaxis()->SetTitle("Residual");
  residuals.SetStats(0);
  for(unsigned bin(0); bin<=num_bins+1; ++bin){
    const double xd(data.GetBinContent(bin));
    const double xm(total.GetBinContent(bin));
    const double sd(fabs(data.GetBinError(bin)));
    const double sm(fabs(total.GetBinError(bin)));
    const double denom(Math::add_in_quadrature(sd, sm));
    double content(0.0);
    if(denom!=0.0){
      content=(xd-xm)/denom;
    }
    residuals.SetBinContent(bin, content);
    residuals.SetBinError(bin, 0.0);
  }

  residuals.SetTitleSize(0.12,"X");
  residuals.SetTitleSize(0.12,"Y");
  residuals.SetLabelSize(0.1,"X");
  residuals.SetLabelSize(0.1,"Y");
  residuals.SetTitleOffset(0.32,"Y");
  residuals.SetMarkerStyle(20);

  TLegend legend(legend_left_, legend_bottom_, legend_right_,
		 1.0-(1.0-legend_top_)/(1.0-vertical_divide_));
  legend.AddEntry(&data, "CMS Data", "lpe");
  for(std::vector<TH1D>::size_type histo(0); histo<mc.size(); ++histo){
    legend.AddEntry(&mc.at(histo), mc_names_.at(histo).c_str(), "f");
  }

  sig=sig+total;

  TCanvas canvas("","", canvas_width_, canvas_height_);
  canvas.Divide(2);
  canvas.cd(1);
  canvas.GetPad(1)->SetPad(0.0, vertical_divide_, 1.0, 1.0);
  canvas.GetPad(1)->SetMargin(0.15, 0.05, 0.0, 0.15/(1.0-vertical_divide_));
  canvas.GetPad(1)->SetLogy(1);
  stack.Draw("hist");
  data.Draw("samee1");
  sig.Draw("samehist");
  legend.Draw("same");
  canvas.cd(2);
  canvas.GetPad(2)->SetPad(0.0, 0.0, 1.0, vertical_divide_);
  canvas.GetPad(2)->SetMargin(0.15, 0.05, 0.15/vertical_divide_, 0.0);
  canvas.GetPad(2)->SetGridy(1);
  residuals.Draw("histp");
  canvas.cd(0);
  canvas.Print(output_name.c_str());
}
