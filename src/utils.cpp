#include "utils.hpp"

void MakeRatioPlot(std::vector<TH1D>& histos, std::vector<std::string>& names,
                   const std::string& out_name){
  if(histos.size()==0 || names.size()!=histos.size()) return;
  const double div(1.0/3.0);
  std::vector<unsigned> colors(0);
  colors.push_back(kRed-4);
  colors.push_back(kCyan-4);
  colors.push_back(kGreen-4);
  colors.push_back(kYellow-4);
  colors.push_back(kMagenta-4);
  colors.push_back(kBlue-4);
  colors.push_back(kOrange-4);
  colors.push_back(kGray);
  std::ostringstream oss("");
  oss << histos.at(0).GetTitle() << ";" << histos.at(0).GetXaxis()->GetTitle() << ";"
      << histos.at(0).GetYaxis()->GetTitle() << std::endl;
  THStack s((std::string("s_")+histos.at(0).GetName()).c_str(),
            oss.str().c_str());
  std::vector<double> edges(histos.at(0).GetNbinsX()+1);
  for(unsigned bin(0); bin<edges.size(); ++bin){
    edges.at(bin)=histos.at(0).GetBinLowEdge(bin+1);
  }
  TH1D* tot(static_cast<TH1D*>(histos.at(0).Clone("abcdefgh")));

  for(int bin(0); bin<=tot->GetNbinsX()+1; ++bin){
    tot->SetBinContent(bin,0.0);
  }
  for(unsigned histo(0); histo<histos.size(); ++histo){
    histos.at(histo).SetTitleSize(0,"X");
    histos.at(histo).SetTitleSize(0.08,"Y");
    //histos.at(histo).SetTitleOffset(0.0,"Y");
  }
  for(unsigned histo(1); histo<histos.size(); ++histo){
    histos.at(histo).SetFillColor(colors.at((histo-1)%colors.size()));
    histos.at(histo).SetLineColor(colors.at((histo-1)%colors.size()));
    s.Add(&histos.at(histo));
    tot->Add(&histos.at(histo));
  }

  double the_min(0.0);
  if(histos.size()>1){
    TH1D temp(histos.at(1));
    the_min=histos.at(1).GetMinimum(0.0);
    for(unsigned h(2); h<histos.size(); ++h){
      double this_min(histos.at(h).GetMinimum(0.0));
      if(this_min<the_min) the_min=this_min;
    }
  }
  const double the_max(std::max(s.GetMaximum(), histos.at(0).GetMaximum()));
  s.SetMinimum(the_min);
  histos.at(0).SetMinimum(the_min);
  s.SetMaximum(the_max);
  histos.at(0).SetMaximum(the_max);

  TH1D rat(histos.at(0));
  rat.SetTitle("");
  rat.GetYaxis()->SetTitle("MC/Data Ratio");
  rat.Divide(tot, &histos.at(0));
  rat.SetStats(0);

  rat.SetTitleSize(0.12,"X");
  rat.SetTitleSize(0.12,"Y");
  rat.SetLabelSize(0.1,"X");
  rat.SetLabelSize(0.1,"Y");
  rat.SetTitleOffset(0.32,"Y");

  TLegend l(0.7, 0.5, 0.95, 1.0-0.15/(1.0-div));
  for(unsigned histo(0); histo<histos.size(); ++histo){
    std::string marker_string(histo==0?"lpe":"f");
    l.AddEntry(&histos.at(histo), names.at(histo).c_str(), marker_string.c_str());
  }

  TCanvas c("c","c");
  c.Divide(2);
  c.cd(1);
  c.GetPad(1)->SetPad(0.0,div,1.0,1.0);
  c.GetPad(1)->SetMargin(0.15,0.05,0.0,0.15/(1.0-div));
  c.GetPad(1)->SetLogy(1);
  s.Draw("hist");
  histos.at(0).Draw("samee1");
  l.Draw("same");
  c.cd(2);
  c.GetPad(2)->SetPad(0.0,0.0,1.0,div);
  c.GetPad(2)->SetMargin(0.15,0.05,0.15/div,0.0);
  c.GetPad(2)->SetGridy(1);
  c.GetPad(2)->SetLogy(1);
  rat.Draw("e1");
  c.cd(0);
  c.Print(out_name.c_str());
}
