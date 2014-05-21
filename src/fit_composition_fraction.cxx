#include "fit_composition_fraction.hpp"

#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFractionFitter.h"
#include "TObjArray.h"
#include "utils.hpp"

int main(int argc, char* argv[]){
  if(argc>1){
    TH1::SetDefaultSumw2();
    std::ifstream file_in(argv[1]);
    float qcd_count[3][2][4], ttbar_count[3][2][4], obs_count[3][2][4];
    get_counts(file_in, qcd_count);
    get_counts(file_in, ttbar_count);
    get_counts(file_in, obs_count);
    file_in.close();

    TH1D qcd_histo("qcd_histo", ";Bin;Events", 24, 0.5, 24.5);
    TH1D ttbar_histo("ttbar_histo", ";Bin;Events", 24, 0.5, 24.5);
    TH1D obs_histo("obs_histo", ";Bin;Events", 24, 0.5, 24.5);

    setup_histo(qcd_histo, qcd_count);
    setup_histo(ttbar_histo, ttbar_count);
    setup_histo(obs_histo, obs_count);

    TObjArray obj_arr(2);
    obj_arr.Add(&qcd_histo);
    obj_arr.Add(&ttbar_histo);

    TFractionFitter fitter(&obs_histo, &obj_arr);
    fitter.Constrain(0, 0.0, 1.0);
    fitter.Constrain(1, 0.0, 1.0);

    int status(fitter.Fit());
    TCanvas c;
    if(status==0){
      TH1F* result(static_cast<TH1F*>(fitter.GetPlot()));
      obs_histo.Draw("Ep");
      result->Draw("same");
      c.Print("crap.pdf");
    }
    
    const int num_bins(1000);
    TH1D scan("scan", ";QCD Fraction;-2*log(L)", num_bins, 0.0, 1.0);
    for(int bin(1); bin<num_bins; ++bin){
      const double x(scan.GetBinCenter(bin));
      const double xp[2]={x,1.0-x};
      scan.SetBinContent(bin, fitter.EvaluateFCN(xp));
    }
    const double the_min(get_minimum(scan));
    for(int bin(1); bin<num_bins; ++bin){
      scan.SetBinContent(bin, scan.GetBinContent(bin)-the_min);
    }
    scan.SetMinimum(0.0);
    scan.SetMaximum(16.0);
    scan.Draw();
    c.Print("crap2.pdf");
  }
}

void get_counts(std::ifstream& file, float count[3][2][4]){
  for(unsigned k(0); k<4; ++k){
    for(unsigned i(0); i<3; ++i){
      for(unsigned j(0); j<2; ++j){
        file >> count[i][j][k];
      }
    }
  }
}

void setup_histo(TH1& h, const float count[3][2][4]){
  h.Sumw2();
  unsigned bin(1);
  for(unsigned k(0); k<4; ++k){
    for(unsigned i(0); i<3; ++i){
      for(unsigned j(0); j<2; ++j, ++bin){
        for(unsigned fill(0); fill<count[i][j][k]; ++fill){
          h.Fill(bin);
        }
      }
    }
  }
}
