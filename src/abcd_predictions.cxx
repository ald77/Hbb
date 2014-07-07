//Produces closure test plots from reduced trees

#include "abcd_predictions.hpp"
#include <iostream>
#include <string>
#include "TLegend.h"
#include "TChain.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "utils.hpp"
#include "style.hpp"
#include "timer.hpp"
#include "math.hpp"

void print_good(double x[3][2][4], double u[3][2][4]){
  for(unsigned sbin_index(0); sbin_index<4; ++sbin_index){
    for(unsigned nb_index(0); nb_index<3; ++nb_index){
      for(unsigned mbb_index(0); mbb_index<2; ++mbb_index){
	std::cout << "$" << fix_width(x[nb_index][mbb_index][sbin_index],4) << "\\pm" << fix_width(u[nb_index][mbb_index][sbin_index],4) << "$ & ";
      }
    }
    std::cout << "\\\\" << std::endl;
  }
}

int main(){
  SetStyle();
  const std::string cut_baseline("passesJSONCut && passesPVCut && passesJet2PtCut && passes2CSVTCut && passesMETSig30Cut && passesMETCleaningCut && passesTriggerCut && passesNumJetsCut && passesMinDeltaPhiCut && passesLeptonVetoCut && passesIsoTrackVetoCut && passesDRCut");
  const std::string cut_2b("num_b_tagged_jets==2");
  const std::string cut_3b("num_b_tagged_jets==3");
  const std::string cut_4b("num_b_tagged_jets>=4");
  const std::string cut_sig("higgs_mass_signal_region");
  const std::string cut_sb("higgs_mass_sideband");
  const std::string cut_sbin1("met_sig>=30.0 && met_sig<50.0");
  const std::string cut_sbin2("met_sig>=50.0 && met_sig<100.0");
  const std::string cut_sbin3("met_sig>=100.0 && met_sig<150.0");
  const std::string cut_sbin4("met_sig>=150.0");
  
  TChain chain("chain","chain");

  chain.Add("reduced_trees/QCD_HT*SyncSkim_mdpscaled.root/reduced_tree");
  /*chain.Add("reduced_trees/T*channel*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTH*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTJets_FullLept*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTJets_SemiLept*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTJets_Hadronic*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTW*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/TTZ*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/W*JetsToLNu_TuneZ2Star*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/WH_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/WW_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/WZ_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/ZH_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/ZZ_*SyncSkim.root/reduced_tree");
  chain.Add("reduced_trees/ZJetsToNuNu*SyncSkim.root/reduced_tree");*/

  //chain.Add("reduced_trees/MET_Run2012*SyncSkim.root/reduced_tree");

  double counts[3][2][4], uncerts[3][2][4];

  Timer timer(24);
  timer.Start();
  for(unsigned short nb_bin(0); nb_bin<3; ++nb_bin){
    std::string cut_nb("");
    switch(nb_bin){
    case 0:
      cut_nb=cut_2b;
      break;
    case 1:
      cut_nb=cut_3b;
      break;
    case 2:
      cut_nb=cut_4b;
      break;
    default:
      break;
    }
    for(unsigned short mbb_bin(0); mbb_bin<2; ++mbb_bin){
      std::string cut_mbb("");
      switch(mbb_bin){
      case 0:
	cut_mbb=cut_sig;
	break;
      case 1:
	cut_mbb=cut_sb;
	break;
      default:
      break;
      }
      for(unsigned short smet_bin(0); smet_bin<4; ++smet_bin){
	std::string cut_smet("");
	switch(smet_bin){
	case 0:
	  cut_smet=cut_sbin1;
	  break;
	case 1:
	  cut_smet=cut_sbin2;
	  break;
	case 2:
	  cut_smet=cut_sbin3;
	  break;
	case 3:
	  cut_smet=cut_sbin4;
	  break;
	default:
	  break;
	}

	const std::string cut_full("full_weight*("+cut_baseline+"&&"+cut_nb+"&&"+cut_mbb+"&&"+cut_smet+")");

	get_count_and_uncertainty(chain, cut_full, counts[nb_bin][mbb_bin][smet_bin], uncerts[nb_bin][mbb_bin][smet_bin]);
	timer.Iterate();
	timer.PrintRemainingTime();
      }
    }
  }
  double indep_counts[3][2][4];
  get_independence_model(counts, indep_counts);
  make_complicated_plot(counts, uncerts);
  print_counts(counts, uncerts);
  print_good(counts, uncerts);
  std::cout << "Independence model:" << std::endl;
  print_matrix(indep_counts);
  print_kappas(counts, uncerts);
}

void sum_on_sbins(double counts_in[3][2][4], double uncerts_in[3][2][4], double counts_out[3][2], double uncerts_out[3][2]){
  for(unsigned nb_index(0); nb_index<3; ++nb_index){
    for(unsigned mbb_index(0); mbb_index<2; ++mbb_index){
      counts_out[nb_index][mbb_index]=0.0;
      uncerts_out[nb_index][mbb_index]=0.0;
      for(unsigned sbin_index(0); sbin_index<4; ++sbin_index){
	counts_out[nb_index][mbb_index]+=counts_in[nb_index][mbb_index][sbin_index];
	uncerts_out[nb_index][mbb_index]+=uncerts_in[nb_index][mbb_index][sbin_index]*uncerts_in[nb_index][mbb_index][sbin_index];
      }
      uncerts_out[nb_index][mbb_index]=sqrt(uncerts_out[nb_index][mbb_index]);
    }
  }
}

void get_kappa_and_uncert(double& kappa, double& uncert, const double& a, const double& b, const double& c, const double& d, const double& ua, const double& ub, const double& uc, const double &ud){
  kappa=a*d/(b*c);
  uncert=sqrt(a*a*b*b*c*c*ud*ud+a*a*b*b*uc*uc*d*d+a*a*ub*ub*c*c*d*d+ua*ua*b*b*c*c*d*d)/(b*b*c*c);
}

void print_kappas(double counts[3][2][4], double uncerts[3][2][4]){
  double counts_sum[3][2], uncerts_sum[3][2];
  sum_on_sbins(counts, uncerts, counts_sum, uncerts_sum);
  double kappa23(0.0), kappa24(0.0), kappa34(0.0);
  double uncert23(0.0), uncert24(0.0), uncert34(0.0);
  get_kappa_and_uncert(kappa23, uncert23, counts_sum[1][0], counts_sum[1][1], counts_sum[0][0], counts_sum[0][1], uncerts_sum[1][0], uncerts_sum[1][1], uncerts_sum[0][0], uncerts_sum[0][1]);
  get_kappa_and_uncert(kappa24, uncert24, counts_sum[2][0], counts_sum[2][1], counts_sum[0][0], counts_sum[0][1], uncerts_sum[2][0], uncerts_sum[2][1], uncerts_sum[0][0], uncerts_sum[0][1]);
  get_kappa_and_uncert(kappa34, uncert34, counts_sum[2][0], counts_sum[2][1], counts_sum[1][0], counts_sum[1][1], uncerts_sum[2][0], uncerts_sum[2][1], uncerts_sum[1][0], uncerts_sum[1][1]);
  std::cout << "Overall: " << std::endl;
  std::cout << kappa23 << "+-" << uncert23 << std::endl;
  std::cout << kappa24 << "+-" << uncert24 << std::endl;
  std::cout << kappa34 << "+-" << uncert34 << std::endl;
  
  for(unsigned sbin_index(0); sbin_index<4; ++sbin_index){
    get_kappa_and_uncert(kappa23, uncert23, counts[1][0][sbin_index], counts[1][1][sbin_index], counts[0][0][sbin_index], counts[0][1][sbin_index], uncerts[1][0][sbin_index], uncerts[1][1][sbin_index], uncerts[0][0][sbin_index], uncerts[0][1][sbin_index]);
    get_kappa_and_uncert(kappa24, uncert24, counts[2][0][sbin_index], counts[2][1][sbin_index], counts[0][0][sbin_index], counts[0][1][sbin_index], uncerts[2][0][sbin_index], uncerts[2][1][sbin_index], uncerts[0][0][sbin_index], uncerts[0][1][sbin_index]);
    get_kappa_and_uncert(kappa34, uncert34, counts[2][0][sbin_index], counts[2][1][sbin_index], counts[1][0][sbin_index], counts[1][1][sbin_index], uncerts[2][0][sbin_index], uncerts[2][1][sbin_index], uncerts[1][0][sbin_index], uncerts[1][1][sbin_index]);
    std::cout << "S-bin: " << sbin_index+1 << std::endl;
    std::cout << kappa23 << "+-" << uncert23 << std::endl;
    std::cout << kappa24 << "+-" << uncert24 << std::endl;
    std::cout << kappa34 << "+-" << uncert34 << std::endl;
  }
}

void get_abcd_prediction_and_uncertainty(double num_count_1, double num_uncert_1,
					 double num_count_2, double num_uncert_2,
					 double den_count, double den_uncert,
					 double& count, double& uncert){
  count=num_count_1*num_count_2/(den_count);
  uncert=std::sqrt(num_count_1*num_count_1*num_count_2*num_count_2*den_uncert*den_uncert+num_count_1*num_count_1*num_uncert_2*num_uncert_2*den_count*den_count+num_uncert_1*num_uncert_1*num_count_2*num_count_2*den_count*den_count)/(den_count*den_count);
}

void get_independence_model(double in[3][2][4], double out[3][2][4]){
  double out1[3], out2[2], out3[4];
  std::vector<double> sum_vec(0);
  for(unsigned index1(0); index1<3; ++index1){
    out1[index1]=0.0;
    sum_vec.clear();
    for(unsigned index2(0); index2<2; ++index2){
      for(unsigned index3(0); index3<4; ++index3){
	sum_vec.push_back(in[index1][index2][index3]);
      }
    }
    out1[index1]=Math::Sum(sum_vec.begin(), sum_vec.end());
  }
  for(unsigned index2(0); index2<2; ++index2){
    out2[index2]=0.0;
    sum_vec.clear();
    for(unsigned index1(0); index1<3; ++index1){
      for(unsigned index3(0); index3<4; ++index3){
	sum_vec.push_back(in[index1][index2][index3]);
      }
    }
    out2[index2]=Math::Sum(sum_vec.begin(), sum_vec.end());
  }
  for(unsigned index3(0); index3<4; ++index3){
    out3[index3]=0.0;
    sum_vec.clear();
    for(unsigned index1(0); index1<3; ++index1){
      for(unsigned index2(0); index2<2; ++index2){
	sum_vec.push_back(in[index1][index2][index3]);
      }
    }
    out3[index3]=Math::Sum(sum_vec.begin(), sum_vec.end());
  }
  sum_vec.clear();
  for(unsigned index3(0); index3<4; ++index3){
    for(unsigned index1(0); index1<3; ++index1){
      for(unsigned index2(0); index2<2; ++index2){
	sum_vec.push_back(in[index1][index2][index3]);
      }
    }
  }
  double sum_sq_inv(Math::Sum(sum_vec.begin(), sum_vec.end()));
  sum_sq_inv=1.0/(sum_sq_inv*sum_sq_inv);
  for(unsigned index3(0); index3<4; ++index3){
    for(unsigned index1(0); index1<3; ++index1){
      for(unsigned index2(0); index2<2; ++index2){
	out[index1][index2][index3]=out1[index1]*out2[index2]*out3[index3]*sum_sq_inv;
      }
    }
  }
}

void print_counts(double counts[3][2][4], double uncerts[3][2][4]){
  print_matrix(counts);
  std::cout << std::endl;
  print_matrix(uncerts);
}

void print_matrix(double mat[3][2][4]){
  for(unsigned sbin_index(0); sbin_index<4; ++sbin_index){
    for(unsigned nb_index(0); nb_index<3; ++nb_index){
      for(unsigned mbb_index(0); mbb_index<2; ++mbb_index){
	std::cout << fix_width(mat[nb_index][mbb_index][sbin_index],4) << ' ';
      }
    }
    std::cout << std::endl;
  }
}

void make_complicated_plot(double counts[3][2][4], double uncerts[3][2][4]){
  double summed_counts[3][2], summed_uncerts[3][2];
  for(unsigned short b_index(0); b_index<3; ++b_index){
    for(unsigned short mbb_index(0); mbb_index<2; ++mbb_index){
      for(unsigned short smet_index(0); smet_index<4; ++smet_index){
	summed_counts[b_index][mbb_index]+=counts[b_index][mbb_index][smet_index];
	summed_uncerts[b_index][mbb_index]+=uncerts[b_index][mbb_index][smet_index]*uncerts[b_index][mbb_index][smet_index];
      }
      summed_uncerts[b_index][mbb_index]=sqrt(summed_uncerts[b_index][mbb_index]);
    }
  }

  TH1D h_closure_test("h_closure_test", "Closure Test;Sample;Events/19.4 fb^{-1}", 15, 0.5, 15.5);
  h_closure_test.GetXaxis()->SetBinLabel(1,"4b,all S");
  h_closure_test.GetXaxis()->SetBinLabel(2,"3b,all S");
  h_closure_test.GetXaxis()->SetBinLabel(3,"2b,all S");
  h_closure_test.GetXaxis()->SetBinLabel(4,"4b,S-bin 1");
  h_closure_test.GetXaxis()->SetBinLabel(5,"3b,S-bin 1");
  h_closure_test.GetXaxis()->SetBinLabel(6,"2b,S-bin 1");
  h_closure_test.GetXaxis()->SetBinLabel(7,"4b,S-bin 2");
  h_closure_test.GetXaxis()->SetBinLabel(8,"3b,S-bin 2");
  h_closure_test.GetXaxis()->SetBinLabel(9,"2b,S-bin 2");
  h_closure_test.GetXaxis()->SetBinLabel(10,"4b,S-bin 3");
  h_closure_test.GetXaxis()->SetBinLabel(11,"3b,S-bin 3");
  h_closure_test.GetXaxis()->SetBinLabel(12,"2b,S-bin 3");
  h_closure_test.GetXaxis()->SetBinLabel(13,"4b,S-bin 4");
  h_closure_test.GetXaxis()->SetBinLabel(14,"3b,S-bin 4");
  h_closure_test.GetXaxis()->SetBinLabel(15,"2b,S-bin 4");

  h_closure_test.SetBinContent(1, summed_counts[2][0]);
  h_closure_test.SetBinContent(2, summed_counts[1][0]);
  h_closure_test.SetBinContent(3, summed_counts[0][0]);
  h_closure_test.SetBinError(1, summed_uncerts[2][0]);
  h_closure_test.SetBinError(2, summed_uncerts[1][0]);
  h_closure_test.SetBinError(3, summed_uncerts[0][0]);

  h_closure_test.SetFillStyle(3003);
  h_closure_test.SetFillColor(1);
  h_closure_test.SetLineColor(1);
  h_closure_test.SetStats(0);

  for(unsigned short b_index(0); b_index<3; ++b_index){
    for(unsigned short smet_index(0); smet_index<4; ++smet_index){
      const unsigned short bin(3*smet_index+6-b_index);
      h_closure_test.SetBinContent(bin, counts[b_index][0][smet_index]);
      h_closure_test.SetBinError(bin, uncerts[b_index][0][smet_index]);
    }
  }

  const double shift(1.0/8.0);
  const double x_2b[10]={1.0-shift, 2.0-shift, 4.0-shift, 5.0-shift, 7.0-shift, 8.0-shift, 10.0-shift, 11.0-shift, 13.0-shift, 14.0-shift};
  const double x_3b[10]={1.0+shift, 3.0-shift, 4.0+shift, 6.0-shift, 7.0+shift, 9.0-shift, 10.0+shift, 12.0-shift, 13.0+shift, 15.0-shift};
  const double x_4b[10]={2.0+shift, 3.0+shift, 5.0+shift, 6.0+shift, 8.0+shift, 9.0+shift, 11.0+shift, 12.0+shift, 14.0+shift, 15.0+shift};
  const double zeroes[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double y_2b[10], y_3b[10], y_4b[10], uncert_2b[10], uncert_3b[10], uncert_4b[10];
  y_2b[0]=summed_counts[0][0]*summed_counts[2][1]/summed_counts[0][1];
  y_2b[1]=summed_counts[0][0]*summed_counts[1][1]/summed_counts[0][1];
  y_3b[0]=summed_counts[1][0]*summed_counts[2][1]/summed_counts[1][1];
  y_3b[1]=summed_counts[1][0]*summed_counts[0][1]/summed_counts[1][1];
  y_4b[0]=summed_counts[2][0]*summed_counts[1][1]/summed_counts[2][1];
  y_4b[1]=summed_counts[2][0]*summed_counts[0][1]/summed_counts[2][1];
  get_abcd_prediction_and_uncertainty(summed_counts[0][0], summed_uncerts[0][0],
				      summed_counts[2][1], summed_uncerts[2][1],
				      summed_counts[0][1], summed_uncerts[0][1],
				      y_2b[0], uncert_2b[0]);
  get_abcd_prediction_and_uncertainty(summed_counts[0][0], summed_uncerts[0][0],
				      summed_counts[1][1], summed_uncerts[1][1],
				      summed_counts[0][1], summed_uncerts[0][1],
				      y_2b[1], uncert_2b[1]);
  get_abcd_prediction_and_uncertainty(summed_counts[1][0], summed_uncerts[1][0],
				      summed_counts[2][1], summed_uncerts[2][1],
				      summed_counts[1][1], summed_uncerts[1][1],
				      y_3b[0], uncert_3b[0]);
  get_abcd_prediction_and_uncertainty(summed_counts[1][0], summed_uncerts[1][0],
				      summed_counts[0][1], summed_uncerts[0][1],
				      summed_counts[1][1], summed_uncerts[1][1],
				      y_3b[1], uncert_3b[1]);
  get_abcd_prediction_and_uncertainty(summed_counts[2][0], summed_uncerts[2][0],
				      summed_counts[1][1], summed_uncerts[1][1],
				      summed_counts[2][1], summed_uncerts[2][1],
				      y_4b[0], uncert_4b[0]);
  get_abcd_prediction_and_uncertainty(summed_counts[2][0], summed_uncerts[2][0],
				      summed_counts[0][1], summed_uncerts[0][1],
				      summed_counts[2][1], summed_uncerts[2][1],
				      y_4b[1], uncert_4b[1]);
  for(unsigned short bin(2); bin<10; ++bin){
    const unsigned short b_bin_2b((bin&1u)?1:2);
    const unsigned short b_bin_3b((bin&1u)?0:2);
    const unsigned short b_bin_4b((bin&1u)?0:1);
    unsigned short smet_bin(0);
    if(bin<4){
      smet_bin=0;
    }else if(bin<6){
      smet_bin=1;
    }else if(bin<8){
      smet_bin=2;
    }else{
      smet_bin=3;
    }
    get_abcd_prediction_and_uncertainty(counts[0][0][smet_bin], uncerts[0][0][smet_bin], counts[b_bin_2b][1][smet_bin], uncerts[b_bin_2b][1][smet_bin], counts[0][1][smet_bin], uncerts[0][1][smet_bin], y_2b[bin], uncert_2b[bin]);
    get_abcd_prediction_and_uncertainty(counts[1][0][smet_bin], uncerts[1][0][smet_bin], counts[b_bin_3b][1][smet_bin], uncerts[b_bin_3b][1][smet_bin], counts[1][1][smet_bin], uncerts[1][1][smet_bin], y_3b[bin], uncert_3b[bin]);
    get_abcd_prediction_and_uncertainty(counts[2][0][smet_bin], uncerts[2][0][smet_bin], counts[b_bin_4b][1][smet_bin], uncerts[b_bin_4b][1][smet_bin], counts[2][1][smet_bin], uncerts[2][1][smet_bin], y_4b[bin], uncert_4b[bin]);
  }
  TGraphErrors g_2b(10, x_2b, y_2b, zeroes, uncert_2b);
  TGraphErrors g_3b(10, x_3b, y_3b, zeroes, uncert_3b);
  TGraphErrors g_4b(10, x_4b, y_4b, zeroes, uncert_4b);
  g_2b.SetMarkerStyle(20); g_2b.SetMarkerColor(2); g_2b.SetLineColor(2);
  g_3b.SetMarkerStyle(20); g_3b.SetMarkerColor(3); g_3b.SetLineColor(3);
  g_4b.SetMarkerStyle(20); g_4b.SetMarkerColor(4); g_4b.SetLineColor(4);

  double max(get_maximum(h_closure_test));
  if(get_maximum(g_2b)>max) max=get_maximum(g_2b);
  if(get_maximum(g_3b)>max) max=get_maximum(g_3b);
  if(get_maximum(g_4b)>max) max=get_maximum(g_4b);
  double min(get_minimum_positive(h_closure_test));
  if(get_minimum_positive(g_2b)<min) min=get_minimum_positive(g_2b);
  if(get_minimum_positive(g_3b)<min) min=get_minimum_positive(g_3b);
  if(get_minimum_positive(g_4b)<min) min=get_minimum_positive(g_4b);
  h_closure_test.SetMinimum(0.9*min); h_closure_test.SetMaximum(1.1*max);
  g_2b.SetMinimum(0.9*min); g_2b.SetMaximum(1.1*max);
  g_3b.SetMinimum(0.9*min); g_3b.SetMaximum(1.1*max);
  g_4b.SetMinimum(0.9*min); g_4b.SetMaximum(1.1*max);

  TCanvas canvas;
  h_closure_test.Draw("e2");
  h_closure_test.Draw("e1same");
  g_2b.Draw("psame");
  g_3b.Draw("psame");
  g_4b.Draw("psame");
  TLegend legend(0.8, 0.85, 1.0, 1.0);
  legend.AddEntry(&h_closure_test, "Observed", "flep");
  legend.AddEntry(&g_2b, "Pred. from 2b", "lep");
  legend.AddEntry(&g_3b, "Pred. from 3b", "lep");
  legend.AddEntry(&g_4b, "Pred. from 4b", "lep");
  legend.Draw("same");
  canvas.SetLogy(1);
  canvas.Print("closure_test.pdf");
}

double add_in_quadrature(const double a, const double b, const double c){
  double x(0.0), y(0.0), z(0.0);
  if(a>=b && a>=c){
    x=a; y=b; z=c;
  }else if(b>=a && b>=c){
    x=b; y=a; z=c;
  }else{
    x=c; y=a; z=b;
  }
  if(x==0.0){
    return 0.0;
  }else{
    const double r1(y/x), r2(z/x);
    return fabs(x)*sqrt(1.0+r1*r1+r2*r2);
  }
}

double add_in_quadrature(const double a, const double b,
			 const double c, const double d){
  double w(0.0), x(0.0), y(0.0), z(0.0);
  if(a>=b && a>=c && a>=d){
    w=a; x=b; y=c; z=d;
  }else if(b>=a && b>=c && b>=d){
    w=b; x=a; y=c; z=d;
  }else if(c>=a && c>=b && c>=d){
    w=c; x=a; y=b; z=d;
  }else{
    w=d; x=a; y=b; z=c;
  }
  if(w==0.0){
    return 0.0;
  }else{
    const double r1(x/w), r2(y/w), r3(z/w);
    return fabs(w)*sqrt(1.0+r1*r1+r2*r2+r3*r3);
  }
}

