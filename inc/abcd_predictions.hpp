#ifndef H_ABCD_PREDICTIONS
#define H_ABCD_PREDICTIONS

#include <cmath>

void get_abcd_prediction_and_uncertainty(double num_count_1, double num_uncert_1,
					 double num_count_2, double num_uncert_2,
					 double den_count, double den_uncert,
					 double& count, double& uncert);

void sub_on_sbins(double counts_in[3][2][4], double uncerts_in[3][2][4], double counts_out[3][2], double uncerts_out[3][2]);
void print_kappas(double counts[3][2][4], double uncerts[3][2][4]);
void get_kappa_and_uncert(double& kappa, double& uncert, const double& a, const double& b, const double& c, const double& d, const double& ua, const double& ub, const double& uc, const double &ud);
void make_complicated_plot(double counts[3][2][4], double uncerts[3][2][4]);
void get_independence_model(double in[3][2][4], double out[3][2][4]);
void print_counts(double counts[3][2][4], double uncerts[3][2][4]);
void print_matrix(double mat[3][2][4]);
void make_simple_plot(double counts[3][2][4], double uncerts[3][2][4]);
double add_in_quadrature(const double a, const double b, const double c);
double add_in_quadrature(const double a, const double b,
			 const double c, const double d);

#endif
