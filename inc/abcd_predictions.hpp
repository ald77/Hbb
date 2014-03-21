#ifndef H_ABCD_PREDICTIONS
#define H_ABCD_PREDICTIONS

#include <cmath>

void get_abcd_prediction_and_uncertainty(double num_count_1, double num_uncert_1,
					 double num_count_2, double num_uncert_2,
					 double den_count, double den_uncert,
					 double& count, double& uncert);

void make_complicated_plot(double counts[3][2][4], double uncerts[3][2][4]);
void make_simple_plot(double counts[3][2][4], double uncerts[3][2][4]);
double add_in_quadrature(const double a, const double b, const double c);
double add_in_quadrature(const double a, const double b,
			 const double c, const double d);

#endif
