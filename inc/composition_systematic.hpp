#ifndef H_COMPOSITION_SYSTEMATIC
#define H_COMPOSITION_SYSTEMATIC

#include <string>
#include "TH1D.h"
#include "TLegend.h"

void get_total(double x[3][2][5]);
void set_content(TH1D& histo,
		 const double qcd[3][2][5],
		 const double other_background[3][2][5],
		 const unsigned short kappa_type,
		 const unsigned short sbin);
void plot(TLegend& l, TH1D& h1, TH1D& h2, TH1D& h3, const std::string& name);

#endif
