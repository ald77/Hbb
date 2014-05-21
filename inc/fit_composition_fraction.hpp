#ifndef H_FIT_COMPOSITION_SYSTEMATIC
#define H_FIT_COMPOSITION_SYSTEMATIC

#include <fstream>
#include "TH1.h"

void get_counts(std::ifstream& file, float count[3][2][4]);
void setup_histo(TH1& h, const float count[3][2][4]);

#endif
