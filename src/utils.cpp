#include "utils.hpp"

#include "TGraph.h"

void add_point(TGraph& graph, const double x, const double y){
  graph.SetPoint(graph.GetN(), x, y);
}
