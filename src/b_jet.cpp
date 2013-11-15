#include <utility>
#include "b_jet.hpp"
#include "TLorentzVector.h"

BJet::BJet(const TLorentzVector vecIn, const double bTagIn):vec(vecIn),bTag(bTagIn){
}

void BJet::SetLorentzVector(const TLorentzVector vecIn){
  vec=vecIn;
}

void BJet::SetBTag(const double bTagIn){
  bTag=bTagIn;
}

TLorentzVector BJet::GetLorentzVector() const{
  return vec;
}

double BJet::GetBTag() const{
  return bTag;
}

bool BJet::operator==(const BJet &jet) const{
  return vec==jet.vec && bTag==jet.bTag;
}

bool BJet::operator!=(const BJet &jet) const{
  return !(*this==jet);
}

bool BJet::operator<(const BJet &jet) const{
  std::pair<double, double> this_one(bTag, -vec.Pt());
  std::pair<double, double> that_one(jet.bTag, -jet.vec.Pt());
  return this_one<that_one;
}

bool BJet::operator>(const BJet &jet) const{
  std::pair<double, double> this_one(bTag, -vec.Pt());
  std::pair<double, double> that_one(jet.bTag, -jet.vec.Pt());
  return this_one>that_one;
}

bool BJet::operator<=(const BJet &jet) const{
  return *this==jet || *this<jet;
}

bool BJet::operator>=(const BJet &jet) const{
  return *this==jet || *this>jet;
}
