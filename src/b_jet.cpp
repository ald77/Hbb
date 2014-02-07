#include <utility>
#include "b_jet.hpp"
#include "TLorentzVector.h"

BJet::BJet(const TLorentzVector vecIn, const double bTagIn, const int jets_AK5PF_indexIn, const int pIdIn, const int momIdIn):vec(vecIn),bTag(bTagIn),jets_AK5PF_index(jets_AK5PF_indexIn),pId(pIdIn),momId(momIdIn){
}

void BJet::SetLorentzVector(const TLorentzVector vecIn){
  vec=vecIn;
}

void BJet::SetBTag(const double bTagIn){
  bTag=bTagIn;
}

void BJet::SetPartonId(const int pIdIn){
  pId=pIdIn;
}

void BJet::SetMotherId(const int momIdIn){
  momId=momIdIn;
}

TLorentzVector BJet::GetLorentzVector() const{
  return vec;
}

double BJet::GetBTag() const{
  return bTag;
}

unsigned int BJet::GetIndex() const{
  return jets_AK5PF_index;
}

int BJet::GetPartonId() const{
  return pId;
}

int BJet::GetMotherId() const{
  return momId;
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
