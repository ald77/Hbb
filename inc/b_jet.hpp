#ifndef H_BJET
#define H_BJET

#include <cfloat>
#include "TLorentzVector.h"

class BJet{
public:
  BJet(const TLorentzVector=TLorentzVector(0.0,0.0,0.0,0.0), const double bTagIn=-DBL_MAX, const unsigned int jets_AK5PF_indexIn=0);

  void SetLorentzVector(const TLorentzVector vecIn);
  void SetBTag(const double bTagIn);
  void SetIndex(const double jets_AK5PF_indexIn);

  TLorentzVector GetLorentzVector() const;
  double GetBTag() const;
  unsigned int GetIndex() const;

  bool operator==(const BJet &jet) const;
  bool operator!=(const BJet &jet) const;
  bool operator<(const BJet &jet) const;
  bool operator>(const BJet &jet) const;
  bool operator<=(const BJet &jet) const;
  bool operator>=(const BJet &jet) const;
private:
  TLorentzVector vec;
  double bTag;
  int jets_AK5PF_index;
};

#endif
