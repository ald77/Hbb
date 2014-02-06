#ifndef H_BJET
#define H_BJET

#include <cfloat>
#include "TLorentzVector.h"

class BJet{
public:
  explicit BJet(const TLorentzVector=TLorentzVector(0.0,0.0,0.0,0.0), const double bTagIn=-DBL_MAX, const int jets_AK5PF_indexIn=0, const int pIdIn=0, const int momIdIn=0);

  void SetLorentzVector(const TLorentzVector vecIn);
  void SetBTag(const double bTagIn);
  void SetPartonId(const int pIdIn);
  void SetMotherId(const int momIdIn);
  void SetIndex(const double jets_AK5PF_indexIn);

  TLorentzVector GetLorentzVector() const;
  double GetBTag() const;
  unsigned int GetIndex() const;
  int GetPartonId() const;
  int GetMotherId() const;

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
  int pId;
  int momId;
};

#endif
