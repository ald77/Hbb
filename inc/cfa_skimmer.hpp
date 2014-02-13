#ifndef H_CFA_SKIMMER
#define H_CFA_SKIMMER

#include <string>
#include "event_handler.hpp"

class CfASkimmer : public EventHandler{
public:
  CfASkimmer(const std::string& in_file_name,
             const bool is_list,
             const double weight_in=1.0);

  void Skim(const std::string& out_file_name,
            const int chargino_mass=-1,
            const int LSP_mass=-1);
};

#endif
