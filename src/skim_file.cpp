//Skims MC backgrounds using the full path to a local .root file (/path/to/local/file.root) or cfA sample name (TTJets_blahblahblah)

#include <string>
#include <unistd.h>
#include "event_handler.hpp"
#include "weights.hpp"

int main(int argc, char *argv[]){
  std::string inFilename("");
  bool isLocal(false);
  int c(0);
  while((c=getopt(argc, argv, "i:l"))!=-1){
    switch(c){
    case 'i':
      inFilename=optarg;
      break;
    case 'l':
      isLocal=true;
      break;
    }
  }

  std::string outFilename("");
  if(!isLocal){
    outFilename="../data/"+inFilename+"_SyncSkim.root";
    inFilename="/net/cms2/cms2r0/cfA/"+inFilename+"/cfA_"+inFilename+"*.root";
  }else{
    std::string baseName(inFilename);
    size_t pos=baseName.find(".root");
    if(pos!=std::string::npos){
      baseName.erase(pos);
    }
    pos=baseName.rfind("/");
    if(pos!=std::string::npos){
      baseName.erase(0,pos+1);
    }
    outFilename="../data/"+baseName+"_SyncSkim.root";
  }

  WeightCalculator w(19399);
  EventHandler eH(inFilename, false, w.GetWeight(inFilename));
  eH.Skim(outFilename);
}
