#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"

int main(int argc, char *argv[]){
  for(int arg(1); arg<argc; ++arg){
    std::cout << "Repairing " << argv[arg] << "\n";
    const std::string subDir("configurableAnalysis");
    const std::string nameA("eventA");
    const std::string nameB("eventB");
    TFile file(argv[arg],"update");
    if(file.IsOpen() && file.IsWritable() && !file.IsZombie()){
      TDirectory *dir=file.GetDirectory(subDir.c_str());
      if(dir==NULL){
	std::clog << "Creating new configurableAnalysis TDirectory.\n";
	dir=file.mkdir(subDir.c_str(),subDir.c_str());
      }else{
	std::clog << "configurableAnalysis TDirectory already exists.\n";
      }
      if(dir!=NULL && dir->IsWritable() && !dir->IsZombie()){
	TKey *keyA(dir->GetKey(nameA.c_str())==NULL?file.GetKey(nameA.c_str()):NULL);
	TKey *keyB(dir->GetKey(nameB.c_str())==NULL?file.GetKey(nameB.c_str()):NULL);
	TObject *objectA(keyA==NULL?NULL:keyA->ReadObj());
	TObject *objectB(keyB==NULL?NULL:keyB->ReadObj());
	TTree *treeA(objectA==NULL?NULL:static_cast<TTree*>(objectA));
	TTree *treeB(objectB==NULL?NULL:static_cast<TTree*>(objectB));
	if(treeA==NULL){
	  std::clog << "Could not find " << nameA << " in root directory of " << argv[arg] << ".\n";
	}else{
	  std::clog << "Cloning " << nameA << ".\n";
	}
	const TTree *eventA(treeA==NULL?NULL:treeA->CloneTree(-1,"fast"));
	if(treeB==NULL){
	  std::clog << "Could not find " << nameB << " in root directory of " << argv[arg] << ".\n";
	}else{
	  std::clog << "Cloning " << nameB << ".\n";
	}
	const TTree *eventB(treeB==NULL?NULL:treeB->CloneTree(-1,"fast"));
	if(dir->cd()){
	  if(eventA!=NULL){
	    dir->Delete((nameA+";*").c_str());
	    if(dir!=NULL && dir->IsWritable() && !dir->IsZombie() && dir->cd()){
	      std::clog << "Writing " << nameA << ".\n";
	      eventA->Write(nameA.c_str(),TObject::kOverwrite);
	    }else{
	      std::cerr << "Error in " << argv[0] << ": configurableAnalysis not writeable for " << nameA << ".\n";
	    }
	  }else{
	    std::clog << "Skipping writing of " << nameA << ".\n";
	  }
	  if(eventB!=NULL){
	    dir->Delete((nameB+";*").c_str());
	    if(dir!=NULL && dir->IsWritable() && !dir->IsZombie() && dir->cd()){
	      std::clog << "Writing " << nameB << ".\n";
	      eventB->Write(nameB.c_str(),TObject::kOverwrite);
	    }else{
	      std::cerr << "Error in " << argv[0] << ": configurableAnalysis not writeable for " << nameB << ".\n";
	    }
	  }else{
	    std::clog << "Skipping writing of " << nameB << ".\n";
	  }
	}else{
	  std::cerr << "Error in " << argv[0] << ": could not cd into configureableAnalysis.\n";
	}
	if(file.cd()){
	  file.cd();
	  file.Delete((nameA+";*").c_str());
	  file.Delete((nameB+";*").c_str());
	  file.Close();
	  file.Open(argv[arg],"update");
	}else{
	  std::cerr << "Error in " << argv[0] << ": could not cd into root directory of " << argv[arg] << ".\n";
	}
      }
    }else{
      std::cerr << "Error in " << argv[0] << ": could not open " << argv[arg] << ".\n";
    }
    if(file.IsOpen()){
      file.Close();
    }
  }
}
