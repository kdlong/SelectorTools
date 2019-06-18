#include <vector>
#include <string>
#include <unordered_map>
#include <TBranch.h>
#include <TTree.h>


struct BranchManager {
  std::vector<TBranch*> branchHolder;
  std::unordered_map<std::string, TBranch*> specificBranch;
  TTree* fChain;

  void SetTree(TTree* fChain_) {
    fChain = fChain_;
  } 

  template<typename T>
  void SetBranch(std::string name, T& holder) {
    branchHolder.push_back(new TBranch());
    fChain->SetBranchAddress(name.c_str(), &holder, &branchHolder.back());
  }

  template<typename T>
  void SetSpecificBranch(std::string name, T& holder) {
    specificBranch[name] = new TBranch();
    fChain->SetBranchAddress(name.c_str(), &holder, &specificBranch[name]);
  }

  void SetEntry(int entry) {
    for(auto& it: branchHolder) {
      it->GetEntry(entry);
    }
  }

  void SetSpecificEntry(int entry, std::string name) {
    specificBranch[name]->GetEntry(entry);
  }
  
  void CleanUp() {
    branchHolder.clear();
    specificBranch.clear();
  }
};



