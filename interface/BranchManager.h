#include <vector>
#include <string>
#include <TBranch.h>
#include <TTree.h>

struct BranchManager {
  std::vector<TBranch*> branchHolder;
  TTree* fChain;

  void SetTree(TTree* fChain_) {
    fChain = fChain_;
  } 
  
  template<typename T>
  void SetBranch(std::string name, T& holder) {
    branchHolder.push_back(new TBranch());
    fChain->SetBranchAddress(name.c_str(), &holder, &branchHolder.back());
  }
  
  void SetEntry(int entry) {
    for(auto& it: branchHolder) {
      it->GetEntry(entry);
    }
  }

  void CleanUp() {
    branchHolder.clear();
  }
};



