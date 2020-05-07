#ifndef BranchManager_h
#define BranchManager_h

#include <TBranch.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <string>
#include <unordered_map>
#include <vector>

struct BranchManager {
    std::vector<TBranch*> branchHolder;
    std::unordered_map<std::string, TBranch*> specificBranch;
    TTree* fChain;

    void SetTree(TTree* fChain_) { fChain = fChain_; }

    template <typename T>
    void SetBranch(std::string name, T& holder) {
        branchHolder.push_back({});
        fChain->SetBranchAddress(name.c_str(), &holder, &branchHolder.back());
    }

    template <typename T>
    void SetSpecificBranch(std::string name, T& holder) {
        if (specificBranch.find(name) != specificBranch.end())
            specificBranch[name] = {};
        fChain->SetBranchAddress(name.c_str(), &holder, &specificBranch[name]);
    }

    void SetEntry(int entry) {
        for (auto& it : branchHolder) {
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

    void CleanSpecificBranch(std::string name) {
        if (specificBranch.find(name) != specificBranch.end())
            specificBranch.at(name) = nullptr;
    }
};

#endif
