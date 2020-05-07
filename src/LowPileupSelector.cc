#include "Analysis/VVAnalysis/interface/LowPileupSelector.h"
#include <TStyle.h>
#include <regex>

void LowPileupSelector::Init(TTree *tree) {
    b.SetTree(tree);
    if (name_.find("we") != std::string::npos ||
        name_.find("wm") != std::string::npos ||
        name_.find("wm") != std::string::npos) {
        isW_ = true;
    } else if (name_.find("DYm50") != std::string::npos) {
        isZ_ = true;
    }
    SelectorBase::Init(tree);
}

void LowPileupSelector::SetBranchesBacon() {
    b.CleanUp();
    b.SetBranch("genVPt", genVPt);
    b.SetBranch("genVPhi", genVPhi);
    b.SetBranch("genVy", genVy);
    b.SetBranch("genVMass", genVMass);
    b.SetBranch("scale1fb", scale1fb);
}

void LowPileupSelector::LoadBranchesBacon(
    Long64_t entry, std::pair<Systematic, std::string> variation) {
    weight = 1;
    b.SetEntry(entry);

    if (isMC_) {
        weight = scale1fb * 1000;
        // weight = genWeight*PUWeight*scale1fb;
    }
    SetComposite();
}

void LowPileupSelector::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();
    if (name_.find("__m") != std::string::npos ||
        name_.find("__tm") != std::string::npos) {
        allChannels_ = {{mp, "mp"}, {mn, "mn"}};
        isE_ = false;
    } else if (name_.find("__e") != std::string::npos ||
               name_.find("__te") != std::string::npos) {
        allChannels_ = {{ep, "ep"}, {en, "en"}};
        isE_ = true;
    }

    InitializeHistogramsFromConfig();
}
