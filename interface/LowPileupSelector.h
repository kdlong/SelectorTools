#ifndef LowPileupSelector_h
#define LowPileupSelector_h

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TSelector.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/VVAnalysis/interface/BranchManager.h"
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/SelectorBase.h"

#include "Analysis/VVAnalysis/interface/helpers.h"

class LowPileupSelector : public SelectorBase {
   public:
    // Derived values
    float weight;
    TLorentzVector pfMetVec;
    bool isW_ = false;
    bool isZ_ = false;
    bool isE_;

    Float_t pfMet;
    Float_t pfMetPhi;
    Float_t genVPt;
    Float_t genVPhi;
    Float_t genVy;
    Float_t genVMass;
    Float_t scale1fb;

    BranchManager b;

    // Readers to access the data (delete the ones you do not need).
    virtual void Init(TTree *tree) override;
    LowPileupSelector(TTree * /*tree*/ = 0) {}
    ~LowPileupSelector() {}
    virtual void SetupNewDirectory() override;

    ClassDefOverride(LowPileupSelector, 0);

   protected:
    virtual void SetBranchesBacon() override;
    void LoadBranchesBacon(
        Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void SetComposite() {}
};

#endif
