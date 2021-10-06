#ifndef WGenSelector_h
#define WGenSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/SelectorTools/interface/ScaleFactor.h"
#include "Analysis/SelectorTools/interface/NanoGenSelectorBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class WGenSelector : public NanoGenSelectorBase {
public :
    // Derived values
    LorentzVector nu;
    LorentzVector wCand;
    LorentzVector wCandMet;
    float mTtrue;
    float mTmet;
    float cenWeight;
    bool doMassVar_ = false;
    bool doMuonVar_ = false;
    
    float ptl_smear;
    float ptl_smear_fill;
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    WGenSelector(TTree * /*tree*/ =0) { }
    ~WGenSelector() { }

    ClassDefOverride(WGenSelector,0);

protected:
    virtual void SetComposite() override;
    virtual void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif



