#ifndef FakeRateSelector_h
#define FakeRateSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "Math/GenVector/VectorUtil.h"
#include <exception>
#include <iostream>
#include <vector>
#include <unordered_map>

// Headers needed by this particular selector
#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include "Analysis/VVAnalysis/interface/ThreeLepSelector.h"
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/BranchManager.h"
#include "Analysis/VVAnalysis/interface/GoodParticle.h"


class FakeRateSelector : public SelectorBase { 
public :
    int sel_MVAStudy, sel_FakeRate;
    


    
    std::unordered_map<int, TH2D*> passingTight2D_map;
    std::unordered_map<int, TH1D*> passingTight1DPt_map;
    std::unordered_map<int, TH1D*> passingTight1DEta_map;
    std::unordered_map<int, TH2D*> passingLoose2D_map;
    std::unordered_map<int, TH1D*> passingLoose1DPt_map;
    std::unordered_map<int, TH1D*> passingLoose1DEta_map;        
                                            

    TTreeReader fReader;
    TTreeReaderValue<Bool_t> HLT_Mu8 = {fReader, "HLT_Mu8"};
    TTreeReaderValue<Bool_t> HLT_Ele8_CaloIdM_TrackIdM_PFJet30 = {fReader, "HLT_Ele8_CaloIdM_TrackIdM_PFJet30"};
    
    double weight = 1.0;

    const float FR_MAX_PT_ = 50;
    const float FR_MAX_ETA_ = 2.5;
                        
    // Readers to access the data (delete the ones you do not need).
    virtual void    SetupNewDirectory() override;
    virtual void    FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void clearValues();
    void fillReco(std::vector<GenPart>& genList, const std::vector<GoodPart>& recoList);
    
    // overloaded or necessary functions
    virtual void    SetBranchesNanoAOD() override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    // Readers to access the data (delete the ones you do not need).
    virtual void    SlaveBegin(TTree *tree) override {return;}
    virtual void    Init(TTree *tree) override; //{return;}
    ///ignore
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override {return;}
    virtual void    SetBranchesUWVV() override {return;}

    
    ThreeLepSelector Analyzer;
    ClassDefOverride(FakeRateSelector,0);
};

#endif
