#ifndef ZSelector_h
#define ZSelector_h

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
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include "Analysis/VVAnalysis/interface/BranchManager.h"

class ZSelector : public SelectorBase {
public :
    ScaleFactor* pileupSF_ = nullptr;
    ScaleFactor* muonSF_ = nullptr;
    ScaleFactor* eIdSF_  = nullptr;
    ScaleFactor* eGsfSF_ = nullptr;
    ScaleFactor* mIdSF_ = nullptr;
    ScaleFactor* mIsoSF_ = nullptr;
    TEfficiency* prefireEff_ = nullptr;

    // Derived values
    bool singleLepton_;
    bool passesLeptonVeto;
    bool passesTrigger;
    UInt_t nCBVIDTightElec;
    UInt_t nCBVIDVetoElec;
    UInt_t nTightIdMuon;
    UInt_t nMediumIdMuon;
    UInt_t nTightIsoMuon;
    UInt_t nLooseIsoMuon;
    Float_t weight;
    Float_t l1Eta;
    Float_t l2Eta;
    Float_t l3Eta;
    Float_t l1Pt;
    Float_t l2Pt;
    Float_t l1Phi;
    Float_t l2Phi;
    Float_t l1Mass;
    Float_t l2Mass;
    Float_t ZMass;
    Float_t Zy;
    Float_t ZEta;
    Float_t ZPt;
    Bool_t l1IsTight;
    Bool_t l2IsTight;

    // Values read from file
    Float_t genWeight;
    Float_t MET;
    Float_t type1_pfMETPhi;
    static const unsigned int N_KEEP_MU_E_ = 15;
    UInt_t nElectron;
    UInt_t nMuon;
    Float_t Electron_pt[N_KEEP_MU_E_];
    Float_t Electron_eta[N_KEEP_MU_E_];
    Float_t Electron_phi[N_KEEP_MU_E_];
    Float_t Electron_mass[N_KEEP_MU_E_];
    Float_t Muon_pt[N_KEEP_MU_E_];
    Float_t Muon_eta[N_KEEP_MU_E_];
    Float_t Muon_phi[N_KEEP_MU_E_];
    Float_t Muon_mass[N_KEEP_MU_E_];
    Int_t Electron_cutBased[N_KEEP_MU_E_];
    Int_t Electron_charge[N_KEEP_MU_E_];
    Int_t Muon_charge[N_KEEP_MU_E_];
    Bool_t Muon_tightId[N_KEEP_MU_E_];
    Bool_t Muon_mediumId[N_KEEP_MU_E_];
    Float_t Muon_pfRelIso04_all[N_KEEP_MU_E_];
    Bool_t Dimuon_Trigger;
    Bool_t Dielectron_Trigger;
    Bool_t SingleMuon_Trigger;
    Bool_t SingleElectron_Trigger;
    Int_t numPU;
    
    TBranch* b_genWeight;
    TBranch* b_MET;
    TBranch* b_type1_pfMETPhi;
    TBranch* b_nElectron;
    TBranch* b_nMuon;
    TBranch* b_Electron_pt;
    TBranch* b_Electron_eta;
    TBranch* b_Electron_phi;
    TBranch* b_Electron_mass;
    TBranch* b_Muon_pt;
    TBranch* b_Muon_eta;
    TBranch* b_Muon_phi;
    TBranch* b_Muon_mass;
    TBranch* b_Electron_cutBased;
    TBranch* b_Electron_charge;
    TBranch* b_Muon_charge;
    TBranch* b_Muon_tightId;
    TBranch* b_Muon_mediumId;
    TBranch* b_Muon_pfRelIso04_all;
    TBranch* b_Dimuon_Trigger;
    TBranch* b_Dielectron_Trigger;
    TBranch* b_SingleMuon_Trigger;
    TBranch* b_SingleElectron_Trigger;
    TBranch* b_numPU;

    // UWVV specific branches
    Float_t nTruePU;

    TBranch* b_nTruePU;
    TBranch* b_nCBVIDTightElec;
    TBranch* b_nCBVIDHLTSafeElec;
    TBranch* b_nWZTightMuon;
    TBranch* b_nWZMediumMuon;
    TBranch* b_l1IsTight;
    TBranch* b_l2IsTight;
    TBranch* b_l1Eta;
    TBranch* b_l2Eta;
    TBranch* b_l1Pt;
    TBranch* b_l2Pt;
    TBranch* b_ZMass;

  BranchManager b;
  
  // Readers to access the data (delete the ones you do not need).
  virtual void    SetScaleFactors() override;
  virtual void    Init(TTree *tree) override;
  ZSelector(TTree * /*tree*/ =0) { }
  ~ZSelector() { }
  virtual void    SetupNewDirectory() override;

    ClassDefOverride(ZSelector,0);

protected:
    bool zlep1IsTight();
    bool zlep2IsTight();
    bool tightZLeptons();
    virtual void    SetBranchesNanoAOD() override;
    virtual void    SetBranchesUWVV() override;
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void ApplyScaleFactors();
    void SetComposite();
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif


