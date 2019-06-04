
#ifndef SameSignLept_h
#define SameSignLept_h

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

class SameSignLept : public SelectorBase {
public :
    ScaleFactor* pileupSF_;
    ScaleFactor* muonSF_;
    ScaleFactor* eIdSF_ ;
    ScaleFactor* eGsfSF_;
    ScaleFactor* mIdSF_;
    ScaleFactor* mIsoSF_;

    bool  returning_value;

    // Derived values
    bool passesLeptonVeto;
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
    Bool_t l1IsTight;
    Bool_t l2IsTight;
    Bool_t l1IsMed;
    Bool_t l2IsMed;

    // Values read from file
    Float_t genWeight;
    Float_t MET;
    Float_t type1_pfMETPhi;
    static const unsigned int N_KEEP_MU_E_ = 10;
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
    UChar_t Muon_pfIsoId[N_KEEP_MU_E_];
    Float_t Muon_pfRelIso04_all[N_KEEP_MU_E_];
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
    TBranch* b_Muon_pfIsoId;
    TBranch* b_Muon_pfRelIso04_all;
    TBranch* b_numPU;

    // Histograms
    TH1D* cutflow_ee_;
    TH1D* cutflow_mm_;
    TH1D* ZMass_ee_;
    TH1D* ZMass_mm_;
    TH1D* ptl1_ee_; 
    TH1D* ptl1_mm_; 
    TH1D* ptl2_ee_; 
    TH1D* ptl2_mm_; 
    TH1D* l1eta_ee_;
    TH1D* l2eta_ee_;
    TH1D* l1eta_mm_;
    TH1D* l2eta_mm_;
    TH1D* l1phi_mm_;
    TH1D* l2phi_mm_;
    TH1D* l1phi_ee_;
    TH1D* l2phi_ee_;
 
    // Readers to access the data (delete the ones you do not need).
    //virtual void    SetScaleFactors() override;
    virtual void    Init(TTree *tree) override;
    SameSignLept(TTree * /*tree*/ =0) { }
    ~SameSignLept() { }
    virtual void    SetupNewDirectory() override;

    ClassDefOverride(SameSignLept,0);

protected:
    bool zlep1IsTight();
    bool zlep2IsTight();
    bool zlep1IsMed();
    bool zlep2IsMed();
    bool tightZLeptons();
    bool medZLeptons();
    bool tight_med_ZLeptons();

    virtual void    SetBranchesNanoAOD() override;
    virtual void    SetBranchesUWVV() override;
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void ApplyScaleFactors();
    void SetMass();
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif


