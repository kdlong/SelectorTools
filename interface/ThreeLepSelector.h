#ifndef ThreeLepSelector_h
#define ThreeLepSelector_h

// #include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include "Analysis/VVAnalysis/interface/BranchManager.h"
#include "Analysis/VVAnalysis/interface/GoodParticle.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVector;

enum PID {PID_MUON = 13, PID_ELECTRON = 11, PID_BJET = 5, PID_CJET = 4, PID_JET};

class ThreeLepSelector : public SelectorBase {
public :
#include "Analysis/VVAnalysis/interface/FourTopScales.h"
    
    /*****************************************/
    /* ____  ____   ___  __  __   ___ __  __ */
    /* || )) || \\ // \\ ||\ ||  //   ||  || */
    /* ||=)  ||_// ||=|| ||\\|| ((    ||==|| */
    /* ||_)) || \\ || || || \||  \\__ ||  || */
    /*****************************************/

    ScaleFactor* pileupSF_;
    ScaleFactor* muonSF_;
    ScaleFactor* eIdSF_ ;
    ScaleFactor* eGsfSF_;
    ScaleFactor* mIdSF_;
    ScaleFactor* mIsoSF_;

    // Common variables
    Float_t genWeight;
    Float_t MET;
    Float_t type1_pfMETPhi;

    int MVAStudy, FourTopMVAEl, FourTopCutBasedEl, FakeRate;
    int yr2016, yr2017, yr2018;

    TTreeReader     fReader, tmpReader;
    
    
    TTreeReaderArray<Float_t>   Electron_pt = {fReader, "Electron_pt"};
    TTreeReaderArray<Float_t>   Electron_eta = {fReader, "Electron_eta"};
    TTreeReaderArray<Float_t>   Electron_phi = {fReader, "Electron_phi"};
    TTreeReaderArray<Float_t>   Electron_mass = {fReader, "Electron_mass"};
    TTreeReaderArray<Int_t>     Electron_cutBased = {fReader, "Electron_cutBased"};
    TTreeReaderArray<Int_t>     Electron_charge = {fReader, "Electron_charge"};
    TTreeReaderArray<Float_t>   Electron_MVA = {fReader, "Electron_MVA"};
    TTreeReaderArray<Float_t>   Electron_miniPFRelIso_all = {fReader, "Electron_miniPFRelIso_all"};
    TTreeReaderArray<Float_t>   Electron_dxy = {fReader, "Electron_dxy"};
    TTreeReaderArray<Float_t>   Electron_dz = {fReader, "Electron_dz"};
    TTreeReaderArray<Float_t>   Electron_sip3d = {fReader, "Electron_sip3d"};
    TTreeReaderArray<Bool_t>    Electron_convVeto = {fReader, "Electron_convVeto"};
    TTreeReaderArray<UChar_t>   Electron_lostHits = {fReader, "Electron_lostHits"};
    TTreeReaderArray<Int_t>     Electron_tightCharge = {fReader, "Electron_tightCharge"};
    TTreeReaderArray<Float_t>   Electron_sieie = {fReader, "Electron_sieie"};
    TTreeReaderArray<Float_t>   Electron_hoe = {fReader, "Electron_hoe"};
    TTreeReaderArray<Float_t>   Electron_deltaEtaSC = {fReader, "Electron_deltaEtaSC"};
    TTreeReaderArray<Float_t>   Electron_eInvMinusPInv = {fReader, "Electron_eInvMinusPInv"};
    TTreeReaderArray<Float_t>   Electron_dr03EcalRecHitSumEt = {fReader, "Electron_dr03EcalRecHitSumEt"};
    TTreeReaderArray<Float_t>   Electron_dr03HcalDepth1TowerSumEt = {fReader, "Electron_dr03HcalDepth1TowerSumEt"};
    TTreeReaderArray<Float_t>   Electron_dr03TkSumPt = {fReader, "Electron_dr03TkSumPt"};
    //TTreeReaderArray<Int_t>     Electron_vidBitmap = {fReader, "Electron_vidBitmap"};
    TTreeReaderArray<Float_t>   Electron_jetRelIso = {fReader, "Electron_jetRelIso"};
    TTreeReaderArray<Float_t>   Electron_eCorr = {fReader, "Electron_eCorr"};

    TTreeReaderArray<Float_t>   Muon_pt  = {fReader, "Muon_pt"};
    TTreeReaderArray<Float_t>   Muon_eta = {fReader, "Muon_eta"};
    TTreeReaderArray<Float_t>   Muon_phi = {fReader, "Muon_phi"};
    TTreeReaderArray<Float_t>   Muon_mass = {fReader, "Muon_mass"};
    TTreeReaderArray<Int_t>     Muon_charge = {fReader, "Muon_charge"};
    TTreeReaderArray<Bool_t>    Muon_tightId = {fReader, "Muon_tightId"};
    TTreeReaderArray<Bool_t>    Muon_mediumId = {fReader, "Muon_mediumId"};
    TTreeReaderArray<UChar_t>   Muon_tkIsoId = {fReader, "Muon_tkIsoId"};
    TTreeReaderArray<Float_t>   Muon_pfRelIso04_all = {fReader, "Muon_pfRelIso04_all"};
    TTreeReaderArray<Float_t>   Muon_miniPFRelIso_all = {fReader, "Muon_miniPFRelIso_all"};
    TTreeReaderArray<Float_t>   Muon_dxy = {fReader, "Muon_dxy"};
    TTreeReaderArray<Float_t>   Muon_dz = {fReader, "Muon_dz"};
    TTreeReaderArray<Float_t>   Muon_sip3d = {fReader, "Muon_sip3d"};
    TTreeReaderArray<Bool_t>    Muon_isGlobal = {fReader, "Muon_isGlobal"};
    TTreeReaderArray<Bool_t>    Muon_isTracker = {fReader, "Muon_isTracker"};
    TTreeReaderArray<Bool_t>    Muon_isPFcand = {fReader, "Muon_isPFcand"};
    TTreeReaderArray<Int_t>     Muon_tightCharge = {fReader, "Muon_tightCharge"};

    TTreeReaderArray<Float_t>   Jet_btagCSVV2 = {fReader, "Jet_btagCSVV2"};
    TTreeReaderArray<Float_t>   Jet_btagDeepB = {fReader, "Jet_btagDeepB"};
    TTreeReaderArray<Float_t>   Jet_eta = {fReader, "Jet_eta"};
    TTreeReaderArray<Float_t>   Jet_phi = {fReader, "Jet_phi"};
    TTreeReaderArray<Float_t>   Jet_pt    = {fReader, "Jet_pt"};
    TTreeReaderArray<Float_t>   Jet_mass = {fReader, "Jet_mass"};
    TTreeReaderArray<Int_t>     Jet_nConstituents = {fReader, "Jet_nConstituents"};
    TTreeReaderArray<Int_t>     Jet_jetId = {fReader, "Jet_jetId"};
    TTreeReaderArray<Int_t>     Jet_hadronFlavour = {fReader, "Jet_hadronFlavour"};
    TTreeReaderArray<Float_t>   Jet_rawFactor = {fReader, "Jet_rawFactor"};

    Int_t     numPU;
    Float_t   Pileup_nTrueInt;
    UInt_t nElectron;
    UInt_t    nMuon;
    UInt_t    nJet;
    UInt_t nGenPart;
    
    Bool_t Flag_goodVertices;
    Bool_t Flag_globalSuperTightHalo2016Filter;
    Bool_t Flag_HBHENoiseFilter;
    Bool_t Flag_HBHENoiseIsoFilter;
    Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter;
    Bool_t Flag_BadPFMuonFilter;
    Bool_t Flag_ecalBadCalibFilter;
    
    Bool_t HLT_DoubleMu8_Mass8_PFHT300;
    Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;
    Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;
    Bool_t HLT_AK8PFJet450;
    Bool_t HLT_PFJet450;

    /// Can be excluded
    TTreeReaderArray<Float_t>   Jet_L1 = {tmpReader, "Jet_L1"};
    TTreeReaderArray<Float_t>   Jet_L2L3 = {tmpReader, "Jet_L2L3"};
    TTreeReaderArray<Int_t> GenPart_pdgId = {tmpReader, "GenPart_pdgId"};
    TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {tmpReader, "GenPart_genPartIdxMother"};
    TTreeReaderArray<Int_t> GenPart_status = {tmpReader, "GenPart_status"};
    
    
    ClassDefOverride(ThreeLepSelector,0);

    /*******************************************************/
    /* __ __  ___  ____  __  ___  ____  __     ____  __    */
    /* || || // \\ || \\ || // \\ || )) ||    ||    (( \   */
    /* \\ // ||=|| ||_// || ||=|| ||=)  ||    ||==   \\    */
    /*  \V/  || || || \\ || || || ||_)) ||__| ||___ \_))   */
    /*******************************************************/
  
    Float_t weight;
    BranchManager b;
    std::vector<GoodPart> goodLeptons;
    std::vector<GoodPart> looseLeptons;
    std::vector<GoodPart> goodJets;
    std::vector<int> jetList;
    std::vector<int> bjetList;
    
    double HT;
    int nJets, nBJets;
    bool passZVeto;
    BTagCalibration calib;
    BTagCalibrationReader btag_reader; // central sys type
    TH2D* Beff_b;
    TH2D* Beff_j;
    ULong64_t event;
    UInt_t lumi;

    TH2D* h_btag_eff_b;
    TH2D* h_btag_eff_c;
    TH2D* h_btag_eff_udsg;

    /************************************************************/
    /* _____ __ __ __  __   ___ ______ __   ___   __  __  __    */
    /* ||    || || ||\ ||  //   | || | ||  // \\  ||\ || (( \   */
    /* ||==  || || ||\\|| ((      ||   || ((   )) ||\\||  \\    */
    /* ||    \\_// || \||  \\__   ||   ||  \\_//  || \|| \_))   */
    /************************************************************/

    void setupMuons();
    void setupElectrons();
    void setupJets();
    void setupChannel();
    void printInfo();
  
    bool isGoodMuon(size_t);
    bool isLooseMuon(size_t);
    bool isGoodJet(size_t);
    bool isGoodBJet(size_t);
    bool isGoodElectron(size_t);
    bool isLooseElectron(size_t);

    bool passMVACut(std::vector<std::vector<double> >, int);
    double mvaInterpolate(double pt, std::vector<double> );
    std::map<int, std::vector<std::vector<double> > >  mvaValues;

    
    bool isLooseMVAElectron(size_t);

    size_t getCloseJetIndex(LorentzVector&, double minDR=10);
    bool doesNotOverlap(size_t);
    bool passFullIso(LorentzVector&, double, double);
    bool doesPassZVeto(GoodPart&, std::vector<GoodPart>&);
    bool passTriggerEmu(size_t);
    double LepRelPt(LorentzVector&, LorentzVector&);
    LorentzVector get4Vector(PID, int);
    bool passFakeableCuts(GoodPart&);
    bool MetFilter();
    float getBtagEffFromFile(double, double, int);
    double getWDecayScaleFactor();
    std::vector<GoodPart>::iterator findJet(std::vector<GoodPart>::iterator&, int);
    std::map<std::string, TTree*> treeMap = {{"tree", nullptr}};
    float bNJets, bnBJets, bHT, bMET, bl1Pt, bl2Pt, blMass, bsphere, bCentral, bShape1, bShape2, bnLeps, bDilepCharge, bnlBJets, bntBJets, bnlLeps;
    float bjMass, bjdr, bj1Pt, bj2Pt, bj3Pt, bj4Pt, bj5Pt, bj6Pt, bj7Pt, bj8Pt, bb1Pt, bb2Pt, bb3Pt, bb4Pt;
    
    //// General Functions
    int getSRBin() const;
    void clearValues();
    void ApplyScaleFactors();

    // Overloaded or necesary functions
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void    SetBranchesNanoAOD() override;
    virtual void    SetupNewDirectory() override;
    // Readers to access the data (delete the ones you do not need).
    virtual void    SetScaleFactors() override;
    virtual void    Init(TTree *tree) override;

    ///ignore
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override {return;}
    virtual void    SetBranchesUWVV() override {return;}
};

#endif
