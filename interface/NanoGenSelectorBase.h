#ifndef NanoGenSelectorBase_h
#define NanoGenSelectorBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <exception>
#include <iostream>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/SelectorTools/interface/ScaleFactor.h"
#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Analysis/SelectorTools/interface/helpers.h"
#include "Analysis/SelectorTools/interface/BranchManager.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"

class NanoGenSelectorBase : public SelectorBase {
public :
    bool isZ_ = false;
    bool altPdf_ = false;
    bool storeCenPdfs_ = false;
    std::vector<int> cenPdfWeightsToStore_ = {};
    bool isMinnlo_ = false;
    int pdfMaxStore_ = 150;
    int pdfCenWeight_ = 0;
    std::string pdfSet_ = "all";
    std::vector<ScaleFactor*> scetlibCorrs_;
    std::unordered_map<std::string, std::vector<int>> minnloPdfMap = {
        {"nnpdf31", {0}},//,1,2,3,4,5,6,7,8}},
        {"nnpdf31cmsw1", {9}},
        {"nnpdf31cmsw2", {10}},
        {"nnpdf31cmsw3", {11}},
        {"nnpdf31cmsw4", {12}},
        {"nnpdf30", {13,15,16,}}, // 14 and 17 are bigger alpha_s vars
        {"ct18", {18}},
        {"ct18z", {18}},
        {"mmht", {19, 20, 21}},
        {"hera", {23,24,25,26,28,29,30}},
    };
    int nScetlibWeights_ = 45;
    PDFWeightsHelper pdfweightshelper_;
    // Derived values
    reco::GenParticleCollection leptons;
    reco::GenParticleCollection bareLeptons;
    reco::GenParticleCollection dressedLeptons;
    reco::GenParticleCollection bornLeptons;
    reco::GenParticleCollection lheLeptons;
    reco::GenParticleCollection preFSRLeptons;
    reco::GenParticleCollection bornNeutrinos;
    reco::GenParticleCollection lheNeutrinos;
    reco::GenParticleCollection preFSRNeutrinos;
    reco::GenParticleCollection fsneutrinos;
    reco::GenParticleCollection neutrinos;
    reco::GenParticleCollection photons;
    std::vector<LorentzVector> jets;
    std::vector<LorentzVector> genjets;
    std::vector<LorentzVector> lhejets;
    LorentzVector genMet;

    int centralWeightIndex_ = 0;
    unsigned int nLeptons_ = 1;
    static const unsigned int MAX_PDF_SETS = 50;
    static const unsigned int N_LHESCALE_WEIGHTS_ = 1000;
    static const unsigned int N_LHEPDF_WEIGHTS_ = 2000;
    static const unsigned int N_LHEPDFAS_WEIGHTS_ = 102;
    static const unsigned int N_MC2HESSIAN_WEIGHTS_ = 60;
    float weight;
    float mVlhe;
    float MV_GEN_;
    float GAMMAV_GEN_;
    bool n3llcorr_ = false;
    bool nnlops_ = false;
    int weightSuppress_ = 0;
    int thweightSuppress_ = 0;
    bool weightSignOnly_ = false;
    bool doTheoryVars_ = false;
    bool doMC2H_ = false;
    bool doPhotons_ = true;
    bool nNeutrinos_ = 0;
    bool doFiducial_ = false;
    bool doBorn_ = false;
    bool doLHE_ = false;
    bool doPreFSR_ = false;
    bool doBareLeptons_ = false;
    float ratio_mass;
    float refWeight = 1;
    float rescaleWeight_ = 1.; //If central weight is modified

    TH1D* mcWeights_;
    TH1D* mcPdfWeights_;
    TH1D* hesPdfWeights_;
    TH1D* scaleWeightsHist_;
    
    double LHEHessianPdfWeight[N_MC2HESSIAN_WEIGHTS_];
    // Values read from file
    TTreeReader     fReader;
    TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
    
    UInt_t nLHEScaleWeight = 1;
    Float_t LHEScaleWeight[N_LHESCALE_WEIGHTS_];
    UInt_t nLHEPdfWeight = 0;
    Float_t LHEPdfWeight[N_LHEPDF_WEIGHTS_];
    Float_t LHEPdfWeights[MAX_PDF_SETS][N_LHEPDF_WEIGHTS_];
    std::array<UInt_t, MAX_PDF_SETS> nLHEPdfWeights = {{0}};
    UInt_t nLHEScaleWeightAltSet1 = 0;
    Float_t LHEScaleWeightAltSet1[N_LHESCALE_WEIGHTS_];
    UInt_t nLHEUnknownWeight = 0;
    Float_t LHEUnknownWeight[100];
    UInt_t nMEParamWeight = 0;
    Float_t MEParamWeight[100];
    UInt_t nLHEUnknownWeightAltSet1 = 0;
    Float_t LHEUnknownWeightAltSet1[100];

    std::array<TBranch*, MAX_PDF_SETS> b_nLHEPdfWeights;
    std::array<TBranch*, MAX_PDF_SETS> b_LHEPdfWeights;
    TBranch* b_nLHEPdfWeight;
    TBranch* b_LHEPdfWeight;
    TBranch* b_nLHEScaleWeightAltSet1;
    TBranch* b_LHEScaleWeightAltSet1;
    TBranch* b_nLHEUnknownWeight;
    TBranch* b_LHEUnknownWeight;
    TBranch* b_nLHEUnknownWeightAltSet1;
    TBranch* b_LHEUnknownWeightAltSet1;

    bool scaleWeights_ = false;
    bool altScaleWeights_ = false;
    std::array<bool, MAX_PDF_SETS> pdfWeights_ = {{false}};
    bool unknownWeights_ = false;
    bool paramWeights_ = false;
    bool paramWeightsUpd_ = false;
    bool unknownWeightsAlt_ = false;

    TTreeReaderValue<UInt_t> nGenDressedLepton = {fReader, "nGenDressedLepton"};
    TTreeReaderArray<Bool_t> GenDressedLepton_hasTauAnc = {fReader, "GenDressedLepton_hasTauAnc"};
    TTreeReaderArray<Float_t> GenDressedLepton_pt = {fReader, "GenDressedLepton_pt"};
    TTreeReaderArray<Float_t> GenDressedLepton_eta = {fReader, "GenDressedLepton_eta"};
    TTreeReaderArray<Float_t> GenDressedLepton_phi = {fReader, "GenDressedLepton_phi"};
    TTreeReaderArray<Float_t> GenDressedLepton_mass = {fReader, "GenDressedLepton_mass"};
    TTreeReaderArray<Int_t> GenDressedLepton_pdgId = {fReader, "GenDressedLepton_pdgId"};
    TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
    TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
    TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
    TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
    TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
    TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
    TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
    TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
    TTreeReaderValue<UInt_t> nLHEPart = {fReader, "nLHEPart"};
    TTreeReaderArray<Float_t> LHEPart_pt = {fReader, "LHEPart_pt"};
    TTreeReaderArray<Float_t> LHEPart_eta = {fReader, "LHEPart_eta"};
    TTreeReaderArray<Float_t> LHEPart_phi = {fReader, "LHEPart_phi"};
    TTreeReaderArray<Float_t> LHEPart_mass = {fReader, "LHEPart_mass"};
    TTreeReaderArray<Int_t> LHEPart_pdgId = {fReader, "LHEPart_pdgId"};
    TTreeReaderValue<UInt_t> nGenJet = {fReader, "nGenJet"};
    TTreeReaderArray<Float_t> GenJet_pt = {fReader, "GenJet_pt"};
    TTreeReaderArray<Float_t> GenJet_eta = {fReader, "GenJet_eta"};
    TTreeReaderArray<Float_t> GenJet_phi = {fReader, "GenJet_phi"};
    TTreeReaderArray<Float_t> GenJet_mass = {fReader, "GenJet_mass"};
    //TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
    //TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
    TTreeReaderValue<Float_t> MET_fiducialGenPt = {fReader, "MET_fiducialGenPt"};
    TTreeReaderValue<Float_t> MET_fiducialGenPhi = {fReader, "MET_fiducialGenPhi"};
    float ht;
    float ptVcorr;
    float yVcorr;
    float mVcorr;
    
    BranchManager b;
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    virtual void    SetScaleFactors() override;
    NanoGenSelectorBase(TTree * /*tree*/ =0) { }
    ~NanoGenSelectorBase() { }
    virtual void    SetupNewDirectory() override;

    ClassDefOverride(NanoGenSelectorBase,0);

protected:
    virtual void    SetBranchesNanoAOD() override;
    virtual void    FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override {}
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void SetComposite() {}
    bool overlapsCollection(const LorentzVector& cand, reco::GenParticleCollection& collection, const float deltaRCut, size_t maxCompare);
    void buildHessian2MCSet();
    reco::GenParticle makeGenParticle(int pdgid, int status, float pt, float eta, float phi, float m);
    double breitWignerWeight(double offset);
};

#endif



