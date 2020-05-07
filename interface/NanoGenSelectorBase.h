#ifndef NanoGenSelectorBase_h
#define NanoGenSelectorBase_h

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/VVAnalysis/interface/BranchManager.h"
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include "Analysis/VVAnalysis/interface/helpers.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"

class NanoGenSelectorBase : public SelectorBase {
   public:
    PDFWeightsHelper pdfweightshelper_;
    // Derived values
    reco::GenParticleCollection leptons;
    reco::GenParticleCollection bareLeptons;
    reco::GenParticleCollection dressedLeptons;
    reco::GenParticleCollection bornLeptons;
    reco::GenParticleCollection lheLeptons;
    reco::GenParticleCollection bornNeutrinos;
    reco::GenParticleCollection lheNeutrinos;
    reco::GenParticleCollection fsneutrinos;
    reco::GenParticleCollection neutrinos;
    reco::GenParticleCollection photons;
    std::vector<LorentzVector> jets;
    LorentzVector genMet;

    int centralWeightIndex_ = 0;
    unsigned int nLeptons_ = 1;
    static const unsigned int N_LHESCALE_WEIGHTS_ = 1000;
    static const unsigned int N_LHEPDF_WEIGHTS_ = 2000;
    static const unsigned int N_LHEPDFAS_WEIGHTS_ = 102;
    static const unsigned int N_MC2HESSIAN_WEIGHTS_ = 60;
    float weight;
    float mVlhe;
    float MV_GEN_;
    float GAMMAV_GEN_;
    bool nnlops_ = false;
    int weightSuppress_ = 0;
    bool weightSignOnly_ = false;
    bool doTheoryVars_ = false;
    bool doMC2H_ = false;
    bool doPhotons_ = true;
    bool doNeutrinos_ = true;
    bool doFiducial_ = false;
    bool doBorn_ = true;
    bool doBareLeptons_ = true;

    float refWeight = 1;

    TH1D* mcWeights_;
    TH1D* mcPdfWeights_;
    TH1D* hesPdfWeights_;
    TH1D* scaleWeights_;

    double LHEHessianPdfWeight[N_MC2HESSIAN_WEIGHTS_];
    // Values read from file
    TTreeReader fReader;
    TTreeReaderValue<Float_t> genWeight = {fReader, "genWeight"};
    TTreeReaderValue<UInt_t> nLHEScaleWeight = {fReader, "nLHEScaleWeight"};
    TTreeReaderArray<Float_t> LHEScaleWeight = {fReader, "LHEScaleWeight"};

    UInt_t nLHEPdfWeight = 0;
    Float_t LHEPdfWeight[N_LHEPDF_WEIGHTS_];
    UInt_t nLHEScaleWeightAltSet1 = 0;
    Float_t LHEScaleWeightAltSet1[N_LHESCALE_WEIGHTS_];
    UInt_t nLHEUnknownWeight = 0;
    Float_t LHEUnknownWeight[100];
    UInt_t nLHEUnknownWeightAltSet1 = 0;
    Float_t LHEUnknownWeightAltSet1[100];

    TBranch* b_nLHEPdfWeight;
    TBranch* b_LHEPdfWeight;
    TBranch* b_nLHEScaleWeightAltSet1;
    TBranch* b_LHEScaleWeightAltSet1;
    TBranch* b_nLHEUnknownWeight;
    TBranch* b_LHEUnknownWeight;
    TBranch* b_nLHEUnknownWeightAltSet1;
    TBranch* b_LHEUnknownWeightAltSet1;

    bool altScaleWeights_ = false;
    bool pdfWeights_ = false;
    bool unknownWeights_ = false;
    bool unknownWeightsAlt_ = false;

    TTreeReaderValue<UInt_t> nGenDressedLepton = {fReader, "nGenDressedLepton"};
    TTreeReaderArray<Bool_t> GenDressedLepton_hasTauAnc = {
        fReader, "GenDressedLepton_hasTauAnc"};
    TTreeReaderArray<Float_t> GenDressedLepton_pt = {fReader,
                                                     "GenDressedLepton_pt"};
    TTreeReaderArray<Float_t> GenDressedLepton_eta = {fReader,
                                                      "GenDressedLepton_eta"};
    TTreeReaderArray<Float_t> GenDressedLepton_phi = {fReader,
                                                      "GenDressedLepton_phi"};
    TTreeReaderArray<Float_t> GenDressedLepton_mass = {fReader,
                                                       "GenDressedLepton_mass"};
    TTreeReaderArray<Int_t> GenDressedLepton_pdgId = {fReader,
                                                      "GenDressedLepton_pdgId"};
    TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
    TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
    TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
    TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
    TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
    TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
    TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
    TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader,
                                                   "GenPart_statusFlags"};
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
    // TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
    // TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
    TTreeReaderValue<Float_t> MET_fiducialGenPt = {fReader,
                                                   "MET_fiducialGenPt"};
    TTreeReaderValue<Float_t> MET_fiducialGenPhi = {fReader,
                                                    "MET_fiducialGenPhi"};
    float ht;
    float ptVlhe;

    BranchManager b;

    // Readers to access the data (delete the ones you do not need).
    virtual void Init(TTree* tree) override;
    NanoGenSelectorBase(TTree* /*tree*/ = 0) {}
    ~NanoGenSelectorBase() {}
    virtual void SetupNewDirectory() override;

    ClassDefOverride(NanoGenSelectorBase, 0);

   protected:
    virtual void SetBranchesNanoAOD() override;
    virtual void FillHistograms(
        Long64_t entry, std::pair<Systematic, std::string> variation) override {
    }
    void LoadBranchesNanoAOD(
        Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void SetComposite() {}
    bool overlapsCollection(const LorentzVector& cand,
                            reco::GenParticleCollection& collection,
                            const float deltaRCut, size_t maxCompare);
    void buildHessian2MCSet();
    reco::GenParticle makeGenParticle(int pdgid, int status, float pt,
                                      float eta, float phi, float m);
    double breitWignerWeight(double offset);
};

#endif
