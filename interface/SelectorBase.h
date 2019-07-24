#ifndef SelectorBase_h
#define SelectorBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TEfficiency.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"

class SelectorBase : public TSelector {
public :
    std::map<std::string, ScaleFactor*> scaleFactors;
    TEfficiency* prefireEff_;
    
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
    enum NtupleType {
        UWVV,
        NanoAOD,
    };

    std::vector<std::string> allChannels_ = {};
    std::vector<std::string> hists1D_ = {};
    std::vector<std::string> hists2D_ = {};

    enum Channel {
        e,
        m,
        ee,
        mm,
        eee,
        eem,
        emm,
        mmm,
        eeee,
        eemm,
        mmee,
        mmmm,
        Inclusive,
        Unknown,
    }; 
    enum Selection {
        tightleptons,
        ZZGenFiducial,
        Wselection,
        Zselection,
        Wselection_Full,
        FakeRateSelectionLoose,
        FakeRateSelectionTight,
        VBSselection_Loose,
        VBSselection_NoZeppenfeld,
        VBSselection_Tight,
        VBSselection_Loose_Full,
        VBSselection_NoZeppenfeld_Full,
        VBSselection_Tight_Full,
        VBSBackgroundControl,
        VBSBackgroundControlATLAS,
        VBSBackgroundControl_Full,
        VBSBackgroundControlLoose,
        VBSBackgroundControlLoose_Full,
        Inclusive2Jet,
        Inclusive2Jet_Full,
        TightWithLooseVeto,
    };

    std::map<std::string, Selection> selectionMap_ = {
      {"tightleptons", tightleptons},
        {"ZZGenFiducial", ZZGenFiducial},
        {"Wselection", Wselection},
        {"Zselection", Zselection},
        {"Wselection_Full", Wselection_Full},
        {"FakeRateSelectionLoose", FakeRateSelectionLoose},
        {"FakeRateSelectionTight", FakeRateSelectionTight},
        {"VBSselection_Loose", VBSselection_Loose},
        {"VBSselection_NoZeppenfeld", VBSselection_NoZeppenfeld},
        {"VBSselection_Tight", VBSselection_Tight},
        {"VBSselection_Loose_Full", VBSselection_Loose_Full},
        {"VBSselection_NoZeppenfeld_Full", VBSselection_NoZeppenfeld_Full},
        {"VBSselection_Tight_Full", VBSselection_Tight_Full},
        {"VBSBackgroundControl", VBSBackgroundControl},
        {"VBSBackgroundControlATLAS", VBSBackgroundControlATLAS},
        {"VBSBackgroundControl_Full", VBSBackgroundControl_Full},
        {"VBSBackgroundControlLoose", VBSBackgroundControlLoose},
        {"VBSBackgroundControlLoose_Full", VBSBackgroundControlLoose_Full},
        {"Inclusive2Jet", Inclusive2Jet},
        {"Inclusive2Jet_Full", Inclusive2Jet_Full},
        {"TightWithLooseVeto", TightWithLooseVeto},
    };

    enum Year {
        yrdefault,
        yr2016,
        yr2017,
        yr2018
    };

    std::map<std::string, Year> yearMap_ = {
        {"default", yrdefault},
        {"2016", yr2016},
        {"2017", yr2017},
        {"2018", yr2018},
    };
    
    std::map<std::string, Channel> channelMap_ = {
        {"e", e},
        {"m", m},
        {"ee", ee},
        {"mm", mm},
        {"eee", eee},
        {"eem", eem},
        {"emm", emm},
        {"mmm", mmm},
        {"eeee", eeee},
        {"eemm", eemm},
        {"mmee", mmee},
        {"mmmm", mmmm},
        {"eeeeGen", eeee},
        {"eemmGen", eemm},
        {"mmeeGen", mmee},
        {"mmmmGen", mmmm},
        {"Inclusive", Inclusive},
    };

    enum Systematic {
        Central,
        jetEnergyScaleUp,
        jetEnergyScaleDown,
        jetEnergyResolutionUp,
        jetEnergyResolutionDown,
        metUnclusteredEnergyUp,
        metUnclusteredEnergyDown,
        muonEfficiencyUp,
        muonEfficiencyDown,
        muonScaleUp,
        muonScaleDown,
        electronEfficiencyUp,
        electronEfficiencyDown,
        electronScaleUp,
        electronScaleDown,
        pileupUp,
        pileupDown,
    }; 


    std::map<Systematic, std::string> variations_ = {{Central, ""}};
    std::map<Systematic, std::string> systematics_ = {};

    TList *currentHistDir_{nullptr};
    TH1D* sumWeightsHist_;

    bool doSystematics_;
    bool addSumweights_;
    bool applyScaleFactors_;
    bool applyPrefiringCorr_;
    
    // Readers to access the data (delete the ones you do not need).
    SelectorBase(TTree * /*tree*/ =0) { }
    virtual ~SelectorBase() { }
    virtual void    SetScaleFactors();
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    // We'll collect pointers to objects from derived classes
    // as they are registered with AddObject, and update them to
    // the new object when a dataset is switched
    std::set<TNamed**> allObjects_;
    // Derived classes override (and call) this to register new objects
    // With AddObject<Type>(localPtr, ...);
    virtual void SetupNewDirectory();

    template<typename T, typename... Args>
    void AddObject(T* &ptr, Args... args) {
        static_assert(std::is_base_of<TNamed, T>::value, "Objects must inheirit from ROOT TNamed to be streamable from PROOF sessions");
        ptr = new T(args...);
        ptr->SetDirectory(gROOT);
        currentHistDir_->Add(ptr);
        allObjects_.insert((TNamed**) &ptr);
    };
    void UpdateDirectory();
    ClassDef(SelectorBase,0);

protected:
    // Maps to the histogram pointers themselves
    std::map<std::string, TH1D*> histMap1D_ = {};
    std::map<std::string, TH2D*> histMap2D_ = {};
    std::map<std::string, TH2D*> weighthists_ = {};
    std::map<std::string, TH3D*> weighthists2D_ {};

    // The histograms for which you also want systematic variations
    std::vector<std::string> systHists_ = {};
    std::vector<std::string> systHists2D_ = {};

    void    SetBranches();
    virtual void    SetBranchesUWVV() { }
    virtual void    SetBranchesNanoAOD() { }
    void    LoadBranches(Long64_t entry, std::pair<Systematic, std::string> variation);
    virtual void    LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) { }
    virtual void    LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { }
    virtual void    FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { }
    std::string name_ = "Unnamed";
    std::string channelName_ = "Unnamed";
    Channel channel_ = Unknown;
    NtupleType ntupleType_ = NanoAOD;
    std::string selectionName_ = "tightleptons";
    Selection selection_ = tightleptons;
    Year year_ = yrdefault;
    bool isMC_;
    float GetPrefiringEfficiencyWeight(std::vector<float>* jetPt, std::vector<float>* jetEta);
    virtual std::string GetNameFromFile() { return ""; }
    void InitializeHistogramFromConfig(std::string name, std::string channel, std::vector<std::string> histData);
    void InitializeHistogramsFromConfig();
    std::vector<std::string> ReadHistDataFromConfig(std::string histDataString);
    std::string getHistName(std::string histName, std::string variationName, std::string channel);
    std::string getHistName(std::string histName, std::string variationName);
    template<typename T, typename... Args>
    void SafeHistFill(std::map<std::string, T*> container, 
            std::string histname, Args... args) {
        if (container[histname] != nullptr)
            container[histname]->Fill(args...);
    };
};

#endif

