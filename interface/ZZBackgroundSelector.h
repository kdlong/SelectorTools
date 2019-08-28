#ifndef ZZBackgroundSelector_h
#define ZZBackgroundSelector_h

#include <vector>
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/ZZSelector.h"

class ZZBackgroundSelector : public ZZSelector {
public :
    ScaleFactor * fakeRate_allE_;
    ScaleFactor * fakeRate_allMu_;
  
    //PPPF represents 3P1F CR and PPFF represents 2P1F CR
    TH1D* Z1MassHistPPPF_;
    TH1D* Z1MassHistPPFF_;
    TH1D* Z2MassHistPPPF_;
    TH1D* Z2MassHistPPFF_;
    //4lMass Histogram
    TH1D* MassHistPPPF_;
    TH1D* MassHistPPFF_;

    TH1D* WeightsHistmmee_;
    TH1D* WeightsHisteemm_;
    
    virtual void    SetupNewDirectory() override;
    virtual void    SlaveBegin(TTree *tree) override;

    ClassDefOverride(ZZBackgroundSelector,0);
private:
    float getl3FakeRate(Long64_t entry);
    float getl4FakeRate(Long64_t entry);
    //In ZZ we only care about 3P1F and 2P2F, contributions are extremely low from other regions like FFFF,PFFF etc
    bool IsPPPFRegion();
    bool IsPPFFRegion();
    void SetZ1Z2Masses();
    float getEventWeight(Long64_t entry);
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif
