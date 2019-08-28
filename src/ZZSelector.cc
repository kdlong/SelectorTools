#include "Analysis/VVAnalysis/interface/ZZSelector.h"
#include "TLorentzVector.h"
#include <boost/algorithm/string.hpp>

void ZZSelector::Init(TTree *tree)
{
    weight_info_ = 0;
    if (isMC_) {
        weight_info_ = GetLheWeightInfo();
    }
    systematics_ = {
        {electronRecoEffUp,"CMS_RecoEff_eUp"},
        {electronRecoEffDown,"CMS_RecoEff_eDown"},
        {electronEfficiencyUp, "CMS_eff_eUp"},
        {electronEfficiencyDown, "CMS_eff_eDown"},
        {muonEfficiencyUp, "CMS_eff_mUp"},
        {muonEfficiencyDown, "CMS_eff_mDown"},
        {pileupUp, "CMS_pileupUp"},
        {pileupDown, "CMS_pileupDown"},
    };

    systHists_ = {
        "yield",
        "Mass",
        "ZMass",
        "ZZPt",
        "Z1Pt",
        "Z2Pt"
        "ZPt",
        "LepPt",
        "Z1lep1_Pt",
        "Z1lep2_Pt",
        "Z2lep1_Pt",
        "Z2lep2_Pt",
        "nTruePU",
    };
    hists1D_ = {
         "yield", "backgroundControlYield","nTruePU","nvtx","ZMass","Z1Mass","Z2Mass","ZZPt",  
         "Z1Pt","Z2Pt","Z1Phi","Z2Phi","dPhiZ1Z2","ZPt","LepPt","LepEta","Lep12Pt","Lep12Eta",
         "Lep34Pt","Lep34Eta","Z1lep1_Eta","Z1lep1_Phi","Z1lep1_Pt","Z1lep1_PdgId","Z1lep2_Eta", 
         "Z1lep2_Phi","Z1lep2_Pt","Z1lep2_PdgId","Z2lep1_Eta","Z2lep1_Phi","Z2lep1_Pt","Z2lep1_PdgId",
         "Z2lep2_Eta","Z2lep2_Phi","Z2lep2_Pt","Z2lep2_PdgId","Mass","nJets",
    };

    ZZSelectorBase::Init(tree);
}

void ZZSelector::SetBranchesUWVV() {
    ZZSelectorBase::SetBranchesUWVV();
    if (isMC_) {
        if (weight_info_ > 0)
            fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
        if ((weight_info_ == 2 || weight_info_ == 3) && doSystematics_)
            fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
    }
    fChain->SetBranchAddress("Mass", &Mass, &b_Mass);
    fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
}

unsigned int ZZSelector::GetLheWeightInfo() {
    std::vector<std::string> noLheWeights = {
        "ggZZ2e2mu", "ggZZ4e", "ggZZ4m", "ggZZ4t", "ggZZ2e2tau","ggZZ2mu2tau"
    };
    std::vector<std::string> scaleAndPdfWeights = {
        "wz3lnu-powheg", "wz3lnu-mg5amcnlo",
        "ZZZ", "WZZ", "WWZ", 
        "zz4l-powheg", "zz4l-amcatnlo",
        "ZZJJTo4L-EWK"
    };
    std::vector<std::string> allLheWeights = {
        //"wzjj-aqgcft", "wzjj-aqgcfm", "wzjj-aqgcfs",
        //"wz-atgc_pt0-200", "wz-atgc_pt200-300",
        //"wz-atgc_pt300"
    };

    if (std::find(noLheWeights.begin(), noLheWeights.end(), name_) != noLheWeights.end())
        return 0;
    if (std::find(scaleAndPdfWeights.begin(), scaleAndPdfWeights.end(), name_) != scaleAndPdfWeights.end())
        return 2;
    if (std::find(allLheWeights.begin(), allLheWeights.end(), name_) != allLheWeights.end())
        return 3;
    return 1;
}
void ZZSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) {

    ZZSelectorBase::LoadBranchesUWVV(entry, variation);

    //b_MtToMET->GetEntry(entry);
    //b_l1Pt->GetEntry(entry);
    //b_l2Pt->GetEntry(entry);
    //b_l3Pt->GetEntry(entry);
    b_Mass->GetEntry(entry);
    b_Pt->GetEntry(entry);
    //std::cout<<"channel in LoadBranches function: "<<channel_<<std::endl;
    if(channel_ == eemm || channel_ == mmee){
      //if(TightZZLeptons()){//i don't think this condition is needed even though it might save time but it messes up sf application for CRs in eemm,mmee states
      if(TightZZLeptons()){
        SetVariables(entry);}
    } 
    if (isMC_) {
        ApplyScaleFactors();
      }
    if (variation.first == Central) {
        if (isMC_ && doSystematics_) {
            if (isMC_ && weight_info_ > 0) {
                b_scaleWeights->GetEntry(entry);
                lheWeights = *scaleWeights;
                if (weight_info_ == 2) {
                    b_pdfWeights->GetEntry(entry);
                    // Only keep NNPDF weights
                    lheWeights.insert(lheWeights.end(), pdfWeights->begin(), 
                        pdfWeights->begin()+std::min(static_cast<size_t>(103), pdfWeights->size()));
                }
                else if (weight_info_ == 3) {
                    b_pdfWeights->GetEntry(entry);
                    lheWeights.insert(lheWeights.end(), pdfWeights->begin(), pdfWeights->end());
                }
            }
        }
      }
    else if (isMC_) {
      //Systematic uncertainties and creating shiftUp and shiftDown histograms
      //Starting with lepton Efficiencies
        if (variation.first == electronEfficiencyUp || variation.first == electronEfficiencyDown || variation.first == electronRecoEffUp || variation.first == electronRecoEffDown || 
            variation.first == muonEfficiencyUp || variation.first == muonEfficiencyDown) {
              ShiftEfficiencies(variation.first);
        }
        else if (variation.first == pileupUp) {
             weight *= pileupSF_->Evaluate1D(nTruePU, ScaleFactor::ShiftUp)/pileupSF_->Evaluate1D(nTruePU);
        }
        else if (variation.first == pileupDown) {
            weight *= pileupSF_->Evaluate1D(nTruePU, ScaleFactor::ShiftDown)/pileupSF_->Evaluate1D(nTruePU);
        }
    }

    auto deltaPhiZZ = [](float phi1, float phi2) {
      float pi = TMath::Pi();
      float dphi = fabs(phi1-phi2);
      if(dphi>pi)
          dphi = 2.0*pi - dphi;
      return dphi;
    };

    dPhiZZ = deltaPhiZZ(Z1Phi,Z2Phi);
}

void ZZSelector::ApplyScaleFactors() {
    //In order to get around the Overflow issue, set the Pt, not ideal.
    //std::cout<<"weight before SF: "<<weight<<std::endl;
    if (channel_ == eeee) {
        float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
        float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
        float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
        float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
        if (eIdSF_ != nullptr and eGapIdSF_ != nullptr ) {
          //In order to get around the Overflow issue, set the Pt, not ideal.
          if(l1IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
            //std::cout<<"weight for GapE: "<<weight<<std::endl;
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
          }
          if(l2IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
          }
          if(l3IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
          }
          if(l4IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
          }
        }
        if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr ) {
          //Applying Electron Reco SFs
          if(l1Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
          }
          if(l2Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
          }
          if(l3Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
          }
          if(l4Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
          }
      }
    }
    else if (channel_ == eemm) {
        float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
        float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
        float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
        float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
        if (eIdSF_ != nullptr and eGapIdSF_ != nullptr ) {
          //In order to get around the Overflow issue, set the Pt, not ideal.
          if(l1IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
          }
          if(l2IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
          }
        }
        //Applying Electron Reco SFs
        if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr ) {
          if(l1Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
          }
          if(l2Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
          }
        }
        if (mIdSF_ != nullptr) {
          weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
          weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
        }
    }
    else if (channel_ == mmee) {
        float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
        float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
        float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
        float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
        if (mIdSF_ != nullptr) {
          weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
          weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
        }
        if (eIdSF_ != nullptr and eGapIdSF_ != nullptr ) {
          if(l3IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
          }
          if(l4IsGap){
            weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
          }
          else{
            weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
          }
        }
        //Applying Electron Reco SFs
        if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr ) {
          if(l3Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
          }
          if(l4Pt < 20){ 
            weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
          }
          else{ 
            weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
          }
        }
    }
    else {
        float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
        float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
        float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
        float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
        if (mIdSF_ != nullptr) {
          weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
          weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
          weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
          weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
        }
    }
    if (pileupSF_ != nullptr) {
      weight *= pileupSF_->Evaluate1D(nTruePU);
    }
}
//Similar to Kenneth's SetShiftedMasses function which i will need later as well
void ZZSelector::SetVariables(Long64_t entry) {
    if(!(e1e2IsZ1(entry))){
      //std::cout<<"e1e2IsZ1 is working"<<std::endl;
      float tempMass = Z1mass;
      Z1mass = Z2mass;
      Z2mass = tempMass;
      float tempPt = Z1pt;
      Z1pt = Z2pt;
      Z2pt = tempPt;
      bool templ1IsTight = l1IsTight;
      l1IsTight = l3IsTight;
      l3IsTight = templ1IsTight;
      bool templ2IsTight = l2IsTight;
      l2IsTight = l4IsTight;
      l4IsTight = templ2IsTight;
      bool templ1IsIso = l1IsIso;
      l1IsIso = l3IsIso;
      l3IsIso = templ1IsIso;
      bool templ2IsIso = l2IsIso;
      l2IsIso = l4IsIso;
      l4IsIso = templ2IsIso;
      bool templ1IsGap = l1IsGap;
      l1IsGap = l3IsGap;
      l3IsGap = templ1IsGap;
      bool templ2IsGap = l2IsGap;
      l2IsGap = l4IsGap;
      l4IsGap = templ2IsGap;
      float templ1Pt = l1Pt;
      l1Pt = l3Pt;
      l3Pt = templ1Pt;
      float templ2Pt = l2Pt;
      l2Pt = l4Pt;
      l4Pt = templ2Pt;
      float templ1Eta = l1Eta;
      l1Eta = l3Eta;
      l3Eta = templ1Eta;
      float templ2Eta = l2Eta;
      l2Eta = l4Eta;
      l4Eta = templ2Eta;
      float templ1Phi = l1Phi;
      l1Phi = l3Phi;
      l3Phi = templ1Phi;
      float templ2Phi = l2Phi;
      l2Phi = l4Phi;
      l4Phi = templ2Phi;
      float templ1SIP3D = l1SIP3D;
      l1SIP3D = l3SIP3D;
      l3SIP3D = templ1SIP3D;
      float templ2SIP3D = l2SIP3D;
      l2SIP3D = l4SIP3D;
      l4SIP3D = templ2SIP3D;
      int templ1PdgId = l1PdgId;
      l1PdgId = l3PdgId;
      l3PdgId = templ1PdgId;
      int templ2PdgId = l2PdgId;
      l2PdgId = l4PdgId;
      l4PdgId = templ2PdgId;
      float templ1Mass = l1Mass;
      l1Mass = l3Mass;
      l3Mass = templ1Mass;
      float templ2Mass = l2Mass;
      l2Mass = l4Mass;
      l4Mass = templ2Mass;
    }
}
void ZZSelector::ShiftEfficiencies(Systematic variation) {
    ScaleFactor::Variation shift = ScaleFactor::Variation::ShiftUp;
    if (variation == electronEfficiencyDown || variation == electronRecoEffDown || variation == muonEfficiencyDown)
        shift = ScaleFactor::Variation::ShiftDown;

    if (channel_ == eeee){ 
      float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
      float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
      float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
      float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
      if (variation == electronRecoEffUp || variation == electronRecoEffDown) {
        //Applying Electron Reco SFs Up/Down for ElectronRecoSyst
        if(l1Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1,shift)/eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        if(l2Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2,shift)/eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        if(l3Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3,shift)/eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3) ;
        }
        if(l4Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4,shift)/eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
      }
      //Applying Electron ID SFs Up/Down for ElectronIDEffSyst
      else if (variation == electronEfficiencyUp || variation == electronEfficiencyDown) {
        if(l1IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift)/eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift)/eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        if(l2IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift)/eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift)/eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        if(l3IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift)/eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift)/eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        if(l4IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift)/eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift)/eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
    }
  }//channel eeee
    else if (channel_ == eemm) {
      //In order to get around the Overflow issue, set the Pt, not ideal.
      float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
      float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
      float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
      float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
      if (variation == electronRecoEffUp || variation == electronRecoEffDown) {
        //Applying Electron Reco SFs Up/Down for ElectronRecoEffSyst
        if(l1Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1,shift)/eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        if(l2Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2,shift)/eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
      }
      //Applying Electron ID SFs Up/Down for ElectronIDEffSyst
      else if (variation == electronEfficiencyUp || variation == electronEfficiencyDown) {
        if(l1IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift)/eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift)/eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        if(l2IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift)/eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift)/eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
      }
      else if (variation == muonEfficiencyUp || variation == muonEfficiencyDown) {
          weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3, shift)/mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
          weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4, shift)/mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
      }
    }
    else if (channel_ == mmee) {
      float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
      float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
      float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
      float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
      if (variation == muonEfficiencyUp || variation == muonEfficiencyDown) {
          weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1, shift)/mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
          weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2, shift)/mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
      }
      else if (variation == electronRecoEffUp || variation == electronRecoEffDown) {
        //Applying Electron Reco SFs Up/Down for ElectronRecoEffSyst 
        if(l3Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3,shift)/eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3) ;
        }
        if(l4Pt < 20){ 
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt,shift)/eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
        }
        else{ 
          weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4,shift)/eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
      }
      //Applying Electron ID SFs Up/Down for ElectronIDEffSyst
      else if (variation == electronEfficiencyUp || variation == electronEfficiencyDown) {
        if(l3IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift)/eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift)/eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        if(l4IsGap){
          weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift)/eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
        else{
          weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift)/eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
      }
    }
    else if (channel_ == mmmm && (variation == muonEfficiencyUp || variation == muonEfficiencyDown)) {
      float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
      float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
      float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
      float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
      weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1, shift)/mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
      //std::cout<<"l1 Shift: "<< mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1, shift)<<std::endl;
      //std::cout<<"l1 sf: "<< mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1)<<std::endl;
      weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2, shift)/mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
      weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3, shift)/mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
      //std::cout<<"l3 Shift: "<< mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3, shift)<<std::endl;
      //std::cout<<"l3 sf: "<< mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3)<<std::endl;
      weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4, shift)/mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
    }
}
bool ZZSelector::PassesZZSelection(){
  if (ZZSelection() && TightZZLeptons())
    return true;
  else
    return false;
}
bool ZZSelector::PassesHZZSelection(){
  if (ZSelection() && TightZZLeptons())
    return true;
  else
    return false;
}
bool ZZSelector::TightZZLeptons() {
    if (tightZ1Leptons() && tightZ2Leptons())
        return true;
    else
        return false;
}
bool ZZSelector::ZZSelection() {
    if ((Z1mass > 60.0 && Z1mass < 120.0) && (Z2mass > 60.0 && Z2mass < 120.0))
        return true;
    else
        return false;
}

//We already require 4 < Z1,Z2 < 120  in the "Loose Skim"
bool ZZSelector::ZSelection() {
    if (Z1mass > 40.0 && Z2mass > 12.0)
        return true;
    else
        return false;
}
bool ZZSelector::Z4lSelection() {
    if (Mass > 80.0 && Mass < 100.0)
        return true;
    else
        return false;
}
//This is no longer needed as I apply SIP < 4.0 on electrons while skimming and we don't need SIP cut on muons with latest HZZID
bool ZZSelector::HZZSIPSelection(){
    if ((l1SIP3D < 4.0 && l2SIP3D < 4.0 && l3SIP3D < 4.0 && l4SIP3D < 4.0))
        return true;
    else
        return false;
}
bool ZZSelector::HZZLowMass() {
    if (Mass > 70.0 && Mass < 110.0)
        return true;
    else
        return false;
}
bool ZZSelector::HZZMediumMass() {
    if (Mass > 138.0 && Mass < 300.0)
        return true;
    else
        return false;
}
bool ZZSelector::TestMuons(){
    if ((Z1mass > 82.0 && Z1mass < 102.0) && (Z2mass < 40.0))
        return true;
    else
        return false;
}

void ZZSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    //bool noBlind = true;
    //Applying the ZZ Selection here
    //std::cout<<"Is fillHistograms working?"<<std::endl;
    //std::cout<<"weight:"<<weight<<std::endl;
    if (!PassesZZSelection()){
        //std::cout<<"It should not be doing this?"<<std::endl;
        return;}
    if ((variation.first == Central || (doaTGC_ && isaTGC_)) && isMC_){
        for (size_t i = 0; i < lheWeights.size(); i++) {
            SafeHistFill(weighthists_, getHistName("yield", variation.second), 1, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthists_, getHistName("Mass", variation.second), Mass, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthists_, getHistName("Pt", variation.second), Pt, i, lheWeights[i]/lheWeights[0]*weight);
        }
      }
    SafeHistFill(histMap1D_, getHistName("yield", variation.second), 1, weight);
    SafeHistFill(histMap1D_, getHistName("Mass", variation.second), Mass,weight);
    SafeHistFill(histMap1D_, getHistName("ZMass", variation.second), Z1mass, weight);
    SafeHistFill(histMap1D_, getHistName("ZMass", variation.second), Z2mass, weight);
    //Making LeptonPt and Eta plots
    SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l1Pt, weight);
    SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l2Pt, weight);
    SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l3Pt, weight);
    SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l4Pt, weight);
    SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l1Eta, weight);
    SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l2Eta, weight);
    SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l3Eta, weight);
    SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l4Eta, weight);
    // Summing 12,34 leptons
    SafeHistFill(histMap1D_, getHistName("Lep12Pt", variation.second), l1Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Lep12Pt", variation.second), l2Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Lep34Pt", variation.second), l3Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Lep34Pt", variation.second), l4Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Lep12Eta", variation.second), l1Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Lep12Eta", variation.second), l2Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Le34Eta", variation.second), l3Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Lep34Eta", variation.second), l4Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Z1Mass", variation.second), Z1mass, weight);
    SafeHistFill(histMap1D_, getHistName("Z2Mass", variation.second), Z2mass, weight);
    SafeHistFill(histMap1D_, getHistName("ZPt", variation.second), Z1pt, weight);
    SafeHistFill(histMap1D_, getHistName("ZPt", variation.second), Z2pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z1Pt", variation.second), Z1pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z2Pt", variation.second), Z2pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z1Phi", variation.second), Z1Phi, weight);
    SafeHistFill(histMap1D_, getHistName("Z2Phi", variation.second), Z2Phi, weight);
    SafeHistFill(histMap1D_, getHistName("dPhiZ1Z2", variation.second), dPhiZZ, weight);
    SafeHistFill(histMap1D_, getHistName("ZZPt", variation.second), Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep1_Pt", variation.second), l1Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep1_Eta", variation.second), l1Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep1_Phi", variation.second), l1Phi, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep1_PdgId", variation.second), l1PdgId, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep2_Pt", variation.second), l2Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep2_Eta", variation.second), l2Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep2_Phi", variation.second), l2Phi, weight);
    SafeHistFill(histMap1D_, getHistName("Z1lep2_PdgId", variation.second), l2PdgId, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep1_Pt", variation.second), l3Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep1_Eta", variation.second), l3Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep1_Phi", variation.second), l3Phi, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep1_PdgId", variation.second), l3PdgId, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep2_Pt", variation.second), l4Pt, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep2_Eta", variation.second), l4Eta, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep2_Phi", variation.second), l4Phi, weight);
    SafeHistFill(histMap1D_, getHistName("Z2lep2_PdgId", variation.second), l4PdgId, weight);
    SafeHistFill(hists2D_, getHistName("Z1lep1_Z1lep2_Pt",variation.second),l1Pt,l2Pt,weight);
    SafeHistFill(hists2D_, getHistName("Z1lep1_Z1lep2_Eta",variation.second),l1Eta,l2Eta,weight);
    SafeHistFill(hists2D_, getHistName("Z1lep1_Z1lep2_Phi",variation.second),l1Phi,l2Phi,weight);
    SafeHistFill(hists2D_, getHistName("Z2lep1_Z2lep2_Pt",variation.second),l3Pt,l4Pt,weight);
    SafeHistFill(hists2D_, getHistName("Z2lep1_Z2lep2_Eta",variation.second),l3Eta,l4Eta,weight);
    SafeHistFill(hists2D_, getHistName("Z2lep1_Z2lep2_Phi",variation.second),l3Phi,l4Phi,weight);
    //2D Z1 vs Z2
    SafeHistFill(hists2D_, getHistName("Z1Mass_Z2Mass",variation.second),Z1mass,Z2mass,weight);

    if (histMap1D_[getHistName("nvtx", variation.second)] != nullptr) {
        b_nvtx->GetEntry(entry);
        histMap1D_[getHistName("nvtx", variation.second)]->Fill(nvtx, weight);
    }
    if (isMC_)
      SafeHistFill(histMap1D_, getHistName("nTruePU", variation.second), nTruePU, weight);
}


void ZZSelector::SetupNewDirectory()
{
    SelectorBase::SetupNewDirectory();
    isaTGC_ = name_.find("atgc") != std::string::npos;
    //std::cout<<"selection in ZZSelector: "<<selection_<<std::endl;
    applyFullSelection_ = (selection_ == ZZselection);
    //std::cout<<applyFullSelection_<<std::endl;
    InitializeHistogramsFromConfig();   
    //std::cout<<"Do Histos get initialized"<<std::endl;
}
