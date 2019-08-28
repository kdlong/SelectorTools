#include "Analysis/VVAnalysis/interface/ZZBackgroundSelector.h"
#include "TLorentzVector.h"
#include <TStyle.h>

void ZZBackgroundSelector::SlaveBegin(TTree * /*tree*/)
{
    //systematics_ = {
    //    {jetEnergyScaleUp, "CMS_scale_jUp"}, 
    //    {jetEnergyScaleDown, "CMS_scale_jDown"}, 
    //    {jetEnergyResolutionUp, "CMS_res_jUp"},
    //    {jetEnergyResolutionDown, "CMS_res_jDown"},
    //};
    //isNonpromptEstimate_ = true;
    //doaQGC_ = false;
    ZZSelector::SlaveBegin(0);
    fakeRate_allE_ = (ScaleFactor *) GetInputList()->FindObject("fakeRate_allE");
    if (fakeRate_allE_ == nullptr ) Abort("Must pass electron fake rate to input list!");
    fakeRate_allMu_ = (ScaleFactor *) GetInputList()->FindObject("fakeRate_allMu");
    if (fakeRate_allMu_ == nullptr ) Abort("Must pass muon fake rate to input list!");
}

void ZZBackgroundSelector::SetupNewDirectory()
{
    ZZSelector::SetupNewDirectory();
   
    // Insure that hist ranges are exactly the same as ZZSelector, just change name
    for (const auto && obj : *currentHistDir_) {
        std::string name = obj->GetName();
        //std::cout<<"histname: " <<name<<std::endl;
        TNamed* named = dynamic_cast<TNamed*>(obj);
        named->SetName(name.insert(name.length()-4, "Fakes_").c_str());
        //std::cout<<"histname in BkgSelector: " <<named<<std::endl;
    } 
      AddObject<TH1D>(Z1MassHistPPPF_, ("Z1Mass_PPPF_"+channelName_).c_str(), "Z1Mass; m_{Z_{1}} [GeV]; Events;", 60, 0, 120);
      AddObject<TH1D>(Z1MassHistPPFF_, ("Z1Mass_PPFF_"+channelName_).c_str(), "Z1Mass; m_{Z_{1}} [GeV]; Events;", 60, 0, 120);
      
      AddObject<TH1D>(Z2MassHistPPPF_, ("Z2Mass_PPPF_"+channelName_).c_str(), "Z2Mass; m_{Z_{2}} [GeV]; Events;", 60, 0, 120);
      AddObject<TH1D>(Z2MassHistPPFF_, ("Z2Mass_PPFF_"+channelName_).c_str(), "Z2Mass; m_{Z_{2}} [GeV]; Events;", 60, 0, 120);

      AddObject<TH1D>(MassHistPPPF_, ("Mass_PPPF_"+channelName_).c_str(), "Mass; m_{4l} [GeV]; Events;", 40, 70, 870);
      AddObject<TH1D>(MassHistPPFF_, ("Mass_PPFF_"+channelName_).c_str(), "Mass; m_{4l} [GeV]; Events;", 40, 70, 870);

      AddObject<TH1D>(WeightsHistmmee_, ("Weights_mmee_"+channelName_).c_str(), "Weight; Event Weight; Events;", 10, -5, 5);
      AddObject<TH1D>(WeightsHisteemm_, ("Weights_eemm_"+channelName_).c_str(), "Weight; Event Weight; Events;", 100, -5, 5);
}
float ZZBackgroundSelector::getEventWeight(Long64_t entry) {
    float evtwgt = 0;

    if(channel_ == mmee && (e1e2IsZ1(entry)))
        WeightsHistmmee_->Fill(1,weight);

    if(channel_ == eemm && !(e1e2IsZ1(entry)))
        WeightsHisteemm_->Fill(1,weight);

    //if ((channel_ == eeee) || (channel_ == mmmm)){
    //}
    //in eemm the e1e2IsZ1 function already sets the correct Z1,Z2 leptons and even the IDs from LoadBranches
    //So only need to take care which is l3 and l4
   // if ((channel_ == eemm) || (channel_ == mmee)){
   //   if(Z2FP()){
   //     float templ3Pt = l3Pt;
   //     l3Pt = l4Pt;
   //     l4Pt = templ3Pt;
   //     float templ3Eta = l3Eta;
   //     l3Eta = l4Eta;
   //     l4Eta = templ3Eta;
   //   }
   // }
    //std::cout<<"Weight in Bkg Seletor getEventWeight function: "<<weight<<std::endl;
      if (IsPPPFRegion()) {
        if (true){
             //std::cout<<"Weight in PPPF: "<<weight<<std::endl;
             Z1MassHistPPPF_->Fill(Z1mass, weight);
             Z2MassHistPPPF_->Fill(Z2mass, weight);
             MassHistPPPF_->Fill(Mass,weight);
         }
         evtwgt = (getl4FakeRate(entry)*weight);
         //WeightsHistPPPF_->Fill(1,evtwgt);
      }
      if (IsPPFFRegion()) {
        if (true) {
            //std::cout<<"Weight in PPFF: "<<weight<<std::endl;
            Z1MassHistPPFF_->Fill(Z1mass, weight);
            Z2MassHistPPFF_->Fill(Z2mass, weight);
            MassHistPPFF_->Fill(Mass,weight);
        }
        evtwgt = ((-1*getl3FakeRate(entry)*getl4FakeRate(entry))*weight);
        //evtwgt = ((getl3FakeRate(entry)*getl4FakeRate(entry))*weight);
        //WeightsHistPPFF_->Fill(1,evtwgt);
    }
    //std::cout<<"eventWeight after Fake Rate: "<<evtwgt<<std::endl;
    return evtwgt;
}

void ZZBackgroundSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) {
    ZZSelector::LoadBranches(entry, variation);
    SetZ1Z2Masses();
    weight = getEventWeight(entry);
}
float ZZBackgroundSelector::getl3FakeRate(Long64_t entry) {
    float pt_fillval = l3Pt < FR_MAX_PT_ ? l3Pt : FR_MAX_PT_ - 0.01;
    float fr = 1;
    if ((channel_ == eeee) || (channel_ == mmee && !(e1e2IsZ1(entry)))){
        fr = fakeRate_allE_->Evaluate2D(pt_fillval, std::abs(l3Eta));
        //fr = 0.03;//avg e fake rate
        //std::cout<<"channel: "<<channel_<<std::endl;
        //std::cout<<"l3 E Fake Rate: "<<fr<<std::endl;
    }
    else if ((channel_ == mmmm) || (channel_ == eemm && (e1e2IsZ1(entry)))){
        fr = fakeRate_allMu_->Evaluate2D(pt_fillval, std::abs(l3Eta));
        //fr = 0.1;//avg mu fake rate
        //std::cout<<"channel: "<<channel_<<std::endl;
        //std::cout<<"l3 Mu Fake Rate: "<<fr<<std::endl;
    }
    else{
      fr=0.;
    }
    return fr/(1-fr); 
}

float ZZBackgroundSelector::getl4FakeRate(Long64_t entry) {
    float pt_fillval = l4Pt < FR_MAX_PT_ ? l4Pt : FR_MAX_PT_ - 0.01;
    float fr = 1;
    if ((channel_ == eeee) || (channel_ == mmee && !(e1e2IsZ1(entry)))){
        fr = fakeRate_allE_->Evaluate2D(pt_fillval, std::abs(l4Eta));
        //fr = 0.03;//avg e fake rate
        //std::cout<<"channel: "<<channel_<<std::endl;
        //std::cout<<"l4 E Fake Rate: "<<fr<<std::endl;
    }
    else if ((channel_ == mmmm) || (channel_ == eemm && (e1e2IsZ1(entry)))){
        fr = fakeRate_allMu_->Evaluate2D(pt_fillval, std::abs(l4Eta));
        //fr = 0.1;//avg mu fake rate
        //std::cout<<"channel: "<<channel_<<std::endl;
        //std::cout<<"l4 Mu Fake Rate: "<<fr<<std::endl;
    }
    else{
      fr=0.;
    }
    return fr/(1-fr); 
}
//Remember that we only build Z1 (Real Z) out of OS-SF tight leptons
bool ZZBackgroundSelector::IsPPPFRegion() {
    return ((tightZ1Leptons() && (Z2PF() || Z2FP())) || (tightZ2Leptons() && (Z1PF() || Z1FP())));
}
//Remember that we only build Z1 (Real Z) out of OS-SF tight leptons
bool ZZBackgroundSelector::IsPPFFRegion() {
    if((channel_ == eeee) || (channel_ == mmmm))
      return ((tightZ1Leptons() && Z2FF()) || (tightZ2Leptons() && Z1FF()) || (Z1FP() && Z2PF()) || (Z1PF() && Z2FP()));
    else
      return ((tightZ1Leptons() && Z2FF()) || (tightZ2Leptons() && Z1FF()));
}
//Remember that we only build Z1 (Real Z) out of OS-SF tight leptons
//This is only for 4e,4mu channels as eemm/mmee channels is already taken care of by SetVariables function which gets
//called in LoadBranches() from ZZSelector. The SetVariables function should work in CRs, SR equally. We only need
//SetZ1Z2Masses function for 4e/4mu in CRs.
void ZZBackgroundSelector::SetZ1Z2Masses() {
    TLorentzVector lepton1;
    lepton1.SetPtEtaPhiM(l1Pt, l1Eta, l1Phi, l1Mass);
    TLorentzVector lepton2;
    lepton2.SetPtEtaPhiM(l2Pt, l2Eta, l2Phi, l2Mass);
    TLorentzVector lepton3;
    lepton3.SetPtEtaPhiM(l3Pt, l3Eta, l3Phi, l3Mass);
    TLorentzVector lepton4;
    lepton4.SetPtEtaPhiM(l4Pt, l4Eta, l4Phi, l4Mass);
    if(tightZ1Leptons() && !tightZ2Leptons()){
      Z1mass = (lepton1+lepton2).M();
      Z2mass = (lepton3+lepton4).M();
      Z1pt = (lepton1+lepton2).Pt();
      Z2pt = (lepton3+lepton4).Pt();
      //In Z2 what is l3 and l4 can change the fake rate a little bit.
      if(Z2FP()){
        float templ3Pt = l3Pt;
        l3Pt = l4Pt;
        l4Pt = templ3Pt;
        float templ3Eta = l3Eta;
        l3Eta = l4Eta;
        l4Eta = templ3Eta;
      }
    }
    else if(tightZ2Leptons() && !tightZ1Leptons()){  
      Z1mass = (lepton3+lepton4).M();
      Z2mass = (lepton1+lepton2).M();
      Z1pt = (lepton3+lepton4).Pt();
      Z2pt = (lepton1+lepton2).Pt();
      //Fakes are l1,l2 from skims, reverse them
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
      //Now we have two fakes identified by l3Pt, l4Pt and l3Eta, l4Eta
      //Further special condition between l3,l4 which one to use for l4fake rate in PPPF region, their IDs still are labeled l1IsTight,l2IsTight
      if(Z1FP()){
        float templ3Pt = l3Pt;
        l3Pt = l4Pt;
        l4Pt = templ3Pt;
        float templ3Eta = l3Eta;
        l3Eta = l4Eta;
        l4Eta = templ3Eta;
      }
        }
    //The last two conditions only matter for TTJets fakes (very small amount)
    else if(Z1FP() && Z2PF()){
      //Make sure I am not making a Z in eemm with an e and mu!
      if ((channel_ == eeee) || (channel_ == mmmm)){
        Z1mass = (lepton2+lepton3).M();
        Z2mass = (lepton1+lepton4).M();
        Z1pt = (lepton2+lepton3).Pt();
        Z2pt = (lepton1+lepton4).Pt();
        //Here the two fakes are l1,l4 and we only need to relabel l1 -> l3
        float templ1Pt = l1Pt;
        l1Pt = l3Pt;
        l3Pt = templ1Pt;
        float templ1Eta = l1Eta;
        l1Eta = l3Eta;
        l3Eta = templ1Eta;}
      else{
      }
    }
    else if(Z1PF() && Z2FP()){
      //Make sure I am not making a Z in eemm with an e and mu!
      if ((channel_ == eeee) || (channel_ == mmmm)){
        Z1mass = (lepton1+lepton4).M();
        Z2mass = (lepton2+lepton3).M();
        Z1pt = (lepton1+lepton4).Pt();
        Z2pt = (lepton2+lepton3).Pt();
        //Here the two fakes are l2,l3 and we only need to relabel l2 -> l4 since this only matters in PPFF region so l3,l4 are interchangeable
        float templ2Pt = l2Pt;
        l2Pt = l4Pt;
        l4Pt = templ2Pt;
        float templ2Eta = l2Eta;
        l2Eta = l4Eta;
        l4Eta = templ2Eta;
      }
      else{
      }
    }
}

