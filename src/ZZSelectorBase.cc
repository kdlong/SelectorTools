#include "Analysis/VVAnalysis/interface/ZZSelectorBase.h"
#include <TStyle.h>
#include <regex>
#include "TParameter.h"

std::string ZZSelectorBase::GetNameFromFile() {
    std::regex expr = std::regex("201[0-9]-[0-9][0-9]-[0-9][0-9]-(.*)-ZZ4l2018");
    std::smatch matches;
    std::string fileName = fChain->GetTree()->GetDirectory()->GetFile()->GetName(); 

    std::regex_search(fileName, matches, expr);
    return std::string(matches.str(1));
}
void ZZSelectorBase::SetScaleFactors()
{
    pileupSF_ = (ScaleFactor *) GetInputList()->FindObject("pileupSF");
    if (pileupSF_ == nullptr ) 
        std::invalid_argument("Must pass pileup weights SF");
    eLowRecoSF_ = (ScaleFactor *) GetInputList()->FindObject("electronLowReco19SF");
    if (eLowRecoSF_ == nullptr ) 
        std::invalid_argument("Must pass electron LowReco18 SF");
    eRecoSF_ = (ScaleFactor *) GetInputList()->FindObject("electronReco19SF");
    if (eRecoSF_ == nullptr ) 
        std::invalid_argument("Must pass electron Reco18 SF");
    eIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronMoriond19SF");
    if (eIdSF_ == nullptr ) 
        std::invalid_argument("Must pass electron Moriond19 SF");
    eGapIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronMoriond19GapSF");
    if (eGapIdSF_ == nullptr ) 
        std::invalid_argument("Must pass electronGap Moriond19 SF");
    mIdSF_ = (ScaleFactor *) GetInputList()->FindObject("muonMoriond19SF");
    if (mIdSF_ == nullptr ) 
        std::invalid_argument("Must pass muon Moriond19 SF");

    //There are L1Prefiring weight and uncertainity in the ZZ UWVV ntuples 
    //prefireEff_ = (TEfficiency*) GetInputList()->FindObject("prefireEfficiencyMap");
    //if (prefireEff_ == nullptr ) 
    //    std::invalid_argument("Must pass prefiring efficiency map");
}

void ZZSelectorBase::Init(TTree *tree)
{

    //allChannels_ = {"eeee","eemm","mmee","mmmm","eee", "eem", "emm", "mmm"};

    SelectorBase::Init(tree);
} 

void ZZSelectorBase::SetBranchesUWVV() {
    if (isMC_){
        fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
        fChain->SetBranchAddress("nTruePU", &nTruePU, &b_nTruePU);
    }
    else {
        //fChain->SetBranchAddress("Flag_duplicateMuonsPass", Flag_duplicateMuonsPass);
        //fChain->SetBranchAddress("Flag_badMuonsPass", Flag_badMuonsPass);
    }
    //std::cout<<"In ZZSelectorBase"<<std::endl;
    //std::cout<<"channelName: "<<channelName_<<std::endl;
    //std::cout<<"enum channel_: "<<channel_<<std::endl;
    //std::cout<<"isMC: "<<isMC_<<std::endl;
    //if (channelName_ == "eeee") {
    if (channel_ == eeee){
        //std::cout<<"enum channel_: "<<channel_<<std::endl;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        //fChain->SetBranchAddress("e1IsFall17isoV2wpHZZ", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        //fChain->SetBranchAddress("e2IsFall17isoV2wpHZZ", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        //fChain->SetBranchAddress("e3IsFall17isoV2wpHZZ", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("e3ZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("e3ZZIsoPass", &l3IsIso, &b_l3IsIso);
        //fChain->SetBranchAddress("e4IsFall17isoV2wpHZZ", &l4IsTight, &b_l4IsTight);
        fChain->SetBranchAddress("e4ZZTightIDNoVtx", &l4IsTight, &b_l4IsTight);
        fChain->SetBranchAddress("e4ZZIsoPass", &l4IsIso, &b_l4IsIso);
        fChain->SetBranchAddress("e1_e2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("e3_e4_Mass", &Z2mass, &b_Z2mass);
        //fChain->SetBranchAddress("e1_e3_Mass", &Zamass, &b_Zamass);
        //fChain->SetBranchAddress("e2_e4_Mass", &Zbmass, &b_Zbmass);
        fChain->SetBranchAddress("e1_e2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("e3_e4_Pt", &Z2pt, &b_Z2pt);
        fChain->SetBranchAddress("e1_e2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("e3_e4_Phi", &Z2Phi, &b_Z2Phi);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("e3Pt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("e4Pt", &l4Pt, &b_l4Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("e3Eta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("e4Eta", &l4Eta, &b_l4Eta);
        fChain->SetBranchAddress("e1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("e2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("e3Phi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("e4Phi", &l4Phi, &b_l4Phi);
        fChain->SetBranchAddress("e1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("e2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("e3SIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("e4SIP3D", &l4SIP3D, &b_l4SIP3D);
        fChain->SetBranchAddress("e1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("e2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("e3PdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("e4PdgId", &l4PdgId, &b_l4PdgId);
        fChain->SetBranchAddress("e1Mass", &l1Mass, &b_l1Mass);
        fChain->SetBranchAddress("e2Mass", &l2Mass, &b_l2Mass);
        fChain->SetBranchAddress("e3Mass", &l3Mass, &b_l3Mass);
        fChain->SetBranchAddress("e4Mass", &l4Mass, &b_l4Mass);
        fChain->SetBranchAddress("e1IsGap", &l1IsGap, &b_l1IsGap);
        fChain->SetBranchAddress("e2IsGap", &l2IsGap, &b_l2IsGap);
        fChain->SetBranchAddress("e3IsGap", &l3IsGap, &b_l3IsGap);
        fChain->SetBranchAddress("e4IsGap", &l4IsGap, &b_l4IsGap);
    }
    //Add 2e2mu channel also but it still needs to differentiate which one is Z1Mass and which one is Z2Mass leptons
    //This is done with a flag at the time of Process for each event on the fly
    else if (channel_ == eemm) {
        //channel_ = eemm;
        //fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight); 
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        //fChain->SetBranchAddress("e1IsFall17isoV2wpHZZ", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        //fChain->SetBranchAddress("e2IsFall17isoV2wpHZZ", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("m1ZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("m1ZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("m2ZZTightIDNoVtx", &l4IsTight, &b_l4IsTight);
        fChain->SetBranchAddress("m2ZZIsoPass", &l4IsIso, &b_l4IsIso);
        fChain->SetBranchAddress("e1_e2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("m1_m2_Mass", &Z2mass, &b_Z2mass);
        fChain->SetBranchAddress("e1_e2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("m1_m2_Pt", &Z2pt, &b_Z2pt);
        fChain->SetBranchAddress("e1_e2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("m1_m2_Phi", &Z2Phi, &b_Z2Phi);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m1Pt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("m2Pt", &l4Pt, &b_l4Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("m1Eta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("m2Eta", &l4Eta, &b_l4Eta);
        fChain->SetBranchAddress("e1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("e2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("e1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("e2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("m1Phi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("m2Phi", &l4Phi, &b_l4Phi);
        fChain->SetBranchAddress("m1SIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("m2SIP3D", &l4SIP3D, &b_l4SIP3D);
        fChain->SetBranchAddress("e1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("e2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("m1PdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("m2PdgId", &l4PdgId, &b_l4PdgId);
        fChain->SetBranchAddress("e1Mass", &l1Mass, &b_l1Mass);
        fChain->SetBranchAddress("e2Mass", &l2Mass, &b_l2Mass);
        fChain->SetBranchAddress("m1Mass", &l3Mass, &b_l3Mass);
        fChain->SetBranchAddress("m2Mass", &l4Mass, &b_l4Mass);
        fChain->SetBranchAddress("e1IsGap", &l1IsGap, &b_l1IsGap);
        fChain->SetBranchAddress("e2IsGap", &l2IsGap, &b_l2IsGap);
        //There are 2 different SFs for electrons depending on whether the electron is from the crack 
        //in ECAL crystals or not but we need a dummy for the muons? Makes life easier later
        fChain->SetBranchAddress("m1IsLoose", &l3IsGap, &b_l3IsGap);
        fChain->SetBranchAddress("m2IsLoose", &l4IsGap, &b_l4IsGap);
    }
    else if (channel_ == mmee) {
        //channel_ = mmee;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        //fChain->SetBranchAddress("e1IsFall17isoV2wpHZZ", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        //fChain->SetBranchAddress("e2IsFall17isoV2wpHZZ", &l2IsTight, &b_l2IsTight); 
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("m1ZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("m1ZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("m2ZZTightIDNoVtx", &l4IsTight, &b_l4IsTight);
        fChain->SetBranchAddress("m2ZZIsoPass", &l4IsIso, &b_l4IsIso);
        fChain->SetBranchAddress("e1_e2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("m1_m2_Mass", &Z2mass, &b_Z2mass);
        fChain->SetBranchAddress("e1_e2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("m1_m2_Pt", &Z2pt, &b_Z2pt);
        fChain->SetBranchAddress("e1_e2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("m1_m2_Phi", &Z2Phi, &b_Z2Phi);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m1Pt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("m2Pt", &l4Pt, &b_l4Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("m1Eta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("m2Eta", &l4Eta, &b_l4Eta);
        fChain->SetBranchAddress("e1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("e2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("e1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("e2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("m1Phi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("m2Phi", &l4Phi, &b_l4Phi);
        fChain->SetBranchAddress("m1SIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("m2SIP3D", &l4SIP3D, &b_l4SIP3D);
        fChain->SetBranchAddress("e1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("e2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("m1PdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("m2PdgId", &l4PdgId, &b_l4PdgId);
        fChain->SetBranchAddress("e1Mass", &l1Mass, &b_l1Mass);
        fChain->SetBranchAddress("e2Mass", &l2Mass, &b_l2Mass);
        fChain->SetBranchAddress("m1Mass", &l3Mass, &b_l3Mass);
        fChain->SetBranchAddress("m2Mass", &l4Mass, &b_l4Mass);
        fChain->SetBranchAddress("e1IsGap", &l1IsGap, &b_l1IsGap);
        fChain->SetBranchAddress("e2IsGap", &l2IsGap, &b_l2IsGap);
        fChain->SetBranchAddress("m1IsLoose", &l3IsGap, &b_l3IsGap);
        fChain->SetBranchAddress("m2IsLoose", &l4IsGap, &b_l4IsGap);
    }
    else if (channel_ == mmmm) {
        //channel_ = mmmm;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        fChain->SetBranchAddress("m1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("m1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        fChain->SetBranchAddress("m2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("m2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("m3ZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("m3ZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("m4ZZTightIDNoVtx", &l4IsTight, &b_l4IsTight);
        fChain->SetBranchAddress("m4ZZIsoPass", &l4IsIso, &b_l4IsIso);
        fChain->SetBranchAddress("m1_m2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("m3_m4_Mass", &Z2mass, &b_Z2mass);
        //fChain->SetBranchAddress("m1_m3_Mass", &Zamass, &b_Zamass);
        //fChain->SetBranchAddress("m2_m4_Mass", &Zbmass, &b_Zbmass);
        fChain->SetBranchAddress("m1_m2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("m3_m4_Pt", &Z2pt, &b_Z2pt);
        fChain->SetBranchAddress("m1_m2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("m3_m4_Phi", &Z2Phi, &b_Z2Phi);
        fChain->SetBranchAddress("m1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("m2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m3Pt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("m4Pt", &l4Pt, &b_l4Pt);
        fChain->SetBranchAddress("m1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("m2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("m3Eta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("m4Eta", &l4Eta, &b_l4Eta);
        fChain->SetBranchAddress("m1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("m2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("m3Phi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("m4Phi", &l4Phi, &b_l4Phi);
        fChain->SetBranchAddress("m1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("m2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("m3SIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("m4SIP3D", &l4SIP3D, &b_l4SIP3D);
        fChain->SetBranchAddress("m1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("m2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("m3PdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("m4PdgId", &l4PdgId, &b_l4PdgId);
        fChain->SetBranchAddress("m1Mass", &l1Mass, &b_l1Mass);
        fChain->SetBranchAddress("m2Mass", &l2Mass, &b_l2Mass);
        fChain->SetBranchAddress("m3Mass", &l3Mass, &b_l3Mass);
        fChain->SetBranchAddress("m4Mass", &l4Mass, &b_l4Mass);
    }
    else if (channel_ == eee) {
        //channel_ = eee;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("e3ZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("e3ZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("e1_e2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("e1_e2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("e1_e2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("e3Pt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("e3Eta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("e1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("e2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("e3Phi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("e1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("e2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("e3SIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("e1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("e2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("e3PdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("e3MtToMET", &l3MtToMET, &b_l3MtToMET);
    }
    else if (channel_ == eem) { 
        //channel_ = eem;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("mZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("mZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("e1_e2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("e1_e2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("e1_e2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("e1_e2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("e1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("e2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("mPt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("mEta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("mPhi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("e1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("e2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("mSIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("e1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("e2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("mPdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("mMtToMET", &l3MtToMET, &b_l3MtToMET);
    }
    else if (channel_ == emm) { 
        //channel_ = emm;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        fChain->SetBranchAddress("eZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("eZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("m1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("m1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        fChain->SetBranchAddress("m2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("m2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("m1_m2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("m1_m2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("m1_m2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("m1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("m2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("m2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("m1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("m2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("ePt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("eEta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("ePhi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("m1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("m2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("eSIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("m1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("m2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("ePdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("eMtToMET", &l3MtToMET, &b_l3MtToMET);
    }
    else if (channel_ == mmm) { 
        //channel_ = mmm;
        fChain->SetBranchAddress("run", &run, &b_run);
        fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
        fChain->SetBranchAddress("evt", &evt, &b_evt);
        //fChain->SetBranchAddress("duplicated", &duplicated, &b_duplicated);
        fChain->SetBranchAddress("m1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("m1ZZIsoPass", &l1IsIso, &b_l1IsIso);
        fChain->SetBranchAddress("m2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("m2ZZIsoPass", &l2IsIso, &b_l2IsIso);
        fChain->SetBranchAddress("m3ZZTightIDNoVtx", &l3IsTight, &b_l3IsTight);
        fChain->SetBranchAddress("m3ZZIsoPass", &l3IsIso, &b_l3IsIso);
        fChain->SetBranchAddress("m1_m2_Mass", &Z1mass, &b_Z1mass);
        fChain->SetBranchAddress("m1_m2_Pt", &Z1pt, &b_Z1pt);
        fChain->SetBranchAddress("m1_m2_Phi", &Z1Phi, &b_Z1Phi);
        fChain->SetBranchAddress("m1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("m2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m3Pt", &l3Pt, &b_l3Pt);
        fChain->SetBranchAddress("m1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("m2Eta", &l2Eta, &b_l2Eta);
        fChain->SetBranchAddress("m3Eta", &l3Eta, &b_l3Eta);
        fChain->SetBranchAddress("m1Phi", &l1Phi, &b_l1Phi);
        fChain->SetBranchAddress("m2Phi", &l2Phi, &b_l2Phi);
        fChain->SetBranchAddress("m3Phi", &l3Phi, &b_l3Phi);
        fChain->SetBranchAddress("m1SIP3D", &l1SIP3D, &b_l1SIP3D);
        fChain->SetBranchAddress("m2SIP3D", &l2SIP3D, &b_l2SIP3D);
        fChain->SetBranchAddress("m3SIP3D", &l3SIP3D, &b_l3SIP3D);
        fChain->SetBranchAddress("m1PdgId", &l1PdgId, &b_l1PdgId);
        fChain->SetBranchAddress("m2PdgId", &l2PdgId, &b_l2PdgId);
        fChain->SetBranchAddress("m3PdgId", &l3PdgId, &b_l3PdgId);
        fChain->SetBranchAddress("m3MtToMET", &l3MtToMET, &b_l3MtToMET);
    }
    else
        throw std::invalid_argument("Invalid channel choice in ZZSelectorBase!");

    fChain->SetBranchAddress("type1_pfMETEt", &type1_pfMETEt, &b_type1_pfMETEt);
}

void ZZSelectorBase::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) {
    throw std::domain_error("NanoAOD ntuples not supported for ZZSelector!");
}

void ZZSelectorBase::SetBranchesNanoAOD() {
    throw std::domain_error("NanoAOD ntuples not supported for ZZSelector!");
}

void ZZSelectorBase::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) {
    weight = 1;
    b_l1Pt->GetEntry(entry);
    b_l2Pt->GetEntry(entry);
    b_l3Pt->GetEntry(entry);
    b_l1Eta->GetEntry(entry);
    b_l2Eta->GetEntry(entry);
    b_l3Eta->GetEntry(entry);
    b_l1Phi->GetEntry(entry);
    b_l2Phi->GetEntry(entry);
    b_l3Phi->GetEntry(entry);
    b_l1SIP3D->GetEntry(entry);
    b_l2SIP3D->GetEntry(entry);
    b_l3SIP3D->GetEntry(entry);
    b_l1PdgId->GetEntry(entry);
    b_l2PdgId->GetEntry(entry);
    b_l3PdgId->GetEntry(entry);
    b_run->GetEntry(entry);
    b_lumi->GetEntry(entry);
    b_evt->GetEntry(entry);
    if (channel_ == eeee || channel_ == eemm || channel_ == mmee || channel_ == mmmm) {
      b_l4Pt->GetEntry(entry);
      b_l4Eta->GetEntry(entry);
      b_l4Phi->GetEntry(entry);
      b_l4SIP3D->GetEntry(entry);
      b_l4PdgId->GetEntry(entry);
      b_l1Mass->GetEntry(entry);
      b_l2Mass->GetEntry(entry);
      b_l3Mass->GetEntry(entry);
      b_l4Mass->GetEntry(entry);
      b_l4IsTight->GetEntry(entry);
      b_l4IsIso->GetEntry(entry);
      //b_l4PVDXY->GetEntry(entry);
      //b_l4PVDZ->GetEntry(entry);
      b_Z2mass->GetEntry(entry);
      b_Z2pt->GetEntry(entry);
      b_Z2Phi->GetEntry(entry);
    }
    b_Z1mass->GetEntry(entry);
    b_Z1pt->GetEntry(entry);
    b_Z1Phi->GetEntry(entry);
    b_type1_pfMETEt->GetEntry(entry); 
    b_l1IsTight->GetEntry(entry);
    b_l1IsIso->GetEntry(entry);
    b_l2IsTight->GetEntry(entry);
    b_l2IsIso->GetEntry(entry);
    b_l3IsTight->GetEntry(entry);
    b_l3IsIso->GetEntry(entry); 
    if (channel_ == eeee || channel_ == eemm || channel_ == mmee){
      b_l1IsGap->GetEntry(entry); 
      b_l2IsGap->GetEntry(entry); 
      b_l3IsGap->GetEntry(entry); 
      b_l4IsGap->GetEntry(entry);
    }
    if (channel_ == eee || channel_ == eem || channel_ == emm || channel_ == mmm) {
      b_l3MtToMET->GetEntry(entry);
    }
    //std::cout<<"IsMC: "<<isMC_<<std::endl;
    if (isMC_) {
        b_genWeight->GetEntry(entry);
        b_nTruePU->GetEntry(entry);
        weight = genWeight;
    }
    if(channel_ == mmee){
      if(e1e2IsZ1(entry))
        weight=0.0;
        //Makes weight 0 if Z1 is ee hence should not go in _mmee histos
    }
    else if(channel_ == eemm){
      if(!(e1e2IsZ1(entry)))
        weight=0.0;
        //Makes weight 0 if Z1 is mm hence should not go in _eemm 
    }
}

//This condition identifies the Z1 candidate
//Required for the 2e2mu state 
//Can I save some kind of flag to identify the order of leptons for each event
//I think this part might need to be added to the ZZSelector and BackgroundSelector and the condition checked at
//process time for each event.
bool ZZSelectorBase::e1e2IsZ1(Long64_t entry){
  //4P Signal region logic where I need to differentiate between two tight pairs and assign Z1 depending on which is closer to mZ
  if(tightZ1Leptons() && tightZ2Leptons()){
    if(fabs(Z1mass-91.1876) < fabs(Z2mass-91.1876))
      return true;
    else
      return false;
  }
  //In CRs it doesn't matter, the tight pair builds the Z and the other pair is X in Z+X.
  else if(tightZ1Leptons() && !tightZ2Leptons()){
    return true;}

  else if(!tightZ1Leptons() && tightZ2Leptons()){
    return false;}
  else
    return false;
}
// Meant to be a wrapper for the tight ID just in case it changes
// To be a function of multiple variables
bool ZZSelectorBase::lep1IsTight() {
    return (l1IsTight && l1IsIso); 
}

bool ZZSelectorBase::lep2IsTight() {
    return (l2IsTight && l2IsIso); 
}
bool ZZSelectorBase::tightSIP() {
    if ((l1SIP3D < 4.0 && l2SIP3D < 4.0 && l3SIP3D < 4.0))
        return true;
    else
        return false;
}
bool ZZSelectorBase::tightZ1Leptons() {
    return lep1IsTight() && lep2IsTight(); 
}
bool ZZSelectorBase::lep3IsTight() {
    return (l3IsTight && l3IsIso);
}
bool ZZSelectorBase::lep4IsTight() {
    return (l4IsTight && l4IsIso);
}
bool ZZSelectorBase::tightZ2Leptons() {
    return lep3IsTight() && lep4IsTight(); 
}
bool ZZSelectorBase::Z1PF(){
  if(lep1IsTight() && !lep2IsTight())
    return true;
  else
    return false;
}
bool ZZSelectorBase::Z1FP(){
  if(lep2IsTight() && !lep1IsTight())
    return true;
  else
    return false;
}
bool ZZSelectorBase::Z1FF(){
  if(!lep1IsTight() && !lep2IsTight())
    return true;
  else
    return false;
}
bool ZZSelectorBase::Z2PF(){
  if(lep3IsTight() && !lep4IsTight())
    return true;
  else
    return false;
}
bool ZZSelectorBase::Z2FP(){
  if(lep4IsTight() && !lep3IsTight())
    return true;
  else
    return false;
}
bool ZZSelectorBase::Z2FF(){
  if(!lep3IsTight() && !lep4IsTight())
    return true;
  else
    return false;
}
