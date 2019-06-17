#include "Analysis/VVAnalysis/interface/ZSelector.h"
#include "TLorentzVector.h"
#include <TStyle.h>
#include <regex>

void ZSelector::Init(TTree *tree)
{
    allChannels_ = {"ee", "mm", "Unknown"};
    hists1D_ = {"CutFlow", "ZMass", "ZEta", "yZ", "ZPt", "ptl1", "etal1", "ptl2", "etal2"};

    b.SetTree(tree);
    SelectorBase::Init(tree);
    
    singleLepton_ = false;
    if (!isMC_ && name_.find("Single") != std::string::npos)
        singleLepton_ = true;
}

void ZSelector::SetScaleFactors() {
    pileupSF_ = (ScaleFactor *) GetInputList()->FindObject("pileupSF");
    if (pileupSF_ == nullptr ) 
        std::invalid_argument("Must pass pileup weights SF");
    eIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronTightIdSF");
    if (eIdSF_ == nullptr ) 
        std::invalid_argument("Must pass electron ID SF");
    eGsfSF_ = (ScaleFactor *) GetInputList()->FindObject("electronGsfSF");
    if (eGsfSF_ == nullptr ) 
        std::invalid_argument("Must pass electron GSF SF");
    mIdSF_ = (ScaleFactor *) GetInputList()->FindObject("muonTightIdSF");
    if (mIdSF_ == nullptr ) 
        std::invalid_argument("Must pass muon ID SF");
    mIsoSF_ = (ScaleFactor *) GetInputList()->FindObject("muonIsoSF");
    if (mIsoSF_ == nullptr ) 
        std::invalid_argument("Must pass muon Iso SF");

    prefireEff_ = (TEfficiency*) GetInputList()->FindObject("prefireEfficiencyMap");
    if (prefireEff_ == nullptr ) 
        std::invalid_argument("Must pass prefiring efficiency map");
}

void ZSelector::SetBranchesUWVV() {
    if (isMC_){
        fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
        fChain->SetBranchAddress("nTruePU", &nTruePU, &b_nTruePU);
    }
    fChain->SetBranchAddress("Mass", &ZMass, &b_ZMass);

    if (channel_ == ee) {
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
    }
    else if (channel_ == mm) {
        fChain->SetBranchAddress("m1IsWZTight", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("m2IsWZTight", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("m1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("m2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("m2Eta", &l2Eta, &b_l2Eta);
    }

    fChain->SetBranchAddress("type1_pfMETEt", &MET, &b_MET);
    fChain->SetBranchAddress("type1_pfMETPhi", &type1_pfMETPhi, &b_type1_pfMETPhi);
    //fChain->SetBranchAddress("nCBVIDTightElec", &nCBVIDTightElec, &b_nCBVIDTightElec);
    //fChain->SetBranchAddress("nCBVIDHLTSafeElec", &nCBVIDHLTSafeElec, &b_nCBVIDHLTSafeElec);
    //fChain->SetBranchAddress("nWZTightMuon", &nWZTightMuon, &b_nWZTightMuon);
    //fChain->SetBranchAddress("nWZMediumMuon", &nWZMediumMuon, &b_nWZMediumMuon);
}

void ZSelector::SetBranchesNanoAOD() {
  b.CleanUp();
  b.SetBranch("nElectron", nElectron);
  b.SetBranch("nMuon", nMuon);
  b.SetBranch("Electron_pt", Electron_pt);
  b.SetBranch("Electron_eta", Electron_eta);
  b.SetBranch("Electron_phi", Electron_phi);
  b.SetBranch("Muon_pt", Muon_pt);
  b.SetBranch("Muon_eta", Muon_eta);
  b.SetBranch("Muon_phi", Muon_phi);
  b.SetBranch("Electron_cutBased", Electron_cutBased);
  b.SetBranch("Muon_tightId", Muon_tightId);
  b.SetBranch("Muon_mediumId", Muon_mediumId);
  b.SetBranch("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
  b.SetBranch("MET_pt", MET);
  b.SetBranch("MET_phi", type1_pfMETPhi);
  b.SetBranch("Electron_charge", Electron_charge);
  b.SetBranch("Muon_charge", Muon_charge);
  b.SetBranch("Electron_mass", Electron_mass);
  b.SetBranch("Muon_mass", Muon_mass);
  //b.SetBranch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", Dimuon_Trigger);
  b.SetBranch("HLT_IsoMu24", SingleMuon_Trigger);
  //b.SetBranch("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", Dielectron_Trigger);
  b.SetBranch("HLT_Ele27_WPTight_Gsf", SingleElectron_Trigger);
  
  if (isMC_) {
    //b.SetBranch("e1GenPt", l1GenPt, );
    //b.SetBranch("e2GenPt", l2GenPt, );
    //b.SetBranch("e3GenPt", l3GenPt, );
    b.SetBranch("genWeight", genWeight);
    b.SetBranch("Pileup_nPU", numPU);
  }
}

void ZSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    weight = 1;
    b.SetEntry(entry);
    if (nElectron > N_KEEP_MU_E_ || nMuon > N_KEEP_MU_E_) {
        std::string message = "Found more electrons or muons than max read number.\n    Found ";
 	message += std::to_string(nElectron);
	message += " electrons.\n    Found ";
        message += std::to_string(nMuon);
	message += " Muons\n  --> Max read number was ";
        message += std::to_string(N_KEEP_MU_E_);
	message += "\nExiting because this can cause problems. Increase N_KEEP_MU_E_ to avoid this error.\n";
        throw std::domain_error(message);
    }

    ZMass = 0;
    l1Pt = 0;
    l2Pt = 0;
    l1Eta = 0;
    l2Eta = 0;
    l1Phi = 0;
    l2Phi = 0;
    l1Mass = 0;
    l2Mass = 0;

    // cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
    nCBVIDTightElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 4);
    nCBVIDVetoElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 1);
    nTightIdMuon = 0;
    nMediumIdMuon = 0;

    for (size_t i = 0; i < nMuon; i++) {
        //nMediumIdMuon += (Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15);
        nMediumIdMuon += (Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15*Muon_pt[i]);
        nTightIdMuon += (Muon_tightId[i] && Muon_pfRelIso04_all[i] < 0.15);
    }

    channel_ = channelMap_[channelName_];
    std::vector<size_t> goodIndices = {};

    if (nMediumIdMuon >= 2) {
        channel_ = mm;
        channelName_ = "mm";
        if (!(Muon_mediumId[0] && Muon_pfRelIso04_all[0] < 0.15
                    && Muon_mediumId[1] && Muon_pfRelIso04_all[1] < 0.15)) {
            for (size_t i = 0; i < nMuon; i++) {
                if (Muon_mediumId[i])
                    goodIndices.push_back(i);
            }
            if (goodIndices.size() < 2) {
                return;
            }
        }
        else 
            goodIndices = {0, 1};
        if (Muon_charge[goodIndices[0]] != Muon_charge[goodIndices[1]]) {
            l1Pt = Muon_pt[goodIndices[0]];
            l2Pt = Muon_pt[goodIndices[1]];
            l1Eta = Muon_eta[goodIndices[0]];
            l2Eta = Muon_eta[goodIndices[1]];
            l1Phi = Muon_phi[goodIndices[0]];
            l2Phi = Muon_phi[goodIndices[1]];
            l1Mass = Muon_mass[goodIndices[0]];
            l2Mass = Muon_mass[goodIndices[1]];
            l1IsTight = (Muon_mediumId[goodIndices[0]] && (Muon_pfRelIso04_all[goodIndices[0]] < 0.15));
            l2IsTight = (Muon_mediumId[goodIndices[1]] && (Muon_pfRelIso04_all[goodIndices[1]] < 0.15));
        }
    }
    else if (nCBVIDTightElec >= 2) {
        channel_ = ee;
        channelName_ = "ee";
        if (!(Electron_cutBased[0] == 4 && Electron_cutBased[1] == 4)) {
            for (size_t i = 0; i < nElectron; i++) {
                if (Electron_cutBased[i] == 4)
                    goodIndices.push_back(i);
            }
            if (goodIndices.size() < 2)
                return;
        }
        else 
            goodIndices = {0, 1};
        if (Electron_charge[goodIndices[0]] != Electron_charge[goodIndices[1]]) {
            l1Pt = Electron_pt[goodIndices[0]];
            l2Pt = Electron_pt[goodIndices[1]];
            l1Eta = Electron_eta[goodIndices[0]];
            l2Eta = Electron_eta[goodIndices[1]];
            l1Phi = Electron_phi[goodIndices[0]];
            l2Phi = Electron_phi[goodIndices[1]];
            l1Mass = Electron_mass[goodIndices[0]];
            l2Mass = Electron_mass[goodIndices[1]];
            l1IsTight = (Electron_cutBased[goodIndices[0]] == 4);
            l2IsTight = (Electron_cutBased[goodIndices[1]] == 4);
        }
    }
    else {
        channel_ = Unknown;
        channelName_ = "Unknown";
    }
    SetComposite();

    if (isMC_) {
        weight = genWeight;
        //b_l1GenPt->GetEntry(entry);
        //b_l2GenPt->GetEntry(entry);
        //b_l3GenPt->GetEntry(entry);
        ApplyScaleFactors();
    }
    else {
        //TODO: add MET filters
        //b_Flag_duplicateMuonsPass->GetEntry(entry);          
        //b_Flag_badMuonsPass->GetEntry(entry);          
    }



    //if (!singleLepton_)
    //    //passesTrigger = (Dimuon_Trigger || SingleMuon_Trigger ||
    //    //        Dielectron_Trigger || SingleElectron_Trigger);
    //    passesTrigger = (Dimuon_Trigger || SingleMuon_Trigger);
    //else
    //    passesTrigger = (!Dimuon_Trigger && SingleMuon_Trigger);
    //    //passesTrigger = ((!Dimuon_Trigger && SingleMuon_Trigger) ||
    //    //        (!Dielectron_Trigger && SingleElectron_Trigger));
    passesTrigger = SingleMuon_Trigger || SingleElectron_Trigger;

    passesLeptonVeto = (std::min(nMediumIdMuon, nLooseIsoMuon) + nCBVIDVetoElec) == 2;
}

void ZSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation){ 
    weight = 1;
    b_l1Pt->GetEntry(entry);
    b_l2Pt->GetEntry(entry);
    b_l1Eta->GetEntry(entry);
    b_l2Eta->GetEntry(entry);

    if (isMC_) {
        b_genWeight->GetEntry(entry);
    }
    b_ZMass->GetEntry(entry);
    b_l1IsTight->GetEntry(entry);
    b_l2IsTight->GetEntry(entry);
    b_MET->GetEntry(entry);
    //b_nCBVIDTightElec->GetEntry(entry);
    //
    passesTrigger = true;
}

void ZSelector::ApplyScaleFactors() {
    if (channel_ == ee) {
        if (eIdSF_ != nullptr) {
	    weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        if (eGsfSF_ != nullptr) {
	    weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
    }
    else if (channel_ == mm) {
        if (mIdSF_ != nullptr) {
	    weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
	}
        if (mIsoSF_ != nullptr) {
	    weight *= mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
	}
    }
    if (pileupSF_ != nullptr) {
        weight *= pileupSF_->Evaluate1D(numPU);
    }
}

void ZSelector::SetComposite() {
    if (l1Pt == 0. || l2Pt == 0.) {
        return;
    }
    TLorentzVector lepton1;
    lepton1.SetPtEtaPhiM(l1Pt, l1Eta, l1Phi, l1Mass);
    TLorentzVector lepton2;
    lepton2.SetPtEtaPhiM(l2Pt, l2Eta, l2Phi, l2Mass);
    auto system = lepton1+lepton2;
    ZMass = system.M();
    ZPt = system.Pt();
    ZEta = system.PseudoRapidity();
    Zy = system.Rapidity();
}

// Meant to be a wrapper for the tight ID just in case it changes
// To be a function of multiple variables
bool ZSelector::zlep1IsTight() {
    return l1IsTight; 
}

bool ZSelector::zlep2IsTight() {
    return l2IsTight; 
}

bool ZSelector::tightZLeptons() {
    return zlep1IsTight() && zlep2IsTight(); 
}

void ZSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    int step = 0;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    if (channel_ != mm && channel_ != ee) 
        return;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    if (!passesTrigger)
        return;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    if (channel_ == ee && (std::abs(l1Eta) > 2.4 || std::abs(l2Eta) > 2.4 ))
        return;
    else if (channel_ == mm && (std::abs(l1Eta) > 2.5 || std::abs(l2Eta) > 2.5 ))
        return;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    if (l1Pt < 25 || l2Pt < 25)
        return;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    if (ZMass > 106.1876 || ZMass < 76.1876)
        return;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    //if (MET > 25)
    //    return;
    //SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    if (!tightZLeptons())
        return;
    SafeHistFill(histMap1D_, getHistName("CutFlow", variation.second), step++, weight);

    SafeHistFill(histMap1D_, getHistName("ZMass", variation.second), ZMass, weight);
    SafeHistFill(histMap1D_, getHistName("ptl1", variation.second), l1Pt, weight);
    SafeHistFill(histMap1D_, getHistName("ptl2", variation.second), l2Pt, weight);
    SafeHistFill(histMap1D_, getHistName("etal1", variation.second), l2Eta, weight);
    SafeHistFill(histMap1D_, getHistName("etal2", variation.second), l2Eta, weight);
    SafeHistFill(histMap1D_, getHistName("ZEta", variation.second), ZEta, weight);
    SafeHistFill(histMap1D_, getHistName("yZ", variation.second), Zy, weight);
    SafeHistFill(histMap1D_, getHistName("ZPt", variation.second), ZPt, weight);
}

void ZSelector::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();

    InitializeHistogramsFromConfig();
}
