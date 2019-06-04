
#include "Analysis/VVAnalysis/interface/SameSignLept.h"
#include "TLorentzVector.h"
#include <TStyle.h>
#include <regex>
#include "TParameter.h"

void SameSignLept::Init(TTree *tree)
{
    SelectorBase::Init(tree);
    if (addSumweights_) {
        TFile* file = fChain->GetTree()->GetDirectory()->GetFile(); 
        TTree* metaInfo = dynamic_cast<TTree*>(file->Get("metaInfo/metaInfo"));
        if (metaInfo == nullptr)
            std::cerr << "WARNING: Failed to add sumWeights histogram" << std::endl;
        else {
            metaInfo->Draw("1>>sumweights", "summedWeights");
        }
    }
}

void SameSignLept::SetBranchesUWVV() {
    throw std::domain_error("UWVV ntuples not defined for Z selector!");
}

void SameSignLept::SetBranchesNanoAOD() {
    fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
    fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
    fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
    fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
    fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
    fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
    fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
    fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
    fChain->SetBranchAddress("Electron_cutBased", &Electron_cutBased, &b_Electron_cutBased);
    fChain->SetBranchAddress("Muon_tightId", &Muon_tightId, &b_Muon_tightId);
    fChain->SetBranchAddress("Muon_mediumId", &Muon_mediumId, &b_Muon_mediumId);
    fChain->SetBranchAddress("Muon_pfIsoId", &Muon_pfIsoId, &b_Muon_pfIsoId);
    fChain->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
    fChain->SetBranchAddress("MET_pt", &MET, &b_MET);
    fChain->SetBranchAddress("MET_phi", &type1_pfMETPhi, &b_type1_pfMETPhi);
    fChain->SetBranchAddress("Electron_charge", &Electron_charge, &b_Electron_charge);
    fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
    fChain->SetBranchAddress("Electron_mass", &Electron_mass, &b_Electron_mass);
    fChain->SetBranchAddress("Muon_mass", &Muon_mass, &b_Muon_mass);
    if (isMC_) {
        fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
        fChain->SetBranchAddress("Pileup_nPU", &numPU, &b_numPU);
    }
}

void SameSignLept::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    weight = 1;
    b_nElectron->GetEntry(entry);
    b_nMuon->GetEntry(entry);
    b_Electron_pt->GetEntry(entry);
    b_Electron_eta->GetEntry(entry);
    b_Electron_phi->GetEntry(entry);
    b_Muon_pt->GetEntry(entry);
    b_Muon_eta->GetEntry(entry);
    b_Muon_phi->GetEntry(entry);
    b_Electron_cutBased->GetEntry(entry);
    b_Muon_tightId->GetEntry(entry);
    b_Muon_mediumId->GetEntry(entry);
    b_Muon_pfRelIso04_all->GetEntry(entry);
    b_Electron_charge->GetEntry(entry);
    b_Muon_charge->GetEntry(entry);
    b_Electron_mass->GetEntry(entry);
    b_Muon_mass->GetEntry(entry);
    b_MET->GetEntry(entry);

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
    // for this, we will need tight electron, medium muon...will this suffice or 
    // do I need something else...?
    nCBVIDTightElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 4);
    nCBVIDVetoElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 1);
    nTightIdMuon = std::count(Muon_tightId, Muon_tightId+nMuon, true);
    nMediumIdMuon = std::count(Muon_mediumId, Muon_mediumId+nMuon, true);
    //nTightIsoMuon = std::count(Muon_pfIsoId, Muon_pfIsoId+nMuon, 4);
    //nLooseIsoMuon = std::count(Muon_pfIsoId, Muon_pfIsoId+nMuon, 1);

    // originally tight muon...let's try medium
    channel_ = channelMap_[channelName_];
    std::vector<size_t> goodIndices = {};
    if (nMediumIdMuon >= 2) {
        channel_ = mm;
        if (!(Muon_mediumId[0] && Muon_mediumId[1])) {
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
	// want same sign so changing from != to ==........
	// additionally changing tight to medium....
	// need to make sure we are calling mediumId elsewhere for muons
        if (Muon_charge[goodIndices[0]] == Muon_charge[goodIndices[1]]) {
            l1Pt = Muon_pt[goodIndices[0]];
            l2Pt = Muon_pt[goodIndices[1]];
            l1Eta = Muon_eta[goodIndices[0]];
            l2Eta = Muon_eta[goodIndices[1]];
            l1Phi = Muon_phi[goodIndices[0]];
            l2Phi = Muon_phi[goodIndices[1]];
            l1Mass = Muon_mass[goodIndices[0]];
            l2Mass = Muon_mass[goodIndices[1]];
            l1IsMed = (Muon_mediumId[goodIndices[0]] && (Muon_pfRelIso04_all[goodIndices[0]] < 0.15));
            l2IsMed = (Muon_mediumId[goodIndices[1]] && (Muon_pfRelIso04_all[goodIndices[1]] < 0.15));
        }
    }
    else if (nCBVIDTightElec >= 2) {
        channel_ = ee;
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
	// wanting same charge. changing to ==
        if (Electron_charge[goodIndices[0]] == Electron_charge[goodIndices[1]]) {
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
    SetMass();

    if (isMC_) {
        b_genWeight->GetEntry(entry);
        //TODO: add scale factors
        //b_numPU->GetEntry(entry);
        //b_l1GenPt->GetEntry(entry);
        //b_l2GenPt->GetEntry(entry);
        //b_l3GenPt->GetEntry(entry);
        //ApplyScaleFactors();
    }
    else {
        //TODO: add MET filters
        //b_Flag_duplicateMuonsPass->GetEntry(entry);          
        //b_Flag_badMuonsPass->GetEntry(entry);          
    }
    // originally has nMediumIdMuon but using medium id...should we change to tight here? 
    passesLeptonVeto = (std::min(nTightIdMuon, nLooseIsoMuon) + nCBVIDVetoElec) == 2;
}

void SameSignLept::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation){ 
    throw std::domain_error("UWVV ntuples not defined for Z selector!");
}

void SameSignLept::ApplyScaleFactors() {
    weight = genWeight;
    if (channel_ == ee) {
        weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
    }
    else if (channel_ == mm) {
        weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
    }
    weight *= pileupSF_->Evaluate1D(numPU);
}

void SameSignLept::SetMass() {
    if (l1Pt == 0. || l2Pt == 0.) {
        return;
    }
    TLorentzVector lepton1;
    lepton1.SetPtEtaPhiM(l1Pt, l1Eta, l1Phi, l1Mass);
    TLorentzVector lepton2;
    lepton2.SetPtEtaPhiM(l2Pt, l2Eta, l2Phi, l2Mass);
    ZMass = (lepton1+lepton2).M();
}

// Meant to be a wrapper for the tight ID just in case it changes
// To be a function of multiple variables

// Modifying a lot since Muons are medium and electrons are tight
bool SameSignLept::zlep1IsTight() {
  if (channel_==ee) {
    returning_value = l1IsTight;
  } else {
    returning_value = false;
  }
  return returning_value;
}
bool SameSignLept::zlep1IsMed() {
  if (channel_==mm) {
    returning_value =  l1IsMed;
  } else {
    returning_value = false;
  }
  return returning_value;
}
bool SameSignLept::zlep2IsTight() {
  if (channel_==ee) {
    returning_value =  l2IsTight;
  } else {
    returning_value = false;
  }
  return returning_value;
}
bool SameSignLept::zlep2IsMed() {
  if (channel_==mm) {
    returning_value= l2IsMed;
  } else {
    returning_value =  false;
  }
  return returning_value;
}
bool SameSignLept::tightZLeptons() {
    return zlep1IsTight() && zlep2IsTight(); 
}
bool SameSignLept::tight_med_ZLeptons() {
  return (zlep1IsTight() && zlep2IsMed()) or (zlep1IsMed() && zlep2IsTight());
}
bool SameSignLept::medZLeptons() {
  return zlep1IsMed() && zlep2IsMed();
}

void SameSignLept::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    cutflow_ee_->Fill(0.,weight);
    cutflow_mm_->Fill(0.,weight);
    if (channel_ == ee)
        cutflow_ee_->Fill(1.,weight);
    else if (channel_ == mm)
        cutflow_mm_->Fill(1.,weight);

    //if (!passesLeptonVeto)
    //    return;

    if (channel_ == ee && (std::abs(l1Eta) > 2.4 || std::abs(l2Eta) > 2.4 ))
        return;
    else if (channel_ == mm && (std::abs(l1Eta) > 2.5 || std::abs(l2Eta) > 2.5 ))
        return;
    if (channel_ == ee)
        cutflow_ee_->Fill(2,weight);
    else if (channel_ == mm)
        cutflow_mm_->Fill(2,weight);

    if (l1Pt < 25 || l2Pt < 20)
        return;
    if (channel_ == ee)
        cutflow_ee_->Fill(3,weight);
    else if (channel_ == mm)
        cutflow_mm_->Fill(3,weight);

    if (ZMass > 101.1876 || ZMass < 81.1876)
        return;

    if (channel_ == ee)
        cutflow_ee_->Fill(4,weight);
    else if (channel_ == mm)
        cutflow_mm_->Fill(4,weight);

    if (MET > 25)
        return;

    if (channel_ == ee)
        cutflow_ee_->Fill(5,weight);
    else if (channel_ == mm)
        cutflow_mm_->Fill(5,weight);

    if (!tightZLeptons())
        return;

    if (channel_ == ee)
        cutflow_ee_->Fill(6,weight);
    else if (channel_ == mm)
        cutflow_mm_->Fill(6,weight);

    if (channel_ == ee) {
        ZMass_ee_->Fill(ZMass, weight);
        ptl1_ee_->Fill(l1Pt, weight);
        ptl2_ee_->Fill(l2Pt, weight);
	l1eta_ee_->Fill(l1Eta, weight);
	l2eta_ee_->Fill(l2Eta, weight);
	l1phi_ee_->Fill(l1Phi, weight);
	l2phi_ee_->Fill(l2Phi, weight);
 }
    else if (channel_ == mm) {
        ZMass_mm_->Fill(ZMass, weight);
        ptl1_mm_->Fill(l1Pt, weight);
        ptl2_mm_->Fill(l2Pt, weight);
	l1eta_mm_->Fill(l1Eta, weight);
	l2eta_mm_->Fill(l2Eta, weight);
	l1phi_mm_->Fill(l1Phi, weight);
	l2phi_mm_->Fill(l2Phi, weight);
 }
    else
        throw std::domain_error("Invalid channel!");
}

void SameSignLept::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();

    AddObject<TH1D>(cutflow_ee_, "cutflow_ee", "Tight leptons; Cut flow", 7, 0, 7);
    AddObject<TH1D>(cutflow_mm_, "cutflow_mm", "Tight leptons; Cut flow", 7, 0, 7);
    AddObject<TH1D>(ZMass_ee_, "ZMass_ee", "Tight leptons; m_{ee} [GeV]", 80, 102, 22);
    AddObject<TH1D>(ZMass_mm_, "ZMass_mm", "Tight leptons; m_{#mu#mu} [GeV]", 80, 102, 22);
    AddObject<TH1D>(ptl1_ee_, "ptl1_ee", "Tight leptons; p_{T}(e_{1}) [GeV]", 100, 0, 200);
    AddObject<TH1D>(ptl1_mm_, "ptl1_mm", "Tight leptons; p_{T}(#mu_{1}) [GeV]", 100, 0, 100);
    AddObject<TH1D>(ptl2_ee_, "ptl2_ee", "Tight leptons; p_{T}(e_{2}) [GeV]", 100, 0, 200);
    AddObject<TH1D>(ptl2_mm_, "ptl2_mm", "Tight leptons; p_{T}(#mu_{2}) [GeV]", 100, 0, 100);

    AddObject<TH1D>(l1eta_mm_, "l1eta_mm", "Tight leptons; eta(#mu_{1})", 100, -5, 5);
    AddObject<TH1D>(l2eta_mm_, "l2eta_mm", "Tight leptons; eta(#mu_{2})", 100, -5, 5);
    AddObject<TH1D>(l1eta_ee_, "l1eta_ee", "Tight leptons; eta(e_{1})", 100, -5, 5);
    AddObject<TH1D>(l2eta_ee_, "l2eta_ee", "Tight leptons; eta(e_{2})", 100, -5, 5);

    AddObject<TH1D>(l1phi_mm_, "l1phi_mm", "Tight leptons; #phi(#mu_{1})", 100, -20, 20);
    AddObject<TH1D>(l2phi_mm_, "l2phi_mm", "Tight leptons; phi(#mu_{2})", 100, -20, 20);
    AddObject<TH1D>(l1phi_ee_, "l1phi_ee", "Tight leptons; phi(e_{1})", 100, -20, 20);
    AddObject<TH1D>(l2phi_ee_, "l2phi_ee", "Tight leptons; phi(e_{2})", 100, -20, 20);

} 
