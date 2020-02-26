#include "Analysis/VVAnalysis/interface/FakeRateSelector.h"
#include <TStyle.h>
#include "TLorentzVector.h"


#define Fill1D(NAME, VALUE_) HistFullFill(histMap1D_, NAME, variation.first, VALUE_, weight);
#define Fill2D(NAME, VALUE1_, VALUE2_) HistFullFill(histMap2D_, NAME, variation.first, VALUE1_, VALUE2_, weight);


typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>LorentzVector;
namespace VectorUtil=ROOT::Math::VectorUtil;

void FakeRateSelector::Init(TTree *tree) {
    allChannels_ = {{e, "e"}, {m, "m"}};
    
    SelectorBase::Init(tree);
    
    TList* slaveClassList = (TList*)GetInputList()->Clone();
    slaveClassList->Add(new TNamed("isSlaveClass", "isSlaveClass"));
    Analyzer.SetInputList(slaveClassList);
    Analyzer.Init(tree);
    Analyzer.SetBranches();
    
}

void FakeRateSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) {
    Analyzer.LoadBranchesNanoAOD(entry,variation);
    channelName_ = "all";
    if(Analyzer.looseLeptons.size() == 1) {
        if(Analyzer.looseLeptons[0].Id() == PID_MUON) channelName_ = "m";
	else channelName_ = "e";
    } 
    channel_ = channelMap_[channelName_];
}

void FakeRateSelector::SetBranchesNanoAOD() {}

void FakeRateSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) {
    // need 3 loose leptons
    size_t N_LEPS_FAKE = 1;
    
    if (Analyzer.looseLeptons.size() != N_LEPS_FAKE)    return;
    int tagLepIdx = -1;
    for(size_t i=0; i < Analyzer.looseLeptons.size(); ++i) {
	auto lep = Analyzer.looseLeptons[i];
	if(Analyzer.passFakeableCuts(lep)) {
	    if(tagLepIdx >= 0) return;
	    tagLepIdx = i;
	}
    }
    LorentzVector met(Analyzer.MET, 0, Analyzer.type1_pfMETPhi, Analyzer.MET);
    LorentzVector tagLep = Analyzer.looseLeptons[tagLepIdx].v;
    double mt = std::sqrt(2.*met.Et()*tagLep.pt()*(1.-cos(VectorUtil::DeltaPhi(tagLep, met))));
        
    if (mt > 20) return;
    if (met.Pt() > 20) return;  // need to fix with skim
    bool hasJet = false;
    for(auto jet: Analyzer.goodJets) {
	if(reco::deltaR(jet.v, tagLep) > 1.) {
	    hasJet = true;
	    break;
	}
    }
    if(!hasJet) return;
    
    float tagPt = tagLep.Pt();
    float tagEta = std::abs(tagLep.Eta());

    // int closeIdx = Analyzer.getCloseJetIndex(tagLep);
    // LorentzVector closeJet  = Analyzer.get4Vector(PID_JET, closeIdx);
    // if(Analyzer.LepRelPt(tagLep, closeJet) > 7.2)
    // 	tagPt *= (1 + std::max(0, Im-I1));
    // else
    // 	tagPt = std::max(tagPt, closeJet.Pt()*I2);
    
    float loose_weight = 1.;
    // float loose_weight = weight;
    // if (channel_ == eee || channel_ == emm) {
    // 	loose_weight /= eIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
    // }p
    // else if (channel_ == eem || channel_ == mmm) {
    // 	loose_weight /= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
    // }
    
    passingLoose2D_map[channel_]->Fill(tagPt, tagEta, loose_weight);
    passingLoose1DPt_map[channel_]->Fill(tagPt, loose_weight);
    passingLoose1DEta_map[channel_]->Fill(tagEta, loose_weight);
    if (Analyzer.goodLeptons.size() == N_LEPS_FAKE) {
    	passingTight2D_map[channel_]->Fill(tagPt, tagEta, weight);
    	passingTight1DPt_map[channel_]->Fill(tagPt, weight);
    	passingTight1DEta_map[channel_]->Fill(tagEta, weight);
    }
}

void FakeRateSelector::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();

    const int nPtBins = 5;
    const int nEtaBins = 3;
    double pt_bins[nPtBins+1] = {10, 15, 25, 35, 50, 70};
    double eta_El_bins[nEtaBins+1] = {0, 0.8, 1.479, 2.5};
    double eta_Mu_bins[nEtaBins+1] = {0, 1.2, 2.1, 2.4};
    std::unordered_map<int, double*> eta_bins= {{e, eta_El_bins}, {m, eta_Mu_bins}};
    for(auto pair : allChannels_) {
	passingTight2D_map[pair.first] = nullptr;
	passingTight1DPt_map[pair.first] = nullptr;
	passingTight1DEta_map[pair.first] = nullptr;
	passingLoose2D_map[pair.first] = nullptr;
	passingLoose1DPt_map[pair.first] = nullptr;
	passingLoose1DEta_map[pair.first] = nullptr;
	
	AddObject<TH2D>(passingTight2D_map[pair.first], ("passingTight2D_"+pair.second).c_str(), "#eta; p_{T} [GeV]", nPtBins, pt_bins, nEtaBins, eta_bins[pair.first]);
	AddObject<TH1D>(passingTight1DPt_map[pair.first], ("passingTight1DPt_"+pair.second).c_str(), "Tight leptons; p_{T} [GeV]", nPtBins, pt_bins);
	AddObject<TH1D>(passingTight1DEta_map[pair.first], ("passingTight1DEta_"+pair.second).c_str(), "Tight leptons; #eta", nEtaBins, eta_bins[pair.first]);
    
	AddObject<TH2D>(passingLoose2D_map[pair.first], ("passingLoose2D_"+pair.second).c_str(), "#eta; p_{T} [GeV]", nPtBins, pt_bins, nEtaBins, eta_bins[pair.first]);
	AddObject<TH1D>(passingLoose1DPt_map[pair.first], ("passingLoose1DPt_"+pair.second).c_str(), "Loose leptons; p_{T} [GeV]", nPtBins, pt_bins);
	AddObject<TH1D>(passingLoose1DEta_map[pair.first], ("passingLoose1DEta_"+pair.second).c_str(), "Loose leptons; #eta", nEtaBins, eta_bins[pair.first]);
    }
}
