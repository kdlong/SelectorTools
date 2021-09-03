#include "Analysis/SelectorTools/interface/WGenSelector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <TStyle.h>
#include <regex>
#include <cmath>
#include <TRandom3.h>
#include "TParameter.h"

void WGenSelector::Init(TTree *tree)
{
    // Don't waste memory on empty e hists
    TParameter<bool>* muOnlyParam = (TParameter<bool>*) GetInputList()->FindObject("muOnly");
    bool muOnly = muOnlyParam != nullptr && muOnlyParam->GetVal();
    allChannels_ = {{mp, "mp"}, {mn, "mn"}}; 
    if (!muOnly) {
        allChannels_.push_back(std::make_pair<Channel, std::string>(ep, "ep"));
        allChannels_.push_back(std::make_pair<Channel, std::string>(en, "en"));
    }

    histMap1D_[{"CutFlow", Unknown, Central}] = {};
    hists1D_ = {"CutFlow", "mWmet", "yWmet", "ptWmet", "mW", "yW", "ptW", "mTtrue", "mTmet",
        "ptl", "etal", "phil", "ptnu", "etanu", "phinu", "MET", "MET_phi",
        "ptj1", "ptj2", "etaj1", "etaj2", "nJets","Ratio_Wmass", 
        "dRlgamma_maxptassoc", "dRlgamma_minassoc", "ptg_closeassoc", "ptg_maxassoc", "nGammaAssoc", 
        "ptgmax_assoc", "ptgmax_assoc",
        "ptl_smear",
    };
    hists2D_ = {"etal_ptl_2D", "etal_ptl_smear_2D"};
    hists3D_ = {"mass_y_pT_3D", };
    systHists3D_ = hists3D_;

    TParameter<bool>* massVar = (TParameter<bool>*) GetInputList()->FindObject("massVar");
    doMassVar_ = massVar != nullptr && massVar->GetVal();
    TParameter<bool>* muonVar = (TParameter<bool>*) GetInputList()->FindObject("muonVar");
    doMuonVar_ = muonVar != nullptr && muonVar->GetVal();
    // Have to read this here as well, otherwise you don't know until after calling NanoGenSelectorBase::Init
    TParameter<bool>* barePart = (TParameter<bool>*) GetInputList()->FindObject("bare");
    doBareLeptons_ = barePart != nullptr && barePart->GetVal();

    if (doMassVar_) {
        systematics_[mWShift50MeVUp] = "mWBWShift50MeVUp";
        systematics_[mWShift50MeVDown] = "mWBWShift50MeVDown";
        systematics_[mWShift100MeVUp] = "mWBWShift100MeVUp";
        systematics_[mWShift100MeVDown] = "mWBWShift100MeVDown";
        if (doBareLeptons_) {
            systematics_[BareLeptons_mWShift100MeVUp] = "bare_mWBWShift100MeVUp";
            systematics_[BareLeptons_mWShift100MeVDown] = "bare_mWBWShift100MeVDown";
        }
    }

    if (doMuonVar_) {
        systematics_[muonScaleUp] = "CMS_scale_mUp";
        systematics_[muonScaleDown] = "CMS_scale_mDown";
        systematics_[muonResolutionUp] = "CMS_res_mUp";
        systematics_[muonResolutionDown] = "CMS_res_mDown";
        if (doBareLeptons_) {
            systematics_[BareLeptons_muonScaleUp] = "bare_CMS_scale_mUp";
            systematics_[BareLeptons_muonScaleDown] = "bare_CMS_scale_mDown";
        }
    }

    systHists_ = {"ptl", "ptl_smear", "yW", "ptW", "mW", "mTrue", "mTmet"};
    systHists2D_ = hists2D_;

    weighthists1D_ = systHists_;
    weighthists2D_ = hists2D_;

    refWeight = 1;
    nLeptons_ = 1;
    nNeutrinos_ = 1;
    doPhotons_ = true;
    // Chose by MC sample
    if (name_.find("nnlops") != std::string::npos) {
        MV_GEN_ = 80398.0;
        GAMMAV_GEN_ = 2088.720;
    }
    else if (name_.find("minnlo") != std::string::npos) {
        MV_GEN_ = 80351.812293789408;
        GAMMAV_GEN_ = 2090.4310808144846;
    }
    else {
        MV_GEN_ = 80419.;
        GAMMAV_GEN_ = 2050;
    }

    NanoGenSelectorBase::Init(tree);

    if (name_.find("N3LLCorr") != std::string::npos) {
        n3llcorr_ = true;
        SetScaleFactors();
        if (scetlibCorrs_.at(0) == nullptr)
            throw std::invalid_argument("Must pass a scalefactor for N3LLCorr file!");
    }
}

void WGenSelector::LoadBranchesNanoAOD(Long64_t entry, SystPair variation) { 
    NanoGenSelectorBase::LoadBranchesNanoAOD(entry, variation);

    // Make sure the central is used for all variations
    // A bit lazy here since I don't actually use the smeared pt for composite quantities
    ptl_smear_fill = ptl_smear;
    if (leptons.size() >= nLeptons_) {
        auto& l = leptons.at(0);
        if (variation.first == Central || variation.first == muonResolutionUp || variation.first == muonResolutionDown) {
            TRandom3 gauss;
            float w = 0.01;
            if (variation.first == muonResolutionUp)
                w += 0.001;
            else if (variation.first == muonResolutionDown)
                w -= 0.001;
            ptl_smear = l.pt()*gauss.Gaus(1, w);
            if (variation.first == Central)
                ptl_smear_fill = ptl_smear;
        }
        else if (variation.first == muonScaleUp || variation.first == BareLeptons_muonScaleUp) {
            leptons.at(0).setP4(makeGenParticle(l.pdgId(), l.status(), l.pt()*1.001, l.eta(), l.phi(), l.mass()).polarP4());
            SetComposite();
            ptl_smear_fill *= 1.001;
        }
        else if (variation.first == muonScaleDown || variation.first == BareLeptons_muonScaleDown) {
            leptons.at(0).setP4(makeGenParticle(l.pdgId(), l.status(), l.pt()*1./1.001, l.eta(), l.phi(), l.mass()).polarP4());
            SetComposite();
            ptl_smear_fill *= 1./1.001;
        }
    }
     
    // This means the observable is corrected differently depending on 
    // the type of leptons. This isn't ideal but oh well
    if (variation.first == Central) {
        cenWeight = weight;
        ptVcorr = wCand.pt();
        mVcorr = wCand.mass();
        ratio_mass = wCand.mass();
    }
    else if (variation.first == PreFSRLeptons) { 
            ptVcorr = wCand.pt();
            mVcorr = wCand.mass();
            yVcorr = wCand.Rapidity();
    }
    else if (variation.first == LHEParticles && systematics_.find(PreFSRLeptons) == std::end(systematics_)) {
        // define at LHE level if it exists
        ptVcorr = wCand.pt();
        mVcorr = wCand.mass();
        yVcorr = wCand.Rapidity();
    }
    else if (variation.first == mWShift10MeVUp)
        weight = cenWeight*breitWignerWeight(10.);
    else if (variation.first == mWShift10MeVDown)
        weight = cenWeight*breitWignerWeight(-10.);
    else if (variation.first == mWShift20MeVUp)
        weight = cenWeight*breitWignerWeight(20.);
    else if (variation.first == mWShift20MeVDown)
        weight = cenWeight*breitWignerWeight(-20.);
    else if (variation.first == mWShift25MeVUp)
        weight = cenWeight*breitWignerWeight(25.);
    else if (variation.first == mWShift25MeVDown)
        weight = cenWeight*breitWignerWeight(-25.);
    else if (variation.first == mWShift50MeVUp) {
        weight = cenWeight*breitWignerWeight(50.);
    }
    else if (variation.first == mWShift50MeVDown)
        weight = cenWeight*breitWignerWeight(-50.);
    else if (variation.first == mWShift100MeVUp || variation.first == BareLeptons_mWShift100MeVUp)
        weight = cenWeight*breitWignerWeight(100.);
    else if (variation.first == mWShift100MeVDown || variation.first == BareLeptons_mWShift100MeVDown)
        weight = cenWeight*breitWignerWeight(-100.);


    if (leptons.size() > 0 && std::abs(leptons.at(0).pdgId()) == 11) {
        if (leptons.at(0).pdgId() > 0) {
            channel_ = en;
            channelName_ = "en";
        }
        else {
            channel_ = ep;
            channelName_ = "ep";
        }
    }
    else if (leptons.size() > 0 && std::abs(leptons.at(0).pdgId()) == 13) {
        if (leptons.at(0).pdgId() > 0) {
            channel_ = mn;
            channelName_ = "mn";
        }
        else {
            channel_ = mp;
            channelName_ = "mp";
        }
    }
    else {
        channel_ = Unknown;
        channelName_ = "Unknown";
        return;
    }

    if (n3llcorr_) {
        //auto* sf = channel_ == mp ? n3llWpSF_ : n3llWmSF_;
        weight *= scetlibCorrs_.at(0)->Evaluate3D(mVcorr, yVcorr, ptVcorr);
    }
}

void WGenSelector::SetComposite() {
    if (leptons.size() == 0) {
        wCandMet = LorentzVector();
        wCand = LorentzVector();
        return;
    }
    else if (neutrinos.size() == 0) {
        wCand = LorentzVector();
        return;
    }
    auto lepP4 = leptons.at(0).polarP4();
    auto compareByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() < b.pt(); };
    auto mt = [] (LorentzVector& l, LorentzVector& v) {
        return std::sqrt(2*l.pt()*v.pt()*(1 - cos(l.phi() - v.phi())));
    };

    auto nup = std::max_element(neutrinos.begin(), neutrinos.end(), compareByPt);
    if (neutrinos.size()) {
        nu = neutrinos.size() > 0 ? nup->polarP4() : LorentzVector();
        wCandMet = lepP4 + genMet;
        wCand = neutrinos.size() > 0 ? lepP4 + nu : LorentzVector();
        mTtrue = mt(lepP4, nu);
    }
    mTmet = mt(lepP4, genMet);
}

void WGenSelector::FillHistograms(Long64_t entry, SystPair variation) { 
    std::string lepType = "";
    FillHistogramsByName(entry, lepType, variation);
}

void WGenSelector::FillHistogramsByName(Long64_t entry, std::string& toAppend, SystPair variation) { 
    int step = 0;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (channel_ != mn && channel_ != en && channel_ != mp && channel_ != ep) 
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (leptons.size() < nLeptons_)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (neutrinos.size() < nLeptons_)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    auto& lep = leptons.at(0);
    if (doFiducial_ && std::abs(lep.eta()) > 2.5)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);


    // Fill smear histograms before

    if (doFiducial_ && lep.pt() < 25)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (variation.first == Central)
        mcWeights_->Fill(weight/std::abs(refWeight));

    if (std::find(theoryVarSysts_.begin(), theoryVarSysts_.end(), variation.first) != theoryVarSysts_.end()) {
        size_t nScaleWeights = nLHEScaleWeight+nLHEScaleWeightAltSet1;
        size_t minimalWeights = nScaleWeights+nMEParamWeight;

        // In the case of CT18, don't store all the CT18Z sets which are mixed in
        size_t nWeights = minimalWeights+pdfMaxStore_;
        if (n3llcorr_)
            nWeights += nScetlibWeights_;
        size_t pdfOffset = nScaleWeights;
        size_t pdfIdx = 0;
        for (size_t i = 0; i < nWeights; i++) {
            float thweight = 1;
            if (i < nLHEScaleWeight)
                thweight = LHEScaleWeight[i];
            else if (i < nScaleWeights)
                thweight = LHEScaleWeightAltSet1[i-nLHEScaleWeight];
            else if (i < nScaleWeights+pdfMaxStore_) {
                while (!pdfWeights_.at(pdfIdx))
                    pdfIdx++;
                thweight = LHEPdfWeights[pdfIdx][i-pdfOffset];
                if (i == pdfOffset+nLHEPdfWeights.at(pdfIdx)-1) {
                    pdfOffset += nLHEPdfWeights.at(pdfIdx++);
                }
            }
            else if (i < minimalWeights+pdfMaxStore_) {
                size_t offset = nScaleWeights+pdfMaxStore_;
                thweight = MEParamWeight[i-offset];
            }
            else {
                int idx = i-minimalWeights-pdfMaxStore_;
                auto* sf = scetlibCorrs_.at(idx);
                float refW = scetlibCorrs_.at(0)->Evaluate3D(mVcorr, yVcorr, ptVcorr);
                thweight = sf->Evaluate3D(mVcorr, yVcorr, ptVcorr)/refW;
            }

            thweight /= rescaleWeight;

            if (((variation.first == ptV0to3 || variation.first == ptV0to3_lhe) && ptVcorr > 3.) ||
                    ((variation.first == ptV3to5 || variation.first == ptV3to5_lhe) && (ptVcorr < 3. || ptVcorr > 5.))  ||
                    ((variation.first == ptV5to7 || variation.first == ptV5to7_lhe) && (ptVcorr < 5. || ptVcorr > 7.)) ||
                    ((variation.first == ptV7to9 || variation.first == ptV7to9_lhe) && (ptVcorr < 7. || ptVcorr > 9.)) ||
                    ((variation.first == ptV9to12 || variation.first == ptV9to12_lhe) && (ptVcorr < 9. || ptVcorr > 12.)) ||
                    ((variation.first == ptV12to15 || variation.first == ptV12to15_lhe) && (ptVcorr < 12. || ptVcorr > 15.)) ||
                    ((variation.first == ptV15to20 || variation.first == ptV15to20_lhe) && (ptVcorr < 15. || ptVcorr > 20.)) ||
                    ((variation.first == ptV20to27 || variation.first == ptV20to27_lhe) && (ptVcorr < 20. || ptVcorr > 27.)) ||
                    ((variation.first == ptV27to40 || variation.first == ptV27to40_lhe) && (ptVcorr < 27. || ptVcorr > 40.)) ||
                    ((variation.first == ptV40toInf || variation.first == ptV40toInf_lhe) && ptVcorr < 40. )) {
                thweight = 1;
            }

            thweight = (!thweightSuppress_ || std::abs(thweight) < thweightSuppress_) ? thweight : (thweight > 0 ? thweightSuppress_ : -1*thweightSuppress_);
            thweight *= weight;
            SafeHistFill(weighthistMap1D_, concatenateNames("mW", toAppend), channel_, variation.first, wCand.mass(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("yW", toAppend), channel_, variation.first, wCand.Rapidity(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptW", toAppend), channel_, variation.first, wCand.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("mWmet", toAppend), channel_, variation.first, wCandMet.mass(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("yWmet", toAppend), channel_, variation.first, wCandMet.Rapidity(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptWmet", toAppend), channel_, variation.first, wCandMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("MET", toAppend), channel_, variation.first, genMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("MET_phi", toAppend), channel_, variation.first, genMet.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptl", toAppend), channel_, variation.first, lep.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptl_smear", toAppend), channel_, variation.first, ptl_smear_fill, i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("etal", toAppend), channel_, variation.first, lep.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("phil", toAppend), channel_, variation.first, lep.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("nJets", toAppend), channel_, variation.first, jets.size(), i, thweight);
            SafeHistFill(weighthistMap2D_, concatenateNames("etal_ptl_2D", toAppend), channel_, variation.first, lep.eta(), lep.pt(), i, thweight);
            SafeHistFill(weighthistMap2D_, concatenateNames("etal_ptl_smear_2D", toAppend), channel_, variation.first, lep.eta(), ptl_smear_fill, i, thweight);
        }
    }

    if (((variation.first == ptV0to3 || variation.first == ptV0to3_lhe) && ptVcorr > 3.) ||
            ((variation.first == ptV3to5 || variation.first == ptV3to5_lhe) && (ptVcorr < 3. || ptVcorr > 5.))  ||
            ((variation.first == ptV5to7 || variation.first == ptV5to7_lhe) && (ptVcorr < 5. || ptVcorr > 7.)) ||
            ((variation.first == ptV7to9 || variation.first == ptV7to9_lhe) && (ptVcorr < 7. || ptVcorr > 9.)) ||
            ((variation.first == ptV9to12 || variation.first == ptV9to12_lhe) && (ptVcorr < 9. || ptVcorr > 12.)) ||
            ((variation.first == ptV12to15 || variation.first == ptV12to15_lhe) && (ptVcorr < 12. || ptVcorr > 15.)) ||
            ((variation.first == ptV15to20 || variation.first == ptV15to20_lhe) && (ptVcorr < 15. || ptVcorr > 20.)) ||
            ((variation.first == ptV20to27 || variation.first == ptV20to27_lhe) && (ptVcorr < 20. || ptVcorr > 27.)) ||
            ((variation.first == ptV27to40 || variation.first == ptV27to40_lhe) && (ptVcorr < 27. || ptVcorr > 40.)) ||
            ((variation.first == ptV40toInf || variation.first == ptV40toInf_lhe) && ptVcorr < 40. )) {
        return;
    }

    SafeHistFill(histMap1D_, concatenateNames("mW", toAppend), channel_, variation.first, wCand.mass(), weight);
    SafeHistFill(histMap1D_, concatenateNames("yW", toAppend), channel_, variation.first, wCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptW", toAppend), channel_, variation.first, wCand.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("mTtrue", toAppend), channel_, variation.first, mTtrue, weight);
    SafeHistFill(histMap1D_, concatenateNames("mTmet", toAppend), channel_, variation.first, mTmet, weight);
    SafeHistFill(histMap1D_, concatenateNames("mWmet", toAppend), channel_, variation.first, wCandMet.mass(), weight);
    SafeHistFill(histMap1D_, concatenateNames("yWmet", toAppend), channel_, variation.first, wCandMet.Rapidity(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptWmet", toAppend), channel_, variation.first, wCandMet.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("MET", toAppend), channel_, variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("MET_phi", toAppend), channel_, variation.first, genMet.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptl", toAppend), channel_, variation.first, lep.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptl_smear", toAppend), channel_, variation.first, ptl_smear_fill, weight);
    SafeHistFill(histMap1D_, concatenateNames("etal", toAppend), channel_, variation.first, lep.eta(), weight);
    SafeHistFill(histMap1D_, concatenateNames("phil", toAppend), channel_, variation.first, lep.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptnu", toAppend), channel_, variation.first, nu.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("etanu", toAppend), channel_, variation.first, nu.eta(), weight);
    SafeHistFill(histMap1D_, concatenateNames("phinu", toAppend), channel_, variation.first, nu.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("nJets", toAppend), channel_, variation.first, jets.size(), weight);
    SafeHistFill(histMap2D_, concatenateNames("etal_ptl_2D", toAppend), channel_, variation.first, lep.eta(), lep.pt(), weight);
    SafeHistFill(histMap2D_, concatenateNames("etal_ptl_smear_2D", toAppend), channel_, variation.first, lep.eta(), ptl_smear_fill, weight);
    SafeHistFill(histMap3D_, "mass_y_pT_3D", channel_, variation.first, wCand.mass(), wCand.Rapidity(), wCand.Pt(), weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, concatenateNames(("ptj"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, concatenateNames(("etaj"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, concatenateNames(("phij"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.phi(), weight);
        }  
    }

   if (variation.first == BareLeptons) {
       //ponerse las pilas, the (variation.first == BareLeptons) refraining the histograms to fill 
        ratio_mass /= wCand.mass();
        // Call the channel central just to avoid appending the name "barelep".
        // These histograms should only be built for the barelepton case, should be understood that they always refer 
        // to the barlep channel implicitly
        SafeHistFill(histMap1D_, "Ratio_Wmass", channel_, Central,  ratio_mass, weight);      
        SafeHistFill(histMap1D_, concatenateNames("nGammaAssoc",toAppend), channel_, Central, photons.size(), weight);

        auto compareByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() < b.pt(); };
        auto compareByDRLead = [lep] (const reco::GenParticle& a, const reco::GenParticle& b) {
            return reco::deltaR(a, lep) < reco::deltaR(b, lep);
        };

        auto gclose = std::min_element(photons.begin(), photons.end(), compareByDRLead);
        auto maxPtg = std::max_element(photons.begin(), photons.end(), compareByPt);
               
        

        SafeHistFill(histMap1D_, "dRlgamma_minassoc", channel_, Central, photons.size() > 0 ? reco::deltaR(*gclose, lep) : 0., weight);
        SafeHistFill(histMap1D_, "dRlgamma_maxptassoc", channel_, Central, photons.size() > 0 ? reco::deltaR(*maxPtg, lep) : 0., weight);
        SafeHistFill(histMap1D_, "ptg_closeassoc", channel_, Central, photons.size() > 0 ? gclose->pt() : 0., weight);
        SafeHistFill(histMap1D_, "ptgmax_assoc", channel_, Central, photons.size() > 0 ? maxPtg->pt() : 0., weight);
    }
}
