#include "Analysis/VVAnalysis/interface/WGenSelector.h"
#include <TRandom3.h>
#include <TStyle.h>
#include <cmath>
#include <regex>
#include "DataFormats/Math/interface/LorentzVector.h"

void WGenSelector::Init(TTree* tree) {
    histMap1D_[{"CutFlow", Unknown, Central}] = {};
    allChannels_ = {{ep, "ep"}, {en, "en"}, {mp, "mp"}, {mn, "mn"}};
    hists1D_ = {
        "CutFlow",
        "mWmet",
        "yWmet",
        "ptWmet",
        "mW",
        "yW",
        "ptW",
        "mTtrue",
        "mTmet",
        "ptl",
        "etal",
        "phil",
        "ptnu",
        "etanu",
        "phinu",
        "MET",
        "MET_phi",
        "ptj1",
        "ptj2",
        "etaj1",
        "etaj2",
        "nJets",
        "dRlgamma_maxptassoc",
        "dRlgamma_minassoc",
        "ptg_closeassoc",
        "ptg_maxassoc",
        "nGammaAssoc",
        "ptgmax_assoc",
        "ptgmax_assoc",
        "ptl_smear",
    };
    hists2D_ = {"etal_ptl_2D", "etal_ptl_smear_2D"};
    doSystematics_ = true;
    systematics_ = {
        {mWShift100MeVUp, "mWShift100MeVUp"},
        {mWShift50MeVUp, "mWShift50MeVUp"},
        {mWShift25MeVUp, "mWShift25MeVUp"},
        {mWShift20MeVUp, "mWShift20MeVUp"},
        {mWShift10MeVUp, "mWShift10MeVUp"},
        {mWShift100MeVDown, "mWShift100MeVDown"},
        {mWShift50MeVDown, "mWShift50MeVDown"},
        {mWShift25MeVDown, "mWShift25MeVDown"},
        {mWShift20MeVDown, "mWShift20MeVDown"},
        {mWShift10MeVUp, "mWShift10MeVUp"},
        {BareLeptons, "barelep"},
        {BornParticles, "born"},
        {LHEParticles, "lhe"},
        {muonScaleUp, "CMS_scale_mUp"},
        {muonScaleDown, "CMS_scale_mDown"},
    };
    systHists_ = hists1D_;
    systHists2D_ = hists2D_;

    weighthists1D_ = {
        "CutFlow", "mW",        "mTmet", "yW",   "ptW",
        "ptl",     "ptl_smear", "etal",  "ptnu", "etanu",
    };
    weighthists2D_ = hists2D_;

    refWeight = 1;
    nLeptons_ = 1;
    doNeutrinos_ = true;
    doPhotons_ = true;

    // Chose by MC sample
    if (name_.find("nnlops") != std::string::npos) {
        MV_GEN_ = 80398.0;
        GAMMAV_GEN_ = 2088.720;
    } else if (name_.find("minnlo") != std::string::npos) {
        MV_GEN_ = 80351.97159;
        GAMMAV_GEN_ = 2084.29889;
    } else {
        MV_GEN_ = 80419.;
        GAMMAV_GEN_ = 2050;
    }

    NanoGenSelectorBase::Init(tree);
}

void WGenSelector::LoadBranchesNanoAOD(Long64_t entry, SystPair variation) {
    NanoGenSelectorBase::LoadBranchesNanoAOD(entry, variation);

    if (leptons.size() >= nLeptons_) {
        auto& l = leptons.at(0);
        if (variation.first == Central) {
            TRandom3 gauss;
            ptl_smear = l.pt() * gauss.Gaus(1, 0.01);
        } else if (variation.first == muonScaleUp) {
            leptons.at(0).setP4(makeGenParticle(l.pdgId(), l.status(),
                                                l.pt() * 1.001, l.eta(),
                                                l.phi(), l.mass())
                                    .polarP4());
            SetComposite();
        } else if (variation.first == muonScaleDown) {
            leptons.at(0).setP4(makeGenParticle(l.pdgId(), l.status(),
                                                l.pt() * 1. / 1.001, l.eta(),
                                                l.phi(), l.mass())
                                    .polarP4());
            SetComposite();
        }
    }

    if (variation.first == Central)
        cenWeight = weight;
    else if (variation.first == LHEParticles) {
        ptVlhe = wCand.pt();
        mVlhe = wCand.mass() * 1000.;
    } else if (variation.first == mWShift10MeVUp)
        weight = cenWeight * breitWignerWeight(10.);
    else if (variation.first == mWShift10MeVDown)
        weight = cenWeight * breitWignerWeight(-10.);
    else if (variation.first == mWShift20MeVUp)
        weight = cenWeight * breitWignerWeight(20.);
    else if (variation.first == mWShift20MeVDown)
        weight = cenWeight * breitWignerWeight(-20.);
    else if (variation.first == mWShift25MeVUp)
        weight = cenWeight * breitWignerWeight(25.);
    else if (variation.first == mWShift25MeVDown)
        weight = cenWeight * breitWignerWeight(-25.);
    else if (variation.first == mWShift50MeVUp)
        weight = cenWeight * breitWignerWeight(50.);
    else if (variation.first == mWShift50MeVDown)
        weight = cenWeight * breitWignerWeight(-50.);
    else if (variation.first == mWShift100MeVUp)
        weight = cenWeight * breitWignerWeight(100.);
    else if (variation.first == mWShift100MeVDown)
        weight = cenWeight * breitWignerWeight(-100.);

    if (leptons.size() > 0 && std::abs(leptons.at(0).pdgId()) == 11) {
        if (leptons.at(0).pdgId() > 0) {
            channel_ = en;
            channelName_ = "en";
        } else {
            channel_ = ep;
            channelName_ = "ep";
        }
    } else if (leptons.size() > 0 && std::abs(leptons.at(0).pdgId()) == 13) {
        if (leptons.at(0).pdgId() > 0) {
            channel_ = mn;
            channelName_ = "mn";
        } else {
            channel_ = mp;
            channelName_ = "mp";
        }
    } else {
        channel_ = Unknown;
        channelName_ = "Unknown";
        return;
    }
}

void WGenSelector::SetComposite() {
    if (leptons.size() == 0) {
        wCandMet = LorentzVector();
        wCand = LorentzVector();
        return;
    } else if (neutrinos.size() == 0) {
        wCand = LorentzVector();
        return;
    }
    auto lepP4 = leptons.at(0).polarP4();
    auto compareByPt = [](const reco::GenParticle& a,
                          const reco::GenParticle& b) {
        return a.pt() < b.pt();
    };
    auto mt = [](LorentzVector& l, LorentzVector& v) {
        return std::sqrt(2 * l.pt() * v.pt() * (1 - cos(l.phi() - v.phi())));
    };

    auto nup =
        std::max_element(neutrinos.begin(), neutrinos.end(), compareByPt);
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

void WGenSelector::FillHistogramsByName(Long64_t entry, std::string& toAppend,
                                        SystPair variation) {
    int step = 0;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_,
                 variation.first, step++, weight);

    if (channel_ != mn && channel_ != en && channel_ != mp && channel_ != ep)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_,
                 variation.first, step++, weight);

    if (leptons.size() < nLeptons_) return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_,
                 variation.first, step++, weight);

    if (neutrinos.size() < nLeptons_) return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_,
                 variation.first, step++, weight);

    auto& lep = leptons.at(0);
    if (doFiducial_ && std::abs(lep.eta()) > 2.5) return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_,
                 variation.first, step++, weight);

    if (variation.first == Central)
        mcWeights_->Fill(weight / std::abs(refWeight));

    if (doFiducial_ && lep.pt() < 25) return;

    float ptl_smear_fill = ptl_smear;
    if (variation.first == muonScaleUp)
        ptl_smear_fill *= 1.001;
    else if (variation.first == muonScaleDown)
        ptl_smear_fill *= 0.999;

    if (std::find(theoryVarSysts_.begin(), theoryVarSysts_.end(),
                  variation.first) != theoryVarSysts_.end()) {
        size_t minimalWeights = *nLHEScaleWeight + nLHEScaleWeightAltSet1 +
                                nLHEUnknownWeight + nLHEUnknownWeightAltSet1;
        size_t nWeights = variation.first == Central
                              ? minimalWeights
                              : minimalWeights + nLHEPdfWeight;

        for (size_t i = 0; i < nWeights; i++) {
            float thweight = 1;
            if (i < *nLHEScaleWeight)
                thweight = LHEScaleWeight[i];
            else if (i < *nLHEScaleWeight + nLHEScaleWeightAltSet1)
                thweight = LHEScaleWeightAltSet1[i - *nLHEScaleWeight];
            else if (i < minimalWeights - nLHEUnknownWeightAltSet1)
                thweight =
                    LHEUnknownWeight[i - minimalWeights + nLHEUnknownWeight +
                                     nLHEUnknownWeightAltSet1];
            else if (i < minimalWeights)
                thweight = LHEUnknownWeight[i - minimalWeights +
                                            nLHEUnknownWeightAltSet1];
            else
                thweight = LHEPdfWeight[i - minimalWeights];

            if (centralWeightIndex_ != -1)
                thweight /= LHEScaleWeight.At(centralWeightIndex_);

            if (((variation.first == ptV0to3 ||
                  variation.first == ptV0to3_lhe) &&
                 ptVlhe > 3.) ||
                ((variation.first == ptV3to5 ||
                  variation.first == ptV3to5_lhe) &&
                 (ptVlhe < 3. || ptVlhe > 5.)) ||
                ((variation.first == ptV5to7 ||
                  variation.first == ptV5to7_lhe) &&
                 (ptVlhe < 5. || ptVlhe > 7.)) ||
                ((variation.first == ptV7to9 ||
                  variation.first == ptV7to9_lhe) &&
                 (ptVlhe < 7. || ptVlhe > 9.)) ||
                ((variation.first == ptV9to12 ||
                  variation.first == ptV9to12_lhe) &&
                 (ptVlhe < 9. || ptVlhe > 12.)) ||
                ((variation.first == ptV12to15 ||
                  variation.first == ptV12to15_lhe) &&
                 (ptVlhe < 12. || ptVlhe > 15.)) ||
                ((variation.first == ptV15to20 ||
                  variation.first == ptV15to20_lhe) &&
                 (ptVlhe < 15. || ptVlhe > 20.)) ||
                ((variation.first == ptV20to27 ||
                  variation.first == ptV20to27_lhe) &&
                 (ptVlhe < 20. || ptVlhe > 27.)) ||
                ((variation.first == ptV27to40 ||
                  variation.first == ptV27to40_lhe) &&
                 (ptVlhe < 27. || ptVlhe > 40.)) ||
                ((variation.first == ptV40toInf ||
                  variation.first == ptV40toInf_lhe) &&
                 ptVlhe < 40.)) {
                thweight = 1;
            }

            thweight *= weight;
            SafeHistFill(weighthistMap1D_, concatenateNames("mW", toAppend),
                         channel_, variation.first, wCand.mass(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("yW", toAppend),
                         channel_, variation.first, wCand.Rapidity(), i,
                         thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptW", toAppend),
                         channel_, variation.first, wCand.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("mWmet", toAppend),
                         channel_, variation.first, wCandMet.mass(), i,
                         thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("yWmet", toAppend),
                         channel_, variation.first, wCandMet.Rapidity(), i,
                         thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptWmet", toAppend),
                         channel_, variation.first, wCandMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("MET", toAppend),
                         channel_, variation.first, genMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_,
                         concatenateNames("MET_phi", toAppend), channel_,
                         variation.first, genMet.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptl", toAppend),
                         channel_, variation.first, lep.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_,
                         concatenateNames("ptl_smear", toAppend), channel_,
                         variation.first, ptl_smear_fill, i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("etal", toAppend),
                         channel_, variation.first, lep.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("phil", toAppend),
                         channel_, variation.first, lep.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("nJets", toAppend),
                         channel_, variation.first, jets.size(), i, thweight);
            SafeHistFill(weighthistMap2D_,
                         concatenateNames("etal_ptl_2D", toAppend), channel_,
                         variation.first, lep.eta(), lep.pt(), i, thweight);
            SafeHistFill(weighthistMap2D_,
                         concatenateNames("etal_ptl_smear_2D", toAppend),
                         channel_, variation.first, lep.eta(), ptl_smear_fill,
                         i, thweight);
        }
    }

    if (((variation.first == ptV0to3 || variation.first == ptV0to3_lhe) &&
         ptVlhe > 3.) ||
        ((variation.first == ptV3to5 || variation.first == ptV3to5_lhe) &&
         (ptVlhe < 3. || ptVlhe > 5.)) ||
        ((variation.first == ptV5to7 || variation.first == ptV5to7_lhe) &&
         (ptVlhe < 5. || ptVlhe > 7.)) ||
        ((variation.first == ptV7to9 || variation.first == ptV7to9_lhe) &&
         (ptVlhe < 7. || ptVlhe > 9.)) ||
        ((variation.first == ptV9to12 || variation.first == ptV9to12_lhe) &&
         (ptVlhe < 9. || ptVlhe > 12.)) ||
        ((variation.first == ptV12to15 || variation.first == ptV12to15_lhe) &&
         (ptVlhe < 12. || ptVlhe > 15.)) ||
        ((variation.first == ptV15to20 || variation.first == ptV15to20_lhe) &&
         (ptVlhe < 15. || ptVlhe > 20.)) ||
        ((variation.first == ptV20to27 || variation.first == ptV20to27_lhe) &&
         (ptVlhe < 20. || ptVlhe > 27.)) ||
        ((variation.first == ptV27to40 || variation.first == ptV27to40_lhe) &&
         (ptVlhe < 27. || ptVlhe > 40.)) ||
        ((variation.first == ptV40toInf || variation.first == ptV40toInf_lhe) &&
         ptVlhe < 40.)) {
        return;
    }

    SafeHistFill(histMap1D_, concatenateNames("mW", toAppend), channel_,
                 variation.first, wCand.mass(), weight);
    SafeHistFill(histMap1D_, concatenateNames("yW", toAppend), channel_,
                 variation.first, wCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptW", toAppend), channel_,
                 variation.first, wCand.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("mTtrue", toAppend), channel_,
                 variation.first, mTtrue, weight);
    SafeHistFill(histMap1D_, concatenateNames("mTmet", toAppend), channel_,
                 variation.first, mTmet, weight);
    SafeHistFill(histMap1D_, concatenateNames("mWmet", toAppend), channel_,
                 variation.first, wCandMet.mass(), weight);
    SafeHistFill(histMap1D_, concatenateNames("yWmet", toAppend), channel_,
                 variation.first, wCandMet.Rapidity(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptWmet", toAppend), channel_,
                 variation.first, wCandMet.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("MET", toAppend), channel_,
                 variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("MET_phi", toAppend), channel_,
                 variation.first, genMet.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptl", toAppend), channel_,
                 variation.first, lep.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptl_smear", toAppend), channel_,
                 variation.first, ptl_smear_fill, weight);
    SafeHistFill(histMap1D_, concatenateNames("etal", toAppend), channel_,
                 variation.first, lep.eta(), weight);
    SafeHistFill(histMap1D_, concatenateNames("phil", toAppend), channel_,
                 variation.first, lep.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptnu", toAppend), channel_,
                 variation.first, nu.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("etanu", toAppend), channel_,
                 variation.first, nu.eta(), weight);
    SafeHistFill(histMap1D_, concatenateNames("phinu", toAppend), channel_,
                 variation.first, nu.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("nJets", toAppend), channel_,
                 variation.first, jets.size(), weight);
    SafeHistFill(histMap2D_, concatenateNames("etal_ptl_2D", toAppend),
                 channel_, variation.first, lep.eta(), lep.pt(), weight);
    SafeHistFill(histMap2D_, concatenateNames("etal_ptl_smear_2D", toAppend),
                 channel_, variation.first, lep.eta(), ptl_smear_fill, weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i) {
            const auto& jet = jets.at(i - 1);
            SafeHistFill(
                histMap1D_,
                concatenateNames(("ptj" + std::to_string(i)).c_str(), toAppend),
                channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_,
                         concatenateNames(("etaj" + std::to_string(i)).c_str(),
                                          toAppend),
                         channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_,
                         concatenateNames(("phij" + std::to_string(i)).c_str(),
                                          toAppend),
                         channel_, variation.first, jet.phi(), weight);
        }
    }

    if (variation.first == BareLeptons) {
        SafeHistFill(histMap1D_, "nGammaAssoc", channel_, variation.first,
                     photons.size(), weight);

        auto compareByPt = [](const reco::GenParticle& a,
                              const reco::GenParticle& b) {
            return a.pt() < b.pt();
        };
        auto compareByDRLead = [lep](const reco::GenParticle& a,
                                     const reco::GenParticle& b) {
            return reco::deltaR(a, lep) < reco::deltaR(b, lep);
        };

        auto gclose =
            std::min_element(photons.begin(), photons.end(), compareByDRLead);
        auto maxPtg =
            std::max_element(photons.begin(), photons.end(), compareByPt);

        SafeHistFill(histMap1D_, "dRlgamma_minassoc", channel_, variation.first,
                     photons.size() > 0 ? reco::deltaR(*gclose, lep) : 0.,
                     weight);
        SafeHistFill(
            histMap1D_, "dRlgamma_maxptassoc", channel_, variation.first,
            photons.size() > 0 ? reco::deltaR(*maxPtg, lep) : 0., weight);
        SafeHistFill(histMap1D_, "ptg_closeassoc", channel_, variation.first,
                     photons.size() > 0 ? gclose->pt() : 0., weight);
        SafeHistFill(histMap1D_, "ptgmax_assoc", channel_, variation.first,
                     photons.size() > 0 ? maxPtg->pt() : 0., weight);
    }
}
