#include "Analysis/SelectorTools/interface/ZGenSelector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <TStyle.h>
#include <regex>
#include <numeric>
#include "TParameter.h"

void ZGenSelector::Init(TTree *tree)
{
    // Don't waste memory on empty e hists
    TParameter<bool>* muOnlyParam = (TParameter<bool>*) GetInputList()->FindObject("muOnly");
    bool muOnly = muOnlyParam != nullptr && muOnlyParam->GetVal();
    allChannels_ = {{mm, "mm"}}; 
    if (!muOnly) {
        allChannels_.push_back(std::make_pair<Channel, std::string>(ee, "ee"));
    }
    // Add CutFlow for Unknown to understand when channels aren't categorized
    histMap1D_[{"CutFlow", Unknown, Central}] = {};
    std::vector<std::string> basehists1D = {"CutFlow", "ZMass", "yZ", "ptZ", "phiZ", "ptl1", "etal1", "phil1", "ptl2", "etal2", "phil2", 
        "ptj1", "ptj2", "ptj3", "etaj1", "etaj2", "etaj3", "phij1", "phij2", "phij3", "nJets",
        "MET", "HT",};
    hists1D_ = basehists1D;
    //std::vector<std::string> partonicChans = {"uu_dd", "uubar_ddbar", "ug_dg", "ubarg_dbarg", "gg", "other"};
    //for (auto& chan : partonicChans) {
    //    for (auto& hist : basehists1D)
    //        hists1D_.push_back(chan + "_" + hist);
    //}
    systHists_ = hists1D_;
    hists3D_ = {"mass_y_pT_3D", };
    systHists3D_ = hists3D_;

    weighthists1D_ = {"CutFlow", "ZMass", "yZ", "ptZ", "phiZ", "ptl1", "etal1", "ptl2", "etal2", 
        "ptj1", "ptj2", "ptj3", "etaj1", "etaj2", "etaj3", "nJets",
        "MET", "HT", };
    nLeptons_ = 2;

    TParameter<bool>* massVar = (TParameter<bool>*) GetInputList()->FindObject("massVar");
    doMassVar_ = massVar != nullptr && massVar->GetVal();

    if (doMassVar_) {
        systematics_[mZShift50MeVUp] = "mZShift50MeVUp";
        systematics_[mZShift50MeVDown] = "mZShift50MeVDown";
        systematics_[mZShift100MeVUp] = "mZShift100MeVUp";
        systematics_[mZShift100MeVDown] = "mZShift100MeVDown";
    }
    doSystematics_ = !systematics_.empty();

    if (name_.find("nnlops") != std::string::npos) {
        MV_GEN_ = 80398.0;
        GAMMAV_GEN_ = 2088.720;
    }
    else if (name_.find("minnlo") != std::string::npos) {
        MV_GEN_ = 91153.509740726733;
        GAMMAV_GEN_ = 2493.2018986110700;
    }
    else {
        MV_GEN_ = 80419.;
        GAMMAV_GEN_ = 2050;
    }

    NanoGenSelectorBase::Init(tree);
    if (name_.find("Corr") != std::string::npos) {
        n3llcorr_ = true;
        SetScaleFactors();
        if (scetlibCorrs_[0] == nullptr)
            throw std::invalid_argument("Must pass a scalefactor for sample with corrections!");
    }
}

void ZGenSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) {
    NanoGenSelectorBase::LoadBranchesNanoAOD(entry, variation);

    if (variation.first == Central) {
        cenWeight = weight;
        ptVcorr = zCand.pt();
        mVcorr = zCand.mass();
        ratio_mass = zCand.mass();
    }
    else if (variation.first == PreFSRLeptons) { 
            ptVcorr = zCand.pt();
            mVcorr = zCand.mass();
            yVcorr = zCand.Rapidity();
    }
    else if (variation.first == LHEParticles && systematics_.find(PreFSRLeptons) == std::end(systematics_)) {
        // define at LHE level if it exists
        ptVcorr = zCand.pt();
        mVcorr = zCand.mass();
        yVcorr = zCand.Rapidity();
    }
    else if (variation.first == mZShift10MeVUp)
        weight = cenWeight*breitWignerWeight(10.);
    else if (variation.first == mZShift10MeVDown)
        weight = cenWeight*breitWignerWeight(-10.);
    else if (variation.first == mZShift20MeVUp)
        weight = cenWeight*breitWignerWeight(20.);
    else if (variation.first == mZShift20MeVDown)
        weight = cenWeight*breitWignerWeight(-20.);
    else if (variation.first == mZShift25MeVUp)
        weight = cenWeight*breitWignerWeight(25.);
    else if (variation.first == mZShift25MeVDown)
        weight = cenWeight*breitWignerWeight(-25.);
    else if (variation.first == mZShift50MeVUp)
        weight = cenWeight*breitWignerWeight(50.);
    else if (variation.first == mZShift50MeVDown)
        weight = cenWeight*breitWignerWeight(-50.);
    else if (variation.first == mZShift100MeVUp)
        weight = cenWeight*breitWignerWeight(100.);
    else if (variation.first == mZShift100MeVDown)
        weight = cenWeight*breitWignerWeight(-100.);

    if (leptons.size() < 2) {
        channel_ = Unknown;
        channelName_ = "Unknown";
        return;
    }
    if (leptons.at(0).pdgId() + leptons.at(1).pdgId() == 0) {
        if (std::abs(leptons.at(0).pdgId()) == 11) {
            channel_ = ee;
            channelName_ = "ee";
        }
        else if (std::abs(leptons.at(0).pdgId()) == 13) {
            channel_ = mm;
            channelName_ = "mm";
        }
    }
    else {
        channel_ = Unknown;
        channelName_ = "Unknown";
    }
    if (n3llcorr_) {
        weight *= scetlibCorrs_[0]->Evaluate3D(mVcorr, yVcorr, ptVcorr);
    }
}

void ZGenSelector::SetComposite() {
    if (leptons.size() >= 2)
        zCand = leptons.at(0).polarP4() + leptons.at(1).polarP4();
}

void ZGenSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    int step = 0;
    int failStep = 0;

    step++;
    if (channel_ != mm && channel_ != ee) 
        failStep = step;

    auto lep1 = leptons.size() > 1 ? leptons.at(0) : reco::GenParticle();
    auto lep2 = leptons.size() > 1 ? leptons.at(1) : reco::GenParticle();

    step++;
    if (zCand.mass() < 50.)
        failStep = step;
    step++;
    if (lep1.pt() < 25. || lep2.pt() < 25.)
        failStep = step;
    step++;
    if ((std::abs(lep1.eta()) > 2.5 || std::abs(lep2.eta()) > 2.5) ||
            (selection_ == ZselectionTight && (std::abs(lep1.eta()) > 2.4 || std::abs(lep2.eta()) > 2.4)))
        failStep = step;
    step++;
    if ((zCand.mass() < 60. || zCand.mass() > 120.) || (selection_ == ZselectionTight && (zCand.mass() < 76.1876 || zCand.mass() > 106.1786)))
        failStep = step;

    for (int j = 0; j < (failStep == 0 ? step : failStep); j++) {
        SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, j, weight);
        //size_t nWeights = *nLHEScaleWeight;
        //for (size_t i = 0; i < nWeights; i++) {
        //    float thweight = LHEScaleWeight[i];
        //    thweight *= weight;
        //    SafeHistFill(weighthistMap1D_, "CutFlow", channel_, variation.first, j, i, thweight);
        //}
    }
    if (doFiducial_ && failStep != 0)
        return;
    if (variation.first == Central)
        mcWeights_->Fill(weight/std::abs(refWeight));

    if (std::find(theoryVarSysts_.begin(), theoryVarSysts_.end(), variation.first) != theoryVarSysts_.end()) {
        size_t nScaleWeights = nLHEScaleWeight+nLHEScaleWeightAltSet1;
        size_t minimalWeights = nLHEScaleWeight+nLHEScaleWeightAltSet1+nLHEUnknownWeight+nLHEUnknownWeightAltSet1;
        size_t allPdfWeights = std::accumulate(nLHEPdfWeights.begin(), nLHEPdfWeights.end(), 1);

        //size_t nWeights = variation.first == Central ? minimalWeights+allPdfWeights : minimalWeights;
        size_t nWeights = minimalWeights+allPdfWeights;
        size_t pdfOffset = nScaleWeights;
        size_t pdfIdx = 0;
        if (n3llcorr_)
            nWeights += nScetlibWeights_;
        for (size_t i = 0; i < nWeights; i++) {
            float thweight = 1;
            if (i < nLHEScaleWeight)
                thweight = LHEScaleWeight[i];
            else if (i < nScaleWeights)
                thweight = LHEScaleWeightAltSet1[i-nLHEScaleWeight];
            else if (i < nScaleWeights+allPdfWeights) {
                thweight = LHEPdfWeights[pdfIdx][i-pdfOffset];
                if (i == pdfOffset+nLHEPdfWeights.at(pdfIdx)-1) {
                    pdfOffset += nLHEPdfWeights.at(pdfIdx++);
                }
            }
            else {
                int idx = i-nScaleWeights-allPdfWeights;
                auto* sf = scetlibCorrs_.at(idx);
                float refW = scetlibCorrs_.at(0)->Evaluate3D(mVcorr, yVcorr, ptVcorr);
                thweight = sf->Evaluate3D(mVcorr, yVcorr, ptVcorr)/refW;
            }

            thweight /= rescaleWeight_;

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

            thweight *= weight;
            SafeHistFill(weighthistMap1D_, "ZMass", channel_, variation.first, zCand.mass(), i, thweight);
            SafeHistFill(weighthistMap1D_, "yZ", channel_, variation.first, zCand.Rapidity(), i, thweight);
            SafeHistFill(weighthistMap1D_, "ptZ", channel_, variation.first, zCand.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, "phiZ", channel_, variation.first, zCand.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, "ptl1", channel_, variation.first, lep1.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, "etal1", channel_, variation.first, lep1.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, "phil1", channel_, variation.first, lep1.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, "ptl2", channel_, variation.first, lep2.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, "etal2", channel_, variation.first, lep2.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, "phil2", channel_, variation.first, lep2.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, "nJets", channel_, variation.first, jets.size(), i, thweight);
            SafeHistFill(weighthistMap1D_, "MET", channel_, variation.first, genMet.pt(), i, thweight);
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

    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);
    SafeHistFill(histMap1D_, "ZMass", channel_, variation.first, zCand.mass(), weight);
    SafeHistFill(histMap1D_, "yZ", channel_, variation.first, zCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, "ptZ", channel_, variation.first, zCand.pt(), weight);
    SafeHistFill(histMap1D_, "phiZ", channel_, variation.first, zCand.phi(), weight);
    SafeHistFill(histMap1D_, "ptl1", channel_, variation.first, lep1.pt(), weight);
    SafeHistFill(histMap1D_, "etal1", channel_, variation.first, lep1.eta(), weight);
    SafeHistFill(histMap1D_, "phil1", channel_, variation.first, lep1.phi(), weight);
    SafeHistFill(histMap1D_, "ptl2", channel_, variation.first, lep2.pt(), weight);
    SafeHistFill(histMap1D_, "etal2", channel_, variation.first, lep2.eta(), weight);
    SafeHistFill(histMap1D_, "phil2", channel_, variation.first, lep2.phi(), weight);
    SafeHistFill(histMap1D_, "nJets", channel_, variation.first, jets.size(), weight);
    SafeHistFill(histMap1D_, "MET", channel_, variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, "HT", channel_, variation.first, ht, weight);
    SafeHistFill(histMap3D_, "mass_y_pT_3D", channel_, variation.first, zCand.mass(), zCand.Rapidity(), zCand.Pt(), weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, ("ptj"+std::to_string(i)).c_str(), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, ("etaj"+std::to_string(i)).c_str(), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, ("phij"+std::to_string(i)).c_str(), channel_, variation.first, jet.phi(), weight);
        }  
    }
    
    // Should check how slow this is. For now it's off 
    return;

    std::string partonicChan = "other";
    if ((*Generator_id1 == 1 && *Generator_id2 == 1) || (*Generator_id1 == 2 && *Generator_id2 == 2))
        partonicChan = "uu_dd";
    else if ((*Generator_id1 == 1 && *Generator_id2 == -1) || (*Generator_id1 == 2 && *Generator_id2 == -2))
        partonicChan = "uubar_ddbar";
    else if (*Generator_id1 == 21 && *Generator_id2 == 21)
        partonicChan = "gg";
    else if ((*Generator_id1 == 1 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == 1) || 
                (*Generator_id1 == 2 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == 2))
        partonicChan = "ug_dg";
    else if ((*Generator_id1 == -1 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == -1) || 
                (*Generator_id1 == -2 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == -2))
        partonicChan = "ubarg_dbarg";
    SafeHistFill(histMap1D_, (partonicChan+"_ZMass").c_str(), channel_, variation.first, zCand.mass(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_yZ").c_str(), channel_, variation.first, zCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_ptZ").c_str(), channel_, variation.first, zCand.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_ptl1").c_str(), channel_, variation.first, lep1.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_etal1").c_str(), channel_, variation.first, lep1.eta(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_phil1").c_str(), channel_, variation.first, lep1.phi(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_ptl2").c_str(), channel_, variation.first, lep2.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_etal2").c_str(), channel_, variation.first, lep2.eta(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_phil2").c_str(), channel_, variation.first, lep2.phi(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_nJets").c_str(), channel_, variation.first, jets.size(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_MET").c_str(), channel_, variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_HT").c_str(), channel_, variation.first, ht, weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, (partonicChan+"_ptj"+std::to_string(i)).c_str(), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, (partonicChan+"_etaj"+std::to_string(i)).c_str(), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, (partonicChan+"_phij"+std::to_string(i)).c_str(), channel_, variation.first, jet.phi(), weight);
        }  
    }
}

