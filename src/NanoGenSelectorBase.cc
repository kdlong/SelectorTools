#include "Analysis/SelectorTools/interface/NanoGenSelectorBase.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TParameter.h"
#include <TStyle.h>
#include <regex>
#include <numeric>

void NanoGenSelectorBase::Init(TTree *tree)
{
    refWeight = 1;
    // Otherwise name isn't read until the base class is called
    TNamed* name = (TNamed *) GetInputList()->FindObject("name");
    bool isMinnlo_ = std::string(name->GetTitle()).find("minnlo") != std::string::npos;

    TParameter<bool>* doTheory = (TParameter<bool>*) GetInputList()->FindObject("theoryUnc");
    doTheoryVars_ = doTheory != nullptr && doTheory->GetVal();
    TParameter<bool>* theoryPrefsr = (TParameter<bool>*) GetInputList()->FindObject("theoryPrefsr");
    bool doTheoryPrefsr = theoryPrefsr != nullptr && theoryPrefsr->GetVal();
    TParameter<bool>* theoryLhe = (TParameter<bool>*) GetInputList()->FindObject("theoryLhe");
    bool doTheoryLhe = theoryLhe != nullptr && theoryLhe->GetVal();
    TParameter<bool>* lhePart = (TParameter<bool>*) GetInputList()->FindObject("lhe");
    doLHE_ = lhePart != nullptr && lhePart->GetVal();
    TParameter<bool>* prefsrPart = (TParameter<bool>*) GetInputList()->FindObject("prefsr");
    doPreFSR_ = prefsrPart != nullptr && prefsrPart->GetVal();
    TParameter<bool>* bornPart = (TParameter<bool>*) GetInputList()->FindObject("born");
    doBorn_ = bornPart != nullptr && bornPart->GetVal();
    TParameter<bool>* barePart = (TParameter<bool>*) GetInputList()->FindObject("bare");
    doBareLeptons_ = barePart != nullptr && barePart->GetVal();

    TParameter<bool>* wSignOnly = (TParameter<bool>*) GetInputList()->FindObject("wSignOnly");
    weightSignOnly_ = wSignOnly != nullptr && wSignOnly->GetVal();
    TParameter<int>* wSuppress = (TParameter<int>*) GetInputList()->FindObject("wSuppress");
    weightSuppress_ = wSuppress != nullptr ? wSuppress->GetVal() : 0;
    TParameter<int>* thwSuppress = (TParameter<int>*) GetInputList()->FindObject("thwSuppress");
    // Default to 10x the nominal
    thweightSuppress_ = thwSuppress != nullptr ? thwSuppress->GetVal() : 10;

    TParameter<int>* storeCenPdfs = (TParameter<int>*) GetInputList()->FindObject("storeCenPdfs");
    storeCenPdfs_ = storeCenPdfs != nullptr ? storeCenPdfs->GetVal() : 0;
    // Only valid for isMinnlo
    storeCenPdfs_ = storeCenPdfs_ && isMinnlo_;
    if (storeCenPdfs_)
        std::cout << "INFO: Also storing central values of main PDF sets\n";

    TNamed* pdfSet = (TNamed*) GetInputList()->FindObject("pdfSet");
    pdfSet_ = pdfSet != nullptr ? pdfSet->GetTitle() : pdfSet_;

    std::cout << "INFO: doLHE = " << doLHE_ << " doPrefsr " << doPreFSR_ << std::endl;
    std::cout << "INFO: doBareLeptons = "<<doBareLeptons_<<"\n";

    if (isMinnlo_ && !weightSuppress_ && !weightSignOnly_) {
        std::cout << "WARNING: You should use wSuppress != 0 or wSignOnly to suppress huge weights in MiNNLO\n";
    }
    if (weightSuppress_ && !weightSignOnly_)
        std::cout << "INFO: wSuppress = " << weightSuppress_ << std::endl;

    if (doBorn_)
        systematics_[BornParticles] = "born";
    if (doBareLeptons_) {
        systematics_[BareLeptons] = "bare";
    }
    if (doPreFSR_)
        systematics_[PreFSRLeptons] = "prefsr";
    if (doLHE_)
        systematics_[LHEParticles] = "lhe";

    if (doTheoryVars_) {
        theoryVarSysts_.insert(theoryVarSysts_.begin(), Central);
    }
    if (doTheoryPrefsr) {
        theoryVarSysts_.insert(theoryVarSysts_.end(), PreFSRLeptons);
    }
    if (doTheoryLhe) {
        theoryVarSysts_.insert(theoryVarSysts_.end(), LHEParticles);
    }
    if (doBareLeptons_) {
        theoryVarSysts_.insert(theoryVarSysts_.end(), BareLeptons);
    }
    TParameter<bool>* doPtVSplit = (TParameter<bool>*) GetInputList()->FindObject("theoryPtV");
    if (doTheoryVars_ && doPtVSplit != nullptr && doPtVSplit->GetVal()) {
        // For now it's based on the LHE kinematics
        TNamed* selection = (TNamed *) GetInputList()->FindObject("selection");
        std::string selectionName = selection == nullptr ? "" : selection->GetTitle();
        SystMap ptvars = {
            {ptV0to3, "ptV0to3"},
            {ptV3to5, "ptV3to5"},
            {ptV5to7, "ptV5to7"},
            {ptV7to9, "ptV7to9"},
            {ptV9to12, "ptV9to12"},
            {ptV12to15, "ptV12to15"},
            {ptV15to20, "ptV15to20"},
            {ptV20to27, "ptV20to27"},
            {ptV27to40, "ptV27to40"},
            {ptV40toInf, "ptV40toInf"},
        };
        for (auto& var : ptvars) {
            systematics_[var.first] = var.second;
            theoryVarSysts_.insert(theoryVarSysts_.end(), var.first);
        }
    }
    doSystematics_ = !systematics_.empty();

    b.SetTree(tree);
    scaleWeights_ = (tree->GetListOfBranches()->FindObject("nLHEScaleWeight") != nullptr);
    altScaleWeights_ = (tree->GetListOfBranches()->FindObject("nLHEScaleWeightAltSet1") != nullptr);
    paramWeights_ = (tree->GetListOfBranches()->FindObject("nMEParamWeight") != nullptr);
    paramWeightsUpd_ = (tree->GetListOfBranches()->FindObject("nLHEReweightingWeightCorrectMass") != nullptr);
    unknownWeights_ = (tree->GetListOfBranches()->FindObject("nLHEUnknownWeight") != nullptr);
    unknownWeightsAlt_ = (tree->GetListOfBranches()->FindObject("nLHEUnknownWeightAltSet1") != nullptr);
    
    bool readAllPdfs = pdfSet_ == "all";
    if (!readAllPdfs && !isMinnlo_)
        std::cerr << "WARNING! Specific pdf selection '" << pdfSet_ << "' is only supported in MiNNLO. No PDF weights will be stored\n";
    else if (isMinnlo_ && pdfSet_ != "all") {
        if (pdfSet_ == "ct18")
            // Don't store CT18Z central value, or central value for your own set
            pdfMaxStore_ = 61+storeCenPdfs_*cenPdfOffset_;
        else if (pdfSet_ == "ct18z")
            pdfCenWeight_ = 61;
    }
    // Should always be the first entry, but CT18 is accidentally 
    // In the same set as CT18Z
    for (size_t i = 0; i < pdfWeights_.size(); i++) {
        if (!readAllPdfs && minnloPdfMap.find(pdfSet_) == std::end(minnloPdfMap))
            std::cerr << "WARNING! PDF set " << pdfSet_ << " is invalid! It will be skipped.\n";
        auto& toRead = minnloPdfMap[pdfSet_];
        bool readPdf = doTheoryVars_ && (readAllPdfs || (isMinnlo_ && std::find(std::begin(toRead), std::end(toRead), i) != std::end(toRead)));
        std::string name = "LHEPdfWeight";
        if (i > 0)
            name += "AltSet" + std::to_string(i);
        if (!readPdf) {
            continue;
        }
        altPdf_ = true;
        if(tree->GetListOfBranches()->FindObject(name.c_str()) != nullptr) {
            std::cout << "INFO: Storing pdf set read from " << name << std::endl;
            pdfWeights_[i] = true;
        }
    }
    // Trigger storing all relevant
    if (storeCenPdfs_) {
        for (auto& entry : minnloPdfMap) {
            // Need to skip it since its central value isn't at index 0 (merged with CT18)
            if (entry.first == "ct18z")
                continue;
            // Just store the first entry, most of the authors are alphas vars
            int i = entry.second.front();
            pdfWeights_[i] = true;
            if (entry.first != pdfSet_)
                cenPdfWeightsToStore_.push_back(i);
        }
    }

    nTheoryWeights_ = pdfSet_ == "all" ? 1000 : 200;

    SelectorBase::Init(tree);

    // Off for now
    doMC2H_ = false && name_.find("cp5") == std::string::npos;
    if (doMC2H_) {
        edm::FileInPath mc2hessianCSV("PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv");
        std::cout << "INFO: Will convert MC PDF set to Hessian with MC2Hessian\n";
        pdfweightshelper_.Init(N_LHEPDF_WEIGHTS_, N_MC2HESSIAN_WEIGHTS_, mc2hessianCSV);
    }
    centralWeightIndex_ = -1;
    // NNLOPSLike is just a config name for one MiNNLO sample
    if (name_.find("nnlops") != std::string::npos && name_.find("nnlopslike") == std::string::npos) {
        if (name_.find("nloOnly") == std::string::npos) {
            centralWeightIndex_ = 9;
            std::cout << "INFO: NNLOPS sample will be weighted by NNLO weight\n";
        }
        else
            std::cout << "INFO: Found NNLOPS sample but not applying weight\n";
    }
    else if (name_.find("ref_nnpdf31") != std::string::npos) {
        centralWeightIndex_ = 0;
        std::cout << "INFO: Sample will be weighted by 1st LHE weight\n";
    }

    doFiducial_ = selection_ != None;
    if (!doFiducial_)
        std::cout << "INFO: No fiducial selection will be applied\n";
    fReader.SetTree(tree);
}

void NanoGenSelectorBase::SetBranchesNanoAOD() {
    if (scaleWeights_) {
        b.SetSpecificBranch("nLHEScaleWeight", nLHEScaleWeight);
        b.SetSpecificBranch("LHEScaleWeight", LHEScaleWeight);
    }
    if (altScaleWeights_) {
        b.SetSpecificBranch("nLHEScaleWeightAltSet1", nLHEScaleWeightAltSet1);
        b.SetSpecificBranch("LHEScaleWeightAltSet1", LHEScaleWeightAltSet1);
    }
    for (size_t i = 0; i < pdfWeights_.size(); i++) {
        if (pdfWeights_.at(i)) {
            std::string name = "LHEPdfWeight";
            if (i > 0)
                name += "AltSet" + std::to_string(i);
            b.SetSpecificBranch("n"+name, nLHEPdfWeights.at(i));
            b.SetSpecificBranch(name, LHEPdfWeights[i]);
        }
        else {
            nLHEPdfWeights.at(i) = 0;
        }
    }
    if (unknownWeights_) {
        b.SetSpecificBranch("nLHEUnknownWeight", nLHEUnknownWeight);
        b.SetSpecificBranch("LHEUnknownWeight", LHEUnknownWeight);
        if (unknownWeightsAlt_) {
            b.SetSpecificBranch("nLHEUnknownWeightAltSet1", nLHEUnknownWeightAltSet1);
            b.SetSpecificBranch("LHEUnknownWeightAltSet1", LHEUnknownWeightAltSet1);
        }
    }
    if (paramWeights_) {
        b.SetSpecificBranch("nMEParamWeight", nMEParamWeight);
        b.SetSpecificBranch("MEParamWeight", MEParamWeight);
    }
    else if (paramWeightsUpd_) {
        b.SetSpecificBranch("nLHEReweightingWeightCorrectMass", nMEParamWeight);
        b.SetSpecificBranch("LHEReweightingWeightCorrectMass", MEParamWeight);
    }
}

void NanoGenSelectorBase::LoadBranchesNanoAOD(Long64_t entry, SystPair variation) { 
    weight = 1;
    if (scaleWeights_) {
        b.SetSpecificEntry(entry, "nLHEScaleWeight");
        b.SetSpecificEntry(entry, "LHEScaleWeight");
    }
    if (altScaleWeights_) {
        b.SetSpecificEntry(entry, "nLHEScaleWeightAltSet1");
        b.SetSpecificEntry(entry, "LHEScaleWeightAltSet1");
    }
    for (size_t i = 0; i < pdfWeights_.size(); i++) {
        if (pdfWeights_.at(i)) {
            std::string name = "LHEPdfWeight";
            if (i > 0)
                name += "AltSet" + std::to_string(i);
            b.SetSpecificEntry(entry, name);
            
            if (storeCenPdfs_ && (std::find(std::begin(cenPdfWeightsToStore_), std::end(cenPdfWeightsToStore_), i)
                    != std::end(cenPdfWeightsToStore_)))
                nLHEPdfWeights[i] = 1;
            else
                b.SetSpecificEntry(entry, "n"+name);
            // A Hack for the fact that CT18 and CT18Z go in the same set
            if (nLHEPdfWeights[i] == 122 && pdfCenWeight_ == 0) {
                nLHEPdfWeights[i] = 61;
            }
        }
    }
    if (unknownWeights_) {
        b.SetSpecificEntry(entry, "LHEUnknownWeight");
        b.SetSpecificEntry(entry, "nLHEUnknownWeight");
        if (unknownWeightsAlt_) {
            b.SetSpecificEntry(entry, "LHEUnknownWeightAltSet1");
            b.SetSpecificEntry(entry, "nLHEUnknownWeightAltSet1");
        }
    }
    if (paramWeights_) {
        b.SetSpecificEntry(entry, "MEParamWeight");
        b.SetSpecificEntry(entry, "nMEParamWeight");
    }
    else if (paramWeightsUpd_) {
        b.SetSpecificEntry(entry, "nLHEReweightingWeightCorrectMass");
        b.SetSpecificEntry(entry, "LHEReweightingWeightCorrectMass");
    }

    if (!scaleWeights_)
        nLHEScaleWeight = 0;
    if (!altScaleWeights_)
        nLHEScaleWeightAltSet1 = 0;
    if (!pdfWeights_.at(0))
        nLHEPdfWeight = {0};
    if (!unknownWeights_)
        nLHEUnknownWeight = 0;
    if (!unknownWeightsAlt_)
        nLHEUnknownWeightAltSet1 = 0;
    if (!(paramWeights_ || paramWeightsUpd_))
        nMEParamWeight = 0;
    fReader.SetLocalEntry(entry);

    channel_ = channelMap_[channelName_];

    // Sort descending
    auto compareMaxByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() > b.pt(); };
    std::vector<unsigned int> idsToKeep = {11, 12, 13, 14};
    if (doPhotons_)
        idsToKeep.push_back(22);

    // This is a bit nasty because it assumes that Central is the first variation, which is NOT guaranteed!
    if (variation.first == Central) {
        bornLeptons.clear();
        dressedLeptons.clear();
        preFSRLeptons.clear();
        bareLeptons.clear();
        lheLeptons.clear();

        bornNeutrinos.clear();
        preFSRNeutrinos.clear();
        fsneutrinos.clear();
        lheNeutrinos.clear();

        jets.clear();
        lhejets.clear();
        for (size_t i = 0; i < *nGenDressedLepton; i++) {
            LorentzVector vec;
            if (GenDressedLepton_hasTauAnc.At(i)) {
                continue;
            }
            dressedLeptons.emplace_back(makeGenParticle(GenDressedLepton_pdgId.At(i), 1, GenDressedLepton_pt.At(i), 
                                        GenDressedLepton_eta.At(i), GenDressedLepton_phi.At(i), GenDressedLepton_mass.At(i)));
        } // No need to sort, they're already pt sorted
        leptons = dressedLeptons;
        // only valid for photos. If there was a radiation, the 746 leptons will have the pre-FSR kinematics.
        // if no radiation, just use the final state leptons
        bool foundStatus746 = false;
        // Do once, only if bare or born leptons are requested
        if (nNeutrinos_ || doBorn_ || doBareLeptons_ || doPreFSR_) {
            for (size_t i = 0; i < *nGenPart; i++) {
                bool fromHardProcessFS = GenPart_status.At(i) == 1 && ((GenPart_statusFlags.At(i) >> 8) & 1);
                bool isHardProcess = (GenPart_statusFlags.At(i) >> 7) & 1;
                bool isPrompt = (GenPart_statusFlags.At(i) >> 0) & 1;
                if (!(isHardProcess || (GenPart_status.At(i) == 1 && isPrompt) || GenPart_status.At(i) == 746))
                    continue;
                if (std::find(idsToKeep.begin(), idsToKeep.end(), std::abs(GenPart_pdgId.At(i))) == idsToKeep.end())
                    continue;

                auto part = makeGenParticle(GenPart_pdgId.At(i), GenPart_status.At(i), GenPart_pt.At(i), 
                        GenPart_eta.At(i), GenPart_phi.At(i), GenPart_mass.At(i));
                if (std::abs(part.pdgId()) == 11 || std::abs(part.pdgId()) == 13) {
                    if (doBareLeptons_ && isPrompt && GenPart_status.At(i) == 1)
                        bareLeptons.emplace_back(part);
                    if (isHardProcess && doBorn_)
                        bornLeptons.emplace_back(part);
                    // Only works with photos!
                    if (doPreFSR_ && ((!foundStatus746 && fromHardProcessFS) || GenPart_status.At(i) == 746)) {
                        if (GenPart_status.At(i) == 746 && !foundStatus746) {
                            preFSRLeptons.clear();
                            foundStatus746 = true;
                        }
                        preFSRLeptons.emplace_back(part);
                    }
                }
                else if (std::abs(part.pdgId()) == 12 || std::abs(part.pdgId()) == 14) {
                    if (GenPart_status.At(i) == 1 && isPrompt)
                        fsneutrinos.emplace_back(part);
                    if (isHardProcess && doBorn_)
                        bornNeutrinos.emplace_back(part);
                    if (doPreFSR_ && ((!foundStatus746 && fromHardProcessFS) || GenPart_status.At(i) == 746)) {
                        if (GenPart_status.At(i) == 746 && !foundStatus746) {
                            preFSRLeptons.clear();
                            foundStatus746 = true;
                        }
                        preFSRNeutrinos.emplace_back(part);
                    }
                }
                else if (std::abs(part.pdgId()) == 22 && isPrompt) {
                    photons.emplace_back(part);
                }
            }
        }
        // Warning! Only really works for the W
        //We need the deltaR for Z case as well, check whether the selection of bare lepton can produce the same distribution , ng_ave < 8
        if (bareLeptons.size() > 0 && doPhotons_) {
            auto& lep = bareLeptons.at(0);
            photons.erase(std::remove_if(photons.begin(), photons.end(), 
                    [lep] (const reco::GenParticle& p) { //std::cout<<"INFO: DeltaR = "<<reco::deltaR(p, lep)<<"\n";
                                                         return reco::deltaR(p, lep) > 0.4; }),  photons.end());
           //for(unsigned int i=0;i<photons.size();++i){std::cout<<"INFO: photons remaining = "<<i<<"\n";}
        }
        if (doLHE_) {
            for (size_t i = 0; i < *nLHEPart; i++) {
                auto part = makeGenParticle(LHEPart_pdgId.At(i), 1, LHEPart_pt.At(i), 
                        LHEPart_eta.At(i), LHEPart_phi.At(i), LHEPart_mass.At(i));

                if (std::abs(part.pdgId()) == 11 || std::abs(part.pdgId()) == 13) {
                    lheLeptons.emplace_back(part);
                }
                else if (std::abs(part.pdgId()) == 12 || std::abs(part.pdgId()) == 14) {
                    lheNeutrinos.emplace_back(part);
                }
                // No dR check here, should add
                else if ((std::abs(part.pdgId()) < 6 || std::abs(part.pdgId()) == 21) && part.pt() > 30) {
                    lhejets.emplace_back(part.polarP4());
                }
            }
        }

        std::sort(lheLeptons.begin(), lheLeptons.end(), compareMaxByPt);
        std::sort(bareLeptons.begin(), bareLeptons.end(), compareMaxByPt);
        std::sort(fsneutrinos.begin(), fsneutrinos.end(), compareMaxByPt);
        std::sort(preFSRLeptons.begin(), preFSRLeptons.end(), compareMaxByPt);
        neutrinos = fsneutrinos;
    }
    else if (variation.first == BareLeptons) {
        leptons = bareLeptons;
        neutrinos = fsneutrinos;
    }
    else if (variation.first == BornParticles) {
        leptons = bornLeptons;
        neutrinos = bornNeutrinos;
    }
    else if (variation.first == PreFSRLeptons) {
        leptons = preFSRLeptons;
        neutrinos = preFSRNeutrinos;
    }
    else if (variation.first == LHEParticles) {
        leptons = lheLeptons;
        neutrinos = lheNeutrinos;
    }
    else if (variation.first == ptV0to3_lhe ||
            variation.first == ptV3to5_lhe || 
            variation.first == ptV5to7_lhe ||
            variation.first == ptV7to9_lhe ||
            variation.first == ptV9to12_lhe ||
            variation.first == ptV12to15_lhe || 
            variation.first == ptV15to20_lhe ||
            variation.first == ptV20to27_lhe ||
            variation.first == ptV27to40_lhe ||
            variation.first == ptV40toInf_lhe ) {
        leptons = lheLeptons;
        neutrinos = lheNeutrinos;
    }
    else {
        leptons = dressedLeptons;
        neutrinos = fsneutrinos;
    }
        
    if (variation.first != LHEParticles) {
        ht = 0;
        jets.clear();
        for (size_t i = 0; i < *nGenJet; i++) {
            auto jet = makeGenParticle(0, 1, GenJet_pt.At(i), 
                    GenJet_eta.At(i), GenJet_phi.At(i), GenJet_mass.At(i));
            // Just in case not using GenJetsNoNu
            if (jet.pt() > 30 && !helpers::overlapsCollection(jet.polarP4(), leptons, 0.4, nLeptons_) &&
                    !helpers::overlapsCollection(jet.polarP4(), neutrinos, 0.4, nNeutrinos_)) {
                ht += jet.pt();
                jets.emplace_back(jet.polarP4());
            }
        } // No need to sort jets, they're already pt sorted

    }
    else {
        jets.clear();
        for (auto& jet : lhejets) {
            if (!helpers::overlapsCollection(jet, leptons, 0.4, nLeptons_))
                jets.emplace_back(jet);
        }
    }

    genMet.SetPt(*MET_fiducialGenPt);
    genMet.SetPhi(*MET_fiducialGenPhi);
    genMet.SetM(0.);
    genMet.SetEta(0.);

    weight = *genWeight;       
    if (weightSignOnly_)
        weight = *genWeight > 0 ? 1 : -1;
    // These options don't really work together
    else if (weightSuppress_ && std::abs(*genWeight) > weightSuppress_) {
        weight = 0;
        *genWeight = 0;
    }

    pdfMaxStore_ = std::min(pdfMaxStore_, std::accumulate(nLHEPdfWeights.begin(), nLHEPdfWeights.end(), 0));

    if (refWeight == 1)
        refWeight = weight;

    if (centralWeightIndex_ != -1 && scaleWeights_) {
        rescaleWeight_ = LHEScaleWeight[centralWeightIndex_];
        weight *= LHEScaleWeight[centralWeightIndex_];
    }
    else if (altPdf_) {
        // Assuming the first set will always have the central values,
        // in general additional sets will be alpha_s variations
        int index = minnloPdfMap[pdfSet_].front();
        float pdfWeight = LHEPdfWeights[index][pdfCenWeight_];
        rescaleWeight_ = pdfWeight;
        weight *= pdfWeight;
    }
    if (doMC2H_)
        buildHessian2MCSet();

    SetComposite();
}

reco::GenParticle NanoGenSelectorBase::makeGenParticle(int pdgid, int status, float pt, float eta, float phi, float m) {
    LorentzVector vec;
    vec.SetPt(pt);
    vec.SetEta(eta);
    vec.SetPhi(phi);
    vec.SetM(m);
    int charge = (pdgid < 0) ? 1: -1;
    if (std::abs(pdgid) ==  12 || std::abs(pdgid) ==  14 || std::abs(pdgid) ==  16 || std::abs(pdgid) ==  22 || std::abs(pdgid) ==  23)
        charge = 0;
    
    auto lep = reco::GenParticle(charge, vec, reco::Particle::Point(), pdgid, status, true);
    return lep;
}

void NanoGenSelectorBase::SetScaleFactors() {
    std::string corrName;
    if (name_.find("minnlo__N3LLCorr") != std::string::npos)
        corrName = "scetlibCorr3D_";
    else if (name_.find("lowpu__N3LLCorr") != std::string::npos)
        corrName = "scetlibCorr3D_MG_";
    else if (name_.find("__MiNNLOCorr") != std::string::npos) {
        corrName = "MGtoMiNNLOCorr3D_";
        nScetlibWeights_ = 1;
    }
    else 
        return;

    corrName += isZ_ ? "Z" : (name_.find("wm") != std::string::npos ? "Wm" : "Wp");

    // Need to check Wm or Wp
    for (int i = 0; i < nScetlibWeights_; i ++) {
        auto objName = corrName+"_var"+std::to_string(i);
        auto* sf = GetInputList()->FindObject(objName.c_str());
        if (sf == nullptr)
            std::invalid_argument("Failed to read ScaleFactor " + objName);
        scetlibCorrs_.push_back(static_cast<ScaleFactor*>(sf));
    }

    if (scetlibCorrs_.at(0) == nullptr && name_.find("N3LLCorr") != std::string::npos) 
        std::invalid_argument("Must pass N3LL correction SF");
}

void NanoGenSelectorBase::buildHessian2MCSet() {
    double pdfWeights[N_LHEPDF_WEIGHTS_];
//    for (size_t i = 0; i < N_LHEPDF_WEIGHTS_; i++) {
//        pdfWeights[i] = LHEPdfWeight[i];
//    }
    pdfweightshelper_.DoMC2Hessian(1., const_cast<const double*>(pdfWeights), LHEHessianPdfWeight);
}

double NanoGenSelectorBase::breitWignerWeight(double offset) {

    double targetMass = MV_GEN_ + offset;
    //double gamma_cen = std::sqrt(MV_GEN_*MV_GEN_*(MV_GEN_*MV_GEN_+GAMMAV_GEN_*GAMMAV_GEN_));
    //double gamma = std::sqrt(targetMass*targetMass*(targetMass*targetMass+GAMMAV_GEN_*GAMMAV_GEN_));
    double s_hat = mVcorr*mVcorr*1000*1000;
    double offshell = s_hat - MV_GEN_*MV_GEN_;
    double offshellOffset = s_hat - targetMass*targetMass;
    double weight = (offshell*offshell + GAMMAV_GEN_*GAMMAV_GEN_*MV_GEN_*MV_GEN_)/
            (offshellOffset*offshellOffset + GAMMAV_GEN_*GAMMAV_GEN_*targetMass*targetMass);
    return weight;
}

void NanoGenSelectorBase::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();
    AddObject<TH1D>(mcWeights_, "genWeights", "gen weights", 200, -10, 10);
    AddObject<TH1D>(mcPdfWeights_, "MCPdfweights", "MC pdf weights", 200, 0, 2);
    AddObject<TH1D>(hesPdfWeights_, "Hesweights", "Hessian pdf weights", 200, 0, 2);
    AddObject<TH1D>(scaleWeightsHist_, "scaleweights", "Scale weights", 200, 0, 2);

    InitializeHistogramsFromConfig();
}
