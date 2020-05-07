
std::string boolStr(bool val) {
    if (val)
        return "true";
    else
        return "false";
}

void ThreeLepSelector::printInfo() {
    int num;

    if (!eventVec[event]) return;
    printf("Event %d, lumi %d\n", (int)event, (int)lumi);
    std::cout << info[event] << std::endl;
    printf(" #: %6s %5s  Id\n", "pt", "eta");
    printf("  : %6.4f %6.4f %6.4f\n", PV_x, PV_y, PV_z);

    num = 0;
    printf("-------leps-----------\n");
    printf("Muon: pt eta id iso tightcharge dz d0 sip ptratio ptrel\n");
    // printf("Elec: pt eta id mva cutbased iso dz d0 sip ptratio ptrel\n");
    // printf("Elec: sieie hoe 1/e-1/pt ecalSumEt hcalSumEt tkSumPt\n");
    for (auto lep : goodLeptons) {
        num++;
        LorentzVector closeJet(get4Vector(PID_JET, getCloseJetIndex(lep.v)));
        printf("%2d: %6.2f %+3.2f %+d", num, lep.Pt(), lep.Eta(),
               -1 * lep.Charge());
        if (lep.Id() == PID_MUON) {
            printf(" %4.2f %d %+5.3f %+5.3f %5.3f, %5.3f %6.3f",
                   Muon_miniPFRelIso_all[lep.index],
                   Muon_tightCharge[lep.index], Muon_dz[lep.index],
                   Muon_dxy[lep.index], Muon_sip3d[lep.index],
                   (Muon_pt[lep.index] / closeJet.Pt()),
                   LepRelPt(lep.v, closeJet));
        } else {
        }
        printf("\n");
    }
    num = 0;
    printf("-----loose-muon-------\n");
    for (auto lep : looseLeptons) {
        if (lep.Id() != PID_MUON) continue;
        int i = lep.index;
        LorentzVector v = lep.v;
        // for(int i = 0; i < (int)nMuon; i++) {
        num++;
        // LorentzVector v = get4Vector(PID_MUON,i);
        int closeIdx = getCloseJetIndex(v);
        LorentzVector closeJet = get4Vector(PID_JET, closeIdx);
        // closeJet += (1-Jet_L2L3[closeIdx])*v;
        closeJet =
            (Jet_L1[closeIdx] * (1 - Jet_rawFactor[closeIdx]) * closeJet - v) *
                Jet_L2L3[closeIdx] +
            v;

        // int i = lep.index;
        printf("%2d: %6.2f %+3.2f %+3.2f %+d", num, Muon_pt[i], Muon_eta[i],
               Muon_phi[i], -1 * Muon_charge[i]);
        printf(" %+5.3f %+5.3f %5.3f | %4.2f %5.3f %6.3f rawpt %6.3f",
               Muon_dz[i], Muon_dxy[i], Muon_sip3d[i], Muon_miniPFRelIso_all[i],
               v.Pt() / closeJet.Pt(), LepRelPt(v, closeJet),
               Jet_rawFactor[getCloseJetIndex(v)]);
        printf("\n");
        printf("    %6.2f %+3.2f %+3.2f | %6.2f %+3.2f %+3.2f | %7.5f\n",
               v.Pt(), v.Eta(), v.Phi(), closeJet.Pt(), closeJet.Eta(),
               closeJet.Phi(), reco::deltaR(closeJet, v));
    }
    // printf("--------all-muon-------\n");
    // for(size_t index = 0; index < nMuon; index++) {
    // 	printf("%2d: %6.2f %+3.2f %+d", (int)index, Muon_pt[index],
    // Muon_eta[index], (int)(-1*Muon_charge[index])); 	printf(" %4.2f %d %+5.3f
    // %+5.3f %5.3f %d", Muon_miniPFRelIso_all[index], Muon_tightCharge[index],
    // 	       Muon_dz[index], Muon_dxy[index], Muon_sip3d[index],
    // Muon_mediumId[index]); 	printf("\n");
    // 	       }

    num = 0;
    printf("-----loose-elec-------\n");
    for (auto lep : looseLeptons) {
        if (lep.Id() != PID_ELECTRON) continue;
        int i = lep.index;
        LorentzVector v = lep.v;
        // for(int i = 0; i < (int)nElectron; i++) {
        num++;
        int closeIdx = getCloseJetIndex(v);
        LorentzVector closeJet = get4Vector(PID_JET, closeIdx);
        // closeJet += (1-Jet_L2L3[closeIdx])*v;
        closeJet =
            (Jet_L1[closeIdx] * (1 - Jet_rawFactor[closeIdx]) * closeJet - v) *
                Jet_L2L3[closeIdx] +
            v;

        // LorentzVector v = get4Vector(PID_ELECTRON, i);
        // int i = lep.index;
        printf("  : %6s %5s %5s %2s %4s %2s %5s %6s %6s %5s %5s %6s\n", "pt",
               "eta", "phi", "ch", "mva", "ID", "iso", "dz", "d0", "sip3d",
               "ratio", "ptrel");
        printf("%2d: %6.2f %+3.2f %+3.2f %+d", num, Electron_pt[i] / CHGPT(i),
               Electron_eta[i], Electron_phi[i], -1 * Electron_charge[i]);
        printf(" %4.2f  %d %5.3f %+5.3f %+5.3f %5.3f %5.3f %6.3f %6.2f\n",
               Electron_MVA[i], Electron_cutBased[i],
               Electron_miniPFRelIso_all[i], Electron_dz[i], Electron_dxy[i],
               Electron_sip3d[i], v.Pt() / closeJet.Pt(), LepRelPt(v, closeJet),
               Electron_pt[i] / CHGPT(i) / Electron_eCorr[i]);
        int DEtaInCut = READFROM(Electron_vidBitmap[i], 6, 3);
        int DPhiInCut = READFROM(Electron_vidBitmap[i], 9, 3);

        printf(" %6s %6s %7s %6s %6s %6s, %5s, %5s\n", "sInIn", "HoverE",
               "1/E-1/p", "ECsum", "HCsum", "TKsum", "DnIn", "DphIn");

        printf(" %6.4f %6.4f %+6.4f %6.4f %6.4f %6.4f, %5d, %5d\n",
               Electron_sieie[i], Electron_hoe[i], Electron_eInvMinusPInv[i],
               Electron_dr03EcalRecHitSumEt[i],
               Electron_dr03HcalDepth1TowerSumEt[i], Electron_dr03TkSumPt[i],
               DEtaInCut, DPhiInCut

        );

        printf(
            " %6s %6s %7s %6s %6s %6s, %5s, %5s\n",
            boolStr(Electron_sieie[i] < 0.031).c_str(),
            boolStr(Electron_hoe[i] < 0.08).c_str(),
            boolStr(abs(Electron_eInvMinusPInv[i]) < 0.01).c_str(),
            boolStr(Electron_dr03EcalRecHitSumEt[i] / Electron_pt[i] *
                        CHGPT(i) <
                    0.45)
                .c_str(),
            boolStr(Electron_dr03HcalDepth1TowerSumEt[i] / Electron_pt[i] *
                        CHGPT(i) <
                    0.25)
                .c_str(),
            boolStr(Electron_dr03TkSumPt[i] / Electron_pt[i] * CHGPT(i) < 0.2)
                .c_str(),
            boolStr(DEtaInCut >= 1).c_str(), boolStr(DPhiInCut >= 3).c_str());
        printf(" %5s %5s %5s %5s %5s %5s %5s %5s %5s\n", "pt", "convV", "missV",
               "tChar", "dz", "d0", "sip3d", "mva", "iso");
        printf(" %5s %5s %5s %5s %5s %5s %5s %5s %5s\n",
               boolStr(Electron_pt[i] / CHGPT(i) > 20).c_str(),
               boolStr(Electron_convVeto[i]).c_str(),
               boolStr(Electron_lostHits[i] == 0).c_str(),
               boolStr(Electron_tightCharge[i] == 2).c_str(),
               boolStr(abs(Electron_dz[i]) < 0.1).c_str(),
               boolStr(abs(Electron_dxy[i]) < 0.05).c_str(),
               boolStr(Electron_sip3d[i] < 4).c_str(),
               boolStr(passMVACut(mvaValues[4], i)).c_str(),
               boolStr(Electron_miniPFRelIso_all[i] < 0.12).c_str());

        printf("\n");
    }

    num = 0;
    printf("-------Jets-----------\n");
    for (auto lep : goodJets) {
        if (lep.passedJetSel()) {
            num++;
            printf("%2d: %6.2f %+4.2f %+4.2f %5.3f %5.3f\n", num, lep.Pt(),
                   lep.Eta(), lep.Phi(), Jet_btagDeepB[lep.index],
                   Jet_btagCSVV2[lep.index]);
        }
    }
    // for(size_t i = 0; i < nJet; i++) {
    // 	if(Jet_btagDeepB[i] < 0.6) continue;
    // 	printf("%2d: %6.2f %+4.2f %+4.2f %6.4f %6.4f %d\n", (int)i, Jet_pt[i],
    // Jet_eta[i], Jet_phi[i], Jet_btagDeepB[i],  Jet_btagCSVV2[i],
    // Jet_jetId[i]);
    // }
    num = 0;
    printf("-------BJets------------\n");
    for (auto lep : goodJets) {
        if (lep.passedBJetSel()) {
            num++;
            printf("%2d: %6.2f %+4.2f %+4.2f %5.3f %5.3f\n", num, lep.Pt(),
                   lep.Eta(), lep.Phi(), Jet_btagDeepB[lep.index],
                   Jet_btagCSVV2[lep.index]);
        }
    }
    printf("-------MET------------\n");
    printf("  : %6.2f \n", MET);

    printf("----------------------\n\n");
}
