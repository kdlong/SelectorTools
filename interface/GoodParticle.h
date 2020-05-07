#ifndef _GOOD_PARTICLE_H_
#define _GOOD_PARTICLE_H_

#include "helpers.h"

// Wrapper for lepton objects so to ignore problems with
// accessing data if particle is muon vs electron

struct GoodPart {
    enum JET_STATUS { NONE = 0, PASS_JET = 1, PASS_BJET = 2, PASS_BOTH = 3 };
    LorentzVector v;
    int pdgId;
    int jetStatus = NONE;
    int index;
    GoodPart() {}

    GoodPart(double pt, double eta, double phi, double m, int pdg = 0)
        : v(pt, eta, phi, m), pdgId(pdg) {}
    GoodPart(LorentzVector v_, int pdg = 0) : v(v_), pdgId(pdg) {}
    double operator[](int index) {
        if (index == 1)
            return v.Px();
        else if (index == 2)
            return v.Py();
        else if (index == 3)
            return v.Pz();
        else
            return -1;
    }

    void SetPdgId(int pdg) { pdgId = pdg; }

    double Pt() { return v.Pt(); }
    double Eta() { return v.Eta(); }
    double Phi() { return v.Phi(); }
    double M() { return v.M(); }
    double E() { return v.E(); }

    int Charge() { return pdgId; }
    int Id() { return std::abs(pdgId); }

    void setPassJetSel() { jetStatus += PASS_JET; }
    void setPassBJetSel() { jetStatus += PASS_BJET; }

    bool passedJetSel() { return jetStatus % 2 == 1; }
    bool passedBJetSel() { return jetStatus / 2 == 1; }
};

// Seems a little less clunky
// Can also add some help functions as time goes on
struct GenPart {
    enum Match { NONE = 0, GEN_ONLY = 1, RECO_ONLY = 2, MATCHED = 3 };
    enum BMatch { bNONE = 0, bGEN_ONLY = 1, bRECO_ONLY = 2, bMATCHED = 3 };

    GoodPart gen;
    GoodPart reco;
    int status = NONE;
    int bstatus = bNONE;

    void SetupGen(double pt, double eta, double phi, double m, int pdg) {
        gen = GoodPart(pt, eta, phi, m, pdg);
        status += GEN_ONLY;
        if (pdg == 5) bstatus += bGEN_ONLY;
    }

    void SetupReco(GoodPart& reco_) {
        reco = reco_;
        status += RECO_ONLY;
        if (reco.passedBJetSel()) bstatus += bRECO_ONLY;
    }

    int Id() {
        if (status % 2 == 1)
            return gId();
        else
            return rId();
    }

    int gId() { return gen.Id(); }
    double gPt() { return gen.Pt(); }
    double gEta() { return gen.Eta(); }
    double gPhi() { return gen.Phi(); }
    double gM() { return gen.M(); }
    LorentzVector gVector() { return gen.v; }

    int rId() { return reco.Id(); }
    double rPt() { return reco.Pt(); }
    double rEta() { return reco.Eta(); }
    double rPhi() { return reco.Phi(); }
    double rM() { return reco.M(); }
    LorentzVector rVector() { return reco.v; }

    bool isMatched() { return status == MATCHED; }
    bool isFaked() { return status == RECO_ONLY; }
    bool noMatched() { return status == GEN_ONLY; }

    bool isBMatched() { return bstatus == bMATCHED; }
    bool isBFaked() { return bstatus == bRECO_ONLY; }
    bool noBMatched() { return bstatus == bGEN_ONLY; }
    bool isGenB() { return bstatus % 2 == 1; }
};

#endif
