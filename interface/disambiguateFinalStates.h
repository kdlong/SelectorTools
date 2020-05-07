/*
 * TSelector deduplicate events in combinatorical ntuple
 * by choosing event with l1_l2_Mass pairing closest to
 * M_PDG(Z) = 91.1876
 *
 * Modified from N. Smith, U. Wisconsin
 *
 * https://github.com/nsmith-/ZHinvAnalysis/blob/master/disambiguateFinalStates.h
 */
#ifndef disambiguateFinalStates_h
#define disambiguateFinalStates_h

#include <TChain.h>
#include <TEntryList.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeFormula.h>

class disambiguateFinalStates : public TSelector {
   public:
    TTree *fChain;

    Float_t Mass;
    ULong64_t evt;
    UInt_t run;

    TBranch *b_Mass;
    TBranch *b_evt;
    TBranch *b_run;

    TEntryList *fBestCandidateEntryList;

    // TODO Add TTreeFormula, fix issue when using chain.Process
    // over multiple files (currently segfaults)
    TTreeFormula *fCutFormula;

    disambiguateFinalStates(TTree * /*tree*/ = 0)
        : fChain(0),
          fBestCandidateEntryList(0),
          fCutFormula(0),
          fCurrentRun(-1),
          fCurrentEvt(-1) {}
    virtual ~disambiguateFinalStates() {
        SafeDelete(fBestCandidateEntryList);
        SafeDelete(fCutFormula);
    }
    virtual Int_t Version() const { return 2; }
    virtual void Begin(TTree *tree);
    virtual void SlaveBegin(TTree *tree);
    virtual void Init(TTree *tree);
    virtual Bool_t Notify();
    virtual Bool_t Process(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
        return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
    }
    virtual void SetOption(const char *option) { fOption = option; }
    virtual void SetObject(TObject *obj) { fObject = obj; }
    virtual void SetInputList(TList *input) { fInput = input; }
    virtual TList *GetOutputList() const { return fOutput; }
    virtual void SlaveTerminate();
    virtual void Terminate();
    void setZCandidateBranchName(const char *name) {
        zCand_name = const_cast<char *>(name);
    }

    ClassDef(disambiguateFinalStates, 0);

   private:
    void findBestEntry();
    UInt_t fCurrentRun;
    ULong64_t fCurrentEvt;
    std::vector<Long64_t> fEntriesToCompare;
    std::vector<Float_t> fEntryDiscriminants;
    char *zCand_name;
};

#endif
