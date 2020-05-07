#include "Analysis/VVAnalysis/interface/FakeRateSelector.h"
#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include "Analysis/VVAnalysis/interface/ZSelector.h"
#include "Analysis/VVAnalysis/interface/LowPileupSelector.h"
#include "Analysis/VVAnalysis/interface/LowPileupZSelector.h"
#include "Analysis/VVAnalysis/interface/LowPileupWSelector.h"
#include "Analysis/VVAnalysis/interface/LowPileupWBackgroundSelector.h"
#include "Analysis/VVAnalysis/interface/ZZGenSelector.h"
#include "Analysis/VVAnalysis/interface/WGenSelector.h"
#include "Analysis/VVAnalysis/interface/ZGenSelector.h"
#include "Analysis/VVAnalysis/interface/NanoGenSelectorBase.h"
#include "Analysis/VVAnalysis/interface/WZSelector.h"
#include "Analysis/VVAnalysis/interface/ThreeLepSelector.h"
#include "Analysis/VVAnalysis/interface/WZSelectorBase.h"
#include "Analysis/VVAnalysis/interface/WZBackgroundSelector.h"
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/disambiguateFinalStates.h"
#include "Analysis/VVAnalysis/interface/disambiguateFinalStatesZZ.h"
#include "Analysis/VVAnalysis/interface/Efficiency.h"


namespace{
  namespace{
    FakeRateSelector pFakeRateSelector;
    WZSelectorBase pWZSelectorBase;
    SelectorBase pSelectorBase;
    ZSelector pZSelector;
    LowPileupZSelector pLowPileupZSelector;
    LowPileupWSelector pLowPileupWSelector;
    LowPileupWSelector pLowPileupWBackgroundSelector;
    LowPileupSelector pLowPileupSelector;
    WZSelector pWZSelector;
    NanoGenSelectorBase pNanoGenSelectorBase;
    ZZGenSelector pZZGenSelector;
    ThreeLepSelector pThreeLepSelector;
    WGenSelector pWGenSelector;
    ZGenSelector pZGenSelector;
    WZBackgroundSelector pWZBackgroundSelector;
    ScaleFactor pScaleFactor;
    disambiguateFinalStates pDisambiguator;
    disambiguateFinalStates pDisambiguatorZZ;
    Efficiency pEfficiency;
  }
}
