#ifndef helpers_h
#define helpers_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

typedef reco::Particle::PolarLorentzVector LorentzVector;

namespace helpers {
template <class T>
bool overlapsCollection(const LorentzVector& cand, T collection,
                        const float deltaRCut, size_t maxCompare) {
    for (size_t i = 0; i < std::min(maxCompare, collection.size()); ++i) {
        if (reco::deltaR(collection[i], cand) < deltaRCut) {
            return true;
        }
    }
    return false;
}

// bool compareByPt (const reco::Candidate& a, const reco::Candidate& b) {
//    return a.pt() > b.pt();
//}
}  // namespace helpers

#endif
