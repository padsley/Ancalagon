// One energy deposit in one BGO crystal (SCNT, see BgoArray.cc), recorded
// by BgoSD::ProcessHits(). crystalID is the real detector-position number
// (1-30, see BgoArray.cc's ComputePlacements()) -- SCNT itself is a single
// shared logical volume placed once per position with copyNo=0 each time
// (only the outer HSNG placement carries the real per-position copy
// number), so BgoSD reads it off HSNG in the touchable history instead of
// off the hit volume's own copy number.
#pragma once

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class BgoHit : public G4VHit {
 public:
  BgoHit() = default;

  void* operator new(size_t);
  void operator delete(void*);

  G4int crystalID = 0;
  G4double edepMeV = 0.0;
  G4ThreeVector position;
  G4double timeNs = 0.0;
  G4int trackID = 0;
  G4String particleName;
};

using BgoHitsCollection = G4THitsCollection<BgoHit>;

extern G4ThreadLocal G4Allocator<BgoHit>* gBgoHitAllocator;

inline void* BgoHit::operator new(size_t) {
  if (!gBgoHitAllocator) gBgoHitAllocator = new G4Allocator<BgoHit>;
  return static_cast<void*>(gBgoHitAllocator->MallocSingle());
}

inline void BgoHit::operator delete(void* hit) {
  gBgoHitAllocator->FreeSingle(static_cast<BgoHit*>(hit));
}
