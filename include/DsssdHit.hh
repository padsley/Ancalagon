// One energy deposit in the DSSSD (a Micron W1, see DetectorConstruction.cc),
// recorded by DsssdSD::ProcessHits(). Position is kept in full (not yet
// mapped to the W1's real 16x16 junction/ohmic strip pitch -- the DSSSD is
// still a single unsegmented silicon slab, see README) so a future strip
// mapping can be added without changing this hit's own fields.
#pragma once

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class DsssdHit : public G4VHit {
 public:
  DsssdHit() = default;

  void* operator new(size_t);
  void operator delete(void*);

  G4double edepMeV = 0.0;
  G4ThreeVector position;
  G4double timeNs = 0.0;
  G4int trackID = 0;
  G4String particleName;
};

using DsssdHitsCollection = G4THitsCollection<DsssdHit>;

extern G4ThreadLocal G4Allocator<DsssdHit>* gDsssdHitAllocator;

inline void* DsssdHit::operator new(size_t) {
  if (!gDsssdHitAllocator) gDsssdHitAllocator = new G4Allocator<DsssdHit>;
  return static_cast<void*>(gDsssdHitAllocator->MallocSingle());
}

inline void DsssdHit::operator delete(void* hit) {
  gDsssdHitAllocator->FreeSingle(static_cast<DsssdHit*>(hit));
}
