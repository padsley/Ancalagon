// Sensitive detector for the BGO array's SCNT (crystal) logical volume,
// see BgoArray.cc. Records one BgoHit per step with nonzero-or-not energy
// deposit (see .cc for why zero-edep steps aren't filtered out here, the
// way most Geant4 examples do).
#pragma once

#include "BgoHit.hh"
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BgoSD : public G4VSensitiveDetector {
 public:
  explicit BgoSD(const G4String& name);
  ~BgoSD() override = default;

  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

 private:
  BgoHitsCollection* fHitsCollection = nullptr;
};
