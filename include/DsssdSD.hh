#pragma once

#include "DsssdHit.hh"
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class DsssdSD : public G4VSensitiveDetector {
 public:
  explicit DsssdSD(const G4String& name);
  ~DsssdSD() override = default;

  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

 private:
  DsssdHitsCollection* fHitsCollection = nullptr;
};
