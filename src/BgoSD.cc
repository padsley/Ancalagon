#include "BgoSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"

BgoSD::BgoSD(const G4String& name) : G4VSensitiveDetector(name) {
  collectionName.insert("BgoHitsCollection");
}

void BgoSD::Initialize(G4HCofThisEvent* hce) {
  fHitsCollection = new BgoHitsCollection(SensitiveDetectorName, collectionName[0]);
  const G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool BgoSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  // Not filtered on edep>0 (unlike most Geant4 examples): this pilot's
  // physics list is transportation-only (see PhysicsList.cc/README), so
  // every step's energy deposit is exactly zero for now regardless of
  // material -- recording every step anyway still gives real hit
  // counts/positions/particle-ID, which is enough to confirm the
  // geometry/SD wiring itself is correct. edep becomes physically
  // meaningful once a real EM/hadronic physics list is added (a
  // separately tracked, not-yet-done item -- see README).
  auto* hit = new BgoHit();
  hit->edepMeV = step->GetTotalEnergyDeposit() / MeV;
  hit->position = step->GetPostStepPoint()->GetPosition();
  hit->timeNs = step->GetPostStepPoint()->GetGlobalTime() / ns;
  hit->trackID = step->GetTrack()->GetTrackID();
  hit->particleName = step->GetTrack()->GetDefinition()->GetParticleName();

  // SCNT (the hit volume, depth 0) is shared: SCNT -> MGOR (1) -> FNGR (2)
  // -> HSNG (3), and only HSNG's own placement carries the real
  // per-position copy number (1-30, see BgoArray.cc) -- see BgoHit.hh's
  // header comment.
  const G4TouchableHandle& touchable = step->GetPreStepPoint()->GetTouchableHandle();
  hit->crystalID = touchable->GetVolume(3)->GetCopyNo();

  fHitsCollection->insert(hit);
  return true;
}
