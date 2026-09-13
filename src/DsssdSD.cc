#include "DsssdSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

DsssdSD::DsssdSD(const G4String& name) : G4VSensitiveDetector(name) {
  collectionName.insert("DsssdHitsCollection");
}

void DsssdSD::Initialize(G4HCofThisEvent* hce) {
  fHitsCollection = new DsssdHitsCollection(SensitiveDetectorName, collectionName[0]);
  const G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool DsssdSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  // See BgoSD::ProcessHits's own comment: not filtered on edep>0, since
  // this pilot's physics list is transportation-only for now.
  auto* hit = new DsssdHit();
  hit->edepMeV = step->GetTotalEnergyDeposit() / MeV;
  hit->position = step->GetPostStepPoint()->GetPosition();
  hit->timeNs = step->GetPostStepPoint()->GetGlobalTime() / ns;
  hit->trackID = step->GetTrack()->GetTrackID();
  hit->particleName = step->GetTrack()->GetDefinition()->GetParticleName();

  fHitsCollection->insert(hit);
  return true;
}
