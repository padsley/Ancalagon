#include "SteppingAction.hh"

#include <cstdio>

#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

void SteppingAction::UserSteppingAction(const G4Step* step) {
  const G4Track* track = step->GetTrack();
  const G4ThreeVector pos = track->GetPosition();
  const G4ThreeVector mom = track->GetMomentum();
  const G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetPhysicalVolume();
  const double edepMeV = step->GetTotalEnergyDeposit() / MeV;

  std::printf("TRAJ %10.4f %10.5f %10.5f %12.6f %12.6f %12.6f %-10s %10.6f\n", pos.z() / cm,
              pos.x() / cm, pos.y() / cm, mom.x() / MeV, mom.y() / MeV, mom.z() / MeV,
              volume ? volume->GetName().c_str() : "OutOfWorld", edepMeV);
}
