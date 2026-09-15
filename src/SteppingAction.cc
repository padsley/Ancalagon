#include "SteppingAction.hh"

#include <cstdio>
#include <cstdlib>

#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

namespace {
// TRAJ_QUIET=1 suppresses the per-step "TRAJ ..." dump below (unset/0
// keeps it, the long-standing default every diagnostic mode used earlier
// relies on). Checked once (std::getenv on every one of potentially
// hundreds of millions of steps in a bulk production run would itself be
// wasteful) via a function-local static -- not a per-instance flag, since
// every SteppingAction in a single process shares the same environment.
bool TrajQuiet() {
  static const bool quiet = [] {
    const char* v = std::getenv("TRAJ_QUIET");
    return v && v[0] == '1';
  }();
  return quiet;
}
}  // namespace

void SteppingAction::UserSteppingAction(const G4Step* step) {
  G4Track* track = step->GetTrack();

  // G4's ion ionisation model (G4ionIonisation/G4BraggIonModel etc.)
  // computes a velocity-dependent "effective charge" for stopping-power
  // purposes (Ziegler-style partial charge-screening at low ion velocity)
  // and calls G4DynamicParticle::SetCharge() with it -- which permanently
  // overwrites GetCharge() for every *other* purpose too, including the
  // magnetic/electric field equation of motion (G4Mag_UsualEqRhs reads
  // the track's current charge fresh at every field evaluation). Left
  // alone, this ion arrives at D1 tracked as ~2.8+ instead of its real,
  // fixed 4+ charge state (confirmed: it drifts from the very first step,
  // in vacuum, before any material is even touched) -- silently wrong
  // rigidity, silently wrong bend angle, no exception or warning anywhere.
  // Force it back to the real charge state every step; a no-op for the
  // plain-proton standalone tests (fNominalChargeIonE defaults to 1, and
  // G4ionIonisation isn't registered for bare protons in the first place).
  if (track->GetParticleDefinition()->GetParticleType() == "nucleus") {
    // GetDynamicParticle() returns const -- the mutation it normally
    // guards against is exactly what a physics model like G4ionIonisation
    // itself does internally (see above), so this is the same class of
    // access, just from user code instead of a process.
    const_cast<G4DynamicParticle*>(track->GetDynamicParticle())->SetCharge(fNominalChargeIonE * eplus);
  }

  if (TrajQuiet()) return;

  const G4ThreeVector pos = track->GetPosition();
  const G4ThreeVector mom = track->GetMomentum();
  const G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetPhysicalVolume();
  const double edepMeV = step->GetTotalEnergyDeposit() / MeV;

  std::printf("TRAJ %10.4f %10.5f %10.5f %12.6f %12.6f %12.6f %-10s %10.6f q=%6.3f m=%12.4f stepLcm=%9.5f\n",
              pos.z() / cm, pos.x() / cm, pos.y() / cm, mom.x() / MeV, mom.y() / MeV, mom.z() / MeV,
              volume ? volume->GetName().c_str() : "OutOfWorld", edepMeV,
              track->GetDynamicParticle()->GetCharge() / eplus,
              track->GetDynamicParticle()->GetMass() / MeV, step->GetStepLength() / cm);
}
