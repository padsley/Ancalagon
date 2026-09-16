#include "PhysicsList.hh"

#include <cstdlib>
#include <cstring>

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Gamma.hh"
#include "G4IonConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4Proton.hh"
#include "G4StepLimiter.hh"

// G4EmStandardPhysics is a G4VPhysicsConstructor -- ordinarily plugged
// into a G4VModularPhysicsList, but perfectly fine to own and delegate to
// directly from a hand-written G4VUserPhysicsList like this one (its own
// ConstructParticle()/ConstructProcess() are just called from ours). It
// constructs and wires up standard EM physics for every particle it
// itself constructs (gamma: photoelectric/Compton/pair/Rayleigh; e-/e+:
// ionisation/bremsstrahlung/multiple scattering; p/ions: ionisation/
// multiple scattering -- including G4GenericIon, so the recoil ion gets
// real energy loss too), which is exactly what BgoSD/DsssdSD (see
// README's "Sensitive detectors") need for their hits' edep to read
// anything but zero.
//
// EM_PHYSICS_LIST=option4 (diagnostic-only env var, defaults to the plain
// G4EmStandardPhysics above if unset) swaps in G4EmStandardPhysics_option4
// instead -- Geant4's own most precise standard EM constructor (finer
// energy-loss/fluctuation binning, more accurate low-energy ion stopping
// and multiple-scattering models; used e.g. by hadron-therapy
// applications where energy-loss straggling accuracy specifically
// matters). Same G4VPhysicsConstructor interface, so nothing else about
// PhysicsList needs to change to use it -- a way to cross-check whether
// this pilot's default choice's straggling is already right, not a
// permanent replacement for it.
namespace {
G4VPhysicsConstructor* MakeEmPhysics() {
  const char* choice = std::getenv("EM_PHYSICS_LIST");
  if (choice && std::strcmp(choice, "option4") == 0) {
    return new G4EmStandardPhysics_option4();
  }
  return new G4EmStandardPhysics();
}
}  // namespace

PhysicsList::PhysicsList() : fEmPhysics(MakeEmPhysics()) {}

PhysicsList::~PhysicsList() { delete fEmPhysics; }

void PhysicsList::ConstructParticle() {
  G4Proton::ProtonDefinition();
  G4Gamma::GammaDefinition();
  G4IonConstructor().ConstructParticle();  // registers G4GenericIon, needed
                                            // for G4IonTable::GetIon() (used
                                            // by --track-chain's recoil ion)
  fEmPhysics->ConstructParticle();  // e-/e+/mu/ions/etc. EM processes need
}

void PhysicsList::ConstructProcess() {
  AddTransportation();
  fEmPhysics->ConstructProcess();

  // G4UserLimits' own MaxAllowedStep is only enforced by a G4StepLimiter
  // process -- G4Transportation alone only limits step length by field
  // curvature/geometry-boundary criteria, which is why a smooth-field
  // container (E1/E2/a dipole/a quad, all one logical volume with no
  // internal geometric boundary) can take one single, large adaptive step
  // across its whole length even with G4UserLimits set on it. Registering
  // this (for every charged particle this pilot ever tracks -- proton,
  // e-/e+ from EM physics, and the generic ion species) is what actually
  // makes DetectorConstruction's FINE_STEP_CM diagnostic knob do anything.
  // A no-op for every particle unless FINE_STEP_CM is set (no G4UserLimits
  // attached to any volume otherwise), so this doesn't change default
  // behavior.
  auto* stepLimiter = new G4StepLimiter();
  auto* particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    if (particle->GetPDGCharge() != 0.0) {
      particle->GetProcessManager()->AddProcess(stepLimiter, -1, -1, 1);
    }
  }
}
