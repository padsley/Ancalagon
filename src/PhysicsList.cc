#include "PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4Gamma.hh"
#include "G4IonConstructor.hh"
#include "G4Proton.hh"

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
PhysicsList::PhysicsList() : fEmPhysics(new G4EmStandardPhysics()) {}

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
}
