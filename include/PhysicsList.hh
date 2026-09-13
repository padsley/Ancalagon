#pragma once

#include "G4VUserPhysicsList.hh"

class G4VPhysicsConstructor;

// Transportation + standard EM physics (G4EmStandardPhysics, delegated
// to rather than reimplemented -- see .cc). The optics chain itself is
// still pure field tracking in vacuum, where EM physics is a no-op, but
// the BGO/DSSSD sensitive detectors (BgoSD/DsssdSD, see README's
// "Sensitive detectors") need real ionisation/EM processes registered
// for their hits' energy deposit to read anything but zero -- this is
// what adds that, without changing anything about the already-validated
// optics tracking.
class PhysicsList : public G4VUserPhysicsList {
 public:
  PhysicsList();
  ~PhysicsList() override;

  void ConstructParticle() override;
  void ConstructProcess() override;

 private:
  G4VPhysicsConstructor* fEmPhysics;
};
