#pragma once

#include "G4UserSteppingAction.hh"

// Dumps (z, x, y, px, py, pz, volume, step-edep) at every step so the
// trajectory through Q1 (or the whole chain) can be inspected/plotted --
// the GEANT4 analogue of GEANT3's TP7.dat ray trace log. Volume name and
// edep were added once real EM physics existed (see README's "Sensitive
// detectors"/"Full chain" caveat) specifically to localize *where* a
// track loses energy/scatters, not just that it does.
class SteppingAction : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* step) override;
};
