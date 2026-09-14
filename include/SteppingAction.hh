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
  // nominalChargeIonE: the charge state (in units of e) the tracked ion is
  // supposed to keep for the whole chain (e.g. 4 for this pilot's 19Ne4+
  // recoil) -- see UserSteppingAction's own comment for why this has to be
  // force-reapplied every step. Irrelevant (default 1) for the plain-
  // proton standalone element tests, which don't hit this bug.
  explicit SteppingAction(double nominalChargeIonE = 1.0) : fNominalChargeIonE(nominalChargeIonE) {}

  void UserSteppingAction(const G4Step* step) override;

 private:
  double fNominalChargeIonE;
};
