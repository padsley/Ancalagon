// Books two ROOT ntuples (TTrees), "Bgo" and "Dsssd" (see EventAction.cc
// for what fills them and README's "Sensitive detectors" for the physics
// context), via Geant4's own G4AnalysisManager -- this writes a genuine
// ROOT file (dragon_hits.root, opened by a plain `root` session with no
// Geant4 dependency) using Geant4's embedded ROOT-format writer, not a
// dependency on an external ROOT install.
#pragma once

#include "G4UserRunAction.hh"

class G4Run;

class RunAction : public G4UserRunAction {
 public:
  RunAction();
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;
};
