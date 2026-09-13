// Prints a one-line per-detector hit summary at the end of each event
// (hit count + total energy deposit) by pulling BgoSD/DsssdSD's own hits
// collections out of the event -- see those classes for why edep reads
// zero for now (no energy-loss physics registered yet).
#pragma once

#include "G4UserEventAction.hh"

class G4Event;

class EventAction : public G4UserEventAction {
 public:
  void EndOfEventAction(const G4Event* event) override;
};
