// Beam-in-gas stopping power, via G4EmCalculator -- i.e. the *same*
// ionization physics already registered in this pilot's own PhysicsList,
// not an independently-sourced formula. Used by PrimaryGeneratorAction to
// find the actual depth into the target gas at which the beam's own
// continuously-degrading energy crosses the configured resonance (see
// ReactionConfig.hh's BKIN card), instead of firing every event from one
// fixed point at one fixed (energy-loss-free) beam energy.
//
// Only usable once the run manager has built its physics tables
// (G4EmCalculator queries those directly) -- i.e. after
// G4RunManager::Initialize(), not from a constructor that runs before it.
#pragma once

#include <vector>

class G4Material;
class G4ParticleDefinition;

class GasStoppingPower {
 public:
  // Builds and caches a kinetic-energy-vs-depth table for ion (Z, A) in
  // `material`, starting at entranceKeMeV and stepping through
  // maxPathLengthCm (or until the ion stops, if sooner). One G4EmCalculator
  // query per substep, all done here at construction -- not per later
  // lookup -- since PrimaryGeneratorAction calls this once (lazily, the
  // first time it's needed) and then reuses it for every event's own,
  // otherwise-cheap, depth lookup.
  GasStoppingPower(int ionZ, int ionA, G4Material* material, double entranceKeMeV,
                    double maxPathLengthCm);

  double EntranceKineticEnergyMeV() const { return fTable.front().keMeV; }
  // The kinetic energy at the far end of maxPathLengthCm (0 if the ion
  // stops before reaching it).
  double ExitKineticEnergyMeV() const { return fTable.back().keMeV; }

  // Depth (cm from the entrance) at which the cached table's kinetic
  // energy first drops to targetKeMeV, by linear interpolation between
  // the two bracketing cached samples. targetKeMeV must be within
  // [ExitKineticEnergyMeV(), EntranceKineticEnergyMeV()]; throws
  // std::runtime_error otherwise (a caller bug -- e.g. handing this a
  // Breit-Wigner sample it should itself have kept inside that range).
  double DepthForKineticEnergyCm(double targetKeMeV) const;

 private:
  struct Sample {
    double depthCm;
    double keMeV;
  };
  std::vector<Sample> fTable;
};
