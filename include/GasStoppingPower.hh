// Beam-in-gas stopping power, via G4EmCalculator -- i.e. the *same*
// ionization physics already registered in this pilot's own PhysicsList,
// not an independently-sourced formula. Used by PrimaryGeneratorAction to
// find the actual depth into the target gas at which the beam's own
// continuously-degrading energy crosses the configured resonance (see
// ReactionConfig.hh's BKIN card), instead of firing every event from one
// fixed point at one fixed (energy-loss-free) beam energy.
//
// Also accumulates Bohr's own energy-loss straggling variance alongside
// the mean energy loss (see StragglingSigmaMeV()): the mean dE/dx curve
// alone says every event of the same species reaches the same depth at
// the same energy, which isn't physical -- real energy loss is a
// statistical process (a finite number of discrete collisions), so real
// ions of the same species/entrance energy actually spread out around
// that mean curve. This is the classical, non-relativistic-collision-
// cross-section-derived Gaussian approximation (Bohr 1915; see e.g. the
// PDG's "Passage of particles through matter" review), valid when the
// number of collisions is large (true here: a dilute, several-cm gas
// path) -- not the more general (and here unnecessary) Landau/Vavilov
// treatment for thin absorbers/few collisions. One simplification against
// full self-consistency with G4's own registered ionisation model: the
// projectile charge in the Bohr formula is taken as the bare Z (fully
// ionized), not G4's own internally-modeled, velocity-dependent effective
// charge (which the mean dE/dx from G4EmCalculator does already account
// for) -- flagged here, not silently assumed exact.
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
  // Builds and caches kinetic-energy- and straggling-variance-vs-depth
  // tables for ion (Z, A) in `material`, starting at entranceKeMeV and
  // stepping through maxPathLengthCm (or until the ion stops, if sooner).
  // One G4EmCalculator query per substep, all done here at construction --
  // not per later lookup -- since PrimaryGeneratorAction calls this once
  // (lazily, the first time it's needed) and then reuses it for every
  // event's own, otherwise-cheap, depth/straggling lookups.
  GasStoppingPower(int ionZ, int ionA, G4Material* material, double entranceKeMeV,
                    double maxPathLengthCm);

  double EntranceKineticEnergyMeV() const { return fTable.front().keMeV; }
  // The kinetic energy at the far end of maxPathLengthCm (0 if the ion
  // stops before reaching it).
  double ExitKineticEnergyMeV() const { return fTable.back().keMeV; }

  // Depth (cm from the entrance) at which the cached table's mean
  // kinetic energy first drops to targetKeMeV, by linear interpolation
  // between the two bracketing cached samples. targetKeMeV must be
  // within [ExitKineticEnergyMeV(), EntranceKineticEnergyMeV()]; throws
  // std::runtime_error otherwise (a caller bug -- e.g. handing this a
  // Breit-Wigner sample it should itself have kept inside that range).
  double DepthForKineticEnergyCm(double targetKeMeV) const;

  // Bohr straggling standard deviation (MeV) accumulated from the
  // entrance to depthCm -- i.e. the real, per-ion spread around the mean
  // kinetic energy DepthForKineticEnergyCm()'s own table assumes at that
  // depth. depthCm is clamped to [0, the table's own max depth].
  double StragglingSigmaMeV(double depthCm) const;

 private:
  struct Sample {
    double depthCm;
    double keMeV;
    double omega2MeV2;  // accumulated Bohr straggling variance, entrance to here
  };
  std::vector<Sample> fTable;
};
