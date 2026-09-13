// Two-body-per-step reaction-vertex kinematics generator, driven entirely
// by a ReactionConfig (see ReactionConfig.hh) instead of being hardcoded
// for one specific reaction -- see reactions/o15ag_19ne.reaction (this
// pilot's own bundled config, 15O(alpha,gamma)19Ne) and README.md's
// "Reaction specification" for the physics and how that file's own
// numbers were derived/cross-checked.
#pragma once

#include <vector>

#include "G4ThreeVector.hh"
#include "ReactionConfig.hh"

class ReactionKinematics {
 public:
  explicit ReactionKinematics(const ReactionConfig& config);

  // Beam kinetic energy (lab frame, target at rest) this resonance needs
  // to be populated -- MeV, total (not per nucleon).
  double BeamKineticEnergyMeV() const;

  struct Event {
    G4ThreeVector recoilMomentumMeV;  // lab frame, MeV/c (ground-state recoil)
    double recoilEnergyMeV = 0.0;     // lab frame, total energy, MeV
    std::vector<G4ThreeVector> gammaMomentaMeV;  // lab frame, MeV/c each (|p|==energy);
                                                  // one per cascade step actually taken
  };

  // Samples one full decay chain, starting at the resonance (level -1) and
  // following ReactionConfig's own BRAT branches (each one weighted by its
  // own percent, renormalized against sibling branches from the same
  // level) down to the ground state (level 0), emitting one gamma per
  // step. Each step's breakup is isotropic in *that step's own parent's*
  // rest frame (no angular correlations between successive gammas are
  // modeled) and boosted to the lab frame with a general (not
  // axis-restricted) relativistic boost, since only the first step's
  // parent (the compound resonance) necessarily moves along +z -- an
  // intermediate excited recoil, after its own first gamma emission,
  // generally does not. Uses G4UniformRand.
  //
  // Scope: the whole chain is treated as instantaneous at the reaction
  // vertex (no separate tracking of an excited intermediate recoil's own,
  // real, but for these lifetimes minute, flight before each subsequent
  // decay); no angular distributions/correlations are modeled anywhere in
  // the chain.
  Event GenerateEvent() const;

 private:
  ReactionConfig fConfig;
  double fBeamMassMeV;
  double fTargetMassMeV;
  double fRecoilGroundMassMeV;

  double LevelMassMeV(int level) const;  // 0 -> ground; >0 -> ground+excitation
};
