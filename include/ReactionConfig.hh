// Plain-text reaction specification, loaded at run time instead of
// hardcoded (see README's "Reaction specification" for the physics this
// pilot's own bundled file, reactions/o15ag_19ne.reaction, encodes, and
// why). Deliberately styled after dragon_2003.ffcards (this project's own
// existing config-file convention): unquoted 4-letter card keyword, then
// that card's own values, one card per line; 'COMM' lines and blank lines
// are ignored; parsing stops at 'SENT' or end of file.
//
// Cards:
//   BEAM  Z  A  massExcessMeV                  -- beam nucleus, ground state
//   TARG  Z  A  massExcessMeV                  -- target nucleus, at rest, ground state
//   RECL  Z  A  massExcessMeV  chargeState      -- recoil nucleus, ground state;
//                                                  chargeState is its production
//                                                  charge (units of e) this
//                                                  separator's fields are tuned for
//   ERES  resonanceEnergyMeV                    -- CM energy above the beam+target
//                                                  threshold at which the
//                                                  compound state forms
//   LEVL  level  excitationMeV  lifetimeSeconds -- one bound excited state of the
//                                                  recoil nucleus below the
//                                                  resonance; `level` is a
//                                                  positive integer, assigned
//                                                  by whoever writes the file
//                                                  (need not be in energy order)
//   BRAT  fromLevel  percent  toLevel           -- one branch of a decay: from
//                                                  `fromLevel` (0 is never a
//                                                  valid fromLevel -- the
//                                                  ground state is stable; -1
//                                                  means the resonance itself,
//                                                  i.e. this branch is one of
//                                                  the resonance's own decay
//                                                  modes) to `toLevel` (0 means
//                                                  the ground state, i.e. this
//                                                  branch ends the cascade),
//                                                  with relative weight
//                                                  `percent` (branches from the
//                                                  same fromLevel are
//                                                  normalized by their own
//                                                  sum, so they don't strictly
//                                                  need to add to 100)
//   SENT                                         -- end of file (optional --
//                                                  parsing also just stops at
//                                                  EOF)
//
// Masses all use the AME mass-excess convention throughout (mass = A*amu +
// massExcessMeV, amu = 931.49432 MeV/c^2), including the target -- unlike
// src/ureact.f's own hemass constant (a literal atomic mass), so every
// particle in this file is specified the same way. Every excited level
// must be reachable from level -1 (the resonance) via BRAT cards, directly
// or through other levels, and every chain must eventually reach level 0
// (the ground state) -- ReactionConfig::Load() does not itself check this;
// ReactionKinematics::GenerateEvent() will loop forever (or crash) on a
// dangling/cyclic level with no path to 0.
#pragma once

#include <map>
#include <string>
#include <vector>

struct ReactionConfig {
  struct ParticleSpec {
    int Z = 0;
    int A = 0;
    double massExcessMeV = 0.0;
  };

  struct Level {
    double excitationMeV = 0.0;
    double lifetimeS = 0.0;
  };

  struct Branch {
    int fromLevel = 0;  // -1 = the resonance
    double percent = 0.0;
    int toLevel = 0;  // 0 = ground state
  };

  ParticleSpec beam;
  ParticleSpec target;
  ParticleSpec recoil;
  int recoilChargeState = 0;
  double resonanceEnergyMeV = 0.0;
  std::map<int, Level> levels;  // keyed by level number, 0/-1 not stored here
  std::vector<Branch> branches;

  // Throws std::runtime_error (with a message naming the file and what's
  // wrong) on a missing file or a malformed/incomplete card.
  static ReactionConfig Load(const std::string& path);
};
