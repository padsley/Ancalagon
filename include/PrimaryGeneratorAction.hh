#pragma once

#include <string>

#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "ReactionKinematics.hh"

// Launches a single particle parallel to z, offset in x/y, from well
// upstream of the first element -- equivalent to a GEANT3 gukine_mitray
// "ip,pkine" ray with theta=phi=0 and a transverse offset.
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  // Proton, at the given total momentum.
  PrimaryGeneratorAction(double x0Cm, double y0Cm, double momentumMeV, double z0Cm = -120.0);
  // An ion (Z, A) with the given net charge state (in units of e -- may
  // differ from Z for a partially-stripped ion) and total momentum. Used
  // for --track-chain: the magnetic and electric elements are tuned for
  // the same design rigidity/kinetic-energy-per-charge, which only a
  // matching ion mass -- not a proton -- can satisfy simultaneously (see
  // README).
  PrimaryGeneratorAction(double x0Cm, double y0Cm, double momentumMeV, double z0Cm, int ionZ,
                         int ionA, int ionChargeState);

  // Fires one real reaction event per call, sampled from `reactionFilePath`
  // (see ReactionConfig.hh for the file format; reactions/o15ag_19ne.reaction
  // is this pilot's own bundled 15O(alpha,gamma)19Ne config) from the
  // vertex (x0,y0,z0)cm: the recoil ion (ground state, at the config's own
  // charge state) plus one gamma per cascade step actually sampled (see
  // ReactionKinematics.hh), all correlated via real relativistic
  // kinematics. Used for --track-reaction. Throws std::runtime_error (from
  // ReactionConfig::Load) if the file can't be read or is malformed --
  // deliberately not caught here, since a broken reaction config should
  // stop the run, not silently fall back to something else.
  PrimaryGeneratorAction(const std::string& reactionFilePath, double x0Cm, double y0Cm,
                         double z0Cm);

  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

 private:
  void GenerateReactionPrimaries(G4Event* event);

  G4ParticleGun* fGun = nullptr;
  // Ion mode only: G4IonTable::GetIon() needs GenericIon to already be
  // constructed, which only happens once the run manager initializes the
  // physics list -- too late to call from this class's own constructor
  // (built before Initialize()), so it's deferred to the first
  // GeneratePrimaries() call instead.
  bool fNeedsIonLookup = false;
  int fIonZ = 0, fIonA = 0, fIonChargeState = 0;

  // Reaction mode only (fGun is unused: a variable number of correlated
  // primaries per event -- one recoil ion + one gamma per cascade step --
  // not a single-particle gun). fReaction owns the kinematics generator
  // built from the config file at construction time; fRecoilZ/A/ChargeState
  // come from that same config's own RECL card, not hardcoded.
  bool fIsReaction = false;
  ReactionKinematics* fReaction = nullptr;
  int fRecoilZ = 0, fRecoilA = 0, fRecoilChargeState = 0;
  double fReactionX0Cm = 0.0, fReactionY0Cm = 0.0, fReactionZ0Cm = 0.0;
};
