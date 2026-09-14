// Generic config-driven reaction-vertex kinematics -- see
// ReactionKinematics.hh for the algorithm/scope and ReactionConfig.hh for
// the file format. The physics specific to *this pilot's own* reaction
// (15O(alpha,gamma)19Ne -- masses, resonance energy, level scheme) lives
// in reactions/o15ag_19ne.reaction's own comments and README.md's
// "Reaction specification", not here.
#include "ReactionKinematics.hh"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

namespace {
constexpr double kAmuMeV = 931.49432;

// General relativistic boost of a 4-vector (eStar, pStar), given in a
// frame moving with velocity `betaVec` (as a fraction of c) relative to
// the lab, into the lab frame. Reduces to the along-Z-only formula used
// by an earlier, single-reaction version of this file when betaVec is
// along Z; needed in general here because only the *first* cascade
// step's parent (the compound resonance) is guaranteed to move along
// +z -- an intermediate excited recoil, after emitting its own first
// gamma, generally recoils in an arbitrary direction.
std::pair<double, G4ThreeVector> BoostToLab(double eStar, const G4ThreeVector& pStar,
                                             const G4ThreeVector& betaVec) {
  const double beta2 = betaVec.mag2();
  if (beta2 < 1.0e-24) return {eStar, pStar};  // parent effectively at rest in the lab
  const double gamma = 1.0 / std::sqrt(1.0 - beta2);
  const double betaDotP = betaVec.dot(pStar);
  const double eLab = gamma * (eStar + betaDotP);
  const G4ThreeVector pLab = pStar + betaVec * (((gamma - 1.0) * betaDotP / beta2) + gamma * eStar);
  return {eLab, pLab};
}
}  // namespace

ReactionKinematics::ReactionKinematics(const ReactionConfig& config)
    : fConfig(config),
      fBeamMassMeV(config.beam.A * kAmuMeV + config.beam.massExcessMeV),
      fTargetMassMeV(config.target.A * kAmuMeV + config.target.massExcessMeV),
      fRecoilGroundMassMeV(config.recoil.A * kAmuMeV + config.recoil.massExcessMeV) {}

double ReactionKinematics::LevelMassMeV(int level) const {
  if (level == 0) return fRecoilGroundMassMeV;
  const auto it = fConfig.levels.find(level);
  if (it == fConfig.levels.end()) {
    throw std::runtime_error("ReactionKinematics: a BRAT card references undefined level " +
                              std::to_string(level) + " (no matching LEVL card)");
  }
  return fRecoilGroundMassMeV + it->second.excitationMeV;
}

double ReactionKinematics::BeamKineticEnergyMeV() const {
  return BeamKineticEnergyMeVForEcm(fConfig.resonanceEnergyMeV);
}

double ReactionKinematics::BeamKineticEnergyMeVForEcm(double ecmAboveThresholdMeV) const {
  const double W = fBeamMassMeV + fTargetMassMeV + ecmAboveThresholdMeV;
  const double eBeam =
      (W * W - fBeamMassMeV * fBeamMassMeV - fTargetMassMeV * fTargetMassMeV) /
      (2.0 * fTargetMassMeV);
  return eBeam - fBeamMassMeV;
}

double ReactionKinematics::SampleBeamKineticEnergyMeV(double keMinMeV, double keMaxMeV) const {
  if (fConfig.resonanceWidthMeV <= 0.0) return BeamKineticEnergyMeV();

  // Direct (inverse-CDF) Breit-Wigner sampling, not naive accept/reject:
  // this reaction's own width can be many orders of magnitude narrower
  // than [keMinMeV, keMaxMeV] (e.g. ~50 eV against ~100s of keV for
  // 15O(alpha,gamma)19Ne's own 4.03 MeV resonance -- see the reaction
  // file's own RWID comment), so a uniform-then-reject sampler over that
  // whole range would need on the order of (range/width) attempts per
  // accepted sample -- here, thousands. tan() inverts the Lorentzian CDF
  // directly; the only rejection needed afterwards is for the rare
  // sample landing outside the physically achievable range at all (the
  // Lorentzian's own heavy tails, or a below-threshold Ecm).
  for (int attempt = 0; attempt < 10000; ++attempt) {
    const double u = G4UniformRand();
    const double ecm =
        fConfig.resonanceEnergyMeV + 0.5 * fConfig.resonanceWidthMeV * std::tan(CLHEP::pi * (u - 0.5));
    if (ecm <= 0.0) continue;  // below the beam+target threshold
    const double keMeV = BeamKineticEnergyMeVForEcm(ecm);
    if (keMeV >= keMinMeV && keMeV <= keMaxMeV) return keMeV;
  }
  throw std::runtime_error(
      "ReactionKinematics::SampleBeamKineticEnergyMeV: no sample landed in [keMinMeV, keMaxMeV] "
      "after 10000 attempts -- check that the target's entrance/exit beam energies actually "
      "straddle the resonance (ERES/RWID vs BKIN and the target's own thickness)");
}

ReactionKinematics::Event ReactionKinematics::GenerateEvent() const {
  return GenerateEvent(BeamKineticEnergyMeV());
}

ReactionKinematics::Event ReactionKinematics::GenerateEvent(double beamKineticEnergyMeV) const {
  const double eBeam = beamKineticEnergyMeV + fBeamMassMeV;
  const double pBeam = std::sqrt(eBeam * eBeam - fBeamMassMeV * fBeamMassMeV);
  // The compound resonance (level -1): formed from beam (along +z) +
  // at-rest target: invariant mass W follows from beam/target 4-momenta
  // directly (reduces to beamMass+targetMass+Ecm when beamKineticEnergyMeV
  // is exactly BeamKineticEnergyMeVForEcm(Ecm), as the two callers above
  // both arrange).
  const double W = std::sqrt(fBeamMassMeV * fBeamMassMeV + fTargetMassMeV * fTargetMassMeV +
                              2.0 * eBeam * fTargetMassMeV);

  double parentMassMeV = W;
  double eParentLab = eBeam + fTargetMassMeV;
  G4ThreeVector pParentLab(0.0, 0.0, pBeam);
  int currentLevel = -1;

  Event ev;
  while (true) {
    double totalPercent = 0.0;
    std::vector<const ReactionConfig::Branch*> options;
    for (const auto& b : fConfig.branches) {
      if (b.fromLevel == currentLevel) {
        options.push_back(&b);
        totalPercent += b.percent;
      }
    }
    if (options.empty()) {
      throw std::runtime_error("ReactionKinematics: no BRAT card has fromLevel=" +
                                std::to_string(currentLevel) +
                                " -- the cascade has nowhere to go from there");
    }

    double r = G4UniformRand() * totalPercent;
    const ReactionConfig::Branch* chosen = options.back();  // guards float rounding at the boundary
    for (const auto* b : options) {
      if (r <= b->percent) {
        chosen = b;
        break;
      }
      r -= b->percent;
    }

    const int daughterLevel = chosen->toLevel;
    const double daughterMassMeV = LevelMassMeV(daughterLevel);

    // Two-body breakup (parent -> daughter + one gamma), isotropic in the
    // parent's own rest frame (no angular correlation with any earlier
    // step in the same chain is modeled).
    const double eGammaStar =
        (parentMassMeV * parentMassMeV - daughterMassMeV * daughterMassMeV) / (2.0 * parentMassMeV);
    const double pStarMag = eGammaStar;
    const double eDaughterStar = parentMassMeV - eGammaStar;

    const double cosTheta = 2.0 * G4UniformRand() - 1.0;
    const double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    const double phi = 2.0 * CLHEP::pi * G4UniformRand();
    const G4ThreeVector pGammaStar(pStarMag * sinTheta * std::cos(phi),
                                    pStarMag * sinTheta * std::sin(phi), pStarMag * cosTheta);
    const G4ThreeVector pDaughterStar = -pGammaStar;

    const G4ThreeVector betaVec = pParentLab / eParentLab;
    const auto [eGammaLab, pGammaLab] = BoostToLab(eGammaStar, pGammaStar, betaVec);
    const auto [eDaughterLab, pDaughterLab] = BoostToLab(eDaughterStar, pDaughterStar, betaVec);
    (void)eGammaLab;  // == |pGammaLab| (massless); not needed separately.

    ev.gammaMomentaMeV.push_back(pGammaLab);

    if (daughterLevel == 0) {
      ev.recoilMomentumMeV = pDaughterLab;
      ev.recoilEnergyMeV = eDaughterLab;
      return ev;
    }

    // Cascade continues: this excited daughter is the next parent.
    currentLevel = daughterLevel;
    parentMassMeV = daughterMassMeV;
    eParentLab = eDaughterLab;
    pParentLab = pDaughterLab;
  }
}
