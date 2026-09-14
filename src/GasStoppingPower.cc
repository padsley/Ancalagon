#include "GasStoppingPower.hh"

#include <cmath>
#include <stdexcept>

#include "G4EmCalculator.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace {
// Fine enough that the table's own linear interpolation error is
// negligible next to the physical uncertainties already in play (the
// resonance width/entrance energy this pilot's own bundled reaction file
// uses are themselves order-of-magnitude/derived estimates, not measured
// beam-tune numbers -- see reactions/o15ag_19ne.reaction's own BKIN/RWID
// comments).
constexpr int kSubsteps = 2000;

// Standard kinematics: the maximum kinetic energy a projectile of mass M,
// kinetic energy keMeV can transfer to a free electron (mass m_e) in one
// elastic collision (see e.g. PDG's "Passage of particles through
// matter").
double TmaxMeV(double keMeV, double massMeV) {
  const double gamma = 1.0 + keMeV / massMeV;
  const double beta2 = 1.0 - 1.0 / (gamma * gamma);
  const double ratio = electron_mass_c2 / massMeV;
  return 2.0 * electron_mass_c2 * beta2 * gamma * gamma /
         (1.0 + 2.0 * gamma * ratio + ratio * ratio);
}
}  // namespace

GasStoppingPower::GasStoppingPower(int ionZ, int ionA, G4Material* material, double entranceKeMeV,
                                    double maxPathLengthCm) {
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(ionZ, ionA, 0.0);
  const double massMeV = ion->GetPDGMass() / MeV;
  const double zSquared = static_cast<double>(ionZ) * ionZ;
  const double electronDensityPerCm3 = material->GetElectronDensity() * cm3;

  G4EmCalculator calc;

  fTable.reserve(kSubsteps + 1);
  const double stepCm = maxPathLengthCm / kSubsteps;
  double depthCm = 0.0;
  double keMeV = entranceKeMeV;
  double omega2MeV2 = 0.0;
  fTable.push_back({depthCm, keMeV, omega2MeV2});
  for (int i = 0; i < kSubsteps; ++i) {
    if (keMeV <= 0.0) {
      fTable.push_back({maxPathLengthCm, 0.0, omega2MeV2});
      break;
    }
    // GetDEDX returns G4's own internal units; dividing by (MeV/cm) reads
    // it off in MeV per cm, the standard G4EmCalculator idiom.
    const double dEdxMeVPerCm = calc.GetDEDX(keMeV * MeV, ion, material) / (MeV / cm);

    // Bohr straggling variance per unit length (see this file's own
    // header comment): d(Omega^2)/dx = 2*pi*re^2*me*c^2*z^2*n_e*Tmax/beta^2.
    const double gamma = 1.0 + keMeV / massMeV;
    const double beta2 = 1.0 - 1.0 / (gamma * gamma);
    const double tmax = TmaxMeV(keMeV, massMeV);
    const double dOmega2dxMeV2PerCm = 2.0 * CLHEP::pi * classic_electr_radius * classic_electr_radius /
                                       (cm * cm) * electron_mass_c2 * zSquared *
                                       electronDensityPerCm3 * tmax / beta2;

    keMeV -= dEdxMeVPerCm * stepCm;
    if (keMeV < 0.0) keMeV = 0.0;
    omega2MeV2 += dOmega2dxMeV2PerCm * stepCm;
    depthCm += stepCm;
    fTable.push_back({depthCm, keMeV, omega2MeV2});
  }
}

double GasStoppingPower::DepthForKineticEnergyCm(double targetKeMeV) const {
  if (targetKeMeV > fTable.front().keMeV || targetKeMeV < fTable.back().keMeV) {
    throw std::runtime_error(
        "GasStoppingPower::DepthForKineticEnergyCm: targetKeMeV is outside "
        "[ExitKineticEnergyMeV(), EntranceKineticEnergyMeV()]");
  }
  for (std::size_t i = 1; i < fTable.size(); ++i) {
    if (fTable[i].keMeV <= targetKeMeV) {
      const Sample& lo = fTable[i - 1];  // higher KE, shallower
      const Sample& hi = fTable[i];      // lower KE, deeper
      if (lo.keMeV == hi.keMeV) return lo.depthCm;
      const double frac = (lo.keMeV - targetKeMeV) / (lo.keMeV - hi.keMeV);
      return lo.depthCm + frac * (hi.depthCm - lo.depthCm);
    }
  }
  return fTable.back().depthCm;  // targetKeMeV == ExitKineticEnergyMeV() exactly
}

double GasStoppingPower::StragglingSigmaMeV(double depthCm) const {
  if (depthCm <= fTable.front().depthCm) return std::sqrt(fTable.front().omega2MeV2);
  if (depthCm >= fTable.back().depthCm) return std::sqrt(fTable.back().omega2MeV2);
  for (std::size_t i = 1; i < fTable.size(); ++i) {
    if (fTable[i].depthCm >= depthCm) {
      const Sample& lo = fTable[i - 1];
      const Sample& hi = fTable[i];
      const double frac = (hi.depthCm == lo.depthCm)
                              ? 0.0
                              : (depthCm - lo.depthCm) / (hi.depthCm - lo.depthCm);
      return std::sqrt(lo.omega2MeV2 + frac * (hi.omega2MeV2 - lo.omega2MeV2));
    }
  }
  return std::sqrt(fTable.back().omega2MeV2);
}
