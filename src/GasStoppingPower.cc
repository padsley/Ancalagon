#include "GasStoppingPower.hh"

#include <stdexcept>

#include "G4EmCalculator.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

namespace {
// Fine enough that the table's own linear interpolation error is
// negligible next to the physical uncertainties already in play (the
// resonance width/entrance energy this pilot's own bundled reaction file
// uses are themselves order-of-magnitude/derived estimates, not measured
// beam-tune numbers -- see reactions/o15ag_19ne.reaction's own BKIN/RWID
// comments).
constexpr int kSubsteps = 2000;
}  // namespace

GasStoppingPower::GasStoppingPower(int ionZ, int ionA, G4Material* material, double entranceKeMeV,
                                    double maxPathLengthCm) {
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(ionZ, ionA, 0.0);
  G4EmCalculator calc;

  fTable.reserve(kSubsteps + 1);
  const double stepCm = maxPathLengthCm / kSubsteps;
  double depthCm = 0.0;
  double keMeV = entranceKeMeV;
  fTable.push_back({depthCm, keMeV});
  for (int i = 0; i < kSubsteps; ++i) {
    if (keMeV <= 0.0) {
      fTable.push_back({maxPathLengthCm, 0.0});
      break;
    }
    // GetDEDX returns G4's own internal units; dividing by (MeV/cm) reads
    // it off in MeV per cm, the standard G4EmCalculator idiom.
    const double dEdxMeVPerCm = calc.GetDEDX(keMeV * MeV, ion, material) / (MeV / cm);
    keMeV -= dEdxMeVPerCm * stepCm;
    if (keMeV < 0.0) keMeV = 0.0;
    depthCm += stepCm;
    fTable.push_back({depthCm, keMeV});
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
