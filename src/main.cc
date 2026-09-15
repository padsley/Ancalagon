#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "G4DynamicParticle.hh"
#include "G4EmCalculator.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VisExecutive.hh"
#include "G4VPhysicalVolume.hh"
#include "MitrayDipoleField.hh"
#include "MitrayEdipoleField.hh"
#include "MitrayQuadrupoleField.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "ReactionConfig.hh"
#include "ReactionKinematics.hh"
#include "SteppingAction.hh"
#include "TargetChamber.hh"

namespace {

// Dumps B(x,y,z) for a quadrupole on a grid, in the element-local
// A-coordinate system, so the CSV can be diffed against the matching
// standalone Fortran probe. No G4RunManager needed -- this exercises only
// the field port.
int RunPoleFieldProbe(const MitrayPoleData& data, const char* outFile, double xlo, double xhi,
                       double ylo, double yhi, double zlo, double zhi) {
  const G4ThreeVector centerCm(0.0, 0.0, 0.0);
  MitrayQuadrupoleField field(data, centerCm);

  std::FILE* f = std::fopen(outFile, "w");
  std::fprintf(f, "x_cm,y_cm,z_cm,bx_T,by_T,bz_T\n");

  const int nx = 5, ny = 5, nz = 21;

  for (int iz = 0; iz < nz; ++iz) {
    const double z = zlo + (zhi - zlo) * iz / (nz - 1);
    for (int ix = 0; ix < nx; ++ix) {
      const double x = xlo + (xhi - xlo) * ix / (nx - 1);
      for (int iy = 0; iy < ny; ++iy) {
        const double y = ylo + (yhi - ylo) * iy / (ny - 1);

        double bfld[3];
        field.FieldInLocalCm(x, y, z, bfld);

        std::fprintf(f, "%10.4f,%10.4f,%10.4f,%16.8E,%16.8E,%16.8E\n", x, y, z, bfld[0], bfld[1],
                     bfld[2]);
      }
    }
  }

  std::fclose(f);
  std::printf("Wrote %s\n", outFile);
  return 0;
}

int RunQ1FieldProbe() {
  return RunPoleFieldProbe(MitrayPoleData::Q1(), "q1_cpp_reference.csv", -3.0, 3.0, -3.0, 3.0,
                            -25.0, 25.0);
}

int RunQ2FieldProbe() {
  return RunPoleFieldProbe(MitrayPoleData::Q2(), "q2_cpp_reference.csv", -4.5, 4.5, -4.5, 4.5,
                            -35.0, 35.0);
}

// Q3 through Q14 combined, grid matching
// validate/mitray_poles_probe_rest.f exactly (including its per-element
// grid extent, derived from each quad's own RAD/L).
int RunQ3ThroughQ14FieldProbe() {
  struct Entry {
    const char* name;
    MitrayPoleData data;
  };
  const Entry entries[] = {
      {"Q3  ", MitrayPoleData::Q3()},   {"Q4  ", MitrayPoleData::Q4()},
      {"Q5  ", MitrayPoleData::Q5()},   {"Q6  ", MitrayPoleData::Q6()},
      {"Q7  ", MitrayPoleData::Q7()},   {"Q8  ", MitrayPoleData::Q8()},
      {"Q9  ", MitrayPoleData::Q9()},   {"Q10 ", MitrayPoleData::Q10()},
      {"Q11 ", MitrayPoleData::Q11()},  {"Q12 ", MitrayPoleData::Q12()},
      {"Q13 ", MitrayPoleData::Q13()},  {"Q14 ", MitrayPoleData::Q14()},
  };

  std::FILE* f = std::fopen("q3_q14_cpp_reference.csv", "w");
  std::fprintf(f, "name,x_cm,y_cm,z_cm,bx_T,by_T,bz_T\n");

  const int nx = 5, ny = 5, nz = 11;
  const G4ThreeVector centerCm(0.0, 0.0, 0.0);

  for (const auto& e : entries) {
    MitrayQuadrupoleField field(e.data, centerCm);
    const double xlo = -0.5 * e.data.RAD, xhi = 0.5 * e.data.RAD;
    const double ylo = xlo, yhi = xhi;
    const double zlo = -1.3 * e.data.L / 2.0 - 5.0, zhi = 1.3 * e.data.L / 2.0 + 5.0;

    for (int iz = 0; iz < nz; ++iz) {
      const double z = zlo + (zhi - zlo) * iz / (nz - 1);
      for (int ix = 0; ix < nx; ++ix) {
        const double x = xlo + (xhi - xlo) * ix / (nx - 1);
        for (int iy = 0; iy < ny; ++iy) {
          const double y = ylo + (yhi - ylo) * iy / (ny - 1);

          double bfld[3];
          field.FieldInLocalCm(x, y, z, bfld);

          std::fprintf(f, "%s,%10.4f,%10.4f,%10.4f,%16.8E,%16.8E,%16.8E\n", e.name, x, y, z,
                       bfld[0], bfld[1], bfld[2]);
        }
      }
    }
  }

  std::fclose(f);
  std::printf("Wrote q3_q14_cpp_reference.csv\n");
  return 0;
}

// Same idea as RunPoleFieldProbe(), for the D1 dipole -- grid matches
// validate/mitray_dipole_probe.f exactly.
int RunDipoleFieldProbe() {
  const MitrayDipoleData d1data = MitrayDipoleData::D1();
  const G4ThreeVector centerCm(0.0, 0.0, 0.0);
  MitrayDipoleField field(d1data, centerCm);

  std::FILE* f = std::fopen("d1_cpp_reference.csv", "w");
  std::fprintf(f, "x_cm,y_cm,z_cm,bx_T,by_T,bz_T\n");

  const int nx = 7, ny = 5, nz = 25;
  const double xlo = -30.0, xhi = 30.0;
  const double ylo = -3.0, yhi = 3.0;
  const double zlo = -60.0, zhi = 60.0;

  for (int iz = 0; iz < nz; ++iz) {
    const double z = zlo + (zhi - zlo) * iz / (nz - 1);
    for (int ix = 0; ix < nx; ++ix) {
      const double x = xlo + (xhi - xlo) * ix / (nx - 1);
      for (int iy = 0; iy < ny; ++iy) {
        const double y = ylo + (yhi - ylo) * iy / (ny - 1);

        double bfld[3];
        field.FieldInLocalCm(x, y, z, bfld);

        std::fprintf(f, "%10.4f,%10.4f,%10.4f,%16.8E,%16.8E,%16.8E\n", x, y, z, bfld[0], bfld[1],
                     bfld[2]);
      }
    }
  }

  std::fclose(f);
  std::printf("Wrote d1_cpp_reference.csv\n");
  return 0;
}

// Same idea as RunDipoleFieldProbe(), for D2 -- grid matches
// validate/mitray_dipole_probe_d2.f exactly.
int RunD2FieldProbe() {
  const MitrayDipoleData d2data = MitrayDipoleData::D2();
  const G4ThreeVector centerCm(0.0, 0.0, 0.0);
  MitrayDipoleField field(d2data, centerCm);

  std::FILE* f = std::fopen("d2_cpp_reference.csv", "w");
  std::fprintf(f, "x_cm,y_cm,z_cm,bx_T,by_T,bz_T\n");

  const int nx = 7, ny = 5, nz = 25;
  const double xlo = -15.0, xhi = 15.0;
  const double ylo = -3.0, yhi = 3.0;
  const double zlo = -55.0, zhi = 55.0;

  for (int iz = 0; iz < nz; ++iz) {
    const double z = zlo + (zhi - zlo) * iz / (nz - 1);
    for (int ix = 0; ix < nx; ++ix) {
      const double x = xlo + (xhi - xlo) * ix / (nx - 1);
      for (int iy = 0; iy < ny; ++iy) {
        const double y = ylo + (yhi - ylo) * iy / (ny - 1);

        double bfld[3];
        field.FieldInLocalCm(x, y, z, bfld);

        std::fprintf(f, "%10.4f,%10.4f,%10.4f,%16.8E,%16.8E,%16.8E\n", x, y, z, bfld[0], bfld[1],
                     bfld[2]);
      }
    }
  }

  std::fclose(f);
  std::printf("Wrote d2_cpp_reference.csv\n");
  return 0;
}

// Same idea as RunPoleFieldProbe(), for an electrostatic deflector --
// dumps E(x,y,z) (kV/cm) so the CSV can be diffed against the matching
// standalone Fortran probe.
int RunEdipoleFieldProbe(const MitrayEdipoleData& data, const char* outFile, double xlo,
                          double xhi, double ylo, double yhi, double zlo, double zhi) {
  const G4ThreeVector centerCm(0.0, 0.0, 0.0);
  MitrayEdipoleField field(data, centerCm);

  std::FILE* f = std::fopen(outFile, "w");
  std::fprintf(f, "x_cm,y_cm,z_cm,ex_kVcm,ey_kVcm,ez_kVcm\n");

  const int nx = 7, ny = 5, nz = 25;

  for (int iz = 0; iz < nz; ++iz) {
    const double z = zlo + (zhi - zlo) * iz / (nz - 1);
    for (int ix = 0; ix < nx; ++ix) {
      const double x = xlo + (xhi - xlo) * ix / (nx - 1);
      for (int iy = 0; iy < ny; ++iy) {
        const double y = ylo + (yhi - ylo) * iy / (ny - 1);

        double efld[3];
        field.FieldInLocalCm(x, y, z, efld);

        std::fprintf(f, "%10.4f,%10.4f,%10.4f,%16.8E,%16.8E,%16.8E\n", x, y, z, efld[0], efld[1],
                     efld[2]);
      }
    }
  }

  std::fclose(f);
  std::printf("Wrote %s\n", outFile);
  return 0;
}

int RunE1FieldProbe() {
  return RunEdipoleFieldProbe(MitrayEdipoleData::E1(), "e1_cpp_reference.csv", -4.0, 4.0, -2.0,
                               2.0, -30.0, 30.0);
}

int RunE2FieldProbe() {
  return RunEdipoleFieldProbe(MitrayEdipoleData::E2(), "e2_cpp_reference.csv", -4.0, 4.0, -2.0,
                               2.0, -30.0, 30.0);
}

// Reaction config file path (see ReactionConfig.hh for the format,
// reactions/o15ag_19ne.reaction for this pilot's own bundled reaction).
// Mirrors src/ureact.f case(20)'s own convention (its LKINE=20 branch
// reads a similar per-run config file named by the INPUT environment
// variable) -- REACTION_INPUT here, not that exact name, to avoid
// colliding with an unrelated pre-existing use of "INPUT" elsewhere.
// Moved above BuildRunManager (it used to sit just above
// RunReactionTracking) so the automatic-retuning code below, which needs
// it too, can call it.
std::string ReactionFilePath() {
  const char* fromEnv = std::getenv("REACTION_INPUT");
  return fromEnv ? std::string(fromEnv) : std::string("o15ag_19ne.reaction");
}

// --- Automatic "Chain" magnetic retuning for real recoil rigidity --------
//
// DetectorConstruction's own fallback (ComputeMagneticRetuneScale(), absent
// a reaction file's RTUN card) only matches a reaction's *idealized*
// (pre-target-energy-loss) recoil rigidity -- reasonable, but not corrected
// for the real rigidity the recoil actually has by the time it reaches D1
// after crossing the target gas/differential-pumping chain (see
// reactions/o15ag_19ne.reaction's own RTUN card comment: that correction
// was originally measured *by hand* -- build, run `--track-reaction 300`,
// average |p| immediately before D1 over many events, divide by charge
// state, divide by D1's native per-charge rigidity). SetUpRetunedChainGeometry()
// below automates exactly that recipe for RunReactionTracking/RunVis's
// "Reaction" branch (the two entry points that track a reaction's actual
// recoil, as opposed to --track-chain's fixed design-orbit ion or
// --track-beam's beam species) -- so it happens by default for whichever
// reaction REACTION_INPUT currently points to, rather than requiring
// someone to remember to hand-measure and paste in an RTUN card.
//
// This measurement needs real G4 tracking through an idealized-scale
// "Chain" geometry first (to see how much rigidity the recoil actually
// loses crossing real material), then a *second*, corrected geometry to
// actually run the requested (visible-to-the-user) events through --
// naively, two separate G4RunManagers. That crashes: G4RunManager's kernel
// chain (G4RunManagerKernel -> G4EventManager -> G4TrackingManager ->
// G4TrackingMessenger) registers UI commands with G4UImanager's singleton
// command tree, and constructing a second G4RunManagerKernel in the same
// process -- even after fully `delete`-ing the first -- segfaults inside
// G4UImanager::AddNewCommand (confirmed empirically; not a documented
// restriction, but a firm one on this Geant4 build regardless). The fix:
// do BOTH phases on the *same* run manager, using G4RunManager's own
// supported "swap in a new G4VUserDetectorConstruction and rebuild
// geometry" mechanism (SetUserInitialization() again + ReinitializeGeometry(true),
// which explicitly cleans up the old solids/logical/physical volumes
// before invoking the new DetectorConstruction's Construct()) instead of
// a second G4RunManager.

// Clears the per-event "have I already recorded this track's D1 entrance
// momentum" bookkeeping (see RigidityCalibrationStepping) -- track IDs are
// reused (reset to 1) across events, so without this a track ID recorded in
// one event would silently suppress recording the same ID's D1 crossing in
// a later event.
class RigidityCalibrationEvent : public G4UserEventAction {
 public:
  explicit RigidityCalibrationEvent(std::set<G4int>* recordedTrackIds)
      : fRecordedTrackIds(recordedTrackIds) {}
  void BeginOfEventAction(const G4Event*) override { fRecordedTrackIds->clear(); }

 private:
  std::set<G4int>* fRecordedTrackIds;
};

// Silent analogue of SteppingAction, purpose-built for the calibration run:
// no printed output, just the recoil's own momentum magnitude the instant
// it enters D1's field container (mirroring the manual recipe's "row
// immediately preceding a D1-volume row" -- vacuum-to-vacuum, so there's no
// material energy loss right at that boundary to distinguish the two). The
// mass>=30000 MeV cut (below any gamma/e-/e+/proton, well below any of this
// pilot's own recoil ions) is the same trick SteppingAction's TRAJ dump
// output already relies on to separate the recoil's own trajectory from its
// correlated gammas'/secondaries'.
class RigidityCalibrationStepping : public G4UserSteppingAction {
 public:
  RigidityCalibrationStepping(double nominalChargeIonE, std::vector<double>* samplesMeV,
                               std::set<G4int>* recordedTrackIds)
      : fNominalChargeIonE(nominalChargeIonE),
        fSamplesMeV(samplesMeV),
        fRecordedTrackIds(recordedTrackIds) {}

  void UserSteppingAction(const G4Step* step) override {
    G4Track* track = step->GetTrack();
    // Same fix as SteppingAction's own -- without it, G4ionIonisation's
    // effective-charge model would drift this ion's tracked charge away
    // from its real, fixed charge state, corrupting the very rigidity this
    // is trying to measure.
    if (track->GetParticleDefinition()->GetParticleType() == "nucleus") {
      const_cast<G4DynamicParticle*>(track->GetDynamicParticle())->SetCharge(fNominalChargeIonE * eplus);
    }
    if (track->GetDynamicParticle()->GetMass() / MeV < 30000.0) return;  // not the recoil

    const G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    if (preVol && preVol->GetName() == "D1" && fRecordedTrackIds->insert(track->GetTrackID()).second) {
      fSamplesMeV->push_back(step->GetPreStepPoint()->GetMomentum().mag() / MeV);
    }
  }

 private:
  double fNominalChargeIonE;
  std::vector<double>* fSamplesMeV;
  std::set<G4int>* fRecordedTrackIds;
};

// Sets up `runManager` (freshly created, nothing registered on it yet) with
// a "Chain" DetectorConstruction/PhysicsList tuned to the real, measured
// magnetic retune scale for whichever reaction REACTION_INPUT currently
// points to -- or, if that reaction file already has its own RTUN card,
// just that value directly (an explicit override always wins, e.g. for a
// hand-tuned figure someone prefers over the automatic measurement). The
// caller still has to register its own primary generator/stepping/run/
// event actions and call Initialize() + BeamOn() itself -- this only
// leaves the geometry/physics user-initializations in place.
//
// When a measurement is needed: phase 1 builds an *idealized*-scale
// geometry (Q1-Q7 need *some* reasonable scale to reach D1 at all -- fine,
// not circular, since a magnetic field does no work, so |p| at D1's own
// entrance is retune-scale-invariant: whatever scale bends Q1-Q7 changes
// the recoil's direction, not its own momentum magnitude, by the time it
// reaches D1, so any reasonable starting scale gives the same
// measurement), fires 300 real reaction events with no printed output, and
// records the recoil's own momentum the instant it enters D1. Phase 2
// swaps in a fresh DetectorConstruction built with that measured scale and
// calls ReinitializeGeometry(true) to actually rebuild the geometry/fields
// with it (see the "Automatic Chain magnetic retuning" comment above for
// why this has to happen on the same run manager, not a second one).
void SetUpRetunedChainGeometry(G4RunManager* runManager) {
  const ReactionConfig cfg = ReactionConfig::Load(ReactionFilePath());

  if (cfg.magneticFieldRetuneScale > 0.0) {
    runManager->SetUserInitialization(new DetectorConstruction("Chain", cfg.magneticFieldRetuneScale));
    runManager->SetUserInitialization(new PhysicsList());
    return;
  }

  // Phase 1: idealized-scale geometry, just to reach D1 for calibration.
  constexpr int kCalibrationEvents = 300;
  runManager->SetUserInitialization(new DetectorConstruction("Chain"));
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserAction(new PrimaryGeneratorAction(ReactionFilePath(), 0.0,
                                                        TargetChamber::BeamApertureWorldYCm(), 0.0));
  std::vector<double> samplesMeV;
  std::set<G4int> recordedTrackIds;
  runManager->SetUserAction(
      new RigidityCalibrationStepping(cfg.recoilChargeState, &samplesMeV, &recordedTrackIds));
  runManager->SetUserAction(new RigidityCalibrationEvent(&recordedTrackIds));
  runManager->Initialize();
  runManager->BeamOn(kCalibrationEvents);

  const MitrayDipoleData d1 = MitrayDipoleData::D1();
  const double nativeRigidityMeVPerCharge = 0.3 * d1.BF * (d1.RB / 100.0) * 1000.0;

  double scale;
  if (samplesMeV.empty()) {
    // No calibration event reached D1 at all -- e.g. a reaction whose
    // idealized-rigidity recoil doesn't even clear the target chamber's own
    // apertures. Fall back to the purely idealized (pre-target-energy-loss)
    // scale rather than dividing by zero.
    std::fprintf(stderr,
                 "WARNING: SetUpRetunedChainGeometry: no calibration event reached D1 "
                 "(0/%d) -- falling back to idealized (pre-target-energy-loss) rigidity.\n",
                 kCalibrationEvents);
    scale = ComputeMagneticRetuneScale(cfg);
  } else {
    double sumP = 0.0;
    for (double p : samplesMeV) sumP += p;
    const double meanP = sumP / samplesMeV.size();
    scale = (meanP / cfg.recoilChargeState) / nativeRigidityMeVPerCharge;
    std::printf(
        "Auto-measured separator retune scale: %.6f (mean |p| at D1 entrance = %.4f MeV/c, "
        "%zu/%d calibration events reached D1)\n",
        scale, meanP, samplesMeV.size(), kCalibrationEvents);
  }

  // Phase 2: rebuild geometry/fields with the corrected scale. The caller
  // registers its own (real) primary generator/stepping/run/event actions
  // next, then calls Initialize() + BeamOn() -- that Initialize() call is
  // what actually applies the ReinitializeGeometry(true) request below.
  runManager->SetUserInitialization(new DetectorConstruction("Chain", scale));
  runManager->ReinitializeGeometry(true);
}

// Builds and initializes a run manager for `element` ("Q1", "D1", "Chain",
// etc.), wiring up the matching primary generator. Shared by the batch
// (--track*) and interactive/visual (--vis) entry points. Uses
// ComputeMagneticRetuneScale's own idealized-or-RTUN-card fallback for
// "Chain" (not the automatic measurement above) -- --track-chain fires a
// fixed design-orbit ion, not any reaction's actual recoil, so there's no
// real per-reaction rigidity to measure here the way there is for
// RunReactionTracking/RunVis's "Reaction" branch.
G4RunManager* BuildRunManager(const std::string& element, double x0Cm, double y0Cm, double pMeV,
                               double z0Cm) {
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  runManager->SetUserInitialization(new DetectorConstruction(element));
  runManager->SetUserInitialization(new PhysicsList());

  if (element == "Chain") {
    // D1/D2 share one design magnetic rigidity (p/q = B*RB = 0.2156 T.m,
    // i.e. p/q = 64.68 MeV/c per unit e); E1/E2 share one design electric
    // "rigidity" (KE/q = E*RB/2 = 0.4721 MeV per unit e) -- see README. A
    // proton can't satisfy both at once; solving m/q = p^2/(2*KE) gives
    // m/q = 4.758 u/e, i.e. (for q=4e) m = 19.03u -- a 19Ne4+ recoil,
    // consistent with this file's own comments about a 19Ne reaction.
    runManager->SetUserAction(new PrimaryGeneratorAction(x0Cm, y0Cm, pMeV, z0Cm, 10, 19, 4));
    // See SteppingAction's own comment: without this, G4ionIonisation's
    // effective-charge model silently drifts this ion's tracked charge
    // away from its real, fixed 4+ state, corrupting every bend downstream.
    runManager->SetUserAction(new SteppingAction(4.0));
  } else {
    runManager->SetUserAction(new PrimaryGeneratorAction(x0Cm, y0Cm, pMeV, z0Cm));
    runManager->SetUserAction(new SteppingAction());
  }
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());

  runManager->Initialize();
  return runManager;
}

// 19Ne recoil stopping power in the target gas, via the *same*
// G4EmCalculator/PhysicsList physics this pilot already uses for both
// real tracking and GasStoppingPower's beam vertex-depth table -- not an
// independently-sourced number. Reported both in MeV/cm and the
// mass-stopping-power convention MeV/(mg/cm^2), over the recoil's own
// actual kinetic-energy range (see README's "Reaction specification":
// ~1.9 MeV at production, degrading as it crosses the target/pumping
// chain). Cross-checks against the recoil's already-documented ~0.197 MeV
// loss crossing the ~9.9cm CELG target cell: dE/dx(~1.7 MeV) * CELG's own
// areal density (density here * 9.9cm) lands right on that figure.
int RunRecoilDedxProbe() {
  auto* runManager = BuildRunManager("Chain", 0.0, 0.0, 258.7, -120.0);
  // BeamOn(0): forces G4 to finish building cuts/couples for every placed
  // material (including TargetGas) -- Initialize() alone isn't enough;
  // G4EmCalculator::GetDEDX() throws G4Exception em0078 ("FindCouple:
  // fail for material") without this.
  runManager->BeamOn(0);
  G4Material* gas = G4Material::GetMaterial(TargetChamber::TargetGasMaterialName());

  // Material-identity check: is "TargetGas" (see TargetChamber.cc's
  // NistNameFor()) really representing molecular gas (H2/He), not, say,
  // atomic/solid hydrogen? Density and the ICRU mean excitation energy
  // (I-value) both distinguish these -- NIST's own "G4_H" entry (what
  // BuildGas() scales from) is defined as gaseous H2 (I=19.2eV, density
  // 8.3748e-05 g/cm3 at NTP -- ICRU37's own H2-gas value, not atomic H's
  // ~14.995eV/solid-H2's ~21.8eV), which is what should print below.
  std::printf("--- TargetGas material identity ---\n");
  std::printf("State: %s\n", gas->GetState() == kStateGas ? "gas" : "NOT gas (unexpected)");
  std::printf("Density: %.6e g/cm^3 (%.6f mg/cm^3)\n", gas->GetDensity() / (g / cm3),
              gas->GetDensity() / (mg / cm3));
  std::printf("Mean excitation energy (I-value): %.4f eV\n",
              gas->GetIonisation()->GetMeanExcitationEnergy() / eV);
  std::printf("Composition: %d element(s)\n", (int)gas->GetNumberOfElements());
  for (size_t i = 0; i < gas->GetNumberOfElements(); ++i) {
    const G4Element* el = gas->GetElement(i);
    std::printf("  %s (Z=%.0f), mass fraction=%.4f, atoms/volume=%.6e /cm^3\n",
                el->GetName().c_str(), el->GetZ(), gas->GetFractionVector()[i],
                gas->GetVecNbOfAtomsPerVolume()[i] / (1.0 / cm3));
  }
  std::printf("\n");

  // dE/dx probe: the reaction's OWN recoil ion (not hardcoded to 19Ne), so
  // this works for whichever reaction REACTION_INPUT points to (e.g.
  // k39pg_40ca.reaction's 40Ca8+ in H2, not just the bundled 19Ne-in-He
  // case).
  const ReactionConfig cfg = ReactionConfig::Load(ReactionFilePath());
  G4ParticleDefinition* ion =
      G4IonTable::GetIonTable()->GetIon(cfg.recoil.Z, cfg.recoil.A, 0.0);
  G4EmCalculator calc;
  const double densityMgPerCm3 = gas->GetDensity() / (mg / cm3);
  std::printf("--- dE/dx for Z=%d A=%d recoil (charge %d+) in TargetGas ---\n", cfg.recoil.Z,
              cfg.recoil.A, cfg.recoilChargeState);
  const double kes[] = {2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.0, 0.8};
  for (double keMeV : kes) {
    const double dEdxMeVPerCm = calc.GetDEDX(keMeV * MeV, ion, gas) / (MeV / cm);
    const double dEdxMeVPerMgCm2 = dEdxMeVPerCm / densityMgPerCm3;
    std::printf("KE=%.2f MeV: dE/dx = %.6f MeV/cm = %.6f MeV/(mg/cm^2)\n", keMeV, dEdxMeVPerCm,
                dEdxMeVPerMgCm2);
  }
  delete runManager;
  return 0;
}

double DefaultZ0Cm(const std::string& element) {
  if (element == "D1") return -140.0;
  if (element == "E1") return -220.0;  // outside E1's +-250cm bounding box
  if (element == "E2") return -270.0;  // outside E2's +-300cm bounding box
  return -120.0;
}

double DefaultMomentumMeV(const std::string& element) {
  if (element == "Chain") return 258.7;   // 4 x 64.68 MeV/c, see BuildRunManager
  return 15000.0;                         // ~150 MeV/u for A=100-ish
}

int RunTracking(const std::string& element, int argc, char** argv) {
  const double x0Cm = argc > 2 ? std::atof(argv[2]) : 1.0;
  const double y0Cm = argc > 3 ? std::atof(argv[3]) : 0.0;
  const double pMeV = argc > 4 ? std::atof(argv[4]) : DefaultMomentumMeV(element);
  // Optional 7th CLI arg overrides the firing z-position -- lets a single
  // element (e.g. Q1, D1) be probed in isolation by starting the design
  // trajectory right at its own entrance instead of always at z=0.
  const double z0Cm = argc > 6 ? std::atof(argv[6]) : DefaultZ0Cm(element);

  auto* runManager = BuildRunManager(element, x0Cm, y0Cm, pMeV, z0Cm);

  const int nEvents = argc > 5 ? std::atoi(argv[5]) : 1;
  runManager->BeamOn(nEvents);

  delete runManager;
  return 0;
}

// Fires real 15O(alpha,gamma)19Ne reaction events (ReactionKinematics.hh)
// from the target's own vertex (world x=0,z=0, y=TargetChamber's own
// BeamApertureWorldYCm() -- the line connecting the gas cell's actual
// EAPG/XAPG entrance/exit apertures, *not* CELL's own geometric center;
// see that function's own comments) through the full chain, instead of
// the fixed design-orbit recoil --track-chain uses.
// Empirically checks a reaction config file's own cascade (see
// ReactionConfig.hh/ReactionKinematics.hh) by sampling many events with
// no Geant4 run manager at all -- ReactionKinematics::GenerateEvent()
// only needs CLHEP (G4UniformRand/G4ThreeVector), not a constructed
// geometry/physics list/run manager -- and tallying how many gammas each
// one emits, so a hand-edited BRAT/LEVL scheme can be sanity-checked
// (branch fractions, that every chain terminates) before spending time on
// full tracking runs.
int RunReactionStats(int argc, char** argv) {
  const std::string path = argc > 2 ? argv[2] : ReactionFilePath();
  const int nEvents = argc > 3 ? std::atoi(argv[3]) : 10000;

  const ReactionConfig config = ReactionConfig::Load(path);
  const ReactionKinematics reaction(config);

  std::printf("Reaction config: %s\n", path.c_str());
  std::printf("Beam kinetic energy: %.4f MeV\n", reaction.BeamKineticEnergyMeV());

  std::map<int, int> gammaCountTally;
  for (int i = 0; i < nEvents; ++i) {
    const ReactionKinematics::Event ev = reaction.GenerateEvent();
    gammaCountTally[static_cast<int>(ev.gammaMomentaMeV.size())]++;
  }
  for (const auto& [nGammas, count] : gammaCountTally) {
    std::printf("  %d gamma(s)/event: %d/%d events (%.2f%%)\n", nGammas, count, nEvents,
                100.0 * count / nEvents);
  }
  return 0;
}

int RunReactionTracking(int argc, char** argv) {
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  SetUpRetunedChainGeometry(runManager);
  runManager->SetUserAction(
      new PrimaryGeneratorAction(ReactionFilePath(), 0.0, TargetChamber::BeamApertureWorldYCm(), 0.0));
  // See SteppingAction's own comment: without this, G4ionIonisation's
  // effective-charge model silently drifts the recoil's tracked charge
  // away from its real, fixed charge state, corrupting every bend
  // downstream.
  runManager->SetUserAction(
      new SteppingAction(ReactionConfig::Load(ReactionFilePath()).recoilChargeState));
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());
  runManager->Initialize();

  const int nEvents = argc > 2 ? std::atoi(argv[2]) : 1;
  runManager->BeamOn(nEvents);

  delete runManager;
  return 0;
}

// Diagnostic: fires the actual *beam* ion (15O, not the recoil) through
// the target region at a given kinetic energy (default: this pilot's own
// design resonance energy, 258.5413 MeV/c), from just upstream of the gas
// cell, along the same world-z beam axis and at
// TargetChamber::BeamApertureWorldYCm() -- to directly measure the beam's
// own real energy loss crossing the gas with full G4 tracking, e.g. to
// re-derive/sanity-check a reaction config file's own BKIN card (see
// README's "Reaction specification" and GasStoppingPower.hh, which do
// this same calculation analytically/via G4EmCalculator for
// PrimaryGeneratorAction's actual per-event use -- this mode is for
// independently cross-checking that against real tracking, not something
// the simulation itself depends on).
int RunBeamThroughTarget(int argc, char** argv) {
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  runManager->SetUserInitialization(new DetectorConstruction("Chain"));
  runManager->SetUserInitialization(new PhysicsList());
  const double beamMomentumMeV = argc > 2 ? std::atof(argv[2]) : 258.5413;  // 15O at Er=0.5036 MeV
  runManager->SetUserAction(new PrimaryGeneratorAction(
      0.0, TargetChamber::BeamApertureWorldYCm(), beamMomentumMeV, -30.0, 8, 15, 8));
  runManager->SetUserAction(new SteppingAction(8.0));
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());
  runManager->Initialize();

  const int nEvents = argc > 3 ? std::atoi(argv[3]) : 1;
  runManager->BeamOn(nEvents);

  delete runManager;
  return 0;
}

// Interactive/visual session: `dragon_g4_pilot --vis [element] [macro]`.
// Opens a viewer (see macros/init_vis.mac -- picks whatever driver is
// available, Qt/OGL/etc.) and either runs the given macro or drops into
// an interactive prompt (/run/beamOn, /vis/... commands by hand).
// `element` "Reaction" fires real 15O(alpha,gamma)19Ne events (see
// RunReactionTracking) through the same full chain "Chain" builds,
// instead of the fixed design-orbit recoil.
int RunVis(int argc, char** argv) {
  const std::string element = argc > 2 ? argv[2] : "Chain";

  G4RunManager* runManager;
  if (element == "Reaction") {
    runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    SetUpRetunedChainGeometry(runManager);
    runManager->SetUserAction(new PrimaryGeneratorAction(ReactionFilePath(), 0.0,
                                                          TargetChamber::BeamApertureWorldYCm(), 0.0));
    runManager->SetUserAction(
        new SteppingAction(ReactionConfig::Load(ReactionFilePath()).recoilChargeState));
    runManager->SetUserAction(new RunAction());
    runManager->SetUserAction(new EventAction());
    runManager->Initialize();
  } else {
    const double pMeV = DefaultMomentumMeV(element);
    const double z0Cm = DefaultZ0Cm(element);
    runManager = BuildRunManager(element, 0.0, 0.0, pMeV, z0Cm);
  }

  auto* visManager = new G4VisExecutive();
  visManager->Initialize();

  G4UImanager* uiManager = G4UImanager::GetUIpointer();

  if (argc > 3) {
    uiManager->ApplyCommand(G4String("/control/execute ") + argv[3]);
  } else {
    auto* ui = new G4UIExecutive(argc, argv);
    uiManager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 1 && std::strcmp(argv[1], "--vis") == 0) {
    return RunVis(argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe") == 0) {
    return RunQ1FieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-q2") == 0) {
    return RunQ2FieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-q3-q14") == 0) {
    return RunQ3ThroughQ14FieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-dipole") == 0) {
    return RunDipoleFieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-dipole-d2") == 0) {
    return RunD2FieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-e1") == 0) {
    return RunE1FieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-e2") == 0) {
    return RunE2FieldProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-q2") == 0) {
    return RunTracking("Q2", argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-dipole") == 0) {
    return RunTracking("D1", argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-e1") == 0) {
    return RunTracking("E1", argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-e2") == 0) {
    return RunTracking("E2", argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-chain") == 0) {
    return RunTracking("Chain", argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-reaction") == 0) {
    return RunReactionTracking(argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--track-beam") == 0) {
    return RunBeamThroughTarget(argc, argv);
  }
  if (argc > 1 && std::strcmp(argv[1], "--probe-recoil-dedx") == 0) {
    return RunRecoilDedxProbe();
  }
  if (argc > 1 && std::strcmp(argv[1], "--reaction-stats") == 0) {
    return RunReactionStats(argc, argv);
  }
  return RunTracking("Q1", argc, argv);
}
