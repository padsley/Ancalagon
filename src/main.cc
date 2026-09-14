#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "G4EmCalculator.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
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

// Builds and initializes a run manager for `element` ("Q1", "D1", "Chain",
// etc.), wiring up the matching primary generator. Shared by the batch
// (--track*) and interactive/visual (--vis) entry points.
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
  G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(10, 19, 0.0);
  G4EmCalculator calc;
  const double densityMgPerCm3 = gas->GetDensity() / (mg / cm3);
  std::printf("TargetGas density: %.6f mg/cm^3\n", densityMgPerCm3);
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

// Reaction config file path (see ReactionConfig.hh for the format,
// reactions/o15ag_19ne.reaction for this pilot's own bundled reaction).
// Mirrors src/ureact.f case(20)'s own convention (its LKINE=20 branch
// reads a similar per-run config file named by the INPUT environment
// variable) -- REACTION_INPUT here, not that exact name, to avoid
// colliding with an unrelated pre-existing use of "INPUT" elsewhere.
std::string ReactionFilePath() {
  const char* fromEnv = std::getenv("REACTION_INPUT");
  return fromEnv ? std::string(fromEnv) : std::string("o15ag_19ne.reaction");
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

  runManager->SetUserInitialization(new DetectorConstruction("Chain"));
  runManager->SetUserInitialization(new PhysicsList());
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
    runManager->SetUserInitialization(new DetectorConstruction("Chain"));
    runManager->SetUserInitialization(new PhysicsList());
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
