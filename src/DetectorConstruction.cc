#include "DetectorConstruction.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "G4Box.hh"
#include "G4ChordFinder.hh"
#include "G4Colour.hh"
#include "G4ClassicalRK4.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4IntegrationDriver.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "BgoArray.hh"
#include "DsssdSD.hh"
#include "MitrayDipoleField.hh"
#include "MitrayEdipoleField.hh"
#include "ReactionConfig.hh"
#include "ReactionKinematics.hh"
#include "RotateAboutY.hh"
#include "TargetChamber.hh"

namespace {

// Colour coding shared by every mode (standalone and chained) so the
// element types are visually distinguishable: quadrupoles blue, dipoles
// red, electrostatic deflectors green, decorative markers get their own
// colours (see ChainMarker's call sites).
const G4VisAttributes kQuadVis(G4Colour(0.3, 0.5, 1.0));
const G4VisAttributes kDipoleVis(G4Colour(1.0, 0.3, 0.3));
const G4VisAttributes kEdipoleVis(G4Colour(0.3, 1.0, 0.3));
const G4VisAttributes kCollimatorVis(G4Colour(0.85, 0.45, 0.2));  // copper

// --- Retuning the separator for whichever reaction is actually loaded,
// not just the bundled 15O(alpha,gamma)19Ne one --------------------------
//
// Every 'POLE'/'DIPO'/'EDIP' card's field-strength value in
// dat/dragon_2014_DSSSD.dat is set for a 19Ne4+ recoil at exactly 258.7
// MeV/c -- the momentum it has *at the reaction vertex*, before crossing
// any of the target gas, differential-pumping chain, or collimators (see
// README's "Reaction specification": that's literally how the 258.7 MeV/c
// design orbit was identified in the first place, from D1/D2/E1/E2's own
// card values). That's a property of the *specific reaction* the file was
// designed around, not of the separator hardware itself -- a different
// reaction config (different masses, ERES, recoil charge state) produces
// a recoil at a completely different rigidity, and the same fixed
// hardware needs a correspondingly different field strength to still bend
// it by its design angles. ComputeMagneticRetuneScale() below is that
// calculation, generalized: it no longer hardcodes one reaction's own
// empirically-measured ratio.
//
// A real separator operator retunes magnet currents/deflector voltages
// for the beam actually delivered, while leaving the hardware (positions,
// apertures, pole/electrode geometry) untouched -- that's what this does:
// scale only each element's field-STRENGTH parameters (verified linear in
// the field formula itself, see MitrayQuadrupoleField::Bpoles'
// grad1..grad5 / MitrayDipoleField's BF-only-appears-via-"db=BF-br" /
// MitrayEdipoleField's edip() "ef" argument), leaving every geometric
// parameter (aperture radii, Z11..Z22 zone boundaries, RB, PHI,
// ALPHA/BETA, the dimensionless Enge C0..C5 fringe coefficients) exactly
// as the real card specifies. This is why the scaling is applied to a
// *copy* of each element's data at the "Chain" call site, never to the
// canonical MitrayPoleData::Q1()/etc. factories themselves -- those
// remain the untouched, bit-exact-validated real card values the
// standalone tracking modes and --probe* field-grid validation still
// depend on.
//
// For a magnetic element, r=p/(qB) at fixed geometric radius r means B
// must scale linearly with p -- hence the returned scale (the ratio of
// the recoil's own rigidity to D1's native per-charge rigidity) applied
// directly to BQD/BHX/BOC/BDC/BDD (quads) and BF (dipoles). An
// electrostatic deflector's design condition instead balances qE=mv^2/RB,
// i.e. KE/q=E*RB/2 (README's "Tracking spot-checks") -- for fixed
// mass/charge, KE (non-relativistic) scales as p^2, so E1/E2's EFF gets
// the *square* of the same ratio (see the "Chain" call sites).

}  // namespace

// This is declared in DetectorConstruction.hh (external linkage, not
// anonymous-namespace-local like everything else here) so main.cc's
// automatic calibration pass (ResolveChainMagneticRetuneScale()) can reuse
// it -- both to bootstrap its own throwaway calibration geometry (Q1-Q7
// need *some* reasonable scale to reach D1 at all -- see that function's
// own comment for why this idealized one is fine for that) and as its own
// fallback if a calibration run somehow has zero events reach D1.
double ComputeMagneticRetuneScale(const ReactionConfig& cfg) {
  // ReactionConfig's own RTUN card (see ReactionConfig.hh) lets a
  // reaction file supply an empirically-measured scale -- the recoil's
  // real, energy-loss-degraded rigidity by the time it reaches D1, not
  // the idealized vertex value -- exactly the level of refinement
  // reactions/o15ag_19ne.reaction's own RTUN card gives (measured: mean
  // 249.4 MeV/c at D1's entrance over 300 --track-reaction events,
  // against a 258.7 MeV/c idealized vertex value for that reaction's
  // 19Ne4+ recoil -- see that card's own comment for the exact
  // reproduction recipe). Without one, fall back to this reaction's own
  // *idealized* (pre-target-energy-loss) recoil rigidity -- a reasonable
  // first cut (it's the same rigor dat/dragon_2014_DSSSD.dat's own card
  // values were originally set to), just not corrected for real energy
  // loss the way a measured RTUN value is. In practice this fallback is
  // now rarely reached directly -- main.cc's "Chain"-building entry
  // points call ResolveChainMagneticRetuneScale() instead, which measures
  // the real, energy-loss-degraded value automatically (see its own
  // comment) whenever a reaction file omits RTUN, rather than settling
  // for this idealized number the way earlier versions of this pilot did.
  if (cfg.magneticFieldRetuneScale > 0.0) return cfg.magneticFieldRetuneScale;

  const ReactionKinematics reaction(cfg);
  constexpr int kSamples = 2000;  // isotropic emission angle -> converges quickly
  double sumRecoilMomentumMeV = 0.0;
  for (int i = 0; i < kSamples; ++i) {
    sumRecoilMomentumMeV += reaction.GenerateEvent().recoilMomentumMeV.mag();
  }
  const double meanRecoilMomentumMeV = sumRecoilMomentumMeV / kSamples;

  // D1's own native per-charge rigidity (p/q = 0.3*B*RB, MeV/c per unit
  // e) -- the fixed hardware constant every field-strength parameter in
  // the chain is ultimately scaled relative to (D1 rather than any other
  // element only because it's the natural "rigidity" reference the
  // README's own validation already uses).
  const MitrayDipoleData d1 = MitrayDipoleData::D1();
  const double nativeRigidityMeVPerCharge = 0.3 * d1.BF * (d1.RB / 100.0) * 1000.0;

  return (meanRecoilMomentumMeV / cfg.recoilChargeState) / nativeRigidityMeVPerCharge;
}

namespace {

// Diagnostic-only per-quad field-strength trim, on top of magneticScale --
// QN_TRIM_SCALE=<scale> for quad N (1-14), default 1.0 (no change). Added
// to generalize the QSLT-focus retrim investigation (Q2_TRIM_SCALE=1.35 --
// see DetectorConstruction's own "Chain" build site comment) to every
// quad, for checking whether the same kind of achromatic-focus residual
// shows up -- and is similarly fixable -- at the other named foci
// (MSLT/FSLT), not just QSLT.
double QuadTrimScale(int n) {
  char name[32];
  std::snprintf(name, sizeof(name), "Q%d_TRIM_SCALE", n);
  const char* env = std::getenv(name);
  return env ? std::atof(env) : 1.0;
}

MitrayPoleData RetunedQuad(MitrayPoleData d, double magneticScale) {
  d.BQD *= magneticScale;
  d.BHX *= magneticScale;
  d.BOC *= magneticScale;
  d.BDC *= magneticScale;
  d.BDD *= magneticScale;
  return d;
}

MitrayDipoleData RetunedDipole(MitrayDipoleData d, double magneticScale) {
  d.BF *= magneticScale;
  return d;
}

// D2-specific empirical trim, on top of the rigidity retuning above.
//
// Tracing an idealized on-axis particle (zero spread, exactly the
// self-consistent rigidity for its own field strength) through the whole
// chain shows every element within its already-documented ~0.1-0.2deg
// per-element precision (D1: 0.13deg, E1: 0.06deg, every quad <=0.23deg)
// -- except D2, whose own bend angle undershoots its 75deg design PHI by
// 2.06deg, 15-25x every other element's residual, with the box-container
// geometry above (worse still, 3.10deg, with the old sphere -- so this
// isn't a container-shape artifact) and unchanged at 1000x tighter
// chord-finder accuracy (so it isn't a numerical-integration artifact
// either). D2's field values are themselves independently bit-exact
// against the real, unmodified GEANT3 Fortran (see README's "Field
// grids" validation) -- the discrepancy isn't a porting bug in the field
// formula, so this doesn't get "fixed" there. The most likely
// explanation is that the underlying MIT-RAYTRACE analytic field model
// (a 3rd-order paraxial expansion) is simply less accurate for D2's much
// more extreme geometry than D1's -- PHI=75deg vs D1's 50deg, and
// critically ALPHA=BETA=29deg vs D1's 5.8deg (a much larger pole-face
// wedge angle) -- an inherent limitation of the 1970s-era model itself,
// not something introduced by this port.
//
// D2's placement is physically fixed (real, immovable hardware), so the
// only lever available is its own field strength -- exactly how a real
// operator would compensate a magnet that's known to under-deliver its
// nominal bend for a given current. kD2ResidualTrim is that compensation:
// empirically found (a couple of iterations of the same idealized-on-axis
// trace this comment describes) to bring D2's own bend angle to within
// 0.001deg of exactly 75deg for the retuned rigidity above -- i.e. it
// corrects D2's own aberration, not a second independent retuning for
// momentum. Verified effective: the same idealized on-axis particle that
// previously died at RC42 (unretuned) or FC3 (retuned, before this trim)
// now reaches the DSSSD cleanly; over 500 real --track-reaction events,
// DSSSD transmission goes from 0/500 to 4/500 -- modest, because QSLT/MSLT
// still clip the bulk of events on real momentum spread (a separate,
// harder problem -- see the session that derived this comment), but a
// real, reproducible, nonzero improvement this trim alone is responsible
// for, not a coincidence of retuning.
//
// This ratio is itself scale-invariant, so it does not need
// re-measuring per reaction: r=p/(qB) at fixed geometric radius is exactly
// preserved under any simultaneous (p, B) rescaling (the trajectory's own
// equation of motion, in arc-length parametrization, is unchanged), so a
// *fractional* field-strength correction that fixes D2's own aberration
// at one rigidity fixes it at any rigidity -- confirmed directly by
// deriving this value at 19Ne's own retuned rigidity and separately
// verifying it reproduces the same <0.01deg residual there.
constexpr double kD2ResidualTrim = 1.02636;

MitrayDipoleData RetunedD2(MitrayDipoleData d, double magneticScale) {
  d.BF *= magneticScale * kD2ResidualTrim;
  return d;
}

// E1/E2's own retune scale -- NOT simply magneticScale^2 in general (see
// the caller's own comment). dat/dragon_2014_DSSSD.dat's whole beamline
// was carded for a 19Ne4+ recoil at 258.7 MeV/c (README's "Reaction
// specification") -- kDesignRecoilMassMeV/kDesignRecoilChargeState are
// that ion's own real ground-state mass (AME mass-excess convention, same
// formula ReactionKinematics uses: A*amu + massExcessMeV, from
// reactions/o15ag_19ne.reaction's own RECL card) and charge state, not
// this pilot's invention. For a recoil of a *different* mass and/or
// charge state than that (e.g. 40Ca, in reactions/k39pg_40ca.reaction),
// the correction below is required: KE/q (what E1/E2's own design
// condition, qE=mv^2/RB, actually constrains) does NOT scale the same way
// magnetic rigidity p/q does when mass differs, since
// KE = p^2/(2m) -- confirmed empirically (see the session that found
// this: a --track-reaction charge-state scan for k39pg_40ca showed E1's
// own bend angle drifting smoothly from -28deg to -9.5deg against its
// 20deg design as the assigned recoil charge state ran 5->16, even though
// every element's *magnetic* retuning was already correct for each
// charge's own rigidity -- only charge 8 (this reaction's actual
// production charge state, coincidentally close to the ~8.4 value that
// makes the old magneticScale^2 shortcut accidentally near-correct) looked
// fine). The correction factor below is exactly 1 when mass/charge match
// the original 19Ne4+ design (reactions/o15ag_19ne.reaction), so this is
// not a behavior change for that reaction -- only for any other one.
constexpr double kDesignRecoilMassMeV = 19 * 931.49432 + 1.7511;  // 19Ne ground state
constexpr int kDesignRecoilChargeState = 4;                       // 19Ne4+

double ComputeElectricRetuneScale(const ReactionConfig& cfg, double magneticScale) {
  const double recoilMassMeV = ReactionKinematics(cfg).RecoilGroundMassMeV();
  return magneticScale * magneticScale * (kDesignRecoilMassMeV / recoilMassMeV) *
         (static_cast<double>(cfg.recoilChargeState) / kDesignRecoilChargeState);
}

MitrayEdipoleData RetunedEdipole(MitrayEdipoleData d, double electricScale) {
  d.EFF *= electricScale;
  return d;
}

// --- Standalone-mode builders (unchanged; theta=0 always) -----------------
//
// name: GEANT volume name ("Q1" or "Q2"). centerXCm/centerZCm: world
// position of the tube's centre, which is also the field's origin (za=0)
// -- GEANT3's ugeo_mpole places the tube exactly there too.
void BuildQuad(G4LogicalVolume* worldLV, G4Material* vacuum, const char* name,
               const MitrayPoleData& data, double centerXCm, double centerZCm) {
  const G4ThreeVector centerCm(centerXCm, 0.0, centerZCm);

  auto* field = new MitrayQuadrupoleField(data, centerCm);
  auto* eqRhs = new G4Mag_UsualEqRhs(field);
  auto* stepper = new G4ClassicalRK4(eqRhs);
  auto* fieldManager = new G4FieldManager(field);
  fieldManager->SetChordFinder(new G4ChordFinder(field, 1.0e-3 * mm, stepper));

  const double radius = field->ApertureRadiusCm() * cm;
  const double halfLength = field->HalfLengthCm() * cm;
  auto* solid = new G4Tubs(name, 0.0, radius, halfLength, 0.0, 360.0 * deg);
  auto* lv = new G4LogicalVolume(solid, vacuum, name);
  lv->SetVisAttributes(kQuadVis);
  new G4PVPlacement(nullptr, G4ThreeVector(centerXCm * cm, 0.0, centerZCm * cm), lv, name, worldLV,
                     false, 0, true);

  lv->SetFieldManager(fieldManager, true);
}

// originXCm/originZCm: world position of D1's field origin (za=0). An
// arbitrary centre point in standalone "D1" mode; box centred on it
// (symmetric margin on both sides -- fine in isolation).
void BuildD1(G4LogicalVolume* worldLV, G4Material* vacuum, double originXCm, double originZCm) {
  const MitrayDipoleData d1data = MitrayDipoleData::D1();
  const G4ThreeVector d1OriginCm(originXCm, 0.0, originZCm);

  auto* field = new MitrayDipoleField(d1data, d1OriginCm);
  auto* eqRhs = new G4Mag_UsualEqRhs(field);
  auto* stepper = new G4ClassicalRK4(eqRhs);
  auto* fieldManager = new G4FieldManager(field);
  fieldManager->SetChordFinder(new G4ChordFinder(field, 1.0e-3 * mm, stepper));

  // Generous bounding box (not a port of ugeo_dipole's exact TRAP wedge --
  // see header comment): big enough in x/z to contain the whole bent
  // trajectory (RB=100cm, PHI=50deg) plus fringe padding, and in y to
  // clear the pole gap.
  auto* d1Solid = new G4Box("D1", 80.0 * cm, 15.0 * cm, 100.0 * cm);
  auto* d1LV = new G4LogicalVolume(d1Solid, vacuum, "D1");
  d1LV->SetVisAttributes(kDipoleVis);
  new G4PVPlacement(nullptr, G4ThreeVector(originXCm * cm, 0.0, originZCm * cm), d1LV, "D1",
                     worldLV, false, 0, true);

  d1LV->SetFieldManager(fieldManager, true);
}

// name: GEANT volume name ("E1" or "E2"). originXCm/originZCm: world
// position of the field's origin (za=0, an arbitrary centre point in
// standalone mode). An electric field does work on the particle, so this
// needs G4EqMagElectricField (not G4Mag_UsualEqRhs) and an 8-variable
// stepper.
void BuildEdipole(G4LogicalVolume* worldLV, G4Material* vacuum, const char* name,
                   const MitrayEdipoleData& data, double originXCm, double originZCm) {
  const G4ThreeVector originCm(originXCm, 0.0, originZCm);

  auto* field = new MitrayEdipoleField(data, originCm);
  auto* eqRhs = new G4EqMagElectricField(field);
  const G4int nvar = 8;
  auto* stepper = new G4ClassicalRK4(eqRhs, nvar);
  auto* fieldManager = new G4FieldManager(field);
  auto* driver = new G4IntegrationDriver<G4ClassicalRK4>(1.0e-3 * mm, stepper,
                                                          stepper->GetNumberOfVariables());
  fieldManager->SetChordFinder(new G4ChordFinder(driver));

  // Generous bounding box (not a port of ugeo_edipol's exact TRAP wedge --
  // see header comment), sized off RB/PHI so it comfortably contains the
  // bend regardless of exactly where za=0 sits relative to it (RB+50cm in
  // both x and z bounds the full bend circle from either end with margin
  // to spare -- e.g. E1: RB=200cm -> +-250cm; E2: RB=250cm -> +-300cm).
  const double halfExtentCm = data.RB + 50.0;
  auto* solid = new G4Box(name, halfExtentCm * cm, 10.0 * cm, halfExtentCm * cm);
  auto* lv = new G4LogicalVolume(solid, vacuum, name);
  lv->SetVisAttributes(kEdipoleVis);
  new G4PVPlacement(nullptr, originCm * cm, lv, name, worldLV, false, 0, true);

  lv->SetFieldManager(fieldManager, true);
}

// --- Full-beamline chain builder -------------------------------------------
//
// Walks dat/dragon_2014_DSSSD.dat's own card sequence -- DRIFT/SHRT advance
// a running (x, z, theta) state exactly as src/mitray_setup.f's
// ugeom_setup does (~lines 606-665); a quadrupole card advances it by
// A+B+L, centred at A+(Z22+L-Z11)/2 from entry (see BuildQuad's comment --
// every quad in this file has Z11==Z22, so that's just L/2). No SHRT in
// this file ever rotates (every angle slot is zero throughout), so the
// only rotation source is a bend at a dipole or electrostatic deflector:
// each one bends the *design* trajectory by PHI along a circular arc of
// radius RB (toward -x, matching this file's sign convention -- verified
// empirically for D1/E1 in the standalone pilots), so its local exit point
// is the standard circle relation (-RB(1-cosPHI), 0, RB sinPHI) and theta
// decreases by PHI. This does NOT replicate GEANT3's own dipole volume
// placement algorithm (ugeom_setup ~line 757-820, the curvature-centre/
// arc-midpoint construction) -- as with the standalone D1 pilot, a
// bending element's field origin here is simply its drift-accumulated
// entry point, which is a documented placement simplification, not a
// physics one (the field itself is still bit-exact validated against the
// real Fortran).
struct ChainState {
  double xCm = 0.0;
  double zCm = 0.0;
  double thetaDeg = 0.0;
};

void Drift(ChainState& s, double lengthCm) {
  const double th = s.thetaDeg * CLHEP::pi / 180.0;
  s.xCm += lengthCm * std::sin(th);
  s.zCm += lengthCm * std::cos(th);
}

void Shift(ChainState& s, double dxCm) {
  const double th = s.thetaDeg * CLHEP::pi / 180.0;
  s.xCm += dxCm * std::cos(th);
  s.zCm -= dxCm * std::sin(th);
}

// Diagnostic-only: DUMP_CHAIN_GEOMETRY=1 prints every named element's own
// reference-axis state (world x/z, local theta) as it's placed -- for
// reconstructing the design axis as an analytic line/arc polyline (for a
// beam-envelope plot, e.g.), without hand-transcribing Drift/Shift lengths
// from source. Unset by default, so omitting it changes nothing.
void DumpChainState(const char* tag, const ChainState& s) {
  if (!std::getenv("DUMP_CHAIN_GEOMETRY")) return;
  std::fprintf(stderr, "DUMP_CHAIN %-8s xCm=%.6f zCm=%.6f thetaDeg=%.6f\n", tag, s.xCm, s.zCm,
               s.thetaDeg);
}

// Placement rotation shared by every rotated volume past a bend -- real
// absorbing geometry (collimators) and the quad/dipole/edipole field
// containers alike, now that none of them fall back to an orientation-
// agnostic G4Orb any more (see ChainQuad's/ChainDipole's own comments for
// why that mattered for correctness, not just looks). Builds
// the inverse of RotateAboutY.hh's local->world active rotation (world =
// R(theta)*local) via CLHEP's rotateAxes(...) taking that rotation's own
// ROWS (same "pass rows for the inverse" trick already used in
// BgoArray.cc/TargetChamber.cc, since R(theta) is orthogonal).
G4RotationMatrix* BeamAxisRotation(double thetaDeg) {
  const double th = thetaDeg * CLHEP::pi / 180.0;
  const double c = std::cos(th), sn = std::sin(th);
  auto* rot = new G4RotationMatrix();
  rot->rotateAxes(G4ThreeVector(c, 0, sn), G4ThreeVector(0, 1, 0), G4ThreeVector(-sn, 0, c));
  return rot;
}

void AttachMagnetic(G4LogicalVolume* lv, G4MagneticField* field) {
  auto* eqRhs = new G4Mag_UsualEqRhs(field);
  auto* stepper = new G4ClassicalRK4(eqRhs);
  auto* fieldManager = new G4FieldManager(field);
  fieldManager->SetChordFinder(new G4ChordFinder(field, 1.0e-3 * mm, stepper));
  lv->SetFieldManager(fieldManager, true);
}

void AttachElectric(G4LogicalVolume* lv, G4ElectricField* field) {
  auto* eqRhs = new G4EqMagElectricField(field);
  auto* stepper = new G4ClassicalRK4(eqRhs, 8);
  auto* fieldManager = new G4FieldManager(field);
  auto* driver = new G4IntegrationDriver<G4ClassicalRK4>(1.0e-3 * mm, stepper,
                                                          stepper->GetNumberOfVariables());
  fieldManager->SetChordFinder(new G4ChordFinder(driver));
  lv->SetFieldManager(fieldManager, true);
}

// Places a quadrupole at the chain's current entry point and advances the
// chain by its pass-through length A+B+L; theta is unchanged (a
// quadrupole doesn't bend the design trajectory). The geometry is always
// the same exact-fit G4Tubs as BuildQuad(), rotated onto the local beam
// axis past a bend via BeamAxisRotation() -- the same helper
// ChainCollimator() already validates. A quad's container used to fall
// back to an orientation-agnostic G4Orb (bounding sphere of the tube)
// once rotated, to sidestep G4PVPlacement's rotation-convention gotcha;
// that was a correctness bug, not just a cosmetic one -- a sphere of
// radius sqrt(RAD^2+halfLength^2) only matches the tube's real
// half-length exactly on-axis, and narrows at any transverse offset, so
// an off-design-orbit trajectory (exactly the population most exposed to
// transmission loss) could exit the field container well short of the
// real magnet's physical extent. A G4Tubs is also a solid of revolution
// about its own axis, so -- unlike the collimator boxes -- there's no
// roll-angle ambiguity to get wrong even if BeamAxisRotation's convention
// were subtly off: only its axis direction (already validated by the
// existing centerXCm/centerZCm sin/cos convention) matters.
//
// maxExtentCm caps the container's half-length (BEFORE that -- the
// natural (L+Z11+Z22)/2 fringe-margin size -- routinely exceeds the real
// gap to this element's neighbour, see below): with G4's strict
// non-overlapping-sibling-volume model, two adjacent field-managed
// containers that overlap don't "both apply" in the overlap region --
// Geant4's navigator resolves that whole region to whichever volume it
// placed first, so the OTHER element's field is silently never applied
// there at all. Empirically, before this cap existed, 13 of this chain's
// 18 quad/dipole/e-dipole elements never registered a single step for the
// design-orbit trajectory (their containers overlapped their upstream
// neighbour's). The elements themselves can't move -- they're the real,
// fixed physical arrangement of the separator -- so each cap here is
// precomputed (offline, from this same chain's own fixed entry/exit
// geometry) as half the centre-to-centre distance to whichever neighbour
// is closest, which guarantees zero overlap between any two adjacent
// containers; see the call sites in Construct() for the actual numbers.
// This does shrink how much of each element's fringe field is
// geometrically captured, but every cap here still clears that element's
// own physical core half-length (data.L/2) with room to spare.
void ChainQuad(ChainState& s, G4LogicalVolume* worldLV, G4Material* vacuum, const char* name,
               const MitrayPoleData& data, double maxExtentCm) {
  DumpChainState(name, s);
  const double th = s.thetaDeg * CLHEP::pi / 180.0;
  const double entryToCentreCm = data.A + (data.Z22 + data.L - data.Z11) / 2.0;
  const double centerXCm = s.xCm + entryToCentreCm * std::sin(th);
  const double centerZCm = s.zCm + entryToCentreCm * std::cos(th);

  auto* field =
      new MitrayQuadrupoleField(data, G4ThreeVector(centerXCm, 0.0, centerZCm), s.thetaDeg);

  const double halfLengthCm = std::min((data.L + data.Z11 + data.Z22) / 2.0, maxExtentCm);
  auto* solid = new G4Tubs(name, 0.0, data.RAD * cm, halfLengthCm * cm, 0.0, 360.0 * deg);
  auto* lv = new G4LogicalVolume(solid, vacuum, name);
  lv->SetVisAttributes(kQuadVis);
  G4RotationMatrix* rot = BeamAxisRotation(s.thetaDeg);
  new G4PVPlacement(rot, G4ThreeVector(centerXCm * cm, 0.0, centerZCm * cm), lv, name, worldLV,
                     false, 0, true);
  AttachMagnetic(lv, field);

  const double passThroughCm = data.A + data.B + data.L;
  s.xCm += passThroughCm * std::sin(th);
  s.zCm += passThroughCm * std::cos(th);
}

// Places a dipole at the chain's current entry point (field origin = the
// entry point, A-axis aligned to the chain's current theta, unaffected by
// anything below -- the physics is keyed off this, not the container),
// then advances the chain to the design-orbit exit point and rotates
// theta by -PHI for everything downstream.
//
// Container: a chord-aligned G4Box, centred on the *arc's own midpoint*
// (at PHI/2), not the entry point -- unlike a quad, a dipole's container
// must be big enough to contain its whole bend chord (2*RB*sin(PHI/2)),
// which for a wide-angle bend like D1's 50 degrees is far bigger than the
// gap to its nearest neighbour (Q2) if centred at the entry point (the
// earlier RB+50cm-at-entry-point formula did exactly that, and is exactly
// why D1 swallowed Q2 -- see ChainQuad's own comment on the general
// problem). Centring at the arc's midpoint instead roughly halves the
// needed reach for the same chord (RB*sqrt(2*(1-cos(PHI/2))) instead of
// the full chord), and -- since the midpoint sits further from both this
// element's neighbours than either endpoint does -- leaves much more of
// containerRadiusCm's own margin for the aperture width and fringe decay
// beyond the idealized zero-width arc. containerRadiusCm is precomputed
// the same way as ChainQuad's maxExtentCm (that formula plus a margin,
// further capped to at most half the centre-to-centre distance to
// whichever neighbour is closest to the midpoint) -- see the call sites.
void ChainDipole(ChainState& s, G4LogicalVolume* worldLV, G4Material* vacuum, const char* name,
                 const MitrayDipoleData& data, double containerRadiusCm) {
  DumpChainState(name, s);
  const double th = s.thetaDeg * CLHEP::pi / 180.0;
  const double originXCm = s.xCm, originZCm = s.zCm;

  auto* field =
      new MitrayDipoleField(data, G4ThreeVector(originXCm, 0.0, originZCm), s.thetaDeg);

  const double phiRad = data.PHI * CLHEP::pi / 180.0;
  const double halfPhiRad = 0.5 * phiRad;
  const double midXLocalCm = -data.RB * (1.0 - std::cos(halfPhiRad));
  const double midZLocalCm = data.RB * std::sin(halfPhiRad);
  const double midXCm = originXCm + std::cos(th) * midXLocalCm + std::sin(th) * midZLocalCm;
  const double midZCm = originZCm - std::sin(th) * midXLocalCm + std::cos(th) * midZLocalCm;

  // Chord-aligned box instead of an orientation-agnostic sphere. The
  // straight chord from entry to exit bisects the bend angle exactly (a
  // circular arc's standard chord/tangent identity -- also confirmed
  // directly from ChainDipole's own exit-point formula below), so its
  // world-frame direction is just the entry theta minus half of PHI;
  // BeamAxisRotation() at that angle is the same helper/convention
  // ChainQuad/ChainCollimator already validate, just evaluated at the
  // bisector instead of the post-bend angle. Half-length along the chord
  // reuses containerRadiusCm as-is (the same precomputed, neighbour-
  // overlap-safe reach the sphere used); the transverse/vertical
  // half-extents come from the field's own aperture (WDIP1/WDIP2, D)
  // plus a flat margin for the arc's sagitta and fringe falloff, each
  // still capped at containerRadiusCm so this can only shrink the sphere's
  // old footprint, never grow it.
  const double chordThetaDeg = s.thetaDeg - data.PHI / 2.0;
  G4RotationMatrix* rot = BeamAxisRotation(chordThetaDeg);
  const double halfWidthCm =
      std::min(std::max(data.WDIP1, data.WDIP2) / 2.0 + 10.0, containerRadiusCm);
  const double halfHeightCm = std::min(data.D / 2.0 + 10.0, containerRadiusCm);
  auto* solid = new G4Box(name, halfWidthCm * cm, halfHeightCm * cm, containerRadiusCm * cm);
  auto* lv = new G4LogicalVolume(solid, vacuum, name);
  lv->SetVisAttributes(kDipoleVis);
  new G4PVPlacement(rot, G4ThreeVector(midXCm * cm, 0.0, midZCm * cm), lv, name, worldLV, false,
                     0, true);
  AttachMagnetic(lv, field);

  const double exitXLocalCm = -data.RB * (1.0 - std::cos(phiRad));
  const double exitZLocalCm = data.RB * std::sin(phiRad);
  s.xCm = originXCm + std::cos(th) * exitXLocalCm + std::sin(th) * exitZLocalCm;
  s.zCm = originZCm - std::sin(th) * exitXLocalCm + std::cos(th) * exitZLocalCm;
  s.thetaDeg -= data.PHI;
  char exitTag[32];
  std::snprintf(exitTag, sizeof(exitTag), "%s_exit", name);
  DumpChainState(exitTag, s);
}

// Same pattern as ChainDipole(), for an electrostatic deflector (also
// bends the design trajectory by PHI along a radius-RB arc -- verified
// empirically in the standalone E1/E2 pilots). Same arc-midpoint
// container centring, same chord-aligned box in place of the old sphere
// (see ChainDipole's comment) -- MitrayEdipoleData has no separate
// WDIP1/WDIP2, D itself is the deflector's own horizontal aperture bound
// (MitrayEdipoleField::FieldInLocalCm's xbmax/xcmax), so it stands in for
// both the transverse and vertical half-extent here.
void ChainEdipole(ChainState& s, G4LogicalVolume* worldLV, G4Material* vacuum, const char* name,
                  const MitrayEdipoleData& data, double containerRadiusCm) {
  DumpChainState(name, s);
  const double th = s.thetaDeg * CLHEP::pi / 180.0;
  const double originXCm = s.xCm, originZCm = s.zCm;

  auto* field =
      new MitrayEdipoleField(data, G4ThreeVector(originXCm, 0.0, originZCm), s.thetaDeg);

  const double phiRad = data.PHI * CLHEP::pi / 180.0;
  const double halfPhiRad = 0.5 * phiRad;
  const double midXLocalCm = -data.RB * (1.0 - std::cos(halfPhiRad));
  const double midZLocalCm = data.RB * std::sin(halfPhiRad);
  const double midXCm = originXCm + std::cos(th) * midXLocalCm + std::sin(th) * midZLocalCm;
  const double midZCm = originZCm - std::sin(th) * midXLocalCm + std::cos(th) * midZLocalCm;

  const double chordThetaDeg = s.thetaDeg - data.PHI / 2.0;
  G4RotationMatrix* rot = BeamAxisRotation(chordThetaDeg);
  const double halfWidthCm = std::min(data.D / 2.0 + 10.0, containerRadiusCm);
  const double halfHeightCm = halfWidthCm;
  auto* solid = new G4Box(name, halfWidthCm * cm, halfHeightCm * cm, containerRadiusCm * cm);
  auto* lv = new G4LogicalVolume(solid, vacuum, name);
  lv->SetVisAttributes(kEdipoleVis);
  new G4PVPlacement(rot, G4ThreeVector(midXCm * cm, 0.0, midZCm * cm), lv, name, worldLV, false,
                     0, true);
  AttachElectric(lv, field);

  const double exitXLocalCm = -data.RB * (1.0 - std::cos(phiRad));
  const double exitZLocalCm = data.RB * std::sin(phiRad);
  s.xCm = originXCm + std::cos(th) * exitXLocalCm + std::sin(th) * exitZLocalCm;
  s.zCm = originZCm - std::sin(th) * exitXLocalCm + std::cos(th) * exitZLocalCm;
  s.thetaDeg -= data.PHI;
  char exitTag[32];
  std::snprintf(exitTag, sizeof(exitTag), "%s_exit", name);
  DumpChainState(exitTag, s);
}

// Places a thin, labeled, purely-decorative wireframe marker box at the
// chain's current position (no state advance) for the two 'TEST' cards
// (TST3, TST4) that mark where MCP0/MCP1 sit -- 'TEST' cards carry no
// data at all in this file (no aperture, no material), so unlike every
// 'RCOL' card (see ChainCollimator below) there is nothing real to port
// here beyond the position itself. Not rotated to the chain's current
// theta, for the same reason the field-container spheres aren't (see
// ChainQuad) -- doubly unnecessary here since it isn't real absorbing
// geometry anyway.
void ChainMarker(ChainState& s, G4LogicalVolume* worldLV, G4Material* vacuum, const char* name,
                  double halfXCm, double halfYCm, double halfThicknessCm,
                  const G4VisAttributes& vis) {
  auto* solid = new G4Box(name, halfXCm * cm, halfYCm * cm, halfThicknessCm * cm);
  auto* lv = new G4LogicalVolume(solid, vacuum, name);
  lv->SetVisAttributes(vis);
  new G4PVPlacement(nullptr, G4ThreeVector(s.xCm * cm, 0.0, s.zCm * cm), lv, name, worldLV, false,
                     0, true);
}

// Real collimator/slit geometry, ported from src/ugeom_mitray.f's
// ugeo_col(pos,irot,data,rname) -- the only aperture cards in this file
// that get real solid geometry: every 'RCOL' card (idata==17), unlike
// 'COLL' cards (idata==13), which src/mitray_setup.f only ever uses for
// acceptance-limiting bookkeeping (line 667-707, feeding jcol/xcol/ycol/
// dxcol/dycol into whichever optical element they bracket) -- ugeo_col
// itself is called only for idata==17 (line 709). Material is copper for
// every one of them (ugstmed.f medium 5, 'COPPER'; the file never varies
// this by rname). Neither RCOL nor COLL cards advance the running
// position (zero-length in ugeom_setup), so this is called at the
// chain's *current* (x,z,theta) with no state change, same as
// ChainMarker.
//
// isCircular is the card's own DATA(1)!=0: an annular copper tube
// (rmin=dataXCm, rmax=dataYCm, half-length halfZCm). Otherwise
// (DATA(1)==0) a rectangular aperture, built exactly as ugeo_col does:
// four identical copper jaw boxes (half-dims magnification*dataXCm,
// magnification*dataYCm, halfZCm -- ugeo_col's own "block size =
// magnification * aperture half-width" rule), one flanking each side of
// the rectangular hole (+-local-x, +-local-y). ugeo_col keys the
// magnification by rname: 20 for QSLT, 10 for MSLT/FSLT, 1 for every
// other rectangular RCOL card in this file (FC1-FC4); circular cards
// ignore it. Because all four jaws share the same box dimensions rather
// than mitered corners, adjacent jaws genuinely overlap at the corners
// whenever magnification*data is much bigger than data itself (true for
// QSLT/MSLT/FSLT) -- exactly how the real file builds it, so overlap
// checking is deliberately off for these four placements only; it stays
// on for placement against everything else in the chain.
//
// Past a bend (theta!=0) this geometry, like the quad/dipole/edipole
// field containers, must track the local beam axis (see BeamAxisRotation's
// own comment), so every placement here carries a real rotation.
// offsetXCm/offsetYCm are the card's own DATA(2)/DATA(3) (always 0 in
// this file, but ported for fidelity to ugeo_col's own GTRMUL offset, not
// hardcoded).
//
// Expect (and this is not a placement bug): with pCheckOverlaps on, a
// handful of these -- 16 in the current "Chain" build, e.g. E1/E2's own
// containers each overlapping one neighbouring collimator -- still report
// overlapping a quad/dipole/edipole field container. That's unavoidable
// given how tightly this beamline is packed (containerRadiusCm is itself
// capped to half the centre-to-centre distance to the *nearest* neighbour,
// so a container reaching that far can still catch a collimator sitting
// off to one side rather than directly ahead); switching those containers
// from spheres to chord-aligned boxes (see ChainDipole's/ChainQuad's own
// comments) only shrinks this set, never grows it. Harmless regardless:
// these containers are vacuum field volumes with no real material
// boundary for a physics process to resolve differently at the overlap.
void ChainCollimator(const ChainState& s, G4LogicalVolume* worldLV, G4Material* copper,
                      const char* name, bool isCircular, double offsetXCm, double offsetYCm,
                      double dataXCm, double dataYCm, double halfZCm, double magnification = 1.0) {
  DumpChainState(name, s);
  const double thetaRad = s.thetaDeg * CLHEP::pi / 180.0;
  G4RotationMatrix* rot = BeamAxisRotation(s.thetaDeg);

  auto worldPos = [&](double localXCm, double localYCm) {
    double wx, wy, wz;
    RotateLocalToWorld(localXCm, localYCm, 0.0, thetaRad, wx, wy, wz);
    return G4ThreeVector((s.xCm + wx) * cm, wy * cm, (s.zCm + wz) * cm);
  };

  if (isCircular) {
    auto* solid = new G4Tubs(name, dataXCm * cm, dataYCm * cm, halfZCm * cm, 0.0, 360.0 * deg);
    auto* lv = new G4LogicalVolume(solid, copper, name);
    lv->SetVisAttributes(kCollimatorVis);
    new G4PVPlacement(rot, worldPos(offsetXCm, offsetYCm), lv, name, worldLV, false, 0, true);
    return;
  }

  const double halfXCm = magnification * dataXCm;
  const double halfYCm = magnification * dataYCm;
  auto* solid = new G4Box(name, halfXCm * cm, halfYCm * cm, halfZCm * cm);

  auto place = [&](const char* suffix, double localXCm, double localYCm) {
    const G4String lvName = G4String(name) + suffix;
    auto* lv = new G4LogicalVolume(solid, copper, lvName);
    lv->SetVisAttributes(kCollimatorVis);
    new G4PVPlacement(rot, worldPos(localXCm, localYCm), lv, lvName, worldLV, false, 0, false);
  };
  place("_X+", offsetXCm + dataXCm + halfXCm, offsetYCm);
  place("_X-", offsetXCm - dataXCm - halfXCm, offsetYCm);
  place("_Y+", offsetXCm, offsetYCm + dataYCm + halfYCm);
  place("_Y-", offsetXCm, offsetYCm - dataYCm - halfYCm);
}

}  // namespace

namespace {
// Diagnostic-only: without a max-step constraint, G4's own adaptive
// stepper can take one huge step across an entire smooth-field container
// (a quad/dipole/edipole's own bounding volume can be tens of cm with no
// intervening geometric boundary) -- see the finding at this call's own
// use sites: for at least Q1/Q2/D1 in "Chain" mode, that single big step
// was not just under-sampled but produced a genuinely different final
// position/momentum than the same trajectory computed with small steps,
// well beyond anything explained by the 1e-3mm chord-accuracy target
// already configured on every field's own G4ChordFinder. FINE_STEP_CM
// caps every logical volume's own max step size (via G4UserLimits, which
// needs PhysicsList.cc's own G4StepLimiter process registration to do
// anything at all); unset by default, so omitting it changes nothing.
void ApplyFineStepIfRequested() {
  const char* fineStepEnv = std::getenv("FINE_STEP_CM");
  if (!fineStepEnv) return;
  const double fineStepCm = std::atof(fineStepEnv);
  auto* limits = new G4UserLimits(fineStepCm * cm);
  for (auto* lv : *G4LogicalVolumeStore::GetInstance()) {
    lv->SetUserLimits(limits);
  }
}
}  // namespace

G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

  // The real beamline folds through a cumulative -180 degrees by Q14 (50
  // at D1 + 20 at E1 + 75 at D2 + 35 at E2), so the chain's world-frame
  // footprint is genuinely large -- roughly x in [-1030, 50] cm, z in
  // [0, 700] cm by Q14 -- not a sign of a bug. Generously oversized so
  // every element (including E2's own +-300cm bounding sphere) fits with
  // margin regardless of exactly where it lands.
  const double worldHalfZ = 16.0 * m;
  const double worldHalfXY = 16.0 * m;
  auto* worldSolid = new G4Box("World", worldHalfXY, worldHalfXY, worldHalfZ);
  auto* worldLV = new G4LogicalVolume(worldSolid, vacuum, "World");
  auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);

  if (fElement == "D1") {
    BuildD1(worldLV, vacuum, 0.0, kD1CenterZCm);
    ApplyFineStepIfRequested();
    return worldPV;
  }

  if (fElement == "Q2") {
    BuildQuad(worldLV, vacuum, "Q2", MitrayPoleData::Q2(), 0.0, kQ1CenterZCm);
    ApplyFineStepIfRequested();
    return worldPV;
  }

  if (fElement == "E1") {
    BuildEdipole(worldLV, vacuum, "E1", MitrayEdipoleData::E1(), 0.0, kE1CenterZCm);
    ApplyFineStepIfRequested();
    return worldPV;
  }

  if (fElement == "E2") {
    BuildEdipole(worldLV, vacuum, "E2", MitrayEdipoleData::E2(), 0.0, kE2CenterZCm);
    ApplyFineStepIfRequested();
    return worldPV;
  }

  if (fElement == "Chain") {
    // Q1 -> Q2 -> D1 -> Q3 -> Q4 -> Q5 -> Q6 -> Q7 -> E1 -> Q8 -> Q9 ->
    // Q10 -> D2 -> Q11 -> Q12 -> E2 -> Q13 -> Q14, spaced and oriented
    // exactly as dat/dragon_2014_DSSSD.dat's real card sequence (line
    // numbers below refer to that file). 'COLL'/'RCOL' cards are aperture
    // cuts only (zero length); 'TEST'/'FCUP' cards net to zero advance in
    // ugeom_setup and are omitted.
    // Whichever reaction config file is actually loaded (same
    // REACTION_INPUT env var convention as main.cc's own
    // ReactionFilePath(), duplicated here in miniature rather than
    // shared, since this is the only other call site) drives two things
    // below: the target chamber's own gas species, and the separator's
    // field-strength retuning.
    const char* reactionInputEnv = std::getenv("REACTION_INPUT");
    const std::string reactionPath =
        reactionInputEnv ? std::string(reactionInputEnv) : std::string("o15ag_19ne.reaction");
    const ReactionConfig reactionConfig = ReactionConfig::Load(reactionPath);

    // Overridable mass slit (MSLT) dispersive (X) half-aperture. The real
    // 2014-tune hardware card value is 0.75cm (see dat/dragon_2014_DSSSD.dat's
    // own 'RCOL' 'MSLT' card) -- MSLT is a physically adjustable slit on the
    // real hardware, not a fixed design constant, and this pilot's own
    // default is deliberately set wide open (2.5cm, the real dat/dragon_2001.dat
    // tune) at padsley's own request (2026-09-16), since MSLT-limited
    // transmission was masking other, more interesting effects (e.g. the
    // gamma-cascade recoil-angle cutoff). Change back to 0.75 (or override
    // via MSLT_HALFGAP_X_CM) to study the real 2014 tune specifically.
    const char* msltHalfGapEnv = std::getenv("MSLT_HALFGAP_X_CM");
    const double msltHalfGapXCm = msltHalfGapEnv ? std::atof(msltHalfGapEnv) : 2.5;

    // Same, for MSLT's non-dispersive (Y) half-aperture -- real 2014-tune
    // value 1.25cm, defaulted here to the same wide-open 2.5cm for the same
    // reason (see above); override via MSLT_HALFGAP_Y_CM.
    const char* msltHalfGapYEnv = std::getenv("MSLT_HALFGAP_Y_CM");
    const double msltHalfGapYCm = msltHalfGapYEnv ? std::atof(msltHalfGapYEnv) : 2.5;

    // Same diagnostic override, for the charge slit's (QSLT) own dispersive
    // (X) half-aperture -- real card value 1.25cm (see the "QSLT" call site
    // below). Defaults to the real value, so omitting it changes nothing.
    const char* qsltHalfGapEnv = std::getenv("QSLT_HALFGAP_X_CM");
    const double qsltHalfGapXCm = qsltHalfGapEnv ? std::atof(qsltHalfGapEnv) : 1.25;

    // Same override again, for QSLT's non-dispersive (Y) half-aperture --
    // real card value 1.25cm too (same call site). Defaults to the real
    // value, so omitting it changes nothing.
    const char* qsltHalfGapYEnv = std::getenv("QSLT_HALFGAP_Y_CM");
    const double qsltHalfGapYCm = qsltHalfGapYEnv ? std::atof(qsltHalfGapYEnv) : 1.25;

    // Target chamber + BGO array sit at the world origin -- the beamline's
    // own 'STRV' start (upstream of Q1), where the reaction actually
    // happens in the real machine. Gas species mirrors src/ugmate_trgt.f's
    // own rule (atarg<1.2 -> H2, else He) -- not hardcoded, so a (p,gamma)
    // reaction config (target=1H) gets a hydrogen target chamber instead
    // of the bundled 15O(alpha,gamma)19Ne config's helium one.
    const TargetChamber::TargetGas gas = reactionConfig.target.A < 1.2
                                              ? TargetChamber::TargetGas::kHydrogen
                                              : TargetChamber::TargetGas::kHelium;
    TargetChamber::Build(worldLV, gas);
    BgoArray::Build(worldLV);
    G4Material* copper = nist->FindOrBuildMaterial("G4_Cu");  // ugstmed.f medium 5

    // fMagneticRetuneScaleOverride (see DetectorConstruction.hh) lets a
    // caller inject an already-resolved scale -- main.cc's reaction-driven
    // entry points always do (SetUpRetunedChainGeometry() there measures
    // the real, energy-loss-degraded value rather than settling for
    // ComputeMagneticRetuneScale's own idealized fallback). <=0 (only
    // reachable if some other call site doesn't resolve one) falls back to
    // that same old behavior. electricScale is NOT simply this squared in
    // general -- see ComputeElectricRetuneScale's own comment.
    const double magneticScale = (fMagneticRetuneScaleOverride > 0.0)
                                      ? fMagneticRetuneScaleOverride
                                      : ComputeMagneticRetuneScale(reactionConfig);
    const double electricScale = ComputeElectricRetuneScale(reactionConfig, magneticScale);

    ChainState s;
    G4VisAttributes tstVis(G4Colour(0.0, 1.0, 1.0));  // cyan: MCP0/MCP1 markers
    tstVis.SetForceWireframe(true);

    // Diagnostic-only per-quad field-strength trims, on top of
    // magneticScale -- QuadTrimScale(n)/QN_TRIM_SCALE, see that function's
    // own comment. Added for probing whether QSLT's own (x|a)
    // achromatic-focus residual (padsley noticed QSLT isn't a clean focus
    // per Hutcheon Fig. 1/Table 2) can be zeroed the same way
    // kD2ResidualTrim above fixes D2's own bend-angle undershoot; then
    // generalized to every quad to check whether the same kind of residual
    // -- and the same kind of fix -- shows up at MSLT/FSLT too. All
    // default to 1.0 (no change from magneticScale alone).
    //
    // Findings, NOT baked in as new defaults (padsley's own call -- these
    // are much bigger corrections than kD2ResidualTrim's 2.6%, so they stay
    // opt-in pending more confidence): Q1_TRIM_SCALE has little effect on
    // QSLT's slope (0.73-0.86 mm/mrad over a 0.90-1.10 scan); Q2_TRIM_SCALE
    // =1.35 alone (Q1 untouched) drives the slope from 0.79 mm/mrad to
    // -0.0005 mm/mrad (confirmed unchanged at 4x finer FINE_STEP_CM, so a
    // real field effect, not integration noise), leaves MSLT's own
    // separate residual essentially unchanged (-2.27 to -2.25 mm/mrad --
    // this fix is properly localized to the Charge focus), and
    // independently improves real DSSSD transmission (500-event
    // o15ag_19ne: 18->108 hits; 1000-event k39pg_40ca: 553->612).

    const MitrayPoleData q1data = RetunedQuad(MitrayPoleData::Q1(), magneticScale * QuadTrimScale(1));
    const double q1EntryToCentreCm = q1data.A + (q1data.Z22 + q1data.L - q1data.Z11) / 2.0;
    s.xCm = 0.0;
    s.zCm = kQ1CenterZCm - q1EntryToCentreCm;  // Q1's own entry point
    s.thetaDeg = 0.0;

    // Sanity-check against the values published in the header/README.
    static_assert(kChainD1EntryZCm > 185.0 && kChainD1EntryZCm < 186.0, "check derivation");

    // The maxExtentCm/containerRadiusCm arguments below (ChainQuad's 5th,
    // ChainDipole/ChainEdipole's 6th) are precomputed offline from this
    // same chain's own fixed geometry: half the centre-to-centre distance
    // to whichever neighbour is closest (arc midpoint for a dipole/
    // e-dipole, see ChainDipole's own comment), so no two adjacent field
    // containers can overlap. See ChainQuad's comment for why this
    // matters (a silent, chain-wide field-masking bug otherwise).
    ChainQuad(s, worldLV, vacuum, "Q1", q1data, 27.5000);         // line 19
    Drift(s, 25.6925);                                   // DF7,  line 32
    ChainQuad(s, worldLV, vacuum, "Q2", RetunedQuad(MitrayPoleData::Q2(), magneticScale * QuadTrimScale(2)), 27.5000);  // line 35
    Drift(s, 26.4);                                       // DF9,  line 46
    ChainCollimator(s, worldLV, copper, "RC9", true, 0, 0, 7.46, 7.62, 26.4);  // line 47
    Drift(s, 26.4);                                       // DF10, line 50
    Drift(s, 11.0075);                                     // DFA1, line 52
    Shift(s, -0.19107);                                    // SH08, line 54
    ChainDipole(s, worldLV, vacuum, "D1", RetunedDipole(MitrayDipoleData::D1(), magneticScale), 58.2879);  // line 57

    Shift(s, 0.19107);                                     // SH09, line 70
    Drift(s, 0.0);                                         // DFA2, line 74
    Drift(s, 30.79);                                       // DF11, line 76
    ChainCollimator(s, worldLV, copper, "QSLT", false, 0, 0, qsltHalfGapXCm, qsltHalfGapYCm, 0.025, 20.0);  // line 77 (charge slit; QSLT_HALFGAP_X_CM/QSLT_HALFGAP_Y_CM override the real 1.25/1.25cm values -- see above)
    Drift(s, 3.9);                                         // DF12, line 80
    Drift(s, 22.65);                                       // DFAA, line 83
    ChainCollimator(s, worldLV, copper, "RC12", true, 0, 0, 4.92, 5.08, 26.55);  // line 84
    Drift(s, 26.55);                                       // DF13, line 87
    Drift(s, 18.32);                                       // DF14, line 90
    ChainQuad(s, worldLV, vacuum, "Q3", RetunedQuad(MitrayPoleData::Q3(), magneticScale * QuadTrimScale(3)), 21.1025);  // line 95

    Drift(s, 0.0);                                         // DF15, line 108
    Drift(s, 16.14);                                       // DF16, line 112
    ChainQuad(s, worldLV, vacuum, "Q4", RetunedQuad(MitrayPoleData::Q4(), magneticScale * QuadTrimScale(4)), 21.1025);  // line 115

    Drift(s, 0.0);                                         // DF17, line 128
    Drift(s, 21.62);                                       // DF18, line 130
    ChainQuad(s, worldLV, vacuum, "Q5", RetunedQuad(MitrayPoleData::Q5(), magneticScale * QuadTrimScale(5)), 27.5000);  // line 135

    Drift(s, 21.62);                                       // DF19, line 148
    ChainQuad(s, worldLV, vacuum, "Q6", RetunedQuad(MitrayPoleData::Q6(), magneticScale * QuadTrimScale(6)), 21.1025);  // line 153

    Drift(s, 0.0);                                         // DF20, line 166
    Drift(s, 16.14);                                       // DF21, line 168
    ChainQuad(s, worldLV, vacuum, "Q7", RetunedQuad(MitrayPoleData::Q7(), magneticScale * QuadTrimScale(7)), 21.1025);  // line 173

    Drift(s, 15.23);                                       // DF22, line 186
    Drift(s, 13.0);                                        // DF23, line 190
    ChainCollimator(s, worldLV, copper, "RC23", true, 0, 0, 7.46, 7.62, 13.0);  // line 191
    Drift(s, 13.0);                                        // DF24, line 194
    Drift(s, 14.97);                                       // DF25, line 198
    ChainCollimator(s, worldLV, copper, "RC25", true, 0, 0, 7.46, 7.62, 14.97);  // line 199
    Drift(s, 14.97);                                       // DF26, line 202
    Drift(s, 0.875);                                       // FD0,  line 204
    Shift(s, -0.0657);                                     // SH23, line 206
    ChainCollimator(s, worldLV, copper, "FC1", false, 0, 0, 5.0, 14.0, 0.875);  // line 207
    Drift(s, 8.875);                                       // FD1,  line 210
    ChainEdipole(s, worldLV, vacuum, "E1", RetunedEdipole(MitrayEdipoleData::E1(), electricScale), 49.8623);  // line 213

    Shift(s, 0.0657);                                      // SH24, line 222
    Drift(s, 8.875);                                       // FD2,  line 226
    ChainCollimator(s, worldLV, copper, "FC2", false, 0, 0, 5.0, 14.0, 0.875);  // line 227
    Drift(s, 0.875);                                       // FD3,  line 230
    Drift(s, 47.625);                                      // DF27, line 232
    ChainCollimator(s, worldLV, copper, "RC27", true, 0, 0, 4.92, 6.96, 47.625);  // line 233
    Drift(s, 47.625);                                      // DF28 (1st), line 236
    ChainCollimator(s, worldLV, copper, "MSLT", false, 0, 0, msltHalfGapXCm, msltHalfGapYCm, 0.025, 10.0);  // line 237 (mass slit; defaults wide open at 2.5x2.5cm -- see msltHalfGapXCm/msltHalfGapYCm above -- override via MSLT_HALFGAP_X_CM/MSLT_HALFGAP_Y_CM)
    Drift(s, 3.6);                                         // DF28 (2nd), line 243
    Drift(s, 22.9);                                        // DFBB, line 245
    ChainCollimator(s, worldLV, copper, "RC28", true, 0, 0, 4.92, 5.08, 26.8);  // line 246
    Drift(s, 26.8);                                        // DF29, line 249
    Drift(s, 13.5425);                                     // DF30, line 251
    ChainCollimator(s, worldLV, copper, "RC30", true, 0, 0, 4.92, 5.08, 13.5425);  // line 252
    Drift(s, 13.5425);                                     // DF31, line 255
    ChainQuad(s, worldLV, vacuum, "Q8", RetunedQuad(MitrayPoleData::Q8(), magneticScale * QuadTrimScale(8)), 27.5000);  // line 260

    Drift(s, 0.0);                                         // DF33, line 273
    Drift(s, 25.695);                                      // DF34, line 275
    ChainQuad(s, worldLV, vacuum, "Q9", RetunedQuad(MitrayPoleData::Q9(), magneticScale * QuadTrimScale(9)), 21.2250);  // line 280

    Drift(s, 15.81);                                       // DF35, line 293
    ChainQuad(s, worldLV, vacuum, "Q10", RetunedQuad(MitrayPoleData::Q10(), magneticScale * QuadTrimScale(10)), 21.2250);  // line 298

    Drift(s, 9.8);                                         // DF37, line 311
    Drift(s, 26.0);                                        // DF38, line 315
    Drift(s, 9.3);                                         // 'DF',  line 319
    ChainDipole(s, worldLV, vacuum, "D2", RetunedD2(MitrayDipoleData::D2(), magneticScale), 52.9418);  // line 324

    Drift(s, 56.076);                                      // DF39, line 341
    Drift(s, 6.025);                                       // DF40, line 346
    ChainCollimator(s, worldLV, copper, "RC40", true, 0, 0, 7.46, 7.62, 24.1);  // line 347
    Drift(s, 6.025);                                       // DF41, line 350
    Drift(s, 24.992);                                      // DF42, line 354
    ChainCollimator(s, worldLV, copper, "RC42", true, 0, 0, 7.46, 7.62, 24.992);  // line 355
    Drift(s, 24.992);                                      // DF43, line 358
    ChainQuad(s, worldLV, vacuum, "Q11", RetunedQuad(MitrayPoleData::Q11(), magneticScale * QuadTrimScale(11)), 21.2250);  // line 363

    Drift(s, 0.0);                                         // DF43 (2nd), line 376
    Drift(s, 15.81);                                       // DF44, line 378
    ChainQuad(s, worldLV, vacuum, "Q12", RetunedQuad(MitrayPoleData::Q12(), magneticScale * QuadTrimScale(12)), 21.2250);  // line 383

    Drift(s, 15.0);                                        // DF46, line 396
    Drift(s, 13.0);                                        // DF47, line 400
    ChainCollimator(s, worldLV, copper, "RC47", true, 0, 0, 7.46, 7.62, 13.0);  // line 401
    Drift(s, 13.0);                                        // DF48, line 404
    Drift(s, 17.175);                                      // DF49, line 408
    ChainCollimator(s, worldLV, copper, "RC49", true, 0, 0, 7.46, 7.62, 17.175);  // line 409
    Drift(s, 17.175);                                      // DF50, line 412
    Drift(s, 0.875);                                       // FD4,  line 414
    Shift(s, -0.089);                                      // SH43, line 416
    ChainCollimator(s, worldLV, copper, "FC3", false, 0, 0, 5.0, 15.0, 0.875);  // line 417
    Drift(s, 8.875);                                       // FD5,  line 420
    ChainEdipole(s, worldLV, vacuum, "E2", RetunedEdipole(MitrayEdipoleData::E2(), electricScale), 85.3127);  // line 423

    Shift(s, 0.089);                                       // SH44, line 432
    Drift(s, 8.875);                                       // FD6,  line 436
    ChainCollimator(s, worldLV, copper, "FC4", false, 0, 0, 5.0, 15.0, 0.875);  // line 437
    Drift(s, 0.875);                                       // FD7,  line 440
    Drift(s, 16.375);                                      // DF51, line 442
    ChainCollimator(s, worldLV, copper, "RC51", true, 0, 0, 7.46, 15.81, 16.375);  // line 443
    Drift(s, 16.375);                                      // DF52, line 446
    Drift(s, 19.775);                                      // DF53, line 450
    ChainCollimator(s, worldLV, copper, "RC53", true, 0, 0, 7.46, 7.62, 19.775);  // line 451
    Drift(s, 19.775);                                      // DF54, line 454
    Drift(s, 12.95);                                       // DF55, line 458
    ChainCollimator(s, worldLV, copper, "RC55", true, 0, 0, 7.46, 7.62, 12.95);  // line 459
    Drift(s, 12.95);                                       // DF56, line 462
    Drift(s, 12.0);                                        // DF51 (2nd), line 464
    ChainQuad(s, worldLV, vacuum, "Q13", RetunedQuad(MitrayPoleData::Q13(), magneticScale * QuadTrimScale(13)), 33.3000);  // line 469

    Drift(s, 19.9);                                        // DF57, line 482
    ChainQuad(s, worldLV, vacuum, "Q14", RetunedQuad(MitrayPoleData::Q14(), magneticScale * QuadTrimScale(14)), 33.3000);  // line 487

    // Past Q14: the beamline's final stretch to the focal-plane detector.
    // No more bending elements ('DIPO'/'EDIP' cards) appear, so theta is
    // fixed at its final value from here on. 'TEST' cards (TST3, TST4)
    // net to zero advance in ugeom_setup, same as before, but per this
    // file's own '15 Jan 2014' comment they mark where MCP0/MCP1 sit --
    // real MCP volumes ('MCPF' cards) never actually appear anywhere in
    // this file, so there is nothing to port there beyond the position
    // marker itself.
    Drift(s, 101.66);                                      // DF58, line 500
    ChainMarker(s, worldLV, vacuum, "TST3_MCP0", 5.0, 5.0, 0.05, tstVis);  // line 501
    ChainCollimator(s, worldLV, copper, "RC58", true, 0, 0, 7.46, 7.62, 58.83465);  // line 503
    Drift(s, 16.0);                                        // DF59, line 505
    ChainCollimator(s, worldLV, copper, "FSLT", false, 0, 0, 2.25, 2.25, 0.025, 10.0);  // line 507 (final slit)
    Drift(s, 43.0);                                        // DF60, line 509
    ChainMarker(s, worldLV, vacuum, "TST4_MCP1", 5.0, 5.0, 0.05, tstVis);  // line 511
    Drift(s, 3.98);                                        // FSB0, line 512
    Drift(s, 7.979);                                       // FSB1, line 514
    ChainCollimator(s, worldLV, copper, "FSB2", true, 0, 0, 4.76, 5.76, 7.979);  // line 516
    Drift(s, 7.979);                                       // FSB3, line 518
    Drift(s, 3.407);                                       // DB0X, line 520

    // 'COLL' 'DSSD' (line 522): the DSSSD focal-plane detector's aperture
    // position. Confirmed by direct search: no real silicon geometry is
    // built anywhere in the GEANT3 source for this ('ugeo_dssd' exists in
    // src/ugeom_mitray.f but is dead code, never called; no 'MCPF' card
    // appears in this file either), so unlike everything else in this
    // pilot this volume has no Fortran source to port at all. Real part
    // identified by the user from Micron Semiconductor Limited's 2024
    // catalogue: a **W1** double-sided silicon strip detector -- 16
    // junction (x) x 16 ohmic (y) strips, 3100um pitch, 49.50mm x 49.50mm
    // active area (catalogue p.33; three W1 packaging variants all share
    // this same active area). Thickness (1mm) and position (centred on
    // the beam axis, at the focal plane -- the running position the
    // file's own 'COLL' 'DSSD' card sits at, i.e. s.xCm/s.zCm here) per
    // the user's own direction, not the catalogue (which offers a wafer
    // thickness range, not one fixed value). Strip segmentation itself
    // isn't modeled -- a single silicon slab, not 16+16 separate strip
    // volumes/readout channels -- see DsssdHit.hh. It is now a real
    // G4VSensitiveDetector (DsssdSD), so steps through it do get recorded
    // as hits (position/time/particle/edep), even without per-strip
    // channel mapping.
    {
      G4NistManager* nist2 = G4NistManager::Instance();
      G4Material* silicon = nist2->FindOrBuildMaterial("G4_Si");
      auto* dsssdSolid = new G4Box("DSSSD", 2.475 * cm, 2.475 * cm, 0.05 * cm);  // 1mm thick
      auto* dsssdLV = new G4LogicalVolume(dsssdSolid, silicon, "DSSSD");
      G4VisAttributes dsssdVis(G4Colour(1.0, 0.0, 1.0));  // magenta
      dsssdLV->SetVisAttributes(dsssdVis);
      new G4PVPlacement(nullptr, G4ThreeVector(s.xCm * cm, 0.0, s.zCm * cm), dsssdLV, "DSSSD",
                         worldLV, false, 0, true);

      auto* dsssdSD = new DsssdSD("DsssdSD");
      G4SDManager::GetSDMpointer()->AddNewDetector(dsssdSD);
      dsssdLV->SetSensitiveDetector(dsssdSD);
    }

    ApplyFineStepIfRequested();
    return worldPV;
  }

  // Standalone Q1 (default).
  BuildQuad(worldLV, vacuum, "Q1", MitrayPoleData::Q1(), 0.0, kQ1CenterZCm);
  ApplyFineStepIfRequested();
  return worldPV;
}
