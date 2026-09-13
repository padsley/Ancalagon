#include "BgoArray.hh"

#include <cmath>
#include <vector>

#include "BgoSD.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Polyhedra.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"

namespace {

// dragon_2003.ffcards values, read by src/ugffgo.f into the COMMON block
// declared in inc/geometry.inc:
constexpr double kSFinger = 5.588;    // FSID(1): scintillator finger side [cm]
constexpr double kZFinger = 7.620;    // FSID(2): scintillator finger length [cm]
constexpr double kAirGap = 0.1270;    // FSID(3): air gap between hexagons [cm]
constexpr double kDAir1 = 0.0355;     // FSID(4): MGO gap/film, side [cm]
constexpr double kDAir2 = 0.3175;     // FSID(5): MGO gap/film, front [cm]
constexpr double kDMtl = 0.0635;      // FSID(6): metal sheet thickness [cm]
constexpr double kWall3 = 0.4978;     // WALL(3) [cm]
constexpr double kBoxWidth = 5.08;    // BGAP [cm]
constexpr double kAprt = 0.4496;      // HOLE [cm]
constexpr double kPmtLength = 2.5;    // PMTR(2) [cm]
constexpr double kPmtRadius = 2.54;   // PMTR(1) [cm] -- mtype_pmt=1 (circular), src/uvinit.f

// src/ugeom_gbox.f lines ~50-53 (ugeo_detector):
constexpr double kHexSmall = kSFinger + 2.0 * kDAir1 + 2.0 * kDMtl + kAirGap;  // hexagon_small_width
constexpr double kHexLarge = 2.0 * kHexSmall / 1.7320508075688772;  // hexagon_large_width (2/sqrt(3))
constexpr double kDepth = kZFinger + kDAir2 + kDMtl;                 // depth

// HSNG (outer housing) envelope in the assembly's own local frame,
// src/ugeom_gbox.f's ugeo_finger (~line 1137-1148): PGON apothem =
// hexagon_small_width/2, half-length (its Z, which becomes world X after
// the front/back rotation) = (d_mtl+d_air(2)+z_finger+pmt_length)/2. Every
// inner layer (FNGR/MGOR/SCNT) and the PMT are built from this same
// envelope, shrunk/offset by the real gap constants above -- see
// BgoArray::Build().
constexpr double kApothemCm = kHexSmall / 2.0;
constexpr double kHalfLengthCm = (kDMtl + kDAir2 + kZFinger + kPmtLength) / 2.0;

struct Placement {
  double xCm, yCm, zCm;
  bool front;  // true: irot_front (local Z axis -> world +X); false: irot_back (-> world -X)
};

// Literal transcription of src/ugeom_gbox.f's ugeo_finger (~lines 1230-1420):
// detectors 1-10 are individually hand-placed (their own x/y/z formulas, all
// using irot_front); detectors 11-30 come from the "Do k = 9,3,-1" / "Do j =
// 1,jm" loop (each (k,j) placing a mirrored back/front pair). Every position is
// placed at world (-x, y, -z) for the x/y/z computed here, matching the
// Fortran's own GSPOS calls exactly.
std::vector<Placement> ComputePlacements() {
  std::vector<Placement> p;
  p.reserve(30);

  // x is *not* reassigned in every one of the Fortran's ten detector blocks
  // (Fortran retains the last value across statements) -- reproduced here by
  // literally not reassigning xBase except where the source does.
  const double xBase = kDMtl / 2.0 + kDAir2 / 2.0 - kPmtLength / 2.0;
  double x, y, z;

  // Detector 1.
  x = xBase - 7.2;
  y = -kAprt - kWall3 - (1.0 / 2.0) * kHexLarge - 0.6;
  z = (5.0 / 2.0) * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 2.
  x = xBase;
  y = -kAprt - kWall3 - (5.0 / 4.0) * kHexLarge - 0.6;
  z = 2.0 * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 3.
  x = xBase - 6.7;
  y = kAprt + kWall3 + (1.0 / 2.0) * kHexLarge + 0.6;
  z = 2.0 * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 4.
  x = xBase;
  y = kAprt + kWall3 + (5.0 / 4.0) * kHexLarge + 0.6;
  z = (3.0 / 2.0) * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 5 (x unchanged from Detector 4, as in the Fortran).
  y = (9.0 / 8.0) * kHexLarge;
  z = kHexSmall / 2.0;
  p.push_back({-x, y, -z, true});

  // Detector 6 (x unchanged).
  y = (9.0 / 8.0) * kHexLarge;
  z = -kHexSmall / 2.0;
  p.push_back({-x, y, -z, true});

  // Detector 7 (x unchanged).
  y = kAprt + kWall3 + (5.0 / 4.0) * kHexLarge + 0.6;
  z = -(3.0 / 2.0) * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 8 (x unchanged).
  y = -kAprt - kWall3 - (5.0 / 4.0) * kHexLarge - 0.6;
  z = -2.0 * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 9.
  x = xBase;
  y = kAprt + kWall3 + (1.0 / 2.0) * kHexLarge + 0.6;
  z = -2.0 * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detector 10.
  x = xBase;
  y = -kAprt - kWall3 - (1.0 / 2.0) * kHexLarge - 0.6;
  z = -(5.0 / 2.0) * kHexSmall;
  p.push_back({-x, y, -z, true});

  // Detectors 11-30.
  for (int k = 9; k >= 3; --k) {
    z = -3.0 / 2.0 * kHexSmall + (k - 3) * kHexSmall / 2.0;
    const int jm = (k == 3 || k == 5 || k == 7 || k == 9) ? 1 : 2;

    for (int j = 1; j <= jm; ++j) {
      const double xAbs = kBoxWidth / 2.0 + (kDepth + kPmtLength) / 2.0;
      double yy = (jm == 1) ? -(3.0 / 8.0) * kHexLarge : -(9.0 / 8.0) * kHexLarge;
      yy += (j - 1) * (3.0 / 2.0) * kHexLarge;

      // "back" (irot_back): x = +xAbs -> placed at world x = -xAbs.
      p.push_back({-xAbs, yy, -z, false});
      // "front" (irot_front): x = -xAbs -> placed at world x = +xAbs.
      p.push_back({xAbs, yy, -z, true});
    }
  }

  return p;
}

}  // namespace

namespace {

// One Z-section (flat top/bottom, no taper) hexagonal PGON, apothem
// (centre-to-flat) `apothemCm` on both faces, half-length `halfLenCm`.
// phiStart=30deg (not 0): a GEANT4/GEANT3 polygon's phiStart is a *vertex*
// angle, but the real hexagons tile along the z_fngr stacking axis via
// FLAT faces (hexagon_small_width is a flat-to-flat measurement) -- a
// 30deg offset puts a flat, not a vertex, on that axis (same reasoning
// for every layer, since they all share the same hexagonal cross-section
// orientation, just shrinking apothem/half-length layer by layer).
G4VSolid* HexPrism(const char* name, double apothemCm, double halfLenCm) {
  const double zPlane[2] = {-halfLenCm * cm, halfLenCm * cm};
  const double rInner[2] = {0.0, 0.0};
  const double rOuter[2] = {apothemCm * cm, apothemCm * cm};
  return new G4Polyhedra(name, 30.0 * deg, 360.0 * deg, 6, 2, zPlane, rInner, rOuter);
}

}  // namespace

void BgoArray::Build(G4LogicalVolume* worldLV) {
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* bgo = nist->FindOrBuildMaterial("G4_BGO");
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  // MgO powder (medium 16 = 'MGO (POWDER)', src/ugmate.f's own material
  // #24: Mg+O 1:1 by atoms, density 1.87 g/cm3 -- a real, checked-in
  // number, much less than solid MgO's ~3.58 g/cm3, i.e. genuinely loose
  // powder packing, not a simplification on this pilot's part).
  auto* mgoPowder = new G4Material("MgOPowder", 1.87 * g / cm3, 2);
  mgoPowder->AddElement(nist->FindOrBuildElement("Mg"), 1);
  mgoPowder->AddElement(nist->FindOrBuildElement("O"), 1);
  // PMT window: src/ugeom_gbox.f's ugeo_pmt uses medium 17 ('GLASS'), but
  // src/ugmate.f's own material table (#25) gives 'GLASS' the exact same
  // a/z/density as material #17 ('SCINTILLATOR', plastic CH) -- a
  // stale-copy-paste artifact in the original file (the same kind already
  // flagged for CMBR's material index in TargetChamber.hh), not a real
  // glass composition. Ported as real NIST glass instead of reproducing
  // that bug.
  G4Material* glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");

  G4VisAttributes bgoVis(G4Colour(1.0, 0.85, 0.0));      // gold/amber: BGO crystal (SCNT)
  G4VisAttributes housingVis(G4Colour(1.0, 1.0, 1.0, 0.05));  // near-invisible: HSNG (air envelope)
  housingVis.SetForceWireframe(true);
  G4VisAttributes fingerVis(G4Colour(0.6, 0.6, 0.6, 0.5));    // grey: FNGR (aluminum can)
  G4VisAttributes reflectorVis(G4Colour(0.9, 0.9, 0.85, 0.4));  // off-white: MGOR (MgO powder)
  G4VisAttributes pmtVis(G4Colour(0.6, 0.8, 0.9, 0.4));       // pale blue: PMT (glass)

  // Real nested stack, src/ugeom_gbox.f's ugeo_finger (~lines 1142-1195)
  // and ugeo_pmt (~lines 1500-1519): HSNG (outer envelope, filled with
  // AIR -- there is no separate solid aluminum "can" wall in the real
  // file at all; FNGR itself, despite its name, is the aluminum housing)
  // > FNGR (aluminum) > MGOR (MgO powder reflector) > SCNT (BGO crystal,
  // the real reaction volume), plus PMT (glass) as HSNG's other direct
  // daughter, both sized/positioned by exactly the same shrink-by-gap
  // formulas the real GSVOLU/GSPOS calls use (each layer's own apothem
  // and half-length derived from the one before it, not independently).
  const double hsngApothemCm = kApothemCm;
  const double hsngHalfLenCm = kHalfLengthCm;
  auto* hsngSolid = HexPrism("HSNG", hsngApothemCm, hsngHalfLenCm);
  auto* hsngLV = new G4LogicalVolume(hsngSolid, air, "HSNG");
  hsngLV->SetVisAttributes(housingVis);

  const double fngrApothemCm = hsngApothemCm - kAirGap / 2.0;
  const double fngrHalfLenCm = hsngHalfLenCm - kPmtLength / 2.0;
  auto* fngrSolid = HexPrism("FNGR", fngrApothemCm, fngrHalfLenCm);
  auto* fngrLV = new G4LogicalVolume(fngrSolid, aluminum, "FNGR");
  fngrLV->SetVisAttributes(fingerVis);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -kPmtLength / 2.0 * cm), fngrLV, "FNGR", hsngLV,
                     false, 0, true);

  const double mgorApothemCm = fngrApothemCm - kDMtl;
  const double mgorHalfLenCm = fngrHalfLenCm - kDMtl / 2.0;
  auto* mgorSolid = HexPrism("MGOR", mgorApothemCm, mgorHalfLenCm);
  auto* mgorLV = new G4LogicalVolume(mgorSolid, mgoPowder, "MGOR");
  mgorLV->SetVisAttributes(reflectorVis);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kDMtl / 2.0 * cm), mgorLV, "MGOR", fngrLV, false, 0,
                     true);

  const double scntApothemCm = mgorApothemCm - kDAir1;
  const double scntHalfLenCm = mgorHalfLenCm - kDAir2 / 2.0;
  auto* scntSolid = HexPrism("SCNT", scntApothemCm, scntHalfLenCm);
  auto* scntLV = new G4LogicalVolume(scntSolid, bgo, "SCNT");
  scntLV->SetVisAttributes(bgoVis);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kDAir2 / 2.0 * cm), scntLV, "SCNT", mgorLV, false, 0,
                     true);

  // Sensitive detector: one BgoSD shared by all 30 crystal positions (all
  // 30 are the same SCNT logical volume, placed 30 times -- see
  // ComputePlacements()/the HSNG placement loop below); BgoSD reads the
  // real per-position ID (1-30) off the ancestor HSNG placement's own
  // copy number, not SCNT's own (always 0) -- see BgoHit.hh.
  auto* bgoSD = new BgoSD("BgoSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(bgoSD);
  scntLV->SetSensitiveDetector(bgoSD);

  // PMT: mtype_pmt==1 (src/uvinit.f), a circular G4Tubs (radius pmt_size,
  // half-length pmt_length/2), HSNG's other direct daughter -- placed at
  // HSNG-local z=depth/2 (depth = z_finger+d_air(2)+d_mtl, the same
  // constant kDepth already used for the crystal's own outer envelope),
  // which puts it flush against FNGR's own far end and HSNG's own front
  // face, with no overlap (verified: FNGR spans [-A/2, A/2-pmt_length] in
  // HSNG's frame, PMT spans exactly [A/2-pmt_length, A/2]).
  auto* pmtSolid = new G4Tubs("PMT", 0.0, kPmtRadius * cm, kPmtLength / 2.0 * cm, 0.0, 360.0 * deg);
  auto* pmtLV = new G4LogicalVolume(pmtSolid, glass, "PMT");
  pmtLV->SetVisAttributes(pmtVis);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, kDepth / 2.0 * cm), pmtLV, "PMT", hsngLV, false, 0,
                     true);

  const std::vector<Placement> placements = ComputePlacements();

  for (std::size_t i = 0; i < placements.size(); ++i) {
    const Placement& pl = placements[i];

    // The assembly's local Z axis (HSNG/FNGR/MGOR/SCNT's shared PGON
    // extrusion axis, and PMT's tube axis) must point along world +X
    // ("front") or -X ("back"). Active rotation by angle th about Y maps
    // local Z-hat to world (sin(th), 0, cos(th)); +90deg gives +X, -90deg
    // gives -X. G4PVPlacement's rotation matrix is the *inverse* of that
    // active rotation (a standard, well-documented Geant4 convention), so
    // the matrix passed here uses the opposite sign of the angle actually
    // wanted.
    auto* rot = new G4RotationMatrix();
    rot->rotateY(pl.front ? -90.0 * deg : 90.0 * deg);

    new G4PVPlacement(rot, G4ThreeVector(pl.xCm * cm, pl.yCm * cm, pl.zCm * cm), hsngLV, "HSNG",
                       worldLV, false, static_cast<int>(i) + 1, true);
  }
}
