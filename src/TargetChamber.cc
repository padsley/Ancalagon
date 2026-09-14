#include "TargetChamber.hh"

#include <utility>
#include <vector>

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4VisAttributes.hh"

namespace {

// --- Gas materials (src/ugmate_trgt.f) -------------------------------------
//
// mtarg (index 27): the full target-pressure gas, ptarg=5/760 atm. The
// real code picks H2 or He by projectile mass (ugmate_trgt.f: atarg<1.2 ->
// H2, else He) -- TargetGas (TargetChamber.hh) is that same switch, ported
// as an explicit enum rather than re-deriving it from atarg here (the
// caller already knows which reaction, and hence which gas, it wants --
// see DetectorConstruction.cc and ReactionKinematics.hh/README's "Reaction
// specification"; this configuration's own reaction, 15O(alpha,gamma)19Ne,
// needs He).
// mcent (index 28): the "central vacuum" gas that fills every EN/EX
// pumping-chain volume, at ptarg*0.002 -- 1/500th of the target pressure,
// same gas species as mtarg.
// Both are approximated the same way: NIST's STP H2/He density scaled by
// the pressure fraction (the file's own reference densities, 8.38e-5 g/cm3
// for H2 and 1.6586e-4 g/cm3 for He, are close but not identical figures
// to the NIST STP ones -- not chased here).
constexpr double kTargetPressureTorr = 5.0;
constexpr double kCentralVacuumFactor = 0.002;

const char* NistNameFor(TargetChamber::TargetGas gas) {
  return gas == TargetChamber::TargetGas::kHydrogen ? "G4_H" : "G4_He";
}

// --- src/ugeom_defin.f (tubetype=0 branch) ---------------------------------
constexpr double kRrms = 6.0;    // outer radius of every EN/EX cone & tube
constexpr double kLen1 = 7.6;
constexpr double kLen2 = 8.7;
constexpr double kLen3 = 10.0;
constexpr double kRilen2 = 0.59, kRiren2 = 0.495;   // EN2C bore radii
constexpr double kRilen3 = 0.75, kRiren3 = 0.65;    // EN3C bore radii
constexpr double kZent1 = -15.6, kZent3 = -36.05, kZent5 = -63.75;
constexpr double kTLrms = 88.0;

// --- src/ugeom_gbox.f's ugeo_detector (box/cell geometry constants) -------
constexpr double kBoxLength = 17.069;
constexpr double kWall2 = 0.3175;       // WALL(2), dragon_2003.ffcards
constexpr double kColLength = 15.24;
constexpr double kColCollarLength = 0.0;  // never assigned in this file's active path
constexpr double kBoxHeight = 20.0;

// CELL's y-offset inside CMBG (ugeo_detector, targtype==0 branch): CELL's
// real position relative to its own box, kept as derived from the original
// GEANT3 model (CELL's own local half-length along this axis, after
// CellRotation(), is 4.208 cm -- moving CELL itself further from CMBG's
// center than this would push it into CMBG's wall, see CMBG's own
// half-height below).
constexpr double kCellYInCmbgCm = 0.5 * kBoxHeight - 4.208 - 0.976;

// CMBR's own y-offset in DETE (ugeo_detector's own DETE placement formula).
// Combined with kCellYInCmbgCm (4.816 cm), this puts CELL's own *geometric
// center* at world y = -2.013 cm -- which looks like it's off the optics'
// y=0 beam axis (see RunAction.hh/main.cc's own older comments on this).
// It isn't a bug: CELL's own center is not the same as the actual
// beam/recoil path through it, which runs through EAPG/XAPG (its own
// entrance/exit apertures, off-center by kApertureLocalZCm -- see below).
// This -6.829 value is tuned, together with kCellYInCmbgCm and
// kApertureLocalZCm, to put THAT path (BeamApertureWorldYCm(), below) on
// the y=0 axis -- UHOL/DHOL (placed at local y = -kBeamHeightCm inside
// CMBR, below) confirm this: they stay pinned to world y=0 regardless of
// this constant's value, and were already on y=0 before anything in this
// file was touched, meaning CELL's center being off-axis was always the
// intended, self-consistent design, not a simplification to fix.
constexpr double kBeamHeightCm = -kBoxHeight / 2.0 + 3.171;

// EAPG/XAPG (CellAndApertures(), below) -- the gas cell's actual entrance/
// exit apertures -- sit at local z = kApertureLocalZCm, not CELL's own
// local z=0 (its geometric center). A primary vertex placed at CELL's
// center (as this pilot's PrimaryGeneratorAction previously was) fires a
// recoil almost straight along local x (world z, see CellRotation()) but
// 2.008 cm below the actual apertures -- it runs into CELL's solid wall
// after ~3.9 cm instead of ever reaching them. See BeamApertureWorldYCm().
constexpr double kApertureLocalZCm = 2.008;

// EAPG/XAPG's own local x=+-kApertureLocalXCm -- i.e. half the total beam
// path length through the gas (EAPG center to XAPG center); see
// BeamPathHalfLengthCm() below.
constexpr double kApertureLocalXCm = 5.315;

// Every gas-filled volume in CellAndApertures() (CELG, EAPG, XAPG) is
// built from BuildGas() with this exact name -- exposed via
// TargetGasMaterialName() so code outside this file (PrimaryGeneratorAction's
// beam-energy-loss depth sampling, see GasStoppingPower) can look the same
// G4Material up post-construction (G4Material::GetMaterial(name)) instead
// of building a second, merely similar one.
constexpr const char* kTargetGasMaterialName = "TargetGas";

G4Material* BuildGas(const char* name, TargetChamber::TargetGas gas, double pressureFractionOfAtm) {
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* stpReference = nist->FindOrBuildMaterial(NistNameFor(gas));
  const double density = stpReference->GetDensity() * pressureFractionOfAtm;
  auto* material = new G4Material(name, density, 1);
  material->AddMaterial(stpReference, 1.0);
  return material;
}

// GEANT3 GSROTM's the1/phi1,the2/phi2,the3/phi3 give the daughter's own
// local X/Y/Z axes as direction cosines (sin(the)cos(phi), sin(the)sin(phi),
// cos(the)) *in the mother's frame* -- i.e. this already is the active
// rotation of the daughter. G4PVPlacement wants the inverse of that (a
// standard, well-documented Geant4 convention, already used the same way in
// BgoArray.cc); since every rotation used here is a pure axis-permutation
// (all direction cosines are 0 or +-1), the inverse is just its transpose.
G4RotationMatrix* CellRotation() {
  // src/ugeom_gbox.f's irot_box: the1=180,phi1=0, the2=270,phi2=180,
  // the3=90,phi3=90 -> local X->mother -Z, local Y->mother +X, local
  // Z->mother +Y. Taken literally this active local->mother rotation A
  // (columns (0,0,-1),(1,0,0),(0,1,0)) has determinant -1 -- an improper
  // transform (reflection), which G4RotationMatrix (unlike GEANT3's own
  // GSROTM, which tolerates it) refuses to represent. CELL is a TRD1 with
  // its DY half-width constant on both faces (shape(3)), i.e. physically
  // symmetric under local Y -> -Y, and every one of its daughters (CELG,
  // EAPG, XAPG) sits at local y=0 -- so flipping which world axis "local
  // Y" maps to changes nothing about the resulting physical placement.
  // Using -A's local-Y column instead gives an equivalent, now-proper
  // (det=+1) rotation A': columns (0,0,-1),(-1,0,0),(0,1,0).
  //
  // CLHEP's rotateAxes(newX,newY,newZ) builds a matrix whose *columns* are
  // (newX,newY,newZ) -- i.e. it builds the active rotation directly, not
  // its inverse -- so to hand G4PVPlacement the inverse of A' (the
  // convention it wants, same as BgoArray.cc), pass A''s rows: that makes
  // the built matrix's columns equal A''s rows, i.e. A'^T = A''s inverse
  // (A' is orthogonal).
  auto* rot = new G4RotationMatrix();
  rot->rotateAxes(G4ThreeVector(0, -1, 0), G4ThreeVector(0, 0, 1), G4ThreeVector(-1, 0, 0));
  return rot;
}

G4RotationMatrix* EapgXapgRotation() {
  // src/ugeom_gbox.f's irot_col: the1=180,phi1=0, the2=90,phi2=90,
  // the3=90,phi3=0 -> tube's local X->CELL -Z, local Y->CELL +Y, local Z
  // (its cylinder axis)->CELL +X. This active rotation (columns
  // (0,0,-1),(0,1,0),(1,0,0)) is already proper (det=+1) as given -- no
  // adjustment needed. EAPG/XAPG are placed as CELL's own daughters (see
  // CellAndApertures()), so this is relative to CELL's local frame only;
  // it is not composed with CellRotation() -- Geant4's mother/daughter
  // hierarchy does that composition automatically, exactly as GEANT3's
  // own volume tree does. Rows passed to rotateAxes for the same reason
  // as CellRotation() above.
  auto* rot = new G4RotationMatrix();
  rot->rotateAxes(G4ThreeVector(0, 0, 1), G4ThreeVector(0, 1, 0), G4ThreeVector(-1, 0, 0));
  return rot;
}

// The real windowless gas cell + its beam apertures, src/ugeom_gbox.f
// lines ~1090-1130. CELL is built here as the mother of CELG/EAPG/XAPG so
// their own GEANT3 local coordinates (relative to CELL) can be transcribed
// directly, with no extra coordinate-composition math needed. `cmbgLV` is
// CMBG's own logical volume (see OuterBoxAndShielding) -- CELL is really
// its daughter (src/ugeom_gbox.f: `gspos('CELL',1,'CMBG',...)`), not a
// direct child of the world; nesting it there instead of flattening (as
// an earlier version of this pilot did, before CMBR/CMBG existed at all)
// avoids CELL and CMBR's own steel shell occupying the same space.
void CellAndApertures(G4LogicalVolume* cmbgLV, TargetChamber::TargetGas gas) {
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  G4Material* targetGas = BuildGas(kTargetGasMaterialName, gas, kTargetPressureTorr / 760.0);

  // CELL: aluminum TRD1 outer shell. shape = (2.115, 6.759, 1.905, 4.208).
  auto* cellSolid = new G4Trd("CELL", 2.115 * cm, 6.759 * cm, 1.905 * cm, 1.905 * cm, 4.208 * cm);
  auto* cellLV = new G4LogicalVolume(cellSolid, aluminum, "CELL");
  cellLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.5)));  // translucent grey

  // CELG: mtarg-gas TRD1 inner shell, the real reaction volume. shape =
  // (1.798, 6.094, 1.588, 3.891), centred on CELL with no rotation.
  auto* celgSolid = new G4Trd("CELG", 1.798 * cm, 6.094 * cm, 1.588 * cm, 1.588 * cm, 3.891 * cm);
  auto* celgLV = new G4LogicalVolume(celgSolid, targetGas, "CELG");
  celgLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.4)));  // pale blue
  new G4PVPlacement(nullptr, G4ThreeVector(), celgLV, "CELG", cellLV, false, 0, true);

  // EAPG/XAPG: mtarg-gas TUBE beam apertures bored through CELL's slanted
  // side wall (tubetype==0 radii: EAPG rmax=0.3cm, XAPG rmax=0.4cm; both
  // half-length 0.5cm), at local x=+-kApertureLocalXCm, z=kApertureLocalZCm,
  // y=0 -- i.e. the real beam/recoil path through this cell is the line
  // (x varies, y=0, z=kApertureLocalZCm), not CELL's own local z=0 (see
  // BeamApertureWorldYCm() below).
  G4RotationMatrix* colRot = EapgXapgRotation();
  auto* eapgSolid = new G4Tubs("EAPG", 0.0, 0.3 * cm, 0.5 * cm, 0.0, 360.0 * deg);
  auto* eapgLV = new G4LogicalVolume(eapgSolid, targetGas, "EAPG");
  eapgLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.4)));
  new G4PVPlacement(colRot, G4ThreeVector(kApertureLocalXCm * cm, 0.0, kApertureLocalZCm * cm),
                     eapgLV, "EAPG", cellLV, false, 0, true);

  auto* xapgSolid = new G4Tubs("XAPG", 0.0, 0.4 * cm, 0.5 * cm, 0.0, 360.0 * deg);
  auto* xapgLV = new G4LogicalVolume(xapgSolid, targetGas, "XAPG");
  xapgLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.4)));
  new G4PVPlacement(colRot, G4ThreeVector(-kApertureLocalXCm * cm, 0.0, kApertureLocalZCm * cm),
                     xapgLV, "XAPG", cellLV, false, 0, true);

  // CELL's placement inside CMBG (see kCellYInCmbgCm above).
  new G4PVPlacement(CellRotation(), G4ThreeVector(0.0, kCellYInCmbgCm * cm, 0.0), cellLV, "CELL",
                     cmbgLV, false, 0, true);
}

// One EN/EX gas TUBE, optionally with a CONE collimator insert centred
// inside it (nullptr solid -> no insert, matching the file's own volumes
// that carry no collimator). Everything here uses material mcent.
void GasTube(G4LogicalVolume* worldLV, G4Material* mcent, const char* name, double halfLenCm,
             double zCm, G4VSolid* insertSolid, const char* insertName) {
  auto* solid = new G4Tubs(name, 0.0, kRrms * cm, halfLenCm * cm, 0.0, 360.0 * deg);
  auto* lv = new G4LogicalVolume(solid, mcent, name);
  lv->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.15)));  // faint pale blue

  if (insertSolid != nullptr) {
    G4NistManager* nist = G4NistManager::Instance();
    auto* insertLV = new G4LogicalVolume(insertSolid, nist->FindOrBuildMaterial("G4_Al"), insertName);
    insertLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.6)));
    new G4PVPlacement(nullptr, G4ThreeVector(), insertLV, insertName, lv, false, 0, true);
  }

  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, zCm * cm), lv, name, worldLV, false, 0, true);
}

// The differential-pumping chain, src/ugeom_trgt_small.f. All z positions
// and half-lengths below are that file's own formulas, evaluated (its
// GEANT3 REALs, not re-derived): entrance side numbered from the target
// outward (EN1G..EN6G), exit side similarly (EX2G..EX10, EX6C).
void DifferentialPumpingChain(G4LogicalVolume* worldLV, TargetChamber::TargetGas gas) {
  G4Material* mcent =
      BuildGas("CentralVacuumGas", gas, (kTargetPressureTorr / 760.0) * kCentralVacuumFactor);
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");

  // EN1G: no collimator insert. zcol per line 89-90.
  const double en1gHalfLen = 0.5 * 0.68;
  const double zcol = -(0.5 * kBoxLength - kWall2 - kColCollarLength + kColLength + en1gHalfLen);
  GasTube(worldLV, mcent, "EN1G", en1gHalfLen, zcol, nullptr, nullptr);

  // EN2G: no collimator insert.
  const double en2gHalfLen = 0.5 * (kZent1 - kZent3 - kLen1 - kLen2);
  const double en2gZ = 0.5 * (kZent1 + kZent3 - kLen1 + kLen2);
  GasTube(worldLV, mcent, "EN2G", en2gHalfLen, en2gZ, nullptr, nullptr);

  // EN3G: EN2C aluminum collimator insert (exact fit, half-length = len2).
  auto* en2c = new G4Cons("EN2C", kRilen2 * cm, kRrms * cm, kRiren2 * cm, kRrms * cm, kLen2 * cm,
                           0.0, 360.0 * deg);
  GasTube(worldLV, mcent, "EN3G", kLen2, kZent3, en2c, "EN2C");

  // EN4G: no collimator insert.
  const double en4gHalfLen = 0.5 * (kZent3 - kZent5 - kLen3 - kLen2);
  const double en4gZ = 0.5 * (kZent3 + kZent5 + kLen3 - kLen2);
  GasTube(worldLV, mcent, "EN4G", en4gHalfLen, en4gZ, nullptr, nullptr);

  // EN5G: EN3C aluminum collimator insert (exact fit, half-length = len3).
  auto* en3c = new G4Cons("EN3C", kRilen3 * cm, kRrms * cm, kRiren3 * cm, kRrms * cm, kLen3 * cm,
                           0.0, 360.0 * deg);
  GasTube(worldLV, mcent, "EN5G", kLen3, kZent5, en3c, "EN3C");

  // EN6G: no collimator insert.
  const double en6gHalfLen = 0.5 * (kTLrms + kZent5 - kLen3);
  const double en6gZ = -kTLrms + en6gHalfLen;
  GasTube(worldLV, mcent, "EN6G", en6gHalfLen, en6gZ, nullptr, nullptr);

  // EX2G: no collimator insert.
  GasTube(worldLV, mcent, "EX2G", 6.7 / 2.0, 26.25, nullptr, nullptr);

  // EX3G: EX2C insert.
  auto* ex2c = new G4Cons("EX2C", (1.81 / 2.0) * cm, kRrms * cm, (1.81 / 2.0) * cm, kRrms * cm,
                           (3.52 / 2.0) * cm, 0.0, 360.0 * deg);
  GasTube(worldLV, mcent, "EX3G", 3.52 / 2.0, 31.44, ex2c, "EX2C");

  // EX4G: no collimator insert.
  GasTube(worldLV, mcent, "EX4G", 4.2 / 2.0, 35.3, nullptr, nullptr);

  // EX5G: EX3C insert.
  auto* ex3c = new G4Cons("EX3C", (1.8 / 2.0) * cm, kRrms * cm, (2.01 / 2.0) * cm, kRrms * cm,
                           (5.5 / 2.0) * cm, 0.0, 360.0 * deg);
  GasTube(worldLV, mcent, "EX5G", 5.5 / 2.0, 40.15, ex3c, "EX3C");

  // EX6G: no collimator insert.
  GasTube(worldLV, mcent, "EX6G", 5.24 / 2.0, 45.48, nullptr, nullptr);

  // EX7G: EX4C insert.
  auto* ex4c = new G4Cons("EX4C", (2.67 / 2.0) * cm, kRrms * cm, (2.67 / 2.0) * cm, kRrms * cm,
                           (4.16 / 2.0) * cm, 0.0, 360.0 * deg);
  GasTube(worldLV, mcent, "EX7G", 4.16 / 2.0, 50.22, ex4c, "EX4C");

  // EX8G: no collimator insert.
  GasTube(worldLV, mcent, "EX8G", 7.85 / 2.0, 56.175, nullptr, nullptr);

  // EX9G: EX5C insert (src/ugeom_trgt_small.f places this with 'MANY', but
  // EX5C's own half-length equals EX9G's -- an exact-fit coincident insert
  // like all the others above, not a genuine overlap; 'ONLY' semantics are
  // fine here).
  auto* ex5c = new G4Cons("EX5C", (2.71 / 2.0) * cm, kRrms * cm, (3.67 / 2.0) * cm, kRrms * cm,
                           (24.9 / 2.0) * cm, 0.0, 360.0 * deg);
  GasTube(worldLV, mcent, "EX9G", 24.9 / 2.0, 72.55, ex5c, "EX5C");

  // EX10: no collimator insert.
  GasTube(worldLV, mcent, "EX10", 3.0 / 2.0, 86.5, nullptr, nullptr);

  // EX6C: placed directly in world with no surrounding gas (src/
  // ugeom_trgt_small.f's own comment: "Since the following collimator
  // extends into Q1, it is placed in the WRLD coordinates with no gas").
  // Real full extent (half-length 22.9/2=11.45cm centred at z=106.05, i.e.
  // z in [94.6, 117.5]) would overlap Q1's own tube (Q1's entry face sits
  // at the beamline's real STRV->Q1 distance, 106.885cm -- see
  // DetectorConstruction.hh's kQ1CenterZCm -- and Q1's radius, 5.3975cm,
  // is smaller than EX6C's outer radius, kRrms=6.0cm, in the z range they'd
  // share). Clipped here to end exactly at Q1's entry face instead of
  // running through it; the trimmed ~5cm sliver is thin aluminum shell at
  // r>~3cm, well outside the beam's own envelope, so this doesn't touch
  // anything trajectory-relevant. rmin at the new end face is linearly
  // interpolated along the same taper the file itself uses (2.725cm at
  // z=94.6 to 3.255cm at z=117.5).
  constexpr double kQ1EntryZCm = 106.885;
  constexpr double kEx6cFullHalfLen = 22.9 / 2.0;         // 11.45
  constexpr double kEx6cFullCentreZ = 106.05;
  constexpr double kEx6cFullStartZ = kEx6cFullCentreZ - kEx6cFullHalfLen;  // 94.6
  constexpr double kEx6cFullEndZ = kEx6cFullCentreZ + kEx6cFullHalfLen;    // 117.5
  constexpr double kEx6cRminStart = 5.45 / 2.0;  // 2.725
  constexpr double kEx6cRminEnd = 6.51 / 2.0;    // 3.255

  const double ex6cHalfLen = (kQ1EntryZCm - kEx6cFullStartZ) / 2.0;
  const double ex6cCentreZ = (kEx6cFullStartZ + kQ1EntryZCm) / 2.0;
  const double clipFraction = (kQ1EntryZCm - kEx6cFullStartZ) / (kEx6cFullEndZ - kEx6cFullStartZ);
  const double ex6cRminAtEntry = kEx6cRminStart + (kEx6cRminEnd - kEx6cRminStart) * clipFraction;

  auto* ex6c = new G4Cons("EX6C", kEx6cRminStart * cm, kRrms * cm, ex6cRminAtEntry * cm, kRrms * cm,
                           ex6cHalfLen * cm, 0.0, 360.0 * deg);
  auto* ex6cLV = new G4LogicalVolume(ex6c, aluminum, "EX6C");
  ex6cLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.6)));
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, ex6cCentreZ * cm), ex6cLV, "EX6C", worldLV,
                     false, 0, true);
}

}  // namespace

// --- CMBR/CMBG outer box + UHOL/PUAI..PUJ2/DHOL/PDAI..PDI4 shielding -------
//
// src/ugeom_gbox.f's ugeo_detector, ~lines 71-732 (tubetype==0 branch,
// what dat/dragon_2014_DSSSD.dat selects; the tubetype==1 'ELSEIF' branch
// from line 734 on is a different target's own hardware, not ported).
// 'DETE' (the mother of this whole region, and of CELL/EN-EX above) is an
// invisible bookkeeping box (GEANT3's own 'SEEN',0) placed in 'WRLD' with
// zero offset -- exactly like CellAndApertures()/DifferentialPumpingChain()
// above, this flattens DETE away and places everything directly in
// worldLV, at the same absolute coordinates DETE's own zero offset would
// give anyway.
//
// CMBR (aluminum -- no: medium 20 is 'STAINLESS STEEL', a real, standard
// GEANT3 material, not the custom material table; TargetChamber.hh's
// previous note calling this "NAI:TL" was a misread of the *medium*
// index against the wrong table -- ugstmed.f's medium 20 maps to
// *material* 26, a standalone standard material, and name_med(20) itself
// already reads 'STAINLESS STEEL') is the outer structural box; CMBG is
// the gas-filled interior cavity (GEANT3 places it with 'MANY' -- this
// port uses an ordinary daughter placement instead, which gives the same
// physical result, a hollow box with a wall(2)-thick shell, without
// GEANT3-specific volume-overlap semantics GEANT4 doesn't have). UHOL/
// DHOL are the beam holes bored straight through that shell (entrance/
// exit); PUAI/PDAI (+ PUBI/PDBI) are the aluminum collimator inserts
// just inside them, in CMBG.
//
// PUA1..PUJ2 (upstream, i.e. entrance side, more negative z) and
// PDA1..PDI4 (downstream/exit side; note there is no separate 'C'
// section and no 'J' section on this side -- the real file combines
// what would be C+D and F+G into one section each on this side, and
// omits J outright, all ported exactly as asymmetric, not "fixed" to
// match the upstream side) are the lead/aluminum/gas collar chain
// bolted onto each face of the box, each lettered section a handful of
// concentric tubes (one, PUG, with a lead BOX instead of a tube as its
// outer layer) sharing one half-length, successively narrower toward
// the beam axis -- see CollarSection/BuildCollarChain below.
//
// SCOPE/CAVEAT: this and the EN/EX differential-pumping chain above
// port two different GEANT3 subroutines (ugeom_gbox.f here vs.
// src/ugeom_trgt_small.f there) that were evidently authored somewhat
// independently and can occupy overlapping z-ranges around the same
// physical region (the EN/EX tubes' own kRrms=6cm envelope is wider
// than most, but not all, of this chain's own layers -- PUG1's lead box
// alone is 15.5x23.5cm). Both are ported exactly as their own source
// gives them; where G4's overlap checker flags the two chains crossing,
// that reflects the real files' own, apparently uncoordinated overlap
// (GEANT3's tracker tolerates arbitrary volume overlaps the way GEANT4's
// strict mother/daughter model does not), not a placement error
// introduced by this port.
namespace {

// One lettered PUx/PDx "collar" section: `layersOuterToInner.size()`
// concentric G4Tubs (each the next's mother -- mirrors how the real
// GSVOLU/GSPOS calls nest PUx2 inside PUx1, PUx3 inside PUx2, etc.), all
// sharing `halfLenCm`, the outermost placed in `motherLV` at local
// z=`zCm`. `namePrefix` + "1", "2", ... names each layer, matching the
// real volume names (PUB1/PUB2/PUB3, etc.) closely enough to identify in
// the geometry tree.
void CollarTubes(G4LogicalVolume* motherLV, const char* namePrefix, double zCm, double halfLenCm,
                  const std::vector<std::pair<G4Material*, double>>& layersOuterToInner,
                  const G4VisAttributes& vis) {
  G4LogicalVolume* parent = motherLV;
  for (std::size_t i = 0; i < layersOuterToInner.size(); ++i) {
    const G4String name = G4String(namePrefix) + std::to_string(i + 1);
    auto* solid = new G4Tubs(name, 0.0, layersOuterToInner[i].second * cm, halfLenCm * cm, 0.0,
                              360.0 * deg);
    auto* lv = new G4LogicalVolume(solid, layersOuterToInner[i].first, name);
    lv->SetVisAttributes(vis);
    const G4ThreeVector pos = (i == 0) ? G4ThreeVector(0.0, 0.0, zCm * cm) : G4ThreeVector();
    new G4PVPlacement(nullptr, pos, lv, name, parent, false, 0, true);
    parent = lv;
  }
}

struct CollarSection {
  const char* name;
  double halfLenCm;
  std::vector<std::pair<G4Material*, double>> layersOuterToInner;  // (material, radiusCm)
};

// Places `sections` back-to-back along +-z starting at `startZCm` (each
// section's own z = signZ*(startZCm + running + halfLen), running then
// advances by 2*halfLen) -- mirrors the real file's own cumulative
// "z = box_length/2 + 2*bp<earlier>_len + ... + bp<this>_len" formulas
// exactly, without hardcoding each one by hand.
void BuildCollarChain(G4LogicalVolume* worldLV, double startZCm, double signZ,
                      const std::vector<CollarSection>& sections, const G4VisAttributes& vis) {
  double runningCm = 0.0;
  for (const auto& sec : sections) {
    const double zCm = signZ * (startZCm + runningCm + sec.halfLenCm);
    CollarTubes(worldLV, sec.name, zCm, sec.halfLenCm, sec.layersOuterToInner, vis);
    runningCm += 2.0 * sec.halfLenCm;
  }
}

// Builds CMBR/CMBG and the whole UHOL..PDI4 shielding chain, returning
// CMBG's own logical volume so CellAndApertures() can nest CELL inside
// it (matching the real file's own hierarchy -- see that function's
// header comment).
G4LogicalVolume* OuterBoxAndShielding(G4LogicalVolume* worldLV, TargetChamber::TargetGas gas) {
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  const double ptargAtm = kTargetPressureTorr / 760.0;
  // ment(1..3)/mex(1..3): src/ugmate_trgt.f's own pressure fractions
  // (0.05/0.036/0.016 x ptarg) for both the entrance and exit collar
  // chains -- identical fractions on both sides, so one gas per fraction
  // is reused for both (see header comment; no physical difference).
  G4Material* gas1 = BuildGas("CollarGas1", gas, ptargAtm * 0.05);
  G4Material* gas2 = BuildGas("CollarGas2", gas, ptargAtm * 0.036);
  G4Material* gas3 = BuildGas("CollarGas3", gas, ptargAtm * 0.016);
  // mbox: CMBG's own interior fill, src/ugmate_trgt.f's 0.056 x ptarg.
  G4Material* boxGas = BuildGas("BoxGas", gas, ptargAtm * 0.056);

  G4VisAttributes steelVis(G4Colour(0.55, 0.55, 0.6, 0.3));
  G4VisAttributes leadVis(G4Colour(0.35, 0.35, 0.4, 0.6));
  G4VisAttributes alVis(G4Colour(0.7, 0.7, 0.7, 0.4));

  // CMBR (steel outer box) + CMBG (gas interior, ordinary daughter -- see
  // header comment on the real file's own 'MANY' there).
  constexpr double kBoxWidthCm = 5.08;  // box_width (BGAP ffcard) -- same global as BgoArray.cc's kBoxWidth
  const double boxHalfXCm = kBoxWidthCm / 2.0;
  auto* cmbrSolid = new G4Box("CMBR", boxHalfXCm * cm, kBoxHeight / 2.0 * cm, kBoxLength / 2.0 * cm);
  auto* cmbrLV = new G4LogicalVolume(cmbrSolid, steel, "CMBR");
  cmbrLV->SetVisAttributes(steelVis);

  auto* cmbgSolid = new G4Box("CMBG", (boxHalfXCm - kWall2) * cm, (kBoxHeight / 2.0 - kWall2 / 2.0) * cm,
                               (kBoxLength / 2.0 - kWall2) * cm);
  auto* cmbgLV = new G4LogicalVolume(cmbgSolid, boxGas, "CMBG");
  cmbgLV->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.05)));
  new G4PVPlacement(nullptr, G4ThreeVector(), cmbgLV, "CMBG", cmbrLV, false, 0, true);

  new G4PVPlacement(nullptr, G4ThreeVector(0.0, kBeamHeightCm * cm, 0.0), cmbrLV, "CMBR", worldLV,
                     false, 0, true);

  // UHOL/DHOL: beam holes through the steel wall itself (CMBR's own
  // daughter, not CMBG's -- see header comment on why these two don't
  // overlap). PUAI/PDAI (+PUBI/PDBI): the aluminum collimator insert
  // just inside, in CMBG.
  constexpr double kBpaInLenCm = 0.472;
  auto uhol = [&](const char* name, G4Material* holeGas, double radiusCm, double zLocalCm) {
    auto* solid = new G4Tubs(name, 0.0, radiusCm * cm, kWall2 / 2.0 * cm, 0.0, 360.0 * deg);
    auto* lv = new G4LogicalVolume(solid, holeGas, name);
    lv->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.1)));
    new G4PVPlacement(nullptr, G4ThreeVector(0.0, -kBeamHeightCm * cm, zLocalCm * cm), lv, name,
                       cmbrLV, false, 0, true);
  };
  // Upstream (entrance, "U") sits at negative z, downstream (exit, "D")
  // at positive z -- same convention as the EN/EX chain above (kZent1
  // etc. negative, EX2G etc. positive).
  uhol("UHOL", gas1, 0.4, -(kBoxLength / 2.0 - kWall2 / 2.0));
  uhol("DHOL", gas1, 0.45, kBoxLength / 2.0 - kWall2 / 2.0);

  auto insertCollimator = [&](const char* namePrefix, G4Material* holeGas, double radiusCm,
                               double zLocalCm) {
    CollarTubes(cmbgLV, namePrefix, zLocalCm, kBpaInLenCm,
                {{aluminum, 1.905}, {holeGas, radiusCm}}, alVis);
  };
  insertCollimator("PUAI", gas1, 0.4, -(kBoxLength / 2.0 - kWall2 - kBpaInLenCm));
  insertCollimator("PDAI", gas1, 0.45, kBoxLength / 2.0 - kWall2 - kBpaInLenCm);

  // Upstream (entrance-side) lettered collar chain, PUA1..PUJ2.
  const std::vector<CollarSection> upstream = {
      {"PUA", 0.472, {{aluminum, 1.905}, {gas1, 0.4}}},
      {"PUB", 1.437, {{lead, 1.353}, {aluminum, 1.035}, {gas1, 0.4}}},
      {"PUC", 0.321, {{lead, 1.353}, {aluminum, 1.035}, {gas2, 0.45}}},
      {"PUD", 0.159, {{lead, 2.88}, {aluminum, 1.035}, {gas2, 0.45}}},
      {"PUE", 2.060, {{lead, 3.20}, {air, 2.53}, {aluminum, 2.09}, {gas2, 0.45}}},
      {"PUF",
       0.499,
       {{lead, 3.2}, {aluminum, 2.53}, {aluminum, 2.09}, {gas3, 1.25}, {aluminum, 1.04},
        {gas3, 0.5}}},
      {"PUH", 0.980, {{aluminum, 5.71}, {gas1, 1.25}, {aluminum, 1.04}, {gas3, 0.5}}},
      {"PUI", 0.585, {{aluminum, 2.53}, {gas1, 1.25}, {aluminum, 1.04}, {gas3, 0.5}}},
      {"PUJ", 0.350, {{aluminum, 1.04}, {gas3, 0.5}}},
  };
  // Chain order is A,B,C,D,E,F,G,H,I,J; PUG alone has a lead BOX (not
  // tube) outer layer, so it's built by hand between two BuildCollarChain
  // calls -- A..F first (running offset now = 2*(A+B+C+D+E+F)), then PUG
  // at that same offset plus its own half-length, then H..J resumed from
  // an offset seeded past both F *and* G (matching the real file's own
  // z formula for PUH1, which includes both 2*bpf_len and 2*bpg_len).
  constexpr double kPugHalfLenCm = 0.476;
  const std::vector<CollarSection> abcdef(upstream.begin(), upstream.begin() + 6);  // A..F
  BuildCollarChain(worldLV, kBoxLength / 2.0, -1.0, abcdef, leadVis);

  const double pugStartCm =
      kBoxLength / 2.0 + 2.0 * (0.472 + 1.437 + 0.321 + 0.159 + 2.060 + 0.499);
  const double pugZCm = -(pugStartCm + kPugHalfLenCm);  // upstream: negative z
  auto* pug1Solid = new G4Box("PUG1", 7.75 * cm, 11.75 * cm, kPugHalfLenCm * cm);
  auto* pug1LV = new G4LogicalVolume(pug1Solid, lead, "PUG1");
  pug1LV->SetVisAttributes(leadVis);
  new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, pugZCm * cm), pug1LV, "PUG1", worldLV, false, 0,
                     true);
  CollarTubes(pug1LV, "PUG2", 0.0, kPugHalfLenCm,
              {{aluminum, 2.53}, {gas1, 1.25}, {aluminum, 1.04}, {gas3, 0.5}}, alVis);

  const double startAfterPugCm = pugStartCm + 2.0 * kPugHalfLenCm;
  const std::vector<CollarSection> hij(upstream.begin() + 6, upstream.end());  // H, I, J
  BuildCollarChain(worldLV, startAfterPugCm, -1.0, hij, leadVis);

  // Downstream (exit-side) chain, PDA1..PDI4 -- asymmetric vs. upstream:
  // no separate 'C' (combined into PDD alongside D), no lead layer in
  // PDB/PDE (aluminum outermost instead), no box section (PDG doesn't
  // exist), no PDJ at all (chain stops after PDI). All exactly as the
  // real file gives it, not harmonized with the upstream side.
  const std::vector<CollarSection> downstream = {
      {"PDA", 0.472, {{aluminum, 1.905}, {gas1, 0.450}}},
      {"PDB", 1.437, {{aluminum, 1.035}, {gas1, 0.450}}},
      {"PDD", 0.321 + 0.159, {{aluminum, 1.035}, {gas2, 0.520}}},  // combines bpc_len+bpd_len
      {"PDE", 2.060, {{aluminum, 2.09}, {gas2, 0.520}}},
      {"PDF",
       0.499 + 0.476,  // combines bpf_len+bpg_len
       {{aluminum, 2.53}, {gas3, 1.25}, {aluminum, 1.04}, {gas3, 0.591}}},
      {"PDH", 0.980, {{aluminum, 5.71}, {gas3, 1.25}, {aluminum, 1.04}, {gas3, 0.591}}},
      {"PDI", 0.585, {{aluminum, 2.53}, {gas3, 1.25}, {aluminum, 1.04}, {gas3, 0.591}}},
  };
  BuildCollarChain(worldLV, kBoxLength / 2.0, 1.0, downstream, alVis);

  return cmbgLV;
}

}  // namespace

void TargetChamber::Build(G4LogicalVolume* worldLV, TargetGas gas) {
  G4LogicalVolume* cmbgLV = OuterBoxAndShielding(worldLV, gas);
  CellAndApertures(cmbgLV, gas);
  DifferentialPumpingChain(worldLV, gas);
}

double TargetChamber::BeamApertureWorldYCm() {
  // See kApertureLocalZCm's own comment: CELL's local z=kApertureLocalZCm
  // maps, through CellRotation(), to world +y (see CellAndApertures()'s
  // EAPG/XAPG placements); CELL's own placement inside CMBG/CMBR
  // (kCellYInCmbgCm + kBeamHeightCm) puts its *center* at this same world y
  // plus kApertureLocalZCm. A primary vertex meant to fire recoils down the
  // gas cell's real (aperture-to-aperture) path belongs at this world y,
  // not at CELL's own center (world y = kCellYInCmbgCm + kBeamHeightCm).
  return kCellYInCmbgCm + kBeamHeightCm + kApertureLocalZCm;
}

double TargetChamber::BeamPathHalfLengthCm() { return kApertureLocalXCm; }

const char* TargetChamber::TargetGasMaterialName() { return kTargetGasMaterialName; }
