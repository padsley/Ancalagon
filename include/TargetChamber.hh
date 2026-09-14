// Port of the DRAGON windowless recirculating gas target, transcribed
// from two GEANT3 source files in full (nothing in either is left out):
//
//   - src/ugeom_gbox.f's ugeo_detector, targtype==0 branch (~lines 1-732):
//     'CMBR' (stainless-steel outer support box) > 'CMBG' (its gas-filled
//     interior) > 'CELL' (aluminum TRD1 gas-cell outer shell) > 'CELG'
//     (mtarg-material TRD1 inner gas volume, the real reaction volume),
//     'EAPG'/'XAPG' (mtarg-material TUBE entrance/exit beam apertures
//     bored through CELL's slanted side wall) -- built by
//     OuterBoxAndShielding()/CellAndApertures() -- plus 'UHOL'/'DHOL' (beam
//     holes through CMBR's own wall), 'PUAI'/'PDAI'+'PUBI'/'PDBI' (aluminum
//     collimator inserts just inside them, in CMBG), and the lettered
//     lead/aluminum/gas collar chains bolted onto each face of the box
//     (PUA1..PUJ2 upstream, PDA1..PDI4 downstream, asymmetric -- see
//     OuterBoxAndShielding() for the exact, non-harmonized layer lists) --
//     all built by OuterBoxAndShielding().
//
//   - src/ugeom_trgt_small.f (the tubetype==0 branch dat/dragon_2014_
//     DSSSD.dat's beamline actually selects): the differential-pumping
//     chain around the cell -- 7 aluminum CONE collimator inserts
//     (EN2C/EN3C/EX2C..EX6C) nested inside 11 stepped-radius mcent-
//     material (low-pressure) gas TUBE volumes (EN1G..EN6G, EX2G..EX10),
//     at the exact Z positions the file derives from inc/uggeom.inc's
//     zent()/len1-3 constants. Built by DifferentialPumpingChain().
//
// The tubetype==1 branch inside ugeo_detector (~line 734 on, its own
// separate hardware for a different target configuration) is not ported
// -- this file's own 'TUBE' ffcard selects tubetype==0.
//
// CMBR's own medium (20) is genuinely 'STAINLESS STEEL', a standard
// GEANT3 material (material index 26, outside the custom material
// table) -- a previous version of this file's own comment mis-read that
// index against the wrong table and called it "NAI:TL"; ugstmed.f's
// medium 20 says 'STAINLESS STEEL' in its own name table, and that's
// what's built here.
//
// SCOPE/CAVEAT: this file and src/ugeom_trgt_small.f (both ported here)
// were evidently authored somewhat independently and can occupy
// overlapping z-ranges around the same physical region -- where G4's own
// overlap checker (every placement here runs with it on) flags the EN/EX
// chain crossing the PU/PD collar chain, that reflects the real files'
// own, apparently uncoordinated overlap (GEANT3's tracker tolerates
// arbitrary volume overlaps in a way GEANT4's strict mother/daughter
// model does not), not a placement error introduced by this port.
//
// COORDINATES: everything here is placed directly in the world frame,
// exactly as src/ugeom_gbox.f's 'DETE' (mother of the whole target region)
// sits in 'WRLD' with zero offset and zero rotation -- so this file's
// local (x,y,z)cm literals equal absolute world (x,y,z)cm.
#pragma once

class G4LogicalVolume;

namespace TargetChamber {

// src/ugmate_trgt.f picks the target gas by projectile mass number
// (`atarg.lt.1.2` -> hydrogen, else helium) -- both real DRAGON target
// gases, selected by whichever reaction is being run (see
// ReactionKinematics.hh/README's "Reaction specification" for the
// specific reaction this pilot's own generator implements, and which gas
// it needs).
enum class TargetGas { kHydrogen, kHelium };

void Build(G4LogicalVolume* worldLV, TargetGas gas);

// World y of the line connecting the EAPG/XAPG beam apertures (the gas
// cell's real entrance/exit) -- i.e. where a primary vertex must sit for a
// recoil fired along +z to have any chance of exiting through the cell's
// actual apertures instead of its solid wall. See TargetChamber.cc's own
// kApertureLocalZCm/BeamApertureWorldYCm() comments.
double BeamApertureWorldYCm();

// Half the beam's total path length through the target gas (EAPG's own
// center to XAPG's, both on the BeamApertureWorldYCm() line) -- i.e. the
// beam enters the gas at world z=-BeamPathHalfLengthCm() and exits at
// world z=+BeamPathHalfLengthCm(), both at world x=0,
// y=BeamApertureWorldYCm(). Used by PrimaryGeneratorAction's beam-energy-
// loss depth sampling (see GasStoppingPower) to know how much gas the
// beam actually crosses.
double BeamPathHalfLengthCm();

// Name of the G4Material every gas-filled volume in this file (CELG, EAPG,
// XAPG) is built from -- look it up with G4Material::GetMaterial(name)
// after Build() has run (which any code calling this needs to have
// happened already anyway, per the same ordering G4RunManager itself
// guarantees: DetectorConstruction::Construct() before any event
// generation), rather than building a second, merely similar material.
const char* TargetGasMaterialName();

}  // namespace TargetChamber
