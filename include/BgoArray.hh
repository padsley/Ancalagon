// Port of the DRAGON BGO gamma-ray array's crystal geometry
// (src/ugeom_gbox.f's ugeo_finger, the SCNT/MGOR/FNGR/HSNG stack, called
// from ugeom.f's ugeo_detector -> ugeo_finger -> ugeo_pmt chain), for the
// real "small pumping tube" target configuration (tubetype=0, what
// dat/dragon_2014_DSSSD.dat's TUBE FFREAD card actually selects).
//
// Each of the 30 detector positions is the real nested stack: HSNG (outer
// envelope, filled with AIR -- there is no separate solid aluminum "can"
// in the real file; FNGR itself is the aluminum housing, despite its
// name) > FNGR (aluminum) > MGOR (MgO powder reflector) > SCNT (BGO
// crystal), plus PMT (glass) as HSNG's other direct daughter. Every
// layer's own apothem/half-length is derived from the one before it by
// the real gap constants (air_gap, d_mtl, d_air), exactly as the real
// GSVOLU/GSPOS calls derive them -- see BgoArray.cc for the exact
// formulas and the one deliberate non-literal choice (real NIST glass
// for the PMT window, in place of the original file's own stale
// copy-paste 'GLASS' material, which reuses the scintillator's plastic
// CH composition -- see BgoArray.cc's own comment). There's no field or
// validated formula here the way the optics elements have, so the
// "correctness" bar for this geometry is faithful transcription of the
// real GSVOLU dimensions and GSPOS positions, which IS done exactly (see
// BgoArray.cc for the literal card-by-card derivation of each of the 30
// positions).
//
// Materials/dimensions come from dragon_2003.ffcards's FSID/WALL/BGAP/
// HOLE/PMTR cards (read by src/ugffgo.f into the COMMON block declared in
// inc/geometry.inc) -- see BgoArray.cc for the exact values and their
// card names.
#pragma once

class G4LogicalVolume;

namespace BgoArray {
// Places all 30 BGO crystals directly in `worldLV` (the real GEANT3 code
// places them directly in 'WRLD' too, via an intermediate 'DETE' mother
// box that itself sits at the world origin with no offset).
void Build(G4LogicalVolume* worldLV);
}  // namespace BgoArray
