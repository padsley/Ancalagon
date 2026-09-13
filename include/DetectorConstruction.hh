#pragma once

#include <string>

#include "G4VUserDetectorConstruction.hh"
#include "MitrayQuadrupoleField.hh"

// Minimal geometry: a vacuum world containing volumes for the
// pilot-ported MITRAY elements, selected by name:
//  - "Q1".."Q14": just that quadrupole, standalone, at theta=0 -- a tube
//    sized exactly as src/ugeom_mitray.f's ugeo_mpole sizes the real
//    GEANT3 volume (radius = RAD, half-length = (L+Z11+Z22)/2), field
//    from MitrayQuadrupoleField. (Only "Q1" and "Q2" are actually wired
//    up as standalone modes at present; Q3-Q14 are only exercised inside
//    "Chain".)
//  - "D1" / "D2": just that dipole, standalone, at theta=0 -- a
//    generously-sized bounding box (NOT a port of ugeo_dipole's exact
//    wedge/TRAP shape -- the field's own zone logic already zeroes B
//    outside the physical region, so an oversized box is physically
//    equivalent for tracking), field from MitrayDipoleField. (Only "D1"
//    is wired up as a standalone mode.)
//  - "E1" / "E2": just that electrostatic deflector, standalone, at
//    theta=0 -- a generously-sized bounding box (NOT a port of
//    ugeo_edipol's exact wedge/TRAP shape, same reasoning as D1), field
//    from MitrayEdipoleField. Uses a G4EqMagElectricField equation of
//    motion (an E field does work on the particle, unlike the
//    purely-magnetic elements) and an 8-variable stepper.
//  - "Chain": the *entire* real beamline -- Q1 -> Q2 -> D1 -> Q3 -> Q4 ->
//    Q5 -> Q6 -> Q7 -> E1 -> Q8 -> Q9 -> Q10 -> D2 -> Q11 -> Q12 -> E2 ->
//    Q13 -> Q14 -- spaced *and oriented* exactly as
//    dat/dragon_2014_DSSSD.dat's real card sequence places them (every
//    'DRIF'/'SHRT'/element card walked in DetectorConstruction.cc, with
//    line-number comments there tying each call back to the file).
//    Elements downstream of a bend are rotated (about world Y) by the
//    accumulated bend angle -- see RotateAboutY.hh and the field
//    classes' `thetaDeg` constructor parameter. GEANT3's own internal
//    algorithm for placing/orienting a dipole's physical volume
//    (src/mitray_setup.f's ugeom_setup, ~line 757-820, using the bend's
//    center-of-curvature and arc-midpoint) is not replicated -- a
//    bending element's field origin here is simply its drift-accumulated
//    entry point. This is a documented placement simplification; it does
//    not affect the already-validated field values, only where exactly
//    an element sits/points in this simplified world relative to
//    GEANT3's own bookkeeping.
class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  explicit DetectorConstruction(std::string element = "Q1") : fElement(std::move(element)) {}

  G4VPhysicalVolume* Construct() override;

  // World-frame z (cm) of the standalone element's entrance/centre --
  // exposed so the primary generator can start particles upstream of it.
  //
  // kQ1CenterZCm is also "Chain" mode's absolute anchor (see
  // DetectorConstruction.cc): 106.885 (the real STRV->Q1 drift, dat/
  // dragon_2014_DSSSD.dat lines 7-14: DF3=75 + DF4=6.5 + DST0=20 + DF5=
  // 5.385) + 12.615 (Q1's own entry-to-centre offset, (Z22+L-Z11)/2 with
  // Z11=Z22=18.89, L=25.23) = 119.5. This makes the target chamber
  // (TargetChamber.cc, built at the world origin = the real STRV/target
  // position) and Q1 sit at their real, mutually-consistent absolute
  // separation, rather than the arbitrary 50.0 placeholder used before the
  // target chamber had real absolute geometry to be consistent with.
  static constexpr double kQ1CenterZCm = 119.5;
  static constexpr double kD1CenterZCm = 50.0;
  static constexpr double kE1CenterZCm = 50.0;
  static constexpr double kE2CenterZCm = 50.0;

  // World-frame z (cm) of D1's field origin in "Chain" mode -- see
  // DetectorConstruction.cc for the derivation from the real beamline
  // file (checked there with a static_assert).
  static constexpr double kChainD1EntryZCm = 185.5;

 private:
  std::string fElement;
};
