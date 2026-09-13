// Port of src/mitray_edipol.f (MIT-RAYTRACE cylindrical electrostatic
// deflector field, as adapted for GEANT3 in the DRAGON simulation) to a
// G4ElectricField.
//
// This is a purely 2D (bend-plane) field model: a cylindrical capacitor of
// bend radius RB, field strength EFF, with an Enge-function fringe falloff
// at each end (same style as the magnetic elements, but here Y only enters
// through a (generally unused -- both E1 and E2 have EC2=EC4=0) shift of
// the fringe longitudinal coordinate S; there is no off-midplane
// finite-difference expansion the way the magnetic dipole/quadrupole have,
// and EY is always 0). Unlike the magnetic dipole, there is no pole-face
// rotation (ALPHA=BETA=0 always) and no curved-face (NSRF) complication --
// mitray_edipol.f has no equivalent of MITRAY_SDIP at all, so this port has
// no scope restriction beyond what src/mitray_edipol.f itself computes.
//
// Field units: src/guefld.f scales the raw EFLD magnitude from
// mitray_edipol.f by 1e-6 to get GEANT3's internal field unit of GeV/cm
// (see src/grkuta.f's docstring) -- i.e. raw EFLD is in kV/cm. XPOS is in
// cm, as for the magnetic elements.
//
// SCOPE: DATA(11)=A (an entrance-position offset, analogous to the
// magnetic dipole's A) is 0 for both E1 and E2 in
// dat/dragon_2014_DSSSD.dat, so this port hardcodes A=0 in the A->B
// position transform rather than threading it through as a parameter;
// porting a future EDIP card with A!=0 would need that restored.
//
// Element parameters below are named after the Fortran DATA(1..40) array
// as read by mitray_setup.f for an 'EDIP' card:
//   DATA(1)=LF1  DATA(2)=LU1  DATA(3)=LF2  DATA(4)=DG   (all unused by the
//                                                         field calc itself)
//   DATA(10)=A(unused) DATA(11)=B(unused) DATA(13)=D DATA(14)=RB
//   DATA(15)=EFF DATA(16)=PHI DATA(17)=EC2 DATA(18)=EC4 DATA(19)=WE
//   DATA(20)=WC (unused)
//   DATA(25..28)=Z11,Z12,Z21,Z22
//   DATA(29..34)=C0..C5 (entrance fringe)  DATA(35..40)=C0..C5 (exit fringe)
//
// NOTE the record layout mitray_setup.f actually reads for 'EDIP' differs
// slightly from the above indices in a few slots (A is DATA(11), B is
// DATA(12), D is DATA(13) per mitray_setup.f's own record boundaries) --
// see MitrayEdipoleData::E1()/E2() for the exact, verified mapping used
// here.
#pragma once

#include "G4ElectricField.hh"
#include "G4ThreeVector.hh"

struct MitrayEdipoleData {
  double D = 0, RB = 0, EFF = 0, PHI = 0;
  double EC2 = 0, EC4 = 0, WE = 0;
  double Z11 = 0, Z12 = 0, Z21 = 0, Z22 = 0;
  double entranceC[6] = {0, 0, 0, 0, 0, 0};  // C0..C5
  double exitC[6] = {0, 0, 0, 0, 0, 0};      // C0..C5

  // E1 electrostatic deflector of the DRAGON separator,
  // dat/dragon_2014_DSSSD.dat
  static MitrayEdipoleData E1();
  // E2 electrostatic deflector -- same file, further downstream.
  static MitrayEdipoleData E2();
};

class MitrayEdipoleField : public G4ElectricField {
 public:
  // worldCenterCm: world-frame position, in cm (plain numbers), of the
  // element's A-axis origin (XA=YA=ZA=0). thetaDeg: rotation of the
  // element's local A-axis frame about the world Y axis, in degrees --
  // see MitrayQuadrupoleField's constructor comment; 0 recovers the
  // original (isolated-element) behavior.
  MitrayEdipoleField(const MitrayEdipoleData& data, const G4ThreeVector& worldCenterCm,
                     double thetaDeg = 0.0);

  // Bfield[0..2] = B (always 0), Bfield[3..5] = E, matching the 6-element
  // contract G4EqMagElectricField expects.
  void GetFieldValue(const G4double point[4], G4double* bField) const override;

  // Literal port of MITRAY_EDIPOL: cm in (A-axis system) / kV/cm out.
  void FieldInLocalCm(double xa, double ya, double za, double efld[3]) const;

  double ApertureHalfGapCm() const { return fData.D / 2.0; }

 private:
  static void Edpp(double d, double rb, double s, const double c[6], double& re, double g[4]);

  MitrayEdipoleData fData;
  G4ThreeVector fCenterCm;
  double fThetaRad;
};
