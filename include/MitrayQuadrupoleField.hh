// Port of src/mitray_poles.f + src/mitray_zone.f (MIT-RAYTRACE multipole
// field, as adapted for GEANT3 in the DRAGON simulation) to a G4MagneticField.
//
// The Fortran routines work in an element-local "A-coordinate system", in
// cm, and return B in Tesla. This class reproduces that arithmetic exactly
// (same formulas, same zone logic) in double precision, then converts to/from
// G4 internal units at the boundary.
//
// Element parameters below are named after the Fortran DATA(1..42) array
// as read by mitray_setup.f for a 'POLE' card:
//   DATA(1)=LF1  DATA(2)=LU1  DATA(3)=LF2   (unused by the field calc itself)
//   DATA(10)=A   DATA(11)=B   DATA(12)=L    DATA(13)=RAD
//   DATA(14)=BQD DATA(15)=BHX DATA(16)=BOC  DATA(17)=BDC  DATA(18)=BDD
//   DATA(19)=Z11 DATA(20)=Z12 DATA(21)=Z21  DATA(22)=Z22
//   DATA(23..28)=C0..C5 (entrance fringe)   DATA(29..34)=C0..C5 (exit fringe)
//   DATA(35..38)=FRH,FRO,FRD,FRDD           DATA(39..42)=DSH,DSO,DSD,DSDD
#pragma once

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

struct MitrayPoleData {
  double LF1 = 0, LU1 = 0, LF2 = 0;
  double A = 0, B = 0, L = 0, RAD = 0;
  double BQD = 0, BHX = 0, BOC = 0, BDC = 0, BDD = 0;
  double Z11 = 0, Z12 = 0, Z21 = 0, Z22 = 0;
  double entranceC[6] = {0, 0, 0, 0, 0, 0};  // C0..C5
  double exitC[6] = {0, 0, 0, 0, 0, 0};      // C0..C5
  double FRH = 0, FRO = 0, FRD = 0, FRDD = 0;
  double DSH = 0, DSO = 0, DSD = 0, DSDD = 0;

  // Q1 quadrupole of the DRAGON separator, dat/dragon_2014_DSSSD.dat
  static MitrayPoleData Q1();
  // Q2 quadrupole -- same file, immediately downstream of Q1. Unlike Q1,
  // Q2 has a small nonzero hexapole component (BHX), so it exercises the
  // multipole terms beyond the pure quadrupole one that Q1's grid never
  // touches.
  static MitrayPoleData Q2();
  // Q3 through Q14: every other quadrupole in the same file's beamline.
  // Every 'POLE' card in dat/dragon_2014_DSSSD.dat has LF1=LU1=LF2=3 and
  // A=B=0, and Z11==Z22 / Z12==Z21 (a "far"/"near" fringe-zone boundary
  // pair, common to entrance and exit) -- FromCard() below is that common
  // shape; Q1()/Q2() above predate it and are left as their own literal
  // transcriptions.
  static MitrayPoleData Q3();
  static MitrayPoleData Q4();
  static MitrayPoleData Q5();
  static MitrayPoleData Q6();
  static MitrayPoleData Q7();
  static MitrayPoleData Q8();
  static MitrayPoleData Q9();
  static MitrayPoleData Q10();
  static MitrayPoleData Q11();
  static MitrayPoleData Q12();
  static MitrayPoleData Q13();
  static MitrayPoleData Q14();

 private:
  static MitrayPoleData FromCard(double L, double RAD, double BQD, double BHX, double BOC,
                                  double BDC, double BDD, double zFar, double zNear,
                                  const double c[6]);
};

class MitrayQuadrupoleField : public G4MagneticField {
 public:
  // worldCenterCm: world-frame position, in cm (plain numbers, not G4
  // quantities), of the element's A-axis origin (XA=YA=ZA=0).
  // thetaDeg: rotation of the element's local A-axis frame about the
  // world Y axis, in degrees -- i.e. how far the beam direction has been
  // bent (by upstream dipoles/electrostatic deflectors) by the time it
  // reaches this element. 0 means the local z axis is aligned with world
  // z (true for an unrotated/first-in-chain element). GetFieldValue()
  // rotates the query point into the local frame and the resulting field
  // back into world coordinates; FieldInLocalCm() itself is unaffected
  // (and thus so is its existing bit-exact validation).
  MitrayQuadrupoleField(const MitrayPoleData& data, const G4ThreeVector& worldCenterCm,
                        double thetaDeg = 0.0);

  void GetFieldValue(const G4double point[4], G4double* bField) const override;

  // Evaluate the field directly in the element-local A-coordinate system,
  // cm in / Tesla out -- the literal port of MITRAY_POLES. Used both by
  // GetFieldValue() and by the standalone validation probe so both paths
  // share one implementation.
  void FieldInLocalCm(double xa, double ya, double za, double bfld[3]) const;

  // Total physical half-length along z of GEANT3's ugeo_mpole bounding
  // volume: (L + Z11 + Z22) / 2, in cm.
  double HalfLengthCm() const { return (fData.L + fData.Z11 + fData.Z22) / 2.0; }
  double ApertureRadiusCm() const { return fData.RAD; }

 private:
  static int Zone(double zb, double zc, double z11, double z12, double z21, double z22);
  static void Bpls(int igp, double d, double s, const double c[6], double& re, double g[6]);
  static void Bpoles(int in, double x, double y, double z, double d, double dh, double dOv,
                      double dd, double ddd, double dsh, double dso, double dsd, double dsdd,
                      const double grad[5], const double c[6], double& bx, double& by,
                      double& bz);

  MitrayPoleData fData;
  G4ThreeVector fCenterCm;
  double fThetaRad;
};
