// Port of src/mitray_dipole.f (MIT-RAYTRACE dipole field, GEANT3 dipole
// coordinate/field transforms) + src/mitray_zone-equivalent dispatch that
// lives inline in MITRAY_DIPOLE itself, to a G4MagneticField.
//
// SCOPE: this ports exactly the code path exercised by the real D1 dipole
// of the DRAGON separator (dat/dragon_2014_DSSSD.dat, 'DIPO' 'D1  '):
//   - MTYP=3 (nonuniform/gradient field, "standard approximation")
//   - IMAP=0 (no external field map -- MITRAY_BDMP/FMAP/DMAP not ported)
//   - flat pole faces, i.e. NSRF=0 for both entrance and exit (all of
//     D1's S2..S8 face-curvature coefficients are zero). The curved-face
//     case (NSRF=1) requires MITRAY_SDIP's iterative closest-point search
//     on the curved effective-field-boundary, which is NOT ported here;
//     the constructor throws if a curved face is supplied.
// MTYP=1,2,5 (uniform-field magnets via MITRAY_BDIP's own field-map/
// standard branches) and MTYP=6 (pretzel magnet) are out of scope.
//
// Coordinate systems (as in the Fortran): A-axis = element-local input
// axis (what GetFieldValue works in), B-axis = entrance pole-face axis,
// C-axis = exit pole-face axis. MITRAY_DIPO_ATOB/BTOC transform position
// A->B->C; MITRAY_DIPO_BBTOBA/BCTOBB transform the field back B->A, C->B.
#pragma once

#include <algorithm>

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

struct MitrayDipoleData {
  double DG = 0;
  int MTYP = 0;
  int IMAP = 0;
  double A = 0, Bp = 0, D = 0, RB = 0, BF = 0;
  double PHI = 0, ALPHA = 0, BETA = 0;
  double NDX = 0, BET1 = 0, GAMA = 0, DELT = 0;
  double Z11 = 0, Z12 = 0, Z21 = 0, Z22 = 0;
  double entranceC[6] = {0, 0, 0, 0, 0, 0};  // C0..C5
  double exitC[6] = {0, 0, 0, 0, 0, 0};      // C0..C5
  double BR1 = 0, BR2 = 0, XCR1 = 0, XCR2 = 0, DELS1 = 0, DELS2 = 0;
  double RCA1 = 0, RCA2 = 0, WDIP1 = 0, WDIP2 = 0;
  double Sraw_entrance[7] = {0, 0, 0, 0, 0, 0, 0};  // DATA(51..57)
  double Sraw_exit[7] = {0, 0, 0, 0, 0, 0, 0};      // DATA(58..64)

  // D1 dipole of the DRAGON separator, dat/dragon_2014_DSSSD.dat
  static MitrayDipoleData D1();
  // D2 dipole -- same file, further downstream (past Q3-Q10 and E1). Same
  // MTYP=3/IMAP=0/flat-face scope as D1.
  static MitrayDipoleData D2();
};

class MitrayDipoleField : public G4MagneticField {
 public:
  // worldCenterCm: world-frame position, in cm (plain numbers), of the
  // element's A-axis origin (XA=YA=ZA=0). thetaDeg: rotation of the
  // element's local A-axis frame about the world Y axis, in degrees --
  // see MitrayQuadrupoleField's constructor comment for the full
  // explanation; 0 recovers the original (isolated-element) behavior.
  MitrayDipoleField(const MitrayDipoleData& data, const G4ThreeVector& worldCenterCm,
                     double thetaDeg = 0.0);

  void GetFieldValue(const G4double point[4], G4double* bField) const override;

  // Literal port of MITRAY_DIPOLE: cm in (A-axis system) / Tesla out.
  void FieldInLocalCm(double xa, double ya, double za, double bfld[3]) const;

  double ApertureHalfGapCm() const { return fData.D / 2.0; }
  double WdipHalfWidthCm() const { return std::max(fData.WDIP1, fData.WDIP2) / 2.0; }

 private:
  // S2..S8 (7 coefficients) computed from the raw DATA(51..57 or 58..64)
  // values, RB and RCA, per lines 326-332 / 411-417 of mitray_dipole.f.
  struct FaceShape {
    double s[7] = {0, 0, 0, 0, 0, 0, 0};  // S2..S8
    bool flat = true;                     // NSRF==0
  };
  static FaceShape MakeFaceShape(const double sraw[7], double rca, double rb);
  static double Zefb(const FaceShape& shape, double x);

  // MITRAY_DIPO_ATOB / _BTOC / _BBTOBA / _BCTOBB
  void AtoB(double xa, double ya, double za, double& xb, double& yb, double& zb) const;
  void BtoC(double xb, double yb, double zb, double& xc, double& yc, double& zc) const;
  void BbToBa(double bxb, double byb, double bzb, double& bxa, double& bya, double& bza) const;
  void BcToBb(double bxc, double byc, double bzc, double& bxb, double& byb, double& bzb) const;

  // MITRAY_SDIP (flat-boundary branch only) + MITRAY_SIJ, folded into one
  // call: returns the fringe-field longitudinal coordinate S for the
  // stencil point offset (iz, jx) steps of DG from the evaluation point
  // (x, z) used for the central (iz=jx=0) call. This is the exact
  // NSRF=0 (flat pole face) reduction of SDIP+SIJ; it is not the general
  // curved-face algorithm.
  double FringeS(const FaceShape& shape, double x, double z, double dels, int iz, int jx) const;

  // MITRAY_NDPP main entry (iz=jx=0) and its MITRAY_NDPPX entry (offset
  // stencil points), merged: given S (already positioned), DR, and the
  // face's C0..C5, returns the scalar field magnitude used to build one
  // of the B1..B12 samples.
  double NdppScalar(double s, double dr, double br, const double c[6]) const;

  // MITRAY_NDIP: field in the local (B- or C-axis) frame, TC=(x,y,z).
  void Ndip(int in, double x, double y, double z, double xcOffset, double zcOffset,
            const FaceShape& shape, const double c[6], double dels, double br, double bfld[3]) const;

  MitrayDipoleData fData;
  G4ThreeVector fCenterCm;
  double fThetaRad;
  FaceShape fEntranceShape;
  FaceShape fExitShape;
};
