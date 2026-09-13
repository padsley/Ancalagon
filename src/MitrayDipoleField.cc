#include "MitrayDipoleField.hh"

#include <cmath>
#include <stdexcept>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "RotateAboutY.hh"

namespace {
// The Fortran consistently uses this divisor (not exactly 180/pi) to
// convert degrees to radians in mitray_dipole.f's coordinate transforms.
constexpr double kD2R = 57.29577951;
// ...except the entrance XC_OFFSET/ZC_OFFSET calculation in
// MITRAY_DIPOLE uses this slightly different (lower-precision) constant --
// a pre-existing inconsistency in the original source, reproduced here
// for bit-exact fidelity rather than "fixed".
constexpr double kD2rEntranceOffset = 57.29578;
}  // namespace

MitrayDipoleData MitrayDipoleData::D1() {
  // Verbatim from the 'DIPO' 'D1  ' card in dat/dragon_2014_DSSSD.dat.
  MitrayDipoleData d;
  d.DG = 2.;
  d.MTYP = 3;
  d.IMAP = 0;
  d.A = 0.;
  d.Bp = 0.;
  d.D = 10.;
  d.RB = 100.;
  d.BF = 0.21561204;
  d.PHI = 50.;
  d.ALPHA = 5.8;
  d.BETA = 5.8;
  d.NDX = 0.;
  d.BET1 = 0.;
  d.GAMA = 0.;
  d.DELT = 0.;
  d.Z11 = 30.;
  d.Z12 = -22.;
  d.Z21 = -22.;
  d.Z22 = 30.;
  const double c[6] = {0.2877, 3.52101, -1.02159, -0.049652, 0.133009, -0.0193801};
  for (int i = 0; i < 6; ++i) {
    d.entranceC[i] = c[i];
    d.exitC[i] = c[i];
  }
  d.BR1 = 0.;
  d.BR2 = 0.;
  d.XCR1 = 0.;
  d.XCR2 = 0.;
  d.DELS1 = 0.;
  d.DELS2 = 0.;
  d.RCA1 = 0.;
  d.RCA2 = 0.;
  d.WDIP1 = 100.;
  d.WDIP2 = 100.;
  // Sraw_entrance / Sraw_exit already zero-initialized -> flat pole faces.
  return d;
}

MitrayDipoleData MitrayDipoleData::D2() {
  // Verbatim from the 'DIPO' 'D2  ' card in dat/dragon_2014_DSSSD.dat.
  MitrayDipoleData d;
  d.DG = 1.;
  d.MTYP = 3;
  d.IMAP = 0;
  d.A = 0.;
  d.Bp = 0.;
  d.D = 12.;
  d.RB = 81.3;
  d.BF = 0.26520546;
  d.PHI = 75.;
  d.ALPHA = 29.;
  d.BETA = 29.;
  d.NDX = 0.;
  d.BET1 = 0.;
  d.GAMA = 0.;
  d.DELT = 0.;
  d.Z11 = 36.;
  d.Z12 = -30.;
  d.Z21 = -30.;
  d.Z22 = 36.;
  const double c[6] = {0.3295, 3.31886, -1.2036, 0.181157, 0.1103868, -0.029513};
  for (int i = 0; i < 6; ++i) {
    d.entranceC[i] = c[i];
    d.exitC[i] = c[i];
  }
  d.BR1 = 0.;
  d.BR2 = 0.;
  d.XCR1 = 0.;
  d.XCR2 = 0.;
  d.DELS1 = 0.;
  d.DELS2 = 0.;
  d.RCA1 = 0.;
  d.RCA2 = 0.;
  d.WDIP1 = 80.;
  d.WDIP2 = 80.;
  // Sraw_entrance / Sraw_exit already zero-initialized -> flat pole faces.
  return d;
}

MitrayDipoleField::FaceShape MitrayDipoleField::MakeFaceShape(const double sraw[7], double rca,
                                                                double rb) {
  FaceShape f;
  const double rb2 = rb * rb, rb3 = rb2 * rb, rb4 = rb3 * rb;
  const double rb5 = rb4 * rb, rb6 = rb5 * rb, rb7 = rb6 * rb;
  f.s[0] = sraw[0] / rb + rca / 2.;                       // S2
  f.s[1] = sraw[1] / rb2;                                 // S3
  f.s[2] = sraw[2] / rb3 + rca * rca * rca / 8.;          // S4
  f.s[3] = sraw[3] / rb4;                                 // S5
  f.s[4] = sraw[4] / rb5 + std::pow(rca, 5) / 16.;        // S6
  f.s[5] = sraw[5] / rb6;                                 // S7
  f.s[6] = sraw[6] / rb7 + std::pow(rca, 7) / 25.6;       // S8
  f.flat = true;
  for (double v : f.s) {
    if (v != 0.) f.flat = false;
  }
  return f;
}

double MitrayDipoleField::Zefb(const FaceShape& shape, double x) {
  const double x2 = x * x, x3 = x2 * x, x4 = x3 * x;
  const double zefb = -(shape.s[0] * x2 + shape.s[1] * x3 + shape.s[2] * x4 + shape.s[3] * x4 * x +
                         shape.s[4] * x4 * x2 + shape.s[5] * x4 * x3 + shape.s[6] * x4 * x4);
  return zefb;
}

MitrayDipoleField::MitrayDipoleField(const MitrayDipoleData& data,
                                      const G4ThreeVector& worldCenterCm, double thetaDeg)
    : fData(data), fCenterCm(worldCenterCm), fThetaRad(thetaDeg * CLHEP::pi / 180.0) {
  if (fData.MTYP != 3) {
    throw std::runtime_error(
        "MitrayDipoleField: only MTYP=3 (nonuniform/gradient, standard approximation) is "
        "ported; got MTYP=" +
        std::to_string(fData.MTYP));
  }
  if (fData.IMAP != 0) {
    throw std::runtime_error("MitrayDipoleField: field-map dipoles (IMAP!=0) are not ported");
  }
  fEntranceShape = MakeFaceShape(fData.Sraw_entrance, fData.RCA1, fData.RB);
  fExitShape = MakeFaceShape(fData.Sraw_exit, fData.RCA2, fData.RB);
  if (!fEntranceShape.flat || !fExitShape.flat) {
    throw std::runtime_error(
        "MitrayDipoleField: curved pole faces (NSRF=1) are not ported -- only flat faces "
        "(as D1 has) are supported");
  }
}

// --- MITRAY_DIPO_ATOB / _BTOC / _BBTOBA / _BCTOBB ------------------------
void MitrayDipoleField::AtoB(double xa, double ya, double za, double& xb, double& yb,
                              double& zb) const {
  const double cosa = std::cos(fData.ALPHA / kD2R);
  const double sina = std::sin(fData.ALPHA / kD2R);
  xb = (fData.A - za) * sina - xa * cosa;
  yb = ya;
  zb = (fData.A - za) * cosa + xa * sina;
  xb -= fData.XCR1;
}

void MitrayDipoleField::BtoC(double xb, double yb, double zb, double& xc, double& yc,
                              double& zc) const {
  const double copab = std::cos((fData.PHI - fData.ALPHA - fData.BETA) / kD2R);
  const double sipab = std::sin((fData.PHI - fData.ALPHA - fData.BETA) / kD2R);
  const double cospb = std::cos((fData.PHI / 2. - fData.BETA) / kD2R);
  const double sinpb = std::sin((fData.PHI / 2. - fData.BETA) / kD2R);
  const double sip2 = std::sin((fData.PHI / 2.) / kD2R);

  const double xt = xb + fData.XCR1;
  const double zt = zb;

  zc = -zt * copab + xt * sipab - 2. * fData.RB * sip2 * cospb;
  xc = -zt * sipab - xt * copab - 2. * fData.RB * sip2 * sinpb;
  xc -= fData.XCR2;
  yc = yb;
}

void MitrayDipoleField::BbToBa(double bxb, double byb, double bzb, double& bxa, double& bya,
                                double& bza) const {
  const double cosa = std::cos(-fData.ALPHA / kD2R);
  const double sina = std::sin(-fData.ALPHA / kD2R);
  bxa = -bzb * sina - bxb * cosa;
  bya = byb;
  bza = -bzb * cosa + bxb * sina;
}

void MitrayDipoleField::BcToBb(double bxc, double byc, double bzc, double& bxb, double& byb,
                                double& bzb) const {
  const double copab = std::cos(-(fData.PHI - fData.ALPHA - fData.BETA) / kD2R);
  const double sipab = std::sin(-(fData.PHI - fData.ALPHA - fData.BETA) / kD2R);
  bzb = -bzc * copab + bxc * sipab;
  bxb = -bzc * sipab - bxc * copab;
  byb = byc;
}

// --- MITRAY_SDIP (flat) + MITRAY_SIJ, folded ------------------------------
double MitrayDipoleField::FringeS(const FaceShape& shape, double x, double z, double dels, int iz,
                                   int jx) const {
  const double ss = z / fData.D + dels;
  if (iz == 0 && jx == 0) return ss;
  const double zo = Zefb(shape, x);
  const double dsd = -(Zefb(shape, x + jx * fData.DG) - zo);
  return ss + (iz * fData.DG + dsd) / fData.D;
}

// --- MITRAY_NDPP (+ MITRAY_NDPPX) -----------------------------------------
double MitrayDipoleField::NdppScalar(double s, double dr, double br, const double c[6]) const {
  const double drr1 = dr / fData.RB, drr2 = drr1 * drr1, drr3 = drr2 * drr1, drr4 = drr3 * drr1;
  double cs = c[0] + s * (c[1] + s * (c[2] + s * (c[3] + s * (c[4] + s * c[5]))));
  if (std::abs(cs) > 70.) cs = std::copysign(70., cs);
  const double e = std::exp(cs);
  const double p0 = 1. + e;
  const double db = fData.BF - br;

  double bfld = 0.;
  if (fData.MTYP == 3) {
    bfld = br +
           (1. - fData.NDX * drr1 + fData.BET1 * drr2 + fData.GAMA * drr3 + fData.DELT * drr4) *
               db / p0;
  } else if (fData.MTYP == 4) {
    bfld = br + (1. / (1. + fData.NDX * drr1)) * db / p0;
  }
  return bfld;
}

// --- MITRAY_NDIP -----------------------------------------------------------
void MitrayDipoleField::Ndip(int in, double x, double y, double z, double xcOffset,
                              double zcOffset, const FaceShape& shape, const double c[6],
                              double dels, double br, double bfld[3]) const {
  const double dx = x - xcOffset, dz = z - zcOffset;
  const double rp = std::sqrt(dx * dx + dz * dz);
  const double dr = rp - fData.RB;
  const double NDX = fData.NDX, BET1 = fData.BET1, GAMA = fData.GAMA, DELT = fData.DELT;
  const double BF = fData.BF, RB = fData.RB;

  if (in == 2) {
    // Uniform (interior) field region -- 3rd-order on/off-midplane
    // expansion in the dipole's field index and shape coefficients.
    const double drr1 = dr / RB, drr2 = drr1 * drr1, drr3 = drr2 * drr1, drr4 = drr3 * drr1;
    if (y == 0.) {
      bfld[0] = 0.;
      bfld[1] = BF * (1. - NDX * drr1 + BET1 * drr2 + GAMA * drr3 + DELT * drr4);
      bfld[2] = 0.;
      return;
    }
    const double yr1 = y / RB, yr2 = yr1 * yr1, yr3 = yr2 * yr1, yr4 = yr3 * yr1;
    const double rr1 = RB / rp, rr2 = rr1 * rr1, rr3 = rr2 * rr1;
    const double brr =
        BF * ((-NDX + 2. * BET1 * drr1 + 3. * GAMA * drr2 + 4. * DELT * drr3) * yr1 -
              (NDX * rr2 + 2. * BET1 * rr1 * (1. - rr1 * drr1) +
               3. * GAMA * (2. + 2. * rr1 * drr1 - rr2 * drr2) +
               4. * DELT * (6. * drr1 + 3. * rr1 * drr2 - rr2 * drr3)) *
                  yr3 / 6.);
    const double by =
        BF * (1. - NDX * drr1 + BET1 * drr2 + GAMA * drr3 + DELT * drr4 -
              .5 * yr2 *
                  (-NDX * rr1 + 2. * BET1 * (1. + rr1 * drr1) +
                   3. * GAMA * drr1 * (2. + rr1 * drr1) + 4. * DELT * drr2 * (3. + rr1 * drr1)) +
              yr4 *
                  (-NDX * rr3 + 2. * BET1 * (rr3 * drr1 - rr2) +
                   3. * GAMA * (4. * rr1 - 2. * rr2 * drr1 + rr3 * drr2) +
                   4. * DELT * (6. + 12. * rr1 * drr1 - 3. * rr2 * drr2 + rr3 * drr3)) /
                  24.);
    bfld[0] = brr * dx / rp;
    bfld[1] = by;
    bfld[2] = brr * dz / rp;
    return;
  }

  // in == 1 (entrance) or in == 3 (exit): fringe field region.
  const double zfb = Zefb(shape, x);
  double drMid = dr;
  if (z > zfb) drMid = std::sqrt(dx * dx + (zfb - zcOffset) * (zfb - zcOffset)) - RB;

  const double sMid = FringeS(shape, x, z, dels, 0, 0);
  const double b0 = NdppScalar(sMid, drMid, br, c);

  if (y == 0.) {
    bfld[0] = 0.;
    bfld[1] = b0;
    bfld[2] = 0.;
    return;
  }

  const double DG = fData.DG;
  double dr1, dr2, dr3, dr4, dr5, dr6, dr7, dr8, dr9, dr10, dr11, dr12;

  if (z > zfb) {
    dr1 = dr2 = dr9 = dr10 = drMid;
    auto shiftedDr = [&](double xp) {
      const double zfbp = Zefb(shape, xp);
      const double dxp = xp - xcOffset;
      return std::sqrt(dxp * dxp + (zfbp - zcOffset) * (zfbp - zcOffset)) - RB;
    };
    dr3 = shiftedDr(x + DG);
    dr5 = dr3;
    dr11 = dr3;
    dr4 = shiftedDr(x - DG);
    dr7 = dr4;
    dr12 = dr4;
    dr6 = shiftedDr(x + 2. * DG);
    dr8 = shiftedDr(x - 2. * DG);
  } else {
    dr1 = std::sqrt(dx * dx + (dz + DG) * (dz + DG)) - RB;
    dr2 = std::sqrt(dx * dx + (dz + 2. * DG) * (dz + 2. * DG)) - RB;
    dr3 = std::sqrt((dx + DG) * (dx + DG) + (dz + DG) * (dz + DG)) - RB;
    dr4 = std::sqrt((dx - DG) * (dx - DG) + (dz + DG) * (dz + DG)) - RB;
    dr5 = std::sqrt((dx + DG) * (dx + DG) + dz * dz) - RB;
    dr6 = std::sqrt((dx + 2. * DG) * (dx + 2. * DG) + dz * dz) - RB;
    dr7 = std::sqrt((dx - DG) * (dx - DG) + dz * dz) - RB;
    dr8 = std::sqrt((dx - 2. * DG) * (dx - 2. * DG) + dz * dz) - RB;
    dr9 = std::sqrt(dx * dx + (dz - DG) * (dz - DG)) - RB;
    dr10 = std::sqrt(dx * dx + (dz - 2. * DG) * (dz - 2. * DG)) - RB;
    dr11 = std::sqrt((dx + DG) * (dx + DG) + (dz - DG) * (dz - DG)) - RB;
    dr12 = std::sqrt((dx - DG) * (dx - DG) + (dz - DG) * (dz - DG)) - RB;
  }

  const double b1 = NdppScalar(FringeS(shape, x, z, dels, 1, 0), dr1, br, c);
  const double b2 = NdppScalar(FringeS(shape, x, z, dels, 2, 0), dr2, br, c);
  const double b3 = NdppScalar(FringeS(shape, x, z, dels, 1, 1), dr3, br, c);
  const double b4 = NdppScalar(FringeS(shape, x, z, dels, 1, -1), dr4, br, c);
  const double b5 = NdppScalar(FringeS(shape, x, z, dels, 0, 1), dr5, br, c);
  const double b6 = NdppScalar(FringeS(shape, x, z, dels, 0, 2), dr6, br, c);
  const double b7 = NdppScalar(FringeS(shape, x, z, dels, 0, -1), dr7, br, c);
  const double b8 = NdppScalar(FringeS(shape, x, z, dels, 0, -2), dr8, br, c);
  const double b9 = NdppScalar(FringeS(shape, x, z, dels, -1, 0), dr9, br, c);
  const double b10 = NdppScalar(FringeS(shape, x, z, dels, -2, 0), dr10, br, c);
  const double b11 = NdppScalar(FringeS(shape, x, z, dels, -1, 1), dr11, br, c);
  const double b12 = NdppScalar(FringeS(shape, x, z, dels, -1, -1), dr12, br, c);

  const double yg1 = y / DG, yg2 = yg1 * yg1, yg3 = yg2 * yg1, yg4 = yg3 * yg1;
  const double bx = yg1 * ((b5 - b7) * 2. / 3. - (b6 - b8) / 12.) +
                     yg3 * ((b5 - b7) / 6. - (b6 - b8) / 12. -
                            (b3 + b11 - b4 - b12 - 2. * b5 + 2. * b7) / 12.);
  const double by = b0 -
                     yg2 * ((b1 + b9 + b5 + b7 - 4. * b0) * 2. / 3. -
                            (b2 + b10 + b6 + b8 - 4. * b0) / 24.) +
                     yg4 * (-(b1 + b9 + b5 + b7 - 4. * b0) / 6. +
                            (b2 + b10 + b6 + b8 - 4. * b0) / 24. +
                            (b3 + b11 + b4 + b12 - 2. * b1 - 2. * b9 - 2. * b5 - 2. * b7 +
                             4. * b0) /
                                12.);
  const double bz = yg1 * ((b1 - b9) * 2. / 3. - (b2 - b10) / 12.) +
                     yg3 * ((b1 - b9) / 6. - (b2 - b10) / 12. -
                            (b3 + b4 - b11 - b12 - 2. * b1 + 2. * b9) / 12.);

  bfld[0] = bx;
  bfld[1] = by;
  bfld[2] = bz;
}

// --- MITRAY_DIPOLE (top level) --------------------------------------------
void MitrayDipoleField::FieldInLocalCm(double xa, double ya, double za, double bfld[3]) const {
  double xb, yb, zb, xc, yc, zc;
  AtoB(xa, ya, za, xb, yb, zb);
  BtoC(xb, yb, zb, xc, yc, zc);

  const double xbmax = fData.WDIP1 / 2., xbmin = -xbmax;
  const double xcmax = fData.WDIP2 / 2., xcmin = -xcmax;

  double local[3];

  if (zb <= fData.Z11 && zb > fData.Z12 && xb >= xbmin && xb <= xbmax) {
    // Entrance fringe field.
    const double xcOffset = fData.RB * std::cos(fData.ALPHA / kD2rEntranceOffset);
    const double zcOffset = -fData.RB * std::sin(fData.ALPHA / kD2rEntranceOffset);
    Ndip(1, xb, yb, zb, xcOffset, zcOffset, fEntranceShape, fData.entranceC, fData.DELS1,
         fData.BR1, local);
    BbToBa(local[0], local[1], local[2], bfld[0], bfld[1], bfld[2]);
    return;
  }

  if (zc > fData.Z21 && zc <= fData.Z22 && xc >= xcmin && xc <= xcmax) {
    // Exit fringe field.
    const double xcOffset = -fData.RB * std::cos(fData.BETA / kD2R);
    const double zcOffset = -fData.RB * std::sin(fData.BETA / kD2R);
    Ndip(3, xc, yc, zc, xcOffset, zcOffset, fExitShape, fData.exitC, fData.DELS2, fData.BR2,
         local);
    double bxb, byb, bzb;
    BcToBb(local[0], local[1], local[2], bxb, byb, bzb);
    BbToBa(bxb, byb, bzb, bfld[0], bfld[1], bfld[2]);
    return;
  }

  if (zb <= fData.Z12 && zc <= fData.Z21) {
    // Uniform interior field region.
    const double xcOffset = -fData.RB * std::cos(fData.BETA / kD2R);
    const double zcOffset = -fData.RB * std::sin(fData.BETA / kD2R);
    Ndip(2, xc, yc, zc, xcOffset, zcOffset, fExitShape, fData.exitC, 0., 0., local);
    double bxb, byb, bzb;
    BcToBb(local[0], local[1], local[2], bxb, byb, bzb);
    BbToBa(bxb, byb, bzb, bfld[0], bfld[1], bfld[2]);
    return;
  }

  // Entrance far field / exit far field / unspecified region: B=0, matching
  // MITRAY_DIPOLE.
  bfld[0] = bfld[1] = bfld[2] = 0.;
}

void MitrayDipoleField::GetFieldValue(const G4double point[4], G4double* bField) const {
  const double dxCm = (point[0] - fCenterCm.x() * cm) / cm;
  const double dyCm = (point[1] - fCenterCm.y() * cm) / cm;
  const double dzCm = (point[2] - fCenterCm.z() * cm) / cm;

  double xa, ya, za;
  RotateWorldToLocal(dxCm, dyCm, dzCm, fThetaRad, xa, ya, za);

  double bfld[3];
  FieldInLocalCm(xa, ya, za, bfld);

  double bx, by, bz;
  RotateLocalToWorld(bfld[0], bfld[1], bfld[2], fThetaRad, bx, by, bz);

  bField[0] = bx * tesla;
  bField[1] = by * tesla;
  bField[2] = bz * tesla;
}
