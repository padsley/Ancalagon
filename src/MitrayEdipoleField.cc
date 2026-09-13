#include "MitrayEdipoleField.hh"

#include <cmath>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "RotateAboutY.hh"

namespace {
// mitray_edipol.f's own B->C position transform uses this constant...
constexpr double kD2rBtoC = 57.29578;
// ...while MITRAY_EDIP_ECTOEA (the C->A field transform) uses this
// different, more precise one. Both are reproduced here exactly as
// written, matching the equivalent inconsistency already noted in
// MitrayDipoleField for the magnetic dipole.
constexpr double kD2rFieldTransform = 57.29577951;
}  // namespace

MitrayEdipoleData MitrayEdipoleData::E1() {
  // Verbatim from the 'EDIP' 'E1 ' card in dat/dragon_2014_DSSSD.dat.
  MitrayEdipoleData d;
  d.D = 10.;
  d.RB = 200.;
  d.EFF = -4.721;
  d.PHI = 20.;
  d.EC2 = 0.;
  d.EC4 = 0.;
  d.WE = 0.;  // mitray_edipol.f defaults WE to 1000*RB when read as 0
  d.Z11 = 15.1658;
  d.Z12 = -14.6504;
  d.Z21 = -14.6504;
  d.Z22 = 15.1658;
  const double c[6] = {0.07901, 3.90918, -0.65329, 1.91401, 0.22838, -0.80791};
  for (int i = 0; i < 6; ++i) {
    d.entranceC[i] = c[i];
    d.exitC[i] = c[i];
  }
  return d;
}

MitrayEdipoleData MitrayEdipoleData::E2() {
  // Verbatim from the 'EDIP' 'E2 ' card in dat/dragon_2014_DSSSD.dat.
  MitrayEdipoleData d;
  d.D = 10.;
  d.RB = 250.;
  d.EFF = -3.77557049;
  d.PHI = 35.;
  d.EC2 = 0.;
  d.EC4 = 0.;
  d.WE = 50.;  // nonzero here, but moot: EC2=EC4=0 for E2 too
  d.Z11 = 15.1658;
  d.Z12 = -14.6504;
  d.Z21 = -14.6504;
  d.Z22 = 15.1658;
  const double c[6] = {0.07901, 3.90918, -0.65329, 1.91401, 0.22838, -0.80791};
  for (int i = 0; i < 6; ++i) {
    d.entranceC[i] = c[i];
    d.exitC[i] = c[i];
  }
  return d;
}

MitrayEdipoleField::MitrayEdipoleField(const MitrayEdipoleData& data,
                                        const G4ThreeVector& worldCenterCm, double thetaDeg)
    : fData(data), fCenterCm(worldCenterCm), fThetaRad(thetaDeg * CLHEP::pi / 180.0) {}

// --- MITRAY_EDPP -----------------------------------------------------------
void MitrayEdipoleField::Edpp(double d, double rb, double s, const double c[6], double& re,
                               double g[4]) {
  const double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s;
  const double cs_raw = c[0] + c[1] * s + c[2] * s2 + c[3] * s3 + c[4] * s4 + c[5] * s5;

  const double rbd = rb / d;
  const double cp1 = (c[1] + 2. * c[2] * s + 3. * c[3] * s2 + 4. * c[4] * s3 + 5. * c[5] * s4) * rbd;
  const double cp2 = (2. * c[2] + 6. * c[3] * s + 12. * c[4] * s2 + 20. * c[5] * s3) * rbd * rbd;
  const double cp3 = (6. * c[3] + 24. * c[4] * s + 60. * c[5] * s2) * rbd * rbd * rbd;
  const double cp4 = (24. * c[4] + 120. * c[5] * s) * rbd * rbd * rbd * rbd;

  double cs = cs_raw;
  if (std::abs(cs) > 70.) cs = std::copysign(70., cs);
  const double e = std::exp(cs);
  re = 1. / (1. + e);
  const double ere = e * re;
  const double ere1 = ere * re, ere2 = ere * ere1, ere3 = ere * ere2, ere4 = ere * ere3;

  const double cp12 = cp1 * cp1, cp13 = cp1 * cp12, cp14 = cp12 * cp12, cp22 = cp2 * cp2;

  g[0] = -cp1 * ere1;  // G1
  g[1] = -(cp2 + cp12) * ere1 + 2. * cp12 * ere2;  // G2
  g[2] = -(cp3 + 3. * cp1 * cp2 + cp13) * ere1 + 6. * (cp1 * cp2 + cp13) * ere2 -
         6. * cp13 * ere3;  // G3
  g[3] = -(cp4 + 4. * cp1 * cp3 + 3. * cp22 + 6. * cp12 * cp2 + cp14) * ere1 +
         (8. * cp1 * cp3 + 36. * cp12 * cp2 + 6. * cp22 + 14. * cp14) * ere2 -
         36. * (cp12 * cp2 + cp14) * ere3 + 24. * cp14 * ere4;  // G4
}

// --- MITRAY_EDIPOL (top level) + MITRAY_EDIP + MITRAY_EDIP_ECTOEA --------
void MitrayEdipoleField::FieldInLocalCm(double xa, double ya, double za, double efld[3]) const {
  const MitrayEdipoleData& p = fData;

  const double xb = -xa, yb = ya, zb = -za;  // A = 0 for E1/E2 (za already A-za with A=0)

  const double copab = std::cos(p.PHI / kD2rBtoC);
  const double sipab = std::sin(p.PHI / kD2rBtoC);
  const double cospb = std::cos((p.PHI / 2.) / kD2rBtoC);
  const double sinpb = std::sin((p.PHI / 2.) / kD2rBtoC);
  const double sip2 = std::sin((p.PHI / 2.) / kD2rBtoC);

  const double zc = -zb * copab + xb * sipab - 2. * p.RB * sip2 * cospb;
  const double xc = -zb * sipab - xb * copab - 2. * p.RB * sip2 * sinpb;
  const double yc = yb;

  const double xbmax = p.D / 2., xbmin = -xbmax;
  const double xcmax = xbmax, xcmin = xbmin;

  auto edip = [&](int in, double x, double y, double z, double xcOffset, double ef,
                  const double c[6], double out[3]) {
    const double dx = x - xcOffset;
    const double rp2 = dx * dx + z * z;

    if (in == 2) {
      out[0] = ef * p.RB * dx / rp2;
      out[1] = 0.;
      out[2] = ef * p.RB * z / rp2;
      return;
    }

    const double rp = std::sqrt(rp2);
    const double dr = rp - p.RB;
    const double sint = z / rp;
    const double cost = std::abs(xcOffset - x) / rp;
    const double theta = std::asin(sint);
    const double we = (p.WE == 0.) ? 1000. * p.RB : p.WE;
    const double s = theta * p.RB / p.D + p.EC2 * y * y / (we * we) +
                     p.EC4 * (y / we) * (y / we) * (y / we) * (y / we);

    double re, g[4];
    Edpp(p.D, p.RB, s, c, re, g);

    const double drr = dr / p.RB, drr2 = drr * drr, drr3 = drr2 * drr, drr4 = drr3 * drr;
    const double efr =
        ef * p.RB * (re - drr2 * g[1] / 2. + drr3 * g[1] / 2. + drr4 * (g[3] - 11. * g[1]) / 24.) /
        rp;
    const double eft = ef * p.RB *
                        (drr * g[0] - drr2 * g[0] / 2. + drr3 * (2. * g[0] - g[2]) / 6. -
                         drr4 * (g[0] - g[2]) / 4.) /
                        rp;

    double ex = efr * cost - eft * sint;
    double ez = efr * sint + eft * cost;
    if (in == 1) ez = -ez;

    out[0] = ex;
    out[1] = 0.;
    out[2] = ez;
  };

  auto ectoea = [&](const double c_axis[3], double a_axis[3]) {
    const double copab2 = std::cos(-p.PHI / kD2rFieldTransform);
    const double sipab2 = std::sin(-p.PHI / kD2rFieldTransform);
    const double ezb = -c_axis[2] * copab2 + c_axis[0] * sipab2;
    const double exb = -c_axis[2] * sipab2 - c_axis[0] * copab2;
    const double eyb = c_axis[1];
    a_axis[0] = exb;
    a_axis[1] = eyb;
    a_axis[2] = ezb;
  };

  if (zb <= p.Z11 && zb > p.Z12 && xb >= xbmin && xb <= xbmax) {
    // Entrance fringe field -- B-axis and A-axis coincide for the E-dipole
    // (pure translation, no rotation), so no transform is needed here.
    edip(1, xb, yb, zb, p.RB, p.EFF, p.entranceC, efld);
    return;
  }

  if (zc > p.Z21 && zc <= p.Z22 && xc >= xcmin && xc <= xcmax) {
    // Exit fringe field.
    double local[3];
    edip(3, xc, yc, zc, -p.RB, -p.EFF, p.exitC, local);
    ectoea(local, efld);
    return;
  }

  if (zb <= p.Z12 && zc <= p.Z21) {
    // Uniform (cylindrical 1/r) field region.
    double local[3];
    edip(2, xc, yc, zc, -p.RB, -p.EFF, p.exitC, local);
    ectoea(local, efld);
    return;
  }

  // Entrance far field / exit far field / unspecified region: E=0.
  efld[0] = efld[1] = efld[2] = 0.;
}

void MitrayEdipoleField::GetFieldValue(const G4double point[4], G4double* bField) const {
  const double dxCm = (point[0] - fCenterCm.x() * cm) / cm;
  const double dyCm = (point[1] - fCenterCm.y() * cm) / cm;
  const double dzCm = (point[2] - fCenterCm.z() * cm) / cm;

  double xa, ya, za;
  RotateWorldToLocal(dxCm, dyCm, dzCm, fThetaRad, xa, ya, za);

  double efld[3];
  FieldInLocalCm(xa, ya, za, efld);

  double ex, ey, ez;
  RotateLocalToWorld(efld[0], efld[1], efld[2], fThetaRad, ex, ey, ez);

  // raw efld is kV/cm (see header); G4EqMagElectricField expects
  // Bfield[0..2]=B, Bfield[3..5]=E in G4 internal units.
  bField[0] = 0.;
  bField[1] = 0.;
  bField[2] = 0.;
  bField[3] = ex * (kilovolt / cm);
  bField[4] = ey * (kilovolt / cm);
  bField[5] = ez * (kilovolt / cm);
}
