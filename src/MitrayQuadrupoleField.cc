#include "MitrayQuadrupoleField.hh"

#include <cmath>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "RotateAboutY.hh"

MitrayPoleData MitrayPoleData::Q1() {
  // Verbatim from the 'POLE' 'Q1  ' card in dat/dragon_2014_DSSSD.dat, read
  // according to the field layout in src/mitray_setup.f.
  MitrayPoleData d;
  d.LF1 = 3.;
  d.LU1 = 3.;
  d.LF2 = 3.;
  d.A = 0.;
  d.B = 0.;
  d.L = 25.23;
  d.RAD = 5.3975;
  d.BQD = -0.09432691;
  d.BHX = 0.;
  d.BOC = 0.;
  d.BDC = 0.;
  d.BDD = 0.;
  d.Z11 = 18.89;
  d.Z12 = -13.494;
  d.Z21 = -13.494;
  d.Z22 = 18.89;
  d.entranceC[0] = 0.295;
  d.entranceC[1] = 6.30221;
  d.entranceC[2] = -3.51059;
  d.entranceC[3] = 0.29528;
  d.entranceC[4] = 1.19866;
  d.entranceC[5] = -0.423408;
  for (int i = 0; i < 6; ++i) d.exitC[i] = d.entranceC[i];
  d.FRH = 0.;
  d.FRO = 0.;
  d.FRD = 0.;
  d.FRDD = 0.;
  d.DSH = 0.;
  d.DSO = 0.;
  d.DSD = 0.;
  d.DSDD = 0.;
  return d;
}

MitrayPoleData MitrayPoleData::Q2() {
  // Verbatim from the 'POLE' 'Q2  ' card in dat/dragon_2014_DSSSD.dat.
  MitrayPoleData d;
  d.LF1 = 3.;
  d.LU1 = 3.;
  d.LF2 = 3.;
  d.A = 0.;
  d.B = 0.;
  d.L = 33.385;
  d.RAD = 7.9375;
  d.BQD = 0.08637239;
  d.BHX = 0.0045675;
  d.BOC = 0.;
  d.BDC = 0.;
  d.BDD = 0.;
  d.Z11 = 23.813;
  d.Z12 = -19.844;
  d.Z21 = -19.844;
  d.Z22 = 23.813;
  d.entranceC[0] = 0.22;
  d.entranceC[1] = 5.367112;
  d.entranceC[2] = -1.99912;
  d.entranceC[3] = 0.911917;
  d.entranceC[4] = -0.663814;
  d.entranceC[5] = 0.348883;
  for (int i = 0; i < 6; ++i) d.exitC[i] = d.entranceC[i];
  d.FRH = 0.;
  d.FRO = 0.;
  d.FRD = 0.;
  d.FRDD = 0.;
  d.DSH = 0.;
  d.DSO = 0.;
  d.DSD = 0.;
  d.DSDD = 0.;
  return d;
}

MitrayPoleData MitrayPoleData::FromCard(double L, double RAD, double BQD, double BHX, double BOC,
                                         double BDC, double BDD, double zFar, double zNear,
                                         const double c[6]) {
  MitrayPoleData d;
  d.LF1 = 3.;
  d.LU1 = 3.;
  d.LF2 = 3.;
  d.A = 0.;
  d.B = 0.;
  d.L = L;
  d.RAD = RAD;
  d.BQD = BQD;
  d.BHX = BHX;
  d.BOC = BOC;
  d.BDC = BDC;
  d.BDD = BDD;
  d.Z11 = zFar;
  d.Z12 = zNear;
  d.Z21 = zNear;
  d.Z22 = zFar;
  for (int i = 0; i < 6; ++i) {
    d.entranceC[i] = c[i];
    d.exitC[i] = c[i];
  }
  return d;
}

// 'POLE' 'Q3' -- pure hexapole (BQD=0, BHX!=0), no fringe falloff (C0..C5
// all 0 in the file).
MitrayPoleData MitrayPoleData::Q3() {
  const double c[6] = {0., 0., 0., 0., 0., 0.};
  return FromCard(18.75, 7.95, 0., 0.01831480, 0., 0., 0., 20., -20., c);
}

// 'POLE' 'Q4  '
MitrayPoleData MitrayPoleData::Q4() {
  const double c[6] = {0.225, 6.22466, -2.38148, 0.2341, -0.72032, 0.72371};
  return FromCard(33.38, 7.9375, 0.07875767, 0., 0., 0., 0., 23.813, -19.844, c);
}

// 'POLE' 'Q5  '
MitrayPoleData MitrayPoleData::Q5() {
  const double c[6] = {0.225, 6.22466, -2.38148, 0.2341, -0.72032, 0.72371};
  return FromCard(33.38, 7.9375, -0.10400018, 0., 0., 0., 0., 23.813, -19.844, c);
}

// 'POLE' 'Q6'
MitrayPoleData MitrayPoleData::Q6() {
  const double c[6] = {0.225, 6.22466, -2.38148, 0.2341, -0.72032, 0.72371};
  return FromCard(33.38, 7.9375, 0.05728922, 0., 0., 0., 0., 23.813, -19.844, c);
}

// 'POLE' 'Q7  ' -- pure hexapole, no fringe falloff (like Q3).
MitrayPoleData MitrayPoleData::Q7() {
  const double c[6] = {0., 0., 0., 0., 0., 0.};
  return FromCard(18.75, 7.95, 0., 0.00385606, 0., 0., 0., 20., -20., c);
}

// 'POLE' 'Q8' -- same L/RAD/Z/fringe shape as Q1.
MitrayPoleData MitrayPoleData::Q8() {
  const double c[6] = {0.295, 6.30221, -3.51059, 0.29528, 1.19866, -0.423408};
  return FromCard(25.23, 5.3975, -0.05091622, 0., 0., 0., 0., 18.89, -13.494, c);
}

// 'POLE' 'Q9' -- same L/RAD/Z/fringe shape as Q4-Q6.
MitrayPoleData MitrayPoleData::Q9() {
  const double c[6] = {0.225, 6.22466, -2.38148, 0.2341, -0.72032, 0.72371};
  return FromCard(33.38, 7.9375, 0.0731129, 0., 0., 0., 0., 23.813, -19.844, c);
}

// 'POLE' 'Q10' -- hexapole (BHX) + decapole (BDC), no pure-quadrupole term,
// no fringe falloff.
MitrayPoleData MitrayPoleData::Q10() {
  const double c[6] = {0., 0., 0., 0., 0., 0.};
  return FromCard(19.9, 8.0, 0., 0.0020298, 0., 0.0007228, 0., 20., -20., c);
}

// 'POLE' 'Q11' -- same L/RAD/Z/fringe shape as Q4-Q6/Q9.
MitrayPoleData MitrayPoleData::Q11() {
  const double c[6] = {0.225, 6.22466, -2.38148, 0.2341, -0.72032, 0.72371};
  return FromCard(33.38, 7.9375, 0.05420925, 0., 0., 0., 0., 23.813, -19.844, c);
}

// 'POLE' 'Q12' -- hexapole + decapole, like Q10.
MitrayPoleData MitrayPoleData::Q12() {
  const double c[6] = {0., 0., 0., 0., 0., 0.};
  return FromCard(19.9, 8.0, 0., 0.015526, 0., 0.005701, 0., 20., -20., c);
}

// 'POLE' 'Q13'
MitrayPoleData MitrayPoleData::Q13() {
  const double c[6] = {0.2535, 5.840314, -3.40247, 1.456423, 1.44575, -0.754832};
  return FromCard(46.7, 7.5, -0.04192126, 0., 0., 0., 0., 20.25, -18.75, c);
}

// 'POLE' 'Q14' -- same L/RAD/Z/fringe shape as Q13.
MitrayPoleData MitrayPoleData::Q14() {
  const double c[6] = {0.2535, 5.840314, -3.40247, 1.456423, 1.44575, -0.754832};
  return FromCard(46.7, 7.5, 0.04687356, 0., 0., 0., 0., 20.25, -18.75, c);
}

MitrayQuadrupoleField::MitrayQuadrupoleField(const MitrayPoleData& data,
                                              const G4ThreeVector& worldCenterCm, double thetaDeg)
    : fData(data), fCenterCm(worldCenterCm), fThetaRad(thetaDeg * CLHEP::pi / 180.0) {}

// --- MITRAY_ZONE ------------------------------------------------------
int MitrayQuadrupoleField::Zone(double zb, double zc, double z11, double z12, double z21,
                                 double z22) {
  if (zb > z11 || zc > z22) return 0;
  if (zb <= z12 && zc <= z21) return 2;

  bool lentr = false, lexit = false;
  if (zb <= z11 && zb > z12) lentr = true;
  if (zc <= z22 && zc > z21) lexit = true;

  if (lentr && lexit) return 4;
  if (lentr) return 1;
  if (lexit) return 3;
  return -1;
}

// --- MITRAY_BPLS --------------------------------------------------------
void MitrayQuadrupoleField::Bpls(int igp, double d, double s, const double c[6], double& re,
                                  double g[6]) {
  const double s2 = s * s, s3 = s2 * s, s4 = s2 * s2, s5 = s4 * s;

  const double cs_raw = c[0] + c[1] * s + c[2] * s2 + c[3] * s3 + c[4] * s4 + c[5] * s5;
  const double cp1 = (c[1] + 2. * c[2] * s + 3. * c[3] * s2 + 4. * c[4] * s3 + 5. * c[5] * s4) / d;
  const double cp2 = (2. * c[2] + 6. * c[3] * s + 12. * c[4] * s2 + 20. * c[5] * s3) / (d * d);
  const double cp3 = (6. * c[3] + 24. * c[4] * s + 60. * c[5] * s2) / (d * d * d);
  const double cp4 = (24. * c[4] + 120. * c[5] * s) / (d * d * d * d);
  const double cp5 = 120. * c[5] / (d * d * d * d * d);

  double cs = cs_raw;
  if (std::abs(cs) > 70.) cs = std::copysign(70., cs);
  const double e = std::exp(cs);
  re = 1. / (1. + e);
  const double ere = e * re;
  const double ere1 = ere * re, ere2 = ere * ere1, ere3 = ere * ere2, ere4 = ere * ere3,
               ere5 = ere * ere4, ere6 = ere * ere5;

  const double cp12 = cp1 * cp1, cp13 = cp1 * cp12, cp14 = cp12 * cp12, cp22 = cp2 * cp2;
  const double cp15 = cp12 * cp13, cp16 = cp13 * cp13, cp23 = cp2 * cp22, cp32 = cp3 * cp3;

  for (int i = 0; i < 6; ++i) g[i] = 0.;
  if (igp == 6) return;

  g[0] = -cp1 * ere1;  // G1
  if (igp == 5) return;

  g[1] = -(cp2 + cp12) * ere1 + 2. * cp12 * ere2;  // G2
  if (igp == 4) return;

  g[2] = -(cp3 + 3. * cp1 * cp2 + cp13) * ere1 + 6. * (cp1 * cp2 + cp13) * ere2 -
         6. * cp13 * ere3;  // G3
  if (igp == 3) return;

  g[3] = -(cp4 + 4. * cp1 * cp3 + 3. * cp22 + 6. * cp12 * cp2 + cp14) * ere1 +
         (8. * cp1 * cp3 + 36. * cp12 * cp2 + 6. * cp22 + 14. * cp14) * ere2 -
         36. * (cp12 * cp2 + cp14) * ere3 + 24. * cp14 * ere4;  // G4
  if (igp != 2) return;

  g[4] = (-cp5 - 5. * cp1 * cp4 - 10. * cp2 * cp3 - 10. * cp12 * cp3 - 15. * cp1 * cp22 -
          10. * cp13 * cp2 - cp15) *
             ere1 +
         (10. * cp1 * cp4 + 20. * cp2 * cp3 + 60. * cp12 * cp3 + 90. * cp1 * cp22 +
          140. * cp13 * cp2 + 30. * cp15) *
             ere2 +
         (-60. * cp12 * cp3 - 90. * cp1 * cp22 - 360. * cp13 * cp2 - 150. * cp15) * ere3 +
         (240. * cp13 * cp2 + 240. * cp15) * ere4 + (-120. * cp15) * ere5;  // G5

  g[5] = (-6. * cp1 * cp5 - 15. * cp2 * cp4 - 15. * cp12 * cp4 - 10. * cp32 -
          60. * cp1 * cp2 * cp3 - 20. * cp13 * cp3 - 15. * cp23 - 45. * cp12 * cp22 -
          15. * cp14 * cp2 - cp16) *
             ere1 +
         (12. * cp1 * cp5 + 30. * cp2 * cp4 + 90. * cp12 * cp4 + 20. * cp32 +
          360. * cp1 * cp2 * cp3 + 280. * cp13 * cp3 + 90. * cp23 + 630. * cp12 * cp22 +
          450. * cp14 * cp2 + 62. * cp16) *
             ere2 +
         (-90. * cp12 * cp4 - 360. * cp1 * cp2 * cp3 - 720. * cp13 * cp3 - 90. * cp23 -
          1620. * cp12 * cp22 - 2250. * cp14 * cp2 - 540. * cp16) *
             ere3 +
         (480. * cp13 * cp3 + 1080. * cp12 * cp22 + 3600. * cp14 * cp2 + 1560. * cp16) * ere4 +
         (-1800. * cp14 * cp2 - 1800. * cp16) * ere5 + 720. * cp16 * ere6;  // G6
}

// --- MITRAY_BPOLES --------------------------------------------------------
void MitrayQuadrupoleField::Bpoles(int in, double x, double y, double z, double d, double dh,
                                    double dOv, double dd, double ddd, double dsh, double dso,
                                    double dsd, double dsdd, const double grad[5],
                                    const double c[6], double& bx, double& by, double& bz) {
  const double x2 = x * x, x3 = x2 * x, x4 = x3 * x, x5 = x4 * x, x6 = x5 * x, x7 = x6 * x;
  const double y2 = y * y, y3 = y2 * y, y4 = y3 * y, y5 = y4 * y, y6 = y5 * y, y7 = y6 * y;
  const double grad1 = grad[0], grad2 = grad[1], grad3 = grad[2], grad4 = grad[3],
               grad5 = grad[4];

  if (in == 2) {
    // Uniform ("hard edge") field region -- no z dependence, no Bz.
    const double b2x = grad1 * y, b2y = grad1 * x;
    const double b3x = grad2 * 2. * x * y, b3y = grad2 * (x2 - y2);
    const double b4x = grad3 * (3. * x2 * y - y3), b4y = grad3 * (x3 - 3. * x * y2);
    const double b5x = grad4 * 4. * (x3 * y - x * y3), b5y = grad4 * (x4 - 6. * x2 * y2 + y4);
    const double b6x = grad5 * (5. * x4 * y - 10. * x2 * y3 + y5),
                 b6y = grad5 * (x5 - 10. * x3 * y2 + 5. * x * y4);
    bx = b2x + b3x + b4x + b5x + b6x;
    by = b2y + b3y + b4y + b5y + b6y;
    bz = 0.;
    return;
  }

  // Entrance (in==1) or exit (in==3) fringe field: Enge-function falloff.
  double re, g[6];

  double s = z / d;
  Bpls(2, d, s, c, re, g);
  const double b2x = grad1 * (re * y - (g[1] / 12.) * (3. * x2 * y + y3) +
                               (g[3] / 384.) * (5. * x4 * y + 6. * x2 * y3 + y5) -
                               (g[5] / 23040.) * (7. * x6 * y + 15. * x4 * y3 + 9. * x2 * y5 + y7));
  const double b2y = grad1 * (re * x - (g[1] / 12.) * (x3 + 3. * x * y2) +
                               (g[3] / 384.) * (x5 + 6. * x3 * y2 + 5. * x * y4) -
                               (g[5] / 23040.) * (x7 + 9. * x5 * y2 + 15. * x3 * y4 + 7. * x * y6));
  const double b2z = grad1 * (g[0] * x * y - (g[2] / 12.) * (x3 * y + x * y3) +
                               (g[4] / 384.) * (x5 * y + 2. * x3 * y3 + x * y5));

  double ss = z / dh + dsh;
  Bpls(3, dh, ss, c, re, g);
  const double b3x = grad2 * (re * 2. * x * y - (g[1] / 48.) * (12. * x3 * y + 4. * x * y3));
  const double b3y =
      grad2 * (re * (x2 - y2) - (g[1] / 48.) * (3. * x4 + 6. * x2 * y2 - 5. * y4));
  const double b3z =
      grad2 * (g[0] * (x2 * y - y3 / 3.) - (g[2] / 48.) * (3. * x4 * y + 2. * x2 * y3 - y5));

  ss = z / dOv + dso;
  Bpls(4, dOv, ss, c, re, g);
  const double b4x = grad3 * (re * (3. * x2 * y - y3) - (g[1] / 80.) * (20. * x4 * y - 4. * y5));
  const double b4y = grad3 * (re * (x3 - 3. * x * y2) - (g[1] / 80.) * (4. * x5 - 20. * x * y4));
  const double b4z = grad3 * g[0] * (x3 * y - x * y3);

  ss = z / dd + dsd;
  Bpls(5, dd, ss, c, re, g);
  const double b5x = grad4 * re * (4. * x3 * y - 4. * x * y3);
  const double b5y = grad4 * re * (x4 - 6. * x2 * y2 + y4);
  const double b5z = grad4 * g[0] * (x4 * y - 2. * x2 * y3 + y5 / 5.);

  ss = z / ddd + dsdd;
  Bpls(6, ddd, ss, c, re, g);
  const double b6x = grad5 * re * (5. * x4 * y - 10. * x2 * y3 + y5);
  const double b6y = grad5 * re * (x5 - 10. * x3 * y2 + 5. * x * y4);
  const double b6z = 0.;

  bx = b2x + b3x + b4x + b5x + b6x;
  by = b2y + b3y + b4y + b5y + b6y;
  bz = b2z + b3z + b4z + b5z + b6z;
}

// --- MITRAY_POLES (top level) --------------------------------------------
void MitrayQuadrupoleField::FieldInLocalCm(double xa, double ya, double za, double bfld[3]) const {
  const MitrayPoleData& p = fData;

  bfld[0] = bfld[1] = bfld[2] = 0.;

  double frh = p.FRH, fro = p.FRO, frd = p.FRD, frdd = p.FRDD;
  if (frh == 0.) frh = 1.;
  if (fro == 0.) fro = 1.;
  if (frd == 0.) frd = 1.;
  if (frdd == 0.) frdd = 1.;

  const double d = 2. * p.RAD;
  const double dh = frh * d, dOv = fro * d, dd = frd * d, ddd = frdd * d;

  const double xb = -xa, yb = ya, zb = p.A - za;
  const double xc = -xb, yc = yb, zc = -zb - p.L;

  const int izone = Zone(zb, zc, p.Z11, p.Z12, p.Z21, p.Z22);
  if (izone == 0 || izone == -1) return;

  const double grad1c = p.BQD / p.RAD;
  const double grad2c = p.BHX / (p.RAD * p.RAD);
  const double grad3c = p.BOC / (p.RAD * p.RAD * p.RAD);
  const double grad4c = p.BDC / (p.RAD * p.RAD * p.RAD * p.RAD);
  const double grad5c = p.BDD / (p.RAD * p.RAD * p.RAD * p.RAD * p.RAD);

  double bx, by, bz;

  if (izone == 2 || izone == 4) {
    const double grad[5] = {grad1c, grad2c, grad3c, grad4c, grad5c};
    Bpoles(2, xc, yc, zc, d, dh, dOv, dd, ddd, p.DSH, p.DSO, p.DSD, p.DSDD, grad, nullptr, bx, by,
           bz);
    if (izone == 2) {
      bfld[0] = bx;
      bfld[1] = by;
      bfld[2] = bz;
      return;
    }
    bfld[0] -= bx;
    bfld[1] -= by;
    bfld[2] -= bz;
  }

  if (izone == 1 || izone == 4) {
    const double grad[5] = {-grad1c, grad2c, -grad3c, grad4c, -grad5c};
    Bpoles(1, xb, yb, zb, d, dh, dOv, dd, ddd, p.DSH, p.DSO, p.DSD, p.DSDD, grad, p.entranceC, bx,
           by, bz);
    if (izone == 1) {
      bfld[0] = -bx;
      bfld[1] = by;
      bfld[2] = -bz;
      return;
    }
    bfld[0] -= bx;
    bfld[1] += by;
    bfld[2] -= bz;
  }

  if (izone == 3 || izone == 4) {
    const double grad[5] = {grad1c, grad2c, grad3c, grad4c, grad5c};
    Bpoles(3, xc, yc, zc, d, dh, dOv, dd, ddd, p.DSH, p.DSO, p.DSD, p.DSDD, grad, p.exitC, bx, by,
           bz);
    if (izone == 3) {
      bfld[0] = bx;
      bfld[1] = by;
      bfld[2] = bz;
      return;
    }
    bfld[0] += bx;
    bfld[1] += by;
    bfld[2] += bz;
  }
}

void MitrayQuadrupoleField::GetFieldValue(const G4double point[4], G4double* bField) const {
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
