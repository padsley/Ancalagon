// Shared helper for chaining MITRAY elements past a bend: every bending
// element (dipole or electrostatic deflector) in dat/dragon_2014_DSSSD.dat
// bends the reference trajectory purely in the horizontal (x-z) plane --
// no element in that file ever applies a 'SHRT' rotation, and no bend has
// a vertical component -- so the whole beamline's accumulated orientation
// reduces to a single scalar angle (rotation about the world Y axis).
//
// World and local frames are related by world = R(theta) * local, with
//   world_x = local_x*cos(theta) + local_z*sin(theta)
//   world_z = -local_x*sin(theta) + local_z*cos(theta)
//   world_y = local_y
// (theta=0 recovers the identity used by every field class before chaining
// was introduced, so this is purely additive -- it doesn't change any
// already-validated theta=0 behavior.)
#pragma once

#include <cmath>

// World-frame displacement (dx, dy, dz) -> local-frame (xa, ya, za).
inline void RotateWorldToLocal(double dx, double dy, double dz, double thetaRad, double& xa,
                                double& ya, double& za) {
  const double c = std::cos(thetaRad), s = std::sin(thetaRad);
  xa = c * dx - s * dz;
  ya = dy;
  za = s * dx + c * dz;
}

// Local-frame field vector (fx, fy, fz) -> world-frame.
inline void RotateLocalToWorld(double fx, double fy, double fz, double thetaRad, double& wx,
                                double& wy, double& wz) {
  const double c = std::cos(thetaRad), s = std::sin(thetaRad);
  wx = c * fx + s * fz;
  wy = fy;
  wz = -s * fx + c * fz;
}
