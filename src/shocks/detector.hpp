#ifndef SHOCKS_DETECTOR_HPP_
#define SHOCKS_DETECTOR_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include <limits>
#include "athena.hpp"

// Restricted nearest-axis temperature-normal port of the ShockIt snapshot.
// local_pass SHA256: 6f5ec8bc89c90ecb1ad2e0135fb8abe7a79f2b43badfc73356fd4ae3f1cc3857
// Its analytic temperature inverse supersedes the plan's older iteration prose.
// The caller owns real synchronized ghosts; no local periodic wrapping occurs.
namespace shocks {

struct DetectorConfig {
  Real gamma=1.4, min_mach=1.1, mach_max=1000, grad_floor=1e-30, div_limit=0;
  int width=1;
  // Zero disables these floors. Pressure floor/alternate criteria are unsupported.
  Real density_floor=0, temperature_floor=0;
};

struct Signal {
  bool valid=false, candidate=false, inside=false;
  Real rho=0, pressure=0, temperature=0, divergence=0;
  Real grad_temperature[3]={}, grad_entropy[3]={}, grad_density[3]={};
  Real normal[3]={};
  Real rho_up=0, rho_down=0, pressure_up=0, pressure_down=0;
  Real temperature_up=0, temperature_down=0;
  Real mach_pressure=std::numeric_limits<Real>::quiet_NaN();
  Real mach_temperature=std::numeric_limits<Real>::quiet_NaN();
};

KOKKOS_INLINE_FUNCTION
bool ConfigValid(const DetectorConfig &c) {
  return Kokkos::isfinite(c.gamma) && c.gamma > 1 && c.width >= 1 &&
      Kokkos::isfinite(c.min_mach) && c.min_mach >= 1 &&
      Kokkos::isfinite(c.mach_max) && c.mach_max > 1 &&
      Kokkos::isfinite(c.grad_floor) && c.grad_floor >= 0 &&
      Kokkos::isfinite(c.div_limit) &&
      Kokkos::isfinite(c.density_floor) && c.density_floor >= 0 &&
      Kokkos::isfinite(c.temperature_floor) && c.temperature_floor >= 0;
}

// Read EOS primitives: IEN is internal energy density, not pressure.
KOKKOS_INLINE_FUNCTION
bool ReadState(const DvceArray5D<Real> &w, int m, int k, int j, int i,
               Real gamma, Real &rho, Real &p, Real &t, Real &s) {
  if (m < 0 || m >= static_cast<int>(w.extent(0)) || k < 0 ||
      k >= static_cast<int>(w.extent(2)) || j < 0 ||
      j >= static_cast<int>(w.extent(3)) || i < 0 ||
      i >= static_cast<int>(w.extent(4)) || w.extent(1) < 5) return false;
  for (int n=0; n<5; ++n) {
    if (!Kokkos::isfinite(w(m,n,k,j,i))) return false;
  }
  rho = w(m,IDN,k,j,i);
  p = (gamma-1)*w(m,IEN,k,j,i);
  if (!(rho > 0) || !(p > 0) || !Kokkos::isfinite(p)) return false;
  t = p/rho;
  s = Kokkos::log(p)-gamma*Kokkos::log(rho);
  return Kokkos::isfinite(t) && Kokkos::isfinite(s);
}

KOKKOS_INLINE_FUNCTION
Real Alignment(const Real a[3], const Real b[3]) {
  Real sa=0, sb=0;
  for (int d=0; d<3; ++d) {
    sa = Kokkos::fmax(sa, Kokkos::fabs(a[d]));
    sb = Kokkos::fmax(sb, Kokkos::fabs(b[d]));
  }
  if (sa == 0) sa = 1;
  if (sb == 0) sb = 1;
  return (a[0]/sa)*(b[0]/sb)+(a[1]/sa)*(b[1]/sb)+(a[2]/sa)*(b[2]/sb);
}

KOKKOS_INLINE_FUNCTION
Real TemperatureMach(Real ratio, const DetectorConfig &c) {
  const Real invalid=std::numeric_limits<Real>::quiet_NaN();
  if (!(ratio > 1) || !Kokkos::isfinite(ratio)) return invalid;
  const Real gp=c.gamma+1, gm=c.gamma-1;
  const Real a=2*c.gamma*gm;
  const Real b=-c.gamma*c.gamma+6*c.gamma-1-ratio*gp*gp;
  const Real cc=-2*gm;
  const Real disc=Kokkos::hypot(b, 2*Kokkos::sqrt(a)*Kokkos::sqrt(-cc));
  const Real x=(-0.5*b+0.5*disc)/a;
  if (!(x > 1) || !Kokkos::isfinite(x)) return invalid;
  const Real cap2=c.mach_max*c.mach_max;
  const Real q=(2*c.gamma*cap2-gm)/gp;
  const Real r=gp*cap2/(gm*cap2+2);
  if (ratio > q/r) return invalid;
  const Real mach=Kokkos::sqrt(x);
  return mach > c.mach_max ? c.mach_max : mach;
}

KOKKOS_INLINE_FUNCTION
Signal DetectCell(const DvceArray5D<Real> &w, int m, int k, int j, int i,
                  Real dx, Real dy, Real dz, bool three_d,
                  const DetectorConfig &c) {
  Signal out;
  if (!ConfigValid(c) || !Kokkos::isfinite(dx) || dx <= 0 ||
      !Kokkos::isfinite(dy) || dy <= 0 || !Kokkos::isfinite(dz) || dz <= 0) {
    return out;
  }
  Real entropy;
  if (!ReadState(w,m,k,j,i,c.gamma,out.rho,out.pressure,
                 out.temperature,entropy)) return out;
  const Real spacing[3]={dx,dy,dz};
  Real grad_pressure[3]={};
  for (int d=0; d<(three_d ? 3 : 2); ++d) {
    const int di=(d == 0), dj=(d == 1), dk=(d == 2);
    Real rp, pp, tp, sp, rm, pm, tm, sm;
    if (!ReadState(w,m,k+dk,j+dj,i+di,c.gamma,rp,pp,tp,sp) ||
        !ReadState(w,m,k-dk,j-dj,i-di,c.gamma,rm,pm,tm,sm)) return out;
    out.grad_temperature[d]=(tp-tm)/(2*spacing[d]);
    out.grad_entropy[d]=(sp-sm)/(2*spacing[d]);
    out.grad_density[d]=(rp-rm)/(2*spacing[d]);
    grad_pressure[d]=(pp-pm)/(2*spacing[d]);
    out.divergence += (w(m,IVX+d,k+dk,j+dj,i+di)-
                       w(m,IVX+d,k-dk,j-dj,i-di))/(2*spacing[d]);
    if (!Kokkos::isfinite(out.grad_temperature[d]) ||
        !Kokkos::isfinite(out.grad_entropy[d]) ||
        !Kokkos::isfinite(out.grad_density[d]) ||
        !Kokkos::isfinite(grad_pressure[d])) return out;
  }
  const Real mag=Kokkos::hypot(Kokkos::hypot(out.grad_temperature[0],
      out.grad_temperature[1]),out.grad_temperature[2]);
  // Reference validates pressure gradients even when temperature defines normals.
  const Real pmag=Kokkos::hypot(Kokkos::hypot(grad_pressure[0],
      grad_pressure[1]),grad_pressure[2]);
  if (!Kokkos::isfinite(mag) || !Kokkos::isfinite(pmag) ||
      !Kokkos::isfinite(out.divergence)) return out;
  for (int d=0; d<3; ++d) out.normal[d]=out.grad_temperature[d]/(mag > 0 ? mag : 1);
  out.candidate=out.divergence < c.div_limit && mag > c.grad_floor &&
      Alignment(out.grad_temperature,out.grad_entropy) > 0 &&
      Alignment(out.grad_temperature,out.grad_density) > 0;
  int axis=0;
  if (Kokkos::fabs(out.normal[1]) > Kokkos::fabs(out.normal[axis])) axis=1;
  if (Kokkos::fabs(out.normal[2]) > Kokkos::fabs(out.normal[axis])) axis=2;
  const int shift=(out.normal[axis] < 0 ? -c.width : c.width);
  const int di=(axis == 0 ? shift : 0), dj=(axis == 1 ? shift : 0);
  const int dk=(axis == 2 ? shift : 0);
  Real rp, pp, tp, sp, rm, pm, tm, sm;
  if (!ReadState(w,m,k+dk,j+dj,i+di,c.gamma,rp,pp,tp,sp) ||
      !ReadState(w,m,k-dk,j-dj,i-di,c.gamma,rm,pm,tm,sm)) return out;
  const bool plus_down=pp >= pm;
  out.rho_down=plus_down ? rp : rm;
  out.rho_up=plus_down ? rm : rp;
  out.pressure_down=plus_down ? pp : pm;
  out.pressure_up=plus_down ? pm : pp;
  out.temperature_down=plus_down ? tp : tm;
  out.temperature_up=plus_down ? tm : tp;
  const Real q=out.pressure_down/out.pressure_up;
  const Real r=out.rho_down/out.rho_up;
  const Real t=out.temperature_down/out.temperature_up;
  // Deviation: invalid touched state flags this sample rather than aborting a cube.
  if (!Kokkos::isfinite(q) || !Kokkos::isfinite(r) || !Kokkos::isfinite(t)) return out;
  if (out.candidate && q > 1) {
    out.mach_pressure=Kokkos::sqrt((q*(c.gamma+1)+(c.gamma-1))/(2*c.gamma));
  }
  if (out.candidate && t > 1) out.mach_temperature=TemperatureMach(t,c);
  if (Kokkos::isinf(out.mach_pressure)) return out;
  const Real expected=q/r;
  const bool triplet=r <= (c.gamma+1)/(c.gamma-1) &&
      Kokkos::fabs(t-expected) <= 0.25*Kokkos::fabs(expected);
  const bool floors=out.rho_up >= c.density_floor &&
      out.temperature_up >= c.temperature_floor;
  out.valid=true;
  out.inside=out.candidate && q > 1 && r > 1 && t > 1 && triplet && floors &&
      (out.mach_pressure >= c.min_mach || out.mach_temperature >= c.min_mach);
  return out;
}

}  // namespace shocks
#endif  // SHOCKS_DETECTOR_HPP_
