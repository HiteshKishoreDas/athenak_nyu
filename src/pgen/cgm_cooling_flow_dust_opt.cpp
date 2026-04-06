#include <iostream>
#include <cmath>
#include <stdexcept>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "pgen.hpp"
#include "units/units.hpp"
#include "utils/profile_reader.hpp"
#include "utils/sn_scheduler.hpp"
#include "utils/random.hpp"
#include "particles/particles.hpp"

constexpr int kMaxDustBins = 32;
constexpr int kMaxDustBinSlots = kMaxDustBins + 2;
constexpr int kMaxDustEdgeSlots = kMaxDustBins + 3;

struct Alims {
  Real amin;
  Real amax;
};

KOKKOS_INLINE_FUNCTION
void SetEquilibriumState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real x1l, Real x1r, Real x2l, Real x2r,
                         Real G, Real r_s, Real rho_s, Real m_g,
                         Real a_g, Real z_g, Real r_m, Real rho_m,
                         Real gm1, const ProfileReader &disk_profile);

KOKKOS_INLINE_FUNCTION
void SetCoolingFlowState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real gm1, const ProfileReader &profile);

KOKKOS_INLINE_FUNCTION
void SetRotation(const DvceArray5D<Real> &u0,
                 int m, int k, int j, int i,
                 Real x1v, Real x2v, Real x3v,
                 Real r_circ, Real v_circ);

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3,
             Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal,
             Real c_outer, Real rho_mean);

void UserSource(Mesh* pm, const Real bdt);
void GravitySource(Mesh* pm, const Real bdt);
void SNSource(Mesh* pm, const Real bdt);
void DensityCeilingSource(Mesh* pm, const Real bdt);
void UserBoundary(Mesh* pm);
void FreeProfile(ParameterInput *pin, Mesh *pm);
void RefinementCondition(MeshBlockPack* pmbp);
void AddDustSource(Mesh *pm, const Real bdt);
void UpdateTurbulenceCache(Mesh *pm);
void CGMHistory(HistoryData *pdata, Mesh *pm);

KOKKOS_INLINE_FUNCTION
void SetScalarState(const DvceArray5D<Real> &u0, int m, int k, int j, int i,
                    int nhydro, int nscalars, Real Z_init, Real Zsol,
                    bool dust_enabled, int ndust, Real D_Z,
                    const DvceArray1D<Real> &dust_fracs) {
  const int nmetal = nhydro;
  const int ndusti = nhydro + 1;
  Real rho = u0(m, IDN, k, j, i);
  if (nscalars > 0) {
    u0(m, nmetal, k, j, i) = Z_init * Zsol * rho;
  }
  for (int n = ndusti; n < nhydro + nscalars; ++n) {
    u0(m, n, k, j, i) = 0.0;
  }
  if (dust_enabled) {
    Real dust_mass = D_Z * Z_init * Zsol * rho;
    for (int id = 0; id < ndust; ++id) {
      u0(m, ndusti + id, k, j, i) = dust_mass * dust_fracs(id);
    }
  }
}

template <typename T>
KOKKOS_INLINE_FUNCTION
Real array_get(const T &a, int i) {
  return a(i);
}

KOKKOS_INLINE_FUNCTION
Real array_get(const Real *a, int i) {
  return a[i];
}

KOKKOS_INLINE_FUNCTION
Real array_get(Real *a, int i) {
  return a[i];
}

template <typename ArrayType, typename T>
KOKKOS_INLINE_FUNCTION
int searchsorted(const ArrayType &a, const int n, const T &x) {
  if (n <= 0) return 0;
  int lo = 0;
  int hi = n;

  if (x <= array_get(a, lo)) return 0;
  if (x > array_get(a, n - 1)) return hi;

  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (array_get(a, mid) < x) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

template <typename EdgeArray>
KOKKOS_INLINE_FUNCTION
int searchsorted_log_edges(const EdgeArray &edges, const int n, const Real x,
                           const Real inv_log_step) {
  if (n <= 0) return 0;
  const Real first = array_get(edges, 0);
  const Real last = array_get(edges, n - 1);
  if (x <= first) return 0;
  if (x > last) return n;

  int idx = static_cast<int>(ceil(log(x / first) * inv_log_step));
  if (idx < 0) idx = 0;
  if (idx > (n - 1)) idx = n - 1;
  return idx;
}

template <typename EdgeArray>
KOKKOS_INLINE_FUNCTION
int searchsorted_padded(const EdgeArray &edges, const int n, const Real x,
                        const Real inv_log_step) {
  if (n <= 0) return 0;
  const Real left_ghost = array_get(edges, 0);
  const Real first_physical = array_get(edges, 1);
  const Real last_physical = array_get(edges, n - 2);
  const Real right_ghost = array_get(edges, n - 1);

  if (x <= left_ghost) return 0;
  if (x > right_ghost) return n;
  if (x <= first_physical) return 1;
  if (x > last_physical) return n - 1;

  int idx = 1 + static_cast<int>(ceil(log(x / first_physical) * inv_log_step));
  if (idx < 1) idx = 1;
  if (idx > (n - 2)) idx = n - 2;
  return idx;
}

KOKKOS_INLINE_FUNCTION
Real square(const Real x) {
  return x * x;
}

KOKKOS_INLINE_FUNCTION
Real cube(const Real x) {
  return x * x * x;
}

KOKKOS_INLINE_FUNCTION
Real quartic(const Real x) {
  const Real x2 = x * x;
  return x2 * x2;
}

KOKKOS_INLINE_FUNCTION
Real dust_k_um(const Real rho_grain) {
  return (4.0 / 3.0) * M_PI * (rho_grain * 1.0e-12);
}

KOKKOS_INLINE_FUNCTION
Real mass_to_num(Real mass_in_bin, Real edges1, Real edges2, const Real rho_grain) {
  const Real k_dust_um = dust_k_um(rho_grain);
  Real dist_change = 4.0 * mass_in_bin / k_dust_um;
  dist_change /= quartic(edges2) - quartic(edges1);
  return dist_change;
}

KOKKOS_INLINE_FUNCTION
Real mass_to_num(Real mass_in_bin, Real edge4_delta, const Real inv_k_dust_um) {
  return mass_in_bin * inv_k_dust_um / edge4_delta;
}

template <typename DistArray, typename EdgeArray>
KOKKOS_INLINE_FUNCTION
Real dist_num_integral(const DistArray &d_arr, const EdgeArray &edges, const int n,
                       Real a1, Real a2, Real rho_grain, Real dust_inv_log_step,
                       bool mass = false) {
  int i1 = searchsorted_padded(edges, n, a1, dust_inv_log_step) - 1;
  int i2 = searchsorted_padded(edges, n, a2, dust_inv_log_step) - 1;

  const Real k_dust_um = dust_k_um(rho_grain);

  Real integral = 0.0;
  if (i1 == i2) {
    if (mass) {
      integral += 0.25 * array_get(d_arr, i1) * k_dust_um *
                  (quartic(a2) - quartic(a1));
    } else {
      integral += array_get(d_arr, i1) * (a2 - a1);
    }
    return integral;
  }

  if (mass) {
    integral += 0.25 * array_get(d_arr, i1) * k_dust_um *
                (quartic(array_get(edges, i1 + 1)) - quartic(a1));
  } else {
    integral += array_get(d_arr, i1) * (array_get(edges, i1 + 1) - a1);
  }

  if (mass) {
    integral += 0.25 * array_get(d_arr, i2) * k_dust_um *
                (quartic(a2) - quartic(array_get(edges, i2)));
  } else {
    integral += array_get(d_arr, i2) * (a2 - array_get(edges, i2));
  }

  if (i2 > i1 + 1) {
    if (mass) {
      for (int i = i1 + 1; i < i2; ++i) {
        integral += 0.25 * array_get(d_arr, i) * k_dust_um *
                    (quartic(array_get(edges, i + 1)) - quartic(array_get(edges, i)));
      }
    } else {
      for (int i = i1 + 1; i < i2; ++i) {
        integral += array_get(d_arr, i) *
                    (array_get(edges, i + 1) - array_get(edges, i));
      }
    }
  }

  return integral;
}

KOKKOS_INLINE_FUNCTION
Real sputtering(Real n_H, Real T, Real Z, Real grain_porosity_factor) {
  Real da_dt = -1.0;
  da_dt *= n_H / 1.0;
  const Real tratio = 2.0e6 / T;
  da_dt /= 1.0 + tratio * tratio * sqrt(tratio);
  da_dt *= Z;
  da_dt *= grain_porosity_factor;
  return da_dt;
}

KOKKOS_INLINE_FUNCTION
Real accretion(Real n_H, Real T, Real Z, Real grain_porosity_factor) {
  Real da_dt = 1.8622e-4;
  da_dt *= n_H / 1.0e3;
  da_dt /= sqrt(T / 10.0);
  da_dt *= Z;
  da_dt *= grain_porosity_factor;
  return da_dt;
}

template <typename EdgeArray>
KOKKOS_INLINE_FUNCTION
void rebin(Real *d_arr, const EdgeArray &edges, const int n, Real shift, Real rho_grain) {
  if (n > kMaxDustEdgeSlots) return;

  const int nbin = n - 1;
  Real mass_in_bin[kMaxDustBinSlots];
  const Real k_dust_um = dust_k_um(rho_grain);
  const Real inv_k_dust_um = 4.0 / k_dust_um;
  for (int i = 0; i < nbin; ++i) mass_in_bin[i] = 0.0;

  int src = 0;
  int dst = 0;
  while (src < nbin && dst < nbin) {
    const Real src_left = array_get(edges, src) + shift;
    const Real src_right = array_get(edges, src + 1) + shift;
    const Real dst_left = array_get(edges, dst);
    const Real dst_right = array_get(edges, dst + 1);
    const Real left = Kokkos::max(src_left, dst_left);
    const Real right = Kokkos::min(src_right, dst_right);
    if (right > left) {
      mass_in_bin[dst] += 0.25 * d_arr[src] * k_dust_um *
                          (quartic(right) - quartic(left));
    }
    if (src_right <= dst_right) {
      ++src;
    } else {
      ++dst;
    }
  }

  for (int i = 1; i < (nbin - 1); ++i) {
    d_arr[i] = mass_to_num(mass_in_bin[i],
                           quartic(array_get(edges, i + 1)) - quartic(array_get(edges, i)),
                           inv_k_dust_um);
  }
  d_arr[0] = mass_in_bin[0];
  d_arr[nbin - 1] = mass_in_bin[nbin - 1];
}

template <typename StateArray, typename FactorArray, typename EdgeArray>
KOKKOS_INLINE_FUNCTION
Real apply_sputter_accrete_rebin(Real *dust_arr, const StateArray &u0,
                                 int m, int k, int j, int i, int ndusti,
                                 int ndust, const FactorArray &mass_to_num_fac,
                                 const EdgeArray &edges_pad, int n_edges,
                                 Real dens, Real temp, Real Z_local, Real Z_solar,
                                 Real grain_porosity_factor, Real rho_grain,
                                 Real dust_a_local, Real dust_b_local,
                                 Real dust_inv_log_step, Real bdt, Real MYR) {
  constexpr Real D_tol = 1.0e-20;

  dust_arr[0] = 0.0;
  dust_arr[ndust + 1] = 0.0;

  Real D_check = 0.0;
  for (int id = 0; id < ndust; ++id) {
    const Real dust_mass = u0(m, ndusti + id, k, j, i);
    dust_arr[id + 1] = dust_mass * mass_to_num_fac(id);
    D_check += dust_mass;
  }
  if (D_check <= D_tol) return 0.0;

  const Real dust_tot_before = dist_num_integral(
      dust_arr, edges_pad, n_edges, dust_a_local, dust_b_local, rho_grain,
      dust_inv_log_step, true);

  Real tot_shift = sputtering(dens, temp, Z_local / Z_solar, grain_porosity_factor);
  tot_shift += accretion(dens, temp, Z_local / Z_solar, grain_porosity_factor);
  tot_shift *= bdt / MYR;

  rebin(dust_arr, edges_pad, n_edges, tot_shift, rho_grain);
  return dust_tot_before;
}

template <typename StateArray, typename FactorArray>
KOKKOS_INLINE_FUNCTION
void finalize_dust_source(const StateArray &u0, Real *dust_arr,
                          int m, int k, int j, int i, int nmetal, int ndusti,
                          int ndust, Real dens, Real dust_tot, Real rate_Z,
                          const FactorArray &num_to_mass_fac) {
  constexpr Real D_tol = 1.0e-20;

  for (int id = 0; id < ndust; ++id) {
    u0(m, ndusti + id, k, j, i) = dust_arr[id + 1] * num_to_mass_fac(id);
  }

  u0(m, nmetal, k, j, i) += rate_Z;
  u0(m, nmetal, k, j, i) = Kokkos::clamp(u0(m, nmetal, k, j, i), 0.0, dens);

  const Real clamp_dust_tot = Kokkos::clamp(dust_tot, 0.0, dens);
  for (int id = 0; id < ndust; ++id) {
    if (dust_tot <= D_tol) {
      u0(m, ndusti + id, k, j, i) = 0.0;
    } else {
      u0(m, ndusti + id, k, j, i) *= clamp_dust_tot / dust_tot;
    }
  }
}

KOKKOS_INLINE_FUNCTION
Real v_grain(Real a, Real M_g, Real n_H, Real T, Real rho_grain) {
  const Real mach_32 = M_g * sqrt(M_g);
  const Real temp_quarter = sqrt(sqrt(T / 1.0e4));
  const Real nH_minus_quarter = 1.0 / sqrt(sqrt(n_H));
  return 1.1 * mach_32 * sqrt(a / 0.1) * temp_quarter * nH_minus_quarter *
         sqrt(rho_grain / 3.5);
}

KOKKOS_INLINE_FUNCTION
Real maxwell_tail_mean(Real v, Real vgr) {
  constexpr Real vgr_tol = 1.0e-20;
  if (vgr <= vgr_tol) return 0.0;

  Real v0 = vgr * sqrt(2.0 / 3.0);
  Real vmean = sqrt(8.0 / M_PI) * v0;
  const Real vv0 = v / v0;
  vmean *= 1.0 + 0.5 * square(vv0);
  vmean *= exp(-0.5 * square(vv0));
  return vmean;
}

KOKKOS_INLINE_FUNCTION
Real maxwell_head_mean(Real v, Real vgr) {
  constexpr Real vgr_tol = 1.0e-20;
  if (vgr <= vgr_tol) return 0.0;

  Real v0 = vgr * sqrt(2.0 / 3.0);
  Real vmean = sqrt(8.0 / M_PI) * v0;
  const Real vv0 = v / v0;
  vmean *= 1.0 + 0.5 * square(vv0);
  vmean *= exp(-0.5 * square(vv0));
  return sqrt(8.0 / M_PI) * v0 - vmean;
}

KOKKOS_INLINE_FUNCTION
Real interaction_rate_pair(Real bdt, Real dist_i, Real dist_j, Real dbin_i, Real dbin_j,
                           Real dbin_i3, Real dbin_j3, Real delta1i, Real delta2i,
                           Real delta3i, Real delta1j, Real delta2j, Real delta3j,
                           Real vgrain_i, Real vol_cc, bool shatter) {
  constexpr Real v_shatt = 0.5 * (1.2 + 2.7);

  Real vproc = 0.0;
  if (shatter) {
    vproc = maxwell_tail_mean(v_shatt, vgrain_i);
  } else {
    const Real pair_size = dbin_i + dbin_j;
    const Real pair_geom = (dbin_i3 + dbin_j3) / cube(pair_size);
    vproc = sqrt(pair_geom) *
            pow(pair_size / (dbin_i * dbin_j), 5.0 / 6.0);
    vproc = maxwell_head_mean(vproc, vgrain_i);
  }

  Real rate = delta3i * delta1j / 3.0 + delta2i * delta2j / 2.0 +
              delta1i * delta3j / 3.0;
  rate *= M_PI * vproc / vol_cc;
  rate *= dist_i * dist_j;
  rate *= 1.0e-3 * bdt;
  return rate;
}

KOKKOS_INLINE_FUNCTION
void calc_interaction(Real bdt, Real intr_arr[][kMaxDustBins], const Real dist[],
                      const Real dbins[], const Real edges[], int nbin,
                      Real M_g, Real rho_gr, Real n_H, Real T, Real vol_cc,
                      bool shatter) {
  const Real F_stick = 10.0;
  const Real gamma_Si = 2.7;
  const Real gamma_C = 1.2;
  const Real gamma_d = 0.5 * (gamma_C + gamma_Si);
  const Real E_Si = 5.4e11;
  const Real E_C = 3.4e10;
  const Real E_d = 0.5 * (E_C + E_Si);
  const Real v_shatt_Si = 2.7;
  const Real v_shatt_C = 1.2;
  const Real v_shatt = 0.5 * (v_shatt_C + v_shatt_Si);

  Real vc = 2.14;
  vc *= F_stick * std::pow(gamma_d, 5.0 / 6.0) /
        std::pow(E_d, 1.0 / 3.0) / std::sqrt(rho_gr);

  for (int i = 0; i < nbin; ++i) {
    const Real dela1i = edges[i + 1] - edges[i];
    const Real dela2i = std::pow(edges[i + 1], 2.0) - std::pow(edges[i], 2.0);
    const Real dela3i = std::pow(edges[i + 1], 3.0) - std::pow(edges[i], 3.0);

    for (int j = 0; j < nbin; ++j) {
      const Real dela1j = edges[j + 1] - edges[j];
      const Real dela2j = std::pow(edges[j + 1], 2.0) - std::pow(edges[j], 2.0);
      const Real dela3j = std::pow(edges[j + 1], 3.0) - std::pow(edges[j], 3.0);

      Real vproc = 0.0;
      if (shatter) {
        vproc = maxwell_tail_mean(v_shatt,
                                  v_grain(dbins[i], M_g, n_H, T, rho_gr));
      } else {
        vproc = (std::pow(dbins[i], 3.0) + std::pow(dbins[j], 3.0)) /
                std::pow(dbins[i] + dbins[j], 3.0);
        vproc = std::sqrt(vproc) *
                std::pow((dbins[i] + dbins[j]) / dbins[i] / dbins[j], 5.0 / 6.0);
        vproc = maxwell_head_mean(vproc,
                                  v_grain(dbins[i], M_g, n_H, T, rho_gr));
      }

      intr_arr[i][j] = dela3i * dela1j / 3.0 + dela2i * dela2j / 2.0 +
                       dela1i * dela3j / 3.0;
      intr_arr[i][j] *= M_PI * vproc / vol_cc;
      intr_arr[i][j] *= dist[i] * dist[j];
      intr_arr[i][j] *= 1.0e-3 * bdt;
    }
  }
}

KOKKOS_INLINE_FUNCTION
void add_coagulate(Real intr_arr[][kMaxDustBins], Real mass_in_bin[],
                   const Real dbins[], const Real edges[], int nbin,
                   Real rho_gr) {
  const Real rho_gr_um = rho_gr * 1.0e-12;
  const Real k_dust_um = (4.0 / 3.0) * M_PI * rho_gr_um;

  for (int i = 0; i < nbin; ++i) {
    for (int j = 0; j < nbin; ++j) {
      const Real ai3 = std::pow(dbins[i], 3.0);
      const Real aj3 = std::pow(dbins[j], 3.0);
      const Real a_coag = std::pow(ai3 + aj3, 1.0 / 3.0);
      const int k = searchsorted(edges, nbin + 1, a_coag) - 1;

      if ((k < nbin) && (k > 0)) {
        mass_in_bin[i] -= intr_arr[i][j] * k_dust_um * ai3;
        mass_in_bin[j] -= intr_arr[i][j] * k_dust_um * aj3;
        mass_in_bin[k] += intr_arr[i][j] * k_dust_um * (ai3 + aj3);
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION
Real sigma_fn(Real M, const Real s) {
  const Real M_inv = 1.0 / M;
  return 0.3 * std::pow(s + M_inv - 0.11, 0.13) / (s + M_inv - 1.0);
}

KOKKOS_INLINE_FUNCTION
Real M_shocked(Real Mrel, Real Mproj, const Real rho_gr, const Real sigma_1,
               const Real s, const Real R, const Real M_1) {
  constexpr Real Mrel_tol = 1.0e-6;
  if (Mrel <= Mrel_tol) return 0.0;

  const Real sigma_r = sigma_fn(Mrel / (1.0 + R), s);
  Real Msh = (1.0 + 2.0 * R) / std::pow(1.0 + R, 9.0 / 16.0) /
             std::pow(sigma_r, 1.0 / 9.0);
  Msh *= std::pow(Mrel / sigma_1 / M_1, 8.0 / 9.0);
  return Msh * Mproj;
}

KOKKOS_INLINE_FUNCTION
Alims a_lim_frag(Real Mfrac, const Real P1, const Real rho_gr, const Real Pv,
                 const Real z_const) {
  const Real rho_gr_um = rho_gr * 1.0e-12;
  Alims afrac_lim{0.0, 0.0};

  afrac_lim.amax = 3.0 / 4.0 / M_PI / rho_gr_um * (1.0 + z_const);
  afrac_lim.amax /= 4.0 * std::pow(z_const, 3.0) * (z_const - 2.0);
  afrac_lim.amax *= Mfrac;
  afrac_lim.amin = afrac_lim.amax * std::pow(P1 / Pv, 1.47);
  return afrac_lim;
}

KOKKOS_INLINE_FUNCTION
Real integrate_frag(Real Anorm, Real a1, Real a2, Real expo, const Real rho_gr) {
  const Real rho_gr_um = rho_gr * 1.0e-12;
  const Real k_dust_um = (4.0 / 3.0) * M_PI * rho_gr_um;
  const Real mex = expo + 4.0;
  return Anorm * k_dust_um * (std::pow(a2, mex) - std::pow(a1, mex)) / mex;
}

template <typename EdgeArray>
KOKKOS_INLINE_FUNCTION
void deposit_frag(Real mass_in_bin[], const EdgeArray &edges, int nbin, Alims alim,
                  Real Anorm, Real expo, const Real rho_gr,
                  const Real dust_inv_log_step) {
  int k = searchsorted_log_edges(edges, nbin + 1, alim.amin, dust_inv_log_step) - 1;
  int l = searchsorted_log_edges(edges, nbin + 1, alim.amax, dust_inv_log_step) - 1;

  if (l >= nbin) {
    Kokkos::abort("deposit_frag overflowed dust-bin bounds");
  }

  if (l < 0) {
    mass_in_bin[0] += integrate_frag(Anorm, alim.amin, alim.amax, expo, rho_gr);
    return;
  }

  Real aleft = alim.amin;
  if (k < 0) {
    mass_in_bin[0] += integrate_frag(Anorm, alim.amin, array_get(edges, 0), expo, rho_gr);
    k = 0;
    aleft = array_get(edges, 0);
  }

  if (k == l) {
    mass_in_bin[k] += integrate_frag(Anorm, aleft, alim.amax, expo, rho_gr);
    return;
  }

  mass_in_bin[k] += integrate_frag(Anorm, aleft, array_get(edges, k + 1), expo, rho_gr);
  for (int i = k + 1; i < l; ++i) {
    mass_in_bin[i] += integrate_frag(Anorm, array_get(edges, i),
                                     array_get(edges, i + 1), expo, rho_gr);
  }
  mass_in_bin[l] += integrate_frag(Anorm, array_get(edges, l), alim.amax, expo, rho_gr);
}

KOKKOS_INLINE_FUNCTION
void add_shatters(Real intr_arr[][kMaxDustBins], Real mass_in_bin[],
                  const Real dbins[], const Real edges[], int nbin, Real M_g,
                  const Real rho_gr, Real n_H, Real T) {
  constexpr Real R = 1.0;
  constexpr Real P1_Si = 3.0e11;
  constexpr Real P1_C = 4.0e10;
  constexpr Real P1 = 0.5 * (P1_Si + P1_C);
  constexpr Real Pv_Si = 5.4e12;
  constexpr Real Pv_C = 5.8e12;
  constexpr Real Pv = 0.5 * (Pv_Si + Pv_C);
  constexpr Real v_shatt_Si = 2.7;
  constexpr Real v_shatt_C = 1.2;
  constexpr Real v_shatt = 0.5 * (v_shatt_Si + v_shatt_C);
  constexpr Real s_Si = 1.2;
  constexpr Real s_C = 1.9;
  constexpr Real s = 0.5 * (s_Si + s_C);
  constexpr Real z_const = 3.4;
  constexpr Real c0_Si = 5.0;
  constexpr Real c0_C = 1.8;
  constexpr Real c0 = 0.5 * (c0_Si + c0_C);
  constexpr Real c0_cgs = c0 * 1.0e5;
  constexpr Real shatter_expo = -3.3;
  constexpr Real mex = shatter_expo + 4.0;

  const Real phi_1 = P1 / rho_gr / c0_cgs / c0_cgs;
  const Real M_1 = 2.0 * phi_1 / (1.0 + std::sqrt(1.0 + 4.0 * s * phi_1));
  const Real sigma_1 = sigma_fn(M_1, s);
  const Real rho_gr_um = rho_gr * 1.0e-12;
  const Real k_dust_um = (4.0 / 3.0) * M_PI * rho_gr_um;
  const Real dust_inv_log_step = 1.0 / std::log(edges[1] / edges[0]);
  constexpr Real Mfrag_tol = 1.0e-20;

  for (int i = 0; i < nbin; ++i) {
    for (int j = 0; j < nbin; ++j) {
      const Real ai3 = std::pow(dbins[i], 3.0);
      const Real aj3 = std::pow(dbins[j], 3.0);
      const Real Mi = k_dust_um * ai3;
      const Real Mj = k_dust_um * aj3;

      const Real v_grain_i = v_grain(dbins[i], M_g, n_H, T, rho_gr);
      const Real vrel = maxwell_tail_mean(v_shatt, v_grain_i);
      const Real Msh = M_shocked(vrel / c0, Mj, rho_gr, sigma_1, s, R, M_1);

      const Real Mfrag = (Msh > 0.5 * Mi) ? Mi : 0.4 * Msh;
      const Real Mleft = Mi - Mfrag;

      if ((Mfrag / Mi) >= Mfrag_tol) {
        const Real a_left = std::pow(Mleft / k_dust_um, 1.0 / 3.0);
        const int k = searchsorted_log_edges(edges, nbin + 1, a_left,
                                             dust_inv_log_step) - 1;
        if ((k < nbin) && (k > 0)) {
          mass_in_bin[i] -= intr_arr[i][j] * Mi;
          mass_in_bin[j] -= intr_arr[i][j] * Mj;
          mass_in_bin[k] += 2.0 * intr_arr[i][j] * Mleft;
        }

        const Alims frag_alim = a_lim_frag(2.0 * Mfrag, P1, rho_gr, Pv, z_const);
        const Real frag_norm = Mfrag * mex / k_dust_um /
                               (std::pow(frag_alim.amax, mex) -
                                std::pow(frag_alim.amin, mex));
        deposit_frag(mass_in_bin, edges, nbin, frag_alim, frag_norm,
                     shatter_expo, rho_gr, dust_inv_log_step);
      }
    }
  }
}

template <typename DistArray, typename BinArray, typename Bin3Array,
          typename EdgeArray, typename DeltaArray>
KOKKOS_INLINE_FUNCTION
void add_coagulate_streamed(Real bdt, Real mass_in_bin[], const DistArray &dist,
                            const BinArray &dbins, const Bin3Array &dbins3,
                            const EdgeArray &edges, const DeltaArray &delta1,
                            const DeltaArray &delta2, const DeltaArray &delta3,
                            int nbin, Real M_g, Real rho_grain, Real n_H, Real T,
                            Real vol_cc, Real dust_inv_log_step) {
  const Real k_dust_um = dust_k_um(rho_grain);

  for (int i = 0; i < nbin; ++i) {
    const Real dbin_i = array_get(dbins, i);
    const Real ai3 = array_get(dbins3, i);
    const Real Mi = k_dust_um * ai3;
    const Real vgrain_i = v_grain(dbin_i, M_g, n_H, T, rho_grain);
    for (int j = 0; j < nbin; ++j) {
      const Real dbin_j = array_get(dbins, j);
      const Real aj3 = array_get(dbins3, j);
      const Real rate = interaction_rate_pair(
          bdt, array_get(dist, i), array_get(dist, j), dbin_i, dbin_j, ai3, aj3,
          array_get(delta1, i), array_get(delta2, i), array_get(delta3, i),
          array_get(delta1, j), array_get(delta2, j), array_get(delta3, j),
          vgrain_i, vol_cc, false);
      const Real a_coag = pow(ai3 + aj3, 1.0 / 3.0);
      const int k = searchsorted_log_edges(edges, nbin + 1, a_coag,
                                           dust_inv_log_step) - 1;
      if ((k < nbin) && (k > 0)) {
        mass_in_bin[i] -= rate * Mi;
        mass_in_bin[j] -= rate * k_dust_um * aj3;
        mass_in_bin[k] += rate * k_dust_um * (ai3 + aj3);
      }
    }
  }
}

template <typename DistArray, typename BinArray, typename Bin3Array,
          typename EdgeArray, typename DeltaArray>
KOKKOS_INLINE_FUNCTION
void add_shatters_streamed(Real bdt, Real mass_in_bin[], const DistArray &dist,
                           const BinArray &dbins, const Bin3Array &dbins3,
                           const EdgeArray &edges, const DeltaArray &delta1,
                           const DeltaArray &delta2, const DeltaArray &delta3,
                           int nbin, Real M_g, const Real rho_grain,
                           Real n_H, Real T, Real vol_cc,
                           Real dust_inv_log_step) {
  constexpr Real R = 1.0;
  constexpr Real P1_Si = 3.0e11;
  constexpr Real P1_C = 4.0e10;
  constexpr Real P1 = 0.5 * (P1_Si + P1_C);
  constexpr Real Pv_Si = 5.4e12;
  constexpr Real Pv_C = 5.8e12;
  constexpr Real Pv = 0.5 * (Pv_Si + Pv_C);
  constexpr Real v_shatt_Si = 2.7;
  constexpr Real v_shatt_C = 1.2;
  constexpr Real v_shatt = 0.5 * (v_shatt_C + v_shatt_Si);
  constexpr Real s_Si = 1.2;
  constexpr Real s_C = 1.9;
  constexpr Real s = 0.5 * (s_Si + s_C);
  constexpr Real z_const = 3.4;
  constexpr Real c0_Si = 5.0;
  constexpr Real c0_C = 1.8;
  constexpr Real c0 = 0.5 * (c0_Si + c0_C);
  constexpr Real c0_cgs = c0 * 1.0e5;
  constexpr Real shatter_expo = -3.3;
  constexpr Real mex = shatter_expo + 4.0;

  const Real phi_1 = P1 / rho_grain / c0_cgs / c0_cgs;
  const Real M_1 = 2.0 * phi_1 / (1.0 + sqrt(1.0 + 4.0 * s * phi_1));
  const Real sigma_1 = sigma_fn(M_1, s);
  const Real k_dust_um = dust_k_um(rho_grain);
  constexpr Real Mfrag_tol = 1.0e-20;

  for (int i = 0; i < nbin; ++i) {
    const Real dbin_i = array_get(dbins, i);
    const Real ai3 = array_get(dbins3, i);
    const Real Mi = k_dust_um * ai3;
    const Real v_grain_i = v_grain(dbin_i, M_g, n_H, T, rho_grain);
    for (int j = 0; j < nbin; ++j) {
      const Real dbin_j = array_get(dbins, j);
      const Real aj3 = array_get(dbins3, j);
      const Real Mj = k_dust_um * aj3;

      const Real rate = interaction_rate_pair(
          bdt, array_get(dist, i), array_get(dist, j), dbin_i, dbin_j, ai3, aj3,
          array_get(delta1, i), array_get(delta2, i), array_get(delta3, i),
          array_get(delta1, j), array_get(delta2, j), array_get(delta3, j),
          v_grain_i, vol_cc, true);
      Real vrel = maxwell_tail_mean(v_shatt, v_grain_i);
      Real Msh = M_shocked(vrel / c0, Mj, rho_grain, sigma_1, s, R, M_1);

      Real Mfrag = (Msh > 0.5 * Mi) ? Mi : 0.4 * Msh;
      Real Mleft = Mi - Mfrag;

      if ((Mfrag / Mi) >= Mfrag_tol) {
        Real a_left = pow(Mleft / k_dust_um, 1.0 / 3.0);
        int k = searchsorted_log_edges(edges, nbin + 1, a_left,
                                       dust_inv_log_step) - 1;
        if ((k < nbin) && (k > 0)) {
          mass_in_bin[i] -= rate * Mi;
          mass_in_bin[j] -= rate * Mj;
          mass_in_bin[k] += 2.0 * rate * Mleft;
        }

        Alims frag_alim = a_lim_frag(2.0 * Mfrag, P1, rho_grain, Pv, z_const);
        Real frag_norm = Mfrag * mex / k_dust_um /
                         (pow(frag_alim.amax, mex) - pow(frag_alim.amin, mex));
        deposit_frag(mass_in_bin, edges, nbin, frag_alim, frag_norm,
                     shatter_expo, rho_grain, dust_inv_log_step);
      }
    }
  }
}

namespace {
  // Constants for gravitational potential
  Real r_scale;
  Real rho_scale;
  Real m_gal;
  Real a_gal;
  Real z_gal;
  Real r_200;
  Real rho_mean;

  // Constants for rotation
  Real r_circ;
  Real v_circ;

  // Constant for metallicity
  Real Z;


  // Constants for SN
  Real r_inj;
  Real e_sn;
  Real m_ej;
  Real sn_ejecta_Z;
  Real sn_ejecta_dust_ratio;

  bool dust_model = false;
  Real D_Z_init = 0.0;
  Real dust_a = 0.0;
  Real dust_b = 0.0;
  Real dust_exp = -3.5;
  Real dust_inv_log_step = 0.0;
  Real rho_gr = 3.5;
  Real grain_porosity = 1.0;
  Real Z_solar = 0.02;
  int Ndust_bins = 0;
  int vturb_est_ncell = 1;
  bool coag_shatt_flag = false;
  Real grain_porosity_factor = 1.0;
  Real vturb_scale = 0.0;
  DvceArray1D<Real> dust_bin_centers;
  DvceArray1D<Real> dust_bin_edges;
  DvceArray1D<Real> dust_bin_edges_padded;
  DvceArray1D<Real> dust_bin_centers_cubed;
  DvceArray1D<Real> dust_delta1;
  DvceArray1D<Real> dust_delta2;
  DvceArray1D<Real> dust_delta3;
  DvceArray1D<Real> dust_delta4;
  DvceArray1D<Real> dust_mass_to_num;
  DvceArray1D<Real> dust_num_to_mass;
  DvceArray1D<Real> dust_bin_mass_fractions;
  HostArray1D<Real> dust_bin_centers_h;
  HostArray1D<Real> dust_bin_edges_h;
  HostArray1D<Real> dust_bin_edges_padded_h;
  HostArray1D<Real> dust_bin_centers_cubed_h;
  HostArray1D<Real> dust_edge4_h;
  HostArray1D<Real> dust_delta1_h;
  HostArray1D<Real> dust_delta2_h;
  HostArray1D<Real> dust_delta3_h;
  HostArray1D<Real> dust_delta4_h;
  HostArray1D<Real> dust_mass_to_num_h;
  HostArray1D<Real> dust_num_to_mass_h;
  DvceArray4D<Real> mach_cache;
  int mach_cache_cycle = -1;

  // Profiles
  ProfileReader profile_reader;
  ProfileReader disk_profile_reader;

  // Refinment condition threshold
  Real ddens_threshold;

  // SN injection persistent buffer
  DvceArray2D<Real> sn_centers_buffer;
  Kokkos::View<int> sn_counter;

  // SN injection flags
  int last_sn_detect_cycle = -1;
  int num_sn_this_cycle = 0;
}
//===========================================================================//
//                               Initialize                                  //
//===========================================================================//

void InitializeDustBins() {
  DvceArray1D<Real> centers("dust_bin_centers", Ndust_bins);
  DvceArray1D<Real> edges("dust_bin_edges", Ndust_bins + 1);
  DvceArray1D<Real> edges_padded("dust_bin_edges_padded", Ndust_bins + 3);
  DvceArray1D<Real> centers_cubed("dust_bin_centers_cubed", Ndust_bins);
  DvceArray1D<Real> delta1("dust_delta1", Ndust_bins);
  DvceArray1D<Real> delta2("dust_delta2", Ndust_bins);
  DvceArray1D<Real> delta3("dust_delta3", Ndust_bins);
  DvceArray1D<Real> delta4("dust_delta4", Ndust_bins);
  DvceArray1D<Real> mass_to_num_fac("dust_mass_to_num", Ndust_bins);
  DvceArray1D<Real> num_to_mass_fac("dust_num_to_mass", Ndust_bins);
  DvceArray1D<Real> mass_fractions("dust_bin_mass_fractions", Ndust_bins);

  auto centers_h = Kokkos::create_mirror_view(centers);
  auto edges_h = Kokkos::create_mirror_view(edges);
  auto edges_padded_h = Kokkos::create_mirror_view(edges_padded);
  HostArray1D<Real> centers_cubed_h("dust_bin_centers_cubed_h", Ndust_bins);
  HostArray1D<Real> edge4_h("dust_edge4_h", Ndust_bins + 3);
  HostArray1D<Real> delta1_h("dust_delta1_h", Ndust_bins);
  HostArray1D<Real> delta2_h("dust_delta2_h", Ndust_bins);
  HostArray1D<Real> delta3_h("dust_delta3_h", Ndust_bins);
  HostArray1D<Real> delta4_h("dust_delta4_h", Ndust_bins);
  HostArray1D<Real> mass_to_num_h("dust_mass_to_num_h", Ndust_bins);
  HostArray1D<Real> num_to_mass_h("dust_num_to_mass_h", Ndust_bins);
  auto mass_fractions_h = Kokkos::create_mirror_view(mass_fractions);

  const Real ratio = std::pow(dust_b / dust_a, 1.0 / Ndust_bins);
  dust_inv_log_step = 1.0 / std::log(ratio);
  Real a_i = dust_a;
  for (int i = 0; i < Ndust_bins + 1; ++i) {
    edges_h(i) = a_i;
    edges_padded_h(i + 1) = a_i;
    a_i *= ratio;
  }
  edges_padded_h(0) = edges_padded_h(1) / ratio / 10.0;
  edges_padded_h(Ndust_bins + 2) = edges_padded_h(Ndust_bins + 1) * ratio * 10.0;
  const Real k_dust_um = dust_k_um(rho_gr);

  for (int i = 0; i < Ndust_bins; ++i) {
    centers_h(i) = 0.5 * (edges_h(i) + edges_h(i + 1));
    centers_cubed_h(i) = cube(centers_h(i));
    const Real e0 = edges_h(i);
    const Real e1 = edges_h(i + 1);
    const Real e0_2 = square(e0);
    const Real e1_2 = square(e1);
    const Real e0_3 = e0_2 * e0;
    const Real e1_3 = e1_2 * e1;
    const Real e0_4 = square(e0_2);
    const Real e1_4 = square(e1_2);
    delta1_h(i) = e1 - e0;
    delta2_h(i) = e1_2 - e0_2;
    delta3_h(i) = e1_3 - e0_3;
    delta4_h(i) = e1_4 - e0_4;
    mass_to_num_h(i) = 4.0 / (k_dust_um * delta4_h(i));
    num_to_mass_h(i) = 0.25 * k_dust_um * delta4_h(i);
  }
  for (int i = 0; i < Ndust_bins + 3; ++i) {
    const Real edge_sq = square(edges_padded_h(i));
    edge4_h(i) = square(edge_sq);
  }

  const Real mex = dust_exp + 4.0;
  Real norm = 0.0;
  for (int i = 0; i < Ndust_bins; ++i) {
    Real frac = 0.0;
    if (std::abs(mex) > 1.0e-12) {
      frac = (std::pow(edges_h(i + 1), mex) - std::pow(edges_h(i), mex)) /
             (std::pow(dust_b, mex) - std::pow(dust_a, mex));
    } else {
      frac = std::log(edges_h(i + 1) / edges_h(i)) / std::log(dust_b / dust_a);
    }
    mass_fractions_h(i) = frac;
    norm += frac;
  }
  for (int i = 0; i < Ndust_bins; ++i) {
    mass_fractions_h(i) /= norm;
  }

  Kokkos::deep_copy(centers, centers_h);
  Kokkos::deep_copy(edges, edges_h);
  Kokkos::deep_copy(edges_padded, edges_padded_h);
  Kokkos::deep_copy(centers_cubed, centers_cubed_h);
  Kokkos::deep_copy(delta1, delta1_h);
  Kokkos::deep_copy(delta2, delta2_h);
  Kokkos::deep_copy(delta3, delta3_h);
  Kokkos::deep_copy(delta4, delta4_h);
  Kokkos::deep_copy(mass_to_num_fac, mass_to_num_h);
  Kokkos::deep_copy(num_to_mass_fac, num_to_mass_h);
  Kokkos::deep_copy(mass_fractions, mass_fractions_h);

  dust_bin_centers = centers;
  dust_bin_edges = edges;
  dust_bin_edges_padded = edges_padded;
  dust_bin_centers_cubed = centers_cubed;
  dust_delta1 = delta1;
  dust_delta2 = delta2;
  dust_delta3 = delta3;
  dust_delta4 = delta4;
  dust_mass_to_num = mass_to_num_fac;
  dust_num_to_mass = num_to_mass_fac;
  dust_bin_mass_fractions = mass_fractions;
  dust_bin_centers_h = centers_h;
  dust_bin_edges_h = edges_h;
  dust_bin_edges_padded_h = edges_padded_h;
  dust_bin_centers_cubed_h = centers_cubed_h;
  dust_edge4_h = edge4_h;
  dust_delta1_h = delta1_h;
  dust_delta2_h = delta2_h;
  dust_delta3_h = delta3_h;
  dust_delta4_h = delta4_h;
  dust_mass_to_num_h = mass_to_num_h;
  dust_num_to_mass_h = num_to_mass_h;
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->phydro == nullptr) {
    throw std::runtime_error("cgm_cooling_flow_dust requires hydro enabled");
  }

  user_srcs_func = UserSource;
  user_bcs_func = UserBoundary;
  user_hist_func = CGMHistory;
  pgen_final_func = FreeProfile;
  user_ref_func = RefinementCondition;

  const int nhydro = pmbp->phydro->nhydro;
  const int nscalars = pmbp->phydro->nscalars;
  if (nscalars < 1) {
    throw std::runtime_error(
        "cgm_cooling_flow_dust requires at least one passive scalar for metallicity");
  }

  if (global_variable::my_rank==0) {
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Units                                         " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Unit Length         : " << 1.0/pmbp->punit->kpc() 
                                          << " kpc"    << std::endl;
    std::cout << "Unit Temperature    : " << 1.0/pmbp->punit->kelvin()
                                          << " K"     << std::endl;     
    std::cout << "Unit Number Density : " << std::pow(pmbp->punit->cm(),3) 
              << " cm^-3" << std::endl;
    std::cout << "Unit Velocity       : " << 1.0/pmbp->punit->km_s() 
                                          << " km/s"  << std::endl;
    std::cout << "Unit Time           : " << 1.0/pmbp->punit->myr() 
                                          << " Myr"   << std::endl;
    std::cout << std::endl;
  }

  // Read in constants
  r_scale   = pin->GetReal("potential", "r_scale");
  rho_scale = pin->GetReal("potential", "rho_scale");
  m_gal     = pin->GetReal("potential", "mass_gal");
  a_gal     = pin->GetReal("potential", "scale_gal");
  z_gal     = pin->GetReal("potential", "z_gal");
  r_200     = pin->GetReal("potential", "r_200");
  rho_mean  = pin->GetReal("potential", "rho_mean");
  r_circ    = pin->GetReal("problem", "r_circ");
  v_circ    = pin->GetReal("problem", "v_circ");
  Z         = pin->GetOrAddReal("problem", "metallicity", 1./3.);

  // Read in SN injection radius and compute energy and mass injection densities
  r_inj = pin->GetReal("SN", "r_inj");
  const Real sphere_vol = (4.0 / 3.0) * M_PI * std::pow(r_inj, 3);
  const Real E_def = 1.0e51;
  const Real M_def = 8.4;
  e_sn = pin->GetOrAddReal("SN", "E_sn", E_def) * pmbp->punit->erg() / sphere_vol;
  m_ej = pin->GetOrAddReal("SN", "M_ej", M_def) * pmbp->punit->msun() / sphere_vol;
  sn_ejecta_Z = pin->GetOrAddReal("problem", "sn_ejecta_Z",
                                  pin->GetOrAddReal("SN", "Z_ej", 0.1));
  sn_ejecta_dust_ratio =
      pin->GetOrAddReal("problem", "sn_ejecta_dust_ratio",
                        pin->GetOrAddReal("SN", "dust_to_metal", 0.0));

  // Read in dust model parameters
  dust_model = pin->GetOrAddBoolean("problem", "dust_model", false);
  if (dust_model) {
    if (nscalars < 2) {
      throw std::runtime_error(
          "dust mode requires at least two passive scalars: metallicity + dust");
    }
    Ndust_bins = nscalars - 1;
    if (Ndust_bins > kMaxDustBins) {
      throw std::runtime_error("Ndust_bins exceeds kMaxDustBins");
    }

    D_Z_init = pin->GetOrAddReal("problem", "D_Z_init", 0.0);
    dust_a = pin->GetReal("problem", "dust_a");
    dust_b = pin->GetReal("problem", "dust_b");
    dust_exp = pin->GetOrAddReal("problem", "dust_exp", -3.5);
    rho_gr = pin->GetReal("problem", "rho_gr");
    grain_porosity = pin->GetOrAddReal("problem", "grain_porosity", 1.0);
    vturb_est_ncell = pin->GetOrAddInteger("problem", "vturb_est_ncell", 1);
    coag_shatt_flag = pin->GetOrAddBoolean("problem", "coag_shatt_flag", false);

    InitializeDustBins();
    grain_porosity_factor = pow(grain_porosity, -2.0 / 3.0);
    if (vturb_est_ncell > 0) {
      const Real stencil_width = static_cast<Real>(2 * vturb_est_ncell + 1);
      const Real n_count = stencil_width * stencil_width * stencil_width;
      vturb_scale = pow(1.0 / n_count, 1.0 / 9.0);
    } else {
      vturb_scale = 0.0;
    }
    if (coag_shatt_flag && vturb_est_ncell > 0) {
      Kokkos::realloc(mach_cache, pmbp->nmb_thispack, indcs.nx3, indcs.nx2, indcs.nx1);
      mach_cache_cycle = -1;
    } else {
      mach_cache = DvceArray4D<Real>();
      mach_cache_cycle = -1;
    }
  } else {
    Ndust_bins = 0;
    D_Z_init = 0.0;
    grain_porosity_factor = 1.0;
    vturb_scale = 0.0;
    dust_bin_centers = DvceArray1D<Real>();
    dust_bin_edges = DvceArray1D<Real>();
    dust_bin_edges_padded = DvceArray1D<Real>();
    dust_bin_centers_cubed = DvceArray1D<Real>();
    dust_delta1 = DvceArray1D<Real>();
    dust_delta2 = DvceArray1D<Real>();
    dust_delta3 = DvceArray1D<Real>();
    dust_delta4 = DvceArray1D<Real>();
    dust_mass_to_num = DvceArray1D<Real>();
    dust_num_to_mass = DvceArray1D<Real>();
    dust_bin_mass_fractions = DvceArray1D<Real>();
    dust_bin_centers_h = HostArray1D<Real>();
    dust_bin_edges_h = HostArray1D<Real>();
    dust_bin_edges_padded_h = HostArray1D<Real>();
    dust_bin_centers_cubed_h = HostArray1D<Real>();
    dust_edge4_h = HostArray1D<Real>();
    dust_delta1_h = HostArray1D<Real>();
    dust_delta2_h = HostArray1D<Real>();
    dust_delta3_h = HostArray1D<Real>();
    dust_delta4_h = HostArray1D<Real>();
    dust_mass_to_num_h = HostArray1D<Real>();
    dust_num_to_mass_h = HostArray1D<Real>();
    mach_cache = DvceArray4D<Real>();
    mach_cache_cycle = -1;
  }

  // Read the density gradient threshold for refinement
  ddens_threshold = pin->GetReal("problem", "ddens_max");

  // Output parameter information
  if (global_variable::my_rank == 0) {
    std::cout << "==============================================" << std::endl;
    std::cout << "Potential Parameters                          " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "r_scale             : " << r_scale    << std::endl;
    std::cout << "rho_scale           : " << rho_scale  << std::endl;
    std::cout << "m_gal               : " << m_gal      << std::endl;
    std::cout << "a_gal               : " << a_gal      << std::endl;
    std::cout << "z_gal               : " << z_gal      << std::endl;
    std::cout << "r_200               : " << r_200      << std::endl;
    std::cout << "rho_mean            : " << rho_mean   << std::endl;
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Other Parameters                              " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "r_circ              : " << r_circ     << std::endl;
    std::cout << "v_circ              : " << v_circ     << std::endl;
    std::cout << "metallicity         : " << Z          << std::endl;
    std::cout << "ddens_threshold     : " << ddens_threshold << std::endl;
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Supernova Parameters                          " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "r_inj               : " << r_inj      << std::endl;
    std::cout << "e_sn                : " << e_sn       << std::endl;
    std::cout << "m_ej                : " << m_ej       << std::endl;
    std::cout << "sn_ejecta_Z         : " << sn_ejecta_Z << std::endl;
    std::cout << "sn_ejecta_dust_ratio: " << sn_ejecta_dust_ratio << std::endl;
    std::cout << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Dust Parameters                               " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "dust_model          : " << dust_model << std::endl;
    if (dust_model) {
      std::cout << "Ndust_bins          : " << Ndust_bins << std::endl;
      std::cout << "D_Z_init            : " << D_Z_init << std::endl;
      std::cout << "dust_a              : " << dust_a << std::endl;
      std::cout << "dust_b              : " << dust_b << std::endl;
      std::cout << "dust_exp            : " << dust_exp << std::endl;
      std::cout << "rho_gr              : " << rho_gr << std::endl;
      std::cout << "grain_porosity      : " << grain_porosity << std::endl;
      std::cout << "vturb_est_ncell     : " << vturb_est_ncell << std::endl;
      std::cout << "coag_shatt_flag     : " << coag_shatt_flag << std::endl;
    }
    std::cout << std::endl;
  }

  // Read the CGM cooling flow profile file
  std::string profile_file = pin->GetString("problem", "profile_file");
  ProfileReaderHost profile_reader_host;
  profile_reader_host.ReadProfiles(profile_file);
  profile_reader = profile_reader_host.CreateDeviceReader();
  if (global_variable::my_rank==0) {
    std::cout << "Successfully loaded CGM profiles from "
              << profile_file << std::endl;
  }

  // Read in the disk profile file
  std::string disk_profile_file = pin->GetString("problem", "disk_profile_file");
  ProfileReaderHost disk_profile_reader_host;
  disk_profile_reader_host.ReadProfiles(disk_profile_file);
  disk_profile_reader = disk_profile_reader_host.CreateDeviceReader();
  if (global_variable::my_rank==0) {
    std::cout << "Successfully loaded disk profiles from "
              << disk_profile_file << std::endl;
  }

  // Count total particles and initialize SN centers buffer
  pmy_mesh_->CountParticles();
  sn_centers_buffer = DvceArray2D<Real>("sn_centers_buffer", 3, pmy_mesh_->nprtcl_total);
  sn_counter = Kokkos::View<int>("sn_counter");
  if (global_variable::my_rank==0) {
    std::cout << "Successfully initialized " << pmy_mesh_->nprtcl_total
              << " particles!" << std::endl;
  }

  if (restart) return;

  // Generate the initial turbulent field
  int nlow   = pin->GetOrAddInteger("problem", "cgm_turb_nlow", 1);
  int nhigh  = pin->GetOrAddInteger("problem", "cgm_turb_nhigh", 8);
  Real expo  = pin->GetOrAddReal("problem", "cgm_turb_expo", 5.0/3.0);
  Real v_rms = pin->GetOrAddReal("problem", "cgm_turb_rms", 0.1);
  Real cgm_turb_xscale = pin->GetOrAddReal("problem", "cgm_turb_xscale", 0.01);
  Real cgm_turb_yscale = pin->GetOrAddReal("problem", "cgm_turb_yscale", 0.01);
  Real cgm_turb_zscale = pin->GetOrAddReal("problem", "cgm_turb_zscale", 0.01);

  // Initialize random state
  RNG_State rstate;
  rstate.idum = -1;

  // Domain size
  Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  Real ly = pmy_mesh_->mesh_size.x2max - pmy_mesh_->mesh_size.x2min;
  Real lz = pmy_mesh_->mesh_size.x3max - pmy_mesh_->mesh_size.x3min;
  Real dkx = 2.0*M_PI/lx;
  Real dky = 2.0*M_PI/ly;
  Real dkz = 2.0*M_PI/lz;

  // Count modes
  int nmodes = 0;
  for (int nkx = -nhigh; nkx <= nhigh; nkx++) {
    for (int nky = -nhigh; nky <= nhigh; nky++) {
      for (int nkz = -nhigh; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) {
          nmodes++;
        }
  }}}

  // Allocate arrays
  DualArray2D<Real> k_modes, aka, akb;
  Kokkos::realloc(k_modes, 3, nmodes);
  Kokkos::realloc(aka, 3, nmodes);
  Kokkos::realloc(akb, 3, nmodes);

  // Generate modes
  int nmode = 0;
  Real total_energy = 0.0;
  for (int nkx = -nhigh; nkx <= nhigh; nkx++) {
    for (int nky = -nhigh; nky <= nhigh; nky++) {
      for (int nkz = -nhigh; nkz <= nhigh; nkz++) {
        if (nkx == 0 && nky == 0 && nkz == 0) continue;
        int nsqr = nkx*nkx + nky*nky + nkz*nkz;
        if (nsqr >= nlow*nlow && nsqr <= nhigh*nhigh) {
          Real kx = dkx*nkx, ky = dky*nky, kz = dkz*nkz;
          Real kiso = sqrt(kx*kx + ky*ky + kz*kz);

          k_modes.h_view(0, nmode) = kx;
          k_modes.h_view(1, nmode) = ky;
          k_modes.h_view(2, nmode) = kz;

          Real norm = 1.0/pow(kiso, (expo+2.0)/2.0);

          Real aval[3], bval[3];
          for (int dir = 0; dir < 3; dir++) {
            aval[dir] = norm * RanGaussianSt(&rstate);
            bval[dir] = norm * RanGaussianSt(&rstate);
          }

          Real k_dirs[3] = {kx, ky, kz};
          Real ka = kx*aval[0] + ky*aval[1] + kz*aval[2];
          Real kb = kx*bval[0] + ky*bval[1] + kz*bval[2];

          for (int dir = 0; dir < 3; dir++) {
            aval[dir] -= k_dirs[dir]*ka/(kiso*kiso);
            bval[dir] -= k_dirs[dir]*kb/(kiso*kiso);

            aka.h_view(dir,nmode) = aval[dir];
            akb.h_view(dir,nmode) = bval[dir];

            total_energy += 0.5*(aval[dir]*aval[dir] + bval[dir]*bval[dir]);
          }
          nmode++;
        }
      }
    }
  }

  Real v_norm = v_rms/sqrt(total_energy);

  k_modes.template modify<HostMemSpace>();
  k_modes.template sync<DevExeSpace>();
  aka.template modify<HostMemSpace>();
  aka.template sync<DevExeSpace>();
  akb.template modify<HostMemSpace>();
  akb.template sync<DevExeSpace>();

  // Capture variables for kernel
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  EOS_Data &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  auto &profile = profile_reader;
  auto &disk_profile = disk_profile_reader;
  auto dust_fracs = dust_bin_mass_fractions;
  const bool dust_enabled = dust_model;
  const int nmetal = nhydro;
  const int ndusti = nhydro + 1;
  const int ndust = Ndust_bins;
  const Real D_Z0 = D_Z_init;

  Real G = pmbp->punit->grav_constant();
  Real r_s = r_scale;
  Real rho_s = rho_scale;
  Real m_g = m_gal;
  Real a_g = a_gal;
  Real z_g = z_gal;
  Real r_m = r_200;
  Real rho_m = rho_mean;

  Real r_c = r_circ;
  Real v_c = v_circ;
  Real Zsol = Z_solar;
  Real Z_ = Z;

  // Use loaded profiles
  par_for("pgen_ic", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x1l =   LeftEdgeX(i-is, nx1, x1min, x1max);
    Real x1r = LeftEdgeX(i+1-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x2l =   LeftEdgeX(j-js, nx2, x2min, x2max);
    Real x2r = LeftEdgeX(j+1-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
    SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
    SetEquilibriumState(u0, m, k, j, i, x1v, x2v, x3v,
                        x1l, x1r, x2l, x2r, G, r_s, rho_s,
                        m_g, a_g, z_g, r_m, rho_m, gm1, disk_profile);

    Real rho = u0(m, IDN, k, j, i);
    u0(m, nmetal, k, j, i) = Z_ * Zsol * rho;
    for (int n = ndusti; n < nhydro + nscalars; ++n) {
      u0(m, n, k, j, i) = 0.0;
    }

    if (dust_enabled) {
      Real dust_mass = D_Z0 * Z_ * Zsol * rho;
      for (int id = 0; id < ndust; ++id) {
        u0(m, ndusti + id, k, j, i) = dust_mass * dust_fracs(id);
      }
    }

    // Compute turbulent velocities by summing Fourier modes
    Real vx = 0.0, vy = 0.0, vz = 0.0;
    for (int n = 0; n < nmodes; n++) {
      Real phase = k_modes.d_view(0,n)*x1v 
	         + k_modes.d_view(1,n)*x2v 
	         + k_modes.d_view(2,n)*x3v;
      Real cos_phase = cos(phase);
      Real sin_phase = sin(phase);

      vx += aka.d_view(0,n)*cos_phase - akb.d_view(0,n)*sin_phase;
      vy += aka.d_view(1,n)*cos_phase - akb.d_view(1,n)*sin_phase;
      vz += aka.d_view(2,n)*cos_phase - akb.d_view(2,n)*sin_phase;
    }

    // Attenuate in the center by 1 - Gaussian
    Real att = 1.0 - exp(-0.5 * ( SQR(x1v)/SQR(cgm_turb_xscale)
                                + SQR(x2v)/SQR(cgm_turb_yscale)
                                + SQR(x3v)/SQR(cgm_turb_zscale)));

    // Normalize to desired RMS velocity
    vx *= v_norm*att; vy *= v_norm*att; vz *= v_norm*att;

    // Add to conserved variables
    Real rho_local = u0(m,IDN,k,j,i);
    Real rho_v1 = u0(m,IM1,k,j,i);
    Real rho_v2 = u0(m,IM2,k,j,i);
    Real rho_v3 = u0(m,IM3,k,j,i);

    u0(m,IEN,k,j,i) += 0.5 * rho_local * (SQR(vx) + SQR(vy) + SQR(vz));
    u0(m,IEN,k,j,i) += rho_v1 * vx + rho_v2 * vy + rho_v3 * vz;
    u0(m,IM1,k,j,i) += rho_local * vx;
    u0(m,IM2,k,j,i) += rho_local * vy;
    u0(m,IM3,k,j,i) += rho_local * vz;
  });

  if (global_variable::my_rank==0) {
    std::cout << "Successfully initialized grid!" << std::endl;
  }

  return;
}

KOKKOS_INLINE_FUNCTION
void SetCoolingFlowState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real gm1, const ProfileReader &profile) {
    // Calculate radius
    Real r = sqrt(x1v*x1v + x2v*x2v + x3v*x3v);
      
    // Get values from profiles via interpolation
    Real rho  = profile.GetDensity(r);
  Real temp = profile.GetTemperature(r);
    Real vr   = profile.GetVelocity(r);
  Real rmin = profile.GetRmin();

    // set vr to go to zero if r < rmin
  if (r < rmin) {
    vr *= (r / rmin);
  }

    // Calculate pressure from temperature
  Real press = rho * temp;
    
    // Set radial velocity components based on position
  Real v1 = 0.0, v2 = 0.0, v3 = 0.0;
  constexpr Real tiny = 1.0e-20;
    if (r > tiny) {  // Avoid division by zero
      // Negative sign accounts for inflowing vr
    v1 = -vr * x1v / r;
    v2 = -vr * x2v / r;
    v3 = -vr * x3v / r;
  }

    // Set state variables
  u0(m, IDN, k, j, i) = rho;
  u0(m, IM1, k, j, i) = rho * v1;
  u0(m, IM2, k, j, i) = rho * v2;
  u0(m, IM3, k, j, i) = rho * v3;
    u0(m, IEN, k, j, i) = press/gm1 + 0.5*rho*(SQR(v1) + SQR(v2) + SQR(v3));
}

KOKKOS_INLINE_FUNCTION
void SetRotation(const DvceArray5D<Real> &u0,
                 int m, int k, int j, int i,
                 Real x1v, Real x2v, Real x3v,
                 Real r_circ, Real v_circ) {
  // Calculate radius
  Real r = sqrt(x1v*x1v + x2v*x2v + x3v*x3v);
  Real R = sqrt(x1v*x1v + x2v*x2v);

  // Calculate azimuthal velocity
  Real vx = 0.0, vy = 0.0, vz = 0.0;
  constexpr Real tiny = 1.0e-20;
  if (r > tiny && R > tiny) {  // Avoid division by zero
    Real v_phi = 0.0;
    Real sin_theta = R / r;

    if (r < r_circ) {
      v_phi = v_circ * sin_theta;
    }
    if (r > r_circ) {
      v_phi = v_circ * sin_theta * r_circ / r;
    }
    
    // Calculate azimuthal velocity components
    vx = -v_phi * x2v / R;
    vy = v_phi * x1v / R;
  }

  // Set state variables
  Real rho = u0(m,IDN,k,j,i);
  Real rho_v1 = u0(m,IM1,k,j,i);
  Real rho_v2 = u0(m,IM2,k,j,i);
  Real rho_v3 = u0(m,IM3,k,j,i);

  u0(m,IEN,k,j,i) += 0.5 * rho * (SQR(vx) + SQR(vy) + SQR(vz));
  u0(m,IEN,k,j,i) += rho_v1 * vx + rho_v2 * vy + rho_v3 * vz;
  u0(m,IM1,k,j,i) += rho * vx;
  u0(m,IM2,k,j,i) += rho * vy;
  u0(m,IM3,k,j,i) += rho * vz;
}

KOKKOS_INLINE_FUNCTION
void SetEquilibriumState(const DvceArray5D<Real> &u0,
                         int m, int k, int j, int i,
                         Real x1v, Real x2v, Real x3v,
                         Real x1l, Real x1r, Real x2l, Real x2r,
                         Real G, Real r_s, Real rho_s, Real m_g,
                         Real a_g, Real z_g, Real r_m, Real rho_m,
                         Real gm1, const ProfileReader &disk_profile) {
    // Calculate radius
    Real R = sqrt(x1v * x1v + x2v * x2v);
    Real R1l = sqrt(x1l * x1l + x2v * x2v);
    Real R1r = sqrt(x1r * x1r + x2v * x2v);
    Real R2l = sqrt(x1v * x1v + x2l * x2l);
    Real R2r = sqrt(x1v * x1v + x2r * x2r);

    // Don't extrapolate past the last entry in the table
  if (R > disk_profile.GetRmax()) return;

    // Calculate Gravitational Potentialsi
    Real c_out = (4.0/3.0) * pow(5 * r_m, 1.5);
    
    Real phi0    = GravPot(x1v, x2v, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi0_1l = GravPot(x1l, x2v, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi0_1r = GravPot(x1r, x2v, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi0_2l = GravPot(x1v, x2l, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi0_2r = GravPot(x1v, x2r, 0.0, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);

  Real phi   = GravPot(x1v, x2v, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi1r = GravPot(x1r, x2v, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi2l = GravPot(x1v, x2l, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi2r = GravPot(x1v, x2r, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);
  Real phi1l = GravPot(x1l, x2v, x3v, G, r_s, rho_s, m_g, a_g, z_g, c_out, rho_m);

  Real f_x1_ = -(phi1r - phi1l) / (x1r - x1l);
  Real f_x2_ = -(phi2r - phi2l) / (x2r - x2l);
    Real dPhi_dR = f_x1_ * x1v / R  + f_x2_ * x2v / R;

    // Compute sound speed squared
  Real temp = disk_profile.GetTemperature(R);
    Real cs2 = temp; // isothermal, no factor of gamma

    // Get densities from profiles via interpolation
    Real rho = disk_profile.GetDensity(R) * exp(-(phi - phi0) / cs2);
    Real rho_1l = disk_profile.GetDensity(R1l) * exp(-(phi1l - phi0_1l) / cs2);
    Real rho_1r = disk_profile.GetDensity(R1r) * exp(-(phi1r - phi0_1r) / cs2);
    Real rho_2l = disk_profile.GetDensity(R2l) * exp(-(phi2l - phi0_2l) / cs2);
    Real rho_2r = disk_profile.GetDensity(R2r) * exp(-(phi2r - phi0_2r) / cs2); 
    Real p_x1_ = temp*(rho_1r - rho_1l) / (x1r - x1l);
    Real p_x2_ = temp*(rho_2r - rho_2l) / (x2r - x2l);

  Real dP_dR_over_rho = 0.0;
  constexpr Real tiny = 1.0e-20;
  if (rho > tiny) {
      dP_dR_over_rho = (p_x1_ * x1v  + p_x2_ * x2v) / (R * rho);
  }

    // Calculate circular velocity
    Real v_phi = sqrt(R*fmax(dP_dR_over_rho - dPhi_dR, 0.0));
    
    // Calculate azimuthal velocity
  Real v1 = 0.0, v2 = 0.0;
    if (R > tiny) {  // Avoid division by zero     
      // Calculate azimuthal velocity components
    v1 = -v_phi * x2v / R;
      v2 =  v_phi * x1v / R;
  }

    // Combine cgm and disk material consistently
    Real rho_cgm  = u0(m, IDN, k, j, i);
  Real mom1_cgm = u0(m, IM1, k, j, i);
  Real mom2_cgm = u0(m, IM2, k, j, i);
  Real mom3_cgm = u0(m, IM3, k, j, i);
    Real E_cgm    = u0(m, IEN, k, j, i);

    Real KE_cgm  = 0.5*(SQR(mom1_cgm) + SQR(mom2_cgm) + SQR(mom3_cgm))/rho_cgm;
  Real Eth_cgm = E_cgm - KE_cgm;

    Real rho_tot  = rho_cgm + rho;
  Real mom1_tot = mom1_cgm + rho * v1;
  Real mom2_tot = mom2_cgm + rho * v2;
  Real mom3_tot = mom3_cgm;

    Real KE_tot = 0.5*(SQR(mom1_tot)+SQR(mom2_tot)+SQR(mom3_tot))/rho_tot;

    // Set state variables
  u0(m, IDN, k, j, i) = rho_tot;
  u0(m, IM1, k, j, i) = mom1_tot;
  u0(m, IM2, k, j, i) = mom2_tot;
  u0(m, IM3, k, j, i) = mom3_tot;
    u0(m, IEN, k, j, i) = Eth_cgm + (rho * temp)/gm1 + KE_tot;
}

//===========================================================================//
//                              Source Terms                                 //
//===========================================================================//

void UserSource(Mesh* pm, const Real bdt) {
  SNSource(pm, bdt);
  GravitySource(pm, bdt);
  if (dust_model) {
    AddDustSource(pm, bdt);
  }
  return;
}

void GravitySource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;

  Real G = pmbp->punit->grav_constant();
  Real r_s = r_scale;
  Real rho_s = rho_scale;
  Real m_g = m_gal;
  Real a_g = a_gal;
  Real z_g = z_gal;
  Real rho_m = rho_mean;
  Real c_out = (4.0/3.0) * pow(5 * r_200, 1.5);

  par_for("gravity_source", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;

    const Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    const Real x1l = LeftEdgeX(i-is,   nx1, x1min, x1max);
    const Real x1r = LeftEdgeX(i+1-is, nx1, x1min, x1max);

    const Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    const Real x2l = LeftEdgeX(j-js,   nx2, x2min, x2max);
    const Real x2r = LeftEdgeX(j+1-js, nx2, x2min, x2max);

    const Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
    const Real x3l = LeftEdgeX(k-ks,   nx3, x3min, x3max);
    const Real x3r = LeftEdgeX(k+1-ks, nx3, x3min, x3max);

    Real phi1l = GravPot(x1l,x2v,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);
    Real phi1r = GravPot(x1r,x2v,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);

    Real phi2l = GravPot(x1v,x2l,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);
    Real phi2r = GravPot(x1v,x2r,x3v,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);

    Real phi3l = GravPot(x1v,x2v,x3l,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);
    Real phi3r = GravPot(x1v,x2v,x3r,G,r_s,rho_s,m_g,a_g,z_g,c_out,rho_m);

    constexpr Real tiny = 1e-20;
    Real f_x1_ = -(phi1r - phi1l) / fmax(x1r - x1l, tiny);
    Real f_x2_ = -(phi2r - phi2l) / fmax(x2r - x2l, tiny);
    Real f_x3_ = -(phi3r - phi3l) / fmax(x3r - x3l, tiny);

    Real density = w0(m, IDN, k, j, i);
    Real src_x1 = bdt * density * f_x1_;
    Real src_x2 = bdt * density * f_x2_;
    Real src_x3 = bdt * density * f_x3_;

    u0(m,IM1,k,j,i) += src_x1;	
    u0(m,IM2,k,j,i) += src_x2;
    u0(m,IM3,k,j,i) += src_x3;
    u0(m,IEN,k,j,i) += (src_x1 * w0(m,IVX,k,j,i) +
                        src_x2 * w0(m,IVY,k,j,i) +
                        src_x3 * w0(m,IVZ,k,j,i));
    
  });

  return;
}

KOKKOS_INLINE_FUNCTION
Real GravPot(Real x1, Real x2, Real x3,
             Real G, Real r_s, Real rho_s,
             Real M_gal, Real a_gal, Real z_gal,
             Real c_outer, Real rho_mean) {
  const Real R2 = fma(x1, x1 , x2*x2);
  const Real R  = sqrt(R2);
  const Real r2 = fma(x3 , x3 , R2);
  const Real r  = sqrt(fmax(r2, 1e-20));

  // NFW component
  Real x = r / r_s;
  Real phi_NFW = -4 * M_PI * G * rho_s * SQR(r_s) * log1p(x) / x;
  
  // Miyamoto-Nagai model
  Real phi_MN = -G * M_gal / sqrt(R2 + SQR(sqrt(fma(x3 , x3 , z_gal*z_gal)) + a_gal));
  
  // Outer component
  Real phi_Outer = 4 * M_PI * G * rho_mean * (c_outer * sqrt(r) + (1.0/6.0) * r2);
  
  // Total potential
  return phi_NFW + phi_MN + phi_Outer;
}

void UpdateTurbulenceCache(Mesh *pm) {
  if (!dust_model || !coag_shatt_flag || vturb_est_ncell <= 0 ||
      mach_cache_cycle == pm->ncycle) {
    return;
  }

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const auto &w0 = pmbp->phydro->w0;

  Real KELVIN = 1.0;
  Real KM_S = 1.0;
  if (pmbp->punit != nullptr) {
    KELVIN = pmbp->punit->temperature_cgs();
    KM_S = pmbp->punit->km_s();
  } else if (global_variable::my_rank == 0) {
    std::cout << "ERROR: <units> block missing..." << std::endl;
    std::exit(1);
  }

  const Real gamma = pmbp->phydro->peos->eos_data.gamma;
  const Real gm1 = gamma - 1.0;
  auto mach_cache_local = mach_cache;
  const int vturb_est_ncell_local = vturb_est_ncell;
  const Real vturb_scale_local = vturb_scale;

  par_for("update_turbulence_cache", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kk = k - ks;
    const int jj = j - js;
    const int ii = i - is;

    const Real dens = w0(m, IDN, k, j, i);
    const Real temp = (w0(m, IEN, k, j, i) * gm1) / dens * KELVIN;

    Real vx_mean = 0.0;
    Real vy_mean = 0.0;
    Real vz_mean = 0.0;
    Real M2x = 0.0;
    Real M2y = 0.0;
    Real M2z = 0.0;
    int n_count = 0;

    for (int ia = -vturb_est_ncell_local; ia <= vturb_est_ncell_local; ++ia) {
      for (int ja = -vturb_est_ncell_local; ja <= vturb_est_ncell_local; ++ja) {
        for (int ka = -vturb_est_ncell_local; ka <= vturb_est_ncell_local; ++ka) {
          const Real vx = w0(m, IVX, k + ka, j + ja, i + ia);
          const Real vy = w0(m, IVY, k + ka, j + ja, i + ia);
          const Real vz = w0(m, IVZ, k + ka, j + ja, i + ia);

          ++n_count;
          const Real deltax = vx - vx_mean;
          const Real deltay = vy - vy_mean;
          const Real deltaz = vz - vz_mean;
          vx_mean += deltax / n_count;
          vy_mean += deltay / n_count;
          vz_mean += deltaz / n_count;
          M2x += deltax * (vx - vx_mean);
          M2y += deltay * (vy - vy_mean);
          M2z += deltaz * (vz - vz_mean);
        }
      }
    }

    Real vturb_g = sqrt((M2x + M2y + M2z) / Real(n_count + 1));
    vturb_g *= vturb_scale_local / KM_S;
    const Real cs_g = sqrt(gamma * temp / KELVIN) * 1.0e-5;
    mach_cache_local(m, kk, jj, ii) = (cs_g > 0.0) ? (vturb_g / cs_g) : 0.0;
  });

  mach_cache_cycle = pm->ncycle;
}

void AddDustSource(Mesh *pm, const Real bdt) {
  if (!dust_model || Ndust_bins <= 0) return;

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int ni = ie - is + 1;
  const int nj = je - js + 1;
  const int nk = ke - ks + 1;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->phydro->u0;
  const auto &w0 = pmbp->phydro->w0;
  const Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  const int nhydro = pmbp->phydro->nhydro;

  Real KELVIN = 1.0;
  Real MYR = 1.0;
  Real SEC = 1.0;
  Real CM3 = 1.0;
  if (pmbp->punit != nullptr) {
    KELVIN = pmbp->punit->temperature_cgs();
    MYR = pmbp->punit->myr();
    SEC = pmbp->punit->s();
    CM3 = cube(pmbp->punit->cm());
  } else if (global_variable::my_rank == 0) {
    std::cout << "ERROR: <units> block missing..." << std::endl;
    std::exit(1);
  }

  const int nmetal = nhydro;
  const int ndusti = nhydro + 1;
  const int ndust = Ndust_bins;
  const int n_edges = ndust + 3;
  const Real dust_a_local = dust_a;
  const Real dust_b_local = dust_b;
  const Real rho_gr_local = rho_gr;
  const Real grain_porosity_factor_local = grain_porosity_factor;
  const Real Z_solar_local = Z_solar;
  const bool coag_shatt = coag_shatt_flag;
  const bool use_mach_cache = coag_shatt && (vturb_est_ncell > 0);
  const Real bdt_sec = bdt / SEC;
  const auto edges_pad = dust_bin_edges_padded;
  const auto edges = dust_bin_edges;
  const auto dbins = dust_bin_centers;
  const auto dbins3 = dust_bin_centers_cubed;
  const auto delta1 = dust_delta1;
  const auto delta2 = dust_delta2;
  const auto delta3 = dust_delta3;
  const auto mass_to_num_fac = dust_mass_to_num;
  const auto num_to_mass_fac = dust_num_to_mass;
  const Real dust_inv_log_step_local = dust_inv_log_step;

  if (use_mach_cache) UpdateTurbulenceCache(pm);
  auto mach_cache_local = mach_cache;

  constexpr Real Z_tol = 1.0e-20;
  constexpr Real D_tol = 1.0e-20;
  constexpr Real M_g_tol = 1.0e-6;

  DvceArray4D<Real> dust_tot_before("dust_tot_before", nmb1 + 1, nk, nj, ni);

  par_for("dust_source_shift", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dens = w0(m, IDN, k, j, i);
    const Real temp = (w0(m, IEN, k, j, i) * gm1) / dens * KELVIN;
    const Real Z_local = Kokkos::max(Z_tol, w0(m, nmetal, k, j, i));
    const int kk = k - ks;
    const int jj = j - js;
    const int ii = i - is;

    Real dust_arr[kMaxDustBinSlots];
    dust_tot_before(m, kk, jj, ii) = apply_sputter_accrete_rebin(
        dust_arr, u0, m, k, j, i, ndusti, ndust, mass_to_num_fac, edges_pad, n_edges,
        dens, temp, Z_local, Z_solar_local, grain_porosity_factor_local, rho_gr_local,
        dust_a_local, dust_b_local, dust_inv_log_step_local, bdt, MYR);
    if (dust_tot_before(m, kk, jj, ii) <= D_tol) return;

    for (int id = 0; id < ndust; ++id) {
      u0(m, ndusti + id, k, j, i) = dust_arr[id + 1] * num_to_mass_fac(id);
    }
  });

  if (coag_shatt) {
    par_for("dust_source_coagshatt", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real dens = w0(m, IDN, k, j, i);
      const Real temp = (w0(m, IEN, k, j, i) * gm1) / dens * KELVIN;
      const Real vol_cc =
          size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3 / CM3;

      const int kk = k - ks;
      const int jj = j - js;
      const int ii = i - is;

      const Real dust_mass_before = dust_tot_before(m, kk, jj, ii);
      if (dust_mass_before <= D_tol) return;

      Real dust_arr[kMaxDustBinSlots];
      dust_arr[0] = 0.0;
      dust_arr[ndust + 1] = 0.0;

      for (int id = 0; id < ndust; ++id) {
        dust_arr[id + 1] = u0(m, ndusti + id, k, j, i) * mass_to_num_fac(id);
      }

      const Real M_g = use_mach_cache ? mach_cache_local(m, kk, jj, ii) : 0.0;
      if (M_g < M_g_tol) return;

      Real mass_in_bin[kMaxDustBins];
      for (int id = 0; id < ndust; ++id) mass_in_bin[id] = 0.0;

      add_coagulate_streamed(bdt_sec, mass_in_bin, dust_arr + 1, dbins, dbins3, edges,
                             delta1, delta2, delta3, ndust, M_g, rho_gr_local,
                             dens, temp, vol_cc, dust_inv_log_step_local);
      add_shatters_streamed(bdt_sec, mass_in_bin, dust_arr + 1, dbins, dbins3, edges,
                            delta1, delta2, delta3, ndust, M_g, rho_gr_local,
                            dens, temp, vol_cc, dust_inv_log_step_local);

      for (int id = 0; id < ndust; ++id) {
        u0(m, ndusti + id, k, j, i) =
            (dust_arr[id + 1] + mass_in_bin[id] * mass_to_num_fac(id)) *
            num_to_mass_fac(id);
      }
    });
  }

  par_for("dust_source_finalize", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real dens = w0(m, IDN, k, j, i);

    const int kk = k - ks;
    const int jj = j - js;
    const int ii = i - is;

    const Real dust_mass_before = dust_tot_before(m, kk, jj, ii);
    if (dust_mass_before <= D_tol) return;

    Real dust_arr[kMaxDustBinSlots];
    dust_arr[0] = 0.0;
    dust_arr[ndust + 1] = 0.0;
    for (int id = 0; id < ndust; ++id) {
      dust_arr[id + 1] = u0(m, ndusti + id, k, j, i) * mass_to_num_fac(id);
    }

    const Real dust_tot = dist_num_integral(
        dust_arr, edges_pad, n_edges, dust_a_local, dust_b_local, rho_gr_local,
        dust_inv_log_step_local, true);
    const Real rate_Z = dust_mass_before - dust_tot;

    finalize_dust_source(u0, dust_arr, m, k, j, i, nmetal, ndusti, ndust,
                         dens, dust_tot, rate_Z, num_to_mass_fac);
  });
}

void SNSource(Mesh* pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  int nhydro = pmbp->phydro->nhydro;
  int nscalars = pmbp->phydro->nscalars;
  auto &u0 = pmbp->phydro->u0;

  Real dr = r_inj;

  if (pm->ncycle != last_sn_detect_cycle) {
    last_sn_detect_cycle = pm->ncycle;

    auto &pr = pmbp->ppart->prtcl_rdata;
    auto &pi = pmbp->ppart->prtcl_idata;
    int npart = pmbp->ppart->nprtcl_thispack;

    auto gids = pmbp->gids;

    Real time = pm->time;
    int nrdata = pmbp->ppart->nrdata;
    Real unit_time = pmbp->punit->time_cgs();

    auto &sn_centers = sn_centers_buffer;
    Kokkos::deep_copy(sn_counter, 0);
    auto d_counter = sn_counter;

    par_for("sn_source", DevExeSpace(), 0, npart-1, KOKKOS_LAMBDA(const int p) {
    
      Real next_sn_time = pr(nrdata-1, p);
    
      if (time > next_sn_time) {
        // Update particle sn tracking
        pi(2, p) += 1;
        int sn_idx = pi(2, p);
        Real par_t_create = pr(nrdata-3, p);
        Real cluster_mass = pr(nrdata-2, p);
        pr(nrdata-1, p) = GetNthSNTime(cluster_mass, par_t_create, unit_time, sn_idx);

        // Register SN center location and adjust for boundaries
        // Warning : if r_inj is large this will lead to unexpected behavior at the start
        int idx = Kokkos::atomic_fetch_add(&d_counter(), 1);
        int m = pi(PGID, p) - gids;

        Real x1min = size.d_view(m).x1min;
        Real x1max = size.d_view(m).x1max;
        Real x2min = size.d_view(m).x2min;
        Real x2max = size.d_view(m).x2max;
        Real x3min = size.d_view(m).x3min;
        Real x3max = size.d_view(m).x3max;

        sn_centers(0, idx) =
            Kokkos::min(Kokkos::max(pr(IPX, p), x1min + dr), x1max - dr);
        sn_centers(1, idx) =
            Kokkos::min(Kokkos::max(pr(IPY, p), x2min + dr), x2max - dr);
        sn_centers(2, idx) =
            Kokkos::min(Kokkos::max(pr(IPZ, p), x3min + dr), x3max - dr);
      }
    });

    DevExeSpace().fence();

    // Get number of SNe on host
    Kokkos::deep_copy(num_sn_this_cycle, d_counter);
  }

  if (num_sn_this_cycle > 0) {
    Real beta = bdt / pm->dt;
    Real e_sn_ = e_sn * beta;
    Real m_ej_ = m_ej * beta;
    Real Z_ej_ = sn_ejecta_Z;
    Real dust_ratio_ej = sn_ejecta_dust_ratio;
    Real dust_ej_ = Z_ej_ * m_ej_ * dust_ratio_ej;
    int nmetal = nhydro;
    int ndusti = nhydro + 1;
    bool dust_enabled = dust_model;
    int ndust = Ndust_bins;
    auto dust_fracs = dust_bin_mass_fractions;
    auto &sn_centers = sn_centers_buffer;
    int num_sn = num_sn_this_cycle;

    par_for("sn_injection", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real x1min = size.d_view(m).x1min;
      Real x1max = size.d_view(m).x1max;
      Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

      Real x2min = size.d_view(m).x2min;
      Real x2max = size.d_view(m).x2max;
      Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

      Real x3min = size.d_view(m).x3min;
      Real x3max = size.d_view(m).x3max;
      Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

      for (int sn = 0; sn < num_sn; ++sn) {
        Real sn_x = sn_centers(0, sn);
        Real sn_y = sn_centers(1, sn);
        Real sn_z = sn_centers(2, sn);
              
        // Calculate distance from SN center
        Real dx = x1v - sn_x;
        Real dy = x2v - sn_y;
        Real dz = x3v - sn_z;
        Real r = sqrt(dx*dx + dy*dy + dz*dz);

        // Inject energy if within injection radius
        if (r <= dr) {
          u0(m, IDN, k, j, i) += m_ej_;
          u0(m, IEN, k, j, i) += e_sn_;
          if (nscalars > 0) {
            u0(m, nmetal, k, j, i) += Z_ej_ * m_ej_ * (1.0 - dust_ratio_ej);
          }
          if (dust_enabled) {
            for (int id = 0; id < ndust; ++id) {
              u0(m, ndusti + id, k, j, i) += dust_ej_ * dust_fracs(id);
            }
          }
        }
      }
    });
  }
}


//===========================================================================//
//                             User Boundary                                 //
//===========================================================================//

void UserBoundary(Mesh* pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;

  int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  int nhydro = pmbp->phydro->nhydro;
  int nscalars = pmbp->phydro->nscalars;
  auto &u0 = pmbp->phydro->u0;
  auto &eos = pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  auto &profile = profile_reader;

  Real r_c = r_circ;
  Real v_c = v_circ;
  Real Zsol = Z_solar;
  Real Z_ = Z;
  int nmetal = nhydro;
  int ndusti = nhydro + 1;
  int ndust = Ndust_bins;
  bool dust_enabled = dust_model;
  auto dust_fracs = dust_bin_mass_fractions;
  Real D_Z0 = D_Z_init;

  // Handle X1 boundaries
  par_for("static_x1", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (n2-1), 0, (ng-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Inner X1 boundary
    if (mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
      Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
      
      SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
      SetScalarState(u0, m, k, j, i, nhydro, nscalars, Z_, Zsol,
                     dust_enabled, ndust, D_Z0, dust_fracs);
    }

    // Outer X1 boundary
    if (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
      int i_out = ie + i + 1;
      Real x1v = CellCenterX(i_out-is, indcs.nx1, x1min, x1max);
      
      SetCoolingFlowState(u0, m, k, j, i_out, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i_out, x1v, x2v, x3v, r_c, v_c);
      SetScalarState(u0, m, k, j, i_out, nhydro, nscalars, Z_, Zsol,
                     dust_enabled, ndust, D_Z0, dust_fracs);
    }
  });

  // Handle X2 boundaries
  par_for("static_x2", DevExeSpace(), 0, nmb1, 0, (n3-1), 0, (ng-1), 0, (n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    // Inner X2 boundary
    if (mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user) {
      Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
      
      SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
      SetScalarState(u0, m, k, j, i, nhydro, nscalars, Z_, Zsol,
                     dust_enabled, ndust, D_Z0, dust_fracs);
    }

    // Outer X2 boundary
    if (mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user) {
      int j_out = je + j + 1;
      Real x2v = CellCenterX(j_out-js, indcs.nx2, x2min, x2max);
      
      SetCoolingFlowState(u0, m, k, j_out, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j_out, i, x1v, x2v, x3v, r_c, v_c);
      SetScalarState(u0, m, k, j_out, i, nhydro, nscalars, Z_, Zsol,
                     dust_enabled, ndust, D_Z0, dust_fracs);
    }
  });

  // Handle X3 boundaries
  par_for("static_x3", DevExeSpace(), 0, nmb1, 0, (ng-1), 0, (n2-1), 0, (n1-1),
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;

    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

    // Inner X3 boundary
    if (mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user) {
      Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
      
      SetCoolingFlowState(u0, m, k, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k, j, i, x1v, x2v, x3v, r_c, v_c);
      SetScalarState(u0, m, k, j, i, nhydro, nscalars, Z_, Zsol,
                     dust_enabled, ndust, D_Z0, dust_fracs);
    }

    // Outer X3 boundary
    if (mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user) {
      int k_out = ke + k + 1;
      Real x3v = CellCenterX(k_out-ks, indcs.nx3, x3min, x3max);
      
      SetCoolingFlowState(u0, m, k_out, j, i, x1v, x2v, x3v, gm1, profile);
      SetRotation(u0, m, k_out, j, i, x1v, x2v, x3v, r_c, v_c);
      SetScalarState(u0, m, k_out, j, i, nhydro, nscalars, Z_, Zsol,
                     dust_enabled, ndust, D_Z0, dust_fracs);
    }
  });
  
}

//===========================================================================//
//                              Refinement                                   //
//===========================================================================//

// Refine region based on density gradient threshold
void RefinementCondition(MeshBlockPack* pmbp) {
  Mesh *pmesh       = pmbp->pmesh;
  int nmb           = pmbp->nmb_thispack;
  int mbs           = pmesh->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &indcs       = pmesh->mb_indcs;
  int &is = indcs.is, nx1 = indcs.nx1;
  int &js = indcs.js, nx2 = indcs.nx2;
  int &ks = indcs.ks, nx3 = indcs.nx3;
  const int nkji = nx3 * nx2 * nx1;
  const int nji  = nx2 * nx1;
  auto &u0       = pmbp->phydro->u0;
  auto &w0       = pmbp->phydro->w0;
  const Real ddens_threshold_local = ddens_threshold;

  par_for_outer("UserRefineCond",DevExeSpace(), 0, 0, 0, (nmb-1),
  KOKKOS_LAMBDA(TeamMember_t tmember, const int m) {

    Real team_ddmax;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tmember, nkji),
    [=](const int idx, Real& ddmax) {
      int k = (idx)/nji;
      int j = (idx - k*nji)/nx1;
      int i = (idx - k*nji - j*nx1) + is;
      j += js;
      k += ks;

      // Calculate density gradient
      Real d2 = (SQR(u0(m,IDN,k,j,i+1) - u0(m,IDN,k,j,i-1))
               + SQR(u0(m,IDN,k,j+1,i) - u0(m,IDN,k,j-1,i))
               + SQR(u0(m,IDN,k+1,j,i) - u0(m,IDN,k-1,j,i)));
      ddmax = fmax((sqrt(d2)/u0(m,IDN,k,j,i)), ddmax);

      // Calculate pressure gradient
      Real p2 = (SQR(w0(m,IEN,k,j,i+1) - w0(m,IEN,k,j,i-1))
               + SQR(w0(m,IEN,k,j+1,i) - w0(m,IEN,k,j-1,i))
	       + SQR(w0(m,IEN,k+1,j,i) - w0(m,IEN,k-1,j,i)));
      ddmax = fmax((sqrt(p2)/w0(m,IEN,k,j,i)), ddmax);
    },Kokkos::Max<Real>(team_ddmax));

    if (team_ddmax > ddens_threshold_local) {refine_flag.d_view(m+mbs) = 1;}
    if (team_ddmax < 0.25*ddens_threshold_local) {refine_flag.d_view(m+mbs) = -1;}

  });

  // sync host and device
  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}

//===========================================================================//
//                              History Output                               //
//===========================================================================//

void CGMHistory(HistoryData *pdata, Mesh *pm) {
  auto &w0 = pm->pmb_pack->phydro->w0;
  auto &size = pm->pmb_pack->pmb->mb_size;
  const EOS_Data &eos_data = pm->pmb_pack->phydro->peos->eos_data;
  const int nhydro = pm->pmb_pack->phydro->nhydro;
  const int nmetal = nhydro;
  const int ndusti = nhydro + 1;
  const int ndust = Ndust_bins;
  const bool dust_enabled = dust_model && ndust > 0;

  const int nhist = dust_enabled ? 4 + ndust : 4;
  if (nhist > NHISTORY_VARIABLES) {
    throw std::runtime_error("CGMHistory requires more history slots than "
                             "NHISTORY_VARIABLES allows");
  }
  pdata->nhist = nhist;
  pdata->label[0] = "U^2";
  pdata->label[1] = "Mcold";
  pdata->label[2] = "Mgas";
  pdata->label[3] = "Zgas";
  if (dust_enabled) {
    for (int id = 0; id < ndust; ++id) {
      pdata->label[4 + id] = "D" + std::to_string(id);
    }
  }

  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int is = indcs.is;
  int js = indcs.js;
  int ks = indcs.ks;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;
  const int nmkji = (pm->pmb_pack->nmb_thispack) * nx3 * nx2 * nx1;
  const int nkji = nx3 * nx2 * nx1;
  const int nji = nx2 * nx1;

  const Real gm1 = eos_data.gamma - 1.0;
  const Real kelvin = (pm->pmb_pack->punit != nullptr)
                          ? pm->pmb_pack->punit->temperature_cgs()
                          : 1.0;
  constexpr Real cold_temp_cutoff = 2.0e4;

  array_sum::GlobalSum sum_this_mb;
  for (int n = 0; n < NREDUCTION_VARIABLES; ++n) {
    sum_this_mb.the_array[n] = 0.0;
  }

  Kokkos::parallel_reduce("CGMHistory", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = (idx) / nkji;
    int k = (idx - m * nkji) / nji;
    int j = (idx - m * nkji - k * nji) / nx1;
    int i = (idx - m * nkji - k * nji - j * nx1) + is;
    k += ks;
    j += js;

    const Real vol = size.d_view(m).dx1 * size.d_view(m).dx2 * size.d_view(m).dx3;
    const Real rho = w0(m, IDN, k, j, i);
    const Real mass = rho * vol;
    const Real temp = (w0(m, IEN, k, j, i) * gm1) / rho * kelvin;

    array_sum::GlobalSum hvars;
    hvars.the_array[0] +=
        (SQR(w0(m, IVX, k, j, i)) + SQR(w0(m, IVY, k, j, i)) + SQR(w0(m, IVZ, k, j, i))) * vol;
    hvars.the_array[1] += (temp < cold_temp_cutoff) ? mass : 0.0;
    hvars.the_array[2] += mass;
    hvars.the_array[3] += mass * w0(m, nmetal, k, j, i);

    if (dust_enabled) {
      for (int id = 0; id < ndust; ++id) {
        hvars.the_array[4 + id] += mass * w0(m, ndusti + id, k, j, i);
      }
    }

    for (int n = nhist; n < NHISTORY_VARIABLES; ++n) {
      hvars.the_array[n] = 0.0;
    }
    mb_sum += hvars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  Kokkos::fence();

  for (int n = 0; n < pdata->nhist; ++n) {
    pdata->hdata[n] = sum_this_mb.the_array[n];
  }
}

//===========================================================================//
//                            Post Main Loop                                 //
//===========================================================================//

void FreeProfile(ParameterInput *pin, Mesh *pm) {
  // Free Kokkos views before Kokkos::finalize is called
  profile_reader = ProfileReader();
  disk_profile_reader = ProfileReader();
  sn_centers_buffer = DvceArray2D<Real>();
  sn_counter = Kokkos::View<int>();
  dust_bin_centers = DvceArray1D<Real>();
  dust_bin_edges = DvceArray1D<Real>();
  dust_bin_edges_padded = DvceArray1D<Real>();
  dust_bin_centers_cubed = DvceArray1D<Real>();
  dust_delta1 = DvceArray1D<Real>();
  dust_delta2 = DvceArray1D<Real>();
  dust_delta3 = DvceArray1D<Real>();
  dust_delta4 = DvceArray1D<Real>();
  dust_mass_to_num = DvceArray1D<Real>();
  dust_num_to_mass = DvceArray1D<Real>();
  dust_bin_mass_fractions = DvceArray1D<Real>();
  dust_bin_centers_h = HostArray1D<Real>();
  dust_bin_edges_h = HostArray1D<Real>();
  dust_bin_edges_padded_h = HostArray1D<Real>();
  dust_bin_centers_cubed_h = HostArray1D<Real>();
  dust_edge4_h = HostArray1D<Real>();
  dust_delta1_h = HostArray1D<Real>();
  dust_delta2_h = HostArray1D<Real>();
  dust_delta3_h = HostArray1D<Real>();
  dust_delta4_h = HostArray1D<Real>();
  dust_mass_to_num_h = HostArray1D<Real>();
  dust_num_to_mass_h = HostArray1D<Real>();
  mach_cache = DvceArray4D<Real>();
  mach_cache_cycle = -1;
}
