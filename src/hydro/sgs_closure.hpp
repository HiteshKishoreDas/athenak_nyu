#ifndef HYDRO_SGS_CLOSURE_HPP_
#define HYDRO_SGS_CLOSURE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sgs_closure.hpp
//! \brief Schmidt & Federrath (2011) hydrodynamic SGS closures, Eqs. (8), (10), (11).

#include <cmath>

#include "athena.hpp"

namespace hydro {
namespace sgs {

struct Coefficients {
  Real c1 = 0.02;
  Real c2 = 0.7;
  Real ceps = 1.5;
  Real ckappa = 0.65;
};

// Inputs must be finite, with rho > 0, K >= 0 and delta > 0. K is energy density,
// and a[i][j] = d(u_i)/d(x_j). The solver checks admissibility before these helpers.
KOKKOS_INLINE_FUNCTION
void Stress(const Real rho, const Real K, const Real delta, const Real a[3][3],
            const Coefficients &c, Real tau[3][3]) {
  if (K == 0.0) {
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        tau[i][j] = 0.0;
      }
    }
    return;
  }

  // Normalize the gradient before forming its quadratic norm. This preserves
  // the direction-dependent nonlinear term for tiny nonzero gradients.
  Real scale = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      scale = fmax(scale, fabs(a[i][j]));
    }
  }
  if (scale == 0.0) {
    // Explicit zero-gradient convention; the nonlinear limit is not unique.
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        tau[i][j] = (i == j) ? -(2.0/3.0)*K : 0.0;
      }
    }
    return;
  }

  Real b[3][3];
  Real q = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      b[i][j] = a[i][j]/scale;
      q += b[i][j]*b[i][j];
    }
  }
  const Real d = (b[0][0] + b[1][1] + b[2][2])/3.0;
  const Real eddy = 2.0*c.c1*delta*sqrt(2.0)*sqrt(rho)*sqrt(K);
  const Real isotropic = -(2.0/3.0)*(1.0-c.c2)*K;
  for (int i=0; i<3; ++i) {
    for (int j=i; j<3; ++j) {
      Real s = 0.5*(b[i][j] + b[j][i]);
      if (i == j) {
        s -= d;
      }
      Real aat = 0.0;
      for (int k=0; k<3; ++k) {
        aat += b[i][k]*b[j][k];
      }
      Real t = (eddy*scale)*s - 2.0*c.c2*K*(aat/q);
      if (i == j) {
        t += isotropic;
      }
      tau[i][j] = t;
      tau[j][i] = t;
    }
  }
}

// Dissipation is energy density per time, deposited into resolved thermal energy.
KOKKOS_INLINE_FUNCTION
Real Dissipation(const Real rho, const Real K, const Real delta,
                 const Coefficients &c) {
  return c.ceps*K*(sqrt(K)/sqrt(rho))/delta;
}

// This coefficient multiplies the gradient of specific SGS energy K/rho.
KOKKOS_INLINE_FUNCTION
Real Diffusivity(const Real rho, const Real K, const Real delta,
                 const Coefficients &c) {
  return c.ckappa*delta*sqrt(rho)*sqrt(K);
}

}  // namespace sgs
}  // namespace hydro
#endif  // HYDRO_SGS_CLOSURE_HPP_
