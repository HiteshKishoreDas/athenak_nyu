// Copyright(C) 2026 Athena code contributors
// Independent Schmidt--Federrath closure checks. No expected value calls Stress().
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "hydro/sgs_closure.hpp"

namespace {
// Frozen before execution: conditioned double-precision component tolerance.
constexpr long double tolerance = 5.0e-13L;
using Matrix = std::array<std::array<long double, 3>, 3>;
struct Case {
  Matrix a{};
  double rho = 1.0, energy = 1.0, delta = 1.0;
  hydro::sgs::Coefficients coefficients;
};

Matrix Oracle(const Case &c) {
  Matrix t{};
  long double q = 0.0L;
  long double norm = 0.0L;
  long double d = c.a[0][0] + c.a[1][1] + c.a[2][2];
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      q += c.a[i][j]*c.a[i][j];
      norm = std::hypot(norm, c.a[i][j]);
    }
  }
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      if (norm == 0.0L) {
        t[i][j] = (i == j) ? -2.0L*c.energy/3.0L : 0.0L;
        continue;
      }
      long double b = 0.0L;
      for (int l=0; l<3; ++l) {
        // Direct published quadratic ratio for moderate inputs. Independent
        // hypot normalization covers norm-squared underflow on arm64, where
        // long double has no extra exponent range.
        if (q == 0.0L) {
          b += (c.a[i][l]/norm)*(c.a[j][l]/norm);
        } else {
          b += c.a[i][l]*c.a[j][l]/q;
        }
      }
      long double strain = (c.a[i][j]+c.a[j][i])/2.0L;
      if (i == j) strain -= d/3.0L;
      t[i][j] = 2.0L*c.coefficients.c1*c.delta*
          std::sqrt(2.0L*c.rho*c.energy)*strain -
          2.0L*c.coefficients.c2*c.energy*b;
      if (i == j) t[i][j] -=
          2.0L*(1.0L-c.coefficients.c2)*c.energy/3.0L;
    }
  }
  return t;
}

long double Production(const Matrix &a, const Matrix &t) {
  long double s = 0.0L;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) s += t[i][j]*(a[i][j]+a[j][i])/2.0L;
  }
  return s;
}

int checks = 0;
long double max_error = 0.0L;
void Check(long double actual, long double expected, long double scale,
           const char *name) {
  ++checks;
  scale = std::fmax(std::fabs(scale), 1.0e-300L);
  long double error = std::fabs(actual-expected)/scale;
  max_error = std::fmax(max_error, error);
  if (!std::isfinite(actual) || error > tolerance) {
    std::cerr << name << ": actual=" << static_cast<double>(actual)
              << " expected=" << static_cast<double>(expected)
              << " scaled_error=" << static_cast<double>(error) << '\n';
    throw std::runtime_error("Independent SGS oracle mismatch");
  }
}

Matrix Rotate(const Matrix &a, const Matrix &r) {
  Matrix b{};
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      for (int p=0; p<3; ++p) {
        for (int q=0; q<3; ++q) b[i][j] += r[i][p]*a[p][q]*r[j][q];
      }
    }
  }
  return b;
}

void Run() {
  static_assert(sizeof(Real) == sizeof(double), "Oracle requires double precision");
  std::vector<Case> cases(8);
  cases[1].energy = 0.0;
  cases[1].a[0][1] = 1.0L;
  for (int i=0; i<3; ++i) cases[2].a[i][i] = -0.25L;
  cases[3].a[0][1] = 1.0L;
  cases[4].a[0][1] = -0.75L;
  cases[4].a[1][0] = 0.75L;
  cases[5].a[0][0] = 0.5L;
  cases[5].a[1][1] = -0.5L;
  cases[6].a[0][0] = 0.2L;
  cases[6].a[1][1] = cases[6].a[2][2] = -0.1L;
  cases[7].a[0][0] = 1.0L;
  // Hand production checks constrain the oracle itself, not just production code.
  Check(Production(cases[2].a, Oracle(cases[2])), 0.5L, 1.0L, "compression");
  Check(Production(cases[3].a, Oracle(cases[3])),
        0.02L*std::sqrt(2.0L), 1.0L, "simple shear");
  Check(Production(cases[4].a, Oracle(cases[4])), 0.0L, 1.0L, "rotation");
  Check(Production(cases[5].a, Oracle(cases[5])),
        0.02L*std::sqrt(2.0L), 1.0L, "traceless strain");
  Check(Production(cases[6].a, Oracle(cases[6])),
        0.0024L*std::sqrt(2.0L)-0.14L, 1.0L, "backscatter");

  std::mt19937_64 random(4492);
  std::uniform_real_distribution<double> gradient(-3.0, 3.0), logstate(-3.0, 3.0);
  for (int n=0; n<1000; ++n) {
    Case c;
    c.rho = std::pow(10.0, logstate(random));
    c.energy = std::pow(10.0, logstate(random));
    c.delta = std::pow(10.0, logstate(random));
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) c.a[i][j] = gradient(random);
    }
    cases.push_back(c);
  }
  for (double c1 : {0.0, 0.02}) {
    for (double c2 : {0.0, 0.7, 1.0}) {
      Case c = cases[6];
      c.coefficients.c1 = c1;
      c.coefficients.c2 = c2;
      c.coefficients.ceps = c.coefficients.ckappa = 0.0;
      cases.push_back(c);
    }
  }
  for (double factor : {1.0e-280, 1.0e-150, 1.0e-20, 1.0e150}) {
    Case c = cases[6];
    for (auto &row : c.a) {
      for (auto &entry : row) entry *= factor;
    }
    cases.push_back(c);
  }
  for (double rho : {1.0e-100, 1.0e100}) {
    for (double energy : {1.0e-100, 1.0e100}) {
      Case c = cases[6];
      c.rho = rho;
      c.energy = energy;
      cases.push_back(c);
    }
  }
  Matrix rotation = {
    {
      {0.36L, -0.48L, 0.8L},
      {0.8L, 0.6L, 0.0L},
      {-0.48L, 0.64L, 0.6L}
    }
  };
  Case rotated = cases[100];
  rotated.a = Rotate(rotated.a, rotation);
  int rotation_index = static_cast<int>(cases.size());
  cases.push_back(rotated);
  Case scaling = cases[6];
  scaling.coefficients.c2 = 0.0;
  int scaling_index = static_cast<int>(cases.size());
  cases.push_back(scaling);
  scaling.rho = 9.0;
  scaling.energy = 4.0;
  scaling.delta = 2.0;
  cases.push_back(scaling);

  int count = static_cast<int>(cases.size());
  Kokkos::View<Real***> a("gradients", count, 3, 3);
  Kokkos::View<Real**> state("states", count, 7), rates("rates", count, 2);
  Kokkos::View<Real***> stress("stress", count, 3, 3);
  auto ah = Kokkos::create_mirror_view(a);
  auto sh = Kokkos::create_mirror_view(state);
  for (int n=0; n<count; ++n) {
    const auto &c = cases[n];
    sh(n,0) = c.rho; sh(n,1) = c.energy; sh(n,2) = c.delta;
    sh(n,3) = c.coefficients.c1; sh(n,4) = c.coefficients.c2;
    sh(n,5) = c.coefficients.ceps; sh(n,6) = c.coefficients.ckappa;
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) ah(n,i,j) = static_cast<Real>(c.a[i][j]);
    }
  }
  Kokkos::deep_copy(a, ah);
  Kokkos::deep_copy(state, sh);
  Kokkos::parallel_for("SGS independent oracle", count, KOKKOS_LAMBDA(int n) {
    Real gradient[3][3], tau[3][3];
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) gradient[i][j] = a(n,i,j);
    }
    hydro::sgs::Coefficients c;
    c.c1 = state(n,3); c.c2 = state(n,4);
    c.ceps = state(n,5); c.ckappa = state(n,6);
    hydro::sgs::Stress(state(n,0), state(n,1), state(n,2), gradient, c, tau);
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) stress(n,i,j) = tau[i][j];
    }
    rates(n,0) = hydro::sgs::Dissipation(state(n,0), state(n,1), state(n,2), c);
    rates(n,1) = hydro::sgs::Diffusivity(state(n,0), state(n,1), state(n,2), c);
  });
  auto th = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), stress);
  auto rh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rates);
  for (int n=0; n<count; ++n) {
    // Expected input is the exact double value supplied to the device.
    Case c = cases[n];
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) c.a[i][j] = ah(n,i,j);
    }
    Matrix expected = Oracle(c), actual{};
    long double scale = c.energy;
    for (const auto &row : expected) {
      for (auto entry : row) scale = std::fmax(scale, std::fabs(entry));
    }
    long double trace = 0.0L;
    for (int i=0; i<3; ++i) {
      trace += th(n,i,i);
      for (int j=0; j<3; ++j) {
        actual[i][j] = th(n,i,j);
        Check(actual[i][j], expected[i][j], scale, "stress component");
        Check(th(n,i,j), th(n,j,i), scale, "symmetry");
      }
    }
    Check(trace, -2.0L*c.energy, scale, "trace");
    Check(Production(c.a, actual), Production(c.a, expected),
          std::fmax(std::fabs(Production(c.a, expected)), scale), "production");
    long double diss = c.coefficients.ceps*
        std::pow(static_cast<long double>(c.energy), 1.5L)/std::sqrt(c.rho)/c.delta;
    long double diff = c.coefficients.ckappa*c.delta*
        std::sqrt(static_cast<long double>(c.rho)*c.energy);
    Check(rh(n,0), diss, std::fabs(diss), "dissipation");
    Check(rh(n,1), diff, std::fabs(diff), "diffusivity");
  }
  Matrix unrotated{};
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) unrotated[i][j] = th(100,i,j);
  }
  Matrix covariant = Rotate(unrotated, rotation);
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) Check(th(rotation_index,i,j), covariant[i][j],
        std::fmax(1.0L, std::fabs(covariant[i][j])), "rotation covariance");
  }
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      long double original = th(scaling_index,i,j);
      long double scaled = th(scaling_index+1,i,j);
      // Remove isotropic pressure to isolate eddy density/energy/width powers.
      if (i == j) { original += 2.0L/3.0L; scaled += 8.0L/3.0L; }
      Check(scaled, 12.0L*original, 4.0L, "eddy scaling");
    }
  }
  Check(rh(scaling_index+1,0), (4.0L/3.0L)*rh(scaling_index,0),
        2.0L, "dissipation scaling");
  Check(rh(scaling_index+1,1), 12.0L*rh(scaling_index,1),
        8.0L, "diffusivity scaling");
  std::cout << std::setprecision(17) << "cases=" << count << " assertions="
            << checks << " max_scaled_error=" << static_cast<double>(max_error)
            << " tolerance=" << static_cast<double>(tolerance) << '\n';
}
}  // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  int result = 0;
  try {
    Run();
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    result = 1;
  }
  Kokkos::finalize();
  return result;
}
