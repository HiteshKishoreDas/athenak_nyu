// Copyright(C) 2026 Athena code contributors
// Evaluate the production closure on filtered states; no solver evolution.
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#include "hydro/sgs_closure.hpp"

int main(int argc, char *argv[]) {
  static_assert(sizeof(Real) == sizeof(double), "double precision required");
  if (argc != 3) {
    std::cerr << "Usage: sgs_filtered_stress "
              << "input-double-records output-double-records\n";
    return 1;
  }
  std::ifstream input(argv[1], std::ios::binary);
  std::ofstream output(argv[2], std::ios::binary);
  if (!input || !output) return 1;
  const hydro::sgs::Coefficients coefficients;
  std::array<double, 12> record;
  while (input.read(reinterpret_cast<char *>(record.data()), sizeof(record))) {
    for (double value : record) {
      if (!std::isfinite(value)) return 1;
    }
    if (record[0] <= 0.0 || record[1] < 0.0 || record[2] <= 0.0) return 1;
    Real gradient[3][3], tau[3][3];
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) gradient[i][j] = record[3+3*i+j];
    }
    hydro::sgs::Stress(record[0], record[1], record[2], gradient, coefficients, tau);
    std::array<double, 10> result{};
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        result[3*i+j] = tau[i][j];
        result[9] += tau[i][j]*gradient[i][j];
      }
    }
    for (double value : result) {
      if (!std::isfinite(value)) return 1;
    }
    output.write(reinterpret_cast<const char *>(result.data()), sizeof(result));
    if (!output) return 1;
  }
  if (!input.eof() || input.gcount() != 0) return 1;
  return 0;
}
