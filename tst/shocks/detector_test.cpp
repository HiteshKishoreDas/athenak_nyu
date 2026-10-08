//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include "shocks/detector.hpp"

namespace {
Real ReadReal(std::istream &in) {
  std::string token;
  if (!(in >> token)) throw std::runtime_error("Incomplete fixture");
  std::size_t end=0;
  const Real value=std::stod(token,&end);
  if (end != token.size()) throw std::runtime_error("Invalid real token");
  return value;
}

int Wrap(int a, int n) {
  return (a%n+n)%n;
}

void PrintSignal(int i, int j, int k, const shocks::Signal &s) {
  std::cout << i << ',' << j << ',' << k << ',' << s.valid << ',';
  std::cout << s.candidate << ',' << s.inside << ',' << s.rho << ',';
  std::cout << s.pressure << ',' << s.temperature << ',' << s.divergence;
  for (int d=0; d<3; ++d) std::cout << ',' << s.grad_temperature[d];
  for (int d=0; d<3; ++d) std::cout << ',' << s.grad_entropy[d];
  for (int d=0; d<3; ++d) std::cout << ',' << s.grad_density[d];
  for (int d=0; d<3; ++d) std::cout << ',' << s.normal[d];
  std::cout << ',' << s.rho_up << ',' << s.rho_down;
  std::cout << ',' << s.pressure_up << ',' << s.pressure_down;
  std::cout << ',' << s.temperature_up << ',' << s.temperature_down;
  std::cout << ',' << s.mach_pressure << ',' << s.mach_temperature << '\n';
}

void RunChecks() {
  DvceArray5D<Real> w("checks",1,5,1,3,3);
  Kokkos::deep_copy(w,1.0);
  DvceArray1D<int> checks("checks",8);
  Kokkos::parallel_for("detector_boundary_checks",1,KOKKOS_LAMBDA(int) {
    shocks::DetectorConfig c;
    const Real cap2=c.mach_max*c.mach_max;
    const Real q=(2*c.gamma*cap2-(c.gamma-1))/(c.gamma+1);
    const Real r=(c.gamma+1)*cap2/((c.gamma-1)*cap2+2);
    const Real cap=q/r;
    checks(0)=Kokkos::isnan(shocks::TemperatureMach(1,c));
    checks(1)=Kokkos::isnan(shocks::TemperatureMach(
        Kokkos::nextafter(cap,std::numeric_limits<Real>::infinity()),c));
    checks(2)=Kokkos::fabs(shocks::TemperatureMach(cap,c)-c.mach_max) < 1e-9;
    checks(3)=!shocks::DetectCell(w,0,0,0,0,1,1,1,false,c).valid;
    c.width=0;
    checks(4)=!shocks::DetectCell(w,0,0,1,1,1,1,1,false,c).valid;
    const Real a[3]={1e-200,1e-200,0}, b[3]={1e200,-1e200,0};
    const Real zero[3]={0,0,0};
    checks(5)=shocks::Alignment(a,b) == 0;
    checks(6)=shocks::Alignment(a,zero) == 0;
    checks(7)=shocks::Alignment(a,a) > 0;
  });
  const auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),checks);
  for (int n=0; n<8; ++n) {
    if (!host(n)) throw std::runtime_error("Device self-check "+std::to_string(n));
  }
  std::cout << "8 device numerical/support self-checks PASS\n";
}

void RunFixture(const char *path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Cannot read fixture");
  int nx, ny, nz;
  if (!(in >> nx >> ny >> nz) || nx < 1 || ny < 1 || nz < 1) {
    throw std::runtime_error("Invalid fixture shape");
  }
  shocks::DetectorConfig c;
  c.gamma=ReadReal(in);
  c.min_mach=ReadReal(in);
  c.mach_max=ReadReal(in);
  c.grad_floor=ReadReal(in);
  c.div_limit=ReadReal(in);
  if (!(in >> c.width)) throw std::runtime_error("Invalid fixture width");
  c.density_floor=ReadReal(in);
  c.temperature_floor=ReadReal(in);
  const Real dx=ReadReal(in), dy=ReadReal(in), dz=ReadReal(in);
  if (!shocks::ConfigValid(c) || c.width > 16) {
    throw std::runtime_error("Unsupported fixture configuration");
  }
  const int ng=c.width;
  const bool three_d=nz > 1;
  const int nk=three_d ? nz+2*ng : 1;
  DvceArray5D<Real> w("fixture",1,5,nk,ny+2*ng,nx+2*ng);
  auto host=Kokkos::create_mirror_view(w);
  HostArray5D<Real> core("core",1,5,nz,ny,nx);
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
        for (int n=0; n<5; ++n) core(0,n,k,j,i)=ReadReal(in);
      }
    }
  }
  // This fixture adapter provides genuine periodic neighbors for reference parity.
  // DetectCell itself never wraps coordinates or fills ghosts.
  for (int k=0; k<nk; ++k) {
    for (int j=0; j<ny+2*ng; ++j) {
      for (int i=0; i<nx+2*ng; ++i) {
        for (int n=0; n<5; ++n) {
          host(0,n,k,j,i)=core(0,n,three_d ? Wrap(k-ng,nz) : 0,
                              Wrap(j-ng,ny),Wrap(i-ng,nx));
        }
      }
    }
  }
  Kokkos::deep_copy(w,host);
  DvceArray1D<shocks::Signal> signals("signals",nx*ny*nz);
  Kokkos::parallel_for("detector_fixture",nx*ny*nz,KOKKOS_LAMBDA(int l) {
    const int k=l%nz, j=(l/nz)%ny, i=l/(nz*ny);
    signals(l)=shocks::DetectCell(w,0,three_d ? k+ng : 0,j+ng,i+ng,
                                dx,dy,dz,three_d,c);
  });
  const auto result=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),signals);
  const auto after=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),w);
  for (std::size_t l=0; l<host.size(); ++l) {
    const Real a=host.data()[l], b=after.data()[l];
    if (!(a == b || (std::isnan(a) && std::isnan(b)))) {
      throw std::runtime_error("Detector mutated primitive state");
    }
  }
  std::cout << "# i,j,k,valid,candidate,inside,rho,pressure,temperature,divergence,"
      "grad_tx,grad_ty,grad_tz,grad_sx,grad_sy,grad_sz,grad_rx,grad_ry,grad_rz,"
      "normal_x,normal_y,normal_z,rho_up,rho_down,pressure_up,pressure_down,"
      "temperature_up,temperature_down,mach_pressure,mach_temperature\n";
  std::cout << std::setprecision(17);
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) PrintSignal(i,j,k,result((i*ny+j)*nz+k));
    }
  }
}
}  // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: shock_detector_test fixture.txt\n";
    return 2;
  }
  Kokkos::initialize();
  int status=0;
  try {
    if (std::string(argv[1]) == "--self-check") {
      RunChecks();
    } else {
      RunFixture(argv[1]);
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    status=1;
  }
  Kokkos::finalize();
  return status;
}
