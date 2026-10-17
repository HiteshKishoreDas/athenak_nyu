//========================================================================================
// Copyright (C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License. See LICENSE for details.
// AthenaK material-curve kinematic validation fixture.
// Prescribed fields are imposed through conserved source updates. Fluid boundary fills,
// primitive conversion, mesh interpolation and the ordinary curve RK2 path still execute.
// Prescribed flows are kinematic; unforced constant-pressure shear is tested separately.
//========================================================================================
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "mesh/mesh_refinement.hpp"
#include "pgen.hpp"

namespace {
struct Field {
  int mode;
  Real velocity[3], acceleration[3], axis[3], rates[3];
  Real frequency, amplitude, wavenumber, signal_speed, density_amplitude;
  Real constant_pressure;
};
Field field;
Real refine_time, coarsen_time;

void SetField(Mesh *pm, Real time) {
  auto *pack = pm->pmb_pack;
  const auto indcs = pm->mb_indcs;
  auto &u = pack->phydro->u0;
  auto &size = pack->pmb->mb_size;
  const auto eos = pack->phydro->peos->eos_data;
  const Field f = field;
  par_for("material_curve_prescribed_field", DevExeSpace(), 0, pack->nmb_thispack-1,
          indcs.ks, indcs.ke, indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const auto s = size.d_view(m);
    Real x[3] = {CellCenterX(i-indcs.is, indcs.nx1, s.x1min, s.x1max),
                 CellCenterX(j-indcs.js, indcs.nx2, s.x2min, s.x2max),
                 CellCenterX(k-indcs.ks, indcs.nx3, s.x3min, s.x3max)};
    Real v[3] = {f.velocity[0], f.velocity[1], f.velocity[2]};
    if (f.mode == 1) {
      for (int d=0; d<3; ++d) v[d] += time*f.acceleration[d];
    } else if (f.mode == 2) {
      for (int d=0; d<3; ++d) v[d] += f.acceleration[d]*sin(f.frequency*time);
    } else if (f.mode == 3) {
      v[0] += f.amplitude*(f.axis[1]*x[2]-f.axis[2]*x[1]);
      v[1] += f.amplitude*(f.axis[2]*x[0]-f.axis[0]*x[2]);
      v[2] += f.amplitude*(f.axis[0]*x[1]-f.axis[1]*x[0]);
    } else if (f.mode == 4) {
      for (int d=0; d<3; ++d) v[d] += f.rates[d]*x[d];
    } else if (f.mode == 5) {
      v[0] += f.amplitude*x[1];
    } else if (f.mode == 6) {
      v[0] += f.amplitude*sin(f.wavenumber*x[1])*cos(f.frequency*time);
    }
    Real rho = 1.0 + f.density_amplitude*sin(f.wavenumber*x[0]);
    Real vmax = fmax(fabs(v[0]), fmax(fabs(v[1]), fabs(v[2])));
    Real cs = f.signal_speed-vmax;
    u(m,IDN,k,j,i) = rho;
    Real energy = 0.0;
    for (int d=0; d<3; ++d) {
      u(m,IM1+d,k,j,i) = rho*v[d];
      energy += 0.5*rho*v[d]*v[d];
    }
    if (eos.is_ideal) {
      u(m,IEN,k,j,i) = (f.constant_pressure > 0.0 ?
          f.constant_pressure/(eos.gamma-1.0) :
          rho*cs*cs/(eos.gamma*(eos.gamma-1.0))) + energy;
    }
  });
}

void Prescribe(Mesh *pm, const Real bdt) {
  // SSPRK2's predictor and completed state both correspond to t+dt. The source
  // callback runs after RKUpdate and before the boundary/primitive refresh.
  SetField(pm, pm->time+pm->dt);
}

void TimedRefinement(MeshBlockPack *pack) {
  auto *mesh = pack->pmesh;
  auto &flags = mesh->pmr->refine_flag;
  const int first = mesh->gids_eachrank[global_variable::my_rank];
  const int decision = mesh->time >= coarsen_time ? -1 :
                       (mesh->time >= refine_time ? 1 : 0);
  for (int m=0; m<pack->nmb_thispack; ++m) flags.h_view(first+m) = decision;
  flags.template modify<HostMemSpace>();
  flags.template sync<DevExeSpace>();
}
}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (pmy_mesh_->pmb_pack->phydro == nullptr) {
    std::cerr << "material_curve_test requires Newtonian hydro\n";
    std::exit(EXIT_FAILURE);
  }
  std::string mode = pin->GetOrAddString("problem", "flow", "translation");
  if (mode=="shock_tube") {
    // Native evolved RH fixture with the existing timed AMR callback; no prescribed source.
    user_srcs=false;
    user_srcs_func=nullptr;
    refine_time=pin->GetOrAddReal("problem","refine_time",1.0e30);
    coarsen_time=pin->GetOrAddReal("problem","coarsen_time",1.0e30);
    if (pmy_mesh_->adaptive) user_ref_func=TimedRefinement;
    ShockTube(pin,restart);
    return;
  }
  field.mode = -1;
  const char *names[] = {"translation", "acceleration", "sinusoidal", "rotation",
                         "strain", "shear", "reversing_shear"};
  for (int i=0; i<7; ++i) {
    if (mode == names[i]) field.mode = i;
  }
  if (mode == "zero") field.mode = 0;
  if (field.mode < 0) {
    std::cerr << "material_curve_test: unknown flow " << mode << '\n';
    std::exit(EXIT_FAILURE);
  }
  for (int d=0; d<3; ++d) {
    std::string suffix = std::to_string(d+1);
    field.velocity[d] = pin->GetOrAddReal("problem", "v"+suffix, 0.0);
    field.acceleration[d] = pin->GetOrAddReal("problem", "a"+suffix, 0.0);
    field.axis[d] = pin->GetOrAddReal("problem", "axis"+suffix, d == 2 ? 1.0 : 0.0);
    field.rates[d] = pin->GetOrAddReal("problem", "s"+suffix, 0.0);
  }
  Real axis_norm = std::sqrt(SQR(field.axis[0])+SQR(field.axis[1])+SQR(field.axis[2]));
  if (!(axis_norm > 0.0)) {
    std::cerr << "material_curve_test: zero rotation axis\n";
    std::exit(EXIT_FAILURE);
  }
  for (int d=0; d<3; ++d) field.axis[d] /= axis_norm;
  field.frequency = pin->GetOrAddReal("problem", "frequency", 1.0);
  field.amplitude = pin->GetOrAddReal("problem", "amplitude", 0.2);
  field.wavenumber = pin->GetOrAddReal("problem", "wavenumber", 1.0);
  field.signal_speed = pin->GetOrAddReal("problem", "signal_speed", 4.0);
  field.density_amplitude = pin->GetOrAddReal("problem", "density_amplitude", 0.0);
  field.constant_pressure = pin->GetOrAddReal("problem", "constant_pressure", 0.0);
  if (field.constant_pressure < 0 || !std::isfinite(field.constant_pressure)) {
    std::cerr << "material_curve_test: constant_pressure must be finite "
              << "and nonnegative\n";
    std::exit(EXIT_FAILURE);
  }
  refine_time = pin->GetOrAddReal("problem", "refine_time", 1.0e30);
  coarsen_time = pin->GetOrAddReal("problem", "coarsen_time", 1.0e30);
  if (user_srcs) user_srcs_func = Prescribe;
  if (pmy_mesh_->adaptive) user_ref_func = TimedRefinement;
  if (!restart) SetField(pmy_mesh_, pmy_mesh_->time);
}
