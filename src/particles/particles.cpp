//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file particles.cpp
//! \brief implementation of Particles class constructor and assorted other functions

#include <iostream>
#include <string>
#include <algorithm>
#include <limits>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "bvals/bvals.hpp"
#include "particles.hpp"

namespace particles {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Particles::Particles(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack) {
  // check this is at least a 2D problem
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle module only works in 2D/3D" <<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // select particle type
  {
    std::string ptype = pin->GetString("particles","particle_type");
    if (ptype.compare("cosmic_ray") == 0) {
      particle_type = ParticleType::cosmic_ray;
    } else if (ptype == "lagrangian_tracer") {
      particle_type = ParticleType::lagrangian_tracer;
    } else if (ptype.compare("lagrangian_mc") == 0) {
      particle_type = ParticleType::lagrangian_mc;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle type = '" << ptype << "' not recognized"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // select pusher algorithm
  {
    std::string ppush = pin->GetString("particles","pusher");
    if (ppush.compare("drift") == 0) {
      pusher = ParticlesPusher::drift;
    } else if (ppush == "lagrangian_tracer") {
      pusher = ParticlesPusher::lagrangian_tracer;
    } else if (ppush.compare("lagrangian_mc") == 0) {
      pusher = ParticlesPusher::lagrangian_mc;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if ((IsLagrangianMC() != (pusher == ParticlesPusher::lagrangian_mc)) ||
      (IsLagrangianTracer() != (pusher == ParticlesPusher::lagrangian_tracer))) {
    std::cerr << "Fluid tracer particle_type and pusher must match" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (IsLagrangianTracer()) {
    auto coord = pmy_pack->pcoord;
    if (!pmy_pack->phydro || pmy_pack->pmhd || pmy_pack->pframe_tracker ||
        coord->is_special_relativistic || coord->is_general_relativistic ||
        coord->is_dynamical_relativistic ||
        pin->GetOrAddString("time","integrator","rk2") != "rk2" ||
        pin->GetOrAddString("time","evolution","dynamic") != "dynamic") {
      std::cerr << "lagrangian_tracer requires dynamic Newtonian hydro with RK2"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int d=0; d<(pmy_pack->pmesh->three_d ? 6 : 4); ++d) {
      auto bc = pmy_pack->pmesh->mesh_bcs[d];
      if (bc != BoundaryFlag::periodic && bc != BoundaryFlag::outflow) {
        std::cerr << "lagrangian_tracer requires periodic/outflow boundaries"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // set dimensions of particle arrays. Note particles only work in 2D/3D
  if (pmy_pack->pmesh->one_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particles only work in 2D/3D, but 1D problem initialized" <<std::endl;
    std::exit(EXIT_FAILURE);
  }
  switch (particle_type) {
    case ParticleType::cosmic_ray:
      {
        Real ppc = pin->GetOrAddReal("particles","ppc",1.0);
        auto &indcs = pmy_pack->pmesh->mb_indcs;
        int ncells = indcs.nx1*indcs.nx2*indcs.nx3;
        Real r_npart = ppc*static_cast<Real>((pmy_pack->nmb_thispack)*ncells);
        nprtcl_thispack = static_cast<int>(r_npart);
        int ndim=4;
        if (pmy_pack->pmesh->three_d) {ndim+=2;}
        nrdata = ndim;
        nidata = 2;
        break;
      }
    case ParticleType::lagrangian_tracer:
    case ParticleType::lagrangian_mc:
      {
        if (pmy_pack->phydro == nullptr && pmy_pack->pmhd == nullptr) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "lagrangian_mc particles require hydro or mhd"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        EquationOfState *peos = (pmy_pack->phydro != nullptr) ?
                                pmy_pack->phydro->peos : pmy_pack->pmhd->peos;
        if (!peos->eos_data.is_ideal) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "lagrangian_mc thermo tracers require an ideal-gas "
                    << "EOS in this implementation" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        nprtcl_thispack = 0;
        nrdata = LMC_NREAL;
        nidata = PSEEDID + 1;
        random_seed = pin->GetOrAddInteger("particles","random_seed",12345);
        next_tracer_tag = pin->GetOrAddInteger("particles","next_tracer_tag",0);
        ParseTracerSeedSchedules(pin);
        if (IsLagrangianMC() && pmy_pack->phydro != nullptr) {
          pmy_pack->phydro->SetSaveUFlxIdn();
        } else if (IsLagrangianMC()) {
          pmy_pack->pmhd->SetSaveUFlxIdn();
        }
        dtnew = std::numeric_limits<Real>::max();
        break;
      }
    default:
      break;
  }
  Kokkos::realloc(prtcl_rdata, nrdata, nprtcl_thispack);
  Kokkos::realloc(prtcl_idata, nidata, nprtcl_thispack);

  ConfigureShockTracking(pin);

  // allocate boundary object
  pbval_part = new ParticlesBoundaryValues(this, pin);
}

//----------------------------------------------------------------------------------------
// destructor

Particles::~Particles() {
}

//----------------------------------------------------------------------------------------
// GetLagrangianMCScalarCount()

int Particles::GetLagrangianMCScalarCount() const {
  if (pmy_pack->phydro != nullptr) return pmy_pack->phydro->nscalars;
  if (pmy_pack->pmhd != nullptr) return pmy_pack->pmhd->nscalars;
  return 0;
}

//----------------------------------------------------------------------------------------
// CreateParticleTags()
// Assigns tags to particles (unique integer).  Note that tracked particles are always
// those with tag numbers less than ntrack.

void Particles::CreateParticleTags(ParameterInput *pin) {
  if (IsFluidTracer()) return;

  std::string assign = pin->GetOrAddString("particles","assign_tag","index_order");

  // tags are assigned sequentially within this rank, starting at 0 with rank=0
  if (assign.compare("index_order") == 0) {
    int tagstart = 0;
    for (int n=1; n<=global_variable::my_rank; ++n) {
      tagstart += pmy_pack->pmesh->nprtcl_eachrank[n-1];
    }

    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = tagstart + p;
    });

  // tags are assigned sequentially across ranks
  } else if (assign.compare("rank_order") == 0) {
    int myrank = global_variable::my_rank;
    int nranks = global_variable::nranks;
    auto &pi = prtcl_idata;
    par_for("ptags",DevExeSpace(),0,(nprtcl_thispack-1),
    KOKKOS_LAMBDA(const int p) {
      pi(PTAG,p) = myrank + nranks*p;
    });

  // tag algorithm not recognized, so quit with error
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Particle tag assignment type = '" << assign << "' not recognized"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace particles
