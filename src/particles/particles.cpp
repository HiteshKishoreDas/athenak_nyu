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
#include <cmath>
#include <sstream>
#include <set>
#include <utility>
#include "shock_version.hpp"

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
    } else if (ppush.compare("lagrangian_mc") == 0) {
      pusher = ParticlesPusher::lagrangian_mc;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Particle pusher must be specified in <particles> block"
                <<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (particle_type == ParticleType::lagrangian_mc &&
      pusher != ParticlesPusher::lagrangian_mc) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "particle_type=lagrangian_mc requires "
              << "pusher=lagrangian_mc" << std::endl;
    std::exit(EXIT_FAILURE);
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
        if (pmy_pack->phydro != nullptr) {
          pmy_pack->phydro->SetSaveUFlxIdn();
        } else {
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

  if (pin->DoesBlockExist("shock_tracking")) {
    shock_enabled = pin->GetOrAddBoolean("shock_tracking", "enabled", false);
    if (!shock_enabled &&
        pin->GetOrAddBoolean("shock_tracking", "history_present", false)) {
      std::cerr << "Shock history cannot be disabled on restart" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if (shock_enabled) {
    if (pusher != ParticlesPusher::lagrangian_mc || !IsLagrangianMC()) {
      std::cerr << "Shock tracking requires lagrangian_mc particles"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pin->GetOrAddString("time", "integrator", "rk2") != "rk2" ||
        (pin->DoesBlockExist("cooling") &&
         pin->GetOrAddBoolean("cooling", "enabled", false))) {
      std::cerr << "Shock MC qualification requires RK2 and cooling disabled"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (auto block=pin->block.begin(); block != pin->block.end(); ++block) {
      if (block->block_name.rfind("output",0) == 0 &&
          pin->DoesParameterExist(block->block_name,"file_type") &&
          pin->GetString(block->block_name,"file_type") == "rst" &&
          pin->GetOrAddBoolean(block->block_name,"single_file_per_node",false)) {
        std::cerr << "Shock MC restart does not support per-node files" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    std::set<int> schedule_ids;
    for (const auto &seed : seed_schedules_) {
      if (!schedule_ids.insert(seed.id).second) {
        std::cerr << "Shock restart requires unique tracer seed schedule IDs"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    auto &c = shock_config;
    c.width = pin->GetOrAddInteger("shock_tracking", "width", 1);
    c.min_mach = pin->GetOrAddReal("shock_tracking", "min_mach", 1.1);
    c.mach_max = pin->GetOrAddReal("shock_tracking", "mach_max", 1000);
    c.grad_floor = pin->GetOrAddReal("shock_tracking", "grad_floor", 1e-30);
    c.div_limit = pin->GetOrAddReal("shock_tracking", "div_limit", 0);
    c.density_floor = pin->GetOrAddReal("shock_tracking", "density_floor", 0);
    c.temperature_floor = pin->GetOrAddReal("shock_tracking", "temperature_floor", 0);
    shock_capacity = pin->GetOrAddInteger("shock_tracking", "event_capacity", 4096);
    shock_lineage = pin->GetString("shock_tracking", "lineage");
    shock_directory = pin->GetOrAddString("shock_tracking", "directory", "shock_events");
    shock_trace = pin->GetOrAddBoolean("shock_tracking", "trace", false);
    bool ok = c.width >= 1 && c.width <= pmy_pack->pmesh->mb_indcs.ng &&
              shock_capacity > 0 && !shock_lineage.empty() &&
              shock_lineage.size() <= 128 && !shock_directory.empty();
    for (char ch : shock_lineage) {
      ok = ok && ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
                  (ch >= '0' && ch <= '9') || ch == '-' || ch == '_');
    }
    for (Real v : {c.min_mach, c.mach_max, c.grad_floor, c.div_limit,
                   c.density_floor, c.temperature_floor}) ok = ok && std::isfinite(v);
    ok = ok && c.min_mach >= 1 && c.mach_max > 1 && c.grad_floor >= 0 &&
         c.div_limit <= 0 && c.density_floor >= 0 && c.temperature_floor >= 0;
    if (!ok) {
      std::cerr << "Invalid shock_tracking configuration: finite floors and "
                << "nonpositive divergence threshold, "
                << "width in [1,nghost], positive capacity and safe lineage required"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    c.gamma = pmy_pack->phydro ? pmy_pack->phydro->peos->eos_data.gamma : 0;
    std::ostringstream contract;
    contract.precision(std::numeric_limits<Real>::max_digits10);
    contract << "v1-entry-upstream-logdistance-mc-massflux-rk2-temperature-"
             << "nearest-reference6f5ec8bc "
             << SHOCK_IMPLEMENTATION_HASH << ' ' << sizeof(Real) << ' '
             << c.gamma << ' ' << c.width << ' ' << c.min_mach
             << ' ' << c.mach_max << ' ' << c.grad_floor << ' ' << c.div_limit
             << ' ' << c.density_floor << ' ' << c.temperature_floor << ' '
             << pin->GetOrAddString("time", "integrator", "rk2") << ' '
             << pin->GetOrAddString("hydro", "reconstruct", "plm") << ' '
             << pin->GetOrAddString("hydro", "rsolver", "llf") << ' '
             << pin->GetOrAddReal("time", "cfl_number", 0.3);
    auto &ms = pmy_pack->pmesh->mesh_size;
    auto &mi = pmy_pack->pmesh->mesh_indcs;
    auto &bi = pmy_pack->pmesh->mb_indcs;
    contract << ' ' << mi.nx1 << ' ' << mi.nx2 << ' ' << mi.nx3 << ' '
             << bi.nx1 << ' ' << bi.nx2 << ' ' << bi.nx3 << ' '
             << ms.x1min << ' ' << ms.x1max << ' ' << ms.x2min << ' ' << ms.x2max
             << ' ' << ms.x3min << ' ' << ms.x3max;
    for (auto bc : pmy_pack->pmesh->mesh_bcs) contract << ' ' << static_cast<int>(bc);
    if (pmy_pack->phydro) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      const char *names[] = {"dfloor","pfloor","tfloor","sfloor"};
      int index = 0;
      for (Real floor : {eos.dfloor,eos.pfloor,eos.tfloor,eos.sfloor}) {
        if (!std::isfinite(floor) || floor < 0) {
          std::cerr << "Shock MC requires finite nonnegative native EOS floors"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        pin->SetRealPrecise("hydro",names[index++],floor);
        contract << ' ' << floor;
      }
      contract << ' ' << pmy_pack->phydro->use_fofc
               << ' ' << pmy_pack->phydro->nscalars;
    }
    contract << " rng " << random_seed << " schedules " << seed_schedules_.size();
    for (const auto &seed : seed_schedules_) {
      // Native GetOrAddReal formats defaults with limited precision. Persist
      // resolved immutable values exactly so authentic restarts retain the contract.
      for (const auto &entry : std::initializer_list<std::pair<const char*,Real>>{
          {"start_time",seed.start_time},{"end_time",seed.end_time},
          {"cadence",seed.cadence},{"x1min",seed.x1min},{"x1max",seed.x1max},
          {"x2min",seed.x2min},{"x2max",seed.x2max},{"x3min",seed.x3min},
          {"x3max",seed.x3max},{"center1",seed.center1},{"center2",seed.center2},
          {"center3",seed.center3},{"radius",seed.radius},
          {"slab_min",seed.slab_min},{"slab_max",seed.slab_max}}) {
        if (!std::isfinite(entry.second)) {
          std::cerr << "Shock tracer schedule contains nonfinite values" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        pin->SetRealPrecise(seed.block_name,entry.first,entry.second);
      }
      if ((seed.has_target_min && !std::isfinite(seed.target_min)) ||
          (seed.has_target_max && !std::isfinite(seed.target_max))) {
        std::cerr << "Shock tracer target contains nonfinite thresholds" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (seed.has_target_min) {
        pin->SetRealPrecise(seed.block_name,"target_min",seed.target_min);
      }
      if (seed.has_target_max) {
        pin->SetRealPrecise(seed.block_name,"target_max",seed.target_max);
      }
      contract << ' ' << seed.id << ' ' << seed.start_time << ' ' << seed.end_time
               << ' ' << seed.cadence << ' ' << seed.count_per_event << ' ' << seed.seed
               << ' ' << static_cast<int>(seed.weight)
               << ' ' << static_cast<int>(seed.region)
               << ' ' << seed.x1min << ' ' << seed.x1max
               << ' ' << seed.x2min << ' ' << seed.x2max
               << ' ' << seed.x3min << ' ' << seed.x3max
               << ' ' << seed.center1 << ' ' << seed.center2 << ' ' << seed.center3
               << ' ' << seed.radius << ' ' << seed.slab_axis << ' ' << seed.slab_min
               << ' ' << seed.slab_max << ' ' << seed.has_target
               << ' ' << seed.target_field.name << ' ' << seed.has_target_min
               << ' ' << seed.has_target_max << ' ' << seed.target_min
               << ' ' << seed.target_max;
    }
    shock_contract = contract.str();
    if (shock_contract.size()+shock_lineage.size()+1 > 16384) {
      std::cerr << "Shock restart contract exceeds supported metadata size" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Kokkos::realloc(shock_state, nprtcl_thispack);
    Kokkos::realloc(shock_events, shock_capacity);
    const auto &bi2 = pmy_pack->pmesh->mb_indcs;
    Kokkos::realloc(shock_signal, pmy_pack->nmb_thispack,
                    bi2.nx3 > 1 ? bi2.nx3 + 2*bi2.ng : 1,
                    bi2.nx2 + 2*bi2.ng, bi2.nx1 + 2*bi2.ng);
    pin->SetBoolean("shock_tracking", "history_present", true);
  }
  if (shock_enabled) ValidateShockPhysics();

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
  if (particle_type == ParticleType::lagrangian_mc) return;

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
