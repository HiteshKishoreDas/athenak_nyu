// Copyright (C) 2026 the AthenaK contributors. BSD-3-Clause.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "material_curves/material_curves.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "srcterms/srcterms.hpp"
#include "srcterms/turb_driver.hpp"
#include "shock_version.hpp"  // NOLINT(build/include_subdir): generated header

namespace material_curves {
void MaterialCurves::ConfigureShocks(ParameterInput *pin) {
  if (!pin->DoesBlockExist("shock_tracking")) return;
  shock_enabled=pin->GetOrAddBoolean("shock_tracking","enabled",false);
  if (!shock_enabled) {
    if (pin->GetOrAddBoolean("shock_tracking","curve_history_present",false)) {
      Fail("curve shock history cannot be disabled on restart");
    }
    return;
  }
  auto hydro=mesh->pmb_pack->phydro;
  if (std::string(DevExeSpace::name())!="Serial" || sizeof(Real)!=8 ||
      !hydro->peos->eos_data.is_ideal || !hydro->peos->eos_data.use_e ||
      hydro->use_fofc || hydro->pvisc || hydro->pcond) {
    Fail("curve shock histories require Serial/double ideal hydro "
         "without FOFC/transport");
  }
  if (mesh->pmb_pack->pturb &&
      !mesh->pmb_pack->pturb->RestartMetadata().StaticFieldsFinite()) {
    Fail("curve shock histories require finite native forcing parameters");
  }
  auto &c=shock_config;
  c.gamma=hydro->peos->eos_data.gamma;
  c.width=pin->GetOrAddInteger("shock_tracking","width",1);
  c.min_mach=pin->GetOrAddReal("shock_tracking","min_mach",1.1);
  c.mach_max=pin->GetOrAddReal("shock_tracking","mach_max",1000);
  c.grad_floor=pin->GetOrAddReal("shock_tracking","grad_floor",1e-30);
  c.div_limit=pin->GetOrAddReal("shock_tracking","div_limit",0);
  c.density_floor=pin->GetOrAddReal("shock_tracking","density_floor",0);
  c.temperature_floor=pin->GetOrAddReal("shock_tracking","temperature_floor",0);
  shock_capacity=pin->GetOrAddInteger("shock_tracking","event_capacity",4096);
  shock_trace=pin->GetOrAddBoolean("shock_tracking","trace",false);
  shock_lineage=pin->GetString("shock_tracking","lineage");
  shock_directory=
      pin->GetOrAddString("shock_tracking","directory","shock_events")+"/curves";
  bool valid=!shock_lineage.empty() && shock_lineage.size()<=128 && shock_capacity>0 &&
             shocks::ConfigValid(c) && c.div_limit<=0 && c.width<=mesh->mb_indcs.ng;
  for (char ch : shock_lineage) {
    valid=valid && ((ch>='a' && ch<='z') || (ch>='A' && ch<='Z') ||
                    (ch>='0' && ch<='9') || ch=='-' || ch=='_');
  }
  if (!valid) Fail("invalid curve shock detector, lineage, capacity or ghost width");
  std::ostringstream contract;
  contract << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "curve-node-history-v1 " << SHOCK_IMPLEMENTATION_HASH << ' '
           << sizeof(ShockNode) << ' ' << c.gamma << ' ' << c.width << ' '
           << c.min_mach << ' ' << c.mach_max << ' ' << c.grad_floor << ' '
           << c.div_limit << ' ' << c.density_floor << ' ' << c.temperature_floor << ' '
           << pin->GetString("time","integrator") << ' '
           << pin->GetOrAddReal("time","cfl_number",0.3) << ' '
           << pin->GetOrAddString("hydro","reconstruct","plm") << ' '
           << pin->GetOrAddString("hydro","rsolver","llf");
  pin->SetRealPrecise("hydro","gamma",c.gamma);
  const auto eos=hydro->peos->eos_data;
  const char *floor_names[]={"dfloor","pfloor","tfloor","sfloor"};
  int index=0;
  for (Real floor : {eos.dfloor,eos.pfloor,eos.tfloor,eos.sfloor}) {
    if (!std::isfinite(floor) || floor<0) {
      Fail("invalid native EOS floor for curve history");
    }
    pin->SetRealPrecise("hydro",floor_names[index++],floor);
    contract << ' ' << floor;
  }
  contract << ' ' << hydro->nscalars;
  for (const char *block : {"hydro","turb_driving","mesh_refinement","problem"}) {
    std::vector<std::pair<std::string,std::string>> settings;
    for (const auto &section : pin->block) {
      if (section.block_name!=block) continue;
      for (const auto &line : section.line) {
        settings.emplace_back(line.param_name,line.param_value);
      }
    }
    std::sort(settings.begin(),settings.end());
    contract << ' ' << block << ' ' << settings.size();
    for (const auto &entry : settings) {
      contract << ' ' << entry.first.size() << ':' << entry.first << ' '
               << entry.second.size() << ':' << entry.second;
    }
  }
  shock_contract=contract.str();
  if (shock_contract.size()+shock_lineage.size()+1>16384) {
    Fail("curve shock contract exceeds bounded metadata");
  }
  pin->SetBoolean("shock_tracking","curve_history_present",true);
}

void MaterialCurves::ObserveShocks(bool births_only) {
  if (!shock_enabled) return;
  auto generation=static_cast<std::uint64_t>(mesh->ncycle);
  if (!births_only && generation==shock_generation) return;
  if (!births_only && shock_generation!=UINT64_MAX && generation!=shock_generation+1) {
    Fail("curve shock observation generation gap");
  }
  std::vector<Node *> nodes;
  for (auto &curve : curves) for (auto &node : curve.nodes) {
    if (!births_only || shock_states.find(node.id)==shock_states.end()) {
      nodes.push_back(&node);
    }
  }
  std::vector<ShockObservation> samples;
  Sample(nodes,&samples);  // Collective even on ranks with no locally owned nodes.
  std::vector<ShockEvent> events;
  const Real low[3]={mesh->mesh_size.x1min,mesh->mesh_size.x2min,mesh->mesh_size.x3min};
  const Real high[3]={mesh->mesh_size.x1max,mesh->mesh_size.x2max,mesh->mesh_size.x3max};
  for (std::size_t p=0; p<nodes.size(); ++p) {
    const auto &node=*nodes[p];
    auto &entry=shock_states[node.id];
    auto &state=entry.encounter;
    auto &o=samples[p].observation;
    if (state.last_generation==UINT64_MAX) {
      state.particle_id=node.id;
      entry.birth=node.birth;
    }
    if (o.time<entry.birth || state.sequence==UINT64_MAX ||
        state.suppressed_reentries==UINT64_MAX) {
      Fail("invalid curve shock birth/sequence");
    }
    for (int d=0; d<3; ++d) {
      long double image=std::floor((static_cast<long double>(node.x[d])-low[d])/
                                  (high[d]-low[d]));
      if (!std::isfinite(image) || image<static_cast<long double>(INT64_MIN) ||
          image>=-static_cast<long double>(INT64_MIN)) {
        Fail("curve shock periodic image exhausted");
      }
      o.image[d]=static_cast<std::int64_t>(image);
      state.image[d]=o.image[d];
    }
    auto phase=state.phase;
    shocks::Event event;
    if (o.valid && o.inside && !std::isfinite(o.mach_temperature)) {
      state.flags |= shocks::temperature_invalid;
    }
    bool emit=shocks::UpdateEncounter(state,o,event);
    if (state.phase==shocks::active && o.inside) {
      if (phase==shocks::outside) {
        entry.first_level=samples[p].level;
        for (int d=0; d<3; ++d) entry.first_spacing[d]=samples[p].spacing[d];
      }
      entry.last_level=samples[p].level;
      for (int d=0; d<3; ++d) entry.last_spacing[d]=samples[p].spacing[d];
    }
    if (emit) {
      if (events.size()>=static_cast<std::size_t>(shock_capacity)) {
        Fail("curve shock event buffer overflow; "
             "replay prior checkpoint with larger capacity");
      }
      ShockEvent out;
      out.encounter=event; out.birth=entry.birth;
      out.first_level=entry.first_level; out.last_level=entry.last_level;
      for (int d=0; d<3; ++d) {
        out.first_spacing[d]=entry.first_spacing[d];
        out.last_spacing[d]=entry.last_spacing[d];
      }
      events.push_back(out);
    }
  }
  if (shock_trace && (!births_only || !nodes.empty())) {
    WriteShockTrace(nodes,samples,births_only);
  }
  if (!births_only) {
    shock_generation=generation;
    shock_time=mesh->time;
    WriteShockChunk(events);
  }
}
}  // namespace material_curves
