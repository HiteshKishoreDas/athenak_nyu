//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file prtcl_thermo_history.cpp
//! \brief shared thermodynamic history output for fluid particles and material curves

#include <sys/stat.h>

#include <algorithm>
#include <cstdint>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "globals.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"
#include "material_curves/material_curves.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "particles/particles.hpp"
#include "particles/tracer_fields.hpp"
#include "shocks/store.hpp"
#include "outputs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

struct FileHeader {
  char magic[32];
  int version;
  int real_size;
  int nfields;
  int names_size;
};



std::string JoinFieldNames(const std::vector<std::string> &names) {
  std::string joined;
  for (std::size_t n=0; n<names.size(); ++n) {
    if (n > 0) joined += "\n";
    joined += names[n];
  }
  return joined;
}

void FatalHistoryOutput(const std::string &msg) {
  std::cout << "### FATAL ERROR in prtcl_thermo_history.cpp" << std::endl
            << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

std::string ThermoFileHeader(const std::vector<std::string> &field_names) {
  const auto names=JoinFieldNames(field_names);
  FileHeader header{};
  std::strncpy(header.magic,"ATHK_PRTCL_THERMO_HISTORY",sizeof(header.magic)-1);
  header.version=2; header.real_size=sizeof(Real);
  header.nfields=field_names.size(); header.names_size=names.size();
  return std::string(reinterpret_cast<const char *>(&header),sizeof(header))+names;
}

std::string ThermoBlockHeader(int records, int nfields, int cycle, Real time) {
  std::string data(16,'\0');
  std::memcpy(data.data(),"ATHKTHPBLK",10);
  const int header[5]={2,records,4,4+nfields,cycle};
  data.append(reinterpret_cast<const char *>(header),sizeof(header));
  data.append(reinterpret_cast<const char *>(&time),sizeof(time));
  return data;
}

std::string ThermoBlock(const std::vector<int> &integers, const std::vector<Real> &reals,
                        int cycle, Real time, int nfields) {
  const std::size_t records=integers.size()/4, width=4+nfields;
  if (integers.size()%4 || records>INT_MAX || reals.size()!=records*width) {
    FatalHistoryOutput("invalid THP record extent");
  }
  auto data=ThermoBlockHeader(static_cast<int>(records),nfields,cycle,time);
  if (!integers.empty()) {
    data.append(reinterpret_cast<const char *>(integers.data()),integers.size()*sizeof(int));
    data.append(reinterpret_cast<const char *>(reals.data()),reals.size()*sizeof(Real));
  }
  return data;
}

void EnsureFileHeader(FILE *file, const std::vector<std::string> &field_names) {
  const std::string names_blob = JoinFieldNames(field_names);
  std::fseek(file, 0, SEEK_END);
  std::int64_t file_size = std::ftell(file);
  if (file_size == 0) {
    const auto data=ThermoFileHeader(field_names);
    if (std::fwrite(data.data(),1,data.size(),file)!=data.size()) {
      FatalHistoryOutput("failed to write particle thermo file header");
    }
    return;
  }

  std::fseek(file, 0, SEEK_SET);
  FileHeader header;
  if (std::fread(&header, sizeof(header), 1, file) != 1) {
    FatalHistoryOutput("existing particle thermo file has a truncated header");
  }
  if (std::strncmp(header.magic, "ATHK_PRTCL_THERMO_HISTORY", 25) != 0) {
    FatalHistoryOutput("existing particle thermo file has an unrecognized header");
  }
  if (header.version != 2) {
    FatalHistoryOutput("existing particle thermo file uses an older schema; "
                       "write to a clean output directory or remove the old file");
  }
  if (header.real_size != static_cast<int>(sizeof(Real)) ||
      header.nfields != static_cast<int>(field_names.size()) ||
      header.names_size != static_cast<int>(names_blob.size())) {
    FatalHistoryOutput("existing particle thermo file schema does not match this output");
  }
  std::string existing_names(header.names_size, '\0');
  if (header.names_size > 0 &&
      std::fread(existing_names.data(), sizeof(char), header.names_size, file) !=
      static_cast<std::size_t>(header.names_size)) {
    FatalHistoryOutput("existing particle thermo file has a truncated schema");
  }
  if (existing_names != names_blob) {
    FatalHistoryOutput("existing particle thermo file variable list does not match this "
                       "output");
  }
  std::fseek(file, 0, SEEK_END);
}


void AdvanceTracerClock(OutputParameters &op, Mesh *pm, ParameterInput *pin) {
  const bool new_cycle=pm->ncycle !=
      pin->GetOrAddInteger(op.block_name,"last_cycle",-1);
  op.file_number++;
  if (op.last_time < 0.0) {
    op.last_time = pm->time;
  } else if (new_cycle && (op.dt == 0.0 ||
             static_cast<float>(pm->time) >=
             static_cast<float>(op.last_time+op.dt))) {
    // Finalize may repeat a scheduled output or sample between cadence thresholds.
    // Advance only a reached threshold, or restart would silently skip the next one.
    op.last_time += op.dt;
  }
  pin->SetInteger(op.block_name, "file_number", op.file_number);
  pin->SetReal(op.block_name, "last_time", op.last_time);
  pin->SetInteger(op.block_name, "last_cycle", pm->ncycle);
}

} // namespace

//----------------------------------------------------------------------------------------
// ctor

ParticleThermoHistoryOutput::ParticleThermoHistoryOutput(ParameterInput *pin, Mesh *pm,
                                                         OutputParameters op) :
    BaseTypeOutput(pin, pm, op),
    npout_thisrank(0),
    npout_total(0),
    tracer_gamma(5.0/3.0) {
  const auto source=pin->GetOrAddString(op.block_name,"source","particles");
  if (source!="particles" && source!="curves" && source!="both") {
    FatalHistoryOutput("source must be particles, curves, or both");
  }
  sample_particles=(source!="curves");
  sample_curves=(source!="particles");
  if (sample_particles && (pm->pmb_pack->ppart == nullptr ||
                          !pm->pmb_pack->ppart->IsFluidTracer())) {
    FatalHistoryOutput("particle thermo source requires particle_type=lagrangian_mc "
                       "or lagrangian_tracer");
  }
  if (sample_curves && !pm->pmb_pack->pcurves) {
    FatalHistoryOutput("curve thermo source requires material_curves/enabled=true");
  }
  if (sample_curves && !pm->pmb_pack->phydro->peos->eos_data.is_ideal) {
    FatalHistoryOutput("curve thermo source requires ideal-gas hydro");
  }
  for (const auto &entry : {std::make_pair(sample_particles,"prtcl_thermo_history"),
                            std::make_pair(sample_curves,"curve_thermo_history")}) {
    if (entry.first) {
      std::error_code error;
      std::filesystem::create_directories(entry.second,error);
      if (error) FatalHistoryOutput("cannot create thermo directory: "+error.message());
    }
  }

  bool has_mhd = (pm->pmb_pack->pmhd != nullptr);
  tracer_gamma = has_mhd ? pin->GetOrAddReal("mhd", "gamma", 5.0/3.0) :
                           pin->GetOrAddReal("hydro", "gamma", 5.0/3.0);
  int nscalars = has_mhd ? pm->pmb_pack->pmhd->nscalars :
                          pm->pmb_pack->phydro->nscalars;
  std::string field_list = pin->DoesParameterExist(op.block_name, "variables") ?
                           pin->GetString(op.block_name, "variables") :
                           pin->GetOrAddString(sample_particles ? "particles" :
                             "material_curves", "track_variables", "default");
  tracer_fields = particles::ParseTracerFieldList(field_list, has_mhd, nscalars,
                                                  op.block_name + "/variables");
  tracer_field_names = particles::TracerFieldNames(tracer_fields);
  for (const auto &field : tracer_fields) sample_shocks |= particles::IsShockField(field.kind);
  if (sample_shocks && ((sample_particles && !pm->pmb_pack->ppart->shock_enabled) ||
                       (sample_curves && !pm->pmb_pack->pcurves->shock_enabled))) {
    FatalHistoryOutput("shock fields require enabled shock tracking for each selected source");
  }
}

//----------------------------------------------------------------------------------------
// LoadOutputData

void ParticleThermoHistoryOutput::LoadOutputData(Mesh *pm) {
  if (sample_curves) {
    pm->pmb_pack->pcurves->ThermoHistory(tracer_fields,curve_thermo_integers,
                                        curve_thermo_reals);
  }
  if (!sample_particles) return;
  particles::Particles *pp = pm->pmb_pack->ppart;
  npout_thisrank = pm->nprtcl_thisrank;
  npout_total = pm->nprtcl_total;
  Kokkos::realloc(outpart_rdata, pp->nrdata, npout_thisrank);
  Kokkos::realloc(outpart_idata, pp->nidata, npout_thisrank);
  Kokkos::realloc(outfield_data, tracer_fields.size(), npout_thisrank);
  if (npout_thisrank <= 0) return;

  Kokkos::deep_copy(outpart_rdata, pp->prtcl_rdata);
  Kokkos::deep_copy(outpart_idata, pp->prtcl_idata);

  auto &indcs = pm->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int ncells1 = indcs.nx1 + 2*indcs.ng;
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb_alloc = std::max(pm->pmb_pack->nmb_thispack, pm->nmb_maxperrank);
  int nfluid = 0, nscalars = 0;
  bool has_mhd = (pm->pmb_pack->pmhd != nullptr);

  if (pm->pmb_pack->phydro != nullptr) {
    nfluid = pm->pmb_pack->phydro->nhydro;
    nscalars = pm->pmb_pack->phydro->nscalars;
  } else {
    nfluid = pm->pmb_pack->pmhd->nmhd;
    nscalars = pm->pmb_pack->pmhd->nscalars;
  }

  HostArray5D<Real> h_w0("prtcl_history_w0", nmb_alloc, nfluid+nscalars,
                         ncells3, ncells2, ncells1);
  HostArray5D<Real> h_bcc("prtcl_history_bcc", nmb_alloc, NMAG,
                          ncells3, ncells2, ncells1);
  if (pm->pmb_pack->phydro != nullptr) {
    Kokkos::deep_copy(h_w0, pm->pmb_pack->phydro->w0);
  } else {
    Kokkos::deep_copy(h_w0, pm->pmb_pack->pmhd->w0);
    Kokkos::deep_copy(h_bcc, pm->pmb_pack->pmhd->bcc0);
  }

  pm->pmb_pack->pmb->mb_size.template sync<HostMemSpace>();
  auto h_size = pm->pmb_pack->pmb->mb_size.h_view;
  int nmb = pm->pmb_pack->nmb_thispack;
  int gids = pm->pmb_pack->gids;

  HostArray4D<shocks::Signal> signals;
  if (sample_shocks) {
    if (pp->shock_generation!=static_cast<std::uint64_t>(pm->ncycle) ||
        pp->shock_time!=pm->time) FatalHistoryOutput("stale particle shock output generation");
    signals=Kokkos::create_mirror_view_and_copy(HostMemSpace(),pp->shock_signal);
  }
  for (int p=0; p<npout_thisrank; ++p) {
    int m = outpart_idata(PGID,p) - gids;
    if (m < 0 || m >= nmb) {
      FatalHistoryOutput("tracer particle references a MeshBlock not owned by this rank");
    }
    int i = static_cast<int>((outpart_rdata(LMCX,p) - h_size(m).x1min)/
                             h_size(m).dx1) + is;
    int j = pm->multi_d ?
      static_cast<int>((outpart_rdata(LMCY,p) - h_size(m).x2min)/h_size(m).dx2) + js : js;
    int k = pm->three_d ?
      static_cast<int>((outpart_rdata(LMCZ,p) - h_size(m).x3min)/h_size(m).dx3) + ks : ks;
    i = std::max(is, std::min(is + indcs.nx1 - 1, i));
    j = std::max(js, std::min(js + indcs.nx2 - 1, j));
    k = std::max(ks, std::min(ks + indcs.nx3 - 1, k));
    for (int n=0; n<static_cast<int>(tracer_fields.size()); ++n) {
      outfield_data(n,p) = particles::IsShockField(tracer_fields[n].kind) ?
          particles::EvaluateShockField(tracer_fields[n].kind,signals(m,k,j,i)) :
          particles::EvaluateTracerFieldHost(tracer_fields[n], h_w0,
                                                              h_bcc, has_mhd,
                                                              tracer_gamma,
                                                              nfluid, m, k, j, i);
    }
  }
}

//----------------------------------------------------------------------------------------
// WriteOutputFile

void ParticleThermoHistoryOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  if (sample_curves) {
    const auto path="curve_thermo_history/"+out_params.file_basename+"."+
        out_params.file_id+".cycle"+std::to_string(pm->ncycle)+shocks::RankName()+
        "of"+std::to_string(global_variable::nranks)+".bits"+
        std::to_string(pm->pmb_pack->pcurves->curve_id_bits)+".thp";
    std::vector<std::string> names={"birth","a","level"};
    names.insert(names.end(),tracer_field_names.begin(),tracer_field_names.end());
    const auto data=ThermoFileHeader(names)+ThermoBlock(curve_thermo_integers,
        curve_thermo_reals,pm->ncycle,pm->time,names.size());
    shocks::CommitFile(path,data,true);
  }
  if (sample_particles) WriteParticleData(pm);

  AdvanceTracerClock(out_params,pm,pin);
}

void ParticleThermoHistoryOutput::WriteParticleData(Mesh *pm) {
  constexpr int int_per_record = 4;   // cycle, tag, seed_id, gid
  int real_per_record = 4 + static_cast<int>(tracer_fields.size()); // time, x, fields

  std::vector<int> local_i(int_per_record*npout_thisrank);
  std::vector<Real> local_r(real_per_record*npout_thisrank);
  for (int p=0; p<npout_thisrank; ++p) {
    local_i[int_per_record*p] = pm->ncycle;
    local_i[int_per_record*p + 1] = outpart_idata(PTAG,p);
    local_i[int_per_record*p + 2] = outpart_idata(PSEEDID,p);
    local_i[int_per_record*p + 3] = outpart_idata(PGID,p);

    local_r[real_per_record*p] = pm->time;
    local_r[real_per_record*p + 1] = outpart_rdata(LMCX,p);
    local_r[real_per_record*p + 2] = outpart_rdata(LMCY,p);
    local_r[real_per_record*p + 3] = outpart_rdata(LMCZ,p);
    for (int n=0; n<static_cast<int>(tracer_fields.size()); ++n) {
      local_r[real_per_record*p + 4 + n] = outfield_data(n,p);
    }
  }

  std::vector<int> all_i;
  std::vector<Real> all_r;
#if MPI_PARALLEL_ENABLED
  std::vector<int> counts(global_variable::nranks), idispl(global_variable::nranks),
                   rcounts(global_variable::nranks), rdispl(global_variable::nranks);
  int local_count = npout_thisrank;
  MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (global_variable::my_rank == 0) {
    int total_i = 0, total_r = 0;
    for (int n=0; n<global_variable::nranks; ++n) {
      idispl[n] = total_i;
      rdispl[n] = total_r;
      total_i += counts[n]*int_per_record;
      total_r += counts[n]*real_per_record;
      rcounts[n] = counts[n]*real_per_record;
      counts[n] *= int_per_record;
    }
    all_i.resize(total_i);
    all_r.resize(total_r);
  }
  int local_i_count = static_cast<int>(local_i.size());
  int local_r_count = static_cast<int>(local_r.size());
  MPI_Gatherv(local_i.data(), local_i_count, MPI_INT,
              all_i.data(), counts.data(), idispl.data(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local_r.data(), local_r_count, MPI_ATHENA_REAL,
              all_r.data(), rcounts.data(), rdispl.data(), MPI_ATHENA_REAL, 0,
              MPI_COMM_WORLD);
#else
  all_i = std::move(local_i);
  all_r = std::move(local_r);
#endif

  if (global_variable::my_rank == 0) {
    std::string fname = "prtcl_thermo_history/" + out_params.file_basename + "." +
                        out_params.file_id + ".thp";
    FILE *file = std::fopen(fname.c_str(), "ab+");
    if (file == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Output file '" << fname << "' could not be opened"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    EnsureFileHeader(file, tracer_field_names);

    const auto header=ThermoBlockHeader(npout_total,tracer_fields.size(),pm->ncycle,pm->time);
    if (std::fwrite(header.data(),1,header.size(),file)!=header.size() ||
        (!all_i.empty() && std::fwrite(all_i.data(),sizeof(int),all_i.size(),file)!=all_i.size()) ||
        (!all_r.empty() &&
         std::fwrite(all_r.data(),sizeof(Real),all_r.size(),file)!=all_r.size())) {
      FatalHistoryOutput("failed to write particle thermo history block");
    }
    std::fclose(file);
  }

}

// Geometry uses the native output scheduler, basename/id and output numbering.
MaterialCurveOutput::MaterialCurveOutput(ParameterInput *pin, Mesh *pm,
                                         OutputParameters op) : BaseTypeOutput(pin,pm,op) {
  if (!pm->pmb_pack->pcurves) FatalHistoryOutput("file_type=curves requires material curves");
  std::error_code error;
  std::filesystem::create_directories("curves",error);
  if (error) FatalHistoryOutput("cannot create curves directory: "+error.message());
}
void MaterialCurveOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  std::ostringstream path;
  path << "curves/" << out_params.file_basename << '.' << out_params.file_id << '.'
       << std::setw(8) << std::setfill('0') << out_params.file_number << ".txt";
  pm->pmb_pack->pcurves->Snapshot(path.str());
  AdvanceTracerClock(out_params,pm,pin);
}
