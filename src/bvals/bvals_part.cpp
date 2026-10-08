//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_part.cpp
//! \brief

#include <climits>
#include <cstdint>
#include <type_traits>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/nghbr_index.hpp"
#include "mesh/mesh.hpp"
#include "particles/particles.hpp"
#include "bvals.hpp"

namespace particles {
static_assert(std::is_trivially_copyable<shocks::EncounterState>::value,
              "Encounter migration requires trivially copyable typed state");

#if MPI_PARALLEL_ENABLED
// Native bytes require the same executable and ABI on every rank.
static int ShockMessageBytes(int count) {
  if (count < 0 || static_cast<size_t>(count) >
      static_cast<size_t>(INT_MAX)/sizeof(shocks::EncounterState)) {
    std::cerr << "FATAL: shock migration MPI byte count exceeds INT_MAX" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return static_cast<int>(count*sizeof(shocks::EncounterState));
}

static void CheckShockBuffers(int count, int nrdata, int nidata) {
  ShockMessageBytes(count);
  if (nrdata <= 0 || nidata <= 0 || count > INT_MAX/nrdata ||
      count > INT_MAX/nidata) {
    std::cerr << "FATAL: shock migration data count exceeds INT_MAX" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
#endif

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::UpdateGID()
//! \brief Updates GID of particles that cross boundary of their parent MeshBlock.  If
//! the new GID is on a different rank, then store in sendlist_buf DvceArray: (1) index of
//! particle in prtcl array, (2) destination GID, and (3) destination rank.

KOKKOS_INLINE_FUNCTION
void UpdateGID(int &newgid, NeighborBlock nghbr, int myrank,
               DvceArray1D<int> send_count,
               DvceArray1D<ParticleLocationData> slist, int p) {
  newgid = nghbr.gid;
#if MPI_PARALLEL_ENABLED
  if (nghbr.rank != myrank) {
    int index = Kokkos::atomic_fetch_add(&send_count(0),1);
    slist(index).prtcl_indx = p;
    slist(index).dest_gid   = nghbr.gid;
    slist(index).dest_rank  = nghbr.rank;
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::SetNewGID()
//! \brief

TaskStatus ParticlesBoundaryValues::SetNewPrtclGID() {
  if (pmy_part->shock_enabled) return SetShockPrtclGID();
  // create local references for variables in kernel
  auto gids = pmy_part->pmy_pack->gids;
  auto &pr = pmy_part->prtcl_rdata;
  auto &pi = pmy_part->prtcl_idata;
  int npart = pmy_part->nprtcl_thispack;
  auto &mbsize = pmy_part->pmy_pack->pmb->mb_size;
  auto &mblev = pmy_part->pmy_pack->pmb->mb_lev;
  auto &meshsize = pmy_part->pmy_pack->pmesh->mesh_size;
  auto myrank = global_variable::my_rank;
  auto &nghbr = pmy_part->pmy_pack->pmb->nghbr;
  bool &multi_d = pmy_part->pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_part->pmy_pack->pmesh->three_d;

  Kokkos::realloc(sendlist, npart);
  send_count.h_view(0) = 0;
  send_count.template modify<HostMemSpace>();
  send_count.template sync<DevExeSpace>();
  auto d_send_count = send_count.d_view;
  auto d_sendlist = sendlist.d_view;
  par_for("part_update",DevExeSpace(),0,(npart-1), KOKKOS_LAMBDA(const int p) {
    int m = pi(PGID,p) - gids;
    int mylevel = mblev.d_view(m);
    Real x1 = pr(IPX,p);
    Real x2 = pr(IPY,p);
    Real x3 = pr(IPZ,p);

    // length of MeshBlock in each direction
    Real lx = (mbsize.d_view(m).x1max - mbsize.d_view(m).x1min);
    Real ly = (mbsize.d_view(m).x2max - mbsize.d_view(m).x2min);
    Real lz = (mbsize.d_view(m).x3max - mbsize.d_view(m).x3min);

    // integer offset of particle relative to center of MeshBlock (-1,0,+1)
    int ix = static_cast<int>((x1 - mbsize.d_view(m).x1min + lx)/lx) - 1;
    int iy = static_cast<int>((x2 - mbsize.d_view(m).x2min + ly)/ly) - 1;
    int iz = static_cast<int>((x3 - mbsize.d_view(m).x3min + lz)/lz) - 1;

    // sublock indices for faces and edges with S/AMR
    int fx = (x1 < 0.5*(mbsize.d_view(m).x1min + mbsize.d_view(m).x1max))? 0 : 1;
    int fy = (x2 < 0.5*(mbsize.d_view(m).x2min + mbsize.d_view(m).x2max))? 0 : 1;
    int fz = (x3 < 0.5*(mbsize.d_view(m).x3min + mbsize.d_view(m).x3max))? 0 : 1;
    fy = multi_d ? fy : 0;
    fz = three_d ? fz : 0;

    // only update particle GID if it has crossed MeshBlock boundary
    if ((abs(ix) + abs(iy) + abs(iz)) != 0) {
      if (iz == 0) {
        if (iy == 0) {
          // x1 face
          int indx = NeighborIndex(ix,0,0,0,0);           // neighbor at same level
          if (nghbr.d_view(m,indx).lev > mylevel) {       // neighbor at finer level
            indx = NeighborIndex(ix,0,0,fy,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}  // neighbor at coarser level
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        } else if (ix == 0) {
          // x2 face
          int indx = NeighborIndex(0,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,0,fx,fz);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        } else {
          // x1x2 edge
          int indx = NeighborIndex(ix,iy,0,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,iy,0,fz,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        }
      } else if (iy == 0) {
        if (ix == 0) {
          // x3 face
          int indx = NeighborIndex(0,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,0,iz,fx,fy);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        } else {
          // x3x1 edge
          int indx = NeighborIndex(ix,0,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(ix,0,iz,fy,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        }
      } else {
        if (ix == 0) {
          // x2x3 edge
          int indx = NeighborIndex(0,iy,iz,0,0);
          if (nghbr.d_view(m,indx).lev > mylevel) {
            indx = NeighborIndex(0,iy,iz,fx,0);
          }
          while (nghbr.d_view(m,indx).gid < 0) {indx++;}
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        } else {
          // corners
          int indx = NeighborIndex(ix,iy,iz,0,0);
          UpdateGID(pi(PGID,p), nghbr.d_view(m,indx), myrank,
                    d_send_count, d_sendlist, p);
        }
      }

      // reset x,y,z positions if particle crosses Mesh boundary using periodic BCs
      if (x1 < meshsize.x1min) {
        pr(IPX,p) += (meshsize.x1max - meshsize.x1min);
      } else if (x1 > meshsize.x1max) {
        pr(IPX,p) -= (meshsize.x1max - meshsize.x1min);
      }
      if (x2 < meshsize.x2min) {
        pr(IPY,p) += (meshsize.x2max - meshsize.x2min);
      } else if (x2 > meshsize.x2max) {
        pr(IPY,p) -= (meshsize.x2max - meshsize.x2min);
      }
      if (x3 < meshsize.x3min) {
        pr(IPZ,p) += (meshsize.x3max - meshsize.x3min);
      } else if (x3 > meshsize.x3max) {
        pr(IPZ,p) -= (meshsize.x3max - meshsize.x3min);
      }
    }
  });
  send_count.template modify<DevExeSpace>();
  send_count.template sync<HostMemSpace>();
  nprtcl_send = send_count.h_view(0);
  Kokkos::resize(sendlist, nprtcl_send);
  // sync sendlist device array with host
  sendlist.template modify<DevExeSpace>();
  sendlist.template sync<HostMemSpace>();

  return TaskStatus::complete;
}

// Uniform-grid shock path preserves the existing disabled particle implementation.
TaskStatus ParticlesBoundaryValues::SetShockPrtclGID() {
  const int npart = pmy_part->nprtcl_thispack;
  if (pmy_part->shock_state.extent(0) != static_cast<size_t>(npart)) {
    std::cerr << "FATAL: particle shock-state extent mismatch" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto *pack = pmy_part->pmy_pack;
  const int gids = pack->gids, nmb = pack->nmb_thispack;
  const bool multi_d = pack->pmesh->multi_d, three_d = pack->pmesh->three_d;
  const auto meshsize = pack->pmesh->mesh_size;
  Kokkos::Array<BoundaryFlag,6> bcs;
  for (int f=0; f<6; ++f) bcs[f] = pack->pmesh->mesh_bcs[f];
  auto pr = pmy_part->prtcl_rdata;
  auto pi = pmy_part->prtcl_idata;
  auto states = pmy_part->shock_state;
  auto sizes = pack->pmb->mb_size.d_view;
  auto levels = pack->pmb->mb_lev.d_view;
  auto neighbors = pack->pmb->nghbr.d_view;
  const int myrank = global_variable::my_rank;
  Kokkos::realloc(sendlist, npart);
  Kokkos::realloc(shock_destroylist, npart);
  Kokkos::realloc(shock_destroy_count, 1);
  send_count.h_view(0) = 0;
  send_count.template modify<HostMemSpace>();
  send_count.template sync<DevExeSpace>();
  shock_destroy_count.h_view(0) = 0;
  shock_destroy_count.template modify<HostMemSpace>();
  shock_destroy_count.template sync<DevExeSpace>();
  auto scount = send_count.d_view;
  auto dcount = shock_destroy_count.d_view;
  auto slist = sendlist.d_view;
  auto dlist = shock_destroylist.d_view;
  par_for("shock_particle_gid",DevExeSpace(),0,npart-1,
  KOKKOS_LAMBDA(const int p) {
    const int m = pi(PGID,p) - gids;
    if (m < 0 || m >= nmb) Kokkos::abort("invalid shock particle MeshBlock");
    Real x[3] = {pr(IPX,p), pr(IPY,p), pr(IPZ,p)};
    Real lower[3] = {meshsize.x1min, meshsize.x2min, meshsize.x3min};
    Real upper[3] = {meshsize.x1max, meshsize.x2max, meshsize.x3max};
    Real blo[3] = {sizes(m).x1min, sizes(m).x2min, sizes(m).x3min};
    Real bhi[3] = {sizes(m).x1max, sizes(m).x2max, sizes(m).x3max};
    const int ndim = three_d ? 3 : (multi_d ? 2 : 1);
    int offset[3] = {};
    bool destroy = false;
    for (int d=0; d<ndim; ++d) {
      if (!Kokkos::isfinite(x[d])) Kokkos::abort("nonfinite shock particle position");
      const Real length = bhi[d] - blo[d];
      if (x[d] < blo[d]-length || x[d] >= bhi[d]+length) {
        Kokkos::abort("shock particle moved beyond adjacent MeshBlock");
      }
      offset[d] = x[d] < blo[d] ? -1 : (x[d] >= bhi[d] ? 1 : 0);
      const int face = x[d] < lower[d] ? 2*d : (x[d] >= upper[d] ? 2*d+1 : -1);
      if (face >= 0 && bcs[face] == BoundaryFlag::outflow) destroy = true;
      if (face >= 0 && bcs[face] != BoundaryFlag::outflow &&
          bcs[face] != BoundaryFlag::periodic) {
        Kokkos::abort("unsupported shock particle physical boundary");
      }
    }
    // Outflow has no neighbor. Emit its evidence on the host before deletion.
    if (destroy) {
      int n = Kokkos::atomic_fetch_add(&dcount(0),1);
      dlist(n) = p;
      return;
    }
    if (offset[0] || offset[1] || offset[2]) {
      const int index = NeighborIndex(offset[0],offset[1],offset[2],0,0);
      const auto neighbor = neighbors(m,index);
      if (neighbor.gid < 0 || neighbor.lev != levels(m)) {
        Kokkos::abort("shock particle requires a uniform adjacent neighbor");
      }
      UpdateGID(pi(PGID,p), neighbor, myrank, scount, slist, p);
    }
    for (int d=0; d<ndim; ++d) {
      const bool low = x[d] < lower[d], high = x[d] >= upper[d];
      if (low || high) {
        if ((low && states(p).image[d] == INT64_MIN) ||
            (high && states(p).image[d] == INT64_MAX)) {
          Kokkos::abort("shock periodic image counter limit reached");
        }
        states(p).image[d] += low ? -1 : 1;
        pr(2*d,p) += low ? upper[d]-lower[d] : lower[d]-upper[d];
        if (pr(2*d,p) < lower[d] || pr(2*d,p) >= upper[d]) {
          Kokkos::abort("shock particle periodic wrap did not enter domain");
        }
      }
    }
  });
  send_count.template modify<DevExeSpace>();
  send_count.template sync<HostMemSpace>();
  shock_destroy_count.template modify<DevExeSpace>();
  shock_destroy_count.template sync<HostMemSpace>();
  nprtcl_send = send_count.h_view(0);
  nprtcl_recv = 0;
  nsends = nrecvs = 0;
  nprtcl_destroy = shock_destroy_count.h_view(0);
  Kokkos::resize(sendlist, nprtcl_send);
  sendlist.template modify<DevExeSpace>();
  sendlist.template sync<HostMemSpace>();
  Kokkos::resize(shock_destroylist, nprtcl_destroy);
  shock_destroylist.template modify<DevExeSpace>();
  shock_destroylist.template sync<HostMemSpace>();
  std::vector<int> exits;
  for (int n=0; n<nprtcl_destroy; ++n) exits.push_back(shock_destroylist.h_view(n));
  std::sort(exits.begin(), exits.end());
  if (!exits.empty()) pmy_part->EmitShockExits(exits);
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::CountSendsAndRecvs()
//! \brief

TaskStatus ParticlesBoundaryValues::CountSendsAndRecvs() {
#if MPI_PARALLEL_ENABLED
  // Sort sendlist on host by destrank.
  namespace KE = Kokkos::Experimental;
  std::sort(KE::begin(sendlist.h_view), KE::end(sendlist.h_view), SortByRank);
  // sync sendlist host array with device.  This results in sorted array on device
  sendlist.template modify<HostMemSpace>();
  sendlist.template sync<DevExeSpace>();

  // load STL::vector of ParticleMessageData with <sendrank, recvrank, nprtcls> for sends
  // from this rank. Length will be nsends; initially this length is unknown
  sends_thisrank.clear();
  if (nprtcl_send > 0) {
    int &myrank = global_variable::my_rank;
    int rank = sendlist.h_view(0).dest_rank;
    int nprtcl = 1;

    for (int n=1; n<nprtcl_send; ++n) {
      if (sendlist.h_view(n).dest_rank == rank) {
        ++nprtcl;
      } else {
        sends_thisrank.emplace_back(ParticleMessageData(myrank,rank,nprtcl));
        rank = sendlist.h_view(n).dest_rank;
        nprtcl = 1;
      }
    }
    sends_thisrank.emplace_back(ParticleMessageData(myrank,rank,nprtcl));
  }
  nsends = sends_thisrank.size();

  // Share number of ranks to send to amongst all ranks
  nsends_eachrank[global_variable::my_rank] = nsends;
  MPI_Allgather(&nsends, 1, MPI_INT, nsends_eachrank.data(), 1, MPI_INT, mpi_comm_part);

  // Now share ParticleMessageData amongst all ranks
  // First create vector of starting indices in full vector
  std::vector<int> nsends_displ;
  nsends_displ.resize(global_variable::nranks);
  nsends_displ[0] = 0;
  for (int n=1; n<(global_variable::nranks); ++n) {
    if (pmy_part->shock_enabled && (nsends_eachrank[n-1] < 0 ||
        nsends_displ[n-1] > INT_MAX-nsends_eachrank[n-1])) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    nsends_displ[n] = nsends_displ[n-1] + nsends_eachrank[n-1];
  }
  if (pmy_part->shock_enabled && (nsends_eachrank.back() < 0 ||
      nsends_displ.back() > INT_MAX-nsends_eachrank.back())) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  int nsends_allranks = nsends_displ[global_variable::nranks - 1] +
                        nsends_eachrank[global_variable::nranks - 1];
  // Load ParticleMessageData on this rank into full vector
  sends_allranks.resize(nsends_allranks, ParticleMessageData(0,0,0));
  for (int n=0; n<nsends_eachrank[global_variable::my_rank]; ++n) {
    sends_allranks[n + nsends_displ[global_variable::my_rank]] = sends_thisrank[n];
  }

  // Share tuples using MPI derived data type for tuple of 3*int
  MPI_Datatype mpi_ituple;
  MPI_Type_contiguous(3, MPI_INT, &mpi_ituple);
  MPI_Type_commit(&mpi_ituple);
  MPI_Allgatherv(MPI_IN_PLACE, nsends_eachrank[global_variable::my_rank],
                   mpi_ituple, sends_allranks.data(), nsends_eachrank.data(),
                   nsends_displ.data(), mpi_ituple, mpi_comm_part);
  MPI_Type_free(&mpi_ituple);
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::InitPrtclRecv()
//! \brief

TaskStatus ParticlesBoundaryValues::InitPrtclRecv() {
#if MPI_PARALLEL_ENABLED
  // load STL::vector of ParticleMessageData with <sendrank,recvrank,nprtcl_recv> for
  // receives // on this rank. Length will be nrecvs, initially this length is unknown
  recvs_thisrank.clear();

  int nsends_allranks = sends_allranks.size();
  for (int n=0; n<nsends_allranks; ++n) {
    if (sends_allranks[n].recvrank == global_variable::my_rank) {
      recvs_thisrank.emplace_back(sends_allranks[n]);
    }
  }
  nrecvs = recvs_thisrank.size();

  // Figure out how many particles will be received from all ranks
  nprtcl_recv=0;
  for (int n=0; n<nrecvs; ++n) {
    if (pmy_part->shock_enabled &&
        (recvs_thisrank[n].nprtcls < 0 ||
         nprtcl_recv > INT_MAX-recvs_thisrank[n].nprtcls)) {
      std::cerr << "FATAL: shock receive count exceeds INT_MAX" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    nprtcl_recv += recvs_thisrank[n].nprtcls;
  }

  if (pmy_part->shock_enabled) {
    CheckShockBuffers(nprtcl_recv, pmy_part->nrdata, pmy_part->nidata);
    Kokkos::realloc(shock_recvbuf, nprtcl_recv);
    shock_recv_req.assign(nrecvs, MPI_REQUEST_NULL);
  }

  // Allocate receive buffer
  Kokkos::realloc(prtcl_rrecvbuf, (pmy_part->nrdata)*nprtcl_recv);
  Kokkos::realloc(prtcl_irecvbuf, (pmy_part->nidata)*nprtcl_recv);

  // Post non-blocking receives
  bool no_errors=true;
  rrecv_req.clear();
  irecv_req.clear();
  for (int n=0; n<nrecvs; ++n) {
    rrecv_req.emplace_back(MPI_REQUEST_NULL);
    irecv_req.emplace_back(MPI_REQUEST_NULL);
  }

  // Init receives for Reals
  int data_start=0;
  for (int n=0; n<nrecvs; ++n) {
    // calculate amount of data to be passed, get pointer to variables
    int data_size = (pmy_part->nrdata)*(recvs_thisrank[n].nprtcls);
    int data_end = data_start + (pmy_part->nrdata)*(recvs_thisrank[n].nprtcls - 1);
    auto recv_ptr = Kokkos::subview(prtcl_rrecvbuf, std::make_pair(data_start, data_end));
    int drank = recvs_thisrank[n].sendrank;
    int tag = 0; // 0 for Reals, 1 for ints

    // Post non-blocking receive
    int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                         mpi_comm_part, &(rrecv_req[n]));
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    data_start += data_size;
  }
  // Init receives for ints
  data_start=0;
  for (int n=0; n<nrecvs; ++n) {
    // calculate amount of data to be passed, get pointer to variables
    int data_size = (pmy_part->nidata)*(recvs_thisrank[n].nprtcls);
    int data_end = data_start + (pmy_part->nidata)*(recvs_thisrank[n].nprtcls - 1);
    auto recv_ptr = Kokkos::subview(prtcl_irecvbuf, std::make_pair(data_start, data_end));
    int drank = recvs_thisrank[n].sendrank;
    int tag = 1; // 0 for Reals, 1 for ints

    // Post non-blocking receive
    int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_INT, drank, tag,
                         mpi_comm_part, &(irecv_req[n]));
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    data_start += data_size;
  }

  if (pmy_part->shock_enabled) {
    int particle_start = 0;
    for (int n=0; n<nrecvs; ++n) {
      const int count = recvs_thisrank[n].nprtcls;
      int ierr = MPI_Irecv(shock_recvbuf.data()+particle_start,
                           ShockMessageBytes(count), MPI_BYTE,
                           recvs_thisrank[n].sendrank, 2, mpi_comm_part,
                           &(shock_recv_req[n]));
      if (ierr != MPI_SUCCESS) no_errors = false;
      particle_start += count;
    }
  }

  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::PackAndSendPrtcls()
//! \brief

TaskStatus ParticlesBoundaryValues::PackAndSendPrtcls() {
#if MPI_PARALLEL_ENABLED
  // Figure out how many particles will be sent from this ranks
  nprtcl_send=0;
  for (int n=0; n<nsends; ++n) {
    if (pmy_part->shock_enabled &&
        (sends_thisrank[n].nprtcls < 0 ||
         nprtcl_send > INT_MAX-sends_thisrank[n].nprtcls)) {
      std::cerr << "FATAL: shock send count exceeds INT_MAX" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    nprtcl_send += sends_thisrank[n].nprtcls;
  }

  bool no_errors=true;
  if (nprtcl_send > 0) {
    if (pmy_part->shock_enabled) {
      CheckShockBuffers(nprtcl_send, pmy_part->nrdata, pmy_part->nidata);
      Kokkos::realloc(shock_sendbuf, nprtcl_send);
      shock_send_req.assign(nsends, MPI_REQUEST_NULL);
    }
    // Allocate send buffer
    Kokkos::realloc(prtcl_rsendbuf, (pmy_part->nrdata)*nprtcl_send);
    Kokkos::realloc(prtcl_isendbuf, (pmy_part->nidata)*nprtcl_send);

    // sendlist on device is already sorted by destrank in CountSendAndRecvs()
    // Use sendlist on device to load particles into send buffer ordered by dest_rank
    int nrdata = pmy_part->nrdata;
    int nidata = pmy_part->nidata;
    auto &pr = pmy_part->prtcl_rdata;
    auto &pi = pmy_part->prtcl_idata;
    auto &rsendbuf = prtcl_rsendbuf;
    auto &isendbuf = prtcl_isendbuf;
    auto d_sendlist = sendlist.d_view;
    const bool shock_enabled = pmy_part->shock_enabled;
    auto states = pmy_part->shock_state;
    auto state_buffer = shock_sendbuf;
    par_for("ppack",DevExeSpace(),0,(nprtcl_send-1), KOKKOS_LAMBDA(const int n) {
      int p = d_sendlist(n).prtcl_indx;
      if (shock_enabled) state_buffer(n) = states(p);
      for (int i=0; i<nidata; ++i) {
        isendbuf(nidata*n + i) = pi(i,p);
      }
      for (int i=0; i<nrdata; ++i) {
        rsendbuf(nrdata*n + i) = pr(i,p);
      }
    });

    // Post non-blocking sends
    Kokkos::fence();
    rsend_req.clear();
    isend_req.clear();
    for (int n=0; n<nsends; ++n) {
      rsend_req.emplace_back(MPI_REQUEST_NULL);
      isend_req.emplace_back(MPI_REQUEST_NULL);
    }

    // Send Reals
    int data_start=0;
    for (int n=0; n<nsends; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = nrdata*(sends_thisrank[n].nprtcls);
      int data_end = data_start + nrdata*(sends_thisrank[n].nprtcls - 1);
      auto send_ptr = Kokkos::subview(prtcl_rsendbuf,std::make_pair(data_start,data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 0; // 0 for Reals, 1 for ints

      // Post non-blocking sends
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                           mpi_comm_part, &(rsend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start += data_size;
    }
    // Send ints
    data_start=0;
    for (int n=0; n<nsends; ++n) {
      // calculate amount of data to be passed, get pointer to variables
      int data_size = nidata*(sends_thisrank[n].nprtcls);
      int data_end = data_start + nidata*(sends_thisrank[n].nprtcls - 1);
      auto send_ptr = Kokkos::subview(prtcl_isendbuf,std::make_pair(data_start,data_end));
      int drank = sends_thisrank[n].recvrank;
      int tag = 1; // 0 for Reals, 1 for ints

      // Post non-blocking sends
      int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_INT, drank, tag,
                           mpi_comm_part, &(isend_req[n]));
      if (ierr != MPI_SUCCESS) {no_errors=false;}
      data_start += data_size;
    }
  }

  if (pmy_part->shock_enabled) {
    int particle_start = 0;
    for (int n=0; n<nsends; ++n) {
      const int count = sends_thisrank[n].nprtcls;
      int ierr = MPI_Isend(shock_sendbuf.data()+particle_start,
                           ShockMessageBytes(count), MPI_BYTE,
                           sends_thisrank[n].recvrank, 2, mpi_comm_part,
                           &(shock_send_req[n]));
      if (ierr != MPI_SUCCESS) no_errors = false;
      particle_start += count;
    }
  }

  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::RecvAndUnpackPrtcls()
//! \brief

TaskStatus ParticlesBoundaryValues::RecvAndUnpackPrtcls() {
  if (pmy_part->shock_enabled) return RecvAndUnpackShockPrtcls();
#if MPI_PARALLEL_ENABLED
  // Sort sendlist on host by index in particle array
  namespace KE = Kokkos::Experimental;
  std::sort(KE::begin(sendlist.h_view), KE::end(sendlist.h_view), SortByIndex);
  // sync sendlist host array with device.  This results in sorted array on device
  sendlist.template modify<HostMemSpace>();
  sendlist.template sync<DevExeSpace>();

  // increase size of particle arrays if needed
  int new_npart = pmy_part->nprtcl_thispack + (nprtcl_recv - nprtcl_send);
  if (nprtcl_recv > nprtcl_send) {
    Kokkos::resize(pmy_part->prtcl_idata, pmy_part->nidata, new_npart);
    Kokkos::resize(pmy_part->prtcl_rdata, pmy_part->nrdata, new_npart);
  }

  // check that particle communications have all completed
  bool bflag = false;
  bool no_errors=true;
  for (int n=0; n<nrecvs; ++n) {
    int test;
    int ierr = MPI_Test(&(rrecv_req[n]), &test, MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (!(static_cast<bool>(test))) {
      bflag = true;
    }
    ierr = MPI_Test(&(irecv_req[n]), &test, MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (!(static_cast<bool>(test))) {
      bflag = true;
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // exit if particle communications have not completed
  if (bflag) {return TaskStatus::incomplete;}

  // unpack particles into positions of sent particles
  if (nprtcl_recv > 0) {
    int nrdata = pmy_part->nrdata;
    int nidata = pmy_part->nidata;
    auto &pr = pmy_part->prtcl_rdata;
    auto &pi = pmy_part->prtcl_idata;
    auto &rrecvbuf = prtcl_rrecvbuf;
    auto &irecvbuf = prtcl_irecvbuf;
    int npart = pmy_part->nprtcl_thispack;
    int nsend = nprtcl_send;
    auto d_sendlist = sendlist.d_view;
    par_for("punpack",DevExeSpace(),0,(nprtcl_recv-1), KOKKOS_LAMBDA(const int n) {
      int p;
      if (n < nsend) {
        p = d_sendlist(n).prtcl_indx; // place particles in holes created by sends
      } else {
        p = npart + (n - nsend);           // place particle at end of arrays
      }
      for (int i=0; i<nidata; ++i) {
        pi(i,p) = irecvbuf(nidata*n + i);
      }
      for (int i=0; i<nrdata; ++i) {
        pr(i,p) = rrecvbuf(nrdata*n + i);
      }
    });
  }

  // At this point have filled npart_recv holes in particle arrays from sends
  // If (nprtcl_recv < nprtcl_send), have to move particles from end of arrays to fill
  // remaining holes
  int nremain = nprtcl_send - nprtcl_recv;
  if (nremain > 0) {
    int &npart = pmy_part->nprtcl_thispack;
    int i_last_hole = nprtcl_send-1;
    int i_next_hole = nprtcl_recv;
    for (int n=1; n<=nremain; ++n) {
      int nend = npart-n;
      if (nend > sendlist.h_view(i_last_hole).prtcl_indx) {
        // copy particle from end into hole
        int next_hole = sendlist.h_view(i_next_hole).prtcl_indx;
        auto rdest = Kokkos::subview(pmy_part->prtcl_rdata, Kokkos::ALL, next_hole);
        auto rsrc  = Kokkos::subview(pmy_part->prtcl_rdata, Kokkos::ALL, nend);
        Kokkos::deep_copy(rdest, rsrc);
        auto idest = Kokkos::subview(pmy_part->prtcl_idata, Kokkos::ALL, next_hole);
        auto isrc  = Kokkos::subview(pmy_part->prtcl_idata, Kokkos::ALL, nend);
        Kokkos::deep_copy(idest, isrc);
        i_next_hole += 1;
      } else {
        // this index contains a hole, so do nothing except find new index of last hole
        i_last_hole -= 1;
      }
    }

    // shrink size of particle data arrays
    Kokkos::resize(pmy_part->prtcl_idata, pmy_part->nidata, new_npart);
    Kokkos::resize(pmy_part->prtcl_rdata, pmy_part->nrdata, new_npart);
  }

  // Update nparticles_thisrank.  Update cost array (use npart_thismb[nmb]?)
  pmy_part->nprtcl_thispack = new_npart;
#endif
  pmy_part->pmy_pack->pmesh->UpdateParticleCounts();
  return TaskStatus::complete;
}

TaskStatus ParticlesBoundaryValues::RecvAndUnpackShockPrtcls() {
#if MPI_PARALLEL_ENABLED
  bool ready = true;
  for (int n=0; n<nrecvs; ++n) {
    MPI_Request *requests[3] = {&rrecv_req[n], &irecv_req[n], &shock_recv_req[n]};
    for (int i=0; i<3; ++i) {
      int done = 0;
      if (MPI_Test(requests[i], &done, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        std::cerr << "FATAL: testing shock particle receive failed" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!done) ready = false;
    }
  }
  if (!ready) return TaskStatus::incomplete;
#endif
  const int npart = pmy_part->nprtcl_thispack;
  const int64_t wide_count = static_cast<int64_t>(npart) + nprtcl_recv -
                            nprtcl_send - nprtcl_destroy;
  if (wide_count < 0 || wide_count > INT_MAX) {
    std::cerr << "FATAL: shock particle count exceeds supported range" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (npart < 0 || nprtcl_send < 0 || nprtcl_recv < 0 || nprtcl_destroy < 0 ||
      static_cast<int64_t>(nprtcl_send)+nprtcl_destroy > npart) {
    std::cerr << "FATAL: invalid shock particle removal count" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int new_npart = static_cast<int>(wide_count);
  if (nprtcl_send || nprtcl_recv || nprtcl_destroy) {
    std::vector<unsigned char> removed(npart, 0);
    for (int n=0; n<nprtcl_send+nprtcl_destroy; ++n) {
      const int p = n < nprtcl_send ? sendlist.h_view(n).prtcl_indx :
                                    shock_destroylist.h_view(n-nprtcl_send);
      if (p < 0 || p >= npart || removed[p]) {
        std::cerr << "FATAL: duplicate/invalid shock particle removal" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      removed[p] = 1;
    }
    const int nsurvive = npart-nprtcl_send-nprtcl_destroy;
    DvceArray1D<int> survivors("shock_survivors", nsurvive);
    auto indices = Kokkos::create_mirror_view(survivors);
    int next = 0;
    for (int p=0; p<npart; ++p) {
      if (!removed[p]) indices(next++) = p;
    }
    Kokkos::deep_copy(survivors, indices);
    const int nrdata = pmy_part->nrdata, nidata = pmy_part->nidata;
    DvceArray2D<Real> new_real("shock_particle_real", nrdata, new_npart);
    DvceArray2D<int> new_int("shock_particle_int", nidata, new_npart);
    DvceArray1D<shocks::EncounterState> new_state("shock_particle_state", new_npart);
    auto old_real = pmy_part->prtcl_rdata;
    auto old_int = pmy_part->prtcl_idata;
    auto old_state = pmy_part->shock_state;
    par_for("shock_compact",DevExeSpace(),0,nsurvive-1,
    KOKKOS_LAMBDA(const int n) {
      const int p = survivors(n);
      for (int i=0; i<nrdata; ++i) new_real(i,n) = old_real(i,p);
      for (int i=0; i<nidata; ++i) new_int(i,n) = old_int(i,p);
      new_state(n) = old_state(p);
    });
#if MPI_PARALLEL_ENABLED
    auto real_buffer = prtcl_rrecvbuf;
    auto int_buffer = prtcl_irecvbuf;
    auto state_buffer = shock_recvbuf;
    par_for("shock_unpack",DevExeSpace(),0,nprtcl_recv-1,
    KOKKOS_LAMBDA(const int n) {
      const int p = nsurvive+n;
      for (int i=0; i<nrdata; ++i) new_real(i,p) = real_buffer(nrdata*n+i);
      for (int i=0; i<nidata; ++i) new_int(i,p) = int_buffer(nidata*n+i);
      new_state(p) = state_buffer(n);
    });
#endif
    Kokkos::fence();
    pmy_part->prtcl_rdata = new_real;
    pmy_part->prtcl_idata = new_int;
    pmy_part->shock_state = new_state;
    pmy_part->nprtcl_thispack = new_npart;
  }
  // Exactly one count collective per completed receive path on every rank.
  pmy_part->pmy_pack->pmesh->UpdateParticleCounts();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::ClearPrtclSend()
//! \brief

TaskStatus ParticlesBoundaryValues::ClearPrtclSend() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  // wait for all non-blocking sends for vars to finish before continuing
  for (int n=0; n<nsends; ++n) {
    int ierr = MPI_Wait(&(rsend_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(isend_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (pmy_part->shock_enabled) {
      ierr = MPI_Wait(&(shock_send_req[n]), MPI_STATUS_IGNORE);
      if (ierr != MPI_SUCCESS) no_errors = false;
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in clearing sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rsend_req.clear();
  isend_req.clear();
  shock_send_req.clear();
#endif
  nsends=0;
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticlesBoundaryValues::ClearPrtclRecv()
//! \brief

TaskStatus ParticlesBoundaryValues::ClearPrtclRecv() {
#if MPI_PARALLEL_ENABLED
  bool no_errors=true;
  // wait for all non-blocking receives to finish before continuing
  for (int n=0; n<nrecvs; ++n) {
    int ierr = MPI_Wait(&(rrecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    ierr = MPI_Wait(&(irecv_req[n]), MPI_STATUS_IGNORE);
    if (ierr != MPI_SUCCESS) {no_errors=false;}
    if (pmy_part->shock_enabled) {
      ierr = MPI_Wait(&(shock_recv_req[n]), MPI_STATUS_IGNORE);
      if (ierr != MPI_SUCCESS) no_errors = false;
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in clearing receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  rrecv_req.clear();
  irecv_req.clear();
  shock_recv_req.clear();
#endif
  nrecvs=0;
  return TaskStatus::complete;
}

} // namespace particles
