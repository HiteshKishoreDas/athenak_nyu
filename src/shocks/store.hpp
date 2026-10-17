// Copyright (C) 2026 the AthenaK contributors
// Licensed under the 3-clause BSD License (the "LICENSE").
#ifndef SHOCKS_STORE_HPP_
#define SHOCKS_STORE_HPP_

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "../athena.hpp"
#include "../globals.hpp"
#include "shocks/encounter.hpp"

namespace shocks {
[[noreturn]] inline void Fail(const std::string &message) {
  std::cerr << "### FATAL shock tracking: " << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

inline uint64_t Digest(uint64_t h, const char *p, size_t count) {
  for (size_t n=0; n<count; ++n) {
    h ^= static_cast<unsigned char>(p[n]);
    h *= UINT64_C(1099511628211);
  }
  return h;
}
inline uint64_t FileDigest(const std::string &path, uint64_t limit = UINT64_MAX) {
  std::ifstream in(path, std::ios::binary);
  if (!in) Fail("cannot read integrity input: " + path);
  uint64_t hash = UINT64_C(14695981039346656037);
  char buffer[65536];
  while (limit > 0) {
    auto amount = static_cast<std::streamsize>(std::min<uint64_t>(limit,sizeof(buffer)));
    in.read(buffer,amount);
    if (!in.gcount()) break;
    hash = Digest(hash,buffer,in.gcount());
    limit -= in.gcount();
  }
  if (in.bad()) Fail("integrity read failed: " + path);
  return hash;
}

inline std::string RankName() {
  return ".rank" + std::to_string(global_variable::my_rank);
}

// Immutable, durable file commit. Existing files permit only byte-identical replay.
// A torn temp is ignored. A torn final, differing continuation, or lost file fails.
inline void CommitFile(const std::string &path, const std::string &data, bool replay) {
  if (std::filesystem::exists(path)) {
    std::ifstream in(path, std::ios::binary);
    std::string old(data.size(), '\0');
    if (!replay || !in || !in.read(old.data(),old.size()) ||
        in.peek() != std::ifstream::traits_type::eof() || old != data) {
      Fail("existing file differs: " + path);
    }
    return;
  }
  const std::string temp = path + ".tmp." + std::to_string(getpid());
  int fd = open(temp.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0664);
  if (fd < 0) Fail("cannot create temporary output: " + temp);
  size_t offset = 0;
  while (offset < data.size()) {
    ssize_t n = write(fd, data.data() + offset, data.size() - offset);
    if (n <= 0) Fail("failed output write: " + temp);
    offset += n;
  }
  if (fsync(fd) || close(fd)) Fail("failed output flush: " + temp);
  if (link(temp.c_str(), path.c_str())) Fail("failed immutable commit: " + path);
  if (unlink(temp.c_str())) Fail("failed temporary cleanup: " + temp);
  auto dir = std::filesystem::path(path).parent_path();
  int dfd = open(dir.empty() ? "." : dir.c_str(), O_RDONLY);
  if (dfd < 0 || fsync(dfd) || close(dfd)) Fail("failed directory flush: " + path);
}

inline void PutObservation(std::ostream &out, const shocks::Observation &o) {
  out << o.time << ',' << o.generation << ',' << o.valid << ',' << o.inside;
  for (Real x : o.x) out << ',' << x;
  for (int64_t x : o.image) out << ',' << x;
  out << ',' << o.rho << ',' << o.pressure << ',' << o.temperature
      << ',' << o.mach_pressure << ',' << o.mach_temperature;
  for (Real x : o.normal) out << ',' << x;
  for (Real x : o.upstream) out << ',' << x;
  for (Real x : o.downstream) out << ',' << x;
}

inline void PutEvent(std::ostream &out, const shocks::Event &e, bool newline = true) {
  out << e.particle_id << ',' << e.sequence << ',' << e.status << ',' << e.flags;
  for (const auto *o : {&e.pre, &e.first, &e.last, &e.post, &e.representative}) {
    out << ',';
    PutObservation(out,*o);
  }
  if (newline) out << '\n';
}

}  // namespace shocks
#endif  // SHOCKS_STORE_HPP_
