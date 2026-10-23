// Copyright (C) 2026 the AthenaK contributors. BSD-3-Clause.
#ifndef OUTPUTS_TRACER_BINARY_HPP_
#define OUTPUTS_TRACER_BINARY_HPP_

#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "parameter_input.hpp"
#include "shock_version.hpp"  // NOLINT(build/include_subdir): generated source identity
#include "shocks/store.hpp"

namespace tracer_output {
inline bool Binary(ParameterInput *pin) {
  const auto format=pin->GetOrAddString("job","tracer_output_format","legacy");
  if (format!="legacy" && format!="binary") {
    shocks::Fail("job/tracer_output_format must be legacy or binary");
  }
  return format=="binary";
}
inline std::string Extension(bool binary) {return binary ? ".tbin" : ".csv";}

// AthenaK convention: a versioned named header and native binary numeric data.
// No C++ structs/padding or integer identities serialized through Real.
class Table {
 public:
  Table(const std::string &kind, const std::string &columns, uint64_t count,
        uint64_t cycle, Real time, const std::string &metadata="") : count_(count) {
    std::istringstream fields(columns);
    std::string field;
    while (std::getline(fields,field,',')) {
      const auto pos=field.find(':');
      if (pos==std::string::npos || pos==0) shocks::Fail("invalid tracer column");
      const auto type=field.substr(pos+1);
      if (type!="u64" && type!="i64" && type!="u32" && type!="i32" &&
          type!="real") shocks::Fail("invalid tracer column type");
      types_.push_back(type);
      width_ += (type=="real" ? sizeof(Real) : (type=="u64" || type=="i64") ? 8 : 4);
    }
    const auto header="kind="+kind+"\ncolumns="+columns+"\nsource="+
        SHOCK_IMPLEMENTATION_HASH+"\nunits=code\n"+metadata;
    if (header.size()>65536 || types_.empty() ||
        count>(std::numeric_limits<size_t>::max()-header.size()-128)/width_) {
      shocks::Fail("tracer binary extent/metadata overflow");
    }
    data_.reserve(header.size()+128+count*width_);
    char magic[32]="Athena tracer binary output";
    data_.append(magic,sizeof(magic));
    Append(uint32_t{1}); Append(uint32_t{0x12345678});
    Append(uint32_t{sizeof(Real)}); Append(static_cast<uint32_t>(header.size()));
    Append(count); Append(cycle); Append(time);
    Append(static_cast<int32_t>(global_variable::my_rank));
    Append(static_cast<int32_t>(global_variable::nranks));
    data_ += header;
  }

  template<typename T> void Put(T value) {
    const std::string type=std::is_same<T,Real>::value ? "real" :
        std::is_integral<T>::value && sizeof(T)==8 ?
        (std::is_signed<T>::value ? "i64" : "u64") :
        std::is_integral<T>::value && sizeof(T)==4 ?
        (std::is_signed<T>::value ? "i32" : "u32") : "unsupported";
    if (values_/types_.size()>=count_ || type!=types_[values_%types_.size()]) {
      shocks::Fail("tracer binary record type/extent mismatch");
    }
    Append(value); ++values_;
  }
  std::string Finish() {
    if (values_%types_.size() || values_/types_.size()!=count_) {
      shocks::Fail("incomplete tracer binary table");
    }
    Append(shocks::Digest(UINT64_C(14695981039346656037),data_.data(),data_.size()));
    return std::move(data_);
  }

 private:
  template<typename T> void Append(T value) {
    data_.append(reinterpret_cast<const char *>(&value),sizeof(value));
  }
  uint64_t count_, values_=0;
  size_t width_=0;
  std::vector<std::string> types_;
  std::string data_;
};

inline std::string ObservationColumns(const std::string &prefix) {
  std::string result;
  for (const auto &field : {"time:real","generation:u64","valid:u32","inside:u32",
      "x1:real","x2:real","x3:real","image1:i64","image2:i64","image3:i64",
      "density:real","pressure:real","temperature:real","mach_pressure:real",
      "mach_temperature:real","normal1:real","normal2:real","normal3:real",
      "upstream_density:real","upstream_pressure:real","upstream_temperature:real",
      "downstream_density:real","downstream_pressure:real","downstream_temperature:real"}) {
    if (!result.empty()) result += ',';
    result += prefix+field;
  }
  return result;
}
inline std::string EventColumns(bool curve) {
  std::string result="particle_id:u64,sequence:u64,status:u32,flags:u32";
  for (const auto *prefix : {"pre.","first.","last.","post.","representative."}) {
    result += ','+ObservationColumns(prefix);
  }
  if (curve) result += ",birth:real,first_level:i32,first_dx1:real,first_dx2:real,"
      "first_dx3:real,last_level:i32,last_dx1:real,last_dx2:real,last_dx3:real";
  return result;
}
inline std::string TraceColumns(bool curve) {
  auto result=std::string("particle_id:u64,sequence:u64,phase:u32,flags:u32,")+
      ObservationColumns("");
  if (curve) result += ",birth:real,level:i32,dx1:real,dx2:real,dx3:real";
  return result;
}
inline void Observation(Table &out, const shocks::Observation &o) {
  out.Put(o.time); out.Put(o.generation);
  out.Put(static_cast<uint32_t>(o.valid)); out.Put(static_cast<uint32_t>(o.inside));
  for (Real x : o.x) out.Put(x);
  for (int64_t x : o.image) out.Put(x);
  for (Real x : {o.rho,o.pressure,o.temperature,o.mach_pressure,o.mach_temperature})
    out.Put(x);
  for (Real x : o.normal) out.Put(x);
  for (Real x : o.upstream) out.Put(x);
  for (Real x : o.downstream) out.Put(x);
}
inline void Event(Table &out, const shocks::Event &e) {
  out.Put(e.particle_id); out.Put(e.sequence);
  out.Put(static_cast<uint32_t>(e.status)); out.Put(e.flags);
  for (const auto *o : {&e.pre,&e.first,&e.last,&e.post,&e.representative})
    Observation(out,*o);
}
}  // namespace tracer_output
#endif  // OUTPUTS_TRACER_BINARY_HPP_
