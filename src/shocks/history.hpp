// Copyright (C) 2026 AthenaK contributors. BSD-3-Clause.
#ifndef SHOCKS_HISTORY_HPP_
#define SHOCKS_HISTORY_HPP_
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include "parameter_input.hpp"
#include "particles/tracer_fields.hpp"
#include "shocks/store.hpp"

namespace shocks {
constexpr uint64_t history_offset=UINT64_C(14695981039346656037);

// Match the native scheduler; an inactive block must never enable tracking.
inline std::string HistoryOutput(ParameterInput *pin, bool curves) {
  std::string found;
  for (const auto &block : pin->block) {
    const auto &name=block.block_name;
    if (name.rfind("output",0)!=0 || !pin->DoesParameterExist(name,"file_type")) continue;
    const auto type=pin->GetString(name,"file_type");
    const bool curve_output=type=="curve" || type=="crv";
    if (!curve_output && type!="prtcl_thermo_history") continue;
    if (pin->GetOrAddString(name,"mode","snapshot")!="encounters") continue;
    const bool active=pin->DoesParameterExist(name,"dcycle") ?
        pin->GetInteger(name,"dcycle")!=0 : pin->GetReal(name,"dt")>0;
    if (!active) continue;
    if (curve_output!=curves) continue;
    if (!pin->DoesParameterExist(name,"variables")) {
      Fail(name+"/variables must explicitly include shock or mach for mode=encounters");
    }
    const bool mhd=pin->DoesBlockExist("mhd");
    const auto physics=mhd ? "mhd" : "hydro";
    const int scalars=pin->DoesParameterExist(physics,"nscalars") ?
        pin->GetInteger(physics,"nscalars") : 0;
    const auto fields=particles::ParseTracerFieldList(pin->GetString(name,"variables"),
        mhd,scalars,name+"/variables");
    const bool requested=std::any_of(fields.begin(),fields.end(),[](const auto &field) {
      return particles::IsShockField(field.kind) || field.kind==particles::TracerFieldKind::mach;
    });
    if (!requested) Fail(name+"/variables must include shock or mach for mode=encounters");
    if (!found.empty()) Fail("only one active encounter output block per carrier is supported");
    found=name;
  }
  return found;
}

struct HistoryRecord {
  uint64_t emission=0;
  Event event;
  Real birth=-1;
  int first_level=-1, last_level=-1;
  Real first_spacing[3]={}, last_spacing[3]={};
};

struct HistoryBatch {
  uint32_t kind=0;
  uint64_t generation=0, previous=UINT64_MAX, previous_file=0, prefix=history_offset;
  Real time=0;
  std::vector<HistoryRecord> records;
};

// One bounded queue per carrier, consumed only by its explicit output block.
class HistoryStore {
 public:
  bool curves=false;
  unsigned bits=0;
  int capacity=0;
  std::string directory, stem, lineage, contract;
  uint64_t last_flush=UINT64_MAX, last_file=0, digest=history_offset;
  std::vector<HistoryRecord> pending;

  void Configure(ParameterInput *pin, const std::string &block, bool curve_nodes,
                 unsigned curve_bits, int cap, const std::string &dir,
                 const std::string &line, const std::string &definition) {
    curves=curve_nodes; bits=curve_bits; capacity=cap; directory=dir;
    lineage=line; contract=definition;
    const auto type=pin->GetString(block,"file_type");
    const auto id=pin->GetOrAddString(block,"id",type);
    if (pin->GetOrAddString(block,"history_id",id)!=id) Fail("encounter output id changed on restart");
    stem=pin->GetOrAddString(block,"history_stem",pin->GetString("job","basename")+"."+id);
    if (stem.empty() || stem.size()>4096 || stem.find('/')!=std::string::npos ||
        stem.find('\\')!=std::string::npos) Fail("invalid encounter output basename/id");
  }
  std::string Path(uint64_t generation, bool suspended=false) const {
    return directory+"/"+stem+(suspended ? ".suspended" : ".cycle")+
        std::to_string(generation)+RankName()+"of"+std::to_string(global_variable::nranks)+
        (curves ? ".crv" : ".thp");
  }
  void Push(const HistoryRecord &record) {
    if (pending.size()>=static_cast<size_t>(capacity)) {
      Fail("unflushed encounter buffer overflow; replay prior checkpoint with larger event_capacity");
    }
    pending.push_back(record);
  }
  std::string Encode(const std::vector<HistoryRecord> &records, uint32_t kind,
                     uint64_t generation, Real time) const {
    if (records.size()>static_cast<size_t>(capacity)) Fail("encounter batch exceeds capacity");
    auto ordered=records;
    std::sort(ordered.begin(),ordered.end(),[](const auto &a,const auto &b) {
      return std::tie(a.emission,a.event.particle_id,a.event.sequence)<
             std::tie(b.emission,b.event.particle_id,b.event.sequence);
    });
    std::string data(16,'\0'); std::memcpy(data.data(),"ATHK_EVENTS",11);
    auto put=[&](const auto &v) {data.append(reinterpret_cast<const char *>(&v),sizeof(v));};
    auto u32=[&](uint32_t v) {put(v);}; auto u64=[&](uint64_t v) {put(v);};
    auto i32=[&](int32_t v) {put(v);};
    auto text=[&](const std::string &v) {u32(v.size());data+=v;};
    u32(1);u32(0x01020304);u32(sizeof(Real));u32(curves);u32(bits);u32(kind);
    u64(generation);u64(last_flush);u64(last_file);u64(digest);u64(ordered.size());
    put(time);i32(global_variable::my_rank);i32(global_variable::nranks);
    text(lineage);text(contract);text(stem);
    std::set<std::pair<uint64_t,uint64_t>> identities;
    for (const auto &r : ordered) {
      if (!identities.emplace(r.event.particle_id,r.event.sequence).second) Fail("duplicate buffered event identity");
      CheckRecord(r,generation,time);
      u64(r.emission);
      const auto id=r.event.particle_id;
      if (id>(curves ? UINT32_MAX : INT32_MAX)) Fail("event ID exceeds signed32 namespace");
      i32(id<=INT32_MAX ? id : static_cast<int64_t>(id)-(INT64_C(1)<<32));
      u64(r.event.sequence);u32(r.event.status);u32(r.event.flags);
      for (const auto *o : {&r.event.pre,&r.event.first,&r.event.last,&r.event.post,
                            &r.event.representative}) {
        put(o->time);u64(o->generation);u32(o->valid);u32(o->inside);
        for (Real v : o->x) put(v);
        for (int64_t v : o->image) put(v);
        for (Real v : {o->rho,o->pressure,o->temperature,o->mach_pressure,o->mach_temperature}) put(v);
        for (const auto *values : {o->normal,o->upstream,o->downstream}) {
          for (int d=0;d<3;++d) put(values[d]);
        }
      }
      put(r.birth);i32(r.first_level);i32(r.last_level);
      for (const auto *values : {r.first_spacing,r.last_spacing}) {
        for (int d=0;d<3;++d) put(values[d]);
      }
    }
    u64(Digest(history_offset,data.data(),data.size()));
    return data;
  }
  HistoryBatch Decode(const std::string &data) const {
    if (data.size()<104 || data.compare(0,16,std::string("ATHK_EVENTS\0\0\0\0\0",16))) {
      Fail("invalid binary encounter header");
    }
    size_t offset=16;
    auto get=[&](auto &value) {
      if (offset>data.size()-8 || sizeof(value)>data.size()-8-offset) Fail("truncated encounter data");
      std::memcpy(&value,data.data()+offset,sizeof(value));offset+=sizeof(value);
    };
    auto u32=[&]() {uint32_t v;get(v);return v;};
    auto u64=[&]() {uint64_t v;get(v);return v;};
    auto i32=[&]() {int32_t v;get(v);return v;};
    auto real=[&]() {Real v;get(v);return v;};
    uint64_t footer;std::memcpy(&footer,data.data()+data.size()-8,8);
    if (footer!=Digest(history_offset,data.data(),data.size()-8) || u32()!=1 ||
        u32()!=0x01020304 || u32()!=sizeof(Real) || u32()!=curves || u32()!=bits) {
      Fail("binary encounter checksum/schema mismatch");
    }
    HistoryBatch b;b.kind=u32();b.generation=u64();b.previous=u64();
    b.previous_file=u64();b.prefix=u64();const auto count=u64();b.time=real();
    const auto rank=i32(),ranks=i32();
    auto text=[&](size_t cap) {
      const auto count=u32();
      if (count>cap || offset>data.size()-8 || count>data.size()-8-offset) Fail("invalid encounter text extent");
      auto v=data.substr(offset,count);offset+=count;return v;
    };
    if (text(128)!=lineage || text(16384)!=contract || text(4096)!=stem ||
        rank!=global_variable::my_rank || ranks!=global_variable::nranks ||
        b.kind>2 || count>static_cast<uint64_t>(capacity) ||
        !std::isfinite(b.time) || b.time<0 ||
        (b.previous!=UINT64_MAX && b.previous>b.generation)) Fail("binary encounter metadata mismatch");
    std::set<std::pair<uint64_t,uint64_t>> identities;
    for (uint64_t n=0;n<count;++n) {
      HistoryRecord r;r.emission=u64();const int32_t id=i32();
      if ((!curves && id<0) || (curves && id==0)) Fail("invalid binary event identity");
      r.event.particle_id=curves ? static_cast<uint32_t>(id) : static_cast<uint64_t>(id);
      r.event.sequence=u64();r.event.status=u32();r.event.flags=u32();
      for (auto *o : {&r.event.pre,&r.event.first,&r.event.last,&r.event.post,&r.event.representative}) {
        o->time=real();o->generation=u64();const auto valid=u32(),inside=u32();
        if (valid>1 || inside>1 || inside>valid) Fail("invalid binary observation flags");
        o->valid=valid;o->inside=inside;
        for (Real &v : o->x) v=real();
        for (int64_t &v : o->image) get(v);
        o->rho=real();o->pressure=real();o->temperature=real();
        o->mach_pressure=real();o->mach_temperature=real();
        for (auto *values : {o->normal,o->upstream,o->downstream}) for (int d=0;d<3;++d) values[d]=real();
      }
      r.birth=real();r.first_level=i32();r.last_level=i32();
      for (auto *values : {r.first_spacing,r.last_spacing}) for (int d=0;d<3;++d) values[d]=real();
      CheckRecord(r,b.generation,b.time);
      if ((b.kind==1 && (r.emission!=b.generation || r.event.status!=partial ||
                        !(r.event.flags&right_censored))) ||
          (b.kind!=1 && (r.event.flags&right_censored) && !(r.event.flags&domain_exit))) {
        Fail("invalid censored event batch kind");
      }
      if (!identities.emplace(r.event.particle_id,r.event.sequence).second) Fail("duplicate pending event identity");
      if (b.kind!=1 && b.previous!=UINT64_MAX && r.emission<=b.previous) Fail("buffered event predates flush boundary");
      b.records.push_back(r);
    }
    if (offset!=data.size()-8) Fail("binary encounter trailing data");
    return b;
  }
  std::string Read(const std::string &path) const {
    std::error_code error;const auto bytes=std::filesystem::file_size(path,error);
    if (error || bytes>22000+static_cast<uint64_t>(capacity)*1100) Fail("missing/oversized encounter file: "+path);
    std::ifstream in(path,std::ios::binary);std::string data(bytes,'\0');
    if (!in.read(data.data(),data.size())) Fail("cannot read encounter file: "+path);
    return data;
  }
  void Verify(bool reject_orphans=false) const {
    std::set<std::string> committed;
    Real later_time=std::numeric_limits<Real>::infinity();
    auto current=last_flush,expected_file=last_file,expected_digest=digest;
    while (current!=UINT64_MAX) {
      committed.insert(std::filesystem::path(Path(current)).filename().string());
      const auto data=Read(Path(current));const auto b=Decode(data);
      if (b.kind!=0 || b.generation!=current || b.time>=later_time ||
          Digest(history_offset,data.data(),data.size())!=expected_file ||
          Digest(b.prefix,data.data(),data.size())!=expected_digest ||
          (b.previous!=UINT64_MAX && b.previous>=current)) Fail("committed encounter prefix corrupted");
      later_time=b.time;current=b.previous;expected_file=b.previous_file;expected_digest=b.prefix;
    }
    if (expected_file!=0 || expected_digest!=history_offset) Fail("missing initial encounter prefix");
    if (reject_orphans && std::filesystem::exists(directory)) {
      const auto prefix=stem+".cycle";
      const auto suffix=RankName()+"of"+std::to_string(global_variable::nranks)+(curves ? ".crv" : ".thp");
      for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        const auto name=entry.path().filename().string();
        if (name.rfind(prefix,0)==0 && name.size()>=suffix.size() &&
            name.compare(name.size()-suffix.size(),suffix.size(),suffix)==0 && !committed.count(name)) {
          Fail("encounter history has a non-prefix committed tail; restart with only the checkpoint prefix in a new directory");
        }
      }
    }
  }
  void Prepare(bool restart) {
    std::error_code error;std::filesystem::create_directories(directory,error);
    if (error) Fail("cannot create binary encounter directory");
    if (restart) Verify(true);
  }
  void Flush(uint64_t generation, Real time) {
    if (last_flush==generation) {
      if (!pending.empty() || FileDigest(Path(generation))!=last_file) Fail("conflicting same-cycle encounter replay");
      return;
    }
    if (last_flush!=UINT64_MAX && generation<last_flush) Fail("encounter flush clock reversed");
    const auto data=Encode(pending,0,generation,time);
    CommitFile(Path(generation),data,true);
    digest=Digest(digest,data.data(),data.size());
    last_file=Digest(history_offset,data.data(),data.size());last_flush=generation;
    pending.clear();
  }
  std::string Save(uint64_t generation, Real time) const {return Encode(pending,2,generation,time);}
  void Restore(const std::string &data, uint64_t generation, Real time) {
    const auto b=Decode(data);
    if (b.kind!=2 || b.generation!=generation || b.time!=time) Fail("pending encounter checkpoint clock mismatch");
    last_flush=b.previous;last_file=b.previous_file;digest=b.prefix;pending=b.records;
    Verify();
  }
  void Suspended(const std::vector<HistoryRecord> &records,uint64_t generation,Real time) const {
    CommitFile(Path(generation,true),Encode(records,1,generation,time),true);
  }
 private:
  void CheckRecord(const HistoryRecord &r,uint64_t generation,Real time) const {
    const auto &e=r.event;
    if (r.emission>generation || e.status!=static_cast<uint32_t>(bool(e.flags&(left_censored|right_censored|domain_exit))) || (e.flags&~1023U) ||
        (e.status==completed && (e.flags&(left_censored|right_censored|domain_exit))) ||
        !e.first.valid || !e.first.inside || !e.last.valid || !e.last.inside ||
        e.last.time<e.first.time || (e.pre.valid && e.pre.time>e.first.time) ||
        (e.post.valid && e.post.time<e.last.time) || (e.status==completed && !e.post.valid)) Fail("invalid buffered encounter evidence");
    for (const auto *o : {&e.pre,&e.first,&e.last,&e.post,&e.representative}) {
      if (!std::isfinite(o->time) || o->time < -1 || o->time>time ||
          (o->time!=-1 && o->generation>r.emission) || (o->valid && o->time<0)) Fail("invalid buffered observation clock");
      for (Real v : o->x) if (!std::isfinite(v)) Fail("nonfinite buffered position");
      if (!o->valid) continue;
      if (o->rho<=0 || o->pressure<=0 || o->temperature<=0 ||
          !std::isfinite(o->rho) || !std::isfinite(o->pressure) || !std::isfinite(o->temperature) ||
          std::isinf(o->mach_pressure) || std::isinf(o->mach_temperature)) Fail("invalid buffered fluid evidence");
    }
    for (const auto *o : {&e.pre,&e.first,&e.last,&e.post,&e.representative}) {
      if (!o->valid) continue;
      for (int d=0;d<3;++d) {
        if (!std::isfinite(o->normal[d]) || !std::isfinite(o->upstream[d]) ||
            !std::isfinite(o->downstream[d]) || o->upstream[d]<=0 || o->downstream[d]<=0) {
          Fail("invalid buffered endpoint evidence");
        }
      }
    }
    if (!curves) {
      if (r.birth!=-1 || r.first_level!=-1 || r.last_level!=-1) Fail("invalid particle event metadata");
      for (int d=0;d<3;++d) if (r.first_spacing[d]!=0 || r.last_spacing[d]!=0) Fail("invalid particle event spacing");
    }
    if (curves) for (int d=0;d<3;++d) {
      if (!std::isfinite(r.first_spacing[d]) || !std::isfinite(r.last_spacing[d]) ||
          r.first_spacing[d]<=0 || r.last_spacing[d]<=0) Fail("invalid buffered mesh spacing");
    }
    if (curves && (!std::isfinite(r.birth) || r.birth<0 || r.birth>e.first.time ||
                   r.first_level<0 || r.last_level<0)) Fail("invalid buffered curve metadata");
  }
};
}  // namespace shocks
#endif
