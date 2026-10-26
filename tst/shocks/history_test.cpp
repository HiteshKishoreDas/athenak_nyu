// Copyright (C) 2026 AthenaK contributors. BSD-3-Clause.
// Protocol/state test only: no live outflow or hardware qualification implied.
#include <cassert>
#include <filesystem>
#include <iostream>
#include "shocks/history.hpp"

int main(int argc,char **argv) {
  global_variable::my_rank=0;global_variable::nranks=1;
  if (argc<2) return 2;
  shocks::HistoryStore store;
  store.directory=argv[1];store.stem="unit";store.lineage="unit";
  store.contract="typed-buffer-unit";store.capacity=2;
  store.Prepare(false);
  auto observe=[](uint64_t generation,Real time,bool inside) {
    shocks::Observation o;o.generation=generation;o.time=time;o.valid=true;o.inside=inside;
    o.rho=o.pressure=o.temperature=1;o.mach_pressure=o.mach_temperature=2;
    o.normal[0]=1;
    for (int d=0;d<3;++d) o.upstream[d]=o.downstream[d]=1;
    return o;
  };
  shocks::HistoryRecord r;r.emission=4;r.event.particle_id=17;r.event.sequence=5;
  r.event.flags=shocks::no_plateau|shocks::cadence_limited;
  r.event.pre=observe(0,0,false);r.event.first=observe(1,.1,true);
  r.event.last=observe(2,.2,true);r.event.post=observe(3,.3,false);
  r.event.representative=r.event.first;
  store.Push(r);
  if (argc>2 && std::string(argv[2])=="duplicate") {
    store.Push(r);store.Save(4,1);return 3;
  }
  auto exit=r;exit.emission=5;exit.event.particle_id=23;
  exit.event.flags|=shocks::right_censored|shocks::domain_exit;
  exit.event.status=shocks::partial;exit.event.post=shocks::Observation{};
  exit.event.post.generation=5;exit.event.post.time=.5;exit.event.post.x[0]=2;
  store.Push(exit);
  if (argc>2 && std::string(argv[2])=="overflow") {store.Push(r);return 3;}
  const auto packet=store.Save(5,1);
  shocks::HistoryStore restored=store;restored.pending.clear();restored.Restore(packet,5,1);
  // No live particle collection is involved: completed records survive deletion.
  assert(restored.pending.size()==2 && restored.pending[1].event.particle_id==23);
  assert(restored.pending[1].event.post.time==.5 && restored.pending[1].emission==5);
  restored.Flush(5,1);const auto digest=restored.digest,last=restored.last_file;
  restored.Flush(5,1);assert(restored.digest==digest && restored.last_file==last);
  auto next=r;next.emission=7;next.event.sequence=6;restored.Push(next);
  const auto pending=restored.Save(8,2);
  shocks::HistoryStore again=restored;again.pending.clear();again.Restore(pending,8,2);
  assert(again.pending.size()==1 && again.pending[0].emission==7);
  again.Flush(10,3);again.Verify();
  if (argc>2 && std::string(argv[2])=="corrupt") {
    std::ofstream out(again.Path(5),std::ios::app|std::ios::binary);out << 'x';out.close();
    again.Verify();return 3;
  }
  auto curve=store;curve.curves=true;curve.bits=0;curve.stem="unitcurve";
  curve.contract="curve-node-history-v4 0 typed-buffer-unit";curve.pending.clear();
  auto node=r;node.event.particle_id=UINT32_MAX;node.event.sequence=UINT64_C(9007199254740993);
  node.birth=0;node.first_level=1;node.last_level=2;
  for (int d=0;d<3;++d) node.first_spacing[d]=node.last_spacing[d]=.1;
  curve.Push(node);const auto typed=curve.Save(4,1);curve.pending.clear();curve.Restore(typed,4,1);
  assert(curve.pending[0].event.particle_id==UINT32_MAX &&
         curve.pending[0].event.sequence==UINT64_C(9007199254740993));
  curve.Flush(4,1);
  std::cout << "BINH01-05 queued destroyed identity, censored exit, pending restart, exact replay/prefix PASS\n";
}
