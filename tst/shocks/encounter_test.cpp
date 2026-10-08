//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include <cassert>
#include <iostream>
#include <limits>
#include "shocks/encounter.hpp"

shocks::Observation Observe(uint64_t n, bool inside, bool valid=true) {
  shocks::Observation o;
  o.generation=n;
  o.time=0.1*n;
  o.valid=valid;
  o.inside=inside;
  o.rho=inside ? 2 : 1;
  o.mach_pressure=2;
  o.normal[0]=1;
  return o;
}
int main(int argc, char **argv) {
  Kokkos::initialize(argc,argv);
  {
    // State logic also executes in a real device kernel, not only host asserts.
    DvceArray1D<int> result("state_checks",1);
    Kokkos::parallel_for("encounter_checks",1,KOKKOS_LAMBDA(int) {
      int checks=0;
      shocks::EncounterState s;
      shocks::Event e;
      for (uint64_t n=0; n<4; ++n) {
        shocks::Observation o; o.time=n; o.generation=n; o.valid=true;
        if (!shocks::UpdateEncounter(s,o,e) && s.phase==shocks::outside) ++checks;
      }
      result(0)=checks;
    });
    auto h=Kokkos::create_mirror_view_and_copy(HostMemSpace(),result);
    assert(h(0)==4);
    shocks::EncounterState s;
    s.particle_id=UINT64_C(9007199254740993);
    s.sequence=UINT64_C(9007199254740993);
    shocks::Event e;
    assert(!shocks::UpdateEncounter(s,Observe(0,false),e));
    assert(!shocks::UpdateEncounter(s,Observe(1,true),e));
    for (uint64_t n=2; n<7; ++n) assert(!shocks::UpdateEncounter(s,Observe(n,true),e));
    assert(!shocks::UpdateEncounter(s,Observe(7,false),e));
    assert(s.phase==shocks::pending_exit);
    // Duplicate generation cannot satisfy the two-observation rearm rule.
    assert(!shocks::UpdateEncounter(s,Observe(7,false),e));
    assert(s.clean==1);
    assert(shocks::UpdateEncounter(s,Observe(8,false),e));
    assert(e.particle_id==UINT64_C(9007199254740993));
    assert(e.sequence==UINT64_C(9007199254740993));
    assert(e.pre.time==0 && e.first.time==0.1 && e.last.generation==6);
    assert(e.post.generation==7 && e.representative.generation==1);
    assert(e.status==shocks::completed && (e.flags & shocks::no_plateau));
    assert(s.sequence==UINT64_C(9007199254740994));
    assert(!shocks::UpdateEncounter(s,Observe(9,true),e));
    assert(!shocks::UpdateEncounter(s,Observe(10,false),e));
    assert(!shocks::UpdateEncounter(s,Observe(11,true),e));
    assert(s.flags & shocks::merged_flicker);
    assert(!shocks::UpdateEncounter(s,Observe(12,false,false),e));
    assert(s.phase==shocks::active && (s.flags & shocks::invalid_observation));
    assert(!shocks::UpdateEncounter(s,Observe(13,false),e));
    assert(!shocks::UpdateEncounter(s,Observe(14,false,false),e));
    assert(!shocks::UpdateEncounter(s,Observe(15,false),e));
    assert(shocks::UpdateEncounter(s,Observe(16,false),e));
    assert(e.post.generation==13 && (e.flags & shocks::merged_flicker));
    shocks::EncounterState left;
    assert(!shocks::UpdateEncounter(left,Observe(0,true),e));
    assert(!shocks::UpdateEncounter(left,Observe(1,false),e));
    assert(shocks::UpdateEncounter(left,Observe(2,false),e));
    assert(e.status==shocks::partial && (e.flags & shocks::left_censored));
    assert(!e.pre.valid && e.pre.time==-1);
    shocks::EncounterState open;
    assert(!shocks::UpdateEncounter(open,Observe(0,true),e));
    shocks::FillEvent(open,e,shocks::right_censored | shocks::domain_exit);
    assert(e.status==shocks::partial && !e.post.valid);
    assert(open.phase==shocks::active && open.sequence==0);
    // EVT09: opposing normals are ambiguous, never silently coherent.
    shocks::EncounterState reversed;
    shocks::UpdateEncounter(reversed,Observe(0,false),e);
    shocks::UpdateEncounter(reversed,Observe(1,true),e);
    auto negative=Observe(2,true); negative.normal[0]=-1;
    shocks::UpdateEncounter(reversed,negative,e);
    assert(reversed.flags & shocks::direction_ambiguous);
    // Downstream numerical lobes are not a second forward material passage.
    shocks::EncounterState echo;
    auto pre=Observe(0,false); pre.rho=2; pre.pressure=4;
    shocks::UpdateEncounter(echo,pre,e);
    auto lobe=Observe(1,true); lobe.upstream[0]=1; lobe.upstream[1]=1;
    lobe.downstream[0]=2; lobe.downstream[1]=4;
    shocks::UpdateEncounter(echo,lobe,e);
    assert(echo.phase==shocks::outside && echo.suppressed_reentries==1);
    lobe.generation=2;
    shocks::UpdateEncounter(echo,lobe,e);
    assert(echo.suppressed_reentries==1);
    pre.generation=3;
    shocks::UpdateEncounter(echo,pre,e);
    shocks::EncounterState saturated=echo;
    saturated.suppressing_reentry=false;
    saturated.suppressed_reentries=UINT64_MAX;
    lobe.generation=5;
    shocks::UpdateEncounter(saturated,lobe,e);
    assert(saturated.suppressed_reentries==UINT64_MAX);
    assert(saturated.flags & shocks::invalid_observation);
    // A new front has the previous downstream state as its own upstream.
    lobe.generation=4; lobe.upstream[0]=2; lobe.upstream[1]=4;
    lobe.downstream[0]=3; lobe.downstream[1]=8;
    shocks::UpdateEncounter(echo,lobe,e);
    assert(echo.phase==shocks::active);
    // Known-seed adversarial sequence: emitted identities strictly increase.
    shocks::EncounterState random;
    uint64_t seed=72941, expected=0;
    for (uint64_t n=0; n<10000; ++n) {
      seed=6364136223846793005ULL*seed+1;
      if (shocks::UpdateEncounter(random,Observe(n,(seed>>61)==0,(seed>>60)!=0),e)) {
        assert(e.sequence==expected++);
      }
    }
    assert(random.sequence==expected);
    std::cout << "EVT01-08, EVT11-14 state checks PASS; seed72941; large typed IDs\n";
  }
  Kokkos::finalize();
}
