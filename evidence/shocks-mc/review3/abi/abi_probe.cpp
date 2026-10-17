#include <cstddef>
#include <iostream>
#include "old_encounter.hpp"
#include "shocks/encounter.hpp"
struct OldHeader {
    char magic[16];
    int version;
    int enabled;
    int nrdata;
    int nidata;
    int nlocal;
    int nschedules;
    std::int64_t next_tag;
  };
struct NewHeader {
  char magic[16];
  int version;
  int enabled;
  int nrdata;
  int nidata;
  int nlocal;
  int nschedules;
  std::int64_t next_tag;
};
static_assert(sizeof(OldHeader)==sizeof(NewHeader));
static_assert(offsetof(OldHeader,magic)==offsetof(NewHeader,magic));
static_assert(offsetof(OldHeader,version)==offsetof(NewHeader,version));
static_assert(offsetof(OldHeader,enabled)==offsetof(NewHeader,enabled));
static_assert(offsetof(OldHeader,nrdata)==offsetof(NewHeader,nrdata));
static_assert(offsetof(OldHeader,nidata)==offsetof(NewHeader,nidata));
static_assert(offsetof(OldHeader,nlocal)==offsetof(NewHeader,nlocal));
static_assert(offsetof(OldHeader,nschedules)==offsetof(NewHeader,nschedules));
static_assert(offsetof(OldHeader,next_tag)==offsetof(NewHeader,next_tag));
static_assert(sizeof(old_shocks::EncounterState)-sizeof(shocks::EncounterState)==sizeof(shocks::Observation));
static_assert(sizeof(old_shocks::Event)==sizeof(shocks::Event));
static_assert(offsetof(old_shocks::EncounterState,particle_id)==offsetof(shocks::EncounterState,particle_id));
static_assert(offsetof(old_shocks::EncounterState,sequence)==offsetof(shocks::EncounterState,sequence));
static_assert(offsetof(old_shocks::EncounterState,last_generation)==offsetof(shocks::EncounterState,last_generation));
static_assert(offsetof(old_shocks::EncounterState,image)==offsetof(shocks::EncounterState,image));
static_assert(offsetof(old_shocks::EncounterState,phase)==offsetof(shocks::EncounterState,phase));
static_assert(offsetof(old_shocks::EncounterState,clean)==offsetof(shocks::EncounterState,clean));
static_assert(offsetof(old_shocks::EncounterState,suppressing_reentry)==offsetof(shocks::EncounterState,suppressing_reentry));
static_assert(offsetof(old_shocks::EncounterState,suppressed_reentries)==offsetof(shocks::EncounterState,suppressed_reentries));
static_assert(offsetof(old_shocks::EncounterState,flags)==offsetof(shocks::EncounterState,flags));
static_assert(offsetof(old_shocks::EncounterState,previous)==offsetof(shocks::EncounterState,previous));
static_assert(offsetof(old_shocks::EncounterState,first)==offsetof(shocks::EncounterState,first));
static_assert(offsetof(old_shocks::EncounterState,last)==offsetof(shocks::EncounterState,last));
static_assert(offsetof(old_shocks::EncounterState,post)==offsetof(shocks::EncounterState,post));
int main() { std::cout << sizeof(NewHeader) << " " << sizeof(shocks::Observation) << " " << sizeof(old_shocks::EncounterState) << " " << sizeof(shocks::EncounterState) << " " << sizeof(shocks::Event) << "\n"; }
