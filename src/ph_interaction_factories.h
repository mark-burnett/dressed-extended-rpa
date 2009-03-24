#ifndef _PH_INTERACTION_FACTORIES_H_
#define _PH_INTERACTION_FACTORIES_H_

#include "Interaction.h"
#include "Modelspace.h"

PHInteraction
build_ph_interaction_from_pp( const PPInteraction &Gpp,
                              const SingleParticleModelspace &spms );

#endif // _PH_INTERACTION_FACTORIES_H_
