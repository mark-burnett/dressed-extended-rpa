#ifndef _NUCLEAR_PANDYA_H_
#define _NUCLEAR_PANDYA_H_

#include "Modelspace.h"
#include "Interaction.h"

double pandya( const PPInteraction &G_pp,
               const SingleParticleModelspace &spms,
               const ParticleHoleState &A,
               const ParticleHoleState &B );

#endif // _NUCLEAR_PANDYA_H_
