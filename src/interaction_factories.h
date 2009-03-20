#ifndef _INTERACTION_FACTORIES_H_
#define _INTERACTION_FACTORIES_H_

#include <string>

#include "Interaction.h"
#include "Modelspace.h"

PPInteraction
build_gmatrix_from_mhj_file( const std::string &filename,
                             const SingleParticleModelspace &spms );

PHInteraction
build_ph_interaction_from_pp( const PPInteraction &Gpp,
                              const SingleParticleModelspace &spms );

#endif // _INTERACTION_FACTORIES_H_
