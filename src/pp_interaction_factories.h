#ifndef _PP_INTERACTION_FACTORIES_H_
#define _PP_INTERACTION_FACTORIES_H_

#include <string>

#include "Interaction.h"
#include "Modelspace.h"

PPInteraction
build_gmatrix_from_mhj_file( const std::string &filename,
                             const SingleParticleModelspace &spms );

#endif // _PP_INTERACTION_FACTORIES_H_
