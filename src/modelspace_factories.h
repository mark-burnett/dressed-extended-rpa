#ifndef _MODELSPACE_FACTORIES_H_
#define _MODELSPACE_FACTORIES_H_

#include <istream>

#include "Modelspace.h"

SingleParticleModelspace
read_sp_modelspace_from_file( const std::string &filename );

ParticleParticleModelspace
build_pp_modelspace_from_sp( const SingleParticleModelspace &spms );

ParticleHoleModelspace
build_ph_modelspace_from_sp( const SingleParticleModelspace &spms );

ParticleHoleModelspace
build_ph_shells_from_sp( const SingleParticleModelspace &spms );

#endif // _MODELSPACE_FACTORIES_H_
