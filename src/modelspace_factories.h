#ifndef _MODELSPACE_FACTORIES_H_
#define _MODELSPACE_FACTORIES_H_

#include <istream>

#include "Modelspace.h"

SingleParticleModelspace
read_sp_modelspace_from_file( const std::string &filename );

// General 2p/ph modelspaces
ParticleParticleModelspace
build_pp_modelspace_from_sp( const SingleParticleModelspace &spms );

ParticleParticleModelspace
build_hh_modelspace_from_sp( const SingleParticleModelspace &spms );

ParticleHoleModelspace
build_ph_modelspace_from_sp( const SingleParticleModelspace &spms );

ParticleHoleModelspace
build_ph_shells_from_sp( const SingleParticleModelspace &spms );

// Self-energy only modelspaces
PPFromSPModelspace
build_ppsp_modelspace_from_sp( const SingleParticleModelspace &spms );

PPFromSPModelspace
build_hhsp_modelspace_from_sp( const SingleParticleModelspace &spms );

#endif // _MODELSPACE_FACTORIES_H_
