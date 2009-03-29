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

// Modelspace used only in self-energy terms
PPFromSP
build_ppsp_modelspace_from_sp( const SingleParticleModelspace &spms );

PPFromSP
build_hhsp_modelspace_from_sp( const SingleParticleModelspace &spms );

PPFromSP
build_phsp_modelspace_from_sp( const SingleParticleModelspace &spms );

PPFromSP
build_hpsp_modelspace_from_sp( const SingleParticleModelspace &spms );

SEModelspace
build_se_modelspace_from_sp( const SingleParticleModelspace &spms );

#endif // _MODELSPACE_FACTORIES_H_
