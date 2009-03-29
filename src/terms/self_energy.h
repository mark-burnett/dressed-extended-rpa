#ifndef _RPA_TERMS_SELF_ENERGY_H_
#define _RPA_TERMS_SELF_ENERGY_H_

#include <boost/foreach.hpp>
#include <boost/bind.hpp>

#include "linalg.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

namespace terms {

util::matrix_t
self_energy( const std::vector< ParticleHoleState > &vec, double E,
             position_t pos, const PPInteraction &Gpp,
             const ParticleParticleModelspace &ppms,
             const ParticleParticleModelspace &hhms,
             const SEModelspace &sems,
             const SingleParticleModelspace &spms );


namespace internal {
double SE_particle_line( const ParticleHoleState &ph, double E,
                         const PPInteraction &Gpp,
                         const ParticleParticleModelspace &ppms,
                         const ParticleParticleModelspace &hhms,
                         const SEModelspace &sems,
                         const SingleParticleModelspace &spms );
double SE_hole_line    ( const ParticleHoleState &ph, double E,
                         const PPInteraction &Gpp,
                         const ParticleParticleModelspace &ppms,
                         const ParticleParticleModelspace &hhms,
                         const SEModelspace &sems,
                         const SingleParticleModelspace &spms );
} // end namespace internal

Term make_self_energy( const PPInteraction &Gpp,
                       const ParticleParticleModelspace &ppms,
                       const ParticleParticleModelspace &hhms,
                       const SEModelspace &sems,
                       const SingleParticleModelspace &spms );

} // end namespace terms

#endif // _RPA_TERMS_SELF_ENERGY_H_
