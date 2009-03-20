#ifndef _RPA_TERMS_LADDER_H_
#define _RPA_TERMS_LADDER_H_

#include "linalg.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

namespace terms {

util::matrix_t
ladder( const std::vector< ParticleHoleState > &vec, double E,
        position_t pos, const PPInteraction &Gpp,
        const ParticleParticleModelspace &ppms,
        const ParticleParticleModelspace &hhms,
        const SingleParticleModelspace &spms );

namespace internal {
double ladder_A_term( const ParticleHoleState &ph1,
                      const ParticleHoleState &ph2, double E,
                      const PPInteraction &Gpp,
                      const ParticleParticleModelspace &ppms,
                      const ParticleParticleModelspace &hhms,
                      const SingleParticleModelspace &spms );
double ladder_B_term( const ParticleHoleState &ph1,
                      const ParticleHoleState &ph2,
                      const PPInteraction &Gpp,
                      const ParticleParticleModelspace &ppms,
                      const ParticleParticleModelspace &hhms,
                      const SingleParticleModelspace &spms );
} // end namespace internal

Term make_ladder( const PPInteraction &Gpp,
                  const ParticleParticleModelspace &ppms,
                  const ParticleParticleModelspace &hhms,
                  const SingleParticleModelspace &spms );

} // end namespace terms

#endif // _RPA_TERMS_LADDER_H_
