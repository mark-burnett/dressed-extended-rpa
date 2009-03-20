#ifndef _RPA_TERMS_SCREENING_H_
#define _RPA_TERMS_SCREENING_H_

#include "linalg.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

namespace terms {

util::matrix_t
screening( const std::vector< ParticleHoleState > &vec, double E,
           position_t pos, const PHInteraction &Gph,
           const SingleParticleModelspace &spms );

namespace internal {
double screening_A_term( const ParticleHoleState &ph1,
                         const ParticleHoleState &ph2, double E,
                         const PHInteraction &Gph,
                         const ParticleHoleModelspace   &phms,
                         const SingleParticleModelspace &spms );
double screening_B_term( const ParticleHoleState &ph1,
                         const ParticleHoleState &ph2,
                         const PHInteraction &Gph,
                         const ParticleHoleModelspace   &phms,
                         const SingleParticleModelspace &spms );
} // end namespace internal

Term make_screening( const PHInteraction &Gph,
                     const ParticleHoleModelspace   &phms,
                     const SingleParticleModelspace &spms );

} // end namespace terms

#endif // _RPA_TERMS_SCREENING_H_
