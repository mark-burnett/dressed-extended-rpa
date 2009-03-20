#ifndef _RPA_TERMS_NON_INTERACTING_H_
#define _RPA_TERMS_NON_INTERACTING_H_

#include "linalg.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

namespace terms {

util::matrix_t
non_interacting( const std::vector< ParticleHoleState > &vec, double E,
                 position_t pos, const SingleParticleModelspace &spms );

Term make_non_interacting( const SingleParticleModelspace &spms );

} // end namespace terms

#endif // _RPA_TERMS_NON_INTERACTING_H_
