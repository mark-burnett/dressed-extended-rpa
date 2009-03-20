#ifndef _RPA_TERMS_STATIC_H_
#define _RPA_TERMS_STATIC_H_

#include "linalg.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

namespace terms {

util::matrix_t
first_order( const std::vector< ParticleHoleState > &vec, double E,
             position_t pos, const PHInteraction &Gph,
             const SingleParticleModelspace &spms );

Term make_first_order( const PHInteraction &Gph );

} // end namespace terms

#endif // _RPA_TERMS_STATIC_H_
