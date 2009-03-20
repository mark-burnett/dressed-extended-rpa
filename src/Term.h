#ifndef _DYANMIC_TERM_H_
#define _DYANMIC_TERM_H_

#include <boost/function.hpp>

#include "linalg.h"
#include "Modelspace.h"

enum position_t { ENUM_A, ENUM_A_STAR, ENUM_B, ENUM_B_STAR };

typedef boost::function<
    util::matrix_t( const std::vector< ParticleHoleState >, double, position_t )
> Term;

#endif // _DYANMIC_TERM_H_
