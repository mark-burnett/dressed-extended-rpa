#ifndef _MATRIX_FACTORY_H_
#define _MATRIX_FACTORY_H_

#include <vector>

#include "linalg.h"
#include "Modelspace.h"
#include "Term.h"

class MatrixFactory {
    public:
        MatrixFactory( const util::matrix_t                   &nstatic_matrix,
                       const std::vector< Term >              &ndynamic_terms,
                       const SingleParticleModelspace         &nspms,
                       const std::vector< ParticleHoleState > &nph_states,
                       int nJ, int nparity, int ntz )
            : static_matrix( nstatic_matrix ), dynamic_terms( ndynamic_terms ),
              ph_states( nph_states ), spms( nspms ),
              J( nJ ), parity( nparity ), tz( ntz ) {
                  assert( J >= 0 );
                  assert( 1 == parity || -1 == parity );
                  assert( 1 >= tz && -1 <= tz );
                  assert( static_matrix.size1() == static_matrix.size2() );
                  assert( 2 * ph_states.size()  == static_matrix.size1() );
              }
        util::matrix_t build( double E ) const;
    private:
        const util::matrix_t                   static_matrix;
        const std::vector< Term >              dynamic_terms;
        const std::vector< ParticleHoleState > ph_states;
        const SingleParticleModelspace         spms;
        int J, parity, tz;
};

util::matrix_t
build_static_rpa_matrix( const std::vector< Term > &terms,
                         const std::vector< ParticleHoleState > &ph_states );
util::matrix_t
build_static_erpa_matrix( const std::vector< Term > &static_terms,
                          const std::vector< Term > &dynamic_terms,
                          const std::vector< ParticleHoleState > &ph_states );

#endif // _MATRIX_FACTORY_H_
