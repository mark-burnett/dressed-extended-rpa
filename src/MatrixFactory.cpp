#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "linalg.h"
#include "Term.h"
#include "MatrixFactory.h"

util::matrix_t
MatrixFactory::build( double E ) const {
    // Copy the static elements
    util::matrix_t result( static_matrix );

    // Generate the dynamic elements
    int size = result.size1() / 2;
    ublas::range first_half( 0, size );
    ublas::range second_half( size, 0 );

    typedef ublas::matrix_range< util::matrix_t > submatrix_t;
    submatrix_t A     ( result, first_half, first_half );
    submatrix_t A_star( result, second_half, second_half );

    BOOST_FOREACH( const Term &t, dynamic_terms ) {
        A      += t( ph_states, E, ENUM_A );
        A_star -= t( ph_states, E, ENUM_A_STAR );
    }
    return result;
}
