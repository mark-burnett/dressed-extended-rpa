/* Linear algebra wrappers
 * 
 * Mark Burnett, November 2008
 */

// <cassert> is required for geev.hpp, but not included in it..
#include <cassert>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include "linalg.h"

namespace lapack = boost::numeric::bindings::lapack;

namespace util {

// --------------------------------------------------------------------
// Linear algebra stuff
// --------------------------------------------------------------------
vector_t eigenvalues( const matrix_t &mat ) {
    vector_t vals(  mat.size1() );
    matrix_t right_vecs( mat.size1(), mat.size2() );
    matrix_t * dummy = 0;
    matrix_t temp(mat);
    lapack::geev( temp, vals, dummy, &right_vecs, lapack::optimal_workspace() );
    return vals;
}

std::pair< vector_t, matrix_t > eig( const matrix_t &mat ) {
    vector_t vals(  mat.size1() );
    matrix_t right_vecs( mat.size1(), mat.size2() );
    matrix_t * dummy = 0;
    matrix_t temp(mat);
    lapack::geev( temp, vals, dummy, &right_vecs, lapack::optimal_workspace() );
    return std::make_pair( vals, right_vecs );
}

} // end namespace util
