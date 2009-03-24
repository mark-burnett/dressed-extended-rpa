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
cvector_t eigenvalues( const cmatrix_t &mat ) {
    cvector_t vals(  mat.size1() );
    cmatrix_t right_vecs( mat.size1(), mat.size2() );
    cmatrix_t * dummy = 0;
    cmatrix_t temp(mat);
    lapack::geev( temp, vals, dummy, &right_vecs, lapack::optimal_workspace() );
    return vals;
}

std::pair< cvector_t, cmatrix_t > eig( const cmatrix_t &mat ) {
    cvector_t vals(  mat.size1() );
    cmatrix_t right_vecs( mat.size1(), mat.size2() );
    cmatrix_t * dummy = 0;
    cmatrix_t temp(mat);
    lapack::geev( temp, vals, dummy, &right_vecs, lapack::optimal_workspace() );
    return std::make_pair( vals, right_vecs );
}

} // end namespace util
