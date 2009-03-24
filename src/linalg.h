/* String manipulation utilities.
 * 
 * Mark Burnett, November 2008
 */

#ifndef _UTIL_LINALG_H_
#define _UTIL_LINALG_H_

#include <complex>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

namespace util {

typedef ublas::vector< double > vector_t;
typedef ublas::matrix< double, ublas::column_major > matrix_t;

typedef std::complex< double > complex_t;
typedef ublas::vector< complex_t > cvector_t;
typedef ublas::matrix< complex_t, ublas::column_major > cmatrix_t;

// --------------------------------------------------------------------
// Linear algebra stuff
// --------------------------------------------------------------------
// double versions
//vector_t eigenvalues( const matrix_t &mat );
//std::pair< vector_t, matrix_t > eig( const matrix_t &mat );

// complex versions
cvector_t eigenvalues( const cmatrix_t &mat );
std::pair< cvector_t, cmatrix_t > eig( const cmatrix_t &mat );

} // end namespace util

#endif // _UTIL_LINALG_H_
