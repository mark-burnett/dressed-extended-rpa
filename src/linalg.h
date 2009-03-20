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

typedef std::complex< double > complex_t;
typedef ublas::vector< complex_t > vector_t;
typedef ublas::matrix< complex_t, ublas::column_major > matrix_t;

// --------------------------------------------------------------------
// Linear algebra stuff
// --------------------------------------------------------------------
vector_t eigenvalues( const matrix_t &mat );
std::pair< vector_t, matrix_t > eig( const matrix_t &mat );

} // end namespace util

#endif // _UTIL_LINALG_H_
