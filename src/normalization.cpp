#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "MatrixFactory.h"
#include "linalg.h"
#include "normalization.h"

namespace ublas = boost::numeric::ublas;

util::cmatrix_t
normalize_rpa_eigenvectors( const util::cmatrix_t &original_vectors ) {
    typedef ublas::vector_range < util::cvector_t > subvector_t;
    typedef ublas::matrix_column< util::cmatrix_t > column_t;

    int size = original_vectors.size2() / 2;
    util::cmatrix_t result( original_vectors );
    // loop over column vectors
    ublas::range forward_range ( 0, size );
    ublas::range backward_range( size, 2*size );
    for ( int i = 0; i < 2 * size; ++i ) {
        util::cvector_t col = column_t( result, i );
        subvector_t forward ( col, forward_range );
        subvector_t backward( col, backward_range );
        double norm = ublas::norm_2( forward ) - ublas::norm_2( backward );
        column_t( result, i ) = col / std::sqrt( norm ); }
    return result; }

    /*
util::cvector_t
make_erpa_eigenvector( const MatrixFactory &mf, double E, double epsilon ) {
    assert( E > epsilon > 0 );
    assert( epsilon < 1 );
    std::vector< std::pair< std::vector< double >, util::cvector_t > > a;
    // NOTE this does 2 evaluations.
    //      if epsilon is small enough we can skip the first
    // evaluate @ E
    //      figure out index of associated eigenvector
    // evaluate @ E - epsilon
    // calculate derivative
    // calculate norm
    // return vector
}
*/
