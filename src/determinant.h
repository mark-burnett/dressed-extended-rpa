#ifndef _DETERMINANT_H_
#define _DETERMINANT_H_
/* This determinant routine was take from forum threads on the internet.
 * I found the same algorithm in at least 2 places with different names, so
 * it's almost impossible to know who gets the real credit.
 */

#include <cassert>

#include <boost/numeric/ublas/lu.hpp>

namespace ublas = boost::numeric::ublas;

template<class M>
double determinant(M const& m) {
    assert( m.size1() == m.size2() );

    // create a working copy of the input
    ublas::matrix<double> mLu(m);
    ublas::permutation_matrix<std::size_t> pivots(m.size1());

    ublas::lu_factorize(mLu, pivots);

    double det = 1.0;

    for (std::size_t i=0; i < pivots.size(); ++i) {
        if (pivots(i) != i)
            det *= -1.0;
        det *= mLu(i,i); }
    return det;
} 

#endif // _DETERMINANT_H_
