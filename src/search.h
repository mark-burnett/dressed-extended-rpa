#ifndef _SEARCH_H_
#define _SEARCH_H_

#include <vector>
#include <boost/tuple/tuple.hpp>

#include "MatrixFactory.h"

boost::tuple< int, int >
split_values( const std::vector< double > &sorted_vals, double E );

int get_num_solutions( const std::vector< double > &left_vals,  double left,
                       const std::vector< double > &right_vals, double right );

std::vector< double >
solve_derpa_eigenvalues( double Emax,
                         const MatrixFactory &mf,
                         const std::vector< double > &asymptotes,
                         double epsilon = 0.01 );
#endif // _SEARCH_H_
