#ifndef _ERPA_FIT_H_
#define _ERPA_FIT_H_

#include <vector>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>

typedef std::vector< double > vector_t;

typedef boost::tuple< vector_t, vector_t > fit_tuple_t;

typedef boost::function< double ( fit_tuple_t ) > fitness_f;


boost::tuple< fit_tuple_t, double >
optimize( const fitness_f &f, const vector_t &pi, const vector_t &ni,
          const double delta, const double error, const int max_iter = 20 );

#endif // _ERPA_FIT_H_
