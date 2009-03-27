#ifndef _INTERVAL_H_
#define _INTERVAL_H_

#include <vector>
#include <boost/numeric/interval.hpp>

typedef boost::numeric::interval< double > interval_t;

std::vector< interval_t >
build_root_intervals( const std::vector< double > &A,
                      const std::vector< double > &B );

std::vector< interval_t >
get_sub_intervals( const std::vector< interval_t > &root_intervals );

std::vector< double >
get_sub_interval_probabilities( const std::vector< interval_t > &root_intervals,
                                const std::vector< interval_t > &sub_intervals);

#endif // _INTERVAL_H_
