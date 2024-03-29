#ifndef _UTIL_FIND_ROOT_H_
#define _UTIL_FIND_ROOT_H_

#include <boost/function.hpp>

namespace util {

double false_position( const boost::function< double (double) > &f,
                       double left_limit, double right_limit,
                       double precision=0.000001, int max_iter=20 );

double false_position( const boost::function< double (double) > &f,
                       double left_limit, double right_limit,
                       double fa, double fb,
                       double precision=0.000001, int max_iter=20 );

} // end namespace util

#endif // _UTIL_FIND_ROOT_H_
