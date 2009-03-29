#include <cassert>
#include <vector>
#include <algorithm>

#include <boost/numeric/conversion/cast.hpp>

#include "intervals.h"

bool intervals_equal( const interval_t &a, const interval_t &b ) {
    return a.upper() == b.upper() && a.lower() == b.lower(); }

// Builds solution intervals by interleaving elements of A and B.
// Both vectors must be sorted, and A[i] < B[i].
std::vector< interval_t >
build_root_intervals( const std::vector< double > &A,
                      const std::vector< double > &B ) {
    assert( A.size() == B.size() );
    std::vector< interval_t > intervals;
    for ( int i = 0; i < boost::numeric_cast<int>(A.size()); ++i ) {
//        std::cout << A[i] << " " << B[i] << std::endl;
        intervals.push_back( interval_t( A[i], B[i] ) ); }
    return intervals;
}

std::vector< interval_t >
get_sub_intervals( const std::vector< interval_t > &root_intervals ) {
    std::vector< double > end_points;
    for ( int r = 0; r < boost::numeric_cast<int>(root_intervals.size()); ++r ){
        end_points.push_back( root_intervals[r].lower() );
        end_points.push_back( root_intervals[r].upper() ); }

    std::sort( end_points.begin(), end_points.end() );

    std::vector< interval_t > sub_intervals;
    for ( int s = 0; s < boost::numeric_cast<int>(end_points.size()) - 1; ++s ){
        sub_intervals.push_back( interval_t( end_points[s], end_points[s+1] ));}
    return sub_intervals; }

std::vector< double >
get_sub_interval_probabilities( const std::vector< interval_t > &root_intervals,
                                const std::vector< interval_t > &sub_intervals){
    using boost::numeric::width;
    using boost::numeric::subset;

    std::vector< double > root_densities;
    for ( int r = 0; r < boost::numeric_cast<int>(root_intervals.size()); ++r ){
        root_densities.push_back( 1/width(root_intervals[r]) ); }

    std::vector< double > probabilities( sub_intervals.size(), 0 );
    for ( int s = 0; s < boost::numeric_cast<int>(sub_intervals.size()); ++s ) {
        for ( int r = 0; r < boost::numeric_cast<int>(root_intervals.size());
                ++r ) {
            if ( subset( sub_intervals[s], root_intervals[r] ) ) {
                probabilities[s] += root_densities[r]
                                  * width( sub_intervals[s] ); } } }
//    std::sort( probabilities.begin(), probabilities.end() );
    return probabilities;
}
