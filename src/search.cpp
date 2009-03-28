#include <cassert>
#include <vector>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>

#include "find_root.h"
#include "determinant.h"
#include "linalg.h"
#include "intervals.h"
#include "MatrixFactory.h"

// Returns the number of elements (above, below) E in vals.
boost::tuple< int, int >
split_values( const std::vector< double > &sorted_vals, double E ) {
    int num_above = sorted_vals.end()
                  - std::upper_bound( sorted_vals.begin(),
                                      sorted_vals.end(), E );
    int num_below = sorted_vals.size() - num_above;

    if ( E == sorted_vals[ sorted_vals.size() - num_above - 1 ] )
        --num_below;
    return boost::make_tuple( num_above, num_below ); }

// Returns the number of self consistent solutions between left and right
// given the sorted eigenvalues for the problem evaluated at left and right.
int get_num_solutions( const std::vector< double > &left_vals,  double left,
                       const std::vector< double > &right_vals, double right ) {
    boost::tuple< int, int > ab_left  = split_values( left_vals,  left );
    boost::tuple< int, int > ab_right = split_values( right_vals, right );

    int     num_solutions =  std::abs( ab_right.get<0>() - ab_left.get<0>() );
    assert( num_solutions == std::abs( ab_left.get<1>()  - ab_right.get<1>() ));

    return num_solutions; }

double base_root_function( double E, const MatrixFactory &mf ) {
    util::matrix_t m = mf.build(E);
    m -= E * ublas::identity_matrix< double >( m.size1() );
    return determinant( m ); }

// This is the main search algorithm for the (D)ERPA.
// Here is a rough outline of the algorithm:
// Evaluate the eigenvalue problem at the boundaries of the region.
// From those eigenvalues,
//      determine the number of solutions in the region.
//      construct root intervals containing each solution
// Using those root intervals,
//      construct sub intervals
//      determine probabilities for the sub intervals
// Check the sub intervals most likely to contain a solution (recursively)
std::vector< double >
solve_region( const MatrixFactory &mf, const interval_t &region ) {
    // Find eigenvalues at left and right boundaries.
    std::vector< double > lower_vals
        = util::sorted_eigenvalues( mf.build( region.lower() ) );
    std::vector< double > upper_vals
        = util::sorted_eigenvalues( mf.build( region.upper() ) );

    // Get basic information about the region.
    int num_solutions = get_num_solutions( lower_vals, region.lower(),
                                           upper_vals, region.upper() );

    // Return an empty vector if there are no solutions in the region:
    if ( 0 == num_solutions )
        return std::vector< double >();

    // Generate the root intervals
    std::vector< interval_t > root_intervals
        = build_root_intervals( lower_vals, upper_vals );

    // Return the answer if there is only one solution in the region.
    // NOTE: We generated the root intervals first, because the root interval
    // is guaranteed to be smaller than the region size.
    // This may not be optimal, since we could provide the root finding
    // algorithm with the y-values at both ends.
    if ( 1 == num_solutions ) {
        boost::function< double(double) > f = boost::bind( base_root_function,
                _1, boost::cref(mf) );
        return std::vector< double >( 1,
                util::false_position( f, root_intervals[0].lower(),
                                         root_intervals[0].upper() ) ); }

    // Loop until we get have a solution for every root interval
    std::vector< double > solutions;
    // We are going to modify root_intervals in the loop, so let's save the size
    for ( int i = 0; i < num_solutions; ++i ) {
        // We have already found all the solutions.
        if ( 0 == root_intervals.size() )
            break;
        // Generate the sub intervals and the associated probabilities
        std::vector< interval_t > sub_intervals
            = get_sub_intervals( root_intervals );
        std::vector< double > prob
            = get_sub_interval_probabilities( root_intervals, sub_intervals );
        // Go through the sub intervals
        for ( int s = 0; s < boost::numeric_cast<int>(sub_intervals.size());
                ++s ) {
            // Locate the sub interval most likely to have a solution.
            int imax = std::max_element( prob.begin(), prob.end() )
                     - prob.begin();
            // Don't bother looking if there is nothing there
            assert( 0 != prob[imax] );
            // Search that location
            std::vector< double > sub_results =
                solve_region( mf, sub_intervals[imax] );
            // If no solution was found, set the probability to 0 and move on.
            if ( 0 == sub_results.size() )
                prob[imax] = 0;
            else {
                // loop over sub results
                for ( int sr = 0;
                        sr < boost::numeric_cast<int>(sub_results.size());
                        ++sr ) {
                    // add the solution to our solution set.
                    solutions.push_back( sub_results[sr] );
                    // identify root interval it belongs to
                    int root_position
                        = get_num_solutions( lower_vals, region.lower(),
                           util::sorted_eigenvalues(mf.build(sub_results[sr])),
                           sub_results[sr] );
                    // delete that root interval
                    root_intervals.erase( root_intervals.begin()
                                        + root_position ); }
                // end this loop (go back to top loop)
                break; } } }
    // If these conditions are not met, we have missed a solution
    assert( num_solutions == boost::numeric_cast<int>(solutions.size()) );
    assert( root_intervals.empty() );
    return solutions; }

// Finds all (D)ERPA solutions up to the next asymptote above Emax.
// epsilon defines how far away from asymptotes to evaluate the problem.
std::vector< double >
solve_derpa_eigenvalues( double Emax,
                         const MatrixFactory &mf,
                         const std::vector< double > &asymptotes,
                         double epsilon ) {
    double lower = 0;
    std::vector< double > results;
    for ( int a = 0; a < boost::numeric_cast<int>(asymptotes.size()); ++a ) {
        interval_t region( lower + epsilon, asymptotes[a] - epsilon );
        // All done.
        if ( lower > Emax )
            break;
        std::vector< double > region_results = solve_region( mf, region );
        results.insert( results.end(), region_results.begin(),
                                       region_results.end() );
        lower = asymptotes[a]; }
    return results;
}
