#include <boost/foreach.hpp>
#include <cassert>
#include <vector>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/interval/io.hpp>

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

std::vector< interval_t >
make_root_intervals( const std::vector< double > &left_vals,  double left,
                     const std::vector< double > &right_vals, double right,
                     const interval_t &region, int num_solutions ) {
    boost::tuple< int, int > ab_left  = split_values( left_vals,  left );
    boost::tuple< int, int > ab_right = split_values( right_vals, right );

    std::vector< double > lower_bounds( num_solutions ),
                          upper_bounds( num_solutions );
    for ( int i = 0; i < num_solutions; ++i ) {
        double lval = right_vals[ ab_right.get<1>() - num_solutions + i ];
        double uval = left_vals [i + ab_left.get<1>()];
//        std::cout << lval << ", " << uval << std::endl;
        lower_bounds[i] = std::max( lval, region.lower() );
        upper_bounds[i] = std::min( uval, region.upper() ); }
    return build_root_intervals( lower_bounds, upper_bounds); }

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
    // If the region is a singleton, leave.
    if ( boost::numeric::singleton( region ) ) {
        std::cout << "region is singleton.  " << region << std::endl;
        return std::vector< double >();
    }
    // Find eigenvalues at left and right boundaries.
    std::vector< double > lower_vals
        = util::sorted_eigenvalues( mf.build( region.lower() ) );
    std::vector< double > upper_vals
        = util::sorted_eigenvalues( mf.build( region.upper() ) );

    // Get basic information about the region.
    int num_solutions = get_num_solutions( lower_vals, region.lower(),
                                           upper_vals, region.upper() );
    std::cout << "(" << region.lower() << ", " << region.upper() << ") -> "
        << num_solutions << std::endl;

    // Return an empty vector if there are no solutions in the region:
    if ( 0 == num_solutions )
        return std::vector< double >();

    // Generate the root intervals
    std::vector< interval_t > root_intervals
        = make_root_intervals( lower_vals, region.lower(),
                               upper_vals, region.upper(),
                               region, num_solutions );

    // Return the answer if there is only one solution in the region.
    // NOTE: We generated the root intervals first, because the root interval
    // is guaranteed to no larger than the region size.
    // This may not be optimal, since we could provide the root finding
    // algorithm with the y-values at both ends.
    if ( 1 == num_solutions ) {
        boost::function< double(double) > f = boost::bind( base_root_function,
                _1, boost::cref(mf) );
        std::vector< double > result( 1,
                util::false_position( f, root_intervals[0].lower(),
                                         root_intervals[0].upper() ) );
        std::cout << "solution found: " << result[0] << std::endl;
        return result; }

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
            std::vector< double > sub_results;
            // Check for the case that we could not narrow down the
            // root finding interval.
            bool bisected = false;
            if ( boost::numeric::equal( sub_intervals[imax], region ) ) {
                bisected = true;
                double lower  = sub_intervals[imax].lower();
                double upper  = sub_intervals[imax].upper();
                double center = ( lower + upper ) / 2;
                std::cout << "Bisecting (" << lower << ", " << center
                    << ", " << upper << ")." << std::endl;
                sub_results = solve_region( mf, interval_t( lower, center ) );
                if ( 0 == sub_results.size() )
                    sub_results = solve_region( mf, interval_t( center, upper));
            } else {
                std::cout << "Solving normally: ("
                    << sub_intervals[imax].lower() << ", "
                    << sub_intervals[imax].upper() << ")" << std::endl;
                sub_results = solve_region( mf, sub_intervals[imax] );
            }
            // If no solution was found, set the probability to 0 and move on.
            if ( 0 == sub_results.size() )
                prob[imax] = 0;
            else {
//                if ( root_intervals.size() < sub_results.size() ) {
//                    if ( bisected )
//                        std::cout << "we bisected." << std::endl;
//                    std::cout << i << "/" << num_solutions <<  " "
//                        << sub_results.size() << std::endl;
//                }
                assert( root_intervals.size() >= sub_results.size() );
                if ( sub_results.size() > 1 ) {
                    for ( int sr = 0;
                            sr < boost::numeric_cast<int>(sub_results.size());
                            ++sr ) {
                        std::cout << "SR: " <<sub_results[sr] << std::endl;
                    }
                }
                // loop over sub results
                for ( int sr = 0;
                        sr < boost::numeric_cast<int>(sub_results.size());
                        ++sr ) {
                    // add the solution to our solution set.
                    double solution = sub_results[sr];
                    solutions.push_back( solution );
//                    std::cout << "solution found: " << sub_results[sr]
//                        << std::endl;
                    // identify root interval it belongs to
                    double epsilon = 0.000001;
                    int root_position
                        = get_num_solutions( lower_vals, region.lower(),
                           util::sorted_eigenvalues(mf.build(solution
                                   - epsilon)),
                           solution - epsilon);
                    std::cout << "rp = " << root_position << std::endl;
                    // delete that root interval
                    assert( 0 != root_intervals.size() );
                    // FIXME don't erase, replace with a tiny interval
                    root_intervals[root_position]
                        = interval_t( solution ); } //, solution ); }
//                    root_intervals.erase( root_intervals.begin()
//                                        + root_position ); }
                // end this loop (go back to top loop)
                break; } } }
    // If these conditions are not met, we have missed a solution
    assert( num_solutions == boost::numeric_cast<int>(solutions.size()) );
//    assert( root_intervals.empty() );
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
        std::cout << "   (" << region.lower() << ", " << region.upper() << ")"
            << std::endl;
        std::vector< double > region_results = solve_region( mf, region );
        results.insert( results.end(), region_results.begin(),
                                       region_results.end() );
        lower = asymptotes[a]; }
    return results;
}
