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

    int     num_solutions =  ab_right.get<1>() - ab_left.get<1>();
    assert( num_solutions == ab_left.get<0>()  - ab_right.get<0>() );

    return num_solutions; }

// Root finding functions
double base_root_function( double E, const MatrixFactory &mf, int index ) {
    std::vector< double > vals = util::sorted_eigenvalues( mf.build(E) );
    return vals[index] - E; }

double root_find_solution( const MatrixFactory &mf, const interval_t &region,
                           const std::vector< double > lower_vals,
                           const std::vector< double > upper_vals,
                           double epsilon ) {
    int index = std::upper_bound( lower_vals.begin(), lower_vals.end(),
            region.lower() ) - lower_vals.begin();
    double flower = lower_vals[index] - region.lower();
    double fupper = upper_vals[index] - region.upper();
    return util::false_position(
            boost::bind( base_root_function, _1, boost::cref(mf),
                index ),
            region.lower(), region.upper(), flower, fupper, epsilon ); }

// This is the main search algorithm for the (D)ERPA.
std::vector< double >
solve_region( const MatrixFactory &mf, const interval_t &region,
              const std::vector< double > lower_vals,
              const std::vector< double > upper_vals,
              double epsilon ) {
    // Determine # solutions
    int num_solutions = get_num_solutions( lower_vals, region.lower(),
                                           upper_vals, region.upper() );
    // If no solutions, return empty.
    if ( 0 == num_solutions ) {
        return std::vector< double >(); }

    // If the interval has no width, but has solutions, return the value.
    if ( std::abs( region.upper() - region.lower() ) < epsilon ) {
        return std::vector< double >( 1, region.lower() ); }

    // If 1 solution, root_find.
    if ( 1 == num_solutions ) {
        std::vector< double > temp( 1,
                root_find_solution( mf, region, lower_vals, upper_vals,
                    epsilon ) );
        return temp; }

    // If > 1 solution, sub-divide region.
    double center = boost::numeric::median( region );
    std::vector< double > center_vals
        = util::sorted_eigenvalues( mf.build( center ) );
    std::vector< double > solutions
        = solve_region( mf, interval_t( region.lower(), center ),
                        lower_vals, center_vals, epsilon );
    { std::vector< double > upper_solutions
        = solve_region( mf, interval_t( center, region.upper() ),
                        center_vals, upper_vals, epsilon );
        solutions.insert( solutions.end(), upper_solutions.begin(),
                upper_solutions.end() ); }
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
        // Evaluate eigenvalues at upper and lower limits.
        std::vector< double > lower_vals
            = util::sorted_eigenvalues( mf.build( region.lower() ) );
        std::vector< double > upper_vals
            = util::sorted_eigenvalues( mf.build( region.upper() ) );
        // Solve inside region
        std::vector< double > region_results = solve_region( mf, region,
                lower_vals, upper_vals, epsilon );
        results.insert( results.end(), region_results.begin(),
                                       region_results.end() );
        lower = asymptotes[a]; }
    return results;
}
