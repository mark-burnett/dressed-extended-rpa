/* Function implementations for angular momentum algebra.
 *
 * Mark Burnett, November 2008
 */

#include <vector>
#include <cmath>

#include "exceptions.h"
//#include "state.h"

#include "angular_momentum.h"

// NOTE this factorial is not threadsafe.
static std::vector< double > factorial_array(2,1); //[1, 1]);
double factorial(int n)
{
    for (int i = factorial_array.size(); i <= n; ++i) {
        size_t size = factorial_array.size();
        factorial_array.push_back( size *
                factorial_array[ size - 1 ] );
    }

    return factorial_array[n];
}

// --------------------------------------------------------------------
// Angular momentum utility functions.
// --------------------------------------------------------------------

// Returns true if j is an integer.
bool
is_integer(double j) {
    return static_cast<double>(static_cast<int>(j)) == j;
}

// Returns true if j is an integer or half integer (0, 0.5, 1, 1.5, ...):
bool
is_half_integer(double j) {
    return static_cast<double>(static_cast<int>(2*j))/2 == j;
}

// Returns true if j is a half integer (0.5, 1.5, 2.5, ...):
bool
is_strict_half_integer(double j) {
    return static_cast<double>( static_cast<int>(j + 0.5) ) == j + 0.5;
}

// Returns true if these angular momenta are triangular.
bool
is_triangular(double ja, double jb, double jc) {
//    return (ja + jb > jc)
//        && (std::abs( ja - jb ) < jc)
//        && ( static_cast<double>(static_cast<int>(ja+jb+jc)) ==
//                  (ja + jb + jc) );
    // The third condition is only necessary for 6j symbols.
    return !(  (std::abs(ja - jb) > jc) ||
               (         ja + jb  < jc) ||
              !( static_cast<double>(static_cast<int>(ja+jb+jc)) ==
                  (ja + jb + jc) ) );
}

// --------------------------------------------------------------------
// Clebsch-Gordan and Wigner-nJ symbols
// --------------------------------------------------------------------

// A cute number, that comes up in a few places...worth pulling out.
double
triangle_coef(double ja, double jb, double jc) {
    return static_cast<double>( factorial( ja + jb - jc)     *
                                factorial( ja - jb + jc)     *
                                factorial(-ja + jb + jc) )   /
           static_cast<double>( factorial( ja + jb + jc + 1) );
}

// Calculates the wigner 3j symbol:
//      (j1, j2, j3)
//      (m1, m2, m3)
double wigner3j( double j1, double j2, double j3,
                 double m1, double m2, double m3  ) {
    // Error checking:

    // Make sure arguments are half integers
    if ( !is_half_integer(j1) || !is_half_integer(j2) || !is_half_integer(j3) ||
         !is_half_integer(m1) || !is_half_integer(m2) || !is_half_integer(m3) )
        throw illegal_angular_momentum();

    // Make sure m's work with given j's
    if ( std::abs(m1) > j1    || std::abs(m2) > j2    || std::abs(m3) > j3    ||
         !is_integer(j1 - m1) || !is_integer(j2 - m2) || !is_integer(j3 - m3) )
        throw illegal_angular_momentum();

    // Angular momentum projections must sum to 0
    if ( 0 != static_cast<int>(m1 + m2 + m3) )
        throw illegal_angular_momentum();

    // Triangular j's
    if ( !is_triangular(j1, j2, j3) )
        return 0;

    // Actual calculation:

    int tmin = static_cast<int>( std::max( std::max(
                    j2 - j3 - m1, j1 - j3 + m2), 0.0 )  );
    int tmax = static_cast<int>( std::min( std::min(
                    j1 - m1, j2 + m2), j1 + j2 - j3  )  ) + 2;

    double sum = 0.0;
    for (int t = tmin; t < tmax; ++t) {
        // FIXME make pretty and with an if statement like wigner6j
        int denom = factorial(t) *
                factorial( m1 - j2 + j3 + t ) *
                factorial(-j1 - m2 + j3 + t ) *
                factorial( j1 + j2 - j3 - t ) *
                factorial( j1 - m1 - t      ) *
                factorial( j2 + m2 - t      ) ;
        if (denom > 0)
            sum += std::pow(-1.0, t) / static_cast<double>(denom);
    }

    // Simple factor out front
    double factor = std::pow(-1.0, j1 - j2 - m3) * std::sqrt(
            triangle_coef(j1, j2, j3) *
            factorial(j1 + m1) * factorial(j1 - m1) *
            factorial(j2 + m2) * factorial(j2 - m2) *
            factorial(j3 + m3) * factorial(j3 - m3) );

    return factor * sum;
}

// Calculates the Clebsch Gordan coefficient:
//  (j1, m1; j2, m2 | J, M)
double
clebsch_gordan(double j1, double m1, double j2, double m2, double J, double M) {
    return std::pow(-1.0, j1 - j2 + M) * std::sqrt( 2*J + 1) *
        wigner3j(j1, j2, J, m1, m2, -M);
}

// Calculates the Wigner 6j symbol:
//      / j1, j2, j3 \   ---
//      \ J1, J2, J3 /
double
wigner6j(double j1, double j2, double j3, double J1, double J2, double J3) {
    // Error checking:

    // Angular momenta are half integers:
    if (!is_half_integer(j1) || !is_half_integer(j2) || !is_half_integer(j3) ||
        !is_half_integer(J1) || !is_half_integer(J2) || !is_half_integer(J3) )
        throw illegal_angular_momentum();

    // Check for non-triangular angular momenta (but return 0 silently):
    if (!is_triangular(j1, j2, j3) || !is_triangular(j1, J2, J3) ||
        !is_triangular(J1, j2, J3) || !is_triangular(J1, J2, j3) )
        return 0.0;

    int tmin = static_cast<int>(std::max( std::max( std::max( std::max(
                j1 + j2 + j3, j1 + J2 + J3), J1 + j2 + J3),
                J1 + J2 + j3), -1.0));
    int tmax = static_cast<int>(std::min( std::min(
                j1 + j2 + J1 + J2, j2 + j3 + J2 + J3), j1 + j3 + J1 + J3)) + 2;

    double sum = 0.0;
    for (int t = tmin; t < tmax; ++t) {
        // check for negative factorial arguments:
        if (!( t - j1 - j2 - j3 < 0 || j1 + j2 + J1 + J2 - t < 0 || 
               t - j1 - J2 - J3 < 0 || j2 + j3 + J2 + J3 - t < 0 || 
               t - J1 - j2 - J3 < 0 || j1 + j3 + J1 + J3 - t < 0 || 
               t - J1 - J2 - j3 < 0 ))
            sum += std::pow(-1.0, t) * factorial(t + 1) / (
                factorial(t - j1 - j2 - j3) * factorial(j1 + j2 + J1 + J2 - t) *
                factorial(t - j1 - J2 - J3) * factorial(j2 + j3 + J2 + J3 - t) *
                factorial(t - J1 - j2 - J3) * factorial(j1 + j3 + J1 + J3 - t) *
                factorial(t - J1 - J2 - j3) );
    }

    double factor = std::sqrt( triangle_coef(j1, j2, j3) *
                               triangle_coef(j1, J2, J3) *
                               triangle_coef(J1, j2, J3) *
                               triangle_coef(J1, J2, j3) );
    return factor * sum;
}
