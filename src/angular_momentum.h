/* Contains routines related to angular momentum algebra, including:
 * Clebsch-Gordan coefficients, Wigner 6-j symbols, etc.
 *
 * Mark Burnett, November 2008
 */

#ifndef _NUCLEAR_ANGULAR_MOMENTUM_H_
#define _NUCLEAR_ANGULAR_MOMENTUM_H_

double factorial(int n);

bool is_integer            (double j);
bool is_half_integer       (double j);
bool is_strict_half_integer(double j);
bool is_triangular         (double j1, double j2, double j3);


double wigner3j(double j1, double j2, double j3,
                double m1, double m2, double m3);
double clebsch_gordan(double j1, double m1, double j2, double m2,
                      double J, double M);
double wigner6j(double j1, double j2, double j3,
                double J1, double J2, double J3);

#endif // _NUCLEAR_ANGULAR_MOMENTUM_H_
