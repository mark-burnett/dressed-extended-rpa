#include <iostream>
#include <cmath>
#include <boost/function.hpp>

#include "exceptions.h"
#include "find_root.h"

namespace util {

// The basic method.
double false_position( const boost::function< double (double) > &f,
                       double left_limit, double right_limit,
                       double precision, int max_iter ) {
    double a = left_limit;
    double b = right_limit;

    double fa = f(a);
    double fb = f(b);

    // fa and fb must have opposite signs to ensure a solution.
    if ( fa * fb > precision ) {
        std::cout << "false_position failure: a = " << a << ", b = " << b
            << ", fa = " << fa << ", fb = " << fb << std::endl;
        throw root_finding_error();
    }

    double c      = a - 1;
    double c_last = a - 1;

    for (int iter = 0; iter < max_iter; ++iter) {
        c = (fb * a - fa * b) / (fb - fa);
        double fc = f(c);

        // Verify both that the answer is close to zero,
        // and that we made a small last step (so it's more stable).
        if ( std::abs(fc) < precision && std::abs(c-c_last) < precision )
            break;

        // Pick a side
        if ( fa * fc < 0 ) {
            b  = c;
            fb = fc;
        } else {
            a  = c;
            fa = fc;
        }
        c_last = c;
    }

    return c;
}

// Method with initial values 
double false_position( const boost::function< double (double) > &f,
                       double left_limit, double right_limit,
                       double fa, double fb,
                       double precision, int max_iter ) {
    double a = left_limit;
    double b = right_limit;

    // fa and fb must have opposite signs to ensure a solution.
    if ( fa * fb > precision ) {
        std::cout << "false_position failure: a = " << a << ", b = " << b
            << ", fa = " << fa << ", fb = " << fb << std::endl;
        throw root_finding_error();
    }

    double c      = a - 1;
    double c_last = a - 1;

    for (int iter = 0; iter < max_iter; ++iter) {
        c = (fb * a - fa * b) / (fb - fa);
        double fc = f(c);
//        std::cout << "f(" << a << ") = " << fa << ", f(" << b << ") = "
//            << fb << ", f(" << c << ") = " << fc << std::endl;

        // Verify both that the answer is close to zero,
        // and that we made a small last step (so it's more stable).
        if ( std::abs(fc) < precision && std::abs(c-c_last) < precision )
            break;

        // Pick a side
        if ( fa * fc < 0 ) {
            b  = c;
            fb = fc;
        } else {
            a  = c;
            fa = fc;
        }
        c_last = c;
    }

    return c;
}

} // end namespace util
