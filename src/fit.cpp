#include <vector>
#include <iostream>

#include <boost/math/special_functions/pow.hpp>

#include "fit.h"

namespace bm = boost::math;

boost::tuple< fit_tuple_t, fit_tuple_t >
find_grad( const fitness_f &f, const fit_tuple_t &current,
           const double x2, const double delta ) {
    // Declare our first and 2nd derivatives
    int psize = current.get<0>().size();
    int nsize = current.get<1>().size();
    vector_t pgrad( psize );
    vector_t ngrad( nsize );
    vector_t plapl( psize );
    vector_t nlapl( nsize );

    // Protons
    for ( int i = 0; i < psize; ++i ) {
        fit_tuple_t plus( current );
        fit_tuple_t minus( current );
        plus.get<0>()[i]  += delta;
        minus.get<0>()[i] -= delta;
        double x2plus  = f( plus );
        double x2minus = f( minus );
        pgrad[i] = ( x2plus - x2minus ) / ( 2 * delta );
        plapl[i] = ( x2plus + x2minus - 2 * x2 ) / bm::pow<2>(delta); }
    // Neutrons
    for ( int i = 0; i < nsize; ++i ) {
        fit_tuple_t plus( current );
        fit_tuple_t minus( current );
        plus.get<1>()[i]  += delta;
        minus.get<1>()[i] -= delta;
        double x2plus  = f( plus );
        double x2minus = f( minus );
        ngrad[i] = ( x2plus - x2minus ) / ( 2 * delta );
        nlapl[i] = ( x2plus + x2minus - 2 * x2 ) / bm::pow<2>(delta); }

    return boost::make_tuple( fit_tuple_t( pgrad, ngrad ),
                              fit_tuple_t( plapl, nlapl ) ); }

vector_t next_vec( const vector_t &current, const vector_t &grad,
                   const vector_t &lapl ) {
    vector_t next( current.size() );
//    std::cout << "next_vec" << std::endl;
//    std::cout << current.size() << " " << grad.size() << " "
//        << lapl.size() << std::endl;
    for ( unsigned int i = 0; i < current.size(); ++i ) {
        next[i] = current[i] - grad[i] / lapl[i]; }
    return next; }

fit_tuple_t find_next( const fit_tuple_t &current, const fit_tuple_t &grad,
                       const fit_tuple_t &lapl ) {
    return boost::make_tuple(
            next_vec( current.get<0>(), grad.get<0>(), lapl.get<0>() ),
            next_vec( current.get<1>(), grad.get<1>(), lapl.get<1>() ) ); }

boost::tuple< fit_tuple_t, double >
optimize( const fitness_f &f, const vector_t &pi, const vector_t &ni,
          const double delta, const double error, const int max_iter ) {
    // initialize
    fit_tuple_t current( pi, ni );
    fit_tuple_t grad;
    fit_tuple_t lapl;
    double x2 = -1;

    // Loop
    for ( int i = 0; i < max_iter; ++i ) {
        // Get the chi-squared for this point
        x2 = f( current );
        std::cout << i << " " << x2 << std::endl;
        // Are we good enough?
        if ( x2 < error )
            break;
        if ( max_iter - 1 == i ) {
            throw "errorz"; }

        // evaluate gradient
        boost::tie( grad, lapl ) = find_grad( f, current, x2, delta );
        // set next guess
        current = find_next( current, grad, lapl ); }

    return boost::make_tuple( current, x2 ); }
