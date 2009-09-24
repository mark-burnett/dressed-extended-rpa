#include <gtest/gtest.h>

#include <vector>

#include "fit.h"

#include <boost/bind.hpp>
#include <boost/math/special_functions/pow.hpp>
namespace bm = boost::math;

double myfit( const fit_tuple_t &guess, const fit_tuple_t &data ) {
    double x2 = 0;
    for ( unsigned int i = 0; i < guess.get<0>().size();  ++i ) {
        x2 += bm::pow<4>( guess.get<0>()[i] - data.get<0>()[i] ); }
    for ( unsigned int i = 0; i < guess.get<1>().size();  ++i ) {
        x2 += bm::pow<4>( guess.get<1>()[i] - data.get<1>()[i] ); }
    return x2; }

void print_fit( const fit_tuple_t &x ) {
    std::cout << "Protons:\n";
    for ( unsigned int i = 0; i < x.get<0>().size();  ++i ) {
        std::cout << x.get<0>()[i] << "\n"; }
    std::cout << "Neutrons:\n";
    for ( unsigned int i = 0; i < x.get<1>().size();  ++i ) {
        std::cout << x.get<1>()[i] << "\n"; }
    std::cout.flush(); }

TEST( Fit, QuadraticSurface ) {
    // Setup target to fit.
    vector_t pt(3);
    vector_t nt(3);
    pt[0] = 3.1; pt[1] = -1.6; pt[2] = 4.72;
    nt[0] = 1.7; nt[1] = -0.5; nt[2] = 0.29;
    fit_tuple_t target( pt, nt );
    
    fitness_f f = boost::bind( myfit, _1, target );

    // Create starting point
    pt[0] = 10;  pt[1] = -2; pt[2] = 1;
    nt[0] = 0.2; nt[1] = 3;  nt[2] = 0;

//    std::cout << "Start point:\n";
//    print_fit( fit_tuple_t( pt, nt ) );
//
//    std::cout << "Target:\n";
//    print_fit( target );

    fit_tuple_t answer;
    double x2;
    boost::tie( answer, x2 ) = optimize( f, pt, nt, 0.001, 1e-10, 20 );

//    std::cout << "End point:\n";
//    print_fit( answer );
}
