#include <cmath>

#include <boost/bind.hpp>

#include "linalg.h"

#include "Modelspace.h"
#include "Interaction.h"

#include "Term.h"
#include "first_order.h"

#include "exceptions.h"

namespace terms {

util::matrix_t
first_order( const std::vector< ParticleHoleState > &vec, double E,
             position_t pos, const PHInteraction &Gph,
             const SingleParticleModelspace &spms ) {
    unsigned int size = vec.size();
    util::matrix_t m( size, size );
    m.clear();

    for ( unsigned int i=0; i < size; ++i ) {
        const ParticleHoleState &ph1 = vec[i];
        double Sa = spms.pfrag[ ph1.ip ][ ph1.ipf ].S
                  * spms.hfrag[ ph1.ih ][ ph1.ihf ].S;
        for ( unsigned int k=i; k < size; ++k ) {
            const ParticleHoleState &ph2 = vec[k];
            ParticleHoleState r_ph2( ph2.ih, ph2.ip, ph2.ihf, ph2.ipf, ph2.J );

            double S = Sa * spms.pfrag[ ph2.ip ][ ph2.ipf ].S
                          * spms.hfrag[ ph2.ih ][ ph2.ihf ].S;

            double value;
            int phase = 1;
            switch ( pos ) {
                case ENUM_A_STAR:
                case ENUM_A: // A and A* are hermitian, but we pretend sym.
                    value = S * Gph( ph1, ph2 );
                    m( i, k ) = value;
                    m( k, i ) = value;
                    break;
                case ENUM_B_STAR:
                case ENUM_B: // B and B* are symmetric
                    phase *= std::pow( -1.0, spms.j[ph2.ip] - spms.j[ph2.ih]
                                           + ph2.J );
                    value = S * Gph( ph1, r_ph2 );
                    m( i, k ) = phase * value;
                    m( k, i ) = phase * value;
                    break;
                default:
                    throw invalid_matrix_position(); } } }
    return m;
    // Dummy code for unused E
    ++E;
}

Term make_first_order( const PHInteraction &Gph,
                       const SingleParticleModelspace &spms ) {
    return boost::bind( first_order, _1, _2, _3,
            boost::cref(Gph), boost::cref(spms) );
}

} // end namespace terms
