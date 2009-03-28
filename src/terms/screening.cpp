#include <cassert>
#include <cmath>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "linalg.h"
#include "angular_momentum.h"
#include "exceptions.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

#include "screening.h"

namespace terms {

util::matrix_t
screening( const std::vector< ParticleHoleState > &vec, double E,
           position_t pos, const PHInteraction &Gph,
           const ParticleHoleModelspace   &phms,
           const SingleParticleModelspace &spms ) {
    int size = vec.size();
    util::matrix_t m( size, size );
    m.clear();

    for ( int i = 0; i < size; ++i ) {
        ParticleHoleState ph1 = vec[i];
        double Sa = spms.pfrag[ph1.ip][ph1.ipf].S
                  * spms.hfrag[ph1.ih][ph1.ihf].S;
        for ( int k = i; k < size; ++k ) {
            ParticleHoleState ph2 = vec[k];
            double S = Sa * spms.pfrag[ph2.ip][ph2.ipf].S
                          * spms.hfrag[ph2.ih][ph2.ihf].S;

            // Phase for A* and B*
            int phase = std::pow( -1.0, spms.j[ ph1.ip ] + spms.j[ ph1.ih ]
                                      + spms.j[ ph2.ip ] + spms.j[ ph2.ih ] );
            // NOTE: Assuming real valued, so A and B are symmetric
            //   (A is normally Hermitian)
            switch ( pos ) {
                case ENUM_A:
                    m( k, i ) = m( i, k ) = S
                      * internal::screening_A_term( ph1, ph2, E,
                                Gph, phms, spms );
                    break;
                case ENUM_A_STAR:
                    m( k, i ) = m( i, k ) = phase * S
                      * internal::screening_A_term( ph1, ph2, -E,
                                Gph, phms, spms );
                    break;
                case ENUM_B:
                    m( k, i ) = m( i, k ) = S
                      * internal::screening_B_term( ph1, ph2, Gph, phms, spms );
                    break;
                case ENUM_B_STAR:
                    m( k, i ) = m( i, k ) = phase * S
                      * internal::screening_B_term( ph1, ph2, Gph, phms, spms );
                    break;
                default:
                    throw invalid_matrix_position(); } } }
    return m;
}

namespace internal {

double screening_A_term( const ParticleHoleState &ph1,
                         const ParticleHoleState &ph2, double E,
                         const PHInteraction &Gph,
                         const ParticleHoleModelspace   &phms,
                         const SingleParticleModelspace &spms ) {
    typedef ParticleHoleState ph_t;

    // Shell indicies
    int ia = ph1.ip;
    int ib = ph1.ih;
    int ic = ph2.ip;
    int id = ph2.ih;

    // Fragment indicies
    int iaf = ph1.ipf;
    int ibf = ph1.ihf;
    int icf = ph2.ipf;
    int idf = ph2.ihf;

    // Useful items
    int J      = ph1.J;
    int tz     = boost::numeric_cast<int>(spms.tz[ia] - spms.tz[ic]);
    int parity = spms.parity[ia] * spms.parity[ic];

    // Quick asserts
    assert( ph1.J == ph2.J );
    assert( spms.parity[ia] * spms.parity[ic] ==
            spms.parity[ib] * spms.parity[id] );
    assert( spms.tz[ia] - spms.tz[ic] == spms.tz[ib] - spms.tz[id] );
    
    int Jpmin = std::max( std::abs( spms.j[ia] - spms.j[ic] ),
                          std::abs( spms.j[ib] - spms.j[id] ) );
    int Jpmax = std::min( spms.j[ia] + spms.j[ic], spms.j[ib] + spms.j[id] );
    assert( Jpmax < boost::numeric_cast<int>(phms[1+tz][(parity+1)/2].size()) );
    assert( Jpmax < boost::numeric_cast<int>(phms[1-tz][(parity+1)/2].size()) );

    double result = 0;
    for ( int Jp = Jpmin; Jp <= Jpmax; ++Jp ) {
        ph_t left(  ia, ic, -1, -1, Jp );
        ph_t right( ib, id, -1, -1, Jp );
        double JpTerm = 0;
        // Forward going terms
        BOOST_FOREACH(ph_t i_ph, phms[1 + tz][(parity+1)/2][Jp]) {
            double Si_ph = spms.pfrag[i_ph.ip][i_ph.ipf].S
                         * spms.hfrag[i_ph.ih][i_ph.ihf].S;
            JpTerm += Si_ph *
                Gph( left, i_ph ) * Gph( i_ph, right ) /
                ( E - ( spms.pfrag[ic][icf].E - spms.hfrag[ib][ibf].E
                        + spms.pfrag[ i_ph.ip ][ i_ph.ipf ].E
                        - spms.pfrag[ i_ph.ih ][ i_ph.ihf ].E ) );
        }
        // Backward going terms
        BOOST_FOREACH(ph_t i_ph, phms[1 - tz][(parity+1)/2][Jp] ) {
            double Si_ph = spms.pfrag[i_ph.ip][i_ph.ipf].S
                         * spms.hfrag[i_ph.ih][i_ph.ihf].S;
            ph_t r_ph( i_ph.ih, i_ph.ip, -1, -1, Jp );
            JpTerm += Si_ph *
                Gph( left, r_ph ) * Gph( r_ph, right ) /
                ( E - ( spms.pfrag[ia][iaf].E - spms.hfrag[id][idf].E
                        + spms.pfrag[ i_ph.ip ][ i_ph.ipf ].E
                        - spms.pfrag[ i_ph.ih ][ i_ph.ihf ].E ) );
        }
        result -= JpTerm * std::pow( -1.0, spms.j[ib] + spms.j[ic] + J + Jp )
                * (  2 * Jp + 1 )
                * wigner6j( spms.j[ia], spms.j[ib], J,
                            spms.j[id], spms.j[ic], Jp );
    }
    return result;
}

double screening_B_term( const ParticleHoleState &ph1,
                         const ParticleHoleState &ph2,
                         const PHInteraction &Gph,
                         const ParticleHoleModelspace   &phms,
                         const SingleParticleModelspace &spms ) {
    typedef ParticleHoleState ph_t;

    // Shell indicies
    int ia = ph1.ip;
    int ib = ph1.ih;
    int ic = ph2.ih;
    int id = ph2.ip;

    // Fragment indicies
    // NOTE: c is now a hole fragment, and d a particle fragment
    int iaf = ph1.ipf;
    int ibf = ph1.ihf;
    int icf = ph2.ihf;
    int idf = ph2.ipf;

    // Useful items
    int J      = ph1.J;
    int tz     = boost::numeric_cast<int>(spms.tz[ia] - spms.tz[ic]);
    int parity = spms.parity[ia] * spms.parity[ic];

    // Quick asserts
    assert( ph1.J == ph2.J );
    assert( spms.parity[ia] * spms.parity[ic] ==
            spms.parity[ib] * spms.parity[id] );
    assert( spms.tz[ia] - spms.tz[ic] == spms.tz[ib] -  spms.tz[id] );

    int Jpmin = std::max( std::abs( spms.j[ia] - spms.j[ic] ),
                          std::abs( spms.j[ib] - spms.j[id] ) );
    int Jpmax = std::min( spms.j[ia] + spms.j[ic], spms.j[ib] + spms.j[id] );
    assert( Jpmax < boost::numeric_cast<int>(phms[1+tz][(parity+1)/2].size()) );

    double result = 0;
    for ( int Jp = Jpmin; Jp <= Jpmax; ++Jp ) {
        ph_t left(  ia, ic, -1, -1, Jp );
        ph_t right( ib, id, -1, -1, Jp );
        double JpTerm = 0;
        BOOST_FOREACH(ph_t i_ph, phms[1 + tz][(parity+1)/2][Jp]) {
            double Si_ph = spms.pfrag[i_ph.ip][i_ph.ipf].S
                         * spms.hfrag[i_ph.ih][i_ph.ihf].S;
            // Forward going terms
            JpTerm += Si_ph *
                Gph( left, i_ph ) * Gph( i_ph, right ) /
                - ( spms.pfrag[id][idf].E - spms.hfrag[ib][ibf].E
                    + spms.pfrag[ i_ph.ip ][ i_ph.ipf ].E
                    - spms.pfrag[ i_ph.ih ][ i_ph.ihf ].E );
        }
        BOOST_FOREACH(ph_t i_ph, phms[1 - tz][(parity+1)/2][Jp] ) {
            double Si_ph = spms.pfrag[i_ph.ip][i_ph.ipf].S
                         * spms.hfrag[i_ph.ih][i_ph.ihf].S;
            // Backward going terms
            ph_t r_ph( i_ph.ih, i_ph.ip, -1, -1, Jp );
            JpTerm += Si_ph *
                Gph( left, r_ph ) * Gph( r_ph, right ) /
                - ( spms.pfrag[ia][iaf].E - spms.hfrag[ic][icf].E
                    + spms.pfrag[ i_ph.ip ][ i_ph.ipf ].E
                    - spms.pfrag[ i_ph.ih ][ i_ph.ihf ].E );
        }
        result -= JpTerm * std::pow( -1.0, spms.j[ib] + spms.j[ic] + J + Jp )
                * (  2 * Jp + 1 )
                * wigner6j( spms.j[ia], spms.j[ic], J,
                            spms.j[id], spms.j[ib], Jp );
    }
    return result;
}

} // end namespace internal

Term make_screening( const PHInteraction &Gph,
                     const ParticleHoleModelspace   &phms,
                     const SingleParticleModelspace &spms ) {
    return boost::bind( screening, _1, _2, _3, boost::cref(Gph),
            boost::cref(phms), boost::cref(spms) );
}

} // end namespace terms
