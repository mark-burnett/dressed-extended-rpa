#include <cmath>
#include <cassert>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "linalg.h"
#include "exceptions.h"
#include "angular_momentum.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"

#include "ladder.h"

namespace terms {

util::matrix_t
ladder( const std::vector< ParticleHoleState > &vec, double E,
        position_t pos, const PPInteraction &Gpp,
        const ParticleParticleModelspace &ppms,
        const ParticleParticleModelspace &hhms,
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
                     * internal::ladder_A_term( ph1, ph2, E, Gpp,
                               ppms, hhms, spms );
                    break;
                case ENUM_A_STAR:
                    m( k, i ) = m( i, k ) = phase * S
                     * internal::ladder_A_term( ph1, ph2, -E, Gpp,
                               ppms, hhms, spms );
                    break;
                case ENUM_B:
                    m( k, i ) = m( i, k ) = S
                     * internal::ladder_B_term( ph1, ph2, Gpp,
                               ppms, hhms, spms );
                    break;
                case ENUM_B_STAR:
                    m( k, i ) = m( i, k ) = phase * S
                     * internal::ladder_B_term( ph1, ph2, Gpp,
                               ppms, hhms, spms );
                    break;
                default:
                    throw invalid_matrix_position(); } } }
    return m;
}

namespace internal {

double ladder_A_term( const ParticleHoleState &ph1,
                      const ParticleHoleState &ph2, double E,
                      const PPInteraction &Gpp,
                      const ParticleParticleModelspace   &ppms,
                      const ParticleParticleModelspace   &hhms,
                      const SingleParticleModelspace &spms ) {
    typedef ParticleParticleState pp_t;
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
    int tz     = boost::numeric_cast<int>(spms.tz[ia] + spms.tz[id]);
    int parity = spms.parity[ia] * spms.parity[id];

    // Quick asserts
    assert( ph1.J == ph2.J );
    assert( spms.parity[ia] * spms.parity[id] ==
            spms.parity[ib] * spms.parity[ic] );
    assert( spms.tz[ia] + spms.tz[id] == spms.tz[ib] + spms.tz[ic] );

    int Jpmin = std::max( std::abs( spms.j[ia] - spms.j[id] ),
                          std::abs( spms.j[ib] - spms.j[ic] ) );
    int Jpmax = std::min( spms.j[ia] + spms.j[id], spms.j[ib] + spms.j[ic] );
    assert( Jpmax < boost::numeric_cast<int>(ppms[1+tz][(parity+1)/2].size()) );

    double result = 0;
    for ( int Jp = Jpmin; Jp <= Jpmax; ++Jp ) {
        pp_t left ( ia, id, -1, -1, Jp );
        pp_t right( ic, ib, -1, -1, Jp );
        double JpTerm = 0;
        // Intermediate terms above Fermi surface
//        assert( 0 != ppms[1 + tz][(parity+1)/2][Jp].size() );
        BOOST_FOREACH(pp_t i_pp, ppms[1 + tz][(parity+1)/2][Jp]) {
            assert( parity == spms.parity[i_pp.ip1]*spms.parity[i_pp.ip2] );
            double Si_pp = spms.pfrag[i_pp.ip1][i_pp.ip1f].S
                         * spms.pfrag[i_pp.ip2][i_pp.ip2f].S;
            JpTerm += Si_pp * Gpp( left, i_pp ) * Gpp( i_pp, right ) /
                ( E - ( - spms.hfrag[ib][ibf].E - spms.hfrag[id][idf].E
                        + spms.pfrag[i_pp.ip1][i_pp.ip1f].E
                        + spms.pfrag[i_pp.ip2][i_pp.ip2f].E ) );
        }
        // Intermediate terms below Fermi surface
//        assert( 0 != hhms[1 + tz][(parity+1)/2][Jp].size() );
        BOOST_FOREACH(pp_t i_hh, hhms[1 + tz][(parity+1)/2][Jp]) {
            assert( parity == spms.parity[i_hh.ip1]*spms.parity[i_hh.ip2] );
            double Si_hh = spms.hfrag[i_hh.ip1][i_hh.ip1f].S
                         * spms.hfrag[i_hh.ip2][i_hh.ip2f].S;
            JpTerm += Si_hh * Gpp( left, i_hh ) * Gpp( i_hh, right ) /
                ( E - ( spms.pfrag[ia][iaf].E + spms.pfrag[ic][icf].E
                        - spms.hfrag[i_hh.ip1][i_hh.ip1f].E
                        - spms.hfrag[i_hh.ip2][i_hh.ip2f].E ) );
        }
        result -= JpTerm * (  2 * Jp + 1 )
                * wigner6j( spms.j[ia], spms.j[ib], J,
                            spms.j[ic], spms.j[id], Jp );
    }
    return result;
}

double ladder_B_term( const ParticleHoleState &ph1,
                      const ParticleHoleState &ph2,
                      const PPInteraction &Gpp,
                      const ParticleParticleModelspace   &ppms,
                      const ParticleParticleModelspace   &hhms,
                      const SingleParticleModelspace &spms ) {
    typedef ParticleParticleState pp_t;
    // Shell indicies
    int ia = ph1.ip;
    int ib = ph1.ih;
    int ic = ph2.ih;
    int id = ph2.ip;

    // Fragment indicies
    int iaf = ph1.ipf;
    int ibf = ph1.ihf;
    int icf = ph2.ihf;
    int idf = ph2.ipf;

    // Useful items
    int J      = ph1.J;
    int tz     = boost::numeric_cast<int>(spms.tz[ia] + spms.tz[id]);
    int parity = spms.parity[ia] * spms.parity[id];

    // Quick asserts
    assert( ph1.J == ph2.J );
    assert( spms.parity[ia] * spms.parity[id] ==
            spms.parity[ib] * spms.parity[ic] );
    assert( spms.tz[ia] + spms.tz[id] == spms.tz[ib] + spms.tz[ic] );

    int Jpmin = std::max( std::abs( spms.j[ia] - spms.j[id] ),
                          std::abs( spms.j[ib] - spms.j[ic] ) );
    int Jpmax = std::min( spms.j[ia] + spms.j[id], spms.j[ib] + spms.j[ic] );
    assert( Jpmax < boost::numeric_cast<int>(ppms[1+tz][(parity+1)/2].size()) );

    double result = 0;
    for ( int Jp = Jpmin; Jp <= Jpmax; ++Jp ) {
        pp_t left ( ia, id, -1, -1, Jp );
        pp_t right( ic, ib, -1, -1, Jp );
        double JpTerm = 0;
        // Intermediate terms above Fermi surface
        BOOST_FOREACH(pp_t i_pp, ppms[1 + tz][(parity+1)/2][Jp]) {
            assert( parity == spms.parity[i_pp.ip1]*spms.parity[i_pp.ip2] );
            double Si_pp = spms.pfrag[i_pp.ip1][i_pp.ip1f].S
                         * spms.pfrag[i_pp.ip2][i_pp.ip2f].S;
            JpTerm += Si_pp * Gpp( left, i_pp ) * Gpp( i_pp, right ) /
                - ( - spms.hfrag[ib][ibf].E - spms.hfrag[ic][icf].E
                    + spms.pfrag[i_pp.ip1][i_pp.ip1f].E
                    + spms.pfrag[i_pp.ip2][i_pp.ip2f].E );
        }
        // Intermediate terms below Fermi surface
        BOOST_FOREACH(pp_t i_hh, hhms[1 + tz][(parity+1)/2][Jp]) {
            assert( parity == spms.parity[i_hh.ip1]*spms.parity[i_hh.ip2] );
            double Si_hh = spms.hfrag[i_hh.ip1][i_hh.ip1f].S
                         * spms.hfrag[i_hh.ip2][i_hh.ip2f].S;
            JpTerm += Si_hh * Gpp( left, i_hh ) * Gpp( i_hh, right ) /
                - (   spms.pfrag[ia][iaf].E + spms.pfrag[id][idf].E
                    - spms.hfrag[i_hh.ip1][i_hh.ip1f].E
                    - spms.hfrag[i_hh.ip2][i_hh.ip2f].E );
        }
        result -= JpTerm * (  2 * Jp + 1 )
                * wigner6j( spms.j[ia], spms.j[ib], J,
                            spms.j[ic], spms.j[id], Jp );
    }
    return result;
}

} // end namespace internal

Term make_ladder( const PPInteraction &Gpp,
                  const ParticleParticleModelspace &ppms,
                  const ParticleParticleModelspace &hhms,
                  const SingleParticleModelspace &spms ) {
    return boost::bind( ladder, _1, _2, _3, boost::cref(Gpp),
            boost::cref(ppms), boost::cref(hhms), boost::cref(spms) );
}

} // end namespace terms
