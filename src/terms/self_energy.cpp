#include <cmath>
#include <cassert>

#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "linalg.h"

#include "Interaction.h"
#include "Modelspace.h"

#include "Term.h"
#include "exceptions.h"

#include "self_energy.h"

namespace terms {

util::matrix_t
self_energy( const std::vector< ParticleHoleState > &vec, double E,
             position_t pos, const PPInteraction &Gpp,
             const ParticleParticleModelspace &ppms,
             const ParticleParticleModelspace &hhms,
             const SEModelspace               &sems,
             const SingleParticleModelspace &spms ) {
    int size = vec.size();
    util::matrix_t m( size, size );
    m.clear();

    if ( ENUM_B == pos or ENUM_B_STAR == pos )
        return m;

    // self energy terms show up only on the diagonal (of course)
    for ( int i=0; i < size; ++i ) {
        ParticleHoleState ph = vec[i];

        switch ( pos ) {
            case ENUM_A:
                m( i, i ) = internal::SE_particle_line( ph, E, Gpp,
                        ppms, hhms, sems, spms )
                    + internal::SE_hole_line( ph, E, Gpp,
                        ppms, hhms, sems, spms );
                break;
            case ENUM_A_STAR:
                // A_STAR phase is always 1 on the diagonal
                m( i, i ) = internal::SE_particle_line( ph, -E, Gpp,
                        ppms, hhms, sems, spms )
                    + internal::SE_hole_line( ph, -E, Gpp,
                        ppms, hhms, sems, spms );
                break;
            default:
                throw invalid_matrix_position();
        }
    }

    return m;
}

namespace internal {

double SE_particle_line( const ParticleHoleState &ph, double E,
                         const PPInteraction &Gpp,
                         const ParticleParticleModelspace &ppms,
                         const ParticleParticleModelspace &hhms,
                         const SEModelspace &sems,
                         const SingleParticleModelspace &spms ) {
    typedef ParticleParticleState pp_t;

    // Shell indicies
    int ia = ph.ip;
    int ib = ph.ih;

    // Fragment indicies
    int iaf = ph.ipf;
    int ibf = ph.ihf;

    int Jpmin = 0;
    int Jpmax = boost::numeric_cast<int>( spms.maxj + spms.j[ia] );
    double result = 0;
    for ( int Jp = Jpmin; Jp <= Jpmax; ++Jp ) {
        double JpTerm = 0;
        // make list of outter left side states
        // loop over outter left states
//        assert( 0 != sems.ph[Jp][ia][iaf].size() );
        if ( Jp < boost::numeric_cast<int>(sems.ph.size()) ) {
            BOOST_FOREACH( pp_t left, sems.ph[Jp][ia][iaf] ) {
            //      make list of inner right side states
            //      loop over inner right side states
                int tz =
                    boost::numeric_cast<int>(spms.tz[left.ip1]
                            + spms.tz[left.ip2]);
                int parity = spms.parity[left.ip1] * spms.parity[left.ip2];
                if ( Jp >= boost::numeric_cast<int>(
                            ppms[1+tz][(parity+1)/2].size()) )
                    continue;
                BOOST_FOREACH( pp_t right, ppms[1+tz][(parity+1)/2][Jp] ) {
            //      add contribution
                    JpTerm += std::pow( Gpp( left, right ), 2 ) /
                        ( E - (   spms.pfrag[right.ip1][right.ip1f].E
                                + spms.pfrag[right.ip2][right.ip2f].E
                                - spms.hfrag[ib][ibf].E
                                - spms.hfrag[left.ip2][left.ip2f].E ) ); } } }
        //  make list of outter right side states
        //  loop over outter right states
        if ( Jp < boost::numeric_cast<int>(sems.pp.size()) ) {
            BOOST_FOREACH( pp_t left, sems.pp[Jp][ia][iaf] ) {
            //      make inner left states
            //      loop over inner left states
                int tz =
                    boost::numeric_cast<int>(spms.tz[left.ip1]
                            + spms.tz[left.ip2]);
                int parity = spms.parity[left.ip1] * spms.parity[left.ip2];
                if ( Jp >= boost::numeric_cast<int>(
                        hhms[1+tz][(parity+1)/2].size()) )
                    continue;
                BOOST_FOREACH( pp_t right, hhms[1+tz][(parity+1)/2][Jp] ) {
            //      add contribution
                    JpTerm += std::pow( Gpp( left, right ), 2 ) /
                        ( spms.pfrag[ia][iaf].E - (
                              spms.hfrag[right.ip1][right.ip1f].E
                            + spms.hfrag[right.ip2][right.ip2f].E
                            - spms.pfrag[left.ip2][left.ip2f].E ) ); } } }

        result += JpTerm * ( 2 * Jp + 1 ); }
    return result / ( 4 * spms.j[ia] + 2 ); }

double SE_hole_line    ( const ParticleHoleState &ph, double E,
                         const PPInteraction &Gpp,
                         const ParticleParticleModelspace &ppms,
                         const ParticleParticleModelspace &hhms,
                         const SEModelspace &sems,
                         const SingleParticleModelspace &spms ) {
    typedef ParticleParticleState pp_t;
    // Shell indicies
    int ia = ph.ip;
    int ib = ph.ih;

    // Fragment indicies
    int iaf = ph.ipf;
    int ibf = ph.ihf;

    int Jpmin = 0;
    int Jpmax = boost::numeric_cast<int>( spms.maxj + spms.j[ib] );
    double result = 0;
    for ( int Jp = Jpmin; Jp <= Jpmax; ++Jp ) {
        double JpTerm = 0;
        // make list of outter left side states
        // loop over outter left states
        if ( Jp < boost::numeric_cast<int>(sems.hh.size()) ) {
            BOOST_FOREACH( pp_t left, sems.hh[Jp][ib][ibf] ) {
            //      make list of inner right side states
            //      loop over inner right side states
                int tz =
                    boost::numeric_cast<int>( spms.tz[left.ip1]
                                            + spms.tz[left.ip2] );
                int parity = spms.parity[left.ip1] * spms.parity[left.ip2];
                if ( Jp >= boost::numeric_cast<int>(
                            ppms[1+tz][(parity+1)/2].size()) )
                    continue;
                BOOST_FOREACH( pp_t right, ppms[1+tz][(parity+1)/2][Jp] ) {
            //      add contribution
                    JpTerm -= std::pow( Gpp( left, right ), 2 ) /
                        ( spms.hfrag[ib][ibf].E - (
                              spms.pfrag[right.ip1][right.ip1f].E
                            + spms.pfrag[right.ip2][right.ip2f].E
                            - spms.hfrag[left.ip2][left.ip2f].E ) ); } } }
        //  make list of outter right side states
        //  loop over outter right states
//        assert( 0 != sems.hp[Jp][ib][ibf].size() );
        if ( Jp < boost::numeric_cast<int>(sems.hp.size()) ) {
            BOOST_FOREACH( pp_t left, sems.hp[Jp][ib][ibf] ) {
            //      make inner left states
            //      loop over inner left states
                int tz =
                    boost::numeric_cast<int>(spms.tz[left.ip1]
                            + spms.tz[left.ip2]);
                int parity = spms.parity[left.ip1] * spms.parity[left.ip2];
                if ( Jp >= boost::numeric_cast<int>(
                            hhms[1+tz][(parity+1)/2].size()) )
                    continue;
                BOOST_FOREACH( pp_t right, hhms[1+tz][(parity+1)/2][Jp] ) {
            //      add contribution
                    JpTerm += std::pow( Gpp( left, right ), 2 ) /
                        ( E - (   spms.pfrag[ia][iaf].E
                                + spms.pfrag[left.ip2][left.ip2f].E
                                - spms.hfrag[right.ip1][right.ip1f].E
                                - spms.hfrag[right.ip2][right.ip2f].E ) ); } } }
        result += JpTerm * ( 2 * Jp + 1 ); }
    return result / ( 4 * spms.j[ib] + 2 ); }

} // end namespace internal

Term make_self_energy( const PPInteraction &Gpp,
                       const ParticleParticleModelspace &ppms,
                       const ParticleParticleModelspace &hhms,
                       const SEModelspace &sems,
                       const SingleParticleModelspace &spms ) {
    return boost::bind( self_energy, _1, _2, _3, boost::cref(Gpp),
            boost::cref(ppms), boost::cref(hhms), boost::cref(sems),
            boost::cref(spms) );
}

} // end namespace terms
