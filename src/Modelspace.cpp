#include <cassert>
#include <iostream>
#include <set>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/foreach.hpp>
#include "Modelspace.h"

// --------------------------------------------------------------------
// SP Modelspace helpers
// --------------------------------------------------------------------

int get_max_pp_J( const SingleParticleModelspace &spms, int tz, int parity ) {
    int J = -1;
    for ( int a = 0; a < spms.size; ++a ) {
        for ( int b = a; b < spms.size;  ++b ) {
            if ( parity != spms.parity[a] * spms.parity[b] ||
                 tz     != spms.tz[a]     + spms.tz[b]     )
                continue;
            int nJ = spms.j[a] + spms.j[b];
            if ( nJ > J )
                J = nJ; } }
    return J; }

int get_max_ph_J( const SingleParticleModelspace &spms, int tz, int parity ) {
    int J = -1;
    for ( int p = 0; p < spms.size; ++p ) {
        for ( int h = 0; h < spms.size; ++h ) {
            // Check for right parity and isospin
            if ( parity != spms.parity[p] * spms.parity[h] ||
                 (  tz   != spms.tz[p]    - spms.tz[h] &&
                   -tz   != spms.tz[p]    - spms.tz[h] ) )
                continue;
            // Check whether we have actual particle and hole fragments here
            if ( spms.pfrag[p].size() < 1 ||
                 spms.hfrag[h].size() < 1 )
                continue;
            int nJ = spms.j[p] + spms.j[h];
            if ( nJ > J )
                J = nJ; } }
    return J; }

void print_sp_state( std::ostream &o, int i,
                     const SingleParticleModelspace &spms ) {
    o << "(" << spms.j[i] << " " << spms.parity[i] << " " << spms.n[i]
        << " " << spms.tz[i] << ")"; }

// --------------------------------------------------------------------
// PP Modelspace helpers
// --------------------------------------------------------------------
double pp_energy( const ParticleParticleState &pp,
                  const SingleParticleModelspace &spms ) {
    return spms.pfrag[pp.ip1][pp.ip1f].E + spms.pfrag[pp.ip2][pp.ip2f].E; }

double hh_energy( const ParticleParticleState &pp,
                  const SingleParticleModelspace &spms ) {
    return spms.hfrag[pp.ip1][pp.ip1f].E + spms.hfrag[pp.ip2][pp.ip2f].E; }

int parity( const ParticleParticleState &pp,
            const SingleParticleModelspace &spms ) {
    return spms.parity[pp.ip1] * spms.parity[pp.ip2]; }

// --------------------------------------------------------------------
// PH Modelspace helpers
// --------------------------------------------------------------------
// Returns the energy of a particle hole state (fragment)
double energy( const ParticleHoleState &ph,
               const SingleParticleModelspace &spms ) {
    return spms.pfrag[ph.ip][ph.ipf].E - spms.hfrag[ph.ih][ph.ihf].E; }

int parity( const ParticleHoleState &ph,
            const SingleParticleModelspace &spms ) {
    return spms.parity[ph.ip] * spms.parity[ph.ih]; }

std::vector< double > get_erpa_asymptotes( int tz, int parity, int J,
                                const ParticleParticleModelspace &ppms,
                                const ParticleParticleModelspace &hhms,
                                const SingleParticleModelspace &spms ) {
    std::set< double > result;
    for ( int pptz = -1; pptz <= 1; ++pptz ) {
        int hhtz = pptz - tz;
        for ( int ppparity = -1; ppparity <= 1; ppparity += 2 ) {
            int hhparity = parity * ppparity;
            int ppJmax = ppms[pptz+1][(ppparity+1)/2].size() - 1;
            for ( int ppJ = 0; ppJ <= ppJmax; ++ppJ ) {
                const std::vector< ParticleParticleState > &pp_states
                    = ppms[1+pptz][(ppparity+1)/2][ppJ];
                int hhJabsmax = boost::numeric_cast<int>(
                           hhms[1+hhtz][(hhparity+1)/2].size() - 1);
                int hhJmin = std::min( std::abs( ppJ - J ), hhJabsmax );
                int hhJmax = std::min( ppJ + J, hhJabsmax );
                for ( int hhJ = hhJmin; hhJ <= hhJmax; ++hhJ ) {
                    const std::vector< ParticleParticleState > &hh_states
                        = hhms[1+hhtz][(hhparity+1)/2][hhJ];
                    BOOST_FOREACH( const ParticleParticleState &pp,
                            pp_states ) {
                        BOOST_FOREACH( const ParticleParticleState &hh,
                                hh_states ) {
                            result.insert( pp_energy(pp, spms)
                                         - hh_energy(hh, spms) ); } } } } } }
    return std::vector< double >( result.begin(), result.end() ); }

void print_ph_modelspace_sizes( std::ostream &o,
                                const ParticleHoleModelspace &phms ) {
    o << " tz  J^parity  size " << std::endl;
    for ( int tz = -1; tz <= 1; ++tz ) {
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            for ( int J = 0; J <
                    boost::numeric_cast<int>(phms[tz+1][(parity+1)/2].size());
                        ++J) { 
                o << tz << " " << J;
                if ( parity > 0 )
                    o << "+ ";
                else
                    o << "- ";
                o << phms[tz+1][(parity+1)/2][J].size() << std::endl; } } } }

void print_ph_modelspace_sizes( std::ostream &o, int tz,
                                const ParticleHoleModelspace &phms ) {
    o << " tz  J^parity  size " << std::endl;
    for ( int parity = -1; parity <= 1; parity += 2 ) {
        for ( int J = 0; J <
                boost::numeric_cast<int>(phms[tz+1][(parity+1)/2].size());
                    ++J) { 
            o << tz << " " << J;
            if ( parity > 0 )
                o << "+ ";
            else
                o << "- ";
            o << phms[tz+1][(parity+1)/2][J].size() << std::endl; } } }

void print_ph_state( std::ostream &o,
                     const ParticleHoleState        &ph,
                     const SingleParticleModelspace &spms ) {
    o << "(";
    print_sp_state( o, ph.ip, spms );
    o << " ";
    print_sp_state( o, ph.ih, spms );
    o << " " << ph.J << ")";
}
