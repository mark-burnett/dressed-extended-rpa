#include <iostream>
#include <boost/numeric/conversion/cast.hpp>
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
                J = nJ;
        }
    }
    return J;
}

int get_max_ph_J( const SingleParticleModelspace &spms, int tz, int parity ) {
    int J = -1;
    for ( int p = 0; p < spms.size; ++p ) {
        for ( int h = 0; h < spms.size; ++h ) {
            // Check for right parity and isospin
            if ( parity != spms.parity[p] * spms.parity[h] ||
                 tz     != spms.tz[p]     - spms.tz[h]     )
                continue;
            // Check whether we have actual particle and hole fragments here
            if ( spms.pfrag[p].size() < 1 ||
                 spms.hfrag[h].size() < 1 )
                continue;
            int nJ = spms.j[p] + spms.j[h];
            if ( nJ > J )
                J = nJ; } }
    return J;
}

void print_sp_state( std::ostream &o, int i,
                     const SingleParticleModelspace &spms ) {
    o << "(" << spms.j[i] << " " << spms.parity[i] << " " << spms.n[i]
        << " " << spms.tz[i] << ")";
}

// --------------------------------------------------------------------
// PP Modelspace helpers
// --------------------------------------------------------------------

// --------------------------------------------------------------------
// PH Modelspace helpers
// --------------------------------------------------------------------

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

void print_ph_state( std::ostream &o,
                     const ParticleHoleState        &ph,
                     const SingleParticleModelspace &spms ) {
    o << "(";
    print_sp_state( o, ph.ip, spms );
    o << " ";
    print_sp_state( o, ph.ih, spms );
    o << " " << ph.J << ")";
}
