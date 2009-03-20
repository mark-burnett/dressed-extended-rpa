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

// --------------------------------------------------------------------
// PP Modelspace helpers
// --------------------------------------------------------------------

// --------------------------------------------------------------------
// PH Modelspace helpers
// --------------------------------------------------------------------

