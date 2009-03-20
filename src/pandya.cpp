/* This file contains routines for building G-matrix objects.
 *
 * Mark Burnett November, 2008
 */

#include <cmath>

#include "Modelspace.h"
#include "Interaction.h"
#include "pandya.h"

#include "angular_momentum.h"

// Performs a Pandya transformation (two particle to particle hole).
double pandya( const PPInteraction &G_pp,
               const SingleParticleModelspace &spms,
               const ParticleHoleState &A,
               const ParticleHoleState &B ) {
    int ia = A.ip;
    int ib = A.ih;
    int ic = B.ip;
    int id = B.ih;
    int J = A.J;

    int Jmin = std::max( std::abs( spms.j[ia] - spms.j[id] ),
                                  std::abs( spms.j[ib] - spms.j[ic] ) );
    int Jmax = std::min( spms.j[ia] + spms.j[id],
                                  spms.j[ic] + spms.j[ib] );
    double elem = 0.0;
    for ( int Jp = Jmin; Jp <= Jmax; ++Jp ) {
        ParticleParticleState pp_A( ia, id, -1, -1, Jp );
        ParticleParticleState pp_B( ib, ic, -1, -1, Jp );

        double phase = std::pow(-1.0, spms.j[ib] + spms.j[ic] + Jp);
        double temp = phase * (2*Jp + 1) * G_pp(pp_A, pp_B) *
            wigner6j( spms.j[ia], spms.j[ib], J,
                      spms.j[ic], spms.j[id], Jp );
        elem += temp;
    }

    return elem;
}
