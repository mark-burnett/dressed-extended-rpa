#ifndef _MODELSPACE_H_
#define _MODELSPACE_H_
/* This modelspace design is about as simple as it gets.
 * Everything is referenced using indicies, and all quantities
 * associated with a state are accessed through the modelspace using the
 * state's index.
 *
 * Examples:
 *  SingleParticleModelspace sp_ms( ... );
 *  int isp    = 3;
 *  int ipfrag = 0;
 *  int ihfrag = 1;
 *
 *  // Quantum numbers:
 *      angular_momentum = sp_ms.j[isp];
 *
 *  // Fragment strength/energy:
 *      particle_spectral_strength = sp_ms.pfrag[isp][iprag].S;
 *      particle_fragment_location = sp_ms.pfrag[isp][iprag].E;
 *      hole_spectral_strength     = sp_ms.hfrag[isp][ihfrag].S;
 *      hole_fragment_location     = sp_ms.hfrag[isp][ihfrag].E;
 * 
 *  ParticleHoleModelspace ph_ms( ... );
 *
 *  // PH & PP ms access format:
 *  int iph = 2;
 *  ParticleHoleState ph = ph_ms[ tz + 1 ][ (parity + 1)/2 ][ J ][ iph ];
 *  
 * Mark Burnett, March 2009
 */

#include <iostream>
#include <vector>

// Constant public data members allow for a simple interface.
struct Fragment {
    Fragment( double nE, double nS )
        : E(nE), S(nS) { }
    double E;
    double S;
};

struct SingleParticleModelspace {
    SingleParticleModelspace(
            const std::vector< double >                  &nj,
            const std::vector< int >                     &nparity,
            const std::vector< int >                     &nn,
            const std::vector< double >                  &ntz,
            const std::vector< std::vector< Fragment > > &npfrag,
            const std::vector< std::vector< Fragment > > &nhfrag,
            double                                        nmaxj )
        : j( nj ), parity( nparity ), n( nn ), tz( ntz ), pfrag( npfrag ),
          hfrag( nhfrag ), maxj( nmaxj ), size( nj.size() ) { }

    const std::vector< double >                  j;
    const std::vector< int >                     parity;
    const std::vector< int >                     n;
    const std::vector< double >                  tz;
    const std::vector< std::vector< Fragment > > pfrag;
    const std::vector< std::vector< Fragment > > hfrag;
    const double                                 maxj;
    const int                                    size;
};

// ip and ih are indecies for sp states in the sp modelspace
// ipf and ihf are fragment indicies for those sp states
struct ParticleHoleState {
    ParticleHoleState( int nip,  int nih, int nipf, int nihf, int nJ )
        : ip( nip ), ih( nih ), ipf( nipf ), ihf( nihf ), J( nJ ) { }
    int ip;
    int ih;
    int ipf;
    int ihf;
    int J;
    bool operator<( const ParticleHoleState &sister ) const {
        if ( ip != sister.ip )
            return ip < sister.ip;
        if ( ih != sister.ih )
            return ih < sister.ih;
        if ( ipf != sister.ipf )
            return ipf < sister.ipf;
        if ( ihf != sister.ihf )
            return ihf < sister.ihf;
        return J < sister.J;
    }
};

// ip1 and ip2 are indecies for sp states in the sp modelspace
// ip1f and ip2f are fragment indicies for those sp states
struct ParticleParticleState {
    ParticleParticleState( int nip1,  int nip2, int nip1f, int nip2f, int nJ )
        : ip1( nip1 ), ip2( nip2 ), ip1f( nip1f ), ip2f( nip2f ), J( nJ ) { }
    int ip1;
    int ip2;
    int ip1f;
    int ip2f;
    int J;
};

// See comment at top of file for usage.
typedef std::vector< std::vector< std::vector<
    std::vector< ParticleHoleState >
> > > ParticleHoleModelspace;

typedef std::vector< std::vector< std::vector<
    std::vector< ParticleParticleState >
> > > ParticleParticleModelspace;

// This type of space is used for the self energy terms in the ERPA
// space[J][ isp ] -> vector< PPState >
typedef std::vector< std::vector<
    std::vector< ParticleParticleState >
> > PPFromSPModelspace;

// Some SP Modelspace functions
int get_max_pp_J( const SingleParticleModelspace &spms, int tz, int parity );
int get_max_ph_J( const SingleParticleModelspace &spms, int tz, int parity );

// Some PH Modelspace functions
double ph_energy( const ParticleHoleState &ph,
                  const SingleParticleModelspace &spms );
std::vector< double > ph_poles( int tz, int parity, int J,
                                const ParticleHoleModelspace &phms,
                                const SingleParticleModelspace &spms );

// IO Functions
void print_ph_modelspace_sizes( std::ostream &o,
                                const ParticleHoleModelspace &phms );

void print_ph_modelspace_sizes( std::ostream &o, int tz,
                                const ParticleHoleModelspace &phms );

void print_sp_state( std::ostream &o, int i,
                     const SingleParticleModelspace &spms );
void print_ph_state( std::ostream &o,
                     const ParticleHoleState        &ph,
                     const SingleParticleModelspace &spms );

#endif // _MODELSPACE_H_
