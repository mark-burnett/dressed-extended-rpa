#include <cassert>
#include <string>
#include <vector>
#include <fstream>

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include "io.h"
#include "Interaction.h"
#include "Modelspace.h"
#include "ph_interaction_factories.h"
#include "modelspace_factories.h"
#include "angular_momentum.h"
#include "pandya.h"
#include "exceptions.h"

namespace ublas = boost::numeric::ublas;

// --------------------------------------------------------------------
// Support functions for the PH Interaction factories
// --------------------------------------------------------------------

// PH interaction elements
typedef std::vector< std::vector< std::vector<
    ublas::symmetric_matrix< double >
> > > PHMatricies;

// Format: indices[tz+1][(parity+1)/2][J]( ip, ih )
typedef std::vector< std::vector< std::vector<
    ublas::matrix< int >
> > > PHIndices;

double ph_interaction_base( const ParticleHoleState    &A,
                            const ParticleHoleState    &B,
                            const PHMatricies             &matricies,
                            const PHIndices               &indices,
                            const SingleParticleModelspace &spms ) {
    int tz     = spms.tz[ A.ip ]     - spms.tz[ A.ih ];
    int parity = spms.parity[ A.ip ] * spms.parity[ A.ih ];

    assert( tz     == spms.tz[ B.ip ]     - spms.tz[ B.ih ] );
    assert( parity == spms.parity[ B.ip ] * spms.parity[ B.ih ] );
    assert( A.J    == B.J );

    int iA = indices[ tz + 1 ][ (parity+1)/2 ][ A.J ]( A.ip, A.ih );
    int iB = indices[ tz + 1 ][ (parity+1)/2 ][ A.J ]( B.ip, B.ih );

    return matricies[ tz + 1 ][ (parity+1)/2 ][ A.J ]( iA, iB );
}

ublas::matrix< int >
make_ph_indices( const std::vector< ParticleHoleState > &shell,
                 const SingleParticleModelspace &spms ) {
    ublas::matrix< int > result( spms.size, spms.size );
    // Initialize the result to -1 (illegal index)
    for ( int i = 0; i < spms.size; ++i ) {
        for ( int j = 0; j < spms.size; ++j ) {
            result( i, j ) = -1; } }
    
    // Loop over shell
    for ( int i = 0; i < boost::numeric_cast<int>(shell.size()); ++i ) {
        int ip = shell[i].ip;
        int ih = shell[i].ih;
        result( ip, ih ) = i; }

    return result;
}

PHIndices
build_ph_indices( const ParticleHoleModelspace   &shells,
                  const SingleParticleModelspace &spms ) {
    PHIndices indices(3); indices.resize(3);
    for ( int tz = -1; tz <= 1; ++tz ) {
        indices[ tz + 1 ].resize(2); // resize parity
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            int Jmax = get_max_ph_J( spms, tz, parity );
            indices[ tz + 1 ][ (parity + 1)/2 ].resize( Jmax + 1 );
            for ( int J = 0; J <= Jmax; ++J ) {
                indices[ tz + 1 ][ (parity + 1)/2 ][ J ]
                    = make_ph_indices( shells[tz+1][(parity+1)/2][J],
                                       spms ); } } }
    return indices;
}

PHMatricies
build_ph_matricies_from_pp( const PPInteraction &Gpp,
                            const PHIndices &indices,
                            const SingleParticleModelspace &spms,
                            const ParticleHoleModelspace &shells ) {
    PHMatricies matricies(3); matricies.resize(3);
    for ( int tz = -1; tz <= 1; ++tz ) {
        matricies[ tz + 1 ].resize(2); // resize parity
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            matricies[ tz + 1 ][ (parity + 1)/2 ].resize( // resize J
                    indices[ tz + 1 ][ (parity + 1)/2 ].size() );
            for ( int J = 0;
                    J < static_cast<int>(matricies[tz+1][(parity+1)/2].size());
                    ++J ) {
                std::vector< ParticleHoleState > shell
                    = shells[ tz + 1 ][ (parity + 1)/2 ][ J ];
                int size = shell.size();
                ublas::symmetric_matrix< double > m( size, size );
                m.clear();
                for ( int i = 0; i < size; ++i ) {
                    for ( int j = 0; j < size; ++j ) {
                        m( i, j ) = pandya( Gpp, spms, shell[i], shell[j] );
                        matricies[ tz + 1 ][ (parity + 1)/2 ][ J ] = m; }} } } }
    return matricies;
}

// --------------------------------------------------------------------
// Actual PH Interaction
// --------------------------------------------------------------------
PHInteraction
build_ph_interaction_from_pp( const PPInteraction &Gpp,
                              const SingleParticleModelspace &spms ) {
    // build shells
    ParticleHoleModelspace shells = build_ph_shells_from_sp( spms );

    // build indicies
    PHIndices indices;
    indices = build_ph_indices( shells, spms );

    // build matricies
    PHMatricies matricies = build_ph_matricies_from_pp( Gpp, indices,
            spms, shells );

    return boost::bind( ph_interaction_base, _1, _2,
            matricies, indices, spms );
}
