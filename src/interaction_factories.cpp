#include <cassert>
#include <string>
#include <vector>
#include <fstream>

#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include "io.h"
#include "Interaction.h"
#include "Modelspace.h"
#include "interaction_factories.h"
#include "modelspace_factories.h"
#include "angular_momentum.h"
#include "pandya.h"
#include "exceptions.h"

namespace ublas = boost::numeric::ublas;

// --------------------------------------------------------------------
// Support functions for the PP Interaction factories
// --------------------------------------------------------------------

// This has the usual format:  matricies[ tz + 1 ][ (parity+1)/2 ][ J ]
typedef std::vector< std::vector< std::vector<
    ublas::symmetric_matrix< double >
> > > PPMatricies;

// FIXME these could be simply indices[J]( ip1, ip2 ) (implied parity & tz)
// Format: indices[tz+1][(parity+1)/2][J]( ip1, ip2 )
typedef std::vector< std::vector< std::vector<
    ublas::symmetric_matrix< int >
> > > PPIndices;

// Holds the sizes of the PPMatricies during construction
typedef std::vector< std::vector< std::vector< int > > > Indexsizes;

// This function is the actual PP workhorse.
double pp_interaction_base( const ParticleParticleState    &A,
                            const ParticleParticleState    &B,
                            const PPMatricies              &matricies,
                            const PPIndices               &indices,
                            const SingleParticleModelspace &spms ) {
    int tz     = spms.tz[ A.ip1 ]     + spms.tz[ A.ip2 ];
    int parity = spms.parity[ A.ip1 ] * spms.parity[ A.ip2 ];

    assert( tz     == spms.tz[ B.ip1 ]     + spms.tz[ B.ip2 ] );
    assert( parity == spms.parity[ B.ip1 ] * spms.parity[ B.ip2 ] );
    assert( A.J    == B.J );

    int phase = 1;
    int iA = indices[ tz + 1 ][ (parity+1)/2 ][ A.J ]( A.ip1, A.ip2 );
    int iB = indices[ tz + 1 ][ (parity+1)/2 ][ B.J ]( B.ip1, B.ip2 );

    // These phases take care of antisymmetrizing the interaction.
    if ( A.ip1 > A.ip2 )
        phase *= std::pow( -1, spms.j[ A.ip1 ] - spms.j[ A.ip2 ] + A.J );
    if ( B.ip1 > B.ip2 )
        phase *= std::pow( -1, spms.j[ B.ip1 ] - spms.j[ B.ip2 ] + B.J );

    return phase * matricies[ tz + 1 ][ (parity+1)/2 ][ A.J ] ( iA, iB );
}

// Returns index that corresponds to the spms index for a state
int get_ms_index_from_line( const std::string &line,
                            const SingleParticleModelspace &spms ) {
    std::vector< std::string > tokens = util::split( line );

    double j   = boost::lexical_cast<double>(   tokens[4] ) / 2;
    int L      = boost::lexical_cast<int>(      tokens[3] );
    int parity = std::pow(-1.0, L);
    int n      = boost::lexical_cast<int>(      tokens[2] );
    double tz  = - boost::lexical_cast<double>( tokens[5] ) / 2;

    for ( int i = 0; i < spms.size; ++i ) {
        if ( spms.j[i]      == j      &&
             spms.n[i]      == n      &&
             spms.parity[i] == parity &&
             spms.tz[i]     == tz        )
            return i; }
    return -1; // state not found -> return illegal index
}

// Builds a vector that maps from the modelspace numbering in the file
// to the indices used in spms
std::vector< int >
build_modelspace_map_from_stream( std::istream &file,
                            const SingleParticleModelspace &spms ) {
    std::string line; // temporary used repeatedly for getline

    // Get the number of shells in the modelspace
    std::getline( file, line );
    int ms_size = boost::lexical_cast< int >(
            util::split( line ).back() );

    // Skip a line
    std::getline( file, line );

    // Read in the shells & construct mapping from one space to other
    std::vector< int > modelspace_map( ms_size );
    for ( int i = 0; i < ms_size; ++i ) {
        std::getline( file, line );
        modelspace_map[i] = get_ms_index_from_line( line, spms ); }

    return modelspace_map;
}

// Makes the PP index lookup for a given tz, parity, J
boost::tuple< ublas::symmetric_matrix< int >, int >
make_pp_indices ( int tz, int parity, int J,
                   const SingleParticleModelspace &spms ) {
    // Initialize the result to -1 (illegal index)
    ublas::symmetric_matrix< int > result( spms.size, spms.size );

    int size = 0; // the number of legal elements we find
    for ( int i = 0; i < spms.size; ++i ) {
        for ( int j = i; j < spms.size; ++j ) {
            if ( is_triangular( spms.j[i], spms.j[j], J )  &&
                 spms.parity[i] * spms.parity[j] == parity &&
                 spms.tz[i]     + spms.tz[j]     == tz        ) {
                result( i, j ) = size;
                ++size; }
            else
                result( i, j ) = -1; } }
    return boost::make_tuple( result, size );
}

// Builds a lookup for finding PP matrix element indices given
// single particle configurations
boost::tuple< PPIndices, Indexsizes >
build_ppindices( const SingleParticleModelspace &spms ) {
    PPIndices ppi(3); ppi.resize(3);
    Indexsizes sizes(3); sizes.resize(3);

    for ( int tz = -1; tz <= 1; ++tz ) {
        ppi  [ tz + 1 ].resize(2);
        sizes[ tz + 1 ].resize(2);
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            int maxJ = get_max_pp_J( spms, tz, parity );
            ppi  [ tz + 1 ][ (parity + 1)/2 ].resize( maxJ + 1 );
            sizes[ tz + 1 ][ (parity + 1)/2 ].resize( maxJ + 1 );
            for ( int J = 0;
                    J < static_cast<int>(ppi[tz + 1][(parity + 1)/2].size());
                    ++J ) {
                boost::tie( ppi  [ tz + 1 ][ (parity + 1)/2 ][ J ],
                            sizes[ tz + 1 ][ (parity + 1)/2 ][ J ] )
                        = make_pp_indices( tz, parity, J, spms ); } } }
    return boost::make_tuple( ppi, sizes );
}

boost::tuple< int, int, int, int, int, int, int, double >
decode_mhj_line( const std::string &line,
                 const std::vector< int > &modelspace_map,
                 const SingleParticleModelspace &spms ) {
    std::vector< std::string > tokens = util::split( line );

    int tz     = - boost::lexical_cast<int>( tokens[0] );
    int parity = std::pow( -1.0, boost::lexical_cast<int>( tokens[1] ) );
    int J      = boost::lexical_cast<int>( tokens[2] ) / 2;

    int ua   = boost::lexical_cast<int>( tokens[3] ) - 1;
    int ub   = boost::lexical_cast<int>( tokens[4] ) - 1;
    int uc   = boost::lexical_cast<int>( tokens[5] ) - 1;
    int ud   = boost::lexical_cast<int>( tokens[6] ) - 1;

    double V = boost::lexical_cast<double>( tokens[7] );

    int ia = modelspace_map[ ua ];
    int ib = modelspace_map[ ub ];
    int ic = modelspace_map[ uc ];
    int id = modelspace_map[ ud ];

    // correct for phase
    if ( ua > ub && ia > ib )
            V *= std::pow( -1, J + spms.j[ia] - spms.j[ib] );
    if ( uc > ud && ic > ib )
            V *= std::pow( -1, J + spms.j[ic] - spms.j[id] );

    // renormalize
    if ( ia == ib )
        V *= std::sqrt(2.0);
    if ( ic == id )
        V *= std::sqrt(2.0);

    return boost::make_tuple( tz, parity, J, ia, ib, ic, id, V );
}

// Reads in the actual matrix elements
PPMatricies
build_matricies_from_stream( std::ifstream &file,
                             const PPIndices &indices,
                             const Indexsizes &sizes,
                             const SingleParticleModelspace &spms,
                             const std::vector< int > &modelspace_map ) {
    // Initialize matricies
    PPMatricies matricies(3); matricies.resize(3);
    for ( int tz = -1; tz <= 1; ++tz ) {
        matricies[ tz + 1 ].resize(2); // resize parity
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            matricies[ tz + 1 ][ (parity + 1)/2 ].resize( // resize J
                    indices[ tz + 1 ][ (parity + 1)/2 ].size() );
            for ( int J = 0;
                    J < static_cast<int>(matricies[tz+1][(parity+1)/2].size());
                    ++J ) {
                // add empty matricies
                int size = sizes[ tz + 1 ][ (parity + 1)/2 ][ J ];
                ublas::symmetric_matrix< double > m( size, size );
                m.clear();
                matricies[ tz + 1 ][ (parity + 1)/2 ][ J ] = m; } } }

    // loop over stream
    std::string line; // temporary used repeatedly for getline
    while ( std::getline( file, line ) ) {
        // decode line
        int tz, parity, J, ia, ib, ic, id;
        double V;
        boost::tie( tz, parity, J, ia, ib, ic, id, V )
            = decode_mhj_line( line, modelspace_map, spms );

        // lookup indices
        int iA = indices[ tz + 1 ][ (parity + 1)/2 ][ J ]( ia, ib );
        int iB = indices[ tz + 1 ][ (parity + 1)/2 ][ J ]( ic, id );


        // add element
        matricies[ tz + 1 ][ (parity + 1)/2 ][ J ]( iA, iB ) = V;
    }

    return matricies;
}

// --------------------------------------------------------------------
// The actual PPInteraction factories:
// --------------------------------------------------------------------

PPInteraction
build_gmatrix_from_mhj_file( const std::string &filename,
                             const SingleParticleModelspace &spms ) {
    std::ifstream file( filename.c_str() );

    std::string line; // temporary used repeatedly for getline

    // Skip the first 8 lines
    for ( int i = 0; i < 8; ++i )
        if ( !std::getline( file, line ) )
            throw file_error();

    // Construct modelspace map and PPIndices
    std::vector< int > modelspace_map
        = build_modelspace_map_from_stream( file, spms );
    PPIndices indices;
    Indexsizes index_sizes;
    boost::tie( indices, index_sizes ) = build_ppindices( spms );

    // Skip 8 more lines
    for ( int i = 0; i < 8; ++i )
        if ( !std::getline( file, line ) )
            throw file_error();

    // Construct PPMatricies
    PPMatricies matricies
        = build_matricies_from_stream( file, indices, index_sizes,
                                       spms, modelspace_map );

    // Bind it all together
    return boost::bind( pp_interaction_base, _1, _2,
            matricies, indices, spms );
}

// --------------------------------------------------------------------
// Support functions for the PH Interaction factories
// --------------------------------------------------------------------

// PH interaction elements
typedef std::vector< std::vector< std::vector<
    ublas::symmetric_matrix< double >
> > > PHMatricies;

// Format: indices[tz+1][(parity+1)/2][J]( ip, ih )
typedef std::vector< std::vector< std::vector<
    ublas::symmetric_matrix< int >
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
    int iB = indices[ tz + 1 ][ (parity+1)/2 ][ B.J ]( B.ip, B.ih );

    return matricies[ tz + 1 ][ (parity+1)/2 ][ A.J ]( iA, iB );
}

ublas::symmetric_matrix< int >
make_ph_indices( int tz, int parity, int J,
                 const SingleParticleModelspace &spms ) {
    // Initialize the result to -1 (illegal index)
    ublas::symmetric_matrix< int > result( spms.size, spms.size );

    int size = 0; // the number of legal elements we find
    for ( int i = 0; i < spms.size; ++i ) {
        for ( int j = i; j < spms.size; ++j ) {
            if ( is_triangular( spms.j[i], spms.j[j], J )  &&
                 spms.parity[i] * spms.parity[j] == parity &&
                 spms.tz[i]     - spms.tz[j]     == tz        ) {
                result( i, j ) = size;
                ++size; }
            else
                result( i, j ) = -1; } }
    return result;
}

PHIndices
build_ph_indices( const SingleParticleModelspace &spms ) {
    PHIndices indices(3); indices.resize(3);
    for ( int tz = -1; tz <= 1; ++tz ) {
        indices[ tz + 1 ].resize(2); // resize parity
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            int Jmax = get_max_ph_J( spms, tz, parity );
            indices[ tz + 1 ][ (parity + 1)/2 ].resize( Jmax + 1 );
            for ( int J = 0; J <= Jmax; ++J ) {
                indices[ tz + 1 ][ (parity + 1)/2 ][ J ]
                    = make_ph_indices( tz, parity, J, spms ); } } }

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
                    for ( int j = i; j < size; ++j ) {
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
    // build indicies
    PHIndices indices;
    indices = build_ph_indices( spms );

    // build shells
    ParticleHoleModelspace shells = build_ph_shells_from_sp( spms );

    // build matricies
    PHMatricies matricies = build_ph_matricies_from_pp( Gpp, indices,
            spms, shells );

    return boost::bind( ph_interaction_base, _1, _2,
            matricies, indices, spms );
}
