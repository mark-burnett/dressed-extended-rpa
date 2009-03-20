#include <string>
#include <list>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "exceptions.h"
#include "io.h"
#include "angular_momentum.h"

#include "Modelspace.h"

SingleParticleModelspace
read_sp_modelspace_from_file( const std::string &filename ) {
    // Result storage
    std::vector< double >                  j;
    std::vector< int >                     parity;
    std::vector< int >                     n;
    std::vector< double >                  tz;
    std::vector< std::vector< Fragment > > pfrag;
    std::vector< std::vector< Fragment > > hfrag;
    double                                 maxj = 0;
    // Loop through lines and put things where they belong
    BOOST_FOREACH( std::string line, util::read_commented_file( filename ) ) {
        std::vector< std::string > sections = util::split(line, ":");
        { // Quantum Numbers (section 0)
            std::vector< std::string > tokens = util::split( sections[0] );
            j.push_back( boost::lexical_cast<double>(  tokens[2] )/2 );
            n.push_back( boost::lexical_cast<int>(     tokens[1] ) );
            tz.push_back( boost::lexical_cast<double>( tokens[0] )/2 );
            char c_parity = boost::lexical_cast<char>( tokens[3] );
            switch (c_parity) {
                case '-':
                    parity.push_back( -1 );
                    break;
                case '+':
                    parity.push_back(  1 );
                    break;
                default:
                    throw file_error();
            }
            if ( j.back() > maxj )
                maxj = j.back();
        } { // Particle fragments (section 1)
            std::vector< std::string > tokens = util::split( sections[1] );
            std::vector< Fragment > frags;
            for( std::vector< std::string >::const_iterator i = tokens.begin();
                    i != tokens.end(); ++i ) {
                if ( i->empty() )
                    break;
                double E = boost::lexical_cast<double>(*i);
                ++i;
                double S = boost::lexical_cast<double>(*i);
                frags.push_back( Fragment( E, S ) );
            }
            pfrag.push_back( frags );
        } { // Hole fragments (section 2)
            std::vector< std::string > tokens = util::split( sections[2] );
            std::vector< Fragment > frags;
            for( std::vector< std::string >::const_iterator i = tokens.begin();
                    i != tokens.end(); ++i ) {
                if ( i->empty() )
                    break;
                double E = boost::lexical_cast<double>(*i);
                ++i;
                double S = boost::lexical_cast<double>(*i);
                frags.push_back( Fragment( E, S ) );
            }
            hfrag.push_back( frags );
        }
    }

    return SingleParticleModelspace( j, parity, n, tz, pfrag, hfrag, maxj );
}

/*
ParticleParticleModelspace
build_pp_modelspace_from_sp( const SingleParticleModelspace &spms ) {
}
*/

ParticleHoleModelspace
build_ph_modelspace_from_sp( const SingleParticleModelspace &spms ) {
    // Outside index is for isospin
    ParticleHoleModelspace phms(3); phms.resize(3);

    // Resize next vectors for indexing parity
    for ( int i = 0; i < 3; ++i ) {
        phms[i].resize(2);
    }

    // Loop over tz, parity, and J
    for ( int tz = -1; tz <= 1; ++tz ) {
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            // Get maximum J for this tz, parity
            int Jmax = get_max_ph_J( spms, tz, parity );
            phms[ tz + 1 ][ (parity + 1)/2 ].resize( Jmax + 1 );
            for ( int J = 0; J <= Jmax; ++J ) {
                for ( int pshell = 0; pshell < spms.size; ++pshell ) {
                    for ( int hshell = 0; hshell < spms.size; ++hshell ) {
                        // Check for valid isospin and parity
                        if (   tz != spms.tz[pshell]     - spms.tz[hshell] ||
                           parity != spms.parity[pshell] * spms.parity[hshell] )
                            continue;
                        // Check for valid angular momenta
                        if ( !is_triangular(spms.j[pshell], spms.j[hshell], J) )
                            continue;
                        // Loop over fragments
                        for ( int ipf = 0;
                            ipf < static_cast<int>(spms.pfrag[pshell].size());
                                ++ipf ) {
                            for ( int ihf = 0;
                            ihf < static_cast<int>(spms.hfrag[hshell].size());
                                ++ihf ) {
                                phms[ tz + 1 ][ (parity + 1)/2 ][J].push_back(
                                        ParticleHoleState( pshell, hshell,
                                            ipf, ihf, J ) ); } } } } } } }
    return phms;
}

// Ignores fragments
ParticleHoleModelspace
build_ph_shells_from_sp( const SingleParticleModelspace &spms ) {
    // Outside index is for isospin
    ParticleHoleModelspace phms(3); phms.resize(3);

    // Resize next vectors for indexing parity
    for ( int i = 0; i < 3; ++i ) {
        phms[i].resize(2);
    }

    // Loop over tz, parity, and J
    for ( int tz = -1; tz <= 1; ++tz ) {
        for ( int parity = -1; parity <= 1; parity += 2 ) {
            // Get maximum J for this tz, parity
            int Jmax = get_max_ph_J( spms, tz, parity );
            phms[ tz + 1 ][ (parity + 1)/2 ].resize( Jmax + 1 );
            for ( int J = 0; J <= Jmax; ++J ) {
                for ( int pshell = 0; pshell < spms.size; ++pshell ) {
                    for ( int hshell = 0; hshell < spms.size; ++hshell ) {
                        // Check for valid isospin and parity
                        if (   tz != spms.tz[pshell]     - spms.tz[hshell] ||
                           parity != spms.parity[pshell] * spms.parity[hshell] )
                            continue;
                        // Check for valid angular momenta
                        if ( !is_triangular(spms.j[pshell], spms.j[hshell], J) )
                            continue;
                        phms[ tz + 1][ (parity+1)/2 ][J].push_back(
                                ParticleHoleState( pshell, hshell,
                                    -1, -1, J ) ); } } } } }
    return phms;
}
