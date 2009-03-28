#include <gtest/gtest.h>

#include <boost/numeric/conversion/cast.hpp>

#include "modelspace_factories.h"

TEST( ModelspaceFactories, UnfragmentedSPModelspace ) {
    SingleParticleModelspace sp_ms =
        read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    EXPECT_EQ( 20, sp_ms.size );

    EXPECT_EQ( 0.5, sp_ms.j[9] );
    EXPECT_EQ(  -1, sp_ms.parity[9] );
    EXPECT_EQ( 0.5, sp_ms.tz[9] );

    EXPECT_EQ( 0, sp_ms.pfrag[3].size() );
    EXPECT_EQ( 1, sp_ms.hfrag[3].size() );

    EXPECT_FLOAT_EQ( -29.688963, sp_ms.hfrag[3][0].E );
    EXPECT_EQ( 1, sp_ms.hfrag[3][0].S );

    EXPECT_EQ( 3.5, sp_ms.maxj );
}

TEST( ModelspaceFactories, FragmentedSPModelspace ) {
    SingleParticleModelspace sp_ms =
        read_sp_modelspace_from_file( "tests/data/frag_modelspace.dat" );
    EXPECT_EQ( 20, sp_ms.size );
    EXPECT_EQ( 1,  sp_ms.pfrag[0].size() );
    EXPECT_EQ( 2,  sp_ms.hfrag[0].size() );
    EXPECT_EQ( 1,  sp_ms.pfrag[3].size() );
    EXPECT_EQ( 2,  sp_ms.hfrag[3].size() );
    EXPECT_EQ( 2,  sp_ms.pfrag[8].size() );
    EXPECT_EQ( 1,  sp_ms.hfrag[8].size() );
    EXPECT_FLOAT_EQ( -47.981882,       sp_ms.hfrag[2][0].E );
    EXPECT_FLOAT_EQ(   0.806225774830, sp_ms.hfrag[2][0].S );
    EXPECT_FLOAT_EQ(   0.316227766107, sp_ms.hfrag[8][0].S );

    EXPECT_FLOAT_EQ( -8.596114,       sp_ms.pfrag[8][1].E );
    EXPECT_FLOAT_EQ(  0.806225774830, sp_ms.pfrag[8][1].S );

    EXPECT_EQ( 3.5, sp_ms.maxj );
}

TEST( ModelspaceFactories, UnfragmentedPHModelspace ) {
    ParticleHoleModelspace phms = build_ph_modelspace_from_sp(
            read_sp_modelspace_from_file(
                "tests/data/ipm_modelspace.dat" ) );
    // 0+ (tz = 0)
    EXPECT_EQ(  4, phms[ 0 + 1 ][ 1 ][ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ( 19, phms[ 0 + 1 ][ 0 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 13, phms[ 0 + 1 ][ 1 ][ 2 ].size() );
}

TEST( ModelspaceFactories, FragmentedPHModelspace ) {
    ParticleHoleModelspace phms = build_ph_modelspace_from_sp(
            read_sp_modelspace_from_file(
                "tests/data/frag_modelspace.dat" ) );
    // 0+ (tz = 0)
    EXPECT_EQ(  68, phms[ 0 + 1 ][ 1 ][ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ( 155, phms[ 0 + 1 ][ 0 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 173, phms[ 0 + 1 ][ 1 ][ 2 ].size() );
}

TEST( ModelspaceFactories, PPModelspace ) {
    ParticleParticleModelspace ppms = build_pp_modelspace_from_sp(
            read_sp_modelspace_from_file(
                "tests/data/ipm_modelspace.dat" ) );
    // Sanity checks
    EXPECT_EQ( 3, ppms.size() );
    EXPECT_EQ( 2, ppms[1].size() );

    // tz = 0 pp states
    // 0+ (tz = 0)
    EXPECT_EQ(  3, ppms[ 0 + 1 ][ (1 + 1)/2 ][ 0 ].size() );
    // 1+ (tz = 0)
    EXPECT_EQ(  8, ppms[ 0 + 1 ][ (1 + 1)/2 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 10, ppms[ 0 + 1 ][ (1 + 1)/2 ][ 2 ].size() );
    // There are no legal J's for - parity in the test space, so:
    // 0- (tz = 0)
    EXPECT_EQ(  0, ppms[ 0 + 1 ][ (1 - 1)/2 ][ 0 ].size() );
    // 1- (tz = 0)             
    EXPECT_EQ(  0, ppms[ 0 + 1 ][ (1 - 1)/2 ][ 1 ].size() );
    // 2- (tz = 0)             
    EXPECT_EQ(  0, ppms[ 0 + 1 ][ (1 - 1)/2 ][ 2 ].size() );

    // tz = -1 pp states
    // 0+ (tz = 0)
    EXPECT_EQ(  3, ppms[ -1 + 1 ][ (1 + 1)/2 ][ 0 ].size() );
    // 1+ (tz = 0)
    EXPECT_EQ(  5, ppms[ -1 + 1 ][ (1 + 1)/2 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ(  5, ppms[ -1 + 1 ][ (1 + 1)/2 ][ 2 ].size() );
    // There are no legal J's for - parity in the test space, so:
    // 0- (tz = 0)
    EXPECT_EQ(  0, ppms[ -1 + 1 ][ (1 - 1)/2 ][ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ(  0, ppms[ -1 + 1 ][ (1 - 1)/2 ][ 1 ].size() );
    // 2- (tz = 0)
    EXPECT_EQ(  0, ppms[ -1 + 1 ][ (1 - 1)/2 ][ 2 ].size() );

    // tz =  1 pp states
    // 0+ (tz = 0)
    EXPECT_EQ(  4, ppms[  1 + 1 ][ (1 + 1)/2 ][ 0 ].size() );
    // 1+ (tz = 0)
    EXPECT_EQ(  7, ppms[  1 + 1 ][ (1 + 1)/2 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ(  8, ppms[  1 + 1 ][ (1 + 1)/2 ][ 2 ].size() );
    // There are no legal J's for - parity in the test space, so:
    // 0- (tz = 0)
    EXPECT_EQ(  0, ppms[  1 + 1 ][ (1 - 1)/2 ][ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ(  0, ppms[  1 + 1 ][ (1 - 1)/2 ][ 1 ].size() );
    // 2- (tz = 0)
    EXPECT_EQ(  0, ppms[  1 + 1 ][ (1 - 1)/2 ][ 2 ].size() );
}

TEST( ModelspaceFactories, HHModelspace ) {
    ParticleParticleModelspace hhms = build_hh_modelspace_from_sp(
            read_sp_modelspace_from_file(
                "tests/data/ipm_modelspace.dat" ) );
    // Sanity checks
    EXPECT_EQ( 3, hhms.size() );
    EXPECT_EQ( 2, hhms[1].size() );

    // tz = 0 hh states
    // 0+ (tz = 0)
    EXPECT_EQ(  8, hhms[ 0 + 1 ][ (1 + 1)/2 ][ 0 ].size() );
    // 1+ (tz = 0)
    EXPECT_EQ( 16, hhms[ 0 + 1 ][ (1 + 1)/2 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 16, hhms[ 0 + 1 ][ (1 + 1)/2 ][ 2 ].size() );
    // 0- (tz = 0)
    EXPECT_EQ(  6, hhms[ 0 + 1 ][ (1 - 1)/2 ][ 0 ].size() );
    // 1- (tz = 0)             
    EXPECT_EQ( 15, hhms[ 0 + 1 ][ (1 - 1)/2 ][ 1 ].size() );
    // 2- (tz = 0)             
    EXPECT_EQ( 14, hhms[ 0 + 1 ][ (1 - 1)/2 ][ 2 ].size() );

    // tz =  1 hh states
    // 0+ (tz = 0)
    EXPECT_EQ(  7, hhms[  1 + 1 ][ (1 + 1)/2 ][ 0 ].size() );
    // 1+ (tz = 0)
    EXPECT_EQ( 11, hhms[  1 + 1 ][ (1 + 1)/2 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ(  9, hhms[  1 + 1 ][ (1 + 1)/2 ][ 2 ].size() );
    // 0- (tz = 0)
    EXPECT_EQ(  3, hhms[  1 + 1 ][ (1 - 1)/2 ][ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ(  7, hhms[  1 + 1 ][ (1 - 1)/2 ][ 1 ].size() );
    // 2- (tz = 0)
    EXPECT_EQ(  6, hhms[  1 + 1 ][ (1 - 1)/2 ][ 2 ].size() );

    // tz = -1 hh states
    // 0+ (tz = 0)
    EXPECT_EQ(  8, hhms[ -1 + 1 ][ (1 + 1)/2 ][ 0 ].size() );
    // 1+ (tz = 0)
    EXPECT_EQ( 12, hhms[ -1 + 1 ][ (1 + 1)/2 ][ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 11, hhms[ -1 + 1 ][ (1 + 1)/2 ][ 2 ].size() );
    // 0- (tz = 0)
    EXPECT_EQ(  3, hhms[ -1 + 1 ][ (1 - 1)/2 ][ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ(  8, hhms[ -1 + 1 ][ (1 - 1)/2 ][ 1 ].size() );
    // 2- (tz = 0)
    EXPECT_EQ(  8, hhms[ -1 + 1 ][ (1 - 1)/2 ][ 2 ].size() );
}

TEST( ModelspaceFactories, PPSPFromSP ) {
    SingleParticleModelspace spms
        = read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    PPFromSPModelspace ppspms = build_ppsp_modelspace_from_sp( spms );

    // Sanity checks
    EXPECT_EQ( 8, ppspms.size() );
    for ( int i = 0; i < boost::numeric_cast<int>(ppspms.size()); ++i ) {
        EXPECT_EQ( 20, ppspms[i].size() ); }

    // This is the f7/2 proton particle, it only has 1 0+/- state
    EXPECT_EQ( 1, ppspms[0][6].size() );
    // These are the common particles, they have 2 0+/- states
    EXPECT_EQ( 2, ppspms[0][7].size() );
    EXPECT_EQ( 2, ppspms[0][8].size() );
    EXPECT_EQ( 2, ppspms[0][9].size() );
    // f7/2 proton particle
    EXPECT_EQ( 3, ppspms[1][6].size() );
    // These are the common particles
    EXPECT_EQ( 6, ppspms[1][7].size() );
    EXPECT_EQ( 5, ppspms[1][8].size() );
    EXPECT_EQ( 4, ppspms[1][9].size() );

    // f7/2 proton particle
    EXPECT_EQ( 5, ppspms[2][6].size() );
    // These are the common particles
    EXPECT_EQ( 7, ppspms[2][7].size() );
    EXPECT_EQ( 7, ppspms[2][8].size() );
    EXPECT_EQ( 4, ppspms[2][9].size() );
}

TEST( ModelspaceFactories, HHSPFromSP ) {
    SingleParticleModelspace spms
        = read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    PPFromSPModelspace hhspms = build_hhsp_modelspace_from_sp( spms );

    // Sanity checks
    EXPECT_EQ( 8, hhspms.size() );
    for ( int i = 0; i < boost::numeric_cast<int>(hhspms.size()); ++i ) {
        EXPECT_EQ( 20, hhspms[i].size() ); }

    // Some J=0 states
    EXPECT_EQ(  6, hhspms[0][0].size() );
    EXPECT_EQ(  4, hhspms[0][1].size() );
    EXPECT_EQ(  6, hhspms[0][2].size() );
    EXPECT_EQ(  2, hhspms[0][3].size() );
    // Some J=1 states
    EXPECT_EQ( 10, hhspms[1][0].size() );
    EXPECT_EQ( 12, hhspms[1][1].size() );
    EXPECT_EQ( 10, hhspms[1][2].size() );
    EXPECT_EQ(  7, hhspms[1][3].size() );
}
