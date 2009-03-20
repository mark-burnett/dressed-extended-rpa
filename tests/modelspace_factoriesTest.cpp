#include <gtest/gtest.h>

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
    EXPECT_EQ(  4, phms[ 0 + 1 ][ 1 ] [ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ( 19, phms[ 0 + 1 ][ 0 ] [ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 13, phms[ 0 + 1 ][ 1 ] [ 2 ].size() );
}

TEST( ModelspaceFactories, FragmentedPHModelspace ) {
    ParticleHoleModelspace phms = build_ph_modelspace_from_sp(
            read_sp_modelspace_from_file(
                "tests/data/frag_modelspace.dat" ) );
    // 0+ (tz = 0)
    EXPECT_EQ(  68, phms[ 0 + 1 ][ 1 ] [ 0 ].size() );
    // 1- (tz = 0)
    EXPECT_EQ( 155, phms[ 0 + 1 ][ 0 ] [ 1 ].size() );
    // 2+ (tz = 0)
    EXPECT_EQ( 173, phms[ 0 + 1 ][ 1 ] [ 2 ].size() );
}
