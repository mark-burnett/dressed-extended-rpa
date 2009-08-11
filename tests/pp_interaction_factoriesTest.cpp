#include <gtest/gtest.h>

#include "pandya.h"
#include "Modelspace.h"
#include "Interaction.h"
#include "modelspace_factories.h"
#include "pp_interaction_factories.h"

#include <cmath>

/* Single particle modelspace map for test interaction file
 * mhj file index - 1 -> my modelspace index
 * 0 -> 0
 * 1 -> 10
 * 2 -> 1
 * 3 -> 11
 * 4 -> 2
 * 5 -> 12
 * 6 -> 3
 * 7 -> 13
 * 8 -> 4
 * 9 -> 14
 * 10 -> 5
 * 11 -> 15
 * 12 -> 6
 * 13 -> 16
 * 14 -> 8
 * 15 -> 19
 * 16 -> 7
 * 17 -> 17
 * 18 -> 9
 * 19 -> 18
 */

TEST( InteractionFactories, PPInteractionFromMHJ ) {
    SingleParticleModelspace spms =
        read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    PPInteraction Gpp = build_gmatrix_from_mhj_file(
            "tests/data/test_interaction.mhj", spms );

    typedef ParticleParticleState pp_t;

    // A few tests that were wrong a long time ago
    EXPECT_FLOAT_EQ( -0.2344912000,
            Gpp( pp_t( 1, 16, -1, -1, 2 ), pp_t( 6, 11, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.3222744927,
            Gpp( pp_t( 1, 16, -1, -1, 3 ), pp_t( 6, 11, -1, -1, 3 ) ) );
    EXPECT_FLOAT_EQ( -0.810887,
            Gpp( pp_t( 10, 11, -1, -1, 2 ), pp_t( 10, 11, -1, -1, 2 ) ) );

    // Check a couple of things that may be failing
    EXPECT_FLOAT_EQ( -0.090404481,
            Gpp( pp_t( 17, 3, -1, -1, 3 ), pp_t( 13, 8, -1, -1, 3 ) ) );
    EXPECT_FLOAT_EQ( -0.022610215,
            Gpp( pp_t( 17, 3, -1, -1, 2 ), pp_t( 13, 8, -1, -1, 2 ) ) );

    // This set was failing.
    EXPECT_FLOAT_EQ( -0.02807734,
            Gpp( pp_t( 7, 8, -1, -1, 1 ), pp_t( 4, 3, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ( -0.25176015,
            Gpp( pp_t( 7, 8, -1, -1, 2 ), pp_t( 4, 3, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.033133,
            Gpp( pp_t( 7, 8, -1, -1, 3 ), pp_t( 4, 3, -1, -1, 3 ) ) );
    EXPECT_FLOAT_EQ( -0.31687135,
            Gpp( pp_t( 7, 8, -1, -1, 4 ), pp_t( 4, 3, -1, -1, 4 ) ) );

    // A few 2+ tests
    EXPECT_FLOAT_EQ( -0.7497247700,
            Gpp( pp_t( 6, 1, -1, -1, 2 ), pp_t( 1, 6, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ(  0.8575491513,
            Gpp( pp_t( 0, 3, -1, -1, 2 ), pp_t( 2, 7, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.3228936712,
            Gpp( pp_t( 2, 7, -1, -1, 2 ), pp_t( 4, 5, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.0447950562,
            Gpp( pp_t( 3, 4, -1, -1, 2 ), pp_t( 3, 4, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ(  0.1354613,
            Gpp( pp_t( 3, 4, -1, -1, 2 ), pp_t( 6, 8, -1, -1, 2 ) ) );

    // Some 2- tests
    EXPECT_FLOAT_EQ(  0.0111985897,
            Gpp( pp_t( 0, 8, -1, -1, 2 ), pp_t( 1, 3, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.0146918729,
            Gpp( pp_t( 2, 3, -1, -1, 2 ), pp_t( 4, 7, -1, -1, 2 ) ) );

    // 0+
    EXPECT_FLOAT_EQ(  0.9258414338,
            Gpp( pp_t( 9, 2, -1, -1, 0 ), pp_t( 5, 0, -1, -1, 0 ) ) );


// This test includes normal hits plus phase changes.
    // This is a test case where J = 0, so there can be no phase difference.
    // 1 12, 3 4, J = 0
    EXPECT_FLOAT_EQ( -0.3966205694,
            Gpp( pp_t( 0, 15, -1, -1, 0 ), pp_t( 1, 11, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( -0.3966205694,
            Gpp( pp_t( 15, 0, -1, -1, 0 ), pp_t( 1, 11, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( -0.3966205694,
            Gpp( pp_t( 15, 0, -1, -1, 0 ), pp_t( 11, 1, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( -0.3966205694,
            Gpp( pp_t( 0, 15, -1, -1, 0 ), pp_t( 11, 1, -1, -1, 0 ) ) );
    // Now J = 1, with just one phase change
    // 3 4, 15 4, J = 1
    EXPECT_FLOAT_EQ( -2.074107922,
            Gpp( pp_t( 1, 11, -1, -1, 1 ), pp_t( 8, 11, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ(  2.074107922,
            Gpp( pp_t( 11, 1, -1, -1, 1 ), pp_t( 8, 11, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ(  2.074107922,
            Gpp( pp_t( 11, 1, -1, -1, 1 ), pp_t( 11, 8, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ( -2.074107922,
            Gpp( pp_t( 1, 11, -1, -1, 1 ), pp_t( 11, 8, -1, -1, 1 ) ) );

    // Now J = 1, Without any phase changes
    // 2 4, 6 10, J = 1
    EXPECT_FLOAT_EQ(  1.979501969,
            Gpp( pp_t( 10, 11, -1, -1, 1 ), pp_t( 12, 14, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ(  1.979501969,
            Gpp( pp_t( 11, 10, -1, -1, 1 ), pp_t( 12, 14, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ(  1.979501969,
            Gpp( pp_t( 11, 10, -1, -1, 1 ), pp_t( 14, 12, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ(  1.979501969,
            Gpp( pp_t( 10, 11, -1, -1, 1 ), pp_t( 14, 12, -1, -1, 1 ) ) );

    // Now J = 1, with two phase changes
    // 1 11, 3 17, J = 1
    EXPECT_FLOAT_EQ(  0.008798113, 
            Gpp( pp_t( 0, 5, -1, -1, 1 ), pp_t( 1, 7, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ( -0.008798113, 
            Gpp( pp_t( 5, 0, -1, -1, 1 ), pp_t( 1, 7, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ(  0.008798113, 
            Gpp( pp_t( 5, 0, -1, -1, 1 ), pp_t( 7, 1, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ( -0.008798113, 
            Gpp( pp_t( 0, 5, -1, -1, 1 ), pp_t( 7, 1, -1, -1, 1 ) ) );

    // Check for sqrt(2) normalization factor
    // Two factors
    EXPECT_FLOAT_EQ( 2 * 0.7004450235,
            Gpp( pp_t( 5, 5, -1, -1, 0 ), pp_t( 8, 8, -1, -1, 0 ) ) );
    // one factor
    EXPECT_FLOAT_EQ( std::sqrt(2.0) * -0.733220442,
            Gpp( pp_t( 6, 6, -1, -1, 6 ), pp_t( 6, 8, -1, -1, 6 ) ) );
}
