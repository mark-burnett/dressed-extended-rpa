#include <gtest/gtest.h>

#include "Modelspace.h"
#include "Interaction.h"
#include "modelspace_factories.h"
#include "pp_interaction_factories.h"
#include "ph_interaction_factories.h"

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

TEST( InteractionFactories, PHFromPP ) {
    SingleParticleModelspace spms =
        read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    PPInteraction Gpp = build_gmatrix_from_mhj_file(
            "tests/data/test_interaction.mhj", spms );
    PHInteraction Gph = build_ph_interaction_from_pp( Gpp, spms );

    typedef ParticleHoleState ph_t;

    EXPECT_FLOAT_EQ( -1.7639774,
            Gph( ph_t( 9, 2, -1, -1, 0 ), ph_t( 18, 12, -1, -1, 0 ) ) );

    EXPECT_FLOAT_EQ(  4.0043716357,
            Gph( ph_t( 0, 10, -1, -1, 1 ), ph_t( 1, 11, -1, -1, 1 ) ) );

    EXPECT_FLOAT_EQ(  1.4165830803,
            Gph( ph_t( 9, 1, -1, -1, 2 ), ph_t( 17, 12, -1, -1, 2 ) ) );

    EXPECT_FLOAT_EQ( -1.4824042775,
            Gph( ph_t( 7, 1, -1, -1, 2 ), ph_t( 17, 12, -1, -1, 2 ) ) );

    EXPECT_FLOAT_EQ(  0.3647625558,
            Gph( ph_t( 9, 1, -1, -1, 2 ), ph_t( 7, 1, -1, -1, 2 ) ) );

    EXPECT_FLOAT_EQ( -0.260148,
            Gph( ph_t( 18, 13, -1, -1, 3 ), ph_t( 8, 3, -1, -1, 3 ) ) );
    
    // These point out a bug in Gpp
    EXPECT_FLOAT_EQ(  0.22191499,
            Gph( ph_t( 9, 1, -1, -1, 2 ), ph_t( 1, 7, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.3647625558,
            Gph( ph_t( 1, 9, -1, -1, 2 ), ph_t( 1, 7, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.22191499,
            Gph( ph_t( 1, 9, -1, -1, 2 ), ph_t( 7, 1, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.44136867,
            Gph( ph_t( 7, 9, -1, -1, 2 ), ph_t( 1, 1, -1, -1, 2 ) ) );

    // Trying to reproduce differences between code versions.
    EXPECT_FLOAT_EQ( -0.70460743,
            Gph( ph_t( 6, 1, -1, -1, 2 ), ph_t( 19, 16, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ(  1.14281230,
            Gph( ph_t( 6, 1, -1, -1, 2 ), ph_t( 16, 19, -1, -1, 2 ) ) );

    // These tests replicate a bug!
    EXPECT_FLOAT_EQ( -0.903618,
            Gph( ph_t( 7, 4, -1, -1, 0 ), ph_t( 3, 8, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( -0.43559521,
            Gph( ph_t( 7, 4, -1, -1, 0 ), ph_t( 8, 3, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( -0.903618,
            Gph( ph_t( 4, 7, -1, -1, 0 ), ph_t( 8, 3, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ( -0.43559521,
            Gph( ph_t( 4, 7, -1, -1, 0 ), ph_t( 3, 8, -1, -1, 0 ) ) );
}
