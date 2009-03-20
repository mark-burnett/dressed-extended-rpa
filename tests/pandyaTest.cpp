#include <gtest/gtest.h>

#include "pandya.h"
#include "Modelspace.h"
#include "Interaction.h"
#include "modelspace_factories.h"
#include "interaction_factories.h"

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

TEST( Pandya, Pandya ) {
    SingleParticleModelspace spms =
        read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    PPInteraction Gpp = build_gmatrix_from_mhj_file(
            "tests/data/test_interaction.mhj", spms );

    typedef ParticleHoleState ph_t;

    EXPECT_FLOAT_EQ(  6.4905633211, pandya( Gpp, spms,
            ph_t( 0, 10, -1, -1, 0 ), ph_t( 1, 11, -1, -1, 0 ) ) );
    EXPECT_FLOAT_EQ(  4.0043716357, pandya( Gpp, spms,
            ph_t( 0, 10, -1, -1, 1 ), ph_t( 1, 11, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ( -5.4604944101, pandya( Gpp, spms,
            ph_t( 0, 1, -1, -1, 1 ), ph_t( 10, 11, -1, -1, 1 ) ) );
    EXPECT_FLOAT_EQ( -0.7428214168, pandya( Gpp, spms,
            ph_t( 6, 1, -1, -1, 3 ), ph_t( 16, 11, -1, -1, 3 ) ) );

    // Transformations that have failed in the past
    EXPECT_FLOAT_EQ(  0.0846002920, pandya( Gpp, spms,
            ph_t( 6, 1, -1, -1, 2 ), ph_t( 6, 1, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ( -0.8093589423, pandya( Gpp, spms,
            ph_t( 6, 1, -1, -1, 2 ), ph_t( 1, 6, -1, -1, 2 ) ) );
    EXPECT_FLOAT_EQ(  1.4165830802, pandya( Gpp, spms,
            ph_t( 9, 1, -1, -1, 2 ), ph_t( 17, 12, -1, -1, 2 ) ) );
}
