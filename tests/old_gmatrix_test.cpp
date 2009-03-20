#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>
#include <algorithm>

#include "nuclear/nucleus.h"
#include "nuclear/NuclearState.h"
#include "nuclear/pandya.h"

#include "rpa/Modelspace.h"
#include "rpa/modelspace_factories.h"
//#include "rpa/gmatrix.h"
#include "rpa/gmatrix_factories.h"

// Declare some global data for these tests
class GMatrix : public testing::Test {
    protected:
        typedef nuclear::sp_state sp;
        typedef nuclear::pp_state pp;
        typedef nuclear::ph_state ph;
        static std::vector< sp > s;
        static nuclear::pp_interaction_t G_pp;
};
nuclear::pp_interaction_t GMatrix::G_pp =
    rpa::factories::build_gmatrix_from_mhj_file(
                "data/Gpp_n3loC_hw12_e3fl.mhj" );

// sp states as labeled in the test interaction file
std::vector< nuclear::sp_state > GMatrix::s = boost::assign::list_of
    ( sp()                 )  // NULL entry (n = -1 will not be represented)
    ( sp(0.5,  1, 0,  0.5) )  // Number:  1
    ( sp(0.5,  1, 0, -0.5) )  // Number:  2
    ( sp(1.5, -1, 0,  0.5) )  // Number:  3
    ( sp(1.5, -1, 0, -0.5) )  // Number:  4
    ( sp(0.5, -1, 0,  0.5) )  // Number:  5
    ( sp(0.5, -1, 0, -0.5) )  // Number:  6
    ( sp(2.5,  1, 0,  0.5) )  // Number:  7
    ( sp(2.5,  1, 0, -0.5) )  // Number:  8
    ( sp(1.5,  1, 0,  0.5) )  // Number:  9
    ( sp(1.5,  1, 0, -0.5) )  // Number: 10
    ( sp(0.5,  1, 1,  0.5) )  // Number: 11
    ( sp(0.5,  1, 1, -0.5) )  // Number: 12
    ( sp(3.5, -1, 0,  0.5) )  // Number: 13
    ( sp(3.5, -1, 0, -0.5) )  // Number: 14
    ( sp(2.5, -1, 0,  0.5) )  // Number: 15
    ( sp(2.5, -1, 0, -0.5) )  // Number: 16
    ( sp(1.5, -1, 1,  0.5) )  // Number: 17
    ( sp(1.5, -1, 1, -0.5) )  // Number: 18
    ( sp(0.5, -1, 1,  0.5) )  // Number: 19
    ( sp(0.5, -1, 1, -0.5) ); // Number: 20

// FIXME This works, but is not the 'best' place for testing pandya
TEST_F(GMatrix, Pandya) {
    EXPECT_FLOAT_EQ(  4.0043716357, nuclear::pandya( G_pp,
                ph( s[1],  s[2],  1 ),
                ph( s[3],  s[4],  1 ) ) );
    EXPECT_FLOAT_EQ(  6.4905633211, nuclear::pandya( G_pp,
                ph( s[1],  s[2],  0 ),
                ph( s[3],  s[4],  0 ) ) );
    EXPECT_FLOAT_EQ( -5.4604944101, nuclear::pandya( G_pp,
                ph( s[1],  s[3],  1 ),
                ph( s[2],  s[4],  1 ) ) );
    EXPECT_FLOAT_EQ( -0.7428214168, nuclear::pandya( G_pp,
                ph( s[13], s[3],  3 ),
                ph( s[14], s[4],  3 ) ) );
    // Tests from failing PH interaction
    EXPECT_FLOAT_EQ( 0.084600292, nuclear::pandya( G_pp,
                ph( s[13], s[3],  2 ),
                ph( s[13], s[3],  2) ) );
    EXPECT_FLOAT_EQ( -0.8093589423, nuclear::pandya( G_pp,
                ph( s[13], s[3],  2 ),
                ph( s[3],  s[13], 2) ) );
    EXPECT_FLOAT_EQ(  1.4165830802, nuclear::pandya( G_pp,
                ph( s[19], s[3],  2 ),
                ph( s[18], s[6],  2) ) );
}

TEST_F(GMatrix, ParticleHole) {
    boost::tuple< rpa::Modelspace, sp_energies_t > temp
        = rpa::factories::build_modelspace_from_file(
            "data/crazyms.dat", nuclear::nucleus( "Ca48" ) );
    rpa::Modelspace ms = temp.get<0>();

    nuclear::ph_interaction_t Gph =
        rpa::factories::build_ph_gmatrix_from_pp( G_pp, ms );
    nuclear::ph_state ph1( s[1],  s[2],  1 );
    nuclear::ph_state ph2( s[3],  s[4],  1 );
    std::cout << ph1 << " :: " << ph2 << std::endl;
    std::cout << "Parity: " << parity(ph1) << ", " << parity(ph2) << std::endl;
    std::cout << "Isospin: " << tz(ph1) << ", " << tz(ph2) << std::endl;
    EXPECT_FLOAT_EQ(  4.0043716357, Gph( ph( s[1],  s[2],  1 ),
                                         ph( s[3],  s[4],  1 ) ) );
//    EXPECT_FLOAT_EQ(  6.4905633211, nuclear::pandya( G_pp,
//                ph( s[1],  s[2],  0 ),
//                ph( s[3],  s[4],  0 ) ) );
//    EXPECT_FLOAT_EQ( -5.4604944101, nuclear::pandya( G_pp,
//                ph( s[1],  s[3],  1 ),
//                ph( s[2],  s[4],  1 ) ) );
//    EXPECT_FLOAT_EQ( -0.7428214168, nuclear::pandya( G_pp,
//                ph( s[13], s[3],  3 ),
//                ph( s[14], s[4],  3 ) ) );
//    // Tests from failing PH interaction
//    EXPECT_FLOAT_EQ( 0.084600292, nuclear::pandya( G_pp,
//                ph( s[13], s[3],  2 ),
//                ph( s[13], s[3],  2) ) );
//    EXPECT_FLOAT_EQ( -0.8093589423, nuclear::pandya( G_pp,
//                ph( s[13], s[3],  2 ),
//                ph( s[3],  s[13], 2) ) );
//    EXPECT_FLOAT_EQ(  1.4165830802, nuclear::pandya( G_pp,
//                ph( s[19], s[3],  2 ),
//                ph( s[18], s[6],  2) ) );
}

// This test contains numbers that appear exactly in the file.
TEST_F(GMatrix, PPStandard) {
    // A couple that were wrong during development.
    EXPECT_FLOAT_EQ( -0.2344912000, G_pp( pp( s[3],  s[14], 2 ),
                                          pp( s[13], s[4],  2 ) ) );
    EXPECT_FLOAT_EQ( -0.3222744927, G_pp( pp( s[3],  s[14], 3 ),
                                          pp( s[13], s[4],  3 ) ) );
    EXPECT_FLOAT_EQ( -0.810887    , G_pp( pp( s[2],  s[4],  2 ),
                                          pp( s[2],  s[4],  2 ) ) );
    // Some 2+ tests
    EXPECT_FLOAT_EQ( -0.7497247700, G_pp( pp( s[13], s[3],  2 ),
                                          pp( s[3],  s[13], 2 ) ) );
    EXPECT_FLOAT_EQ(  0.8575491513, G_pp( pp( s[1],  s[7],  2 ),
                                          pp( s[5],  s[17], 2 ) ) );
    EXPECT_FLOAT_EQ( -0.3228936712, G_pp( pp( s[5],  s[17], 2 ),
                                          pp( s[9],  s[11], 2 ) ) );
    EXPECT_FLOAT_EQ( -0.0447950562, G_pp( pp( s[7],  s[9],  2 ),
                                          pp( s[7],  s[9],  2 ) ) );
    // Some 2- tests
    EXPECT_FLOAT_EQ(  0.0111985897, G_pp( pp( s[1],  s[15], 2 ),
                                          pp( s[3],  s[7],  2 ) ) );
    EXPECT_FLOAT_EQ( -0.0146918729, G_pp( pp( s[5],  s[7],  2 ),
                                          pp( s[9],  s[17], 2 ) ) );
    EXPECT_FLOAT_EQ( -0.0146918729, G_pp( pp( s[5],  s[7],  2 ),
                                          pp( s[9],  s[17], 2 ) ) );
}

// This test includes normal hits plus phase changes.
TEST_F(GMatrix, PPPhase) {
    // This is a test case where J = 0, so there can be no phase difference.
    // 1 12, 3 4, J = 0
    EXPECT_FLOAT_EQ( -0.3966205694, G_pp( pp( s[1],  s[12], 0 ),
                                          pp( s[3],  s[4],  0 ) ) );
    EXPECT_FLOAT_EQ( -0.3966205694, G_pp( pp( s[12], s[1],  0 ),
                                          pp( s[3],  s[4],  0 ) ) );
    EXPECT_FLOAT_EQ( -0.3966205694, G_pp( pp( s[12], s[1],  0 ),
                                          pp( s[4],  s[3],  0 ) ) );
    EXPECT_FLOAT_EQ( -0.3966205694, G_pp( pp( s[1],  s[12], 0 ),
                                          pp( s[4],  s[3],  0 ) ) );
    // Now J = 1, with just one phase change
    // 3 4, 15 4, J = 1
    EXPECT_FLOAT_EQ( -2.074107922,  G_pp( pp( s[3],  s[4],  1 ),
                                          pp( s[15], s[4],  1 ) ) );
    EXPECT_FLOAT_EQ(  2.074107922,  G_pp( pp( s[4],  s[3],  1 ),
                                          pp( s[15], s[4],  1 ) ) );
    EXPECT_FLOAT_EQ(  2.074107922,  G_pp( pp( s[4],  s[3],  1 ),
                                          pp( s[4],  s[15], 1 ) ) );
    EXPECT_FLOAT_EQ( -2.074107922,  G_pp( pp( s[3],  s[4],  1 ),
                                          pp( s[4],  s[15], 1 ) ) );
    // Now J = 1, Without any phase changes
    // 2 4, 6 10, J = 1
    EXPECT_FLOAT_EQ(  1.979501969,  G_pp( pp( s[2],  s[4],  1 ),
                                          pp( s[6],  s[10], 1 ) ) );
    EXPECT_FLOAT_EQ(  1.979501969,  G_pp( pp( s[4],  s[2],  1 ),
                                          pp( s[6],  s[10], 1 ) ) );
    EXPECT_FLOAT_EQ(  1.979501969,  G_pp( pp( s[4],  s[2],  1 ),
                                          pp( s[10], s[6],  1 ) ) );
    EXPECT_FLOAT_EQ(  1.979501969,  G_pp( pp( s[2],  s[4],  1 ),
                                          pp( s[10], s[6],  1 ) ) );
    // Now J = 1, with two phase changes
    // 1 11, 3 17, J = 1
    EXPECT_FLOAT_EQ(  0.008798113,  G_pp( pp( s[1],  s[11], 1 ),
                                          pp( s[3],  s[17], 1 ) ) );
    EXPECT_FLOAT_EQ( -0.008798113,  G_pp( pp( s[11], s[1],  1 ),
                                          pp( s[3],  s[17], 1 ) ) );
    EXPECT_FLOAT_EQ(  0.008798113,  G_pp( pp( s[11], s[1],  1 ),
                                          pp( s[17], s[3],  1 ) ) );
    EXPECT_FLOAT_EQ( -0.008798113,  G_pp( pp( s[1],  s[11], 1 ),
                                          pp( s[17], s[3],  1 ) ) );
}

// Test for the extra sqrt(2) factor needed based on the mhj file
TEST_F(GMatrix, PPNormalization) {
    // Two factors
    EXPECT_FLOAT_EQ( 2 * 0.7004450235, G_pp( pp( s[11], s[11], 0 ),
                                             pp( s[15], s[15], 0 ) ) );
    EXPECT_FLOAT_EQ( std::sqrt(2.0) * -0.733220442,
                                       G_pp( pp( s[13], s[13], 6 ),
                                             pp( s[13], s[15], 6 ) ) );
}

TEST_F(GMatrix, PPMissedStates) {
    // Angular momentum mismatch
    EXPECT_THROW( G_pp( pp( s[8],  s[3],  2 ),
                        pp( s[3],  s[8],  3 ) ),
                map_key_missing );
    // Parity mismatch
    EXPECT_THROW( G_pp( pp( s[4],   s[8],  1 ),
                        pp( s[10],  s[8],  1 ) ),
                map_key_missing );
    // Isospin mismatch
    EXPECT_THROW( G_pp( pp( s[1],  s[5],  0 ),
                        pp( s[2],  s[6],  0 ) ),
                map_key_missing );
}

// Test for blocks that do not exist.
// NOTE: in principle it could be nice if these returned 0
TEST_F(GMatrix, PPMissedBlocks) {
    sp big_s(6.5, 1, 0, 0.5);
    EXPECT_THROW( G_pp( pp( big_s, big_s, 10 ),
                        pp( big_s, big_s, 10 ) ),
                map_key_missing );
    EXPECT_THROW( G_pp( pp( big_s, big_s,  8 ),
                        pp( big_s, big_s,  8 ) ),
                map_key_missing );
}

TEST_F(GMatrix, PHFunctionalTest) {
    sp_list states( s.begin() + 1, s.end() );
    rpa::Modelspace ms = rpa::factories::build_modelspace_from_states( states,
            nuclear::nucleus( "Ca48" ) );
    nuclear::ph_interaction_t G_ph =
        rpa::factories::build_ph_gmatrix_from_pp( G_pp, ms );

    // These were failing before.
//    EXPECT_FLOAT_EQ( -0.8093589423, G_ph( ph( s[13], s[3],  2 ),
//                                          ph( s[13], s[3],  2 ) ) );
//    EXPECT_FLOAT_EQ( -0.8093589423, G_ph( ph( s[13], s[3],  2 ),
//                                          ph( s[3],  s[13], 2 ) ) );
    EXPECT_FLOAT_EQ(  1.4165830803, G_ph( ph( s[19], s[3],  2 ),
                                          ph( s[18], s[6],  2 ) ) );
    // Test it out on some 2+ states
    // Grab some states

    ph A2 = ph( s[19], s[3],  2 );
    ph B2 = ph( s[18], s[6],  2 );
    ph C2 = ph( s[17], s[3],  2 );
    ph D2 = ph( s[13], s[3],  2 );

    // Some hits, correct values calculated from old Python code
//    EXPECT_FLOAT_EQ(  1.41658308026, G_ph( A2, B2 ) );
    EXPECT_FLOAT_EQ( -1.48240427748, G_ph( C2, B2 ) );
    EXPECT_FLOAT_EQ(  0.40459437798, G_ph( A2, A2 ) );
    EXPECT_FLOAT_EQ(  0.36476255578, G_ph( A2, C2 ) );
    EXPECT_FLOAT_EQ(  0.0846003    , G_ph( D2, D2 ) );

    // Now test some 3- states
    ph A3 = ph( sp(0.5, -1, 1, -0.5), sp(2.5,  1, 0, -0.5), 3 );
    ph B3 = ph( sp(2.5, -1, 0,  0.5), sp(2.5,  1, 0,  0.5), 3 );
    ph C3 = ph( sp(1.5, -1, 1,  0.5), sp(1.5, -1, 0,  0.5), 3 );

    EXPECT_FLOAT_EQ(  1.41658308026, G_ph( A2, B2 ) );
    EXPECT_FLOAT_EQ(  1.41658308026, G_ph( A2, B2 ) );
    EXPECT_FLOAT_EQ(  1.41658308026, G_ph( A2, B2 ) );
    EXPECT_FLOAT_EQ(  1.41658308026, G_ph( A2, B2 ) );

    EXPECT_FLOAT_EQ( -0.260148     , G_ph( A3, B3 ) );
}
