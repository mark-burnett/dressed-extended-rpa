/* Unit-testing for the simple angular momentum functions.
 *
 * Mark Burnett, November 2008
 */
#include <gtest/gtest.h>

#include "exceptions.h"
#include "angular_momentum.h"

TEST(AngularMomentum, Predicates) {
    // is_integer
    EXPECT_TRUE(  is_integer( 3.0) );
    EXPECT_TRUE(  is_integer(-1.0) );
    EXPECT_TRUE(  is_integer( 0.0) );

    EXPECT_FALSE( is_integer( 0.5) );
    EXPECT_FALSE( is_integer(-1.1) );
    EXPECT_FALSE( is_integer( 7.2) );

    // is_half_integer
    EXPECT_TRUE(  is_half_integer( 3.0) );
    EXPECT_TRUE(  is_half_integer(-1.5) );
    EXPECT_TRUE(  is_half_integer( 0.0) );

    EXPECT_FALSE( is_half_integer( 0.6) );
    EXPECT_FALSE( is_half_integer(-1.1) );
    EXPECT_FALSE( is_half_integer( 7.2) );

    // is_strict_half_integer
    EXPECT_TRUE(  is_strict_half_integer( 3.5) );
    EXPECT_TRUE(  is_strict_half_integer(-1.5) );
    EXPECT_TRUE(  is_strict_half_integer( 0.5) );

    EXPECT_FALSE( is_strict_half_integer( 0.6) );
    EXPECT_FALSE( is_strict_half_integer(-1.1) );
    EXPECT_FALSE( is_strict_half_integer( 7.2) );

    EXPECT_FALSE( is_strict_half_integer( 3.0) );
    EXPECT_FALSE( is_strict_half_integer(-1.0) );
    EXPECT_FALSE( is_strict_half_integer( 0.0) );

    // is_triangular
    EXPECT_TRUE(  is_triangular(  1,    1,    0   ) );
    EXPECT_TRUE(  is_triangular(  1,    1,    1   ) );
    EXPECT_TRUE(  is_triangular(  1,    1,    2   ) );

    EXPECT_TRUE(  is_triangular(  0.5,  1.5,  1   ) );
    EXPECT_TRUE(  is_triangular(  1,    2.5,  1.5 ) );
    EXPECT_TRUE(  is_triangular(  0,    0,    0   ) );

    EXPECT_FALSE( is_triangular(  1,    1,   -1   ) );
    EXPECT_FALSE( is_triangular(  1,    1,    3   ) );

    EXPECT_FALSE( is_triangular(  0.5,  1.5,  0.5 ) );
    EXPECT_FALSE( is_triangular(  0.5,  1,    0   ) );
    EXPECT_FALSE( is_triangular(  0,    1.5,  1   ) );

    EXPECT_FALSE( is_triangular(  1.5, -1.5,  0   ) );
    EXPECT_FALSE( is_triangular( -1.5,  1.5,  0   ) );
    EXPECT_FALSE( is_triangular(  1.5,  0,   -1.5 ) );
}

TEST(AngularMomentum, Wigner3J) {
    // Legal evaluations
    // We use EXPECT_FLOAT_EQ to give a slightly wider tolerance
    EXPECT_FLOAT_EQ( 0.186989398002, wigner3j( 6,    4,    2,
                                                        0,    0,    0));
    EXPECT_FLOAT_EQ( 0.387298334621, wigner3j( 2,    0.5,  1.5,
                                                        1,   -0.5, -0.5));
    EXPECT_FLOAT_EQ(-0.218217890236, wigner3j( 1.5,  2.5,  4,
                                                       -0.5,  0.5,  0));
    EXPECT_FLOAT_EQ(-0.4472135955,   wigner3j( 2,    2,    0,
                                                        1,   -1,    0));
    EXPECT_FLOAT_EQ( 0.36514837167,  wigner3j( 2,    1,    1,
                                                        0,    0,    0));
    EXPECT_EQ(       0.5,            wigner3j( 1.5,  0,    1.5,
                                                        0.5,  0,   -0.5));

    // Non-triangular J's (return 0)
    EXPECT_EQ(0, wigner3j( 1,  1, 7, 0,  0,  0));

    // Illegal J's (not half integers, or negative)
    EXPECT_THROW(wigner3j( 0.1,  1, 1, 1, 0, 1 ),
                 illegal_angular_momentum);
    EXPECT_THROW(wigner3j( 1,   -1, 0, 0, 0, 0 ),
                 illegal_angular_momentum);

    // Illegal m's (do not match their j's)
    EXPECT_THROW(wigner3j( 6,     4,   2,    0.5,  0,    0   ),
                 illegal_angular_momentum);
    EXPECT_THROW(wigner3j(  2,    1,    2,   3,   -1,   -2   ),
                 illegal_angular_momentum);
    EXPECT_THROW(wigner3j(  1.5,  2.5,  4,   0,    0.5, -0.5 ),
                 illegal_angular_momentum);
    EXPECT_THROW(wigner3j( 6,     4,   2,   -1,    1,    1   ),
                 illegal_angular_momentum);
}

TEST(AngularMomentum, ClebschGordan) {
    // Illegal calls
    EXPECT_THROW(clebsch_gordan(1, 0.5, 1, 0, 0.5, 0.5),
                 illegal_angular_momentum);
    EXPECT_THROW(clebsch_gordan(1, 1,   1, 0, 0,   0),
                 illegal_angular_momentum);

    // Legal calls
    EXPECT_FLOAT_EQ(0.57735026919,
            clebsch_gordan( 0.5,  0.5, 1,    0,   0.5,  0.5));
    EXPECT_FLOAT_EQ(0.816496580928,
            clebsch_gordan( 0.5,  0.5, 1,    0,   1.5,  0.5));
    EXPECT_FLOAT_EQ(0.774596669241,
            clebsch_gordan( 2,    1,   0.5, -0.5, 1.5,  0.5));
    EXPECT_FLOAT_EQ(0.894427191,
            clebsch_gordan( 2,    1,   0.5,  0.5, 2.5,  1.5));
}

TEST(AngularMomentum, Wigner6J) {
    // Illegal cals
    EXPECT_THROW(wigner6j(1, 2,   0.4, 1, 2, 0.5),
            illegal_angular_momentum);
    EXPECT_THROW(wigner6j(1, 1.9, 0.5, 1, 2, 0.5),
            illegal_angular_momentum);

    // Silent 0 returns
    EXPECT_EQ(0, wigner6j(1,    1,   3, 1,   1,   2 ));
    EXPECT_EQ(0, wigner6j(0.5,  0.5, 1, 0.5, 4.5, 1 ));
    EXPECT_EQ(0, wigner6j(1,   -1,   0, 1,   1,   0 ));

    // Normal evaluations
    EXPECT_FLOAT_EQ( 0.0436435780472,
            wigner6j( 1,   2,   3,   2,   1,   2   ));
    EXPECT_FLOAT_EQ( 0.166666666667,
            wigner6j( 1,   1,   1,   1,   2,   1   ));
    EXPECT_FLOAT_EQ(-0.333333333333,
            wigner6j( 1,   1,   1,   0,   1,   1   ));
    EXPECT_FLOAT_EQ( 0.142857142857,
            wigner6j( 1,   2,   3,   5,   2,   3   ));
    EXPECT_FLOAT_EQ( 0.0690065559342,
            wigner6j( 1.5, 3.5, 2,   3.5, 2.5, 1   ));
}
