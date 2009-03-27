#include <gtest/gtest.h>

#include <vector>
#include <boost/assign/list_of.hpp>

#include "intervals.h"

TEST( Intervals, BuildRootIntervals ) {
    std::vector< double > left = boost::assign::list_of
        (1) (2) (6) (9);
    std::vector< double > right = boost::assign::list_of
        (4) (7) (8) (12);

    std::vector< interval_t > root = build_root_intervals( left, right );

    using boost::numeric::equal;
    EXPECT_EQ( 4, root.size() );
    EXPECT_TRUE( equal( interval_t( 1,  4 ), root[0] ) );
    EXPECT_TRUE( equal( interval_t( 2,  7 ), root[1] ) );
    EXPECT_TRUE( equal( interval_t( 6,  8 ), root[2] ) );
    EXPECT_TRUE( equal( interval_t( 9, 12 ), root[3] ) );
}

TEST( Intervals, GetSubIntervals ) {
    std::vector< interval_t > root;
    root.push_back( interval_t( 1, 4 ) );
    root.push_back( interval_t( 2, 7 ) );
    root.push_back( interval_t( 6, 8 ) );
    root.push_back( interval_t( 9, 12 ) );

    std::vector< interval_t > sub = get_sub_intervals( root );
    EXPECT_EQ( 7, sub.size() );

    using boost::numeric::equal;
    EXPECT_TRUE( equal( interval_t( 1,  2 ), sub[0] ) );
    EXPECT_TRUE( equal( interval_t( 2,  4 ), sub[1] ) );
    EXPECT_TRUE( equal( interval_t( 4,  6 ), sub[2] ) );
    EXPECT_TRUE( equal( interval_t( 6,  7 ), sub[3] ) );
    EXPECT_TRUE( equal( interval_t( 7,  8 ), sub[4] ) );
    EXPECT_TRUE( equal( interval_t( 8,  9 ), sub[5] ) );
    EXPECT_TRUE( equal( interval_t( 9, 12 ), sub[6] ) );
}

TEST( Intervals, GetProbabilities ) {
    std::vector< interval_t > root;
    root.push_back( interval_t( 1, 4 ) );
    root.push_back( interval_t( 2, 7 ) );
    root.push_back( interval_t( 6, 8 ) );
    root.push_back( interval_t( 9, 12 ) );

    std::vector< interval_t > sub = get_sub_intervals( root );
    std::vector< double > prob = get_sub_interval_probabilities( root, sub );

    EXPECT_DOUBLE_EQ( 1 * ( 1.0/3 ),         prob[0] );
    EXPECT_DOUBLE_EQ( 2 * ( 1.0/3 + 1.0/5 ), prob[1] );
    EXPECT_DOUBLE_EQ( 2 * ( 1.0/5 ),         prob[2] );
    EXPECT_DOUBLE_EQ( 1 * ( 1.0/5 + 1.0/2 ), prob[3] );
    EXPECT_DOUBLE_EQ( 1 * ( 1.0/2 ),         prob[4] );
    EXPECT_DOUBLE_EQ( 1 * ( 0 ),             prob[5] );
    EXPECT_DOUBLE_EQ( 3 * ( 1.0/3 ),         prob[6] );
}
