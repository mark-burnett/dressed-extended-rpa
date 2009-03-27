#include <gtest/gtest.h>

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/assign/list_of.hpp>

#include "search.h"

TEST( Search, SplitValues ) {
    typedef boost::tuple< int, int > tup_t;
    std::vector< double > vals = boost::assign::list_of
        (-3.0) (-2.0) (-1.0) (0.1) (1.0) (2.0) (2.9);

    EXPECT_EQ( tup_t( 4, 3 ), split_values( vals,  0   ) );
    EXPECT_EQ( tup_t( 2, 4 ), split_values( vals,  1   ) );
    EXPECT_EQ( tup_t( 7, 0 ), split_values( vals, -5   ) );
    EXPECT_EQ( tup_t( 6, 0 ), split_values( vals, -3   ) );
    EXPECT_EQ( tup_t( 0, 7 ), split_values( vals,  3   ) );
    EXPECT_EQ( tup_t( 0, 6 ), split_values( vals,  2.9 ) );
}

TEST( Search, GetNumSolutions ) {
    std::vector< double > left_vals  = boost::assign::list_of
        (3) (7) (10);
    std::vector< double > right_vals = boost::assign::list_of
        (1) (3) (8);

    EXPECT_EQ( 1, get_num_solutions( left_vals,   4,
                                     right_vals,  5 ) );
    EXPECT_EQ( 1, get_num_solutions( left_vals,   0,
                                     right_vals,  2 ) );
    EXPECT_EQ( 2, get_num_solutions( left_vals,   0,
                                     right_vals,  4 ) );
    EXPECT_EQ( 0, get_num_solutions( left_vals,  -1,
                                     right_vals,  0 ) );
    EXPECT_EQ( 0, get_num_solutions( left_vals,  11,
                                     right_vals, 12 ) );
}
