#include <gtest/gtest.h>

#include <cmath>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include "exceptions.h"
#include "find_root.h"

using namespace std;
using namespace util;
using namespace boost::lambda;

TEST(FalsePosition, Success) {
    EXPECT_FLOAT_EQ( 0.5671433083,
            false_position( _1 * bind((double(*)(double)) exp, _1) - 1,
                -1, 1, 1e-10, 200 ) );
    EXPECT_FLOAT_EQ( 0.8041330975,
            false_position( 11 * bind((double(*)(double,int)) pow, _1, 11) - 1,
                 0, 1, 1e-10, 200 ) );
    EXPECT_FLOAT_EQ( 2.0945515532,
            false_position( bind((double(*)(double,int)) pow, _1, 3) - 2*_1 - 5,
                 2, 3, 1e-10, 200 ) );
}

TEST(FalsePosition, OutOfBounds) {
    EXPECT_THROW( false_position( _1 * bind((double(*)(double)) exp, _1) - 1,
                                  -1, 0, 1e-10, 200 ), root_finding_error );
}
