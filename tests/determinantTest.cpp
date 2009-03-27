#include <gtest/gtest.h>

#include <boost/numeric/ublas/matrix.hpp>

#include "determinant.h"

namespace ublas = boost::numeric::ublas;

TEST( Determinant, Simple ) {
    ublas::matrix< int > m1( 2, 2 );
    m1( 0, 0 ) = 3;  m1( 0, 1 ) = 4;
    m1( 1, 0 ) = 1;  m1( 1, 1 ) = 2;

    EXPECT_EQ( 2, determinant( m1 ) );

    ublas::matrix< int > m2( 2, 2 );
    m2( 0, 0 ) = 2;  m2( 0, 1 ) = 4;
    m2( 1, 0 ) = 1;  m2( 1, 1 ) = 2;

    EXPECT_EQ( 0, determinant( m2 ) );
}

TEST( Determinant, Eigenvalues ) {
    ublas::matrix< int > m3( 2, 2 );
    m3( 0, 0 ) = 3;  m3( 0, 1 ) = 4;
    m3( 1, 0 ) = 2;  m3( 1, 1 ) = 1;

    EXPECT_EQ( 0, determinant( m3 - 5 * ublas::identity_matrix<int>(2) ) );
    EXPECT_EQ( 0, determinant( m3 + 1 * ublas::identity_matrix<int>(2) ) );
}
