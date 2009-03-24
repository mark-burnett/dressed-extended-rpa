#include <vector>
#include <complex>
#include <gtest/gtest.h>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "linalg.h"

#include "Modelspace.h"
#include "Interaction.h"
#include "Term.h"
#include "MatrixFactory.h"

#include "modelspace_factories.h"
#include "pp_interaction_factories.h"
#include "ph_interaction_factories.h"
#include "term_factories.h"


double real( const util::complex_t &a ) {
    return std::real( a );
}

TEST( DRPA, FunctionalTest ) {
    // Build all the crap we use:
    SingleParticleModelspace spms
        = read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    ParticleHoleModelspace phms = build_ph_modelspace_from_sp( spms );
    PPInteraction Gpp
        = build_gmatrix_from_mhj_file( "tests/data/test_interaction.mhj", spms);
    PHInteraction Gph = build_ph_interaction_from_pp( Gpp, spms );
    std::vector< Term > terms = build_rpa_terms( Gph, spms );

    // Build some matricies and check some eigenvalues

    // 0-
    {
        int tz     =  0;
        int parity = -1;
        int J      =  0;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 16.819696, real_vals[ vals.size() / 2 ] ) << "0-";
    }

    // 0+
    {
        int tz     =  0;
        int parity =  1;
        int J      =  0;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 27.696529, real_vals[ vals.size() / 2 ] ) << "0+";
    }
    // 1-
    {
        int tz     =  0;
        int parity = -1;
        int J      =  1;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 14.122624, real_vals[ vals.size() / 2 + 1 ] ) << "1-";
    }
    // 1+
    {
        int tz     =  0;
        int parity =  1;
        int J      =  1;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ(  5.769120, real_vals[ vals.size() / 2 ] ) << "1+";
    }
    // 2-
    {
        int tz     =  0;
        int parity = -1;
        int J      =  2;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 16.382633, real_vals[ vals.size() / 2 ] ) << "2-";
    }

    // 2+
    {
        int tz     =  0;
        int parity =  1;
        int J      =  2;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 1.5833673, real_vals[ vals.size() / 2 ] ) << "2+";
        EXPECT_FLOAT_EQ( 4.9466329, real_vals[ vals.size() / 2 + 1 ] ) << "2+";
    }
    // 3-
    {
        int tz     =  0;
        int parity = -1;
        int J      =  3;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 11.084711, real_vals[ vals.size() / 2 ] ) << "3-";
    }
    // 3+
    {
        int tz     =  0;
        int parity =  1;
        int J      =  3;
        util::cvector_t vals = util::eigenvalues(
              util::cmatrix_t(build_static_matrix( terms,
                      phms[tz+1][(parity+1)/2][J] ) ) );
        std::vector< double > real_vals( vals.size() );
        real_vals.resize( vals.size() );
        std::transform( vals.begin(), vals.end(),
                        real_vals.begin(), real );

        std::sort( real_vals.begin(), real_vals.end() );

        EXPECT_FLOAT_EQ( 2.098350, real_vals[ vals.size() / 2 ] ) << "3+";
    }
}
