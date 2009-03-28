#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "linalg.h"
#include "Term.h"
#include "MatrixFactory.h"

util::matrix_t
MatrixFactory::build( double E ) const {
    // Copy the static elements
    util::matrix_t result( static_matrix );

    // Generate the dynamic elements
    int size = result.size1() / 2;
    ublas::range first_half( 0, size );
    ublas::range second_half( size, 2*size );

    typedef ublas::matrix_range< util::matrix_t > submatrix_t;
    submatrix_t A     ( result, first_half,  first_half );
    submatrix_t A_star( result, second_half, second_half );

    BOOST_FOREACH( const Term &t, dynamic_terms ) {
        A      += t( ph_states, E, ENUM_A );
        A_star -= t( ph_states, E, ENUM_A_STAR );
    }
    return result;
}

util::matrix_t
build_static_rpa_matrix( const std::vector< Term > &terms,
                         const std::vector< ParticleHoleState > &ph_states ) {
    int size = ph_states.size();

    util::matrix_t m( 2 * size, 2 * size );
    m.clear();

    typedef ublas::matrix_range< util::matrix_t > submatrix_t;
    ublas::range first_half( 0, size );
    ublas::range second_half( size, 2*size );

    submatrix_t A     ( m, first_half,  first_half );
    submatrix_t A_star( m, second_half, second_half );
    submatrix_t B     ( m, second_half, first_half );
    submatrix_t B_star( m, first_half,  second_half );

    BOOST_FOREACH( const Term &t, terms ) {
        A      += t( ph_states, 0, ENUM_A );
        B      += t( ph_states, 0, ENUM_B );
        A_star -= t( ph_states, 0, ENUM_A_STAR );
        B_star -= t( ph_states, 0, ENUM_B_STAR ); }
    return m; }

util::matrix_t
build_static_erpa_matrix( const std::vector< Term > &static_terms,
                          const std::vector< Term > &dynamic_terms,
                          const std::vector< ParticleHoleState > &ph_states ) {
    int size = ph_states.size();

    util::matrix_t m( 2 * size, 2 * size );
    m.clear();

    typedef ublas::matrix_range< util::matrix_t > submatrix_t;
    ublas::range first_half( 0, size );
    ublas::range second_half( size, 2*size );

    submatrix_t A     ( m, first_half,  first_half );
    submatrix_t A_star( m, second_half, second_half );
    submatrix_t B     ( m, second_half, first_half );
    submatrix_t B_star( m, first_half,  second_half );

    BOOST_FOREACH( const Term &t, static_terms ) {
        A      += t( ph_states, 0, ENUM_A );
        B      += t( ph_states, 0, ENUM_B );
        A_star -= t( ph_states, 0, ENUM_A_STAR );
        B_star -= t( ph_states, 0, ENUM_B_STAR ); }

    BOOST_FOREACH( const Term &t, dynamic_terms ) {
        B      += t( ph_states, 0, ENUM_B );
        B_star -= t( ph_states, 0, ENUM_B_STAR ); }
    return m; }
