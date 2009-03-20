#include <boost/bind.hpp>

#include "linalg.h"

#include "Modelspace.h"
#include "Interaction.h"

#include "Term.h"
#include "non_interacting.h"

#include "exceptions.h"

namespace terms {

util::matrix_t
non_interacting( const std::vector< ParticleHoleState > &vec, double E,
                 position_t pos, const SingleParticleModelspace &spms ) {
    int size = vec.size();
    util::matrix_t m( size, size );
    m.clear();

    if ( ENUM_B == pos || ENUM_B_STAR == pos )
        return m;

    for ( int i=0; i < size; ++i ) {
        int ip  = vec[i].ip;
        int ih  = vec[i].ih;
        int ipf = vec[i].ipf;
        int ihf = vec[i].ihf;
        switch ( pos ) {
            case ENUM_A:
            case ENUM_A_STAR:
                    m( i, i ) = spms.pfrag[ip][ipf].E - spms.hfrag[ih][ihf].E;
                break;
            default:
                throw invalid_matrix_position(); } }
    return m;
    // Dummy code for unused E
    ++E;
}

Term make_non_interacting( const SingleParticleModelspace &spms ) {
    return boost::bind( non_interacting, _1, _2, _3, spms );
}

} // end namespace terms
