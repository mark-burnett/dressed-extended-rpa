#include <vector>

#include "Interaction.h"
#include "Modelspace.h"
#include "Term.h"

#include "terms/non_interacting.h"
#include "terms/first_order.h"

#include "terms/screening.h"
#include "terms/ladder.h"
#include "terms/self_energy.h"

std::vector< Term > build_rpa_terms( const PHInteraction &Gph,
                                     const SingleParticleModelspace &spms ) {
    std::vector< Term > tvec;

    tvec.push_back( terms::make_non_interacting( spms ) );
    tvec.push_back( terms::make_first_order( Gph, spms ) );

    return tvec;
}

std::vector< Term > build_dynamic_erpa_terms(
                                    const PHInteraction &Gph,
                                    const PPInteraction &Gpp,
                                    const ParticleHoleModelspace     &phms,
                                    const ParticleParticleModelspace &ppms,
                                    const ParticleParticleModelspace &hhms,
                                    const SEModelspace               &sems,
                                    const SingleParticleModelspace   &spms ) {
    std::vector< Term > tvec;

    tvec.push_back( terms::make_screening( Gph, phms, spms ) );
    tvec.push_back( terms::make_ladder( Gpp, ppms, hhms, spms ) );
    tvec.push_back( terms::make_self_energy( Gpp, ppms, hhms, sems, spms ) );
    return tvec;
    // Dummy code
    std::vector< Term > fvec;
    fvec.push_back( terms::make_screening( Gph, phms, spms ) );
    fvec.push_back( terms::make_ladder( Gpp, ppms, hhms, spms ) );
    fvec.push_back( terms::make_self_energy( Gpp, ppms, hhms, sems, spms ) );
}

std::vector< Term > build_dynamic_derpa_terms(
                                    const PHInteraction &Gph,
                                    const PPInteraction &Gpp,
                                    const ParticleHoleModelspace     &phms,
                                    const ParticleParticleModelspace &ppms,
                                    const ParticleParticleModelspace &hhms,
//                                    const PPFromSPModelspace         &ppspms,
//                                    const PPFromSPModelspace         &hhspms,
                                    const SingleParticleModelspace   &spms ) {
    std::vector< Term > tvec;

    tvec.push_back( terms::make_screening( Gph, phms, spms ) );
    tvec.push_back( terms::make_ladder( Gpp, ppms, hhms, spms ) );
// NOTE: no self energy is used for DERPA -- it's already in the "dressing"
//    tvec.push_back( terms::make_self_energy( Gpp, ppms, hhms,
//                                             ppspms, hhspms, spms ) );
    return tvec;
}
