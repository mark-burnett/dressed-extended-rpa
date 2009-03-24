#ifndef _TERM_FACTORIES_H_
#define _TERM_FACTORIES_H_

#include "Interaction.h"
#include "Modelspace.h"
#include "Term.h"

std::vector< Term > build_rpa_terms( const PHInteraction &Gph,
                                     const SingleParticleModelspace &spms );

std::vector< Term > build_static_erpa_terms(
                                    const PHInteraction &Gph,
                                    const PPInteraction &Gpp,
                                    const ParticleHoleModelspace     &phms,
                                    const ParticleParticleModelspace &ppms,
                                    const ParticleParticleModelspace &hhms,
                                    const SingleParticleModelspace   &spms );

std::vector< Term > build_dynamic_derpa_terms(
                                    const PHInteraction &Gph,
                                    const PPInteraction &Gpp,
                                    const ParticleHoleModelspace     &phms,
                                    const ParticleParticleModelspace &ppms,
                                    const ParticleParticleModelspace &hhms,
//                                    const PPFromSPModelspace         &ppspms,
//                                    const PPFromSPModelspace         &hhspms,
                                    const SingleParticleModelspace   &spms );

#endif // _TERM_FACTORIES_H_
