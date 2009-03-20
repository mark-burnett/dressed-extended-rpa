#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include <boost/function.hpp>

#include "Modelspace.h"

typedef boost::function< double ( const ParticleParticleState &,
                                  const ParticleParticleState & )
> PPInteraction;

typedef boost::function< double ( const ParticleHoleState &,
                                  const ParticleHoleState & )
> PHInteraction;

#endif // _INTERACTION_H_
