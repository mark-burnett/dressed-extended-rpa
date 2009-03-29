#include <gtest/gtest.h>

#include <vector>

#include "Modelspace.h"
#include "modelspace_factories.h"

TEST( Modelspace, Poles ) {
    // Build the modelspaces
    SingleParticleModelspace spms =
        read_sp_modelspace_from_file( "tests/data/ipm_modelspace.dat" );
    ParticleHoleModelspace phms = build_ph_modelspace_from_sp( spms );

    // Just check the 0+
    std::vector< double > poles = get_ph_poles( 0, 1, 0, phms, spms );

    EXPECT_EQ( 4, poles.size() );
    EXPECT_FLOAT_EQ( 40.854398, poles[0] );
    EXPECT_FLOAT_EQ( 40.935378, poles[1] );
    EXPECT_FLOAT_EQ( 42.117380, poles[2] );
    EXPECT_FLOAT_EQ( 42.420139, poles[3] );
}
