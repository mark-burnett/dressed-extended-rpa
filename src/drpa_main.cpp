#include <iostream>
#include <fstream>
#include <string>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/cmdline.hpp>

#include "Modelspace.h"
#include "Interaction.h"
#include "MatrixFactory.h"

#include "modelspace_factories.h"
#include "ph_interaction_factories.h"
#include "pp_interaction_factories.h"
#include "term_factories.h"

namespace po = boost::program_options;

int main( int argc, char *argv[] ) {

    // Options parsed by the command line
    po::options_description cmdline_desc("RPA program options");
    cmdline_desc.add_options()
        ("help", "This help message.")
        ("config", po::value<std::string>(), "Name of the configuration file.");
    po::variables_map cmdline_vm;
    po::store( po::parse_command_line( argc, argv, cmdline_desc ), cmdline_vm );
    po::notify( cmdline_vm );

    // Make sure we got the options we wanted
    if ( cmdline_vm.count("help") ) {
        std::cout << cmdline_desc << std::endl;
        return 1;
    }
    if ( !cmdline_vm.count("config") ) {
        std::cerr << "Configuration file not specified.\n";
        return 1;
    }

    std::cout << "Reading configuration from '"
        << cmdline_vm["config"].as<std::string>() << "'." << std::endl;

    // Options stored in the configuration file
    po::options_description config_desc("Configuration file options.");
    config_desc.add_options()
        ("interaction_file", po::value<std::string>(), "Interaction filename.")
        ("modelspace_file",  po::value<std::string>(), "Modelspace filename.")
        ("output_file",      po::value<std::string>(), "Full output filename.");
    po::variables_map config_vm;
    {
        std::ifstream cfile(cmdline_vm["config"].as<std::string>().c_str());
        po::store( po::parse_config_file( cfile,
                    config_desc ), config_vm );
    }
    po::notify( config_vm );

    // Make sure we got the options we needed
    if ( !config_vm.count("interaction_file") ) {
        std::cerr << "Interaction file not specified in config file.\n";
        return 1;
    }
    if ( !config_vm.count("modelspace_file") ) {
        std::cerr << "Modelspace file not specified in config file.\n";
        return 1;
    }
    if ( !config_vm.count("output_file") ) {
        std::cerr << "Output file not specified in config file.\n";
        return 1;
    }

    // Instantiate the object graph for the calculation.
    // Modelspaces
    std::cout << "Building modelspaces." << std::endl;
    SingleParticleModelspace spms =
            read_sp_modelspace_from_file(
                config_vm["modelspace_file"].as<std::string>() );
    ParticleHoleModelspace phms = build_ph_modelspace_from_sp( spms );
    /*
    ParticleParticleModelspace ppms   = build_pp_modelspace_from_sp( spms );
    ParticleParticleModelspace hhms   = build_hh_modelspace_from_sp( spms );
    PPFromSPModelspace         ppspms = build_ppsp_modelspace_from_sp( spms );
    PPFromSPModelspace         hhspms = build_hhsp_modelspace_from_sp( spms );
    */
    std::cout << "Modelspaces built.  PH modelspace sizes:" << std::endl;
    print_ph_modelspace_sizes( std::cout, phms );

    // Particle-hole interaction
    std::cout << "Building interaction objects." << std::endl;
    PPInteraction Gpp = build_gmatrix_from_mhj_file(
            config_vm["interaction_file"].as<std::string>(), spms );
    PHInteraction Gph = build_ph_interaction_from_pp( Gpp, spms );
    std::cout << "Finished building interactions." << std::endl;

    // Build terms
    // loop/build matricies

    /*
    // Actual calculation
    std::cout << "Building RPA calculation object." << std::endl;
    rpa::Calculation myRPA = rpa::factories::build_standard_rpa(
            Gph, ms, sp_energies );
    std::cout << "RPA calculation object built." << std::endl;

    // Perform the calculation
    typedef std::pair< Jpi_t, boost::tuple< util::vector_t,
                                            util::matrix_t > > result_t;
    std::list< result_t > results = myRPA.execute( 0 );

    // Write output
    std::cout << "Calculation complete.  Writing results." << std::endl;
    {
        std::ofstream outfile(
                config_vm["output_file"].as<std::string>().c_str() );

        BOOST_FOREACH( const result_t &result, results ) {
            const Jpi_t &Jpi = result.first;
            const boost::tuple< util::vector_t, util::matrix_t > &r =
                result.second;
            outfile << "(" << Jpi.first << ", " << Jpi.second << ")\n"
                    << r.get<0>() << "\n"
                    << r.get<1>() << "\n\n";
        }
    }
    */

    return 0;
}
