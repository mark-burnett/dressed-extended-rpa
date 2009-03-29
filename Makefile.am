AUTOMAKE_OPTIONS = subdir-objects

INCLUDES = $(PROJECT_INCLUDES) $(BOOST_NUMERIC_BINDINGS_CFLAGS)

# --- Convenience libraries ---
noinst_LTLIBRARIES      = src/libderpa.la
src_libderpa_la_SOURCES = src/angular_momentum.cpp src/io.cpp src/pandya.cpp\
						  src/linalg.cpp\
						  src/find_root.cpp\
						  src/Modelspace.cpp\
						  src/modelspace_factories.cpp\
						  src/pp_interaction_factories.cpp\
						  src/ph_interaction_factories.cpp\
						  src/MatrixFactory.cpp\
						  src/intervals.cpp\
						  src/search.cpp\
						  src/terms/non_interacting.cpp\
						  src/terms/first_order.cpp\
						  src/terms/screening.cpp\
						  src/terms/ladder.cpp\
						  src/terms/self_energy.cpp\
						  src/term_factories.cpp

# --- Main programs ---
bin_PROGRAMS     = bin/drpa bin/erpa bin/plot_eigenvalues

bin_drpa_SOURCES = src/drpa_main.cpp
bin_drpa_LDADD   = src/libderpa.la

bin_erpa_SOURCES = src/erpa_main.cpp
bin_erpa_LDADD   = src/libderpa.la

bin_plot_eigenvalues_SOURCES = src/plot_eigenvalues.cpp
bin_plot_eigenvalues_LDADD   = src/libderpa.la


# --- Tests ---
check_PROGRAMS   = bin/test
bin_test_SOURCES = tests/test_main.cpp\
				   tests/angular_momentumTest.cpp\
				   tests/find_rootTest.cpp\
				   tests/modelspace_factoriesTest.cpp\
				   tests/ModelspaceTest.cpp\
				   tests/pp_interaction_factoriesTest.cpp\
				   tests/ph_interaction_factoriesTest.cpp\
				   tests/pandyaTest.cpp\
				   tests/drpaTest.cpp\
				   tests/determinantTest.cpp\
				   tests/intervalsTest.cpp\
				   tests/searchTest.cpp
bin_test_LDADD   = src/libderpa.la
LIBS             = $(GTEST_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(LAPACK_LIBS)

TESTS = bin/test

# --- Cleanup ---
MAINTAINERCLEANFILES = configure config.log config.status\
					   aclocal.m4 autom4te.cache Makefile.in\
					   config/config.guess config/config.sub config/depcomp\
					   config/install-sh config/ltmain.sh config/missing
