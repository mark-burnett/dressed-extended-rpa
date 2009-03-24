AUTOMAKE_OPTIONS = subdir-objects

INCLUDES = $(PROJECT_INCLUDES) $(BOOST_NUMERIC_BINDINGS_CFLAGS)

# --- Convenience libraries ---
noinst_LTLIBRARIES     = src/liberpa.la
src_liberpa_la_SOURCES = src/angular_momentum.cpp src/io.cpp src/pandya.cpp\
						 src/linalg.cpp\
						 src/Modelspace.cpp\
						 src/modelspace_factories.cpp\
						 src/pp_interaction_factories.cpp\
						 src/ph_interaction_factories.cpp\
						 src/MatrixFactory.cpp\
						 src/terms/non_interacting.cpp\
						 src/terms/first_order.cpp\
						 src/terms/screening.cpp\
						 src/terms/ladder.cpp\
						 src/terms/self_energy.cpp\
						 src/term_factories.cpp

# --- Main programs ---
bin_PROGRAMS     = bin/drpa
bin_drpa_SOURCES = src/drpa_main.cpp
bin_drpa_LDADD   = src/liberpa.la


# --- Tests ---
check_PROGRAMS   = bin/test
bin_test_SOURCES = tests/test_main.cpp\
				   tests/angular_momentumTest.cpp\
				   tests/modelspace_factoriesTest.cpp\
				   tests/pp_interaction_factoriesTest.cpp\
				   tests/ph_interaction_factoriesTest.cpp\
				   tests/pandyaTest.cpp\
				   tests/drpaTest.cpp
bin_test_LDADD   = src/liberpa.la
LIBS             = $(GTEST_LIBS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(LAPACK_LIBS)

TESTS = bin/test

# --- Cleanup ---
MAINTAINERCLEANFILES = configure config.log config.status\
					   aclocal.m4 autom4te.cache Makefile.in\
					   config/config.guess config/config.sub config/depcomp\
					   config/install-sh config/ltmain.sh config/missing
