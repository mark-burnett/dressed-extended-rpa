AUTOMAKE_OPTIONS = subdir-objects

INCLUDES = $(PROJECT_INCLUDES)

# --- Convenience libraries ---
noinst_LTLIBRARIES     = src/liberpa.la
src_liberpa_la_SOURCES = src/angular_momentum.cpp src/io.cpp src/pandya.cpp\
						 src/Modelspace.cpp\
						 src/modelspace_factories.cpp\
						 src/interaction_factories.cpp\
						 src/MatrixFactory.cpp\
						 src/terms/non_interacting.cpp\
						 src/terms/first_order.cpp\
						 src/terms/screening.cpp\
						 src/terms/ladder.cpp\
						 src/terms/self_energy.cpp\
						 src/term_factories.cpp

# --- Main programs ---


# --- Tests ---
check_PROGRAMS   = bin/test
bin_test_SOURCES = tests/test_main.cpp\
				   tests/modelspace_factoriesTest.cpp\
				   tests/interaction_factoriesTest.cpp\
				   tests/pandyaTest.cpp
bin_test_LDADD   = src/liberpa.la
LIBS        = $(GTEST_LIBS)

TESTS = bin/test

# --- Cleanup ---
MAINTAINERCLEANFILES = configure config.log config.status\
					   aclocal.m4 autom4te.cache Makefile.in\
					   config/config.guess config/config.sub config/depcomp\
					   config/install-sh config/ltmain.sh config/missing
