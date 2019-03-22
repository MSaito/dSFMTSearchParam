noinst_PROGRAMS = dSFMTdc calc_parity calc_equidist calc_fixpoint \
test_linearity

dSFMTdc_SOURCES = dSFMTsearch.hpp dSFMTdc.cpp Annihilate.h Annihilate.cpp \
AlgorithmDSFMTEquidistribution.hpp w128.hpp calc_fixpoint.cpp mpicontrol.hpp
calc_equidist_SOURCES = dSFMTsearch.hpp calc_equidist.cpp Annihilate.cpp \
AlgorithmDSFMTEquidistribution.hpp
test_linearity_SOURCES = dSFMTsearch.hpp test_linearity.cpp Annihilate.cpp
calc_fixpoint_SOURCES = dSFMTsearch.hpp calc_fixpoint.cpp Annihilate.cpp \
printBinary.cpp printBinary.h calc_fixpoint.h calc_fixpoint_main.cpp
calc_parity_SOURCES = dSFMTsearch.hpp calc_parity.cpp Annihilate.cpp \
printBinary.cpp printBinary.h calc_fixpoint.h calc_parity.h \
calc_parity_main.cpp

dSFMTdc_CXXFLAGS = -DNO_MPI
LDADD = -lMTToolBox
AM_CXXFLAGS = -Wall -O2 -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
EXTRA_DIST = memo.txt
