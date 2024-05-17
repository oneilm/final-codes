# Makefile for final-codes
#
# This is the only makefile; there are no makefiles in subdirectories. Users
# should not need to edit this makefile (doing so would make it hard to stay up
# to date with repo version). Rather in order to change OS/environment-specific
# compilers and flags, create the file make.inc, which overrides the defaults
# below (which are for ubunutu linux/gcc system).


# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran


# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy
CFLAGS= -fPIC -O3 -march=native -funroll-loops -std=c99
CXXFLAGS= -std=c++11 -fPIC -O3 -march=native -funroll-loops


# For your OS, override the above by placing make variables in make.inc
-include make.inc

.PHONY: usage hermexps clean

default: usage


usage:
	@echo "------------------------------------------------------------------------"
	@echo "Makefile for final-codes, specify which test to run:"
	@echo ""
	@echo "-- QUADRATURE --"
	@echo "  adapgaus            adaptive gauss-legendre quadrature"
	@echo "  gaussq              netlib code for several Gaussian quadratures"
	@echo "  alpertsqrt          gauss-trapezoidal rules for sqrt singularities"
	@echo ""
	@echo "-- SPECIAL FUNCTIONS --"
	@echo "  hermexps            code for hermite functions, nodes, quad weights"
	@echo "  hermexps3d          TODO : code for 3d hermite functions and quads"
	@echo "  jacexps             TODO"
	@echo "  legearc             compute legendre nodes in arclength on a curve"
	@echo ""
	@echo "-- ALGORITHMS --"
	@echo "  fftw_wrap           TODO"
	@echo "  rsortanyn           sorting routine for real arrays"
	@echo ""
	@echo "-- UTILITIES --"
	@echo "  clean               remove all object and executable files"
	@echo "  quaplot             old-school gnuplot 2D plotting"
	@echo "  pplot               python output plotting routines"
	@echo "  prini_omp           TODO : OpenMP compatible prini routines"
	@echo ""
	@echo "------------------------------------------------------------------------"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@


#-- QUADRATURE --
adapgaus: src/adapgaus.o src/prini.o src/legeexps.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_adapgaus.f -o build/int2 $^
	(cd build; ./int2)

gaussq: src/gaussq.o src/prini.o src/legeexps.o src/adapgaus.o src/adapgaus_quad.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_gaussq.f90 -o build/int2 $^
	(cd build; ./int2)

alpertsqrt: src/alpertsqrt.o src/prini.o src/adapgaus_quad.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_alpertsqrt.f90 -o build/int2 $^
	(cd build; ./int2)


#-- SPECIAL FUNCTIONS --
hermexps: src/hermexps.o src/prini.o src/pplot.o src/gammanew_eval.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_hermexps.f90 -o build/int2 $^
	(cd build; ./int2)

hermexps3d: src/hermexps3d.o src/hermexps.o src/prini.o src/pplot.o src/gammanew_eval.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_hermexps3d.f90 -o build/int2 $^
	(cd build; ./int2)

legearc: src/legearc.o src/prini.o src/legeexps.o src/quaplot.o src/adapgaus.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_legearc.f -o build/int2 $^
	(cd build; ./int2)




#-- ALGORITHMS --
fftw_wrap: src/prini.o

rsortanyn: src/prini.o src/rsortanyn.o src/quaplot.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_rsortanyn.f -o build/int2 $^
	(cd build; ./int2)



#-- UTILITIES --
quaplot: src/quaplot.o src/prini.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_quaplot.f -o build/int2 $^
	(cd build; ./int2)

pplot: src/pplot.o src/prini.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_pplot.f90 -o build/int2 $^
	(cd build; ./int2)



clean:
	rm -f src/*.o
	rm -f testing/*.o
