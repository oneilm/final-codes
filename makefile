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
	@echo "-- SPECIAL FUNCTIONS --"
	@echo "  hermexps            code for hermite functions, nodes, quad weights"
	@echo ""
	@echo "-- UTILITIES --"
	@echo "  clean               remove all object and executable files"
	@echo "  pplot               python output plotting routines"
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



#-- SPECIAL FUNCTIONS --
hermexps: src/hermexps.o src/prini.o src/pplot.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_hermexps.f90 -o build/int2 $^
	(cd build; ./int2)


#-- UTILITIES --
pplot: src/pplot.o src/prini.o
	rm -f build/*
	$(FC) $(FFLAGS) testing/test_pplot.f90 -o build/int2 $^
	(cd build; ./int2)



clean:
	rm -f src/*.o
	rm -f testing/*.o
