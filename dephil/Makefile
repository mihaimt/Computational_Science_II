# Makefile
##########  in 'make all', grav_force, grav_force_lvl1 and grav_force_lvl2 are exculded, since grav_force_lvlx takes over their function

# Definitions

CC = gcc
FC = gfortran
OPT = -O2
DEBUG = -g -Wall
FFLAGS = -funderscoring -fno-range-check -fmessage-length=0
CFLAGS = -fno-leading-underscore 
DEPS = module_grav_IO.o module_grav_precalcs.o module_grav_mass.o module_grav_parameters.o
IDIR =
LIBS = 
FOBJ_LEVEL_0 = archive/grav_force.o
FOBJ_LEVEL_1 = archive/grav_force_lvl1.o
FOBJ_LEVEL_2 = archive/grav_force_lvl2.o
FOBJ_LEVEL_X = archive/grav_force_lvlx.o
FOBJ_OSC_F = archive/oscillation_force.o
FOBJ_PURE = grav_pure_lvlx.o
COBJ = bin/c_invsqrt64.o
EXE_LEVEL_0 = grav_force
EXE_LEVEL_1 = grav_force_lvl1
EXE_LEVEL_2 = grav_force_lvl2
EXE_LEVEL_X = grav_force_lvlx
EXE_OSC_F = oscillation_force
EXE_PURE = grav_pure_lvlx


%.o: %.c
	$(CC) -c $(DEBUG) $(OPT) -o $@ $< $(CFLAGS)

%.o: %.f90
	$(FC) -c $(DEBUG) $(OPT) -o $@ $< $(FFLAGS)

default: $(EXE_OSC_M) $(EXE_LEVEL_X)

all: $(EXE_LEVEL_X) $(EXE_OSC_F) $(EXE_OSC_M)

$(EXE_LEVEL_0): $(FOBJ_LEVEL_0) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_LEVEL_1): $(FOBJ_LEVEL_1) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_LEVEL_2): $(FOBJ_LEVEL_2) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_LEVEL_X): $(FOBJ_LEVEL_X) $(DEPS) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_OSC_F): $(FOBJ_OSC_F) $(DEPS)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_PURE): $(FOBJ_PURE) $(DEPS)
	$(FC) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *.mod *~ core
	cd archive/; rm -f *.o *.mod *~ core
	cd bin/; rm -f *.o *.mod *~ core

cleanall: clean
	rm -f $(EXE_LEVEL_0) $(EXE_LEVEL_1) $(EXE_LEVEL_2) $(EXE_LEVEL_X) $(EXE_OSC_F) $(EXE_PURE)


# Dependencies

$(FOBJ_LEVEL_X) : archive/$(EXE_LEVEL_X).f90 module_grav_IO.o module_grav_precalcs.o module_grav_mass.o module_grav_parameters.o

$(FOBJ_OSC_F) : archive/$(EXE_OSC_F).f90 module_grav_IO.o module_grav_precalcs.o module_grav_mass.o module_grav_parameters.o

$(FOBJ_PURE) : $(EXE_PURE).f90 module_grav_IO.o module_grav_precalcs.o module_grav_mass.o module_grav_parameters.o
