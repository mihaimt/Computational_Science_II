CC = gcc
FC = gfortran
OPT = -O1
DEBUG = -g -Wall
FFLAGS = -funderscoring -fno-range-check -fmessage-length=0
CFLAGS = -fno-leading-underscore 
DEPS =
IDIR =
LIBS = 
FOBJ_LEVEL_0 = grav_force.o
FOBJ_LEVEL_1 = grav_force_lvl1.o
FOBJ_LEVEL_2 = grav_force_lvl2.o
COBJ = c_invsqrt64.o
EXE_LEVEL_0 = grav_force
EXE_LEVEL_1 = grav_force_lvl1
EXE_LEVEL_2 = grav_force_lvl2


%.o: %.c $(DEPS)
	$(CC) -c $(DEBUG) $(OPT) -o $@ $< $(CFLAGS)

%.o: %.f90 $(DEPS)
	$(FC) -c $(DEBUG) $(OPT) -o $@ $< $(FFLAGS)

$(EXE_LEVEL_0): $(FOBJ_LEVEL_0) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_LEVEL_1): $(FOBJ_LEVEL_1) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

$(EXE_LEVEL_2): $(FOBJ_LEVEL_2) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

all: $(EXE_LEVEL_0) $(EXE_LEVEL_1) $(EXE_LEVEL_2)

.PHONY: clean

clean:	cleanobj
	rm -f $(EXE_LEVEL_0) $(EXE_LEVEL_1) $(EXE_LEVEL_2)

cleanobj:
	rm -f *.o *~ core
