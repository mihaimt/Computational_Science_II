CC = gcc
FC = gfortran
OPT = -O1
DEBUG = -g -Wall
FFLAGS = -funderscoring -fno-range-check -fmessage-length=0
CFLAGS = -fno-leading-underscore 
DEPS =
IDIR =
LIBS = 
FOBJ = grav_force.o 
COBJ = c_invsqrt64.o
EXE = grav_force


%.o: %.c $(DEPS)
	$(CC) -c $(DEBUG) $(OPT) -o $@ $< $(CFLAGS)

%.o: %.f90 $(DEPS)
	$(FC) -c $(DEBUG) $(OPT) -o $@ $< $(FFLAGS)

$(EXE): $(FOBJ) $(COBJ)
	$(FC) -o $@ $^ $(LIBS)

.PHONY: clean

clean:	cleanobj
	rm -f $(EXE)

cleanobj:
	rm -f *.o *~ core
