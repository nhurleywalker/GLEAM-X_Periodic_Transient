CC = gcc

# If necessary, edit the following to ensure that the CSPICE library and header file (SpiceUsr.h) are locatable
LDLIBS = -l:cspice.a -lm
CFLAGS = -Wall -Wextra -march=native -Ofast

TARGETS = bc_corr \
		  de430.bsp

all: $(TARGETS)

bc_corr: bc_corr.o barycentre.o vec.o
#predictions: predictions.o barycentre.o vec.o

de430.bsp:
	wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/$@

clean:
	$(RM) $(TARGETS) *.o
