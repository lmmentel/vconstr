PROG = $(HOME)/bin/rcvxc
OBJSNOVEC = \
cbrain.o\
rcvxc.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
FFLAGS = -O2
LFLAGS =

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
