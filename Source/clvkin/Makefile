PROG = $(HOME)/bin/cvkin
OBJSNOVEC = \
cvkin.o\
vbrain.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
FFLAGS = -O2

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf -bmaxdata:419430400 $(FFLAGS) -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
