PROG = $(HOME)/bin/ensnergy
OBJSNOVEC = \
ebrain.o\
ensnergy.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
FFLAGS = -O

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
