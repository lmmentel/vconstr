PROG = $(HOME)/bin/spinexc
OBJSNOVEC = \
clcexc.o\
clcrho.o\
clcvxc.o\
ebrain.o\
exc.o\
wrt.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
FFLAGS = -O
LFLAGS =

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
