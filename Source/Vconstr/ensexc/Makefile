PROG = $(HOME)/bin/ensexc
OBJSNOVEC = \
clcrho.o\
clcvxc.o\
ebrain.o\
ensexc.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
FFLAGS = -O

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf $(FFLAGS) -bmaxdata:419430400 -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
