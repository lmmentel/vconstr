PROG = $(HOME)/bin/excdns
OBJSNOVEC = \
clcexc.o\
clcrho.o\
clcvxc.o\
cledet.o\
ebrain.o\
exc.o\
wrt.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
FFLAGS = -O

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf -bmaxdata:419430400 -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
