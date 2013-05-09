PROG = $(HOME)/bin/dsfun
OBJSNOVEC = \
clcens.o\
clwxci.o\
lrchxc.o\
vlchxc.o\
grdrho.o\
dvksmt.o\
dbrain.o\
dsfun.o\
brains.o\
sntial.o\
vntial.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/atmol/libatm.a
#LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/Vconstr/lib/xclib.a $(HOME)/MYRTA/ATMOL/atmol/lib/libatm.a
FFLAGS = -O2
LFLAGS =

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf -bmaxdata:419430400 $(FFLAGS) -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
