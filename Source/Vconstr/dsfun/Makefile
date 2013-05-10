
FC=ifort

PROG = ../bin/dsfun
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

LIBS =  ../lib/xclib.a 
FFLAGS = -O2 -autodouble -fixed
LFLAGS =

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					$(FC)  $(FFLAGS) -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					$(FC) $(FFLAGS) -c $*.f
