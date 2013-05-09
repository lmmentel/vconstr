PROG = $(HOME)/bin/scf
OBJSNOVEC = \
sbrain.o\
clcexc.o\
clcmod.o\
clcmof.o\
clrhod.o\
clrhof.o\
clcvxc.o\
cledet.o\
clsorb.o\
clvmat.o\
clcwxi.o\
prepvs.o\
clvkli.o\
clcdvx.o\
scf.o\
spavxc.o\
vntial.o\
wrt.o

LIBS = $(HOME)/Vconstr/lib/dnslib.a $(HOME)/atmol/libatm.a $(HOME)/Vconstr/lib/xclib.a
FFLAGS = -O

$(PROG) :  $(OBJSNOVEC) $(LIBS)
					xlf $(FFLAGS) -o $(PROG) *.o $(LIBS) 

$(OBJSNOVEC):;
					xlf $(FFLAGS) -c $*.f
