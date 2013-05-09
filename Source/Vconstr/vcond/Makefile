PROG = ${HOME}/bin/vcond
OBJSNOVEC = \
cond.o\
p1dens.o\
rdsamo.o\
vbrain.o\
vcndex.o\
unpack.o

LIBS = $(HOME)/atmol_dist/build/lib/libatm.a $(HOME)/Vconstr/lib/dnslib.a
FFLAGS = -O
LFLAGS = -O

$(PROG): $(OBJSNOVEC) $(LIBS)
	xlf $(LFLAGS) -o $(PROG) *.o $(LIBS)

$(OBJSNOVEC):;
	xlf $(FFLAGS) -c $*.f

