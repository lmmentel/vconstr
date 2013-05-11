.f.a:
#
# Makefile for dnstyfun-library
#
FC=ifort

PROG=../xclib
OBJSNOVEC = \
$(PROG).a(corlsd.o)\
$(PROG).a(gcor.o)\
$(PROG).a(lyp.o)\
$(PROG).a(pw91c.o)\
$(PROG).a(pw91c1.o)\
$(PROG).a(pw91x.o)\
$(PROG).a(pw91x1.o)\
$(PROG).a(xcbeck.o)\
$(PROG).a(xcener.o)\
$(PROG).a(xcpot.o)\
$(PROG).a(xcprep.o)\
$(PROG).a(xcvpot.o)\
$(PROG).a(xcvwnd.o)\
$(PROG).a(xcvwne.o)\
$(PROG).a(xcvwnf.o)\
$(PROG).a(xcvwnp.o)\
$(PROG).a(xcprdw.o)

FFLAGS = -O2 -autodouble -fixed 

$(PROG).a:: $(OBJSNOVEC) $(OBJSVEC);
				ar -ru $(PROG).a *.o

$(OBJSNOVEC):;
				$(FC) $(FFLAGS) -c $*.f

