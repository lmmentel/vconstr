.f.a:
#
# Makefile for dnstyfun-library
#

PROG=$(HOME)/Vconstr/lib/xclib
OBJSNOVEC = \
$(PROG).a(corlsd.o)\
$(PROG).a(gcor.o)\
$(PROG).a(hcth.o)\
$(PROG).a(lyp.o)\
$(PROG).a(pw91c.o)\
$(PROG).a(pw91c1.o)\
$(PROG).a(pw91x.o)\
$(PROG).a(pw91x1.o)\
$(PROG).a(vvv.o)\
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

FFLAGS = -O2

$(PROG).a:: $(OBJSNOVEC) $(OBJSVEC);
				ar -ru $(PROG).a *.o

$(OBJSNOVEC):;
				xlf $(FFLAGS) -c $*.f

