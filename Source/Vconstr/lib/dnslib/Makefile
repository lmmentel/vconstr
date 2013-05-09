.f.a:
#
# Makefile for dnstyfun-library
#

PROG=$(HOME)/Vconstr/lib/dnslib
OBJSNOVEC = \
$(PROG).a(aovl.o)\
$(PROG).a(aovls.o)\
$(PROG).a(aovlsn.o)\
$(PROG).a(bd.o)\
$(PROG).a(calch.o)\
$(PROG).a(calcni.o)\
$(PROG).a(ccint.o)\
$(PROG).a(clcdip.o)\
$(PROG).a(clcpci.o)\
$(PROG).a(clcpot.o)\
$(PROG).a(clcvhr.o)\
$(PROG).a(clcvls.o)\
$(PROG).a(clnhmx.o)\
$(PROG).a(clvkin.o)\
$(PROG).a(cpclr.o)\
$(PROG).a(dmp.o)\
$(PROG).a(epsbvw.o)\
$(PROG).a(erep.o)\
$(PROG).a(exchan.o)\
$(PROG).a(fcalch.o)\
$(PROG).a(fccint.o)\
$(PROG).a(ftrans.o)\
$(PROG).a(frcctr.o)\
$(PROG).a(getgss.o)\
$(PROG).a(getdmp.o)\
$(PROG).a(hmatmo.o)\
$(PROG).a(init.o)\
$(PROG).a(matrix.o)\
$(PROG).a(rddmp.o)\
$(PROG).a(rdhfmo.o)\
$(PROG).a(rdsmo.o)\
$(PROG).a(rdmat.o)\
$(PROG).a(sym.o)\
$(PROG).a(valgrd.o)\
$(PROG).a(vksmat.o)\
$(PROG).a(wmtrx.o)\
$(PROG).a(wrtkso.o)\
$(PROG).a(wrtorb.o)\
$(PROG).a(ylms.o)

FFLAGS = -O2

$(PROG).a:: $(OBJSNOVEC) $(OBJSVEC);
				ar -ru $(PROG).a *.o

$(OBJSNOVEC):;
				xlf $(FFLAGS) -c $*.f
