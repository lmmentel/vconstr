      subroutine aovlv(x0,y0,z0,valao)
c**********************************************************************
c
c calculate values of all ao's (gaussians) in point (x0,y0,z0)
c
c !! for the functions d(x*x), d(y*y) and d(z*z) an extra 
c   normalization factor of 1/sqrt(3) must be added
c
c !! extra normalization factors for d(x*x-y*y) = N*(d(x*x)-d(y*y))
c    and d(3z*z-R*R) = M(2d(z*z) - d(x*x) - d(y*y))   are :
c       N = sqrt(3/4)      M = 1/2
c
c**********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(explim = 300.0d0)
c
c** /gaus3/ is read from dumpfile in routines getgss 
c
      common/gaus3/zeta(10,160),coef(10,160),cx(100),cy(100),cz(100),
     1 zz(100),itype(160),kad(160),icentr(160),mbase(160),ntrm(160),
     2 title(10),acc1,acc2,imin(4),imax(4),nbasis,natoms,iblock,maxblk,
     3 irun,ngroup,ncentg,nosym,cgx(12),cgy(12),cgz(12),idddd,norbd,
     4 nnbase(160)      
c
      dimension valao(nbasis)
      valao(1:nbasis)=0.d0
c
      do 20 igroup = 1,ngroup
c
c transform coordinate (x0,y0,z0) from main coordinate frame to
c   coordinate (x,y,z) in a frame with the atomic orbital center 
c   as an origin
c
        ic = icentr(igroup)
        x  = x0-cx(ic)
        y  = y0-cy(ic)
        z  = z0-cz(ic)
        r2 = x**2+y**2+z**2
c
        if (idddd.eq.1) print*,'WARNING : combine option not valid'
c
        iaobas = mbase(igroup)
c
c** do for every gaussian in the group 'igroup'
c
        do 10 itrm = 1,ntrm(igroup)
          expon = r2*zeta(itrm,igroup)
          if (expon.lt.explim) then
            prb   = coef(itrm,igroup)*exp(-expon)
            goto (1,2,3) itype(igroup)
c
c** s type atomic orbitals
c
    1 valao(iaobas+1) = valao(iaobas+1)+prb
      goto 20
c
c** p type atomic orbitals
c
    2 valao(iaobas+1) = valao(iaobas+1)+prb*x
      valao(iaobas+2) = valao(iaobas+2)+prb*y
      valao(iaobas+3) = valao(iaobas+3)+prb*z
      goto 20
c
c** d type atomic orbitals
c
    3 valao(iaobas+1) = valao(iaobas+1)+prb*x*y
      valao(iaobas+2) = valao(iaobas+2)+prb*x*z
      valao(iaobas+3) = valao(iaobas+3)+prb*y*z
c
c*** calculate x2, y2, and z2 values
c
      prb = prb/sqrt(3.0d0)
      valao(iaobas+4) = valao(iaobas+4)+prb*x*x
      valao(iaobas+5) = valao(iaobas+5)+prb*y*y
      valao(iaobas+6) = valao(iaobas+6)+prb*z*z
c
      endif
c
   10   continue
   20 continue
c
      return
      end
c
      subroutine aovlw(x0,y0,z0,valao)
c*******************************************************************
c ao-function values. ao's restored from atmol4 dumpfile.
c for ordering and normalizing coeficients, see routine gssinw.
c*******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(explim = 300d0)
c
c** /gaus4/ is read from dumpfile in routines getgss an gssinw
c
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq
c
      dimension valao(nbasis)
      valao(1:nbasis)=0.d0
c
c** calculate normalizing coefficients
c
      sqrt3 = sqrt(3.d0)
      sqrt34 = 0.5d0*sqrt3
c
      xd0 = 0.5d0
      xd1 = sqrt3
      xdmin1 = sqrt3
      xd2 = sqrt34
      xdmin2 = sqrt3
      xf0 = 0.5d0
      xf1 = sqrt(1.d0/24.0d0)
      xfmin1 = xf1
      xf2 = sqrt(1.d0/60.0d0)
      xfmin2 = xf2
      xf3 = sqrt(1.d0/360.0d0)
      xfmin3 = xf3
      xg0 = 1.d0/8.0d0
      xg1 = sqrt(7.0d0/100.0d0)
      xgmin1 = xg1
      xg2 = sqrt(1.d0/720.0d0)
      xgmin2 = 2.0d0*xg2
      xg3 = sqrt(1.d0/2520.0d0)
      xgmin3 = xg3
      xg4 = sqrt(1.d0/20160.0d0)
      xgmin4 = xg4
c
c** loop over all cgto's (or groups)
c
      do 120 igroup = 1,ngrp
c
c** calculate coordinates of (x0,y0,z0) in local coordinate frame with
c** cgto-center as origin.
c
        x = x0-cart(icen3(igroup)+1)
        y = y0-cart(icen3(igroup)+2)
        z = z0-cart(icen3(igroup)+3)
        x2 = x*x
        y2 = y*y
        z2 = z*z
        r2 = x2+y2+z2
        iaobas = norb(igroup)
        do 110 itrm = 1,ntrm(igroup)
          izeta = ilifc(igroup)+itrm
          expon = r2*zeta(izeta)
          if (expon.lt.explim) then
            prb = ctran(izeta)*exp(-expon)
c
c** first check if this is a special cgto with nquant=1.  this is a
c** s orbital : s=d(x*x)+d(y*y)+d(z*z),  which can be added to
c** simulate a 6-d set instead of a 5-d set.
c
      if (nquant(igroup).eq.1)  then
        valao(iaobas+1) = valao(iaobas+1)+prb*r2/sqrt(5.d0)
      goto 110
      endif
c
      goto (1,2,3,4,5) lquant(igroup)+1
c
c** s type atomic orbitals
    1 valao(iaobas+1) = valao(iaobas+1)+prb
      goto 110
c
c** p type atomic orbitals
    2 valao(iaobas+1) = valao(iaobas+1)+prb*z
      valao(iaobas+2) = valao(iaobas+2)+prb*x
      valao(iaobas+3) = valao(iaobas+3)+prb*y
      goto 110
c
c** d type atomic orbitals
    3 valao(iaobas+1) = valao(iaobas+1)+prb*(3*z2-r2)*xd0
      valao(iaobas+2) = valao(iaobas+2)+prb*(x*z)*xd1
      valao(iaobas+3) = valao(iaobas+3)+prb*(y*z)*xdmin1
      valao(iaobas+4) = valao(iaobas+4)+prb*(x2-y2)*xd2
      valao(iaobas+5) = valao(iaobas+5)+prb*(x*y)*xdmin2
      goto 110
c
c** f type atomic orbitals
    4 valao(iaobas+1) = valao(iaobas+1)+prb*z*(5*z2-3*r2)*xf0
      valao(iaobas+2) = valao(iaobas+2)+prb*3*x*(5*z2-r2)*xf1
      valao(iaobas+3) = valao(iaobas+3)+prb*3*y*(5*z2-r2)*xfmin1
      valao(iaobas+4) = valao(iaobas+4)+prb*15*z*(x2-y2)*xf2
      valao(iaobas+5) = valao(iaobas+5)+prb*30*x*y*z*xfmin2
      valao(iaobas+6) = valao(iaobas+6)+prb*15*x*(x2-3*y2)*xf3
      valao(iaobas+7) = valao(iaobas+7)+prb*15*y*(3*x2-y2)*xfmin3
      goto 110
c
c** g type atomic orbitals
    5 valao(iaobas+1) = valao(iaobas+1)+
     +                 prb*(35*z2*z2-30*z2*r2+3*r2*r2)*xg0
      valao(iaobas+2) = valao(iaobas+2)+
     +                 prb*5*x*z*(7*z2-3*r2)*xg1
      valao(iaobas+3) = valao(iaobas+3)+
     +                 prb*5*y*z*(7*z2-3*r2)*xgmin1
      valao(iaobas+4) = valao(iaobas+4)+
     +                 prb*15*(7*z2-r2)*(x2-y2)*xg2
      valao(iaobas+5) = valao(iaobas+5)+
     +                 prb*15*(7*z2-r2)*(x*y)*xgmin2
      valao(iaobas+6) = valao(iaobas+6)+
     +                 prb*105*(x2-3*y2)*z*x*xg3
      valao(iaobas+7) = valao(iaobas+7)+
     +                 prb*105*(3*x2-y2)*z*y*xgmin3
      valao(iaobas+8) = valao(iaobas+8)+
     +                 prb*105*(x2*x2-6*x2*y2+y2*y2)*xg4
      valao(iaobas+9) = valao(iaobas+9)+
     + prb*4*105*x*y*(x2-y2)*xgmin4
c
      endif
c
  110   continue
  120 continue
c
      return
      end
