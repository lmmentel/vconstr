      subroutine getgss(iunit,nbas,atmol4)
c******************************************************************
c
c This routine determines if the gaussians are stored on an ATMOL3
c dumpfile or an ATMOL4 dumpfile and it calls the appropriate routine 
c to restored the gaussians from dumpfile.
c
c******************************************************************
      implicit real*8  (a-h,o-z),integer  (i-n)
c
      logical atmol4
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq  

c
      write(6,'(''I approach secini'')')
      call secini(1,iunit)
      write(6,'(''I passed secini'')')
      call secget(191,1,ibloc)
      write(6,'(''I passed secget'')')
      ibl = ibloc+9+1-nipw()
      call search(ibl,iunit)
      write(6,'(''I passed iunit'')')
      call find(iunit)
      write(6,'(''I passed find'')')
      call get(zeta,nw)
      write(6,'(''I passed get'')')
      atmol4=.false.
      if (nw.eq.511) atmol4=.true.
c
      if (atmol4) then 
        call gssinw(iunit,nbas)
        write(6,'(''I passed gssinw'')')
      else
        call gssinv(iunit,nbas)
        write(6,'(''I passed gssinw'')')
      endif
c
      return
      end
c
      subroutine gssinv(iunit,nbas)
c *****************************************************************
c
c  get contracted gaussians from 'atmol3' dumpfile
c
c  !! these functions are not normalized. an additional factor of
c     2**(-1/4) * pi**(-5/8) must be used. (see (old) manual, page 4.53)
c
c *******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/gaus3/zeta(10,160),coef(10,160),cx(100),cy(100),cz(100),
     1 zz(100),itype(160),kad(160),icentr(160),mbase(160),ntrm(160),
     2 title(10),acc1,acc2,imin(4),imax(4),nbasis,natoms,iblock,maxblk,
     3 irun,ngroup,ncentg,nosym,cgx(12),cgy(12),cgz(12),idddd,norbd,
     4 nnbase(160)
c
c** factor = 2**(-1/4) * pi**(-5/8) ***
c
      factor = 0.41117228237d0
      n4626 = 2*10*160+400+5*160/nipw()+12+16/nipw()+
     + 3*12+2/nipw()+160/nipw()
c
      call secget(191,1,ibloc)
      call search(ibloc,iunit)
      call reads(zeta,n4626,iunit)
      if (idddd.eq.0) then
        nbas = nbasis
      else
        nbas = norbd
      end if
      do i = 1,ngroup
        do j = 1,ntrm(i)
          coef(j,i) = coef(j,i)*factor
        enddo
      enddo
      write(6,'(/,''  restored from integv dumpfile :'',i4,
     +'' contracted gaussians'')') nbas
c
      return
      end
c
      subroutine gssinw(iunit,nbas)
c
c** restore contracted gaussians from atmol4 dumpfile.
c
c ordering and normalizing coefficients of ao's :
c
c p : p(z)    n = 1.0
c     p(x)    n = 1.0
c     p(y)    n = 1.0
c
c d : d(3z*z-r*r)     n = 1/2
c     d(xz)           n = sqrt(3)
c     d(yz)           n = sqrt(3)
c     d(x*x-y*y)      n = sqrt(3/4)
c     d(xy)           n = sqrt(3)
c
c f : z*(5z*z-3*r*r)                   n = 1/2
c     3*x(5*z*z-r*r)                   n = sqrt(1/24)
c     3*y(5*z*z-r*r)                   n = sqrt(1/24)
c     15*z(x*x-y*y)                    n = sqr(1/60)
c     30xyz                            n = sqrt(1/60)
c     15(x*x*x-3*y*y*x)                n = sqrt(1/360)
c     15(3*y*x*x-y*y*y)                n = sqrt(1/360)
c
c g : 35z*z*z*z - 30z*z*r*r + 3r*r*r*r     n = 1/8
c     5(7z*z*z - 3z*r*r)*x                 n = sqrt(7/100)
c     5(7z*z*z - 3z*r*r)*y                 n = sqrt(7/100)
c     15(7z*z-r*r)(x*x-y*y)                n = sqrt(1/720)
c     15(7z*z-r*r)(x*y)                    n = sqrt(1/180)
c     105(x*x-3*y*y)zx                     n = sqrt(1/2520)
c     105(3*x*x-y*y)zy                     n = sqrt(1/2520)
c     105(x*x*x*x-6*x*x*y*y+y*y*y*y)       n = sqrt(1/20160)
c     105(4*x*x*x*y-4*x*y*y*y)             n = sqrt(1/20160)
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq
c
c*** fac = 2**(-1/4) * pi**(-5/8) ***
c
      fac = 0.41117228237d0
      n16 = 16/nipw()
      n7018 = 2*1300+10*340/nipw()+2/nipw()+10+600
     *      +2*200+2/nipw()+5
c
      call secget(191,1,ibloc)
      call search(ibloc,iunit)
      call reads(zeta,n16,iunit)
      call reads(ctran,n7018,iunit)
      do 60 igrp = 1,ngrp
        do 50 j = 1,ntrm(igrp)
          ctran(ilifc(igrp)+j) = ctran(ilifc(igrp)+j)*fac*
     + sqrt(2*lquant(igrp)+1.0d0)
   50   continue
   60 continue
      nbas = nbasis
      write(6,'(/,''  restored from atmol4 dumpfile :'',i4,
     +'' contracted gaussians'')') nbas
c
      return
      end
