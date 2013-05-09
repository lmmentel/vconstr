      subroutine aovlsv(x0,y0,z0,valao,gradx,grady,gradz,grads,vn)
c**********************************************************************
c
c calculate values, gradients and Laplacian of all aO's (gaussians) 
c   in point (x0,y0,z0)
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
      parameter(explim = 500.d0)
c
c** /gaus3/ is read from dumpfile in routines getgss an gssinv
c
      common/gaus3/zeta(10,160),coef(10,160),cx(100),cy(100),cz(100),
     1 zz(100),itype(160),kad(160),icentr(160),mbase(160),ntrm(160),
     2 title(10),acc1,acc2,imin(4),imax(4),nbasis,natoms,iblock,maxblk,
     3 irun,ngroup,ncentg,nosym,cgx(12),cgy(12),cgz(12),idddd,norbd,
     4 nnbase(160)      
c
      dimension valao(nbasis),gradx(nbasis),grady(nbasis),
     +  gradz(nbasis),grads(nbasis)
c
      valao(1:nbasis)=0.d0
      gradx(1:nbasis)=0.d0
      grady(1:nbasis)=0.d0
      gradz(1:nbasis)=0.d0
      grads(1:nbasis)=0.d0
c
      sqrt3  = sqrt(3.d0)
      sqrt3i = 1.d0/sqrt3
c
c** calculate nuclear potential
c
      vn = 0.0d0
      do i = 1,natoms
        x2=(x0-cx(i))**2
        y2=(y0-cy(i))**2
        z2=(z0-cz(i))**2
        rr=sqrt(x2+y2+z2)
        if (rr.gt.1.0e-75) then
          vn = vn+zz(i)/rr
        else
          vn = vn+1.0d75
        endif
      enddo
c
      do 20 igroup = 1,ngroup
c
c** transform coordinate (x0,y0,z0) from main coordinate frame to
c** coordinate (x,y,z) in a frame with the atomic orbital center 
c** as an origin
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
          dprb  = 2*zeta(itrm,igroup)*prb
          ar    = 2*zeta(itrm,igroup)*r2
          ax    = -2*zeta(itrm,igroup)*x
          ay    = -2*zeta(itrm,igroup)*y
          az    = -2*zeta(itrm,igroup)*z
          goto (1,2,3) itype(igroup)
c
c** s type atomic orbitals
c
    1 valao(iaobas+1) = valao(iaobas+1)+prb
      gradx(iaobas+1) = gradx(iaobas+1)+ax*prb
      grady(iaobas+1) = grady(iaobas+1)+ay*prb
      gradz(iaobas+1) = gradz(iaobas+1)+az*prb
      grads(iaobas+1) = grads(iaobas+1)+dprb*(ar-3.d0)
      goto 10
c
c** p type atomic orbitals
c
    2 valao(iaobas+1) = valao(iaobas+1)+prb*x
      valao(iaobas+2) = valao(iaobas+2)+prb*y
      valao(iaobas+3) = valao(iaobas+3)+prb*z
      prbx = prb*x
      prby = prb*y
      prbz = prb*z
      gradx(iaobas+1) = gradx(iaobas+1)+ax*prbx+prb
      grady(iaobas+1) = grady(iaobas+1)+ay*prbx
      gradz(iaobas+1) = gradz(iaobas+1)+az*prbx
      gradx(iaobas+2) = gradx(iaobas+2)+ax*prby
      grady(iaobas+2) = grady(iaobas+2)+ay*prby+prb
      gradz(iaobas+2) = gradz(iaobas+2)+az*prby
      gradx(iaobas+3) = gradx(iaobas+3)+ax*prbz
      grady(iaobas+3) = grady(iaobas+3)+ay*prbz
      gradz(iaobas+3) = gradz(iaobas+3)+az*prbz+prb
      dprx = dprb*x
      dpry = dprb*y
      dprz = dprb*z
      grads(iaobas+1) = grads(iaobas+1)+dprx*(ar-5.d0)
      grads(iaobas+2) = grads(iaobas+2)+dpry*(ar-5.d0)
      grads(iaobas+3) = grads(iaobas+3)+dprz*(ar-5.d0)
      goto 10
c
c** d type atomic orbitals
c
    3 valao(iaobas+1) = valao(iaobas+1)+prb*x*y
      valao(iaobas+2) = valao(iaobas+2)+prb*x*z
      valao(iaobas+3) = valao(iaobas+3)+prb*y*z
c
      prbx  = prb*x
      prby  = prb*y
      prbz  = prb*z
      prbxx = prbx*x
      prbxy = prbx*y
      prbxz = prbx*z
      prbyy = prby*y
      prbyz = prby*z
      prbzz = prbz*z
      gradx(iaobas+1) = gradx(iaobas+1)+ax*prbxy+prby
      grady(iaobas+1) = grady(iaobas+1)+ay*prbxy+prbx
      gradz(iaobas+1) = gradz(iaobas+1)+az*prbxy
      gradx(iaobas+2) = gradx(iaobas+2)+ax*prbxz+prbz
      grady(iaobas+2) = grady(iaobas+2)+ay*prbxz
      gradz(iaobas+2) = gradz(iaobas+2)+az*prbxz+prbx
      gradx(iaobas+3) = gradx(iaobas+3)+ax*prbyz
      grady(iaobas+3) = grady(iaobas+3)+ay*prbyz+prbz
      gradz(iaobas+3) = gradz(iaobas+3)+az*prbyz+prby
c
      dprx  = dprb*x
      dpry  = dprb*y
      dprz  = dprb*z
      dprxx = dprx*x
      dprxy = dprx*y
      dprxz = dprx*z
      dpryy = dpry*y
      dpryz = dpry*z
      dprzz = dprz*z
      grads(iaobas+1) = grads(iaobas+1)+dprxy*(ar-7.d0)
      grads(iaobas+2) = grads(iaobas+2)+dprxz*(ar-7.d0)
      grads(iaobas+3) = grads(iaobas+3)+dpryz*(ar-7.d0)
c
c** calculate x2, y2, and z2 values
c
      prb = prb*sqrt3i
      valao(iaobas+4) = valao(iaobas+4)+prb*x*x
      valao(iaobas+5) = valao(iaobas+5)+prb*y*y
      valao(iaobas+6) = valao(iaobas+6)+prb*z*z
c
      prbx  = prbx*sqrt3i
      prby  = prby*sqrt3i
      prbz  = prbz*sqrt3i
      prbxx = prbxx*sqrt3i
      prbyy = prbyy*sqrt3i
      prbzz = prbzz*sqrt3i
      gradx(iaobas+4) = gradx(iaobas+4)+ax*prbxx+2.d0*prbx
      grady(iaobas+4) = grady(iaobas+4)+ay*prbxx
      gradz(iaobas+4) = gradz(iaobas+4)+az*prbxx
      gradx(iaobas+5) = gradx(iaobas+5)+ax* prbyy
      grady(iaobas+5) = grady(iaobas+5)+ay*prbyy+2.d0*prby
      gradz(iaobas+5) = gradz(iaobas+5)+az*prbyy
      gradx(iaobas+6) = gradx(iaobas+6)+ax*prbzz
      grady(iaobas+6) = grady(iaobas+6)+ay*prbzz
      gradz(iaobas+6) = gradz(iaobas+6)+az*prbzz+2.d0*prbz
c
      dprb  = dprb*sqrt3i
      dprx  = dprx*sqrt3i
      dpry  = dpry*sqrt3i
      dprz  = dprz*sqrt3i
      dprxx = dprxx*sqrt3i
      dpryy = dpryy*sqrt3i
      dprzz = dprzz*sqrt3i
      grads(iaobas+4) = grads(iaobas+4)+dprxx*(ar-7.d0)+2.d0*prb
      grads(iaobas+5) = grads(iaobas+5)+dpryy*(ar-7.d0)+2.d0*prb
      grads(iaobas+6) = grads(iaobas+6)+dprzz*(ar-7.d0)+2.d0*prb
c
          endif
   10   continue
   20 continue
c
      return
      end
c
      subroutine aovlsw(x0,y0,z0,valao,gradx,grady,gradz,grads,vn)
c**********************************************************************
c
c calculate the Laplacian of all ao's in point (x0,y0,z0), and return 
c   them in the array grads
c
c *******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(nmxl = 6)
      parameter(explim = 500.d0)
c
c** /gaus4/ is read from dumpfile in routines getgss an gssinw
c
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq
c
      dimension valao(nbasis),gradx(nbasis),grady(nbasis),
     +  gradz(nbasis),grads(nbasis)
      dimension cosphi(nmxl),sinphi(nmxl),
     +  pnl((nmxl+1)*(nmxl+2)/2),dpnl((nmxl+1)*(nmxl+2)/2),
     +  d2pnl((nmxl+1)*(nmxl+2)/2),tmat(9),dtmat(27),
     +  fodsph(3),sodsph(6),dao(3),d2ao(6),newm(9)
c
      pi = acos(-1.d0)
c
      valao(1:nbasis)=0.d0
      gradx(1:nbasis)=0.d0
      grady(1:nbasis)=0.d0
      gradz(1:nbasis)=0.d0
      grads(1:nbasis)=0.d0
c
c** calculate nuclear potential
c
      vn = 0.0d0
      do i = 1,ncen
        cx=cart(3*(i-1)+1)
        cy=cart(3*(i-1)+2)
        cz=cart(3*(i-1)+3)
        x2=(x0-cx)**2
        y2=(y0-cy)**2
        z2=(z0-cz)**2
        rr=sqrt(x2+y2+z2)
        if (rr.gt.1.0e-75) then
          vn = vn+zz(i)/rr
        else
          vn = vn+1.0d75
        endif
      enddo
c
c** loop over all cgto's (or groups)
c
      do 120 igroup = 1,ngrp
c
c** calculate coordinates of (x0,y0,z0) in local coordinate frame with
c** cgto-center as origin.
c
        x = x0 - cart(icen3(igroup)+1)
        y = y0 - cart(icen3(igroup)+2)
        z = z0 - cart(icen3(igroup)+3)
c
c** prepare
c
        lx=max(lquant(igroup),1)
        call prcrds(lx,x,y,z,r,costht,sintht,cosphi,sinphi)
        call capsph(lx,costht,sintht,pnl,dpnl,d2pnl)
        call trfdrv(r,costht,sintht,cosphi(1),sinphi(1),tmat,dtmat)
c
        r2=x**2+y**2+z**2
        iaobas=norb(igroup)
c
        ll=lquant(igroup)
        if (ll+1.gt.5) then
          print*,'ERROR, H-type functions.'
          stop
        endif
c
        mnmbr=2*ll+1
        call trnsfm(newm,mnmbr)
c
        do 110 itrm = 1,ntrm(igroup)
          izeta=ilifc(igroup)+itrm
          expon=r2*zeta(izeta)
          if (expon.lt.explim) then
            gssfn=ctran(izeta)*exp(-expon)
            alfrs=2*expon
c
c** First check if this is a special CGTO with NQUANT=1.  This is a
c** S orbital : S=D(x*x)+D(y*y)+D(z*z),  which can be added to
c** simulate a 6-d set instead of a 5-d set.
c
          if (nquant(igroup).eq.1) print*,'WARNING : special CGTO'
c
          if (ll.eq.0) then
            prb  =gssfn
            dprb =-2.d0*zeta(izeta)*r*gssfn
            d2prb=2.d0*zeta(izeta)*(alfrs-1.d0)*gssfn
          elseif (ll.eq.1) then
            prb  =r*gssfn
            dprb =(1.d0-alfrs)*gssfn
            d2prb=2.d0*zeta(izeta)*r*(alfrs-3.d0)*gssfn
          else
            prb  =r**ll*gssfn
            dprb =r**(ll-1)*(ll-alfrs)*gssfn
            d2prb=r**(ll-2)*(ll*(ll-1)-alfrs*(2*ll+1)+
     + alfrs**2)*gssfn
          endif
          do m=1,mnmbr
            mm=m-ll-1
            lll=((ll+1)*(ll+2))/2-ll+abs(mm)
            call fsdsph(prb,dprb,d2prb,pnl(lll),dpnl(lll),d2pnl(lll),
     + lx,mm,cosphi,sinphi,tmat,dtmat,fodsph,sodsph,ao,dao,d2ao)
            nm=newm(m)
            valao(iaobas+nm) = valao(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*ao
            gradx(iaobas+nm) = gradx(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*dao(1)
            grady(iaobas+nm) = grady(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*dao(2)
            gradz(iaobas+nm) = gradz(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*dao(3)
            grads(iaobas+nm) = grads(iaobas+nm) + 
     + sqrt(4*pi/mnmbr)*(d2ao(1)+d2ao(4)+d2ao(6))
          enddo
        endif
c
  110   continue
  120 continue
c
       
      return
      end
c
      subroutine trnsfm(newm,mnmbr)
c***********************************************************************
c
c  transform the value of m for the new subroutine aovaln, so that the 
c    order of the orbitals is identical to that in the old subroutine aovalw
c
c**********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension newm(mnmbr)
      newm(1:mnmbr)=0
c
      if (mnmbr.eq.1) then
c** S type atomic orbital
        newm(1)=1
      else if (mnmbr.eq.3) then
c** P type atomic orbital
        newm(1)=3
        newm(2)=1
        newm(3)=2
      else if (mnmbr.eq.5) then
c** D type atomic orbital
        newm(1)=5
        newm(2)=3
        newm(3)=1
        newm(4)=2
        newm(5)=4
      else if (mnmbr.eq.7) then
c** F type atomic orbital
        newm(1)=7
        newm(2)=5
        newm(3)=3
        newm(4)=1
        newm(5)=2
        newm(6)=4
        newm(7)=6
      else if (mnmbr.eq.9) then
c** G type atomic orbital
        newm(1)=9
        newm(2)=7
        newm(3)=5
        newm(4)=3
        newm(5)=1
        newm(6)=2
        newm(7)=4
        newm(8)=6
        newm(9)=8
      else
        write(6,'(/,'' ERROR in subroutine trfsfm. L='',i2)')mnmbr
        stop
      endif
      return
      end
