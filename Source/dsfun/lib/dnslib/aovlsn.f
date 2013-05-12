      subroutine aovlnv(x0,y0,z0,valao,gradx,grady,gradz,gradxx,
     + gradxy,gradxz,gradyy,gradyz,gradzz)
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
      parameter(explim = 250.d0)
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
     + gradz(nbasis),gradxx(nbasis),gradxy(nbasis),gradxz(nbasis),
     + gradyy(nbasis),gradyz(nbasis),gradzz(nbasis)
c
      valao(1:nbasis)=0.d0
      gradx(1:nbasis)=0.d0
      grady(1:nbasis)=0.d0
      gradz(1:nbasis)=0.d0
      gradxx(1:nbasis)=0.d0
      gradxy(1:nbasis)=0.d0
      gradxz(1:nbasis)=0.d0
      gradyy(1:nbasis)=0.d0
      gradyz(1:nbasis)=0.d0
      gradzz(1:nbasis)=0.d0
c
      sqrt3  = sqrt(3.0d0)
      sqrt3i = 1.0d0/sqrt3
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
        x2 = x**2
        y2 = y**2
        z2 = z**2
        r2 = x2+y2+z2
     
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
            dprb  = 2.d0*zeta(itrm,igroup)*prb
            axx   = 2.d0*zeta(itrm,igroup)*x2
            ayy   = 2.d0*zeta(itrm,igroup)*y2
            azz   = 2.d0*zeta(itrm,igroup)*z2
            ax    = 2.d0*zeta(itrm,igroup)*x
            ay    = 2.d0*zeta(itrm,igroup)*y
            az    = 2.d0*zeta(itrm,igroup)*z
            goto (1,2,3) itype(igroup)
c
c** s type atomic orbitals
c
    1 valao(iaobas+1)  = valao(iaobas+1)+prb
      gradx(iaobas+1)  = gradx(iaobas+1)-ax*prb
      grady(iaobas+1)  = grady(iaobas+1)-ay*prb
      gradz(iaobas+1)  = gradz(iaobas+1)-az*prb
      gradxx(iaobas+1) = gradxx(iaobas+1)+dprb*(axx-1.d0)
      gradxy(iaobas+1) = gradxy(iaobas+1)+ax*ay*prb
      gradxz(iaobas+1) = gradxz(iaobas+1)+ax*az*prb
      gradyy(iaobas+1) = gradyy(iaobas+1)+dprb*(ayy-1.d0)
      gradyz(iaobas+1) = gradyz(iaobas+1)+ay*az*prb
      gradzz(iaobas+1) = gradzz(iaobas+1)+dprb*(azz-1.d0)
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
      gradx(iaobas+1) = gradx(iaobas+1)-ax*prbx+prb
      grady(iaobas+1) = grady(iaobas+1)-ay*prbx
      gradz(iaobas+1) = gradz(iaobas+1)-az*prbx
      gradx(iaobas+2) = gradx(iaobas+2)-ax*prby
      grady(iaobas+2) = grady(iaobas+2)-ay*prby+prb
      gradz(iaobas+2) = gradz(iaobas+2)-az*prby
      gradx(iaobas+3) = gradx(iaobas+3)-ax*prbz
      grady(iaobas+3) = grady(iaobas+3)-ay*prbz
      gradz(iaobas+3) = gradz(iaobas+3)-az*prbz+prb
      dprx = dprb*x
      dpry = dprb*y
      dprz = dprb*z
      gradxx(iaobas+1) = gradxx(iaobas+1)+dprx*(axx-3.d0)
      gradxy(iaobas+1) = gradxy(iaobas+1)+dpry*(axx-1.d0)
      gradxz(iaobas+1) = gradxz(iaobas+1)+dprz*(axx-1.d0)
      gradyy(iaobas+1) = gradyy(iaobas+1)+dprx*(ayy-1.d0)
      gradyz(iaobas+1) = gradyz(iaobas+1)+prbx*ay*az
      gradzz(iaobas+1) = gradzz(iaobas+1)+dprx*(azz-1.d0)
      gradxx(iaobas+2) = gradxx(iaobas+1)+dpry*(axx-1.d0)
      gradxy(iaobas+2) = gradxy(iaobas+1)+dprx*(ayy-1.d0)
      gradxz(iaobas+2) = gradxz(iaobas+1)+ax*prby*az
      gradyy(iaobas+2) = gradyy(iaobas+1)+dpry*(ayy-3.d0)
      gradyz(iaobas+2) = gradyz(iaobas+1)+dprz*(ayy-1.d0)
      gradzz(iaobas+2) = gradzz(iaobas+1)+dpry*(azz-1.d0)
      gradxx(iaobas+3) = gradxx(iaobas+1)+dprz*(axx-1.d0)
      gradxy(iaobas+3) = gradxy(iaobas+1)+ax*ay*prbz
      gradxz(iaobas+3) = gradxz(iaobas+1)+dprx*(azz-1.d0)
      gradyy(iaobas+3) = gradyy(iaobas+1)+dprz*(ayy-1.d0)
      gradyz(iaobas+3) = gradyz(iaobas+1)+dpry*(azz-1.d0)
      gradzz(iaobas+3) = gradzz(iaobas+1)+dprz*(azz-3.d0)
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
      gradx(iaobas+1) = gradx(iaobas+1)-ax*prbxy+prby
      grady(iaobas+1) = grady(iaobas+1)-ay*prbxy+prbx
      gradz(iaobas+1) = gradz(iaobas+1)-az*prbxy
      gradx(iaobas+2) = gradx(iaobas+2)-ax*prbxz+prbz
      grady(iaobas+2) = grady(iaobas+2)-ay*prbxz
      gradz(iaobas+2) = gradz(iaobas+2)-az*prbxz+prbx
      gradx(iaobas+3) = gradx(iaobas+3)-ax*prbyz
      grady(iaobas+3) = grady(iaobas+3)-ay*prbyz+prbz
      gradz(iaobas+3) = gradz(iaobas+3)-az*prbyz+prby
c
      dprx  = dprb*x
      dpry  = dprb*y
      dprz  = dprb*z
      dprxy = dprx*y
      dprxz = dprx*z
      dpryz = dpry*z
      gradxx(iaobas+1) = gradxx(iaobas+1)+dprxy*(axx-3.d0)
      gradxy(iaobas+1) = gradxy(iaobas+1)+prb*(1.d0+axx*ayy-axx-ayy)
      gradxz(iaobas+1) = gradxz(iaobas+1)+dpryz*(axx-1.d0)
      gradyy(iaobas+1) = gradyy(iaobas+1)+dprxy*(ayy-3.d0)
      gradyz(iaobas+1) = gradyz(iaobas+1)+dprxz*(ayy-1.d0)
      gradzz(iaobas+1) = gradzz(iaobas+1)+dprxy*(azz-1.d0)
      gradxx(iaobas+2) = gradxx(iaobas+2)+dprxz*(axx-3.d0)
      gradxy(iaobas+2) = gradxy(iaobas+2)+dpryz*(axx-1.d0)
      gradxz(iaobas+2) = gradxz(iaobas+2)+prb*(1.d0+axx*azz-axx-azz)
      gradyy(iaobas+2) = gradyy(iaobas+2)+dprxz*(ayy-1.d0)
      gradyz(iaobas+2) = gradyz(iaobas+2)+dprxy*(azz-1.d0)
      gradzz(iaobas+2) = gradzz(iaobas+2)+dprxz*(azz-3.d0)
      gradxx(iaobas+3) = gradxx(iaobas+3)+dpryz*(axx-1.d0)
      gradxy(iaobas+3) = gradxy(iaobas+3)+dprxz*(ayy-1.d0)
      gradxz(iaobas+3) = gradxz(iaobas+3)+dprxy*(azz-1.d0)
      gradyy(iaobas+3) = gradyy(iaobas+3)+dpryz*(ayy-3.d0)
      gradyz(iaobas+3) = gradyz(iaobas+3)+prb*(1.d0+ayy*azz-ayy-azz)
      gradzz(iaobas+3) = gradzz(iaobas+3)+dpryz*(azz-3.d0)
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
      gradx(iaobas+4) = gradx(iaobas+4)-ax*prbxx+2.d0*prbx
      grady(iaobas+4) = grady(iaobas+4)-ay*prbxx
      gradz(iaobas+4) = gradz(iaobas+4)-az*prbxx
      gradx(iaobas+5) = gradx(iaobas+5)-ax*prbyy
      grady(iaobas+5) = grady(iaobas+5)-ay*prbyy+2.d0*prby
      gradz(iaobas+5) = gradz(iaobas+5)-az*prbyy
      gradx(iaobas+6) = gradx(iaobas+6)-ax*prbzz
      grady(iaobas+6) = grady(iaobas+6)-ay*prbzz
      gradz(iaobas+6) = gradz(iaobas+6)-az*prbzz+2.d0*prbz
c
      dprxy = dprxy*sqrt3i
      dprxz = dprxz*sqrt3i
      dpryz = dpryz*sqrt3i
      gradxx(iaobas+4) = gradxx(iaobas+4)+prb*(2.d0+axx*(axx-5.d0))
      gradxy(iaobas+4) = gradxy(iaobas+4)+dprxy*(axx-2.d0)
      gradxz(iaobas+4) = gradxz(iaobas+4)+dprxz*(axx-2.d0)
      gradyy(iaobas+4) = gradyy(iaobas+4)+prb*axx*(ayy-1.d0)
      gradyz(iaobas+4) = gradyz(iaobas+4)+dpryz*axx
      gradzz(iaobas+4) = gradzz(iaobas+4)+prb*axx*(azz-1.d0)
      gradxx(iaobas+5) = gradxx(iaobas+5)+prb*ayy*(axx-1.d0)
      gradxy(iaobas+5) = gradxy(iaobas+5)+dprxy*(ayy-2.d0)
      gradxz(iaobas+5) = gradxz(iaobas+5)+dprxz*ayy
      gradyy(iaobas+5) = gradyy(iaobas+5)+prb*(2.d0+ayy*(ayy-5.d0))
      gradyz(iaobas+5) = gradyz(iaobas+5)+dpryz*(ayy-2.d0)
      gradzz(iaobas+5) = gradzz(iaobas+5)+prb*ayy*(azz-1.d0)
      gradxx(iaobas+6) = gradxx(iaobas+6)+prb*azz*(axx-1.d0)
      gradxy(iaobas+6) = gradxy(iaobas+6)+dprxy*azz
      gradxz(iaobas+6) = gradxz(iaobas+6)+dprxz*(azz-2.d0)
      gradyy(iaobas+6) = gradyy(iaobas+6)+prb*azz*(ayy-1.d0)
      gradyz(iaobas+6) = gradyz(iaobas+6)+dpryz*(azz-2.d0)
      gradzz(iaobas+6) = gradzz(iaobas+6)+prb*(2.d0+azz*(azz-5.d0))
          endif
   10   continue
   20 continue
c
      return
      end
c
      subroutine aovlnw(x0,y0,z0,valao,gradx,grady,gradz,gradxx,
     + gradxy,gradxz,gradyy,gradyz,gradzz)
c**********************************************************************
c
c calculate the Laplacian of all ao's in point (x0,y0,z0), and return 
c   them in the array grads
c
c *******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(nmxl = 6)
      parameter(explim = 250.d0)
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
     + gradz(nbasis),gradxx(nbasis),gradxy(nbasis),gradxz(nbasis),
     + gradyy(nbasis),gradyz(nbasis),gradzz(nbasis)
      dimension cosphi(nmxl),sinphi(nmxl),
     +  pnl((nmxl+1)*(nmxl+2)/2),dpnl((nmxl+1)*(nmxl+2)/2),
     +  d2pnl((nmxl+1)*(nmxl+2)/2),tmat(9),dtmat(27),
     +  fodsph(3),sodsph(6),dao(3),d2ao(6),newm(9)
c
      valao(1:nbasis)=0.d0
      gradx(1:nbasis)=0.d0
      grady(1:nbasis)=0.d0
      gradz(1:nbasis)=0.d0
      gradxx(1:nbasis)=0.d0
      gradxy(1:nbasis)=0.d0
      gradxz(1:nbasis)=0.d0
      gradyy(1:nbasis)=0.d0
      gradyz(1:nbasis)=0.d0
      gradzz(1:nbasis)=0.d0
c
      pi = acos(-1.d0)
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
            gssfn=exp(-expon)*ctran(izeta)
            alfrs=2.d0*zeta(izeta)*r**2
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
              gradxx(iaobas+nm) = gradxx(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*d2ao(1)
              gradxy(iaobas+nm) = gradxy(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*d2ao(2)
              gradxz(iaobas+nm) = gradxz(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*d2ao(3)
              gradyy(iaobas+nm) = gradyy(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*d2ao(4)
              gradyz(iaobas+nm) = gradyz(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*d2ao(5)
              gradzz(iaobas+nm) = gradzz(iaobas+nm) +
     + sqrt(4*pi/mnmbr)*d2ao(6)
            enddo
          endif
c
  110   continue
  120 continue
c
      return
      end
