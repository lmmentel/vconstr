      subroutine prcrds(lx,x,y,z,r,costht,sintht,cosphi,sinphi)
c ======================================================================
c
c purpose: prepare (spherical-) coordinate tables
c
c input  : lx,x,y,z
c output : r,costht,sintht,cosphi(lx),sinphi(lx)
c
c *=====================================================================
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension cosphi(lx),sinphi(lx)
c
      parameter(epsz = 1.d-25)
      parameter(zero = 0.d0)
      parameter(one  = 1.d0)
      parameter(eps2 = 1.d-6)
c
      rxy=sqrt(x*x+y*y)
      r=sqrt(x*x+y*y+z*z)
      if (r.lt.epsz) then
        x=epsz
        y=epsz
        z=epsz
        rxy=sqrt(x*x+y*y)
        r=sqrt(x*x+y*y+z*z)
      endif
      costht=z/r
      sintht=rxy/r
C
C TRANSFORMATION IS SINGULAR FOR THETA = ZERO.
C GRADIENT SHOULD BE A SMOOTH FUNCTION -> USE SMALL THETA
C
      if (abs(sintht).lt.eps2) sintht=eps2
      if (rxy.lt.epsz) then
        cosphi(1)=one
        sinphi(1)=zero
      else
        cosphi(1)=x/rxy
        sinphi(1)=y/rxy
      endif
c
      do 10 m=2,lx
        cosphi(m)=cosphi(m-1)*cosphi(1)-
     + sinphi(m-1)*sinphi(1)
        sinphi(m)=sinphi(m-1)*cosphi(1)+
     + cosphi(m-1)*sinphi(1)
 10   continue
c
      return
      end
c
      subroutine capsph(lx,costht,sintht,pnl,dpnl,d2pnl)
c ======================================================================
c
c purpose: calculation of the theta dependent part of spherical 
c
c                            m
c          harmonics c(l,m)*P (cos(theta)), and its first derivative
c                            l 
c
c          with respect to theta
c
c input  : lx     - maximum order of the spherical harmonics.
c          np,npx - number of points.
c cosht,sintht(np)
c
c output : pnl,dpnl,d2pnl
c          the value of the theta dependent part of the normalized
c          real spherical harmonics, and the first and second derivative
c          with respect to theta
c
c          y(l,m) is found in column (l+1)**2-l+m, so
c            y(0, 0) is found in column 1.
c            y(1, 0) is found in column 2.
c            y(1, 1) is found in column 3.
c            y(2, 0) is found in column 4.
c            y(2, 1) is found in column 5.
c            y(2, 2) is found in column 6.
c
c remarks:
c
c * NEGATIVE VALUES OF M INCLUDED FOR CONVENIENCE.
c * LEGENDRE POLYMIAL AND DERIVATIVES ARE MULTIPLIED BY NORMALIZATION
c   FACTOR.
c * THE LEGENDRE POLYNOMIALS FOUND WITH THE REcURSIVE FORMULAS DIFFERS
c   A FACTOR (-1)**M WITH THE DEFINITION USED IN WEISSBLUTH.
c   THIS IS COMPENSATED AT THE END.
c
c ----------
c see the book 'numerical recipes' by william h. press et. al.
c we take
c
c    l             l
c   y (theta,phi)=p (cos(theta)) for m.eq.0
c    m             0
c
c    l             l
c   y (theta,phi)=p (cos(theta))cos(m*phi) for m.gt.0
c    m             m
c
c    l             l
c   y (theta,phi)=p  (cos(theta))sin(-m*phi) for m.lt.0
c    m             -m
c
c we calculate p(m,m)(x) (x=cos(theta)) as
c
c   if (m.eq.0) p(m,m)(x)=p(0,0)(x)=1
c
c   if (m.gt.0) p(m,m)(x)=-(2m-1)sqrt(1-x**2)*p(m-1,m-1)(x)
c
c we calculate p(m+1,m)(x) as
c
c   p(m+1,m)(x)=x(2m+1)p(m,m)(x).
c
c we calculate p((m+2,m+3,..,l),m)(x) as
c
c             x(2l-1)p(l-1,m)(x)-(l+m-1)p(l-2,m)(x)
c   p(l,m)(x)=-------------------------------------
c                              l-m
c
c we calculates the norms using
c
c      1
c     /   l            2                  2   (l+m):
c     / (p (cos(theta)) d(cos(theta)) = -----*------
c     /   m                             2*l+1 (l-m):
c   -1
c
c   (n: denotes the factorial on n).
c
c   if (m.eq.0) then
c
c      2*pi
c      /
c      / 1 d(phi) = 2*pi
c      /
c     0
c
c   if (m.gt.0) then
c
c     2*pi                    2*pi
c     /    2                  /    2
c     / sin (m*phi) d(phi) =  / cos (m*phi) d(phi) = pi
c     /                       /
c    0                       0
c
c ======================================================================
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(zero = 0.0d0)
      parameter(one  = 1.0d0)
      parameter(two  = 2.0d0)
      parameter(four = 4.0d0)
c
      dimension pnl((lx+1)*(lx+2)/2),dpnl((lx+1)*(lx+2)/2),
     +  d2pnl((lx+1)*(lx+2)/2)
c
      ilm(l,m) = ( ((l+1)*(l+2))/2 - l + m  )
c
c ======================================================================
c
c     *** plm(0,0).
c
      pi = acos(-1.0d0)
      twopi =two*pi
      fourpi=four*pi
      j00=ilm(0,0)
      pnl(j00) = one
c
      dpnl(j00)  = zero
      d2pnl(j00) = zero
c
c     *** plm(m,m), m=1,2,..,lx.
c
      do 240 m=1,lx
        jm1m1=ilm(m-1,m-1)
        jmm=ilm(m,m)
        pnl(jmm)=pnl(jm1m1)*(2*m-1)*sintht
        dpnl(jmm)=(2*m-1)*costht*pnl(jm1m1)+
     + (2*m-1)*sintht*dpnl(jm1m1)
        d2pnl(jmm)=-(2*m-1)*(sintht*pnl(jm1m1)-
     + 2*costht*dpnl(jm1m1)-sintht*d2pnl(jm1m1))
  240 continue
c
c     *** plm(m+1,m), m=0,1,..,lx-1.
c
      do 340 m=0,lx-1
        jmm=ilm(m,m)
        jm1m=ilm(m+1,m)
        pnl(jm1m)=costht*(2*m+1)*pnl(jmm)
        dpnl(jm1m)=-(2*m+1)*sintht*pnl(jmm)+
     + costht*(2*m+1)*dpnl(jmm)
        d2pnl(jm1m)=(2*M+1)*(-costht*pnl(jmm)-
     + 2*sintht*dpnl(jmm)+costht*d2pnl(jmm))
  340 continue

c
c     *** plm(m+2,m+3,..,lx,m), m=0,1,..,lx-2.
c
      do 450 m=0,lx-2
        do 440 l=m+2,lx
          jlm=ilm(l,m)
          jl1m=ilm(l-1,m)
          jl2m=ilm(l-2,m)
          pnl(jlm)=(costht*(2*l-1)*pnl(jl1m)-
     + (l+m-1)*pnl(jl2m))/(l-m)
          dpnl(jlm)=(-(2*l-1)*sintht*pnl(jl1m)+
     + costht*(2*l-1)*dpnl(jl1m)-(l+m-1)*dpnl(jl2m))/(l-m)
          d2pnl(jlm)=(-(2*l-1)*costht*pnl(jl1m)-
     + 2*(2*l-1)*sintht*dpnl(jl1m)+
     + (2*l-1)*costht*d2pnl(jl1m)-(l+m-1)*d2pnl(jl2m))/(l-m)
  440   continue
  450 continue
c
c     *** normalize.
c
      do 520 l=0,lx
        jl0=ilm(l,0)
        ylmnrm=fourpi/(2*l+1)
        fctr=one/sqrt(ylmnrm)
        pnl(jl0)=fctr*pnl(jl0)
        dpnl(jl0)=fctr*dpnl(jl0)
        d2pnl(jl0)=fctr*d2pnl(jl0)
  520 continue
c
      do 550 l=1,lx
        ylmnrm=twopi/(2*l+1)
        do 540 m=1,l
          jlm=ilm(l,m)
          ylmnrm=ylmnrm*(l+m)*(l-m+1)
          fctr=one/sqrt(ylmnrm)
          pnl(jlm) = fctr*pnl(jlm)
          dpnl(jlm)=fctr*dpnl(jlm)
          d2pnl(jlm)=fctr*d2pnl(jlm)
  540   continue
  550 continue
c
      return
      end
c
      subroutine fsdsph(r,dr,d2r,th,dth,d2th,lx,m,cosphi,sinphi,tmat,
     +             dtmat,fodsph,sodsph,zlm,fodcar,sodcar)
c ======================================================================
c
c purpose: Calculate the derivatives of the spherical harmonics in 
c          cartesian coordinates
c
c input  : r,dr,d2r     - radial part and derivatives
c          th,dth,d2th  - theta dependent part and dervatives
c          l,lx         - angular quantum number, maximum
c          m            - magnetic quantum number
c          sinphi(lx),cosphi(lx)
c          tmat(3,3)    - transformation matrix
c          dtmat(3,3,3) - derivatives of transformation matrix
c output : zlm          - real spherical harmonic
c          fodcar(3)    - first order derivative in cartesian coord.
c          sodcar(6)    - second order derivative
c          fodsph(3)    - first order derivative in spherical coord.
c          sodsph(6)    - second order derivative
c
c *=====================================================================
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(zero = 0.0d0)
      parameter(one  = 1.0d0)
c
      dimension sinphi(lx),cosphi(lx),tmat(3,3),dtmat(3,3,3),
     +  fodsph(3),sodsph(6),fodcar(3),sodcar(6)
c
      iij(i,j)= ( ((i-1)*(6-i))/2 + j )
c
c     construct derivatives in spherical coordinates
c
      ph=one
      dph=zero
      if (m.lt.0) then
        ph =sinphi(-m)
        dph=-m*cosphi(-m)
      elseif (m.gt.0) then
        ph =cosphi(m)
        dph=-m*sinphi(m)
      endif
c
      zlm=r*th*ph
c
      fodsph(1)=dr*th*ph
      fodsph(2)=r*dth*ph
      fodsph(3)=r*th*dph
c
      d2ph=-m**2*ph
      sodsph(1)=d2r*th*ph
      sodsph(2)=dr*dth*ph
      sodsph(3)=dr*th*dph
      sodsph(4)=r*d2th*ph
      sodsph(5)=r*dth*dph
      sodsph(6)=r*th*d2ph
c
c     transform first order derivative to cartesian coord.
c
      fodcar(1:3)=0.d0
      do j=1,3
        do i=1,3
          if (i+j.ne.9) then
            fodcar(i)=fodcar(i)+tmat(i,j)*fodsph(j)
          endif
        enddo
      enddo
c
c     transform second order derivative to cartesian coord.
c
      sodcar(1:6)=0.d0
      do i=1,3
        do j=i,3
          ic=iij(i,j)
          do k=1,3
            do l=1,3
              if (k.gt.l) then
                is=iij(l,k)
              else
                is=iij(k,l)
              endif
                sodcar(ic)=sodcar(ic)+
     + tmat(i,k)*tmat(j,l)*sodsph(is)+
     + tmat(i,k)*dtmat(j,l,k)*fodsph(l)
            enddo
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine trfdrv(rad,cost,sint,cosp,sinp,tmat,dtmat)
c ======================================================================
c
c purpose:  transform the gradient in spherical cooordinates with
c           matrix t to cartesian coordinates.
c           the derivatives of this matrix with respect to r,theta, 
c           and phi are also calculated, since they are needed for the 
c           transformation of the second order derivatives.
c
c output :  tmat(3,3),dtmat(3,3,3)
c
c *=====================================================================
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero = 0.d0)
      parameter(one  = 1.d0)
      parameter(epsz = 1.d-25)
c
      dimension tmat(3,3),dtmat(3,3,3)
c
      rinv=one/rad
      costht=cost
      sintht=sint
      if (abs(sintht).lt.epsz) then
        csctht=zero
      else
        csctht=one/sintht
      endif
      cottht=costht*csctht
      cosphi=cosp
      sinphi=sinp
c
c     t
c
      tmat(1,1)=cosphi*sintht
      tmat(1,2)=cosphi*costht*rinv
      tmat(1,3)=-(csctht*sinphi*rinv)
      tmat(2,1)=sinphi*sintht
      tmat(2,2)=costht*sinphi*rinv
      tmat(2,3)=cosphi*csctht*rinv
      tmat(3,1)=costht
      tmat(3,2)=-(sintht*rinv)
      tmat(3,3)=zero
c
c
c     dt/dr
c
      dtmat(1,1,1)=zero
      dtmat(1,2,1)=-(cosphi*costht*rinv**2)
      dtmat(1,3,1)=csctht*sinphi*rinv**2
      dtmat(2,1,1)=zero
      dtmat(2,2,1)=-(costht*sinphi*rinv**2)
      dtmat(2,3,1)=-(cosphi*csctht*rinv**2)
      dtmat(3,1,1)=zero
      dtmat(3,2,1)=sintht*rinv**2
      dtmat(3,3,1)=zero
c
c     dt/dtheta
c
      dtmat(1,1,2)=cosphi*costht
      dtmat(1,2,2)=-(cosphi*sintht*rinv)
      dtmat(1,3,2)=cottht*csctht*sinphi*rinv
      dtmat(2,1,2)=costht*sinphi
      dtmat(2,2,2)=-(sinphi*sintht*rinv)
      dtmat(2,3,2)=-(cosphi*cottht*csctht*rinv)
      dtmat(3,1,2)=-sintht
      dtmat(3,2,2)=-(costht*rinv)
      dtmat(3,3,2)=zero
c
c     dt/dphi
c
      dtmat(1,1,3)=-(sinphi*sintht)
      dtmat(1,2,3)=-(costht*sinphi*rinv)
      dtmat(1,3,3)=-(cosphi*csctht*rinv)
      dtmat(2,1,3)=cosphi*sintht
      dtmat(2,2,3)=cosphi*costht*rinv
      dtmat(2,3,3)=-(csctht*sinphi*rinv)
      dtmat(3,1,3)=zero
      dtmat(3,2,3)=zero
      dtmat(3,3,3)=zero
c
      return
      end
