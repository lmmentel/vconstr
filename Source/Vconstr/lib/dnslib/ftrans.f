      subroutine spherw
c***********************************************************************
c calculation of transformation matrices for the transformation
c of cartesian to spherical gaussians.
c
c matrices calculated for s,p,d,f and g functions
c
c the atmol4 (or integw) ordering of the basis functions is used
c (m=0, m=1, m=-1, m=2, m=-2, ...etc).
c
c it is assumed that the cartesian components are ordered as :
c
c s : s
c p : x, y, z
c d : x2, y2, z2, xy, xz, yz       (x2 is x squared)
c f : x3, y3, z3, x2y, x2z, xy2, y2z, xz2, yz2, xyz
c g : x4, y4, z4, x3y, x3z, xy3, y3z, xz3, yz3, x2y2, x2z2, y2z2,
c     x2yz, xy2z, xyz2
c
c all cartesian components belonging to the same group must have the
c same normalization coefficient n. the value n must be chosen in
c such a way that the functions with highest powers of z (p(z), d(z2),
c f(z3) and g(z4)) are normalized to 1.
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero   = 0.0d0)
      parameter(half   = 0.5d0)
      parameter(one    = 1.0d0)
      parameter(onhalf = 1.5d0)
      parameter(three  = 3.0d0)
      parameter(four   = 4.0d0)
      parameter(five   = 5.0d0)
      parameter(six    = 6.0d0)
      parameter(eight  = 8.0d0)
      parameter(fiftn  = 15.0d0)
      parameter(fourty = 40.0d0)
      common/spheri/s(1,1),p(3,3),d(5,6),f(7,10),g(9,15),istrt(0:4)
c
      p(1:3,1:3)=zero
      d(1:5,1:6)=zero
      f(1:7,1:10)=zero
      g(1:9,1:15)=zero
c
      istrt(0) = 1
      istrt(1) = 1+1
      istrt(2) = 1+9+1
      istrt(3) = 1+9+30+1
      istrt(4) = 1+9+30+70+1
c
      s(1,1) = one
c
c*** p(0) = p(z) ;  p(1) = p(x) ;  p(-1) = p(y)

      p(1,3) = one
      p(2,1) = one
      p(3,2) = one
c
c*** d(0) = z2-0.5x2-0.5y2
c

      d(1,1) = -half
      d(1,2) = -half
      d(1,3) = one
c
c*** d(1) = sqrt(3) xz ;  d(-1) = sqrt(3) yz
c
      d(2,5) = sqrt(three)
      d(3,6) = sqrt(three)
c
c*** d(2) = 0.5 sqrt(3) (x2-y2)
c
      d(4,1) = half*sqrt(three)
      d(4,2) = -half*sqrt(three)
c
c*** d(-2) = sqrt(3) xy
c
      d(5,4) = sqrt(three)
c
c*** f(0) = 1/2 (5z3 - 3r2z)
c
      f(1,3) = one
      f(1,5) = -onhalf
      f(1,7) = -onhalf
c
c*** f(1) = 1/sqrt(24) (15xz2-3xr2)
c
      f(2,1) = -sqrt(three/eight)
      f(2,6) = -sqrt(three/eight)
      f(2,8) =  sqrt(six)
c
c*** f(-1) = 1/sqrt(24) (15yz2-3yr2)
c
      f(3,2) = -sqrt(three/eight)
      f(3,4) = -sqrt(three/eight)
      f(3,9) =  sqrt(six)
c
c*** f(2) = 1/sqrt(60) (15zx2-15zy2)
c
      f(4,5) = sqrt(fiftn/four)
      f(4,7) = -sqrt(fiftn/four)
c
c*** f(-2) = 1/sqrt(60) (30xyz)
c
      f(5,10) = sqrt(fiftn)
c
c*** f(3) = 1/sqrt(360) (1.0)-3y2x)
c
      f(6,1) = sqrt(five/eight)
      f(6,6) = -sqrt(one/fourty)
c
c*** f(-3) = 1/sqrt(360) (3yx2-15y3)
c
      f(7,4) = sqrt(one/fourty)
      f(7,2) = -sqrt(five/eight)
c
      return
      end
c
      subroutine spherv(idddd)
c***********************************************************************
c transformation matrix for the transfo of cartesian d functions
c ("6 d set" : d(x2),d(y2),d(z2),d(xy),d(xz),d(yz) ) to spherical
c d functions ("5 d set" : d(xy),d(xz),d(yz),d(x2-y2),d(3z2-r2) )
c
c it is assumed that all six cartesian d functions have the same
c normalization coefficient n. with this value of n, the functions
c d(xy), d(xz) and d(yz) must be normalized to unity.
c
c transformation coefficients are stored in array d of common/spheri/
c
c if no transformation from 6-d set to 5-d set is necessary (idddd = 0
c in this case) then the extra normalization coefficients for the
c functions d(x2), d(y2) and d(z2) are stored.
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero  = 0.0d0)
      parameter(one   = 1.0d0)
      parameter(two   = 2.0d0)
      parameter(three = 3.0d0)
      parameter(four  = 4.0d0)
      dimension d5(5,6)
      common/spheri/s(1,1),p(3,3),d(6,6),dummy(199),istrt(0:4)
      equivalence(d5(1,1),d(1,1))
c
      p(1:3,1:3)=zero
      d(1:6,1:6)=zero
c
      istrt(0) = 1
      istrt(1) = 1+1
      istrt(2) = 1+9+1
c
      s(1,1) = one
      p(1,1) = one
      p(2,2) = one
      p(3,3) = one
c
      if (idddd.eq.0) then
c
c store normalization coefficients
c
c*** d(xy), d(xz), d(yz), d(x2), d(y2), d(z2)
c
        d(1,4) = one
        d(2,5) = one
        d(3,6) = one
        d(4,1) = one/sqrt(three)
        d(5,2) = one/sqrt(three)
        d(6,3) = one/sqrt(three)
      else
c
c store transformation matrix. used array d5 (equivalenced with
c array d in common/spheri/ )
c
c*** d(xy), d(xz), d(yz), d(x2-y2), d(3z2-r2)
c
        d5(1,4) = one
        d5(2,5) = one
        d5(3,6) = one
        d5(4,1) = sqrt(one/four)
        d5(4,2) = -sqrt(one/four)
        d5(5,1) = -one/(two*sqrt(three))
        d5(5,2) = -one/(two*sqrt(three))
        d5(5,3) =  one/sqrt(three)
      endif
c
      return
      end
c
      subroutine carsph(xint,lqi,lqj,ncomp,lqntmx)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension xint(15,15),ncomp(0:lqntmx)
      dimension xtran(9,9),dummy(245)
      common/spheri/s(1,1),p(3,3),d(5,6),f(7,10),g(9,15),istrt(0:4)
      equivalence(dummy(1),s(1,1))
c
      ncari = ncomp(lqi)
      ncarj = ncomp(lqj)
      nsphi = 2*lqi+1
      nsphj = 2*lqj+1
c
      call mmm(dummy(istrt(lqi)),nsphi,ncari,xint(1,1),
     + dummy(istrt(lqj)),nsphj,ncarj,xtran(1,1))
c
      xint(1:nsphi,1:nsphj)=xtran(1:nsphi,1:nsphj)
c
      return
      end
c
      subroutine cart6d(xint,ncomp,lqi,lqj,lqntmx)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension xint(15,15),ncomp(0:lqntmx)
      dimension xtran(9,9),dummy(245)
      common/spheri/s(1,1),p(3,3),d(6,6),dddum(199),istrt(0:4)
      equivalence(dummy(1),s(1,1))
c
      ncari = ncomp(lqi)
      ncarj = ncomp(lqj)
      nsphi = 2*lqi+1
      nsphj = 2*lqj+1
      if (lqi.eq.2) nsphi = 6
      if (lqj.eq.2) nsphj = 6
c
      call mmm(dummy(istrt(lqi)),nsphi,ncari,xint(1,1),
     + dummy(istrt(lqj)),nsphj,ncarj,xtran(1,1))
c
      xint(1:nsphi,1:nsphj)=xtran(1:nsphi,1:nsphj)
c
      return
      end
c
      subroutine mmm(a,na1,na2,b,c,nc1,nc2,x)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero = 0.0d0)
      dimension a(na1,na2),b(15,15),c(nc1,nc2),x(9,9)
c
      do ia1 = 1,na1
        do ic1 = 1,nc1
          x(ia1,ic1) = zero
          do ia2 = 1,na2
            do ic2 = 1,nc2
              x(ia1,ic1) = x(ia1,ic1)+a(ia1,ia2)*b(ia2,ic2)*
     + c(ic1,ic2)
            enddo
          enddo
        enddo
      enddo
c
      return
      end
