      subroutine pw91x1 (np, lnlpot, d, s, u, v, nilrho, ex, vx)
c
c ======================================================================
c
c purpose: gga91 exchange for a spin-unpolarized electronic system
c
c input  : np - nr. of points
c          lnlpot - .true.: calculate potential as well
c          d(np) - density
c          s(np) -  abs(grad d)/(2*kf*d)
c          u(np) -   (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
c          v(np) - (laplacian d)/(d*(2*kf)**2)
c          nilrho(np) - as usual
c output : ex(np) - exchange energy per electron
c          vx(np) -
c
c remark * u,v input, and vx output, only when lnlpot = .true.
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    lnlpot, nilrho
c
      parameter (zero = 0d0,  one = 1d0,  two = 2d0,  three = 3d0,
     +           tol  = 1d-5)
c
      dimension d(np), s(np), u(np), v(np), ex(np), vx(np), nilrho(np)
c
      data a1,a2,a3,a4/0.19645d0,0.27430d0,0.15084d0,100.d0/
      data ax,a,b1/-0.7385588d0,7.7956d0,0.004d0/
      data thrd,thrd4/0.333333333333d0,1.33333333333d0/
c
c ======================================================================
c
      do 100 ip = 1, np
c
        if (nilrho(ip)) then
          ex(ip) = zero
          if (lnlpot) vx(ip) = zero
          goto 100
        endif
c
        fac = ax*d(ip)**thrd
        s2 = s(ip)*s(ip)
        s3 = s2*s(ip)
        s4 = s2*s2
        p0 = one  / sqrt (one + a * a * s2)
c
        tmp = a*s(ip) + one/p0
        if (abs (tmp-one) .gt. tol) then
          p1 = log (tmp)
        else
          tmp = tmp - one
          p1 = tmp - tmp**2/two + tmp**3/three
        endif
c
        p2 = exp(-a4*s2)
        p3 = one / (one + a1*s(ip)*p1 + b1*s4)
        p4 = one + a1 * s(ip) * p1 + (a2-a3*p2) * s2
        f = p3*p4
        ex(ip) = fac*f
c-pph
        ex(ip) = ex(ip) - fac
c
c local exchange option
c ex = fac
c
        if (lnlpot) then
c
c energy done. now the potential:
c
          p5 = b1*s2-(a2-a3*p2)
          p6 = a1*s(ip)*(p1+a*s(ip)*p0)
          p7 = two * (a2-a3*p2)+2.d0*a3*a4*s2*p2-4.d0*b1*s2*f
          fs = p3*(p3*p5*p6+p7)
          p8 = two * s(ip)*(b1-a3*a4*p2)
          p9 = a1*p1+a*a1*s(ip)*p0*(three-a*a*s2*p0*p0)
          p10 = 4.d0*a3*a4*s(ip)*p2*(two-a4*s2)-8.d0*b1*s(ip)*f-
     +          4.d0*b1*s3*fs
          p11 = -p3*p3*(a1*p1+a*a1*s(ip)*p0+4.d0*b1*s3)
          fss = p3*p3*(p5*p9+p6*p8) + two*p3*p5*p6*p11+p3*p10+p7*p11
          vx(ip) = fac*(thrd4*f-(u(ip)-thrd4*s3)*fss-v(ip)*fs)
c-pph
          vx(ip) = vx(ip) - fac*thrd4
        endif
 100  continue
c
c
      return
      end
