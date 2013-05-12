      subroutine xcbeck (calcv,  nspin,  nx, n, rho,   drho, d2rho,
     +                   rhotrd, nilrho, e,  v)
c
c ======================================================================
c
c purpose: evaluate the energy density and (optionally) the potential
c          according to beckes gradient correction formula.
c          a.d.becke
c
c input  : calcv  - (logical): calculate the potential (if false: only
c                   the energy density)
c          nspin  - nr. of independent spins
c          nx     - size of array (max nr. of points)
c          n      - actual nr. of points for which the functions have
c                   to be evaluated.
c          rho    - charge density (for all spins)
c          drho   - partial (cartesian) derivatives of rho.
c          d2rho  - partial second derivatives
c          rhotrd - rho**1/3 (for all spins)
c          nilrho - (logicals) function result must be set to zero
c                   (the density is in fact zero, irrespective of the
c                   value contained in array rho)
c
c in-out : e,v    - energy and potential (the latter spin-dependent)
c                   the values computed in this routine are added to
c                   the input values.
c
c scratch: b6term,f,grad,rho43,rho53,sinhix,sqrt1x,term1,term2,term3,
c          tlapl,x,xsq  - (n)
c
c remark * if nspin=1, then rho should be half the total charge density.
c
c called from: xc.
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    calcv, nilrho
c
      parameter (zero = 0d0,       one   = 1d0,
     +           two  = 2d0,       three = 3d0,
     +           r4o3 = 4d0/3d0,   tol   = 1d-5,
     +           b    = 0.0042d0,  b6    = 6d0*b)
c
      dimension  rho(n,nspin),    drho(n,3,nspin), d2rho(n,6,nspin),
     +           rhotrd(n,nspin), nilrho(n,nspin), e(n), v(nx,nspin)
      dimension  b6term(n), f(n),      grad(n),  rho43(n), rho53(n),
     +           sinhix(n), sqrt1x(n), term1(n), term2(n), term3(n),
     +           tlapl(n),  x(n),      xsq(n)
c
c ======================================================================
c
      do 200 is = 1, nspin
c
        do 10 k = 1, n
          rho43(k) = rho(k,is) * rhotrd(k,is)
          rho53(k) = rho43(k) * rhotrd(k,is)
          grad(k) = sqrt (drho(k,1,is)**2 +drho(k,2,is)**2 +
     +                    drho(k,3,is)**2)
          x(k) = grad(k) / rho43(k)
          xsq(k) = x(k) **2
          sqrt1x(k) = sqrt (one + xsq(k))
c
          tmp = x(k) + sqrt1x(k)
          if (abs (tmp-one) .gt. tol) then
            sinhix(k) = log (tmp)
          else
            tmp = tmp - one
            sinhix(k) = tmp - tmp**2/two + tmp**3/three
          endif
c
          b6term(k) = b6 * xsq(k) / sqrt1x(k)
          f(k) = one / (one + b6*x(k)*sinhix(k))
   10   continue
c
c
        factor = -b
        if (nspin.eq.1) factor = two*factor
        do 20 k = 1, n
          if (.not.nilrho(k,is))
     +      e(k) = e(k) +  factor * rho43(k) * xsq(k) * f(k)
   20   continue
c
        if (.not.calcv) goto 200
c
        do 100 k = 1, n
          term1(k) = drho(k,1,is) * (drho(k,1,is)*d2rho(k,1,is) + two*
     +      (drho(k,2,is)*d2rho(k,2,is) + drho(k,3,is)*d2rho(k,3,is))) +
     +       drho(k,2,is) * (drho(k,2,is)*d2rho(k,4,is) +
     +                       two*drho(k,3,is)*d2rho(k,5,is)) +
     +       drho(k,3,is)**2 *d2rho(k,6,is)
  100   continue
        do 105 k = 1, n
          if (term1(k) .ne.zero) term1(k) = term1(k) / grad(k)
  105   continue
c
        do 110 k = 1, n
          term2(k) = -r4o3 * rho53(k) * x(k)*xsq(k) +
     +                term1(k) / rho43(k)
  110   continue
        do 120 k = 1, n
          term1(k) = b6*f(k)*term2(k)
  120   continue
        do 130 k = 1, n
          term2(k) = one / (one+xsq(k)) +
     +               two*f(k) * (two - b6term(k))
  130   continue
c
        do 140 k = 1, n
          term3(k) = term1(k) * ( (one+two*f(k)) *sinhix(k) +
     +               term2(k) *x(k)/ sqrt1x(k) )
  140   continue
        do 150 k = 1, n
          term1(k) = r4o3 * rho53(k) * xsq(k)
          tlapl(k) = d2rho(k,1,is) + d2rho(k,4,is) + d2rho(k,6,is)
          term2(k) = tlapl(k) * (one + f(k) - f(k)*b6term(k))
  150   continue
c
        do 160 k = 1, n
          if (.not.nilrho(k,is))
     +      v(k,is) = v(k,is)  - b * f(k) / rho43(k) *
     +                (term1(k)-term2(k)+term3(k))
  160   continue
c
  200 continue
c
c
      return
      end
