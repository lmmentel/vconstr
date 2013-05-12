      subroutine xcprdw (lcalcv, nspin,  nx,     n,
     +                   rho,    drho,   rhotrd, relrho,
     +                   sum,    dsum,   d2sum,  sumtrd, sumsxt,
     +                   nilrho, e, v)
c
c ======================================================================
c
c purpose: evaluation of the gradient correlation correction function
c          according to: j.p. perdew, phys.rev.b 33 (1986) 8822.
c          formulas taken (rewritten) from: liangyou fan, ph.d.thesis,
c          calcary, alberta, canada (1992), page 62-71.
c
c input  : lcalcv - (logical): calculate the potential (if false: only
c                   the energy density)
c          nspin  - nr. of independent spins
c          nx     - size of arrays as declared in the calling program.
c          n      - nr. of points for which the items are computed
c          rho    - density, for all spins
c          drho   - partial derivatives of rho, with respet to xyz,
c                   for spinup, down
c          rhotrd - rho**(1/3), for spinup, down
c          relrho - rho(spin)/rho(total)
c          sum,dsum,... - similarly, for the total charge density
c          nilrho - (logical) output functions must be set to zero,
c                   because the charge density is in fact zero, irres-
c                   pective of the value contained in array rho.
c
c in-out : e,v    - energy density function, and (spin-dependent) pot.
c                   values computed here are added to the input values
c
c scratch: crho,dcdrho,denomi,dinv,gradsq,pade,phi,phitrm,
c          vindep,temp*,etemp - (n)
c          vspin        - (n,nspin)
c
c remark * all densities are assumed definite positive.
c        * if nspin=1, then rho should be half the total density.
c        * a few function calls for vector operations are used to avoid
c          compiler "optimization" re-arrangemens in the order of
c          factors, which might yield underflow in extreme cases.
c
c called from: xc.
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    lcalcv,  nilrho
c
      parameter (alpha = 0.023266d0,  beta  = 7.389d-6,
     +           gamma = 8.723d0,     delta = 0.472d0,     ff = 0.11d0,
     +           c0    = 0.001667d0,  c1    = 0.002568d0,  cinf = c0+c1,
     +           g     = 1.745d0 *ff*cinf)
c
c Constants
      parameter (explim = 200.d0)
c
      parameter (zero = 0d0,  one = 1d0,  two = 2d0,  three = 3d0,
     +           r1o3 = 1d0/3d0,  r4o3 = 4d0/3d0,
     +           r5o3 = 5d0/3d0,  r11o3 = 11d0/3d0, r5o6 = 5d0/6d0,
     +           r7o6 = 7d0/6d0,  r1e4  = 1d4)
c
      dimension  rho(n,nspin),      drho(n,3,nspin),
     +           rhotrd(n,nspin),   relrho(n,nspin),
     +           sum(n), dsum(n,3), d2sum(n,6),   sumtrd(n),
     +           sumsxt(n),         nilrho(n)
      dimension  e(n),      v(nx,nspin)
      dimension  crho(n),   dcdrho(n), denomi(n), dinv(n),
     +           etemp(n),  gradsq(n), pade(n),   phi(n),
     +           phitrm(n), temp1(n),  temp2(n),
     +           vindep(n), vspin(n,nspin)
c
c ======================================================================
c
      pi=acos(-1.0d0)
      fourpi=4.0d0*pi
c
      crs3 = three/fourpi
      crs = crs3 ** r1o3
      r2h1o3 = two ** r1o3
c
c pade=rational function (second term in c(rho))
c denomi=1/denominator of the rational function
c
      a2 = alpha *crs
      a1 = beta *crs**2
      b2 = gamma *crs
      b1 = delta *crs**2
      b0 = r1e4*beta *crs3
c
c store 1/rho**4/3 in temp2
c
      do 10 k = 1, n
        temp1(k) = c1*sum(k) + (a2*sumtrd(k) + a1) *sumtrd(k)
        denomi(k) = one/ (sum(k) + (b2*sumtrd(k) +b1) * sumtrd(k) + b0)
   10 continue
      call vmuvw (n, temp1, denomi, pade)
c
c abs(grad)**2, and x ( =abs(grad)/rho**4/3 )
c
      do 20 k = 1, n
        gradsq(k) = dsum(k,1)**2 + dsum(k,2)**2 + dsum(k,3)**2
        temp2(k) = one / sum(k) ** r4o3
        temp1(k) = sqrt (gradsq(k)) * temp2(k)
   20 continue
c
c c(rho) = c0+rational function
c phi=g*grad/(c(rho)*rho**7/6) = g*(grad/rho**4/3) *rho**1/6 /c(rho)
c
      do 30 k = 1, n
        crho(k) = c0 + pade(k)
        phi(k) = g * temp1(k) * sumsxt(k) / crho(k)
   30 continue
c
c d**-1 (polarization function)
c
      if (nspin.eq.1) then
        do 40 k = 1, n
          dinv(k) = one
   40   continue
      else
        do 50 k = 1, n
          if (relrho(k,1) .gt.zero) then
            dinv(k) = relrho(k,1) **r5o3
          else
            dinv(k) = zero
          endif
          if (relrho(k,2) .gt.zero) then
            dinv(k) = dinv(k) + relrho(k,2) **r5o3
          endif
   50   continue
      endif
c
      do 70 k = 1, n
        if (phi(k) .lt.explim) then
          temp1(k) = exp (-phi(k))
        else
          temp1(k) = zero
        endif
   70 continue
c
c energy density, apart from a factor grad(rho)**2
c
      do 80 k = 1, n
        if (.not.nilrho(k)) then
          if (nspin.eq.2) dinv(k) = one / (r2h1o3 * sqrt (dinv(k)))
          etemp(k) = dinv(k) * temp1(k) * crho(k) * temp2(k)
        else
          etemp(k) = zero
        endif
   80 continue
c
c ----------------------------------------------------------------------
c potential
c
      if (.not.lcalcv) goto 500
c
c spin-dependent part of the potential
c
      if (nspin.eq.1) then
        do 140 k = 1, n
          vspin(k,1) = zero
  140   continue
c
      else
c
        do 150 k = 1, n
          temp1(k) = rhotrd(k,1) **2
          temp2(k) = rhotrd(k,nspin) **2
          vspin(k,1) = -r5o6 * (temp1(k)-temp2(k)) /
     +             (temp1(k)*rho(k,1) + temp2(k)*rho(k,nspin))
          vspin(k,nspin) = -vspin(k,1)
  150   continue
c
        do 200 is = 1, nspin
          js = 3-is
c
c temp1=grad(other spin).dot.grad(grad)
c
          do 160 k = 1, n
            temp1(k) = drho(k,1,js)*dsum(k,1) + drho(k,2,js)*dsum(k,2) +
     +                 drho(k,3,js)*dsum(k,3)
  160     continue
c
          do 170 k = 1, n
            temp2(k) = (one-phi(k)) * relrho(k,js) *gradsq(k) -
     +                 (two-phi(k)) * temp1(k)
  170     continue
c
          do 180 k = 1, n
            vspin(k,is) = vspin(k,is) * temp2(k)
  180     continue
c
  200   continue
c
      endif
c
c d/drho (c(rho))
c
      const = three * r1e4*beta * crs**2
      do 300 k = 1, n
        temp1(k) = (alpha*sumtrd(k) + two*beta*crs)  *sumtrd(k)
        temp2(k) = (gamma*sumtrd(k) + two*delta*crs) *sumtrd(k) + const
  300 continue
c
      do 310 k = 1, n
        dcdrho(k) = -crs /(three*sum(k)) * denomi(k) *
     +              ( temp1(k) - temp2(k)*pade(k) )
  310 continue
c
c phitrm=(4/3 - 11phi/3 + 7phi**2/6) / rho
c
      do 320 k = 1, n
        phitrm(k) = ((r7o6*phi(k) - r11o3) * phi(k) + r4o3) /sum(k)
  320 continue
c
c laplacian
c grad*grad(grad)
c
      do 330 k = 1, n
        temp1(k) = d2sum(k,1)+d2sum(k,4)+d2sum(k,6)
c
        temp2(k) = dsum(k,1) * (dsum(k,1)*d2sum(k,1) + two*
     +    (dsum(k,2)*d2sum(k,2) + dsum(k,3)*d2sum(k,3))) +
     +    dsum(k,2) * (dsum(k,2)*d2sum(k,4) + two*dsum(k,3)*d2sum(k,5))+
     +    dsum(k,3)**2 *d2sum(k,6)
  330 continue
      do 340 k = 1, n
        if (gradsq(k) .gt.zero) temp2(k) = temp2(k) / gradsq(k)
  340 continue
c
c
c accumulate the spin-independent terms
c
      do 350 k = 1, n
        vindep(k) = (two-phi(k)) * temp1(k) +
     +              phi(k)*(phi(k)-three) * temp2(k) -
     +              gradsq(k) * ( phitrm(k) +
     +                (phi(k)**2-phi(k)-one) * dcdrho(k)/crho(k) )
  350 continue
c
c
c total potential
c
      do 410 is = 1, nspin
        do 400 k = 1, n
          v(k,is) = v(k,is) - etemp(k) * (vindep(k) + vspin(k,is))
  400   continue
  410 continue
c
c
c energy
c
  500 do 510 k = 1, n
        e(k) = e(k) +  etemp(k) * gradsq(k)
  510 continue
c
c
      return
      end
