      subroutine lyp (calcv,  nspin, nx, n, rho,   drho,   d2rho,
     +                relrho, sum,   dsum,  d2sum, sumtrd, nilrho,
     +                elyp,   vlyp)
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    calcv, nilrho
c
      dimension  rho    (n,nspin),     drho   (n,3,nspin),
     +           d2rho  (n,6,nspin),   relrho (n,nspin),
     +           nilrho (n,3),
     +           sum    (n),           dsum   (n,3),
     +           d2sum  (n,6),         sumtrd (n),
     +           elyp   (nx),          vlyp   (nx,nspin)
c
c ======================================================================
c
c purpose: Lee-Yang-Parr energy and (optionally) the potential
c
c input  : calcv  - (logical) compute the potential
c          nspin  - 1 or 2: nr of independent spins
c          nx     - leading dimension of target arrays (ELYP, VLYP)
c                   as declared in the caller
c          n      - vector length
c          rho    - spin densities. NB: for NSPIN=1 the first column
c                   contains HALF the total density, the second column
c                   is not used
c          drho   - partial first derivatives of RHO
c          d2rho  - partial second derivatives of RHO
c          relrho - relative spin densities: spin density/total density
c          sum    - total density
c          dsum   - partial first derivatives of SUM
c          d2sum  - partial second derivatives
c          sumtrd - SUM ** 1/3
c          nilrho - logicals to flag that results should not be computed
c                   for indicated elements of the vector
c
c in-out : elyp   - energy. MUST BE INITIALIZED ON INPUT: values are
c                   added here, not assigned
c          vlyp   - potentials. MUST ALSO BE INITIALIZED
c
c remark * potential formulas from C.Lee, W.Yang, R.G.Parr,
c          Phys.Rev.B 37,2 (1988) p785
c        * energy formula from B.G.Johnson et.al. J.Chem.Phys. 98,7
c          (1993) p5612 in order to use expression (for energy) that
c          does not involve second partial derivatives
c
c *=====================================================================
c
      parameter (r7     = 7d0, 
     +           r0p3   = 0.3d0,      r11    = 11d0,
     +           r36    = 36d0,       r47    = 47d0,
     +           r1o9   = 1d0/9d0,    r2o3   = 2d0/3d0,
     +           r4o3   = 4d0/3d0,    r5o3   = 5d0/3d0,
     +           r8o3   = 8d0/3d0,    r11o3  = 11d0/3d0,
     +           explim = 200.d0)
      parameter (ca     = 0.04918d0,  cb   = 0.132d0,
     +           cc     = 0.2533d0,   cd   = 0.349d0)
c
      dimension  gamd (3),    gamdd (3),
     +           ed   (3),    edd   (3),
     +           dd   (3),    ddd   (3),
     +           fd   (3),    fdd   (3),
     +           rd   (3),    rdd   (3),
     +           gd   (3)
c
      pi=acos(-1.0d0)
c
c
c ======================================================================
c
c     ------------------
c     computational part
c     ------------------
c
      test = 2.d0 * (cc / explim)**3
c
      fac2 = 2.d0 ** r2o3
      cfac = r0p3 * (3.d0 * pi**2) ** r2o3
      c1 = 2.d0**r11o3 * cfac * ca * cb
c
      do 100 i = 1, n
        if (nilrho(i,3)) goto 100
c
        suminv = 1.d0 / sum(i)
        rho43 = sum(i) * sumtrd(i)
c
        reldif = relrho(i,nspin) - relrho(i,1)
        relfac = relrho(i,1) * relrho(i,nspin)
c
        if (nspin.eq.1) then
          gam = 1.d0
          gama = 0.d0
          gamb = 0.d0
          do 5 j = 1, 3
            gamd(j) = 0.d0
            gamdd(j) = 0.d0
    5     continue
        else
c
          gam = 2.d0 * (1.d0 - relrho(i,1)**2 - relrho(i,nspin)**2)
          gama =  4.d0 * relrho(i,nspin) * reldif * suminv
          gamb = -4.d0 * relrho(i,1) * reldif * suminv
c
          if (calcv) then
            jj = 1
            do 10 j = 1, 3
              gamd(j) = 4.d0 * reldif * (relrho(i,nspin)*drho(i,j,1) -
     +                  relrho(i,1)*drho(i,j,nspin)) * suminv
              gamdd(j) = (relrho(i,nspin)*drho(i,j,1) - relrho(i,1)*
     +                   drho(i,j,nspin)) * 4.d0*suminv**2 *
     +                   (drho(i,j,nspin) - drho(i,j,1) -
     +                   3.d0*reldif*dsum(i,j)) + 4.d0*suminv*reldif *
     +                   (relrho(i,nspin)*d2rho(i,jj,1) -
     +                   relrho(i,1)*d2rho(i,jj,nspin))
              jj = min (jj+3, 6)
   10       continue
          endif
c
        endif
c
c
        if (calcv) rnabl = d2sum(i,1) + d2sum(i,4) + d2sum(i,6)
        grad2 = dsum(i,1)**2 + dsum(i,2)**2 + dsum(i,3)**2
c
        if (nilrho(i,1)) then
          r83a = 0.d0
          rnabla = 0.d0
          grada2 = 0.d0
        else
          r83a = rho(i,1)**r8o3
          if (calcv) rnabla = d2rho(i,1,1) + d2rho(i,4,1) + d2rho(i,6,1)
          grada2 = drho(i,1,1)**2 + drho(i,2,1)**2 + drho(i,3,1)**2
        endif
c
        if (nilrho(i,nspin)) then
          r83b = 0.d0
          rnablb = 0.d0
          gradb2 = 0.d0
        else
          r83b = rho(i,nspin)**r8o3
          if (calcv)
     +     rnablb = d2rho(i,1,nspin) +d2rho(i,4,nspin) +d2rho(i,6,nspin)
          gradb2 = drho(i,1,nspin)**2 + drho(i,2,nspin)**2 +
     +             drho(i,3,nspin)**2
        endif
c
        d = 1.d0 / (1.d0 + cd/sumtrd(i))
        da = d**2 * cd / (3.d0*sum(i)*sumtrd(i))
        db = da
c
        if (calcv) then
          jj = 1
          do 20 j = 1, 3
            dd(j) = da * dsum(i,j)
            ddd(j) = da * ( d2sum(i,jj) + dsum(i,j)**2 *
     +               (cd*d/sumtrd(i) - 2.d0) * r2o3*suminv )
            jj = min (jj+3, 6)
   20     continue
        endif
c
        f = gam * d
        fa = gama * d  +  gam * da
        fb = gamb * d  +  gam * db
c
        if (calcv) then
          do 30 j = 1, 3
            fd(j) = gamd(j) * d + gam * dd(j)
            fdd(j) = gamdd(j)*d + 2.d0*gamd(j)*dd(j) + gam*ddd(j)
   30     continue
        endif
c
        r = 1.d0 / sum(i)**r5o3
        ra = -r5o3 * r * suminv
        rb = ra
c
        if (calcv) then
          jj = 1
          do 40 j = 1, 3
            rd(j) = ra * dsum(i,j)
            rdd(j) = ra * (d2sum(i,jj) - r8o3 * dsum(i,j)**2 * suminv)
            jj = min (jj+3, 6)
   40     continue
        endif
c
c
        elyp(i) = elyp(i) - ca * f * sum(i)
c
        if (calcv) then
          vlyp(i,1) = vlyp(i,1) - ca * (fa*sum(i) + f)
          if (nspin.eq.2) vlyp(i,2) = vlyp(i,2) - ca * (fb*sum(i) + f)
        endif
c
c       ------------------------------------------------------------
c       remaining terms depend on exponential factor which falls off
c       rapidly as rho goes to zero
c       ------------------------------------------------------------
c
        if (sum(i) .lt. test) goto 100
c
        arg = cc / sumtrd(i)
        e = exp (-arg)
        ea = cc * e * suminv / (3.d0 * sumtrd(i))
        eb = ea
c
        if (calcv) then
          jj = 1
          do 50 j = 1, 3
            ed(j) = ea * dsum(i,j)
            edd(j) = ea * ( d2sum(i,jj) + dsum(i,j)**2 *
     +               (cc/sumtrd(i) - 4.d0) / (3.d0*sum(i)) )
            jj = min (jj+3, 6)
   50     continue
        endif
c
        g = f * e * r
        ga = fa * e * r  +  f * ea * r  +  f * e * ra
        gb = fb * e * r  +  f * eb * r  +  f * e * rb
c
        if (calcv) then
          gnabl = 0.d0
          do 60 j = 1, 3
            gd(j) = fd(j) * e * r  +  f * ed(j) * r  +  f * e * rd(j)
            gnabl = gnabl + fdd(j)*e*r + f*edd(j)*r + f*e*rdd(j) + 
     +              2.d0*(fd(j)*ed(j)*r+fd(j)*e*rd(j)+f*ed(j)*rd(j))
   60     continue
        endif
c
c
        t83 = r83a + r83b
c
c
        temp = suminv**2
c
        gaa = grada2 * temp
        gbb = gradb2 * temp
c
        gab = 0.d0
        if (.not.(nilrho(i,1) .or. nilrho(i,nspin)))
     +  gab = (drho(i,1,1)*drho(i,1,nspin) +
     +         drho(i,2,1)*drho(i,2,nspin) +
     +         drho(i,3,1)*drho(i,3,nspin)) * temp
c
c
        term83 = c1 * relfac * t83 / sum(i)**r8o3
c
        elyp(i) = elyp(i) - e * (d/sumtrd(i)) * rho43 * term83
c
        del = cc  +  cd * d
        d11 = del -  r11 * sumtrd(i)
c
        flfac = -ca * cb * e * d
        flaa = r1o9 * relfac * (sumtrd(i) - 3.d0*del - relrho(i,1) *
     +         d11) - relrho(i,nspin)**2 * sumtrd(i)
        flbb = r1o9 * relfac * (sumtrd(i) - 3.d0*del - relrho(i,nspin)*
     +         d11) - relrho(i,1)**2 * sumtrd(i)
        flab = r1o9 * relfac * (r47*sumtrd(i) - r7*del) -
     +         r4o3 * sumtrd(i)
c
        elyp(i) = elyp(i) + flfac * (flaa*gaa + flab*gab + flbb*gbb)
c
c       ---------
c       potential
c       ---------
c
        if (.not.calcv) goto 100
c
c
        t1 = 2.d0 * fac2 * ca * cb * cfac
        t2 = ca * cb / 4.d0
        tgd = gd(1)*dsum(i,1) + gd(2)*dsum(i,2) + gd(3)*dsum(i,3)
        t3 = 4.d0*g*rnabl
        tsrg = rnabl*sum(i) - grad2
        t4 = ca * cb / r36
c
        term1 = ga*t83 + r8o3*g*r83a/rho(i,1)
        term2 = sum(i)*gnabl + 4.d0*tgd
        term3 = 3.d0*rho(i,1)*gnabl
        trg = drho(i,1,1)*gd(1) + drho(i,2,1)*gd(2) + drho(i,3,1)*gd(3)
        term4 = rho(i,1)*rnabla + rho(i,nspin)*rnablb
        term5 = grada2 + gradb2
c
        vlyp(i,1) = vlyp(i,1) -
     +              t1 * term1 -
     +              t2 * (term2 + t3 + ga*tsrg) -
     +              t4 * ( term3 + 4.d0*trg + 4.d0*g*rnabla +
     +                     3.d0*ga*term4 + ga*term5 )
c
        if (nspin.eq.1) goto 100
c
        term1 = gb*t83 + r8o3*g*r83b/rho(i,nspin)
        term3 = 3.d0*rho(i,nspin)*gnabl
        trg = drho(i,1,nspin)*gd(1) + drho(i,2,nspin)*gd(2) +
     +        drho(i,3,nspin)*gd(3)
c
        vlyp(i,nspin) = vlyp(i,nspin) -
     +                  t1 * term1 -
     +                  t2 * (term2 + t3 + gb*tsrg) -
     +                  t4 * ( term3 + 4.d0*trg + 4.d0*g*rnablb +
     +                         3.d0*gb*term4 + gb*term5 )
c
  100 continue
c
c
      return
      end
