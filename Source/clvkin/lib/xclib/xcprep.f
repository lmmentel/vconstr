      subroutine xcprep (calc2d, nspin, nx, n, den,    dden,   d2den,
     +                   nilrho, rho,   drho,  d2rho,  rhotrd, relrho,
     +                   sum,    dsum,  d2sum, sumtrd, sumsxt)
c
c ======================================================================
c
c purpose: preparation for various xc routines: compute rho**(1/3),
c          rho**(1/6), ...
c
c input  : calc2d - (logical) flag for calculation of second derivatives
c          nspin  - nr. of independent spins
c          nx     - declaration size of some arrays
c          n      - nr. of points
c          den    - input charge density (nspin columns)
c          dden   -       partial derivatives (cartesian)
c          d2den  -       partial second derivatives
c
c output : nilrho - (logical) 0.d0 charge density
c                   n.b.: if rho is small, it is set to 1.d0, to
c                   avoid numerical problems in formulas
c          rho    - density, adapted from input-values when almost zero
c                   also, for nspin=1, the rho is half the input value
c                   (input is then the total density, rho is the spin-
c                   density, even if nspin=1)
c          drho   - partial derivatives, adapted
c          d2rho  - partial second derivatives, identical to input
c                   if nspin=2, half the input for nspin=1
c          rhotrd - rho**(1/3), for all spins
c          relrho - relative density: rho(spin)/rho
c          sum    - total density (sum of all spins)
c          dsum   - derivatives
c          d2sum  - second derivatives
c          sumtrd - sum**(1/3)
c          sumsxt - sum**(1/6)
c
c called from: xc.
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    calc2d, nilrho
c
      parameter (small = 1.d-25)
c
      parameter (r1o3 = 1d0/3d0,  r1o6 = 1d0/6d0)
c
      dimension  den(nx,nspin),      dden(nx,3,nspin),
     +           d2den(nx,6,nspin),  rho(n,nspin),
     +           drho(n,3,nspin),    d2rho(n,6,nspin),
     +           rhotrd(n,nspin),    relrho(n,nspin),
     +           nilrho(n,3),        sum(n),
     +           dsum(n,3),          d2sum(n,6),
     +           sumtrd(n),          sumsxt(n)
c
c ======================================================================
c
      if (nspin.eq.1) then
        call vmucv (n, den, 5.d-1, rho)
        do 1 ic = 1, 3
          call vmucv (n, dden(1,ic,1), 5.d-1, drho(1,ic,1))
    1   continue
c
        if (calc2d) then
          do 2 ic = 1, 6
            call vmucv (n, d2den(1,ic,1), 5.d-1, d2rho(1,ic,1))
    2     continue
        endif
c
      else
        do 5 is = 1, nspin
          call vcuv (n, den(1,is), rho(1,is))
          do 3 ic = 1, 3
            call vcuv (n, dden(1,ic,is), drho(1,ic,is))
    3     continue
c
          if (calc2d) then
            do 4 ic = 1, 6
              call vcuv (n, d2den(1,ic,is), d2rho(1,ic,is))
    4       continue
          endif
c
    5   continue
      endif
c
c
      sum(1:n) = 0.d0
      nilrho(1:n,1:3) = .false.
c
      do 20 is = 1, nspin
        do 10 k = 1, n
          if (den(k,is) .gt.small) then
            sum(k) = sum(k) + den(k,is)
          else
            nilrho(k,is) = .true.
            rho(k,is) = 1.d0
            drho(k,1,is) = 0.d0
            drho(k,2,is) = 0.d0
            drho(k,3,is) = 0.d0
          endif
   10   continue
   20 continue
c
c
      do 30 ic = 1, 3
        call vauvw (n, drho(1,ic,1), drho(1,ic,nspin), dsum(1,ic))
   30 continue
c
      if (nspin.eq.1) then
        do 40 k = 1, n
          nilrho(k,3) = nilrho(k,1)
   40   continue
c
        if (calc2d) then
          do 50 ic = 1, 6
            call vcuv (n, d2den(1,ic,1), d2sum(1,ic))
   50     continue
        endif
c
      else
        do 60 k = 1, n
          nilrho(k,3) = nilrho(k,1) .and. nilrho(k,nspin)
   60   continue
c
        if (calc2d) then
          do 70 ic = 1, 6
            call vauvw (n, d2den(1,ic,1), d2den(1,ic,2), d2sum(1,ic))
   70     continue
        endif
      endif
c
c
      do 90 k = 1, n
        if (nilrho(k,3)) sum(k) = 1.d0
   90 continue
c
      do 100 k = 1, n
        sumsxt(k) = sum(k) ** r1o6
        sumtrd(k) = sumsxt(k) ** 2
  100 continue
c
      if (nspin.eq.1) then
        factor = 5.d-1 ** r1o3
        do 120 k = 1, n
          relrho(k,1) = 5.d-1
          rhotrd(k,1) = factor * sumtrd(k)
  120   continue
c
      else
        do 140 k = 1, n
          relrho(k,1) = rho(k,1) / sum(k)
          relrho(k,2) = rho(k,2) / sum(k)
          rhotrd(k,1) = rho(k,1) ** r1o3
          rhotrd(k,2) = rho(k,2) ** r1o3
  140   continue
      endif
c
      do 160 is = 1, nspin
        do 150 k = 1, n
          if (nilrho(k,is)) relrho(k,is) = 0.d0
  150   continue
  160 continue
c
c
      return
      end
