      subroutine xcpot (nspin, n, nrow, rho, loclxc, xcpar, xpot, 
     + cpot)
c
c ======================================================================
c
c purpose: calculation of exchange correlation potential, due to the
c          density, in a number of points.
c          either the (classical) x-alpha method is chosen, or
c          the method due to vosko, wilk and nusair, with optionally
c          a correction according to stoll, pavlidou and preuss.
c
c input  : nspin  - nr. of independent spins
c          n      - nr. of points
c          nrow   - row dimension of array rho and pot, as specified
c                   in the calling routine
c          rho    - the charge density in n points, nspin components
c          loclxc - option parameter for the type of exchange correl.
c                   0: x-alpha
c                   1,2: vosko-wilk-nusair, vwn+stoll
c          xcpar  - x-alpha parameter (loclxc=0)
c
c output : pot    - exchange correlation potential function, nspin
c                   components
c
c scratch: w      - (n,12)
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    stoll
c
      parameter (xconst = -1.240700981799d0,
     +           restrc =  0.7937005259841d0,
     +           third  = 1d0 / 3d0,
     +           onehlf = 1.5d0)
c
      dimension  rho(nrow,nspin), xpot(nrow,nspin), cpot(nrow,nspin)
      dimension  w(n,12)
c
c ======================================================================
c
c the kohn-sham exchange potential
c
      if (nspin.eq.1) then
        do i = 1, n
          xpot(i,1) = xconst * restrc * rho(i,1)**third
        enddo
      else
        do i = 1, n
          xpot(i,1) = xconst * rho(i,1)**third
          xpot(i,2) = xconst * rho(i,2)**third
        enddo
      endif
c
c correlation part
c choose model: x-alpha, or vosko-wilk-nusair
c
      if (loclxc.eq.0) then
c
c x-alpha
c
        do is = 1, nspin
          xpot(1:n,is) = onehlf*xcpar*xpot(1:n,is)
        enddo
c
      elseif ((loclxc.ne.1).and.(loclxc.ne.2)) then
        write(6,'(''ERROR; unrecognized loclxc option in xcpot'')')
        stop
      else 
c
c vosko, wilk, nusair
c the result is returned (from xcvwnp) in the first 2 columns of w
c
        stoll = mod (loclxc,2) .eq.0
c
        call xcvwnp (nspin,  n, nrow, rho, stoll,
     +               w(1,1), w(1,3), w(1,4), w(1,5), w(1,6),
     +               w(1,7), w(1,8), w(1,9), w(1,10))
c
        do is = 1, nspin
          cpot(1:n,is) = cpot(1:n,is) + w(1:n,is)
        enddo
c
      endif
c
c
      return
      end
