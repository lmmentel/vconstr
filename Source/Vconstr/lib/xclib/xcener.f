      subroutine xcener (nspin, n, nrow, rho, loclxc, xcpar, ex, ec)
c
c ======================================================================
c
c purpose: calculation of exchange correlation energy, due to the
c          density in a number of points.
c          optionally the (classical) x-alpha method is chosen, or
c          the method due to vosko, wilk and nusair, with optionally
c          a correction according to stoll, pavlidou and preuss.
c          returned is the energy function, including a factor rho.
c
c input  : nspin  - nr. of independent spins
c          n      - nr. of points
c          nrow   - row dimension of array rho, as specified in the
c                   calling routine.
c          rho    - the charge density in n points, nspin components
c          loclxc - option parameter for the type of exchange correl.
c                   0: x-alpha
c                   1,2: vosko-wilk-nusair, vwn+stoll
c          xcpar  - x-alpha parameter
c
c output : e      - energy function in each point
c
c scratch: w      - (n,9)
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical   stoll
c
      parameter (xconst = -0.9305257363491d0,
     +           restrc =  0.7937005259841d0,
     +           fthird = 4d0 / 3d0,
     +           onehlf = 1.5d0)
c
      dimension  ex(n), ec(n), rho(nrow,nspin)
      dimension  w(n,9)
c
c ======================================================================
c
c the kohn-sham exchange energy
c
      if (nspin.eq.1) then
        do 10 i = 1, n
          ex(i) = xconst * restrc * rho(i,1)**fthird
   10   continue
c
      else
        do 20 i = 1, n
          ex(i) = xconst * (rho(i,1)**fthird + rho(i,2)**fthird)
   20   continue
      endif
c
c correlation part
c choose model:  x-alpha, or vosko-wilk-nusair
c
      if (loclxc.eq.0) then
c
c x-alpha
c
        ex(1:n) = onehlf*xcpar*ex(1:n)
c
      else
c
c vosko,wilk,nusarir
c
c the result, from xcvwne, is returned in the first column of
c array w. it still has to be multiplied by the (total)
c electronic charge density, rho.
c
        stoll = mod (loclxc,2) .eq.0
c
        call xcvwne (nspin,  n, nrow, rho,   stoll,
     +               w(1,1), w(1,2), w(1,3), w(1,4), w(1,5),
     +               w(1,6), w(1,7))
c
        if (nspin.eq.1) then
          ec(1:n) = rho(1:n,1)*w(1:n,1)
        else
          ec(1:n) = (rho(1:n,1)+rho(1:n,2))*w(1:n,1)
        endif
c
      endif
c
c
      return
      end
