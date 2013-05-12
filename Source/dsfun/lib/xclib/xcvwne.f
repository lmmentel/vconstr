      subroutine xcvwne (nspin,  n,      nrow,   rho,    stoll,  e,
     +                   rhosxt, rhotrd, zeta,   fzeta,  gzeta,  w)
c
c ======================================================================
c
c purpose: calculation of correlation energy function, according to
c          vosko, wilk and nusair (can.j.ph. 58 (1980) 1200), with
c          optionally the correction according to stoll, pavlidou
c          and preuss  (theor. chim. acta 149 (1978) 143).
c
c input  : nspin - nr. of independent spins
c          n     - vector length
c          nrow  - row dimension of array rho, as specified in the
c                  calling routine
c          rho   - density in n points, nspin components
c          stoll - (logical) stoll correction applied
c
c output : e     - energy function
c
c scratch: rhosxt,rhotrd,zeta,fzeta,gzeta(n), w(n,3)
c
c remark * zeta is the spin polarization parameter:
c          zeta = (rho1-rho2) / (rho1+rho2)
c        * several functions of zeta are used:
c          fzeta = (1+zeta)**(4/3) + (1-zeta)**(4/3) - 2)/(2**(4/3)-2)
c          gzeta = zeta**4 * fzeta
c          hzeta = 9/8 * (2**(4/3)-2) * ( fzeta-gzeta )
c          hzeta is stored, when needed, in fzeta.
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    stoll
c
      parameter (zero  = 0d0,      one  = 1d0,    two    = 2d0,
     +           sixth = one/6d0,  half = 0.5d0,  fthird = 4d0/3d0,
     +           eps   = 1d-70)
c 2 ** -(1/6)
      parameter (factr1 = 0.8908987181403d0,
     +           factr2 = factr1*factr1)
c 1/ (2**(4/3) - 2)
      parameter (c1zeta = 1.923661050932d0,
     +           c2zeta = 9d0 / (8d0 * c1zeta) )
c
      dimension  e(n),            fzeta(n),  gzeta(n),
     +           rho(nrow,nspin), rhotrd(n), rhosxt(n),
     +           w(n,3),          zeta(n)
c
      data cap,       cbp,        ccp,     x0p,       bp,      cp, qp
     +    /0.0310907d0, -0.0311676d0, 1.24742d0, -0.104980d0, 3.72744d0,
     +     12.9352d0,   6.15199d0/,
     +     caf,        cbf,       ccf,     x0f,       bf,      cf, qf
     +    /0.01554535d0, -0.144601d0, 3.37666d0, -0.325000d0, 7.06042d0,
     +     18.0578d0,    4.73093d0/,
     +     caalf,       cbalf,        ccalf,    x0alf,  balf, calf, qalf
     +    /-0.01688685d0, -0.000414034d0, 0.317708d0, -0.00475840d0,
     +     1.13107d0,   13.0045d0,     7.12311d0/
c
c ======================================================================
c
      if (nspin.eq.1) then
c
c spin unpolarized
c
        do i = 1, n
          rhosxt(i) = rho(i,1)**sixth
          rhotrd(i) = rhosxt(i) * rhosxt(i)
        enddo
c
        call xcvwnf (n,   rhosxt, rhotrd, cap, cbp, ccp,
     +               x0p, bp,     cp,     qp,
     +               e,   w(1,1), w(1,2))
c
c stoll correction
c
        if (stoll) then
          do i = 1, n
            rhosxt(i) = rhosxt(i) * factr1
            rhotrd(i) = rhotrd(i) * factr2
          enddo
c
          call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,
     +                 w(1,3), w(1,1), w(1,2))
c
          do i = 1, n
            e(i) = e(i) - w(i,3)
          enddo
        endif
c
      else
c
c spin polarized
c
        do i = 1, n
          rhosxt(i) = (rho(i,1) + rho(i,2)) ** sixth
          rhotrd(i) = rhosxt(i) * rhosxt(i)
c
          zeta(i) = rho(i,1) + rho(i,2)
          if (zeta(i) .ne. zero)
     +      zeta(i) = (rho(i,1) - rho(i,2)) / (rho(i,1) + rho(i,2))
c
          w(i,1) = max (eps, one+zeta(i))
          w(i,2) = max (eps, one-zeta(i))
        enddo
c
        do i = 1, n
          fzeta(i) = ( w(i,1)**fthird + w(i,2)**fthird - two ) *c1zeta
          gzeta(i) = zeta(i)**4 * fzeta(i)
        enddo
c
        call xcvwnf (n,   rhosxt, rhotrd, cap, cbp, ccp,
     +               x0p, bp,     cp,     qp,
     +               e,   w(1,1), w(1,2))
c
        call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +               x0f,    bf,     cf,     qf,
     +               w(1,3), w(1,1), w(1,2))
c
        do i = 1, n
          e(i) = e(i) + gzeta(i) * (w(i,3) - e(i))
          fzeta(i) = c2zeta * (fzeta(i) - gzeta(i))
        enddo
c
        call xcvwnf (n,      rhosxt, rhotrd, caalf, cbalf, ccalf,
     +               x0alf,  balf,   calf,   qalf,
     +               w(1,3), w(1,1), w(1,2))
c
        e(1:n) = e(1:n) + fzeta(1:n) * w(1:n,3)
c
c stoll correction
c
        if (stoll) then
          do 120 i = 1, n
            rhosxt(i) = rho(i,1)**sixth
            rhotrd(i) = rhosxt(i) * rhosxt(i)
  120     continue
c
          call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,
     +                 w(1,3), w(1,1), w(1,2))
c
          do 200 i = 1, n
            e(i) = e(i) - half * (one+zeta(i)) * w(i,3)
            rhosxt(i) = rho(i,2)**sixth
            rhotrd(i) = rhosxt(i) * rhosxt(i)
  200     continue
c
          call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,
     +                 w(1,3), w(1,1), w(1,2))
c
          do 300 i = 1, n
            e(i) = e(i) - half * (one-zeta(i)) * w(i,3)
  300     continue
        endif
c
      endif
c
c
      return
      end
