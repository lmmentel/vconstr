      subroutine xcvwnp (nspin,  n, nrow, rho, stoll, p,
     +                   rhosxt, rhotrd, zeta, fzeta, gzeta, hzeta,
     +                   dzeta,  w)
c
c ======================================================================
c
c purpose: calculation correlation potential function, according to
c          vosko, wilk and nusair (can.j.ph. 58 (1980) 1200), with
c          optionally the correction according to stoll, pavlidou
c          and preuss (theor. chim. acta 149 (1978) 143).
c
c input  : nspin  - nr. of independent spins
c          n      - vector length
c          nrow   - row dimension of array rho, as specified in the
c                   calling routine
c          rho    - density in n points, nspin components
c
c output : p      - potential function, for each spin
c
c scratch: rhosxt - (n)
c          rhotrd - (n)
c          zeta   - (n)
c          fzeta  - (n)
c          gzeta  - (n)
c          hzeta  - (n)
c          dzeta  - (n)
c          w      - (n,3)
c
c remark * zeta is the spin polarization parameter:
c             zeta = (rho1-rho2) / (rho1+rho2)
c        * several functions of zeta are used:
c          fzeta = (1+zeta)**(4/3) + (1-zeta)**(4/3) - 2)/(2**(4/3)-2)
c          gzeta = zeta**4 * fzeta
c          hzeta = 9/8 * (2**(4/3)-2) * ( fzeta-gzeta )
c          and the derivates of h and g with respect to zeta
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    stoll
c
      parameter (fthird = 4d0 / 3d0,
     +           one    = 1d0,
     +           sixth  = 1d0 / 6d0,
     +           third  = 1d0 / 3d0,
     +           two    = 2d0,
     +           eps    = 1d-70)
c 2 ** -(1/6)
      parameter (factr1 = 0.8908987181403d0,
     +           factr2 = factr1*factr1)
c 1/ (2**(4/3) - 2)
      parameter (c1zeta = 1.923661050932d0,
     +           c2zeta = 9d0 / (8d0 * c1zeta) )
c
      dimension  rho(nrow,nspin), rhotrd(n), rhosxt(n),
     +           p(n,nspin),
     +           zeta(n), fzeta(n), gzeta(n), hzeta(n), dzeta(n),
     +           w(n,3)
c
      data cap,       cbp,        ccp,     x0p,       bp,      cp, qp
     +    /0.0310907d0, -0.0311676d0, 1.24742d0, -0.104980d0, 3.72744d0,
     +     12.9352d0,   6.15199d0/,
     +     caf,        cbf,       ccf,     x0f,       bf,      cf, qf
     +    /0.01554535d0, -0.144601d0, 3.37666d0, -0.325000d0, 7.06042d0,
     +     18.0578d0,    4.73093d0/,
     +     caalf,       cbalf,        ccalf,    x0alf,  balf, calf, qalf
     +    /-0.01688685d0, -0.000414034d0, 0.317708d0, -0.00475840d0,
     +     1.13107d0,    13.0045d0,     7.12311d0/
c
c ======================================================================
c
c spin restricted
c
      if (nspin.eq.1) then
c
        do 10 i = 1, n
          rhosxt(i) = rho(i,1)**sixth
          rhotrd(i) = rhosxt(i) **2
   10   continue
c
        call xcvwnf (n,   rhosxt, rhotrd, cap, cbp, ccp,
     +               x0p, bp,     cp,     qp,
     +               p,   w(1,1), w(1,2))
c
        call xcvwnd (n,      rhosxt, rhotrd, cap, cbp, ccp,
     +               x0p,    bp,     cp,     qp,  w(1,3))
c
        do 20 i = 1, n
          p(i,1) = p(i,1) - w(i,3)
   20   continue
c
c stoll correction
c
        if (stoll) then
          do 30 i = 1, n
            rhosxt(i) = factr1 * rhosxt(i)
            rhotrd(i) = factr2 * rhotrd(i)
   30     continue
c
          call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,
     +                 w(1,3), w(1,1), w(1,2))
c
          do 40 i = 1, n
            p(i,1) = p(i,1) - w(i,3)
   40     continue
c
          call xcvwnd (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,  w(1,3))
c
          do 50 i = 1, n
            p(i,1) = p(i,1) + w(i,3)
   50     continue
        endif
c
      else
c
c spin polarized
c
        do 80 i = 1, n
          rhosxt(i) = ( rho(i,1) + rho(i,2) )**sixth
          rhotrd(i) = rhosxt(i) * rhosxt(i)
c
          zeta(i) = (rho(i,1) - rho(i,2)) / (rho(i,1) + rho(i,2))
c
          w(i,1) = max (eps, one+zeta(i))
          w(i,2) = max (eps, one-zeta(i))
   80   continue
c
        do 100 i = 1, n
          dzeta(i) = w(i,1)**fthird + w(i,2)**fthird
          fzeta(i) = ( dzeta(i) - two ) * c1zeta
          gzeta(i) = zeta(i)**4 * fzeta(i)
  100   continue
c
        call xcvwnd (n,      rhosxt, rhotrd, cap, cbp, ccp,
     +               x0p,    bp,     cp,     qp,  p)
c
        call xcvwnd (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +               x0f,    bf,     cf,     qf,  w(1,3))
c
        do 130 i = 1, n
          p(i,1) = -p(i,1) + gzeta(i) * (p(i,1) - w(i,3))
          hzeta(i) = c2zeta * (fzeta(i) - gzeta(i))
  130   continue
c
        call xcvwnd (n,      rhosxt, rhotrd, caalf, cbalf, ccalf,
     +               x0alf,  balf,   calf,   qalf,  w(1,3))
c
        do 140 i = 1, n
          p(i,1) = p(i,1) - hzeta(i) * w(i,3)
  140   continue
c
c derivatives of the zeta functions h and g, first the
c derivative of f, stored into d, then the derivative of g
c into f, and the derivative of h into d. so g and h are
c preserved, and f and d contain their derivatives
c
        do 160 i = 1, n
          p(i,2) = p(i,1)
c
          w(i,1) = max (eps, one+zeta(i))
          w(i,2) = max (eps, one-zeta(i))
c
          dzeta(i) = c1zeta * fthird * ( w(i,1)**third -
     +                                   w(i,2)**third )
          fzeta(i) = zeta(i)**3 * (4 * fzeta(i) + zeta(i) * dzeta(i))
          dzeta(i) = c2zeta * (dzeta(i) - fzeta(i))
  160   continue
c
        call xcvwnf (n,      rhosxt, rhotrd, cap, cbp, ccp,
     +               x0p,    bp,     cp,     qp,
     +               w(1,3), w(1,1), w(1,2))
c
        do 180 i = 1, n
          p(i,1) = p(i,1) + w(i,3) *
     +             (one - gzeta(i) - (-zeta(i)+one) * fzeta(i))
          p(i,2) = p(i,2) + w(i,3) *
     +             (one - gzeta(i) - (-zeta(i)-one) * fzeta(i))
  180   continue
c
        call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +               x0f,    bf,     cf,     qf,
     +               w(1,3), w(1,1), w(1,2))
c
        do 210 i = 1, n
          p(i,1) = p(i,1) + w(i,3) *
     +             (gzeta(i) + (-zeta(i)+one) * fzeta(i))
          p(i,2) = p(i,2) + w(i,3) *
     +             (gzeta(i) + (-zeta(i)-one) * fzeta(i))
  210   continue
c
        call xcvwnf (n,      rhosxt, rhotrd, caalf, cbalf, ccalf,
     +               x0alf,  balf,   calf,   qalf,
     +               w(1,3), w(1,1), w(1,2))
c
        do 220 i = 1, n
          p(i,1) = p(i,1) + w(i,3) *
     +             (hzeta(i) + (-zeta(i)+one) * dzeta(i))
          p(i,2) = p(i,2) + w(i,3) *
     +             (hzeta(i) + (-zeta(i)-one) * dzeta(i))
  220   continue
c
c stoll correction
c
        if (stoll) then
          do 250 i = 1, n
            rhosxt(i) = rho(i,1)**sixth
            rhotrd(i) = rhosxt(i) * rhosxt(i)
  250     continue
c
          call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,
     +                 w(1,3), w(1,1), w(1,2))
c
          do 260 i = 1, n
            p(i,1) = p(i,1) - w(i,3)
  260     continue
c
          call xcvwnd (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,  w(1,3))
c
          do 270 i = 1, n
            p(i,1) = p(i,1) + w(i,3)
  270     continue
c
          do 280 i = 1, n
            rhosxt(i) = rho(i,2)**sixth
            rhotrd(i) = rhosxt(i) * rhosxt(i)
  280     continue
c
          call xcvwnf (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,
     +                 w(1,3), w(1,1), w(1,2))
c
          do 290 i = 1, n
            p(i,2) = p(i,2) - w(i,3)
  290     continue
c
          call xcvwnd (n,      rhosxt, rhotrd, caf, cbf, ccf,
     +                 x0f,    bf,     cf,     qf,  w(1,3))
c
          do 300 i = 1, n
            p(i,2) = p(i,2) + w(i,3)
  300     continue
        endif
c
      endif
c
c
      return
      end
