      subroutine xcvwnd (n,      rhosxt, rhotrd, ca, cb, cc,
     +                   x0,     b,      c,      q,
     +                   result)
c
c ======================================================================
c
c purpose: evaluation of a function, for a vector of arguments, occur-
c          ring in the expression for the exchange correlation poten-
c          tial, due to vosko, wilk and nusair (can. j. phys. 58
c          (1980) 1200).
c          this function is the derivative of the function, calculated
c          in xcvwnf, with respect to the charge density.
c
c input  : n      - vector length
c          rhosxt - rho**(1/6),  rho is the charge density
c          rhotrd - rho**(1/3)
c          ca..q  - parameters occurring in the function
c
c output : result - vector of function values
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
c gamma = (3/(4*pi)) ** (1/3)
c
      parameter (gamma  = 0.6203504908994d0,
     +           sqrtga = 0.7876233178997d0,
     +           half   = 0.5d0,  third = 1.d0/3d0)
c
      dimension  rhosxt(n), rhotrd(n), result(n)
c
c ======================================================================
c
      do 10 i = 1, n
        result(i) = cc * q * sqrtga * rhosxt(i) /
     +             ( q*q*rhotrd(i) + ( 2.d0*sqrtga + b*rhosxt(i) )**2 )
   10 continue
c
      do 20 i = 1, n
        result(i) = 1.d0 - result(i) - cb * sqrtga /
     +              (sqrtga - x0 * rhosxt(i))
   20 continue
c
      do 30 i   = 1, n
        result(i) = result(i) + (cb-1.d0) *
     +              (gamma + half * b * sqrtga * rhosxt(i)) /
     +              (gamma + b * sqrtga * rhosxt(i) + c * rhotrd(i))
   30 continue
c
      result(1:n) = third*ca * result(1:n)
c
      return
      end
