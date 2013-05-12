      subroutine xcvwnf (n,      rhosxt, rhotrd, ca, cb, cc,
     +                   x0,     b,      c,      q,
     +                   result, w1,     w2)
c
c ======================================================================
c
c purpose: evaluation of a function, for a vector of arguments, occur-
c          ring in the expression for the exchange correlation poten-
c          tial, due to vosko, wilk and nusair (can. j. phys. 58
c          (1980) 1200).
c
c input  : n      - vector length
c          rhosxt - rho**(1/6),  rho is the charge density
c          rhotrd - rho**(1/3)
c          ca..q  - parameters occurring in the function
c
c output : result - vector of function values
c
c scratch: w1,w2  - (n)
c
c remark * the function is given by equation (4.4) in the paper by vwn
c          for computational purpose it is rewritten as
c
c          f = log(g) - 2*cb*log (sqrt(g)-x0*rhosxt) + (cb-1)*
c
c                 log (g+b*sqrt(g)*rhosxt+c*rhotrd) +
c
c                 cc * atan ( (q*rhosxt)/(2*sqrt(g)+b*rhosxt) )
c
c          g = (3/(4*pi)) ** (1/3)
c
c          in this form only two logarithm evaluations are needed, and
c          possible numerical problems for rho --> 0.0 are avoided.
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
c gamma = (3/(4*pi)) ** (1/3)
c
      parameter (alogga = -0.4774706527671d0,
     +           gamma  =  0.6203504908994d0,
     +           sqrtga =  0.7876233178997d0,
     +           one    =  1.0d0,
     +           two    =  2.0d0)
c
      dimension  result(n), rhotrd(n), rhosxt(n),
     +           w1(n),     w2(n)
c
c ======================================================================
c
      do 10 k = 1, n
        w1(k) = sqrtga - x0 * rhosxt(k)
        result(k) = log (w1(k))
   10 continue
c
      do 20 k = 1, n
        result(k) = alogga - two * cb * result(k)
        w1(k) = gamma + b * sqrtga * rhosxt(k) + c * rhotrd(k)
        w2(k) = log (w1(k))
   20 continue
c
      do 30 k = 1, n
        result(k) = result(k) + (cb - one) * w2(k)
        w1(k) = q * rhosxt(k) / ( two * sqrtga + b * rhosxt(k) )
   30 continue
c
      do 40 k = 1, n
        result(k) = ca * ( result(k) + cc * atan (w1(k)) )
   40 continue
c
c
      return
      end
