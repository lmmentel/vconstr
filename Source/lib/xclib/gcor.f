      subroutine gcor (lnlpot, np, a, a1, b1, b2, b3, b4, p, rs, nilrho,
     +                 gg, ggrs)
c
c ======================================================================
c
c input  : lnlpot, np, a, a1, b1, b2, b3, b4, p, rs(np), nilrho(np)
c
c output : gg(np), ggrs(np)
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    lnlpot, nilrho
c
      parameter (one = 1d0,   two = 2d0,  three = 3d0,
     +           half = 0.5d0,  tol = 1d-5)
c
      dimension  rs(np), gg(np), ggrs(np), nilrho(np)
c
c ======================================================================
c
      do 100 ip = 1, np
        if (nilrho(ip)) goto 100
c
        q0 = -two * a * (one + a1 * rs(ip))
        rs12 = sqrt (rs(ip))
        rs32 = rs12 * rs(ip)
        rsp = rs(ip) ** p
        q1 = two * a * (b1*rs12 + b2*rs(ip) + b3*rs32 + b4*rs(ip)*rsp)
c
        tmp = one/q1
        if (abs(tmp) .gt. tol) then
          q2 = log (one + tmp)
        else
          q2 = ((tmp/three - half)*tmp + one) * tmp
        endif
c
        gg(ip) = q0 * q2
c
        if (lnlpot) then
          q3 = a * (b1/rs12 + two*b2 + three*b3*rs12 +
     +              two*b4*(p+one)*rsp)
          ggrs(ip) = -two * a * a1 * q2 -
     +               q0 * q3 / (q1*(q1+one))
        endif
c
  100 continue
c
c
      return
      end
