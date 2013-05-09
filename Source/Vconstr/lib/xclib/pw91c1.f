      subroutine pw91c1(np,lnlpot,rs,zet,t,uu,vv,ww,nilrho,ec,h,
     + dvcup,dvcdn)
c
c ======================================================================
c
c purpose: gga91 correlation energy and potential
c
c input  : np
c          lnlpot - .true.: calculate potential
c          rs     - seitz radius
c          zet    - relative spin polarization
c          t      - abs(grad d)/(d*2.*ks*g)
c          uu     - (grad d)*grad(abs(grad d))/(d**2 * (2*ks*g)**3)
c          vv     - (laplacian d)/(d * (2*ks*g)**2)
c          ww     -  (grad d)*(grad zet)/(d * (2*ks*g)**2
c          nilrho - logical array determining wether non-local
c                   correction can and for numerical
c                   reasons should be neglected
c
c output : h(np) - nonlocal part of correlation energy per electron
c          dvcup(np), dvcdn(np) - potential
c
c scratch: work  - (np,6)
c          ec    - (np,3): ec(*,1) lsda correlation per electron
c                             ec(*,2) derivative w.r.t. rs
c                             ec(*,3) derivative w.r.t. zet
c
c remark * uu, vv, ww are input, and dvcup, dvcdn are output only when
c          lnlpot = .true.
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical    lnlpot,nilrho
c
      dimension  rs    (np),      zet    (np),   t  (np),
     +           uu    (np),      vv     (np),   ww (np),
     +           h     (np),      nilrho (np),
     +           dvcup (np),      dvcdn  (np),   ec (np,3)
      dimension  work  (np,6)
c
      parameter (zero  = 0d0,  one = 1d0,    two = 2d0,
     +           three = 3d0,  eps = 1d-20,  tol = 1d-5)
c
      data xnu,cc0,cx,alf/15.75592d0,0.004235d0,-0.001667212d0,0.09d0/
      data c1,c2,c3,c4/0.002568d0,0.023266d0,7.389d-6,8.723d0/
      data c5,c6,a4/0.472d0,7.389d-2,100.d0/
      data thrdm,thrd2/-0.333333333333d0,0.666666666667d0/
c
c ======================================================================
c-pph
      pi = acos(-one)
      call corlsd (np, lnlpot, rs, zet, nilrho,
     +             work(1,1), work(1,2), work(1,3), work(1,4),
     +             work(1,5), work(1,6), ec(1,1), ec(1,2), ec(1,3))
c
      do 1000 ip = 1, np
c
        if (nilrho(ip)) then
          h(ip) = zero
          if (lnlpot) then
            dvcup(ip) = zero
            dvcdn(ip) = zero
          endif
          goto 1000
        endif
c
CCC     if (abs(zet(ip)).gt.0.99d0) zet(ip) = zet(ip)*0.99d0
c
        fk = 1.91915829d0/rs(ip)
        fk = max(eps,fk)
        sk = sqrt(4.d0*fk/pi)
        tmp1 = zero
        if (zet(ip) .gt.-one) tmp1 = (one+zet(ip)) ** thrd2
        tmp2 = zero
        if (zet(ip) .lt. one) tmp2 = (one-zet(ip)) ** thrd2
        g = (tmp1+tmp2) / two
c
        bet = xnu*cc0
        delt = two * alf / bet
        g3 = g**3
        g4 = g3*g
        pon = -delt*ec(ip,1)/(g3*bet)
c
c GtV 96.03.05: numerical precision of exponential
CCC     b = delt / (exp(pon) - one)
        b = delt
        if (b.ne.zero) then
          if (abs(pon) .gt. tol) then
            tmp = exp(pon) - one
          else
            tmp = pon + (pon**2 / two) + (pon**3 / 6d0)
          endif
          b = b / tmp
        endif
c GtV 96.03.05: end
c
        b2 = b*b
        t2 = t(ip) **2
        t4 = t2 **2
        t6 = t4*t2
        rs2 = rs(ip) **2
        rs3 = rs2 * rs(ip)
        q4 = one + b*t2
        q5 = one + b*t2 + b2*t4
        q6 = c1 + c2*rs(ip) + c3*rs2
        q7 = one + c4*rs(ip) + c5*rs2 + c6*rs3
        cc = -cx + q6/q7
        r0 = (sk/fk) **2
        r1 = a4 * r0 * g4
        coeff = cc - cc0 - 3d0 * cx / 7d0
        r2 = xnu * coeff * g3
        r3 = exp (-r1*t2)
c
c GtV 96.03.05: numerical precision of logarithm
CCC     h0 = g3 * (bet/delt) * log (one + delt*q4*t2/q5)
        arg = delt * q4 * t2
        if (arg.ne.zero) arg = arg / q5
        if (abs(arg) .gt. tol) then
          tmp = log (one + arg)
        else
          tmp = arg - (arg**2 / two) + (arg**3 / three)
        endif
        h0 = g3 * (bet/delt) * tmp
c GtV: end
c
        h1 = r3 * r2 * t2
        h(ip) = h0 + h1
c
        if (lnlpot) then
c
c local correlation option:
c h = 0.0d0
c energy done. now the potential:
c
          ccrs = (c2+ two*c3*rs(ip))/q7 -
     +           q6*(c4+ two*c5*rs(ip)+3d0*c6*rs2)/q7**2
          rsthrd = rs(ip)/3d0
          r4 = rsthrd*ccrs/coeff
c
          tmp1 = zero
          if (zet(ip) .gt.-one) tmp1 = (one+zet(ip)) ** thrdm
          tmp2 = zero
          if (zet(ip) .lt. one) tmp2 = (one-zet(ip)) ** thrdm
          gz = (tmp1-tmp2) / 3d0
c
c -pph
c
          b2fac = delt*b + b2
          bg = -3d0 * ec(ip,1) * b2fac / (bet*g4)
          bec = b2fac / (bet*g3)
          q8 = q5**2 + delt*q4*q5*t2
          q9 = one + two*b*t2
          h0b = -bet * g3 * b * t6 * (two + b*t2) / q8
          h0rs = -rsthrd * h0b * bec * ec(ip,2)
          fact0 =  two*delt-6d0*b
          fact1 = q5*q9 + q4*q9**2
          h0bt =  two*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
          h0rst = rsthrd*t2*h0bt*bec*ec(ip,2)
          h0z = 3d0*gz*h0/g + h0b*(bg*gz+bec*ec(ip,3))
          h0t =  two*bet*g3*q9/q8
          h0zt = 3d0*gz*h0t/g+h0bt*(bg*gz+bec*ec(ip,3))
          fact2 = q4*q5+b*t2*(q4*q9+q5)
          fact3 =  two*b*q5*q9+delt*fact2
          h0tt = 4d0*bet*g3*t(ip)*(2.d0*b/q8-(q9*fact3/q8)/q8)
          h1rs = r3*r2*t2*(-r4+r1*t2/3d0)
          fact4 =  two-r1*t2
          h1rst = r3*r2*t2*(2.d0*r4*(one -r1*t2)-thrd2*r1*t2*fact4)
          h1z = gz*r3*r2*t2*(3.d0-4.d0*r1*t2)/g
          h1t =  two*r3*r2*(one -r1*t2)
          h1zt =  two*gz*r3*r2*(3.d0-11d0*r1*t2+4d0*r1*r1*t4)/g
          h1tt = 4d0*r3*r2*r1*t(ip)*(-two+r1*t2)
          hrs = h0rs+h1rs
          hrst = h0rst+h1rst
          ht = h0t+h1t
          htt = h0tt+h1tt
          hz = h0z+h1z
          hzt = h0zt+h1zt
          comm = h(ip)+hrs+hrst+t2*ht/6d0 + 7d0*t2*t(ip)*htt/6d0
          pref = hz-gz*t2*ht/g
          fact5 = gz*(two*ht + t(ip)*htt) / g
          comm = comm - pref*zet(ip) -
     +              uu(ip)*htt - vv(ip)*ht - ww(ip)*(hzt-fact5)
          dvcup(ip) = comm + pref
          dvcdn(ip) = comm - pref
        endif
 1000 continue
c
c
      return
      end
