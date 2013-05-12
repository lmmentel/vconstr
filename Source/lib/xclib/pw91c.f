      subroutine pw91c  (lnlpot, nspin, npdim, np, dens, grdrho, d2rho,
     +                   nilrho, eprdew, xcpot)
c
c ======================================================================
c
c purpose: calculate Perdew's correction
c
c input  : nspin -  =1: spin restricted
c          npdim -
c          np - nr. of points
c          dens(npdim, nspin)- density
c          grdrho(npdim,3,nspin) - gradient of the density
c          d2rho(npdim,6,nspin) - second derivatives
c          nilrho - (logical)
c
c output : eprdew(npdim) - correction to exchange energy in all points
c          xcpot(npdim,npsin) -
c
c scratch: rs(np), zeta(np), t(np), uu(np), vv(np), ww(np), aaa(np),
c          dv(np,2), ec(np,3)
c
c
c ======================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical lnlpot, nilrho
c
      parameter(half =1.0d0/2.0d0, sixth = 1.0d0/6.0d0)
      parameter(twthrd = 2.0d0/3.0d0)
      parameter(third =1.0d0/3.0d0)
      parameter(sevsix=7.0d0/6.0d0)
      parameter(eps = 1.d-70,  zetax = 0.9999d0)
c
      dimension dens(npdim,nspin)
      dimension grdrho(npdim,3,nspin)
      dimension d2rho(npdim,6,nspin)
      dimension eprdew(npdim), xcpot(npdim,nspin),
     +          dr(3), d2r(6)
      dimension rs(np), zeta(np), t(np), uu(np), vv(np), ww(np),
     +          aaa(np), dv(np, 2), nilrho(np,3),
     +          ec(np,3)
c
c ----------------------------------------------------------------------
c
      iij(i,j)= ( ((i-1)*(6-i))/2 + j )
c
c ======================================================================
c
      pi=acos(-1.d0)
      fourpi=4.d0*pi
c
      ckf = (3 * pi**2) ** third
      cks = sqrt (4*ckf/pi)
      eprdew(1:np) = 0.d0
c
      do 1110 i = 1, np
        if (nilrho(i,3)) goto 1110
c
        do 30 j = 1, 3
          if (nilrho(i,1)) then
            dr(j) = grdrho(i,j,nspin)
          elseif (nilrho(i,nspin)) then
            dr(j) = grdrho(i,j,1)
          else
            dr(j) = grdrho(i,j,1)
            if (nspin.eq.2) dr(j) = dr(j) + grdrho(i,j,nspin)
          endif
 30     continue
c
        grd = max (eps, sqrt(dr(1)**2+dr(2)**2+dr(3)**2))
c
        if (lnlpot) then
c
          do 50 j = 1, 6
            if (nilrho(i,1)) then
              d2r(j) = d2rho(i,j,nspin)
            elseif (nilrho(i,nspin)) then
              d2r(j) = d2rho(i,j,1)
            else
              d2r(j) = d2rho(i,j,1)
              if (nspin.eq.2) d2r(j) = d2r(j) + d2rho(i,j,nspin)
            endif
 50       continue
c
          rlapla = d2r(1) + d2r(4) + d2r(6)
c
        endif
c
c
        if (nspin.eq.1) then
          zeta(i) = 0.d0
          rho = dens(i,1)
        else
c
          if (nilrho(i,1)) then
            zeta(i) = -1.d0
            rho = dens(i,nspin)
          elseif (nilrho(i,nspin)) then
            zeta(i) = 1.d0
            rho = dens(i,1)
          else
            zeta(i)=(dens(i,1)-dens(i,2))/
     +           (dens(i,1)+dens(i,2))
            rho = dens(i,1) + dens(i,2)
          endif
          if (zeta(i) .gt. zetax) zeta(i) = zetax
          if (zeta(i) .lt.-zetax) zeta(i) = -zetax
c
        endif
c
        tmp1 = max (0.d0, 1.d0+zeta(i))
        tmp2 = max (0.d0, 1.d0-zeta(i))
        g = half * (tmp1**twthrd + tmp2**twthrd)
        t(i) = grd/(2 * g * cks * rho ** sevsix)
        rs(i)=(3.d0/(fourpi*rho))**third
c
        if (lnlpot) then
          rks = max (eps, 2*(3*rho/pi)**sixth)
          tmp = 0.d0
          do 100 k = 1, 3
            do 90 l = 1, 3
              if (k.gt.l) then
                iii = iij(l,k)
              else
                iii = iij(k,l)
              endif
              tmp = tmp + dr(k)*dr(l)*d2r(iii)
 90         continue
 100      continue
          tmp = tmp/grd
          uu(i) = tmp/(rho**2*(2*rks*g)**3)
          vv(i) = rlapla/(rho*(2*rks*g)**2)
c
          ww(i) = 0.d0
          if (nspin.eq.2) then
            tmp = 0.d0
            do 200 j = 1, 3
              if (nilrho(i,1)) then
                tmp = tmp + dr(j)*(-grdrho(i,j,2))/rho -
     +                zeta(i)*dr(j)**2/rho
              elseif (nilrho(i,nspin)) then
                tmp = tmp + dr(j)*grdrho(i,j,1)/rho -
     +                zeta(i)*dr(j)**2/rho
              else
                tmp = tmp + dr(j)*(grdrho(i,j,1)-grdrho(i,j,2))/rho -
     +                zeta(i)*dr(j)**2/rho
              endif
 200        continue
            ww(i) = tmp/(rho*(2*rks*g)**2)
          endif
        endif
c
 1110 continue
c
      call pw91c1 (np, lnlpot, rs, zeta, t, uu, vv, ww,
     +             nilrho(1,3), ec, aaa, dv(1,1), dv(1,2))
c
      do 2000 i = 1, np
        if (nspin.eq.1) then
          rho = dens(i,1)
        else
          if (nilrho(i,1)) then
            rho = dens(i,2)
          elseif (nilrho(i,2)) then
            rho = dens(i,1)
          else
            rho = dens(i,1) + dens(i,2)
          endif
        endif
c
        eprdew(i)= rho*aaa(i)
c
        if (lnlpot) then
          do 300 ispin = 1, nspin
            xcpot(i, ispin) = xcpot(i,ispin) + dv(i,ispin)
 300      continue
        endif
c
 2000 continue
c
c
      return
      end
