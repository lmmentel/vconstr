      subroutine pw91x  (lnlpot, nspin, npdim, np, dens, grdrho, d2rho,
     +                   nilrho, enlx,  xcpot)
c
c ======================================================================
c
c purpose: calculate Perdew's 91 exchange correction
c          J.P. Perdew
c
c input  : nspin -  =1: spin restricted
c          npdim -
c          np - nr. of points
c          dens(npdim, nspin)- density
c          grdrho() - gradient of density in cartesian coordinates
c
c output : enlx(npdim) - correction to exchange energy in all points
c          xcpot(npdim,nspin) - non-local correction is added
c
c scratch: rhotmp, s, u, v, aaa, bbb
c
c ======================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical lnlpot, nilrho
c
      parameter (zero = 0d0,  one = 1d0)
      parameter (third = 1d0/3d0)
      parameter (frthrd = 4d0/3d0)
c
      dimension dens(npdim,nspin)
      dimension grdrho(npdim,3,nspin), d2rho(npdim,6,nspin)
      dimension enlx(npdim), xcpot(npdim,nspin)
      dimension rhotmp(np), s(np), u(np), v(np), aaa(np), bbb(np),
     +          nilrho(np,nspin)
c
c ----------------------------------------------------------------------
c
      iij(i,j)= ( ((i-1)*(6-i))/2 + j )
c
c ======================================================================
c
      pi = acos(-one)
      ckf = (3d0*pi*pi)**third
c
      enlx(1:np) = 0.d0
c
      do 1120 ispin = 1, nspin
c
        do 1110 i = 1, np
          if (nilrho(i,ispin)) goto 1110
c
          grd = sqrt (grdrho(i,1,ispin)**2+
     +                grdrho(i,2,ispin)**2+
     +                grdrho(i,3,ispin)**2)
          rhotmp(i) = dens(i,ispin)*nspin
          grdtmp = grd * nspin
c
          if (lnlpot) then
             rlapla = (d2rho(i,1,ispin) + d2rho(i,4,ispin) +
     +                 d2rho(i,6,ispin) )*nspin
          endif
c
          rho43 = rhotmp(i) ** frthrd
          rkf = ckf * rhotmp(i) ** third
          s(i) = grdtmp / (rho43 * 2d0 * ckf)
c
          if (lnlpot) then
            tmp = zero
            do 100 k = 1, 3
              do 50 l = 1, 3
                if (k.gt.l) then
                  iii = iij(l,k)
                else
                  iii = iij(k,l)
                endif
                tmp = tmp + grdrho(i,k,ispin)*grdrho(i,l,ispin)*
     +                d2rho(i,iii,ispin)
 50           continue
 100        continue
            tmp = tmp*nspin**3/grdtmp
            u(i) = tmp/(rhotmp(i)**2*(2*rkf)**3)
            v(i) = rlapla/(rhotmp(i)*(2*rkf)**2)
          endif
c
 1110   continue
c
c
        call pw91x1 (np, lnlpot, rhotmp, s, u, v, nilrho(1,ispin),
     +               aaa, bbb)
c
        do 1115 i = 1, np
          enlx(i) = enlx(i) + rhotmp(i)*aaa(i)/nspin
 1115   continue
c
        if (lnlpot) then
          do 1116 i = 1, np
            xcpot(i, ispin) = xcpot(i,ispin) + bbb(i)
 1116     continue
        endif
c
 1120 continue
c
c
      return
      end
