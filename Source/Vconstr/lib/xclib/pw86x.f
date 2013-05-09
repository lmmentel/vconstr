      subroutine pw86x (lnlpot, nspin, npdim, np, dens, grdrho, d2rho,
     +                  nilrho, enlx,  xcpot)
c
c ======================================================================
c
c purpose: calculate perdew''s ''86 exchange correction
c
c          Perdew and Yue, Phys. Rev. B(33), 1986, page 8800
c
c input  : lnlpot - true: calculate also pot.
c                   false: energy density only
c          nspin -  =1: spin restricted
c          npdim - max. nr. of points
c          np - nr. of points
c          dens(npdim, nspin)- density
c          grdrho() - gradient of density in cartesian coordinates
c          d2rho() - second order derivatives (only necessary for pot.)
c          nilrho() - as usual
c
c output:  enlx(npdim) - nonlocal correction to exchange energy added
c                        in all points
c          xcpot(npdim,nspin) - non-local correction to pot. is added
c                               in case lnlpot = .true.
c
c *=====================================================================
c
      implicit   integer (i-n)
      implicit   double precision (a-h,o-z)
c
      logical lnlpot, nilrho
c
      parameter (three = 3.0d0, four = 4.0d0,
     +           ssmall = 1.0d-20)
      parameter (third =1.0d0/3.0d0)
c
      dimension dens(npdim,nspin)
      dimension grdrho(npdim,3,nspin), d2rho(npdim,6,nspin)
      dimension enlx(npdim), xcpot(npdim,nspin)
      dimension nilrho(np, nspin)
c
      dimension grdcsi(3)
c
      iij(i,j) = ( ((i-1)*(6-i))/2 + j )
c
c ======================================================================
c
      ppp = .0864d0
      pm = 1.d0/15
      pb = 14.d0
      pc = .2d0
      pi = acos(-1.d0)
      cx = three/four*(three/pi)**(third)
      d1 = 4*(3*pi**2)**(2*third)
      ckf = (3*pi*pi)**third
c
      do 1120 ispin = 1, nspin
        do 1110 ip = 1, np
          if (nilrho(ip,ispin)) goto 1110
          grd=sqrt(grdrho(ip,1,ispin)**2+
     +             grdrho(ip,2,ispin)**2+
     +             grdrho(ip,3,ispin)**2)
          rhotmp = dens(ip,ispin)
          rho43 = rhotmp**(4*third)
          rho53 = rhotmp**(5*third)
          rho13 = rhotmp**(third)
c
c s = 0 leads to problems so use very small s instead
c
          ss = max(grd/(rho43*2*ckf)*nspin**(-third),ssmall)
c
c fff = nonlocal part of enhancement factor in its spin restricted form
c
          fff = (1.d0+ppp*ss**2/pm+pb*ss**4+pc*ss**6 )**pm - 1.d0
c
c energy density
c
          enlx(ip) = enlx(ip) - cx*nspin**(4*third)*rho43*fff/nspin
c
          if (lnlpot) then
c
c potential
c
            rlapla = (d2rho(ip,1,ispin) + d2rho(ip,4,ispin) +
     +               d2rho(ip,6,ispin) )
            sum = 0.d0
            do 10 k = 1, 3
              sum = sum + grdrho(ip,k,ispin)**2
 10         continue
            eta = rlapla/(d1*rho53)*nspin**(-2*third)
            sum = 0.d0
            do 40 k = 1, 3
              grdcsi(k) = 0.d0
              do 30 l = 1, 3
                if (k.gt.l) then
                  iii = iij(l,k)
                else
                  iii = iij(k,l)
                endif
                if (grd.gt.ssmall) then
c
c grad rho = 0 ignored
c
                  grdcsi(k) = grdcsi(k) +
     +                       grdrho(ip,l,ispin)*d2rho(ip,iii,ispin)/grd
                endif
 30           continue
 40         continue
            sum = 0.d0
            do 45 k = 1, 3
              sum = sum + grdrho(ip,k,ispin)*grdcsi(k)
 45         continue
            uu = sum/(24*pi**2*rhotmp**3*nspin)
c
c derivatives of enhancement factor
c
            dfff = pm*(2*ppp*ss/pm + 4*pb*ss**3 + 6*pc*ss**5)*
     -             (1.d0+ppp*ss**2/pm+pb*ss**4+pc*ss**6)**(-1.d0+pm)
            ddfff =
     +      (-1.d0 + pm)*pm*(2*ppp*ss/pm + 4*pb*ss**3 + 6*pc*ss**5)**2*
     -      (1.d0 + ppp*ss**2/pm + pb*ss**4 + pc*ss**6)**(-2 + pm) +
     -      pm*(2*ppp/pm + 12*pb*ss**2 + 30*pc*ss**4)*
     -      (1.d0 + ppp*ss**2/pm + pb*ss**4 + pc*ss**6)**(-1.d0 + pm)
c
            vxlda = -(three/pi)**third*nspin**third*rho13
                    xcpot(ip,ispin) = xcpot(ip,ispin) +
     +              vxlda*( fff - three/four*eta/ss*dfff -
     +              three/four*(uu-four/three*ss**(three))*
     +              (-dfff/ss**2 + ddfff/ss))
          endif
c
 1110   continue
c
 1120 continue
c
c
      return
c ======================================================================
      end
