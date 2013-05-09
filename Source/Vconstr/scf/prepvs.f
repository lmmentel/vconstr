      subroutine prepvs(vxc,rho,drho,nspin,npnt)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(third  = 1.d0/3.d0, fthird = 4.d0/3.d0)
      parameter(xconst = -0.9305257363491d0,
     +          restrc =  0.7937005259841d0)
      parameter(b = 0.0042d0, b6 = 6d0*b)
c
      dimension vxc(npnt,nspin)
      dimension rho(npnt,nspin),drho(npnt,3,nspin)
c
      vxc(1:npnt,1:nspin)=0.d0
c
      if (nspin.eq.1) then
        vxc(1:npnt,1) = xconst*restrc*rho(1:npnt,1)**third
      else
        vxc(1:npnt,1) = xconst*(rho(1:npnt,1)**third+
     + rho(1:npnt,2)**third)
      endif
c
      factor = -b
      if (nspin.eq.1) then
        factor = 2.d0*factor
        xfactr = restrc
      endif
c
      do is=1,nspin
c
        do ipnt=1,npnt
c
          vxc(ipnt,is) = 2.d0*xfactr*xconst*rho(ipnt,is)**third
          rho13=rho(ipnt,is)**third
          rho43=rho(ipnt,is)**fthird
          grad=sqrt(drho(ipnt,1,is)**2+drho(ipnt,2,is)**2+
     + drho(ipnt,3,is)**2)
          x=grad/rho43
          xsq=x**2
          sqrt1x=sqrt(1.d0+xsq)
          sinhix=log(x+sqrt1x)
          f=1.d0/(1.d0+b6*x*sinhix)
          vB=2.d0*factor*rho13*xsq*f
c
          vxc(ipnt,is) = vxc(ipnt,is)+vB
c
        enddo
      enddo
c
      return
      end
