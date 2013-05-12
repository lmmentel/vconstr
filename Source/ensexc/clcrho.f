      subroutine clcrho(rho,drho,d2rho,nspin,norb,npnt,npntmx,pmo,
     + vmopao,grid,weight,atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical atmol4
      dimension rho(npntmx,nspin),drho(npntmx,3,nspin),
     + d2rho(npntmx,6,nspin)
      dimension pmo(norb*(norb+1)/2,nspin),vmopao(norb*norb)
      dimension grid(npntmx,3),weight(npnt)
      logical ltrian
      dimension valao(norb),gradx(norb),grady(norb),
     +  gradz(norb),gradxx(norb),gradxy(norb),gradxz(norb),
     +  gradyy(norb),gradyz(norb),gradzz(norb)
      dimension valmo(norb),gxmo(norb),gymo(norb),gzmo(norb),
     + gxxmo(norb),gxymo(norb),gxzmo(norb),gyymo(norb),gyzmo(norb),
     + gzzmo(norb)
      dimension check(10,nspin)
c
      ltrian = .false.

      check(1:10,1:nspin) = 0.d0
      do 20 ipnt=1,npnt
        m=(ipnt-1)*norb+1
        x=grid(ipnt,1)
        y=grid(ipnt,2)
        z=grid(ipnt,3)
        w=weight(ipnt)
        if (atmol4) then
          call aovlnw(x,y,z,valao,gradx,grady,gradz,gradxx,gradxy,
     + gradxz,gradyy,gradyz,gradzz)
        else
          call aovlnv(x,y,z,valao,gradx,grady,gradz,gradxx,gradxy,
     + gradxz,gradyy,gradyz,gradzz)
        endif
        call vecmat(valao,vmopao,norb,valmo)
        call vecmat(gradx,vmopao,norb,gxmo)
        call vecmat(grady,vmopao,norb,gymo)
        call vecmat(gradz,vmopao,norb,gzmo)
        call vecmat(gradxx,vmopao,norb,gxxmo)
        call vecmat(gradxy,vmopao,norb,gxymo)
        call vecmat(gradxz,vmopao,norb,gxzmo)
        call vecmat(gradyy,vmopao,norb,gyymo)
        call vecmat(gradyz,vmopao,norb,gyzmo)
        call vecmat(gradzz,vmopao,norb,gzzmo)
        do 10 is=1,nspin
          rho(ipnt,is)=valmat(pmo(1,is),norb,valmo,ltrian)
          call valgrd(dxao,dyao,dzao,pmo(1,is),norb,valmo,
     + gxmo,gymo,gzmo)
          drho(ipnt,1,is)=dxao
          drho(ipnt,2,is)=dyao
          drho(ipnt,3,is)=dzao
          call val2grd(dxx,dxy,dxz,dyy,dyz,dzz,pmo(1,is),norb,
     + valmo,gxmo,gymo,gzmo,gxxmo,gxymo,gxzmo,gyymo,gyzmo,gzzmo)
          d2rho(ipnt,1,is)=dxx
          d2rho(ipnt,2,is)=dxy
          d2rho(ipnt,3,is)=dxz
          d2rho(ipnt,4,is)=dyy
          d2rho(ipnt,5,is)=dyz
          d2rho(ipnt,6,is)=dzz
          check(1,is)=check(1,is)+w*rho(ipnt,is)
          check(2,is)=check(2,is)-w*x*drho(ipnt,1,is)
          check(3,is)=check(3,is)-w*y*drho(ipnt,2,is)
          check(4,is)=check(4,is)-w*z*drho(ipnt,3,is)
          check(5,is)=check(5,is)+w*x**2*d2rho(ipnt,1,is)/2
          check(6,is)=check(6,is)+w*x*y*d2rho(ipnt,2,is)
          check(7,is)=check(7,is)+w*x*z*d2rho(ipnt,3,is)
          check(8,is)=check(8,is)+w*y**2*d2rho(ipnt,4,is)/2
          check(9,is)=check(9,is)+w*y*z*d2rho(ipnt,5,is)
          check(10,is)=check(10,is)+w*z**2*d2rho(ipnt,6,is)/2
   10   continue
   20 continue
      close(11)
      if (nspin.eq.2) then
        nel=nint(check(1,1))
        check(1:10,1)=check(1:10,1)-nel
        write(6,'(/'' alpa  '',10(2x,g8.2))')(check(j,1),j=1,10)
        nel=nint(check(1,2))
        check(1:10,2)=check(1:10,2)-nel
        write(6,'('' beta  '',10(2x,g8.2))/')(check(j,2),j=1,10)
      else
        nel=nint(check(1,1))
        check(1:10,1)=check(1:10,1)-nel
        write(6,'(/''  nel   '',10(2x,g8.2))/')(check(j,1),j=1,10)
      endif
c 
      return
      end
