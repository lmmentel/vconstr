      subroutine clcrho(rho,drho,d2rho,norb,npnt,npntmx,nspin,rnel,
     + pmo,vmopao,grid,weight,atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical atmol4,ltrian
      dimension rho(npntmx,nspin),drho(npntmx,3,nspin),
     + d2rho(npntmx,6,nspin),pmo(norb*(norb+1)/2,nspin),
     + vmopao(norb*norb)
      dimension grid(npntmx,3),weight(npnt)
      dimension valao(norb),gradx(norb),grady(norb),
     +  gradz(norb),gradxx(norb),gradxy(norb),gradxz(norb),
     +  gradyy(norb),gradyz(norb),gradzz(norb)
      dimension valmo(norb),gxmo(norb),gymo(norb),gzmo(norb),
     + gxxmo(norb),gxymo(norb),gxzmo(norb),gyymo(norb),
     + gyzmo(norb),gzzmo(norb)
c
      ltrian = .false.
c
      chkval=-rnel
      chkgx=-rnel
      chkgy=-rnel
      chkgz=-rnel
      chkgxx=-2*rnel
      chkgxy=-rnel
      chkgxz=-rnel
      chkgyy=-2*rnel
      chkgyz=-rnel
      chkgzz=-2*rnel
      do ip=1,npnt
        x=grid(ip,1)
        y=grid(ip,2)
        z=grid(ip,3)
        w=weight(ip)
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
        do is=1,nspin
          rho(ip,is)=valmat(pmo(1,is),norb,valmo,ltrian)
          call valgrd(dxao,dyao,dzao,pmo(1,is),norb,valmo,
     + gxmo,gymo,gzmo)
          drho(ip,1,is)=dxao
          drho(ip,2,is)=dyao
          drho(ip,3,is)=dzao
          call val2grd(dxx,dxy,dxz,dyy,dyz,dzz,pmo(1,is),norb,
     + valmo,gxmo,gymo,gzmo,gxxmo,gxymo,gxzmo,gyymo,gyzmo,gzzmo)
          d2rho(ip,1,is)=dxx
          d2rho(ip,2,is)=dxy
          d2rho(ip,3,is)=dxz
          d2rho(ip,4,is)=dyy
          d2rho(ip,5,is)=dyz
          d2rho(ip,6,is)=dzz
          chkval=chkval+w*rho(ip,is)
          chkgx=chkgx-w*x*drho(ip,1,is)
          chkgy=chkgy-w*y*drho(ip,2,is)
          chkgz=chkgz-w*z*drho(ip,3,is)
          chkgxx=chkgxx+w*x**2*d2rho(ip,1,is)
          chkgxy=chkgxy+w*x*y*d2rho(ip,2,is)
          chkgxz=chkgxz+w*x*z*d2rho(ip,3,is)
          chkgyy=chkgyy+w*y**2*d2rho(ip,4,is)
          chkgyz=chkgyz+w*y*z*d2rho(ip,5,is)
          chkgzz=chkgzz+w*z**2*d2rho(ip,6,is)
        enddo
      enddo
      close(11)
      write(6,'(/,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2)') 
     + chkval,chkgx,chkgy,chkgz,chkgxx,chkgxy,chkgxz,chkgyy,chkgyz,
     + chkgzz
c
      return
      end
