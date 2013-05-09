      subroutine clcrho(rho,drho,d2rho,norb,npnt,npntmx,rnel,pmo,
     + vmopao,grid,weight,atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical atmol4,ltrian
      dimension rho(npntmx,1),drho(npntmx,3,1),d2rho(npntmx,6,1),
     + pmo(norb*(norb+1)/2),vmopao(norb*norb)
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
      do 20 ip=1,npnt
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
        rho(ip,1)=valmat(pmo,norb,valmo,ltrian)
        call valgrd(dxao,dyao,dzao,pmo,norb,valmo,gxmo,gymo,gzmo)
        drho(ip,1,1)=dxao
        drho(ip,2,1)=dyao
        drho(ip,3,1)=dzao
        call val2grd(dxx,dxy,dxz,dyy,dyz,dzz,pmo,norb,valmo,
     + gxmo,gymo,gzmo,gxxmo,gxymo,gxzmo,gyymo,gyzmo,gzzmo)
        d2rho(ip,1,1)=dxx
        d2rho(ip,2,1)=dxy
        d2rho(ip,3,1)=dxz
        d2rho(ip,4,1)=dyy
        d2rho(ip,5,1)=dyz
        d2rho(ip,6,1)=dzz
        chkval=chkval+w*rho(ip,1)
        chkgx=chkgx-w*x*drho(ip,1,1)
        chkgy=chkgy-w*y*drho(ip,2,1)
        chkgz=chkgz-w*z*drho(ip,3,1)
        chkgxx=chkgxx+w*x**2*d2rho(ip,1,1)
        chkgxy=chkgxy+w*x*y*d2rho(ip,2,1)
        chkgxz=chkgxz+w*x*z*d2rho(ip,3,1)
        chkgyy=chkgyy+w*y**2*d2rho(ip,4,1)
        chkgyz=chkgyz+w*y*z*d2rho(ip,5,1)
        chkgzz=chkgzz+w*z**2*d2rho(ip,6,1)
   20 continue
      close(11)
      write(6,'(/,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2)') 
     + chkval,chkgx,chkgy,chkgz,chkgxx,chkgxy,chkgxz,chkgyy,chkgyz,
     + chkgzz
c
      return
      end
