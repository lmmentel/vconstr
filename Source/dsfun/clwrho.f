      subroutine clwrho(norb,npnt,npntmx,pmo,valmo,grid,weight)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical ltrian
      dimension pmo(norb*(norb+1)/2),valmo(npnt*norb)
      dimension grid(npntmx,3),weight(npnt)
      dimension gxmo(norb),gymo(norb),gzmo(norb),gxxmo(norb),
     + gxymo(norb),gxzmo(norb),gyymo(norb),gyzmo(norb),gzzmo(norb)
c
      ltrian = .false.
c
      open(11,file='hfvls')
      rewind(11)
c
      chkval=0.d0
      chkgx=0.d0
      chkgy=0.d0
      chkgz=0.d0
      chkgxx=0.d0
      chkgxy=0.d0
      chkgxz=0.d0
      chkgyy=0.d0
      chkgyz=0.d0
      chkgzz=0.d0
      do ipnt=1,npnt
        m=(ipnt-1)*norb+1
        x=grid(ipnt,1)
        y=grid(ipnt,2)
        z=grid(ipnt,3)
        w=weight(ipnt)
        read(11,*) (gxmo(i),i=1,norb)
        read(11,*) (gymo(i),i=1,norb)
        read(11,*) (gzmo(i),i=1,norb)
        read(11,*) (gxxmo(i),i=1,norb)
        read(11,*) (gxymo(i),i=1,norb)
        read(11,*) (gxzmo(i),i=1,norb)
        read(11,*) (gyymo(i),i=1,norb)
        read(11,*) (gyzmo(i),i=1,norb)
        read(11,*) (gzzmo(i),i=1,norb)
        rho=valmat(pmo,norb,valmo(m),ltrian)
        call valgrd(dxao,dyao,dzao,pmo,norb,valmo(m),gxmo,gymo,gzmo)
        call val2grd(dxx,dxy,dxz,dyy,dyz,dzz,pmo,norb,valmo(m),
     + gxmo,gymo,gzmo,gxxmo,gxymo,gxzmo,gyymo,gyzmo,gzzmo)
        chkval=chkval+w*rho
        chkgx=chkgx-w*x*dxao
        chkgy=chkgy-w*y*dyao
        chkgz=chkgz-w*z*dzao
        chkgxx=chkgxx+w*x**2*dxx
        chkgxy=chkgxy+w*x*y*dxy
        chkgxz=chkgxz+w*x*z*dxz
        chkgyy=chkgyy+w*y**2*dyy
        chkgyz=chkgyz+w*y*z*dyz
        chkgzz=chkgzz+w*z**2*dzz 
        write(13,*)rho/2,dxao/2,dyao/2,dzao/2,dxx/2,dxy/2,dxz/2,dyy/2,
     + dyz/2,dzz/2
      enddo
      close(11)
      write(6,'(g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2,g8.2)') 
     + chkval,chkgx,chkgy,chkgz,chkgxx/2,chkgxy,chkgxz,chkgyy/2,chkgyz,
     + chkgzz/2
c
      return
      end
