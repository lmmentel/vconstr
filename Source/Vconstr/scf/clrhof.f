      subroutine clrhof(rho,drho,d2rho,nspin,norb,npnt,npntmx,pmo,
     + valmo,grid,weight)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical ltrian
      dimension rho(npntmx,nspin),drho(npntmx,3,nspin),
     + d2rho(npntmx,6,nspin)
      dimension pmo(norb*(norb+1)/2,nspin),valmo(npnt*norb)
      dimension grid(npntmx,3),weight(npnt)
      dimension gxmo(norb),gymo(norb),gzmo(norb),gxxmo(norb),
     + gxymo(norb),gxzmo(norb),gyymo(norb),gyzmo(norb),gzzmo(norb),
     + check(10,nspin)
c
      ltrian = .false.
c
      open(11,file='hfvls')
      rewind(11)
c
      check(1:10,1:nspin) = 0.d0
      do 20 ipnt=1,npnt
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
        do 10 is=1,nspin
          rho(ipnt,is)=valmat(pmo(1,is),norb,valmo(m),ltrian)
          call valgrd(dxao,dyao,dzao,pmo(1,is),norb,valmo(m),
     + gxmo,gymo,gzmo)
          drho(ipnt,1,is)=dxao
          drho(ipnt,2,is)=dyao
          drho(ipnt,3,is)=dzao
          call val2grd(dxx,dxy,dxz,dyy,dyz,dzz,pmo(1,is),norb,
     + valmo(m),gxmo,gymo,gzmo,gxxmo,gxymo,gxzmo,gyymo,gyzmo,gzzmo)
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
