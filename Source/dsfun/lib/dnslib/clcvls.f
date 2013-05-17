      subroutine clcvls(dns,ddns,dsdns,vnuc,valmo,grdmo,npnt,npntmx,
     + norb,pnomo,vmopao,grid,weight,rnel)
      use basisModule
      use gridAOvaluesModule
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps    = 1.d-4)
c
      parameter(one    = 1.0d0)
      parameter(half   = 0.5d0)
      parameter(three  = 3.0d0)
c
      logical ltrian,lorb
      dimension dns(npnt),ddns(npntmx,3),dsdns(npnt),
     + vnuc(npnt),valmo(npnt*norb),grdmo(npntmx*norb,4)
      dimension pnomo(norb*(norb+1)/2),vmopao(norb*norb),
     + grid(npntmx,3),weight(npnt)
      dimension chkvl(norb),chkgx(norb),chkgy(norb),
     + chkgz(norb),chknb(norb)
      dimension valao(norb),gradx(norb),grady(norb),
     +  gradz(norb),grads(norb)
clmm..add system type to hold system information
      type(basisType)  :: basis

clmm..read basis set information
      call newBasis(basis)
c
      ltrian = .false.
c
      tmp  = 0.d0
      tmpr = 0.d0
      tmps = 0.d0
      tmpf = 0.d0
      vonw = 0.d0
      chkvl(1:norb) = 0.d0
      chkgx(1:norb) = 0.d0
      chkgy(1:norb) = 0.d0
      chkgz(1:norb) = 0.d0
      chknb(1:norb) = 0.d0
c
clmm..temrorary set number of points to 1
      npnt = 1 
      write(*,*) 'before the loop in -clcvls-' 
      write(*,*) 'number of points = ', npnt
      do 100 ipnt=1,npnt
        k=(ipnt-1)*norb
        m=k+1
        x=grid(ipnt,1)
        y=grid(ipnt,2)
        z=grid(ipnt,3)
        rs=x*x+y*y+z*z
        w=weight(ipnt)
clmm        if (atmol4) then
clmm          call aovlsw(x,y,z,valao,gradx,grady,gradz,grads,vn)
clmm        else
clmm          call aovlsv(x,y,z,valao,gradx,grady,gradz,grads,vn)
clmm        endif
clmm..call new routine for calculating values, gradients and laplacian
clmm..using gamess-us basis set format
        call AOvalueAtPoint(basis,x,y,z,valao,gradx,grady,gradz,
     & grads,vn)
        vnuc(ipnt)=vn
        call vecmat(valao,vmopao,norb,valmo(m))
        call vecmat(gradx,vmopao,norb,grdmo(m,1))
        call vecmat(grady,vmopao,norb,grdmo(m,2))
        call vecmat(gradz,vmopao,norb,grdmo(m,3))
        call vecmat(grads,vmopao,norb,grdmo(m,4))
        do i=1,norb
          chkvl(i)=chkvl(i)+w*valmo(k+i)*valmo(k+i)
          chkgx(i)=chkgx(i)-w*valmo(k+i)*x*grdmo(k+i,1)
          chkgy(i)=chkgy(i)-w*valmo(k+i)*y*grdmo(k+i,2)
          chkgz(i)=chkgz(i)-w*valmo(k+i)*z*grdmo(k+i,3)
          chknb(i)=chknb(i)+w*rs*(grdmo(k+i,1)*grdmo(k+i,1)+
     + grdmo(k+i,2)*grdmo(k+i,2)+grdmo(k+i,3)*grdmo(k+i,3)+
     + valmo(k+i)*grdmo(k+i,4))
        enddo
        dn=valmat(pnomo,norb,valmo(m),ltrian)
        call valgrd(dx,dy,dz,pnomo,norb,valmo(m),
     + grdmo(m,1),grdmo(m,2),grdmo(m,3))
        ds=valsgrd(pnomo,valmo(m),grdmo(m,1),grdmo(m,2),
     + grdmo(m,3),grdmo(m,4),norb)
        tmp=tmp+dn*w
        tmpr=tmpr+(x*dx+y*dy+z*dz)*w
        tmps=tmps+w*ds*(x*x+y*y+z*z)
        tmpf=tmpf+dn*(dx+dy+dz)*w
        if (dn.gt.1.d-20) vonw=vonw+w*(dx*dx+dy*dy+dz*dz)/(8*dn)
        dns(ipnt)=dn
        ddns(ipnt,1)=dx
        ddns(ipnt,2)=dy
        ddns(ipnt,3)=dz
        dsdns(ipnt)=ds
  100 continue
      
      write(*,*) 'after the loop in -clcvls-' 
c
      write(6,'(/'' density : '',g10.4,'' integration errors '',
     + 4(g10.4))')rnel,tmp-rnel,-tmpr/3-rnel,tmps/6-rnel,tmpf
      if (abs(tmp-rnel).gt.1.d-3) write(6,'(''WARNING; '',
     + ''seems like the basis orbitals are symmetry adapted'')')
      write(6,'('' Von Weizacker kinetic energy : '',g12.6)')
     + vonw
c
      write(6,'('' Orbital Errors upon Integration''/)')
      do 30 i=1,norb
        t1=abs(chkvl(i)-one)
        t2=abs(chkgx(i)-half)
        t3=abs(chkgy(i)-half)
        t4=abs(chkgz(i)-half)
        t5=abs(chknb(i)-three)
        if ((t1.gt.eps).or.(t2.gt.eps).or.(t3.gt.eps).or.(t4.gt.eps).or.
     + (t5.gt.eps)) then
          write(6,'(i3,5(1x,g8.2))')i,t1,t2,t3,t4,t5
          lorb=.false.
        endif
   30 continue
      if (lorb) write(6,'('' No errors greater than '',g12.2)')eps
c
clmm..deallocate arrays with basis set information
      call deleteBasis(basis)
      return
      end
