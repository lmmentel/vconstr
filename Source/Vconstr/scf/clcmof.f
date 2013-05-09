      subroutine clcmof(valmo,npnt,npntmx,norb,vmopao,grid,weight,
     + atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(eps  = 1.d-6)
c
      logical lorb,atmol4
      dimension valmo(npnt*norb)
      dimension vmopao(norb*norb)
      dimension grid(npntmx,3),weight(npnt)
      dimension chkvl(norb),chkgx(norb),chkgy(norb),
     + chkgz(norb),chkgxx(norb),chkgxy(norb),chkgxz(norb),
     + chkgyy(norb),chkgyz(norb),chkgzz(norb)
      dimension valao(norb),gradx(norb),grady(norb),
     +  gradz(norb),gradxx(norb),gradxy(norb),gradxz(norb),
     +  gradyy(norb),gradyz(norb),gradzz(norb)
      dimension gxmo(norb),gymo(norb),gzmo(norb),gxxmo(norb),
     + gxymo(norb),gxzmo(norb),gyymo(norb),gyzmo(norb),gzzmo(norb)
c
      chkvl(1:norb)  = 0.d0
      chkgx(1:norb)  = 0.d0
      chkgy(1:norb)  = 0.d0
      chkgz(1:norb)  = 0.d0
      chkgxx(1:norb) = 0.d0
      chkgxy(1:norb) = 0.d0
      chkgxz(1:norb) = 0.d0
      chkgyy(1:norb) = 0.d0
      chkgyz(1:norb) = 0.d0
      chkgzz(1:norb) = 0.d0
c
      open(11,file='hfvls')
      rewind(11)
      do 100 ipnt=1,npnt
        k=(ipnt-1)*norb
        m=k+1
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
        do i=1,norb
          chkvl(i)=chkvl(i)+w*valao(i)**2
          chkgx(i)=chkgx(i)-w*valao(i)*x*gradx(i)
          chkgy(i)=chkgy(i)-w*valao(i)*y*grady(i)
          chkgz(i)=chkgz(i)-w*valao(i)*z*gradz(i)
          chkgxx(i)=chkgxx(i)+
     + w*x**2*(valao(i)*gradxx(i)+gradx(i)**2)
          chkgxy(i)=chkgxy(i)+
     + w*x*y*(valao(i)*gradxy(i)+gradx(i)*grady(i))
          chkgxz(i)=chkgxz(i)+
     + w*x*z*(valao(i)*gradxz(i)+gradx(i)*gradz(i))
          chkgyy(i)=chkgyy(i)+
     + w*y**2*(valao(i)*gradyy(i)+grady(i)**2)
          chkgyz(i)=chkgyz(i)+
     + w*y*z*(valao(i)*gradyz(i)+grady(i)*gradz(i))
          chkgzz(i)=chkgzz(i)+
     + w*z**2*(valao(i)*gradzz(i)+gradz(i)**2)
        enddo
        call vecmat(valao,vmopao,norb,valmo(m))
        call vecmat(gradx,vmopao,norb,gxmo)
        call vecmat(grady,vmopao,norb,gymo)
        call vecmat(gradz,vmopao,norb,gzmo)
        call vecmat(gradxx,vmopao,norb,gxxmo)
        call vecmat(gradxy,vmopao,norb,gxymo)
        call vecmat(gradxz,vmopao,norb,gxzmo)
        call vecmat(gradyy,vmopao,norb,gyymo)
        call vecmat(gradyz,vmopao,norb,gyzmo)
        call vecmat(gradzz,vmopao,norb,gzzmo)
        write(11,*) (gxmo(i),i=1,norb)
        write(11,*) (gymo(i),i=1,norb)
        write(11,*) (gzmo(i),i=1,norb)
        write(11,*) (gxxmo(i),i=1,norb)
        write(11,*) (gxymo(i),i=1,norb)
        write(11,*) (gxzmo(i),i=1,norb)
        write(11,*) (gyymo(i),i=1,norb)
        write(11,*) (gyzmo(i),i=1,norb)
        write(11,*) (gzzmo(i),i=1,norb)
  100 continue
      close(11)
c
      write(6,'(/,'' HF-MO values calculated and stored on disk'')')
c
      lorb=.true.
      write(6,'(/'' Orbital errors upon integration'')')
      do i=1,norb
        t1=abs(chkvl(i)-1.d0)
        t2=abs(chkgx(i)-5.d-1)
        t3=abs(chkgy(i)-5.d-1)
        t4=abs(chkgz(i)-5.d-1)
        t5=abs(chkgxx(i)-1.d0)
        t6=abs(chkgxy(i)-5.d-1)
        t7=abs(chkgxz(i)-5.d-1)
        t8=abs(chkgyy(i)-1.d0)
        t9=abs(chkgyz(i)-5.d-1)
        t10=abs(chkgzz(i)-1.d0)
        if ((t1.gt.eps).or.(t2.gt.eps).or.(t3.gt.eps).or.(t4.gt.eps).or.
     + (t5.gt.eps).or.(t6.gt.eps).or.(t7.gt.eps).or.(t8.gt.eps).or.
     + (t9.gt.eps).or.(t10.gt.eps)) then
          write(6,'(i3,10(1x,g8.2))')i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
          lorb=.false.
        endif
      enddo
      if (lorb) write(6,'('' No errors greater than '',g12.2)')eps
c
      return
      end
