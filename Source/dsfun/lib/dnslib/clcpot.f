      subroutine clcpot(norb,nmos,npnt,vksmo,hksmat,occmo,vhartr,vnuc,
     + fxyz,valmo,grdsmo,grid)
c***********************************************************************
c
c Construct the KS-potential from the orbitals
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension vksmo(norb*norb),hksmat(norb*norb),occmo(nmos),
     + vhartr(npnt),vnuc(npnt),fxyz(3),grid(npnt,3)
      dimension valmo(norb*npnt),grdsmo(norb*npnt)
c
      open(91,file='vxact.dat')
      rewind(91)
      do ip=1,npnt
        k=(ip-1)*norb
        vs=0.d0
        rho=0.d0
        do i=1,nmos
          iimo=(i-1)*norb
          tmpd=0.d0
          tmps=0.d0
          do j=1,norb
            tmpd=tmpd+vksmo(iimo+j)*valmo(k+j)
            tmps=tmps+vksmo(iimo+j)*grdsmo(k+j)
          enddo
          ei=hksmat((i-1)*norb+i)
          vs=vs+occmo(i)*tmpd*(tmps/2+ei*tmpd)
          rho=rho+occmo(i)*tmpd*tmpd
        enddo
        vd=0.d0
        do i=1,3
          vd=vd+fxyz(i)*grid(ip,i)
        enddo
        vst=vs/rho-vhartr(ip)+vnuc(ip)+vd
        write(91,*) vst
      enddo
      close(91)
c
      return
      end
