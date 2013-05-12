      subroutine clvmat(vmat,npnt,npntmx,norb,vxc,valmo,weight)
c***********************************************************************
c
c Calculate matrix elements of potential vxc in HF-MO basis. 
c The matrix elements are returned in array vmat
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter (nspin = 2)
c
      dimension vmat(norb*(norb+1)/2,nspin)
      dimension vxc(npntmx,nspin),valmo(norb*npnt),weight(npnt)
c
      vmat(1:norb*(norb+1)/2,1:nspin)=0.d0
c
      do ipnt=1,npnt
        k=(ipnt-1)*norb
        w=weight(ipnt)
        x1=w*vxc(ipnt,1)
        x2=w*vxc(ipnt,2)
        do i=1,norb
          ii=i*(i-1)/2
          vi=valmo(k+i)
          y1=vi*x1
          y2=vi*x2
          do j=1,i
            vj=valmo(k+j)
            vmat(ii+j,1)=vmat(ii+j,1)+vj*y1
            vmat(ii+j,2)=vmat(ii+j,2)+vj*y2
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine clvmts(vmat,vxc,valmo,weight,norbsm,nsym,npnt,npntmx,
     + norb)
c***********************************************************************
c
c Calculate matrix elements of potential vxc in HF-MO basis. 
c The matrix elements are returned in array vmat
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter (nspin=2)
c
      dimension vmat(norb*(norb+1)/2,nspin)
      dimension vxc(npntmx,nspin),valmo(norb*npnt),weight(npnt),
     + norbsm(1:norb,0:nsym)
c
      vmat(1:norb*(norb+1)/2,1:nspin) = 0.d0
c
      do isym=1,nsym
        do ipnt=1,npnt
          k=(ipnt-1)*norb
          x1=weight(ipnt)*vxc(ipnt,1)
          x2=weight(ipnt)*vxc(ipnt,2)
          do i=1,norbsm(isym,0)
            iorb=norbsm(i,isym)
            ii=iorb*(iorb-1)/2
            y1=x1*valmo(k+iorb)
            y2=x2*valmo(k+iorb)
            do j=1,i
              jorb=norbsm(j,isym)
              vmat(ii+jorb,1)=vmat(ii+jorb,1)+y1*valmo(k+jorb)
              vmat(ii+jorb,2)=vmat(ii+jorb,2)+y2*valmo(k+jorb)
            enddo
          enddo
        enddo
      enddo
c
      return
      end
