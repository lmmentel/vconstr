      subroutine vksmat(vmat,npnt,norb,vxc,valmo,weight)
c***********************************************************************
c
c Calculate matrix elements of potential vxc in HF-MO basis. 
c The matrix elements are returned in array vmat
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension vmat(norb*(norb+1)/2)
      dimension vxc(npnt),valmo(norb*npnt),weight(npnt)
c
      vmat(1:norb*(norb+1)/2) = 0.d0
c
      do ipnt=1,npnt
        k=(ipnt-1)*norb
        x=weight(ipnt)*vxc(ipnt)
        do i=1,norb
          ii=i*(i-1)/2
          y=x*valmo(k+i)
          do j=1,i
            vmat(ii+j)=vmat(ii+j)+y*valmo(k+j)
          enddo
        enddo
      enddo
c
      return
      end
c
      subroutine vksmts(vmat,vxc,valmo,weight,norbsm,nsym,npnt,norb)
c***********************************************************************
c
c Calculate matrix elements of potential vxc in HF-MO basis. 
c The matrix elements are returned in array vmat
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension vmat(norb*(norb+1)/2)
      dimension vxc(npnt),valmo(norb*npnt),weight(npnt),
     + norbsm(1:norb,0:nsym)
c
      vmat(1:norb*(norb+1)/2) = 0.d0
c
      do isym=1,nsym
        do ipnt=1,npnt
          k=(ipnt-1)*norb
          x=weight(ipnt)*vxc(ipnt)
          do i=1,norbsm(isym,0)
            iorb=norbsm(i,isym)
            ii=iorb*(iorb-1)/2
            y=x*valmo(k+iorb)
            do j=1,i
              jorb=norbsm(j,isym)
              vmat(ii+jorb)=vmat(ii+jorb)+y*valmo(k+jorb)
            enddo
          enddo
        enddo
      enddo
c
      return
      end
