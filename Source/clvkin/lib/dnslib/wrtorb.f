      subroutine wrtorb(vecks,occmo,ityp,ndvcr,nmos,norb,npnt,npold,
     + valmo)
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension vecks(norb*norb),occmo(nmos),ityp(nmos),valmo(norb*npnt)
      dimension orbdns(ndvcr)
c
      open(55,file='orb.plt')
      rewind(55)
      open(56,file='orb.dat')
      rewind(56)
      do ip=1,npnt
        orbdns(1:ndvcr)=0.d0
        k=(ip-1)*norb
        do i=1,nmos
          it=ityp(i)
          if (it.ge.0) then
            tmp=0.d0
            ii=(i-1)*norb
            do j=1,norb
              tmp=tmp+vecks(ii+j)*valmo(k+j)
            enddo
            orbdns(it)=orbdns(it)+occmo(i)*tmp*tmp
          endif
        enddo
        if (ip.ge.npold) write(55,'(10(f12.6,2x))')
     + (orbdns(i),i=1,ndvcr)
        write(56,*)(orbdns(i),i=1,ndvcr)
      enddo
      close(55)
      close(56)
c
      return
      end

