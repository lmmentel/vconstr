      subroutine clsorb(orbdns,vksmo,occmo,ityp,nspin,nmos,nmomx,
     + norb,npnt,valmo)
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension orbdns(npnt,nmomx,nspin)
      dimension vksmo(norb*norb,nspin),occmo(nmomx,nspin),
     + valmo(norb*npnt),ityp(nmomx,nspin)
c
      orbdns(1:npnt,1:nmos,1:nspin)=0.d0
      do is=1,nspin
        do i=1,nmos
          l=ityp(i,is)
          if (l.ne.0) then
            ii=(i-1)*norb
            do ip=1,npnt
              k=(ip-1)*norb
              tmp=0.d0
              do j=1,norb
                tmp=tmp+valmo(k+j)*vksmo(ii+j,is)
              enddo
              orbdns(ip,l,is)=orbdns(ip,l,is)+occmo(i,is)*tmp*tmp
            enddo
          endif
        enddo
      enddo
c
      return
      end
