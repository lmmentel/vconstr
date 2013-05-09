      subroutine clcdvx(wvxc,wovs,rho,orbdns,weight,occtyp,npnt,npntmx,
     + ndvx)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(eps=1.0d-12)
c
      dimension wvxc(ndvx)
      dimension wovs(ndvx),occtyp(ndvx)
      dimension rho(npnt),orbdns(npntmx,ndvx),weight(npnt)
      dimension sltmat(ndvx,ndvx)
      dimension eresp(ndvx)
      dimension am(ndvx,ndvx)
c
      sltmat(1:ndvx,1:ndvx)=0.d0
c
c** first diagonal
c
      do i=1,ndvx
        do j=1,i
          do ip=1,npnt
            if (rho(ip).gt.eps) then
              sltmat(i,j)=sltmat(i,j)-
     + orbdns(ip,i)*orbdns(ip,j)*weight(ip)/rho(ip)
            elseif (j.eq.ndvx) then
              sltmat(i,j)=sltmat(i,j)-orbdns(ip,i)*weight(ip)
            endif
          enddo
        enddo
      enddo
c
c** convert to complete matrix
c
      do i=1,ndvx
        do j=1,i-1
          sltmat(j,i)=sltmat(i,j)
        enddo
      enddo
      do i=1,ndvx
      do j=1,ndvx
      am(i,j)=-sltmat(i,j)
      enddo
      enddo
c
c** add occupation to diagonal elements
c
      do i=1,ndvx
        sltmat(i,i)=sltmat(i,i)+occtyp(i)
      enddo
c
c** add occupation for highest occupied orbitals
c
      sltmat(1:ndvx,ndvx)=sltmat(1:ndvx,ndvx)+occtyp(1:ndvx)
c
c** invert the matrix sltmat => sltmat
c
      tol=1.d-20
      call minvr(sltmat,tol,det,ier,ndvx)
      if (ier.eq.1) then
        write(6,'(/,'' ERROR : Unable to invert Slater matrix'')')
        stop
      endif
c
      wvxc(1:ndvx)=0.d0
      do i=1,ndvx
        do j=1,ndvx
          wvxc(i)=wvxc(i)+sltmat(i,j)*wovs(j)
        enddo
      enddo
c
      do i=1,ndvx
      eresp(i)=0.d0
      do j=1,ndvx
      eresp(i)=eresp(i)+wvxc(j)*am(i,j)/occtyp(i)
      enddo
      enddo
      write(6,'(/''  erv  : '',10(1x,f11.6))')
     +(eresp(i),i=1,ndvx)
      return
      end
