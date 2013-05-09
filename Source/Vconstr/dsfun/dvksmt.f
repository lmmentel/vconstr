      subroutine dvksmt(dvxc,npnt,nmos,ndvcr,rho,orbdnst,Avec,weight,
     + occmo,ityp)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(fourth=1.0d0/4.0d0)
      parameter(eps=1.0d-12)
c
      dimension dvxc(nmos)
      dimension ityp(nmos)
      dimension occmo(nmos)
      dimension rho(npnt),orbdnst(npnt*ndvcr),Avec(npnt),weight(npnt)
      dimension sltmat(ndvcr,ndvcr),gdrmat(ndvcr),occtyp(ndvcr)
c
      sltmat(1:ndvcr,1:ndvcr)=0.d0
c
c** first diagonal
c
      do ipnt=1,npnt
        k=(ipnt-1)*ndvcr
        do i=1,ndvcr
          do j=1,i
            if (rho(ipnt).gt.eps) then
              sltmat(i,j)=sltmat(i,j)+
     + orbdnst(k+i)*orbdnst(k+j)*weight(ipnt)/rho(ipnt)
            elseif (abs(rho(ipnt)-orbdnst(k+i)).lt.eps) then
              sltmat(i,j)=sltmat(i,j)+
     + orbdnst(k+j)*weight(ipnt)
            endif
          enddo
        enddo
      enddo
c
c** convert to complete matrix
c
      do i=1,ndvcr
        do j=1,i-1
          sltmat(j,i)=sltmat(i,j)
        enddo
      enddo
c
c** subtract occupation from diagonal elements
c
      occtyp(1:ndvcr)=0.d0
      do i=1,nmos
        l=ityp(i)
        if (l.ne.0) occtyp(l)=occtyp(l)+occmo(i)
      enddo
      do i=1,ndvcr
        sltmat(i,i)=sltmat(i,i)-occtyp(i)
      enddo
c
c** subtract one for highest occupied orbitals
c
      sltmat(1:ndvcr,ndvcr)=sltmat(1:ndvcr,ndvcr)-occtyp(1:ndvcr)
c
c** invert the matrix sltmat => sltmat
c
      tol=1.d-20
      call minvr(sltmat,tol,det,ier,ndvcr)
      if (ier.eq.1) then
        write(6,'(/,'' ERROR : density response function could not'',
     +        '' be inverted.'')')
        stop
      endif
c
      gdrmat(1:ndvcr)=0.d0
      do ipnt=1,npnt
        k=(ipnt-1)*ndvcr
        do i=1,ndvcr
          if (rho(ipnt).gt.eps) then
            gdrmat(i)=gdrmat(i)-fourth*orbdnst(k+i)*
     + Avec(ipnt)*weight(ipnt)/rho(ipnt)
          elseif (i.eq.ndvcr) then
            gdrmat(i)=gdrmat(i)-fourth*Avec(ipnt)*weight(ipnt)
          endif
        enddo
      enddo
c
      dvxc(1:ndvcr)=0.d0
      do i=1,ndvcr
        do j=1,ndvcr
          dvxc(i)=dvxc(i)+sltmat(i,j)*gdrmat(j)
        enddo
      enddo
c
      return
      end
