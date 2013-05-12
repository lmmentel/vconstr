      subroutine cledet(occmo,pdndn,vdnmo,vmopao,rhotot,grid,weight,
     + rnel,ndets,nmos,norb,npnt,intpnt,vxc,exc,atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-6)
c
      logical atmol4
      dimension occmo(nmos)
      dimension pdndn(norb*(norb+1)/2)
      dimension vdnmo(norb*norb),vmopao(norb*norb)
      dimension rhotot(npnt),grid(npnt,3),weight(npnt)
      dimension vxc(npnt,7),exc(npnt,7)
c
      dimension pdnmo(norb*(norb+1)/2)
      dimension rho(npnt,1),drho(npnt,3,1),d2rho(npnt,6,1)
      dimension vxcens(npnt,7),excens(npnt,7)
      dimension cdet(ndets),ifromo(ndets),detmat(ndets,ndets),
     + ensmat(ndets,ndets)
c
      ifr=1
c
      rn=0.d0 
      occorb=2.d0
c
      exc(1:npnt,1:7)=0.d0
      vxc(1:npnt,1:7)=0.d0
c
      do imo=1,nmos
        if ((abs(occmo(imo)).gt.eps).and.
     + (abs(occmo(imo)-occorb).gt.eps)) then
          ifromo(ifr)=imo
          ifr=ifr+1
          rn=rn+occmo(imo)
        endif
      enddo
      write(6,'(//'' Energy determined from '',i4,
     + '' degenerate determinants.''/)')ndets
      nl=nint(rn)/2
      detmat(1:ndets,1:ndets)=0.d0
      detmat(1:ndets,1:nl)=occorb
c     
c** each row is an occupation vector
c
      imat=2
      do imo=1,nl
        do jmo=nl+1,ndets
          detmat(imat,imo)=0.d0
          detmat(imat,jmo)=occorb
          imat=imat+1
        enddo
      enddo
      ensmat(1:ndets,1:ndets)=detmat(1:ndets,1:ndets)
c
c** invert the matrix ensmat => ensmat
c
      tol=1.d-20
      call minvr(ensmat,tol,det,ier,ndets)
      if (ier.eq.1) then  
        write(6,'(/,'' ERROR : ensembles could not be determined.'')')
        stop 
      endif
c
      do i=1,ndets
        cdet(i)=0.d0
        do j=1,ndets
          jmo=ifromo(j)
c
c** the coefficients in front of each occmoupation vector correspond to a row
c** vector, as each row in the matrix detmat was an occupation vector
c
          cdet(i)=cdet(i)+occmo(ifromo(j))*ensmat(j,i)
          pdndn(jmo*(jmo+1)/2)=detmat(i,j)
        enddo
        write(6,'(/'' '',i2,''   '',f6.4)')i,cdet(i)
        call tmtdag(pdndn,nmos,pdnmo,norb,vdnmo,5.d-1)
        call wmtrx('####  ks-density matrix in mo basis  ####',
     + pdnmo,norb,1.d-2)
        call clcrho(rho,drho,d2rho,norb,npnt,npnt,rnel,pdnmo,
     + vmopao,grid,weight,atmol4)
        call clcvxc(npnt,npnt,rho,drho,d2rho,vxcens,excens)
        call clcexc(intpnt,npnt,excens,weight)
        do j=1,7
          exc(1:npnt,j)=exc(1:npnt,j)+cdet(i)*excens(1:npnt,j)
          vxc(1:npnt,j)=vxc(1:npnt,j)+ 
     + cdet(i)*rho(1:npnt,1)*vxcens(1:npnt,j)/rhotot(1:npnt)
        enddo
      enddo
      write(6,'(//''  Final calculation ''/)')
      call clcexc(intpnt,npnt,exc,weight)
c
      return
      end
