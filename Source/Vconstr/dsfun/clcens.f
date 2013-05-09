      subroutine clcens(vksmo,pksks,occmo,occorb,eone,norb,nmos,nfomos,
     + iint)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps    = 1.0d-8)
c
      dimension vksmo(norb*norb),pksks(norb*(norb+1)/2),occmo(nmos)
      dimension detmat(nfomos,nfomos),ensmat(nfomos,nfomos),
     + ifromo(nfomos)
      dimension scrtc(norb*(norb+1)/2)
c
c** ensemble representation
c
      ifr=1
      c = 5.d-1
      exks=0.d0
      ecks=0.d0
c
      rn=0.d0
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
     + '' degenerate determinants.''/)')nfomos
      nl=nint(rn)/2
      detmat(1:nfomos,1:nfomos)=0.d0
      detmat(1:nfomos,1:nl)=occorb
c
c** each row is an occupation vector
c
      imat=2
      do imo=1,nl
        do jmo=nl+1,nfomos
          detmat(imat,imo)=0.d0
          detmat(imat,jmo)=occorb
          imat=imat+1
        enddo
      enddo
      ensmat(1:nfomos,1:nfomos)=detmat(1:nfomos,1:nfomos)
c
c** invert the matrix ensmat => ensmat
c
      tol=1.d-20
      call minvr(ensmat,tol,det,ier,nfomos)
      if (ier.eq.1) then
        write(6,'(/,'' ERROR : ensembles could not be determined.'')')
        stop
      endif
c
      do i=1,nfomos
        cdet=0.d0
        do j=1,nfomos
          jmo=ifromo(j)
c
c** the coefficients in front of each occupation vector correspond to a row
c** vector, as each row in the matrix detmat was an occupation vector
c
          cdet=cdet+occmo(ifromo(j))*ensmat(j,i)
          pksks(jmo*(jmo+1)/2)=detmat(i,j)
        enddo
        call tmtdag(pksks,nmos,scrtc,norb,vksmo,c)
        call erep(scrtc,norb,iint,ectmp,extmp)
        write(6,'('' '',i2,''   '',f6.4,'' '',f8.4,''   '',f8.4,
     + ''   '',3(15f3.0,/))')i,cdet,ectmp,extmp,
     + (pksks(k*(k+1)/2),k=1,nmos)
        exks=exks+cdet*extmp
        ecks=ecks+cdet*ectmp
      enddo
c
      write(6,'(//)')
      write(6,'('' KS coulomb repulsion energy            :''
     +,f16.10)')ecks
      write(6,'('' KS exchange energy                     :''
     +,f16.10)')exks
      write(6,'('' KS energy expectation value            :''
     +,f16.10)')eone+ecks+exks
c
      return
      end
