      subroutine cledet(occmo,pksks,vksmo,vmopao,rhotot,grid,weight,
     + valmo,ndets,nmos,norb,npnt,intpnt,iint,vxc,exc,lb88,lpw91x,
     + lprdwc,lpw91c,llyp,ldisk,atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-6)
c
      logical lb88,lpw91x,lprdwc,lpw91c,llyp,ldisk,atmol4
      dimension occmo(nmos)
      dimension pksks(norb*(norb+1)/2)
      dimension vksmo(norb*norb),vmopao(norb*norb)
      dimension rhotot(npnt),grid(npnt,3),weight(npnt)
      dimension valmo(npnt*norb)
      dimension vxc(npnt),exc(npnt,4)
c
      dimension pksmo(norb*(norb+1)/2)
      dimension rho(npnt,1),drho(npnt,3,1),d2rho(npnt,6,1)
      dimension vxcens(npnt),excens(npnt,4)
      dimension cdet(ndets),ifromo(ndets),detmat(ndets,ndets),
     + ensmat(ndets,ndets)
c
      nspin=1
      npntmx=npnt
c
      ifr=1
c
      rn=0.d0
      occorb=2.d0
c
      exc(1:npnt,1:4)=0.d0
      vxc(1:npnt)=0.d0
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
      ehks=0.d0
      exks=0.d0
      do i=1,ndets
        cdet(i)=0.d0
        do j=1,ndets
          jmo=ifromo(j)
c
c** the coefficients in front of each occmoupation vector correspond to a row
c** vector, as each row in the matrix detmat was an occupation vector
c
          cdet(i)=cdet(i)+occmo(ifromo(j))*ensmat(j,i)
          pksks(jmo*(jmo+1)/2)=detmat(i,j)
        enddo
        write(6,'(/'' '',i2,''   '',f6.4)')i,cdet(i)
        call tmtdag(pksks,nmos,pksmo,norb,vksmo,5.d-1)
        call wmtrx('####  ks-density matrix in mo basis  ####',
     + pksmo,norb,1.d-2)
        if (ldisk) then
          call clrhof(rho,drho,d2rho,nspin,norb,npnt,npntmx,
     + pksmo,valmo,grid,weight)
        else
          call clrhod(rho,drho,d2rho,nspin,norb,npnt,npntmx,
     + pksmo,vmopao,valmo,grid,weight,atmol4)
        endif
        call clcvxc(nspin,npntmx,npnt,rho,drho,d2rho,
     + vxcens,excens,lb88,lpw91x,lprdwc,lpw91c,llyp)
        call erep(pksmo,norb,iint,ectmp,extmp)
        call clcexc(exld,exnl,ecld,ecnl,intpnt,npntmx,excens,weight)
        call prtres(ectmp,exld,exnl,extmp,ecld,ecnl)
        do j=1,4
          exc(1:npnt,j)=exc(1:npnt,j)+cdet(i)*excens(1:npnt,j)
        enddo
        vxc(1:npnt)=vxc(1:npnt)+
     + cdet(i)*rho(1:npnt,1)*vxcens(1:npnt)/rhotot(1:npnt)
        ehks=ehks+cdet(i)*ectmp
        exks=exks+cdet(i)*extmp
      enddo
      write(6,'(//''  Final calculation ''/)')
      call clcexc(exld,exnl,ecld,ecnl,intpnt,npntmx,exc,weight)
      call prtres(ehks,exld,exnl,exks,ecld,ecnl)
c
      return
      end
c
      subroutine prtres(ehks,exld,exnl,exks,ecld,ecnl)
      implicit real*8 (a-h,o-z),integer(i-n)

      write(6,'('' KS coulomb repulsion energy            :'',  
     +f16.10)') ehks
      write(6,'(/'' KS local density exchange energy       :'',
     +f16.10)')exld
      write(6,'('' KS non-local exchange energy           :'',
     +f16.10)')exnl
      write(6,'('' KS total exchange energy               :'',
     +f16.10)')exld+exnl
      write(6,'('' KS orbital-exchange energy             :'',
     +f16.10)')exks
c       
      write(6,'(/'' KS local density correlation energy    :'',
     +f16.10)')ecld
      write(6,'('' KS non-local correlation energy        :'',
     +f16.10)')ecnl
      write(6,'('' KS total correlation energy            :'',
     +f16.10)')ecld+ecnl
c
      return
      end
