      subroutine ebrain(norb,rnel,idmp,ismo,isks,isao,atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-8)
      parameter(nfomx = 2)
c
      logical atmol4,ltrian
      dimension pksks(norb*(norb+1)/2),pksmo(norb*(norb+1)/2),
     + pksimo(norb*(norb+1)/2,nfomx),vksao(norb*norb),
     + vksmo(norb*norb),vaopao(norb*norb),vmoao(norb*norb),
     + vmopao(norb*norb),scrtc(norb*norb),
     + valmo(norb),valao(norb),occmo(norb)
      dimension ensrho(nfomx),cdet(nfomx),ifromo(nfomx),
     + detmat(nfomx,nfomx),ensmat(nfomx,nfomx)
c
      dimension pkstmp(norb*(norb+1)/2)
c
c** read dumpfile 
c
      call getdmp(idmp,ismo,vmoao,occmo,norb)
      rnel=sum(occmo)
c
c** vmoao contains the orbitals in (symmetry adapted or primitive) ao basis
c
      if (isao.gt.0) then
        call getdmp(idmp,isao,vaopao,occmo,norb)
c
c** transform MO orbitals from symmetry adapted to primitive ao basis 
c**  vmopao = vmoao*vaopao
c
        call matml(vmoao,vaopao,vmopao,norb)
      else
        vmopao(1:norb*norb)=vmoao(1:norb*norb)
      endif
c
c** get Kohn-Sham Orbitals in ao basis
c
      call getdmp(idmp,isks,vksao,occmo,norb)
c
c** transform natural orbitals from ao to mo basis
c
      scrtc(1:norb*norb)=vmoao(1:norb*norb)
c
c** invert the matrix scrtc => scrtc
c
      tol = 1.0d-7
      call minvr(scrtc,tol,det,ier,norb)
      if (ier.ne.0) then
        write(6,'(/,'' ERROR : transformation from ao to '',
     +       ''mo basis is singular'')')
        stop
      endif
c
c** transform Kohn-Sham orbitals from ao to mo basis
c**   vksmo = vksao*scrtc
c
      call matml(vksao,scrtc,vksmo,norb)
c
      mvpr = 0
      sno = 0.d0
      pksks(1:norb*(norb+1)/2)=0.d0
      do i=1,norb
        ii=i*(i+1)/2
        sno = sno+occmo(i)
        if (occmo(i).gt.1.d-6) mvpr=mvpr+1
        pksks(ii)=occmo(i)
      enddo
      if (abs(sno-rnel).gt.eps) then
        write(6,'(/''ERROR; sum of natorb occupation numbers :'',
     + f8.3)') sno
        write(6,'(''       number of electrons              :'',
     + f8.3)') rnel
        stop
      endif
c
      do imo=norb,1,-1
        if (occmo(imo).gt.1.d-6) then
          nmos=imo
          exit
        endif
      enddo
      call tmtdag(pksks,nmos,pksmo,norb,vksmo,5.d-1)
c
c      call prvec(vksmo,norb,mvpr,' Kohn-Sham orbitals in mo-basis ')
      call wmtrx('####  ks-density matrix in mo basis  ####',
     + pksmo,norb,1.d-2)
c
c** ensemble representation
c
      ifr=1
      nfomos=0
c
      rn=0.d0
      occorb=2.d0
c
      do imo=1,nmos
        if ((abs(occmo(imo)-occorb).gt.eps).and.
     + (abs(occmo(imo)).gt.eps)) nfomos=nfomos+1
      enddo
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
      pkstmp(1:norb*(norb+1)/2)=0.d0
      do i=1,nfomos
        cdet(i)=0.d0
        do j=1,nfomos
          jmo=ifromo(j)
c
c** the coefficients in front of each occmoupation vector correspond to a row
c** vector, as each row in the matrix detmat was an occupation vector
c
          cdet(i)=cdet(i)+occmo(ifromo(j))*ensmat(j,i)
          pksks(jmo*(jmo+1)/2)=detmat(i,j)
        enddo
        write(6,'(/'' '',i2,''   '',f6.4)')i,cdet(i)
        call tmtdag(pksks,nmos,pksimo(1,i),norb,vksmo,5.d-1)
        call wmtrx('####  ks-density matrix in mo basis  ####',
     + pksimo(1,i),norb,1.d-2)
        pkstmp(1:norb*(norb+1)/2)=pkstmp(1:norb*(norb+1)/2)+
     + cdet(i)*pksimo(1:norb*(norb+1)/2,i)
      enddo
      call wmtrx('####  ks-density matrix in mo basis  ####',
     + pkstmp,norb,1.d-2)
c
      ltrian = .false.
c
      open(99,file='points')
      rewind 99
      read(99,*)npnt,npold
      write(6,'(/,'' number of gridpoints in numerical'',
     + '' integration :'',i5)')npold-1
c
      open(20,file='vhartr1')
      rewind(20)
      open(21,file='vhartr2')
      rewind(21)
      open(22,file='vhartr.ks')
      rewind(22)
c
      open(23,file='vcond1')
      rewind(23)
      open(24,file='vcond2')
      rewind(24)
      open(25,file='vcond.ks')
      rewind(25)
c
      do ipnt=1,npnt
        read(99,*) x,y,z,w
        if (atmol4) then
          call aovlw(x,y,z,valao)
        else
          call aovlv(x,y,z,valao)
        endif
        call vecmat(valao,vmopao,norb,valmo)
        rho=valmat(pksmo,norb,valmo,ltrian)
c
        tmp=0.d0
        do ifo=1,nfomos
          ensrho(ifo)=valmat(pksimo(1,ifo),norb,valmo,ltrian)
          tmp=tmp+cdet(ifo)*ensrho(ifo)
        enddo
c
        read(20,*) vH1
        read(21,*) vH2
        vHtmp=(cdet(1)*ensrho(1)*vH1+cdet(2)*ensrho(2)*vH2)/rho
        write(22,*) vHtmp
c
        read(23,*) vc1
        read(24,*) vc2
        vctmp=(cdet(1)*ensrho(1)*vc1+cdet(2)*ensrho(2)*vc2)/rho
        write(25,*) vctmp
c
      enddo
c
      write(6,'(//)')
c
      return
      end
