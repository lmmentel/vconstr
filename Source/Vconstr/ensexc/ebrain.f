      subroutine ebrain(norb,npnt,intpnt,rnel,idmp,ismo,isks,isao,
     + atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-8)
      parameter(nfomx = 2)
      parameter(nspin = 1)
c
      logical atmol4
      dimension pksks(norb*(norb+1)/2),pksmo(norb*(norb+1)/2),
     + pksimo(norb*(norb+1)/2,nfomx),vksao(norb*norb),
     + vksmo(norb*norb),vaopao(norb*norb),vmoao(norb*norb),
     + vmopao(norb*norb),scrtc(norb*norb),occmo(norb)
      dimension cdet(nfomx),ifromo(nfomx),detmat(nfomx,nfomx),
     + ensmat(nfomx,nfomx)
      dimension exc(npnt,8),eexc(npnt,8),evxc(npnt,7,nspin)
      dimension erho(npnt,nspin),edrho(npnt,3,nspin),
     + ed2rho(npnt,6,nspin)
      dimension grid(npnt,3),weight(npnt)      
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
      do ip=1,npnt
        read(99,*) (grid(ip,i),i=1,3),weight(ip)
      enddo
c
      exld=0.d0
      ecld=0.d0
      exb=0.d0
      expw=0.d0
      ecp=0.d0
      ecpw=0.d0
      eclyp=0.d0
      ehcth=0.d0
c
      exc(1:npnt,1:8)=0.d0
c
      do ifo=1,nfomos
        eexc(1:npnt,1:8)=0.d0
        call clcrho(erho,edrho,ed2rho,nspin,norb,npnt,npnt,
     + pksimo(1,ifo),vmopao,grid,weight,atmol4)
        call clcvxc(nspin,npnt,npnt,erho,edrho,ed2rho,evxc,eexc)
        do ip=1,npnt
          exld=exld+cdet(ifo)*eexc(ip,1)*weight(ip)
          ecld=ecld+cdet(ifo)*eexc(ip,2)*weight(ip)
          exb=exb+cdet(ifo)*eexc(ip,3)*weight(ip)
          expw=expw+cdet(ifo)*eexc(ip,4)*weight(ip)
          ecp=ecp+cdet(ifo)*eexc(ip,5)*weight(ip)
          ecpw=ecpw+cdet(ifo)*eexc(ip,6)*weight(ip)
          eclyp=eclyp+cdet(ifo)*eexc(ip,7)*weight(ip)
c
c*** HCTH-Functional
c
          rha=5.d-1*erho(ip,1)
          drxa=5.d-1*edrho(ip,1,1)
          drya=5.d-1*edrho(ip,2,1)
          drza=5.d-1*edrho(ip,3,1)
          rhb=5.d-1*erho(ip,1)
          drxb=5.d-1*edrho(ip,1,1)
          dryb=5.d-1*edrho(ip,2,1)
          drzb=5.d-1*edrho(ip,3,1)
          za=sqrt(drxa**2+drya**2+drza**2)
          zb=sqrt(drxb**2+dryb**2+drzb**2)
          zab=sqrt(drxa*drxb+drya*dryb+drza*drzb)
          call hcth(dfdra,dfdza,dfdrb,dfdzb,dfdzab,rha,rhb,
     + za,zb,zab,.true.,totalF_xc)
          ehcth=ehcth+cdet(ifo)*totalF_xc*weight(ip)
          eexc(ip,8)=totalF_xc
c
        enddo
c
      write(6,'(/''** Exchange energies **'')')
      write(6,'(6x,''LDA    :'',12x,f10.6)')exld
      write(6,'(6x,''Becke  :'',f10.6,2x,f10.6)')exb,exld+exb
      write(6,'(6x,''PW91   :'',f10.6,2x,f10.6)')expw,exld+expw
c
      write(6,'(/''** Correlation energies **'')')
      write(6,'(6x,''LDA    :'',12x,f10.6)')ecld
      write(6,'(6x,''LYP    :'',12x,f10.6)')eclyp
      write(6,'(6x,''Perdew :'',f10.6,2x,f10.6)')ecp,ecld+ecp
      write(6,'(6x,''PW91   :'',f10.6,2x,f10.6)')ecpw,ecld+ecpw
c         
      write(6,'(/6x,''hcthxc :'',12x,f10.6)')ehcth
c
        do k=1,8
          exc(1:npnt,k)=exc(1:npnt,k)+
     + cdet(ifo)*eexc(1:npnt,k)/erho(1:npnt,1)
        enddo
      enddo
c
      exc(1:npnt,3)=exc(1:npnt,1)+exc(1:npnt,3)
      exc(1:npnt,4)=exc(1:npnt,1)+exc(1:npnt,4)
      exc(1:npnt,5)=exc(1:npnt,2)+exc(1:npnt,5)
      exc(1:npnt,6)=exc(1:npnt,2)+exc(1:npnt,6)
c
      open(84,file='fexc.dat') 
      rewind(84)
      do ip=intpnt,npnt
        write(84,*)exc(ip,1),exc(ip,2),exc(ip,3),exc(ip,4),
     + exc(ip,5),exc(ip,6),exc(ip,7)
      enddo
c
      write(6,'(/''** Exchange energies **'')')
      write(6,'(6x,''LDA    :'',12x,f10.6)')exld
      write(6,'(6x,''Becke  :'',f10.6,2x,f10.6)')exb,exld+exb
      write(6,'(6x,''PW91   :'',f10.6,2x,f10.6)')expw,exld+expw
c
      write(6,'(/''** Correlation energies **'')')
      write(6,'(6x,''LDA    :'',12x,f10.6)')ecld
      write(6,'(6x,''LYP    :'',12x,f10.6)')eclyp
      write(6,'(6x,''Perdew :'',f10.6,2x,f10.6)')ecp,ecld+ecp
      write(6,'(6x,''PW91   :'',f10.6,2x,f10.6)')ecpw,ecld+ecpw
c
      write(6,'(/6x,''hcthxc :'',12x,f10.6)')ehcth
c
      write(6,'(//)')
c
      return
      end
