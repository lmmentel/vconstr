      subroutine brains(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmn,
     + crrmx,thresh,dqmax,dvdmp,scfdmp,kpens,info,nppr,nvpr,npnt,npold,
     + norb,nmos,nmomx,itrx,ibcens,iscens,iint,idmp,ismo,isno,isao,
     + isoe,isks,lsym,lintsm,lrdocc,lrfun,nmosa,nmosb,atmol4)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.0d-8)
      parameter(fparin = 1.01d0)
c
      logical lsym,lintsm,lrdocc,lrfun,lincpr,atmol4
      dimension info(nppr)
      dimension occmo(nmomx),soccmo(nmomx,2)
c
      logical ltrian,lcens,lvhart
      dimension fxyz(3),dxyz(3)
      dimension ii1(8),jj1(8),ll1(8),kk1(8)
      dimension iorder(nmomx),ityp(nmomx),nmoss(2)
      dimension iold(norb),inew(norb)
      dimension qorb(nmomx),qorbld(kpens*nmomx),chkval(nmomx),
     + tsorb(nmomx)
      dimension pksks(norb*(norb+1)/2),pmo(norb*(norb+1)/2),
     + pmoao(norb*(norb+1)/2),pnomo(norb*(norb+1)/2),
     + pnoao(norb*(norb+1)/2),pksmo(norb*(norb+1)/2),
     + pksao(norb*(norb+1)/2),vxcmat(norb*(norb+1)/2),
     + vhrmat(norb*(norb+1)/2),hmatx(norb*(norb+1)/2),
     + vnmat(norb*(norb+1)/2),tsmat(norb*(norb+1)/2),
     + pksiks(norb*(norb+1)/2),pksimo(norb*(norb+1)/2),
     + pksiao(norb*(norb+1)/2),pvksmo(norb*(norb+1)/2)
      dimension pnono(norb*(norb+1)/2),vinomo((norb-nmosa)*norb),
     +pnomvo(norb*(norb+1)/2),pnovo((norb-nmosa)*(norb-nmosa+1)/2)
      dimension vksmo(norb*norb),vmoao(norb*norb),hksmat(norb*norb),
     +vmopao(norb*norb),pnona(norb*(norb+1)/2),
     +pnonb(norb*(norb+1)/2),pnova((norb-nmosa)*(norb-nmosa+1)/2),
     +vinoma((norb-nmosa)*norb),pnomva(norb*(norb+1)/2),
     +pnovb((norb-nmosb)*(norb-nmosb+1)/2),
     +vinomb((norb-nmosb)*norb),pnomvb(norb*(norb+1)/2)
      dimension pnomoa(norb*norb),pnomob(norb*norb),gaa(norb*norb),
     + gbb(norb*norb),gab(norb*norb),vnomoa(norb*norb),
     +vnoaoa(norb*norb), 
     +vnomob(norb*norb),occnoa(norb),occnob(norb),vnoma(norb*norb),
     +vnomb(norb*norb),occna(norb),occnb(norb),sab(norb*norb)
     +,vnom(norb*norb),vnomo(norb*norb),occno(norb),occn(norb) 
      dimension vnomos(norb*norb),vnoms(norb*norb),occns(norb),
     +occnos(norb),pnons(norb*(norb+1)/2),gabs(norb*norb)
      dimension rho(npnt),vxc(npnt),vnuc(npnt)
      dimension eorb(norb),iorbsm(norb)
      dimension grid(npnt,3),weight(npnt)
      dimension dns(npnt),ddns(npnt,3),dsdns(npnt)
      dimension valmo(norb*npnt)
      dimension occtyp(nmos),anr(nmos),stp(nmos) 
      dimension vck(npnt),vresp(npnt),erv(norb),orbdnst(npnt*nmos)
      dimension dymo(784),valdo(28*npnt),ia(28),anorm(28)
      dimension skd(nmos,28),oc(28),vip(28),evip(nmos)
      dimension psnomo(norb*(norb+1)/2),srho(npnt,2)
      dimension dnss(npnt),ddnss(npnt,3),dsdnss(npnt),
     +dnst(npnt),ddnst(npnt,3),dsdnst(npnt),sdns(npnt,2),
     +sddns(npnt,3,2),svxc(npnt,2),spksmo(norb*(norb+1)/2,2),
     +spksks(norb*(norb+1)/2,2),vhrma1(norb*(norb+1)/2)
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
      allocatable grdmo(:,:),norbsm(:,:),drho(:,:),am(:,:),
     + bm(:,:),cm(:,:),dsrho(:),vhartr(:),vcond(:),
     +pkd(:,:),apkd(:,:)
c
      ltrian = .false.
      lcens  = .false.
      lvhart = .false.
c
      norbtr=norb*(norb+1)/2
      iorbsm(1:norb)=0
c
      kkk=1
      parin = 10.d0
      epse  = 2.d-4
      epsq  = 1.d-6
c
      scfdpi = 1.d0-scfdmp
c
      npntmx=npnt
      intpnt=npold-1
      npold=1
c
      occorb=2.d0
c
      allocate(grdmo(norb*npnt,4),stat=ialloc)
      if (ialloc.ne.0) then
        write(6,'(/'' No memory allocated for derivativess'')')
        stop
      endif
c
      do ipnt=1,npnt
        read(99,*) grid(ipnt,1),grid(ipnt,2),grid(ipnt,3),weight(ipnt)
      enddo
      close(99)
c
      open(15,file='vxc.plt')
      rewind(15)
      open(25,file='dns.plt')
      rewind(25)
      open(16,file='ddvxc')
      rewind(16)
      open(17,file='vxcab')
      rewind(17)
c
c** get one electron matrix (on mo basis) from dumpfile
c
c     here should be a call to gamess-us file read
c      call hmatmo(idmp,isoe,hmatx,norb)
      if (lsym) then
c
c** determine symmetry of the orbitals, using the orbitals itself if
c** lintsm, else using the hamiltonian matrix
c
        if (lintsm) then
          call orbsym(iorbsm,nsym,norb,thresh,iint)
        else
          call hsym(iorbsm,nsym,hmatx,norb,thresh)
        endif
        write(6,'(/'' Symmetry ordering of Kohn-Sham vectors'')')
        write(6,'(10(/,20i4))')(iorbsm(i),i=1,norb)
        allocate(norbsm(1:norb,0:nsym),stat=ialloc)
        if (ialloc.ne.0) then
          write(6,'(/'' Unable to allocate symmetry specification'')')
          lsym=.false.
        endif
        norbsm(1:norb,0:nsym)=0
        do isym=1,nsym
          item=0
          do iorb=1,norb
            if (iorbsm(iorb).eq.isym) then
              item=item+1
              norbsm(item,isym)=iorb
            endif
          enddo
          norbsm(isym,0)=item
        enddo
      endif
c
c** read kinetic and electron-nuclear attraction integrals from dumpfile,
c** section 192, and store in the vector vnmat
c
clmm      call rdmat(idmp,norb,tsmat,vnmat,enuc)
clmm      write(6,'(/''  nuclear repulsion energy : '',f16.8,/)')enuc
c
c** read dumpfile : vmopao - molecular orbitals in primitive ao basis
c**                 vmoao  - molecular orbitals in ao basis
c**                 pmo    - mo density matrix in mo basis
c**                 pnomo  - ci density matrix in mo basis
c
      psnomo(1:norbtr)=0.d0
clmm      call rddmp(vmopao,vmoao,pmo,psnomo,norb,nvpr,rnel1,idmp,ismo,
clmm     + isao,4)
clmm      call rddmp(vmopao,vmoao,pmo,pnomo,norb,nvpr,rnel,idmp,ismo,isao,
clmm     + isno)
c      call getdmp(idmp,isno,vmoao,occno,norb)
c            do i=1,norb
c      write(6,'(''MO numbers   : '',i5,(1x,f11.5))') i,occno(i)
c      write(6,'(''vnoao   : '',14(1x,f11.5))')
c     + (vmoao((i-1)*norb+j),j=1,norb)
c      enddo
c      call clcvhr(pnomo,norb,iint,vhrmat)
c
c** set elements of vectors, corresponding to the diagonal elements of a
c** triangular matrix (norb,norb), equal to the occupation of the HF-MOs
c
      pksks(1:norbtr)=0.d0
      if (lrdocc) then
        rksnel=sum(pmo)
        if (abs(rnel-rksnel).gt.1.d-6) then
          write(6,'(/''ERROR; sum of KS occupation : '',f8.2)')rksnel
          write(6,'(''       sum of no occupation : '',f8.2)')rnel
          stop
        endif
        pksks(1:norbtr)=pmo(1:norbtr)
        do imo=norb,1,-1
          jmo=imo*(imo+1)/2
          if (pksks(jmo).gt.eps) then
            nmos=imo
            if (nmos.gt.nmomx) then
              write(6,'(''ERROR; nmos > nmomx'')')
              write(6,'(''  nmos  = '',i3,/,''  nmomx = '',i3)')
     + nmos,nmomx
              stop
            endif
            exit
          endif
        enddo
        occmo(1:nmomx)=0.d0
        do imo=1,nmos
          jmo=imo*(imo+1)/2
          occmo(imo)=pksks(jmo)
        enddo
      else
        rksnel=sum(occmo)
        if (abs(rnel-rksnel).gt.1.d-6) then
          write(6,'(/''ERROR; sum of KS occupation : '',f8.2)')rksnel
          write(6,'(''       sum of no occupation : '',f8.2)')rnel
          stop
        endif
        do imo = 1,nmos
          pksks(imo*(imo+1)/2) = occmo(imo)
        enddo
        rnelo=sum(occmo)
        if (abs(rnelo-rnel).gt.eps) then
          write(6,'(/''ERROR; sum of Kohn-Sham occupation numbers :'',
     + f8.3)') rnelo
          stop
        endif
      endif
      pksmo(1:norbtr)=pksks(1:norbtr)
      rnel=sum(occmo)
      if (nvpr.lt.nmos) nvpr=nmos
      call clcvhr(pksmo,norb,iint,vhrmat)
c
c** dipole : transform KS-density matrix from mo to ao basis
c     
      call tmtdag(pksmo,norb,pksao,norb,vmoao,5.d-1)
      call clcdip(dxyz,pksao,idmp,norb)
      write(6,'(/'' initial dipole moment'')')
      write(6,'(''  dx = '',f14.10,''    dy = '',f14.10, 
     + ''    dz = '',f14.10)')dxyz(1),dxyz(2),dxyz(3)
c    
c      if (lrfun) then
        allocate(drho(npnt,3),stat=ialloc)
        if (ialloc.ne.0) then
          write(6,'(/'' No memory allocated for density gradients'')')
          stop
        endif
        allocate(dsrho(npnt),stat=ialloc)
        if (ialloc.ne.0) then
          write(6,'(/'' No memory allocated for Laplacian of '',
     + ''density'')')
          stop
        endif
c      endif
c
      call clcvls(dnst,ddnst,dsdnst,vnuc,valmo,grdmo,npnt,npntmx,
     + norb,pnomo,vmopao,grid,weight,rnel,atmol4)
      call clcvls(dnss,ddnss,dsdnss,vnuc,valmo,grdmo,npnt,npntmx,
     + norb,psnomo,vmopao,grid,weight,rnel1,atmol4)
      nmoss(1:2)=0
      do i=1,nmomx
      do j=1,2
      soccmo(i,j)=0.d0
      enddo
      enddo
      do i=1,nmos
      soccmo(i,1)=1.0d0
      nmoss(1)=nmoss(1)+1
      if (occmo(i).lt.1.5d0) then
      soccmo(i,2)=0.d0
      else
      soccmo(i,2)=1.0d0
      nmoss(2)=nmoss(2)+1
      endif
      enddo
      do isp=1,2
      do i=1,norbtr
      spksks(i,isp)=0.d0
      enddo
      do i=1,nmoss(isp)
      spksks(i*(i+1)/2,isp)=soccmo(i,isp)
      enddo
      do i=1,norbtr
      spksmo(i,isp)=spksks(i,isp)
      enddo
      enddo 
      open(76,file='vnuc.dat')
      rewind(76)
      do ip=1,npnt
        write(76,*) vnuc(ip)
      enddo
c
      nroot=28
      na1=10
      ma1=28
      i=1
      open(88,file='dag.dat')
      do id=1,nroot
      do irow=1,na1
      if(irow.eq.10) then
      read(88,*) dymo(i)
      i=i+1
      else
      read (88,*) dymo(i),dymo(i+1),dymo(i+2)
      i=i+3
      endif
      enddo
      enddo
      close(88)
      open(77,file='dymo.dat')
      rewind(77)
      do id=1,nroot
      do irow=1,28
      i=(id-1)*ma1+irow
      write(77,*) dymo(i)
      enddo
      enddo
      close(77)
      open(66,file='iag.dat')
      rewind(66)
      do i=1,ma1
      read(66,*) ia(i)
      enddo
      close(66) 
      open(89,file='ocag.dat')
      rewind(89)
      do i=1,nroot
      read(89,*) oc(i)
      enddo
      close(89)
      do i=1,npnt
      do j=1,nroot
      ij=(i-1)*nroot+j
      valdo(ij)=0.0d0
      enddo
      enddo
      do ipnt=1,npnt
      k=(ipnt-1)*nroot
      j=(ipnt-1)*norb
      do id=1,nroot
      ic=(id-1)*ma1
      do im=1,ma1
      imo=ia(im)
      valdo(k+id)=valdo(k+id)+valmo(j+imo)*dymo(ic+im)
      enddo
      enddo
      enddo
      do i=1,nroot
      anorm(i)=0.0d0
      enddo
      do ipnt=1,npnt
      k=(ipnt-1)*nroot
      do id=1,nroot
      anorm(id)=anorm(id)+valdo(k+id)*valdo(k+id)*weight(ipnt)
      enddo
      enddo
      write(6,'(/"Dyson norm :",10(1x,f11.6))')
     +(anorm(i),i=1,nroot) 
c
      open(11,file='rhcc.dat')
      rewind(11)
      do ipnt=npold,npnt
        write(11,*) dns(ipnt),ddns(ipnt,1),ddns(ipnt,2),ddns(ipnt,3),
     + dsdns(ipnt)
      enddo
      close(11)
c
      do ipnt=1,npnt
c      dnss(ipnt)=0.d0
      sdns(ipnt,1)=0.5d0*(dnst(ipnt)+dnss(ipnt))
      sdns(ipnt,2)=0.5d0*(dnst(ipnt)-dnss(ipnt))
      do i=1,3
c      ddnss(ipnt,i)=0.d0
      sddns(ipnt,i,1)=0.5d0*(ddnst(ipnt,i)+ddnss(ipnt,i))
      sddns(ipnt,i,2)=0.5d0*(ddnst(ipnt,i)-ddnss(ipnt,i))
      enddo
      enddo
      open(90,file='vntial.dat',status='old',err=222)
      write(6,'(/'' initial potential read from file''/)')
      rewind(90)
      do ip=1,npnt
        read(90,*) vxc(ip)
      enddo
      close(90)
      goto 223
  222 do i=1,2
      do ipnt=1,npnt
      dns(ipnt)=sdns(ipnt,i)
      do j=1,3
      ddns(ipnt,j)=sddns(ipnt,j,i)
      enddo
      enddo 
      if (i.eq.2) alpha=1.336d0
      call sntial(npnt,dns,ddns,vxc,alpha,beta,gamma)
      do ipnt=1,npnt
      svxc(ipnt,i)=vxc(ipnt)
      enddo
      enddo
  223 continue
c
      open(36,file='vhartr.dat',status='old',err=333)
      rewind(36)
      allocate(vhartr(npnt),stat=ialloc)
      if (ialloc.ne.0) then
        write(6,'(/'' No memory allocated for vhartr'')')
      else
        write(6,'(/'' Allocate memeory for vhartr and read from '',
     + ''file'')')
        do ip=1,npnt
          read(36,*) vhartr(ip)
        enddo
        lvhart=.true.   
        if (lsym) then
          call vksmts(vhrmat,vhartr,valmo,weight,norbsm,nsym,intpnt,
     + norb)
        else
          call vksmat(vhrmat,intpnt,norb,vhartr,valmo,weight)
        endif
      endif
      close(36) 
  333 continue
      open(57,file='vcond.dat',status='old',err=332)
      rewind(57)
      allocate(vcond(npnt),stat=ialloc)
      if (ialloc.ne.0) then
        write(6,'(/'' No memory allocated for vcond'')')
      else
        write(6,'(/'' Allocate memory for vcond and read from '',
     + ''file'')')
        do ip=1,npnt
          read(57,*) vcond(ip)
        enddo
        close(57)
      endif
  332 continue
c
c** In the following cycles, pksmo will converge to the KS-density matrix
c
      itrens=itrx+1
      ippr=1
c
      do 100 itr=0,itrx
        itrens=itrens+1
      if (itr.ne.0) then
      do isp=1,2
      do i=1,norbtr
      pksmo(i)=spksmo(i,isp)
      enddo
      if (isp.eq.1) then
      call clcvhr(pksmo,norb,iint,vhrmat)
      else
      call clcvhr(pksmo,norb,iint,vhrma1)
      endif
      enddo
      do i=1,norbtr
      vhrmat(i)=vhrmat(i)+vhrma1(i)
      enddo
      endif
      do isp=1,2
      nmos=nmoss(isp)
      do i=1,nmos
      occmo(i)=soccmo(i,isp)
      enddo
      do i=1,norbtr
      pksmo(i)=spksmo(i,isp)
      pksks(i)=spksks(i,isp)
      enddo
c
        write(6,'(//,'' ----- itr '',i4,'' -----'')')itr
c
        tstr = 0.d0
        tester = 0.d0
c
        if (itr.ne.0) then
          if (lrfun) then
            call grdrho(orbdnst,rho,drho,dsrho,ityp,nmos,ndvcr,norb,
     + npnt,npntmx,pksmo,vksmo,occmo,valmo,grdmo,grid,weight,rnel)
            call lrchxc(vxc,rho,orbdnst,drho,dsrho,dns,ddns,dsdns,
     + weight,occmo,dvdmp,dvmax,drhmx,idvmax,idrhmx,df,ityp,intpnt,
     + npnt,npntmx,nmos,ndvcr)
          else
      do ipnt=1,npnt
      vxc(ipnt)=svxc(ipnt,isp)
      rho(ipnt)=srho(ipnt,isp)
      dns(ipnt)=sdns(ipnt,isp)
      enddo      
      call vlchxc(vxc,rho,dns,dvmax,drhmx,idvmax,idrhmx,df,
     + crrmn,crrmx,npnt,intpnt)
      do ipnt=1,npnt
      svxc(ipnt,isp)=vxc(ipnt)
      enddo
      endif
      else
      do ipnt=1,npnt
      vxc(ipnt)=svxc(ipnt,isp)
      dns(ipnt)=sdns(ipnt,isp)
      rho(ipnt)=dns(ipnt)
      enddo 
        endif
c
c** calculate matrix elements vxcmat in mo basis
c
        if (lsym) then
          call vksmts(vxcmat,vxc,valmo,weight,norbsm,nsym,intpnt,norb)
        else
          call vksmat(vxcmat,intpnt,norb,vxc,valmo,weight)
        endif
c
c        if (.not.lvhart) then
c          call clcvhr(pksmo,norb,iint,vhrmat)
c        endif
c
c** add elements to h-matrix
c
        hksmat(1:norbtr)=hmatx(1:norbtr)+vhrmat(1:norbtr)+
     + vxcmat(1:norbtr)
c
c** expand matrix hksmat in triangular form to full n*n matrix
c
        call sqrmat(hksmat,norb)
c
c** diagonalize matrix hksmat and return the eigenvalues on the main 
c** diagonal of hksmat. vksmo returns the eigenvectors in mo basis.
c
        call diaglz(hksmat,vksmo,norb)
        do iorb=1,norb
          eorb(iorb)=hksmat((iorb-1)*norb+iorb)
        enddo
        nmos=nmoss(isp)
c
        if (itr.ne.0) qorbld(1:nmos)=qorb(1:nmos)       
c
        do imo=1,nmos
          iorder(imo)=imo
        enddo
        do imo=1,nmos-1
          emin=eorb(imo)
          do jmo=imo+1,nmos
            if (eorb(jmo).lt.emin) then
              eorb(imo)=eorb(jmo)
              eorb(jmo)=emin
              emin=eorb(imo)
              itmp=iorder(imo)
              iorder(imo)=iorder(jmo)
              iorder(jmo)=itmp
            endif
          enddo
        enddo
        do imo=nmos,1,-1
          if (occmo(iorder(imo)).gt.eps) then
            homo=eorb(imo)
            exit
          endif
        enddo
        do imo=nmos+1,norb
          if ((eorb(imo)-homo).lt.eps) then
            write(6,'(''WARNING; Spectrum is incorrect '',f8.4)')
     + eorb(imo)
            if (nmos.lt.nmomx) then
              insmo=1
              enew=eorb(imo)
              do jmo=nmos,1,-1
                if (enew.gt.eorb(jmo)) then
                  insmo=jmo+1
                  exit
                endif
              enddo
              do jmo=nmos,insmo,-1
                eorb(jmo+1)=eorb(jmo)
                iorder(jmo+1)=iorder(jmo)
              enddo
              eorb(insmo)=enew
              iorder(insmo)=imo
              nmos=nmos+1
            else
              write(6,'(/''ERROR; nmos > nmomx'')')
              write(6,'(''  nmos  = '',i3,/,''  nmomx = '',i3)')
     + nmos,nmomx
              stop
            endif
          endif
        enddo
        do imo=1,nmos
          jmo=iorder(imo)
          qorb(imo)=occmo(jmo)
        enddo
        do imo=nmos,1,-1
          if (qorb(imo).gt.eps) then
            nmopr=imo
            exit
          endif
        enddo
c
        elumo=1.d3
        do iorb=nmos+1,norb
          if (eorb(iorb).lt.elumo) elumo=eorb(iorb)
        enddo
c
        write(6,'(/)')
        if (lsym) then
c
c** delete symmetry contaminating elements from Kohn-Sham vectors
c
          smermx=-1.d0
          do i=1,norb
            xnorm=0.d0
            ii=(i-1)*norb
            do j=1,norb
              ij=ii+j
              if (iorbsm(i).eq.iorbsm(j)) then
                xnorm=xnorm+vksmo(ij)*vksmo(ij)
              else
                vksmo(ij)=0.d0
              endif
            enddo
            smerr=1.d0-xnorm
            if (smerr.gt.smermx) smermx=smerr
            if (smerr.gt.1.d-6) then
              write(6,'(''WARNING Orbital '',i3,'' has '',g8.2,
     + ''% admixture of wrong symmetry'')')i,100*abs(xnorm-1.d0)
            endif
            if (xnorm.lt.1.d-2) then
              write(6,'(''ERROR in normalization. Orbital '',
     + i3,'' has norm '',f5.2)')i,xnorm
              stop
            endif
c
c** normalize KS vector
c
            xnorm = 1.d0/sqrt(xnorm)
            do j=1,norb
              vksmo(ii+j)=vksmo(ii+j)*xnorm
            enddo
          enddo
          if (abs(smermx).gt.1.d-6) 
     + write(6,'(''  largest symmetry error : '',g12.4)')smermx
          write(6,'(''  symm  : '',10(4x,i3,5x))')
     + (iorbsm(iorder(i)),i=1,nmopr)
        endif
c
c** check the KS orbital density
c
        chkval(1:nmomx)=0.d0
        do ipnt = 1,intpnt
          k=(ipnt-1)*norb
          do imo=1,nmopr
            tmp=0.d0
            i=iorder(imo)
            ii=(i-1)*norb
            do j=1,norb
              tmp=tmp+vksmo(ii+j)*valmo(k+j)
            enddo
            chkval(imo)=chkval(imo)+tmp*tmp*weight(ipnt)
          enddo
        enddo
c
        write(6,'(''  integ : '',10(1x,g11.4))')
     + (chkval(i)-1.d0,i=1,nmopr)
        write(6,'(''  occ   : '',10(1x,f11.5))')
     + (qorb(i),i=1,nmopr)
        write(6,'(''  ev    : '',10(1x,f11.5))')
     + (eorb(i),i=1,nmopr)
c
c** calculate kinetic energy from the one electron energy of the Hartree-Fock 
c** h-matrix (h=T+Vnuc) minus the electron-nuclear atraction energy (Vnuc)
c
        c = 0.5d0
c
c** transform matrix KS-density pksks from KS-orbital basis to HF-MO basis by
c** pvksmo = vksmo pksks vksmo+; pvksmo(kl)=vksmo(ki)*pksks(ij)*vksmo(lj) sum 
c** over i,j. pksks was already initialized with the occupation number input
c
        call tmtdag(pksks,nmos,pvksmo,norb,vksmo,c)
c
c** transform KS-density matrix from mo to ao basis
c
        call tmtdag(pvksmo,norb,pksao,norb,vmoao,c)
c
        evnks=0.d0
        etsks=0.d0
        do iorb=1,norbtr
          evnks=evnks+vnmat(iorb)*pksao(iorb)
          etsks=etsks+tsmat(iorb)*pksao(iorb)
        enddo
        call erep(pvksmo,norb,iint,ehks,exks)
c
c** calculate KS density
c
        do ipnt=1,npnt
          m=(ipnt-1)*norb+1
          pks=valmat(pvksmo,norb,valmo(m),ltrian)
          tstr=tstr+abs(pks-rho(ipnt))*weight(ipnt)
          rho(ipnt)=pks
        enddo
c
c** calculate tester (difference in density)
c
        do ipnt=1,intpnt
          tester=tester+abs(dns(ipnt)-rho(ipnt))*weight(ipnt)
        enddo
c
        write(6,'(/''  error :'',f14.8,'' ('',g14.8,'')'')')tester,tstr
        write(6,'(''  ts    : '',f12.6)')etsks
        write(6,'(''  en    : '',f12.6)')evnks
        write(6,'(''  eh    : '',f12.6)')ehks
        write(6,'(''  ex    : '',f12.6)')exks
        write(6,'(''  total : '',f12.6)')etsks+evnks+ehks+exks+enuc
        write(6,'(/''  elumo : '',f8.4)')elumo
c
c** dipole
c
        call clcdip(dxyz,pksao,idmp,norb)
        write(6,'(/''  dx = '',f14.10,''    dy = '',f14.10,
     + ''    dz = '',f14.10)')dxyz(1),dxyz(2),dxyz(3)
c
        it=1
        ityp(1:nmos)=0
        do imo=1,nmos
          occimo=occmo(imo)
          if (occimo.gt.epsq) then
            if (ityp(imo).eq.0) then
              ityp(imo)=it
              eimo=hksmat((imo-1)*norb+imo)
              do jmo=imo+1,nmos
                occjmo=occmo(jmo)
                if (occjmo.gt.epsq) then
                  ejmo=hksmat((jmo-1)*norb+jmo)
                  if ((abs(eimo-ejmo).lt.epse).and.
     + (abs(occimo-occjmo).lt.epsq)) ityp(jmo)=ityp(imo)
                endif
              enddo
              it=it+1
            endif
          endif
        enddo
        ndvcr=it-1
c
       write(6,'(/''  stype : '',10(4x,i3,5x))')
     + (ityp(iorder(i)),i=1,nmopr)
c
c** determine possible ensemble
c
        if ((itr.eq.ibcens).and.(kpens.ne.0)) then
          write(6,'(/'' ** start ensdet option **'')')
          itrens=iscens
          lcens=.true.
        endif
c
        if (lcens.and.(itrens.eq.iscens)) then
c
c** calculate new occupation pattern
c
          if (itr.ne.0) then 
            lincpr=.true.
            do imo=1,nmos
              if (abs(qorb(imo)-qorbld(imo)).gt.epsq) then
                lincpr=.false.
                exit
              endif
            enddo
            if (lincpr) then
              parin=fparin*parin
            else
              parin=parin/fparin
            endif
            write(6,'('' parin  :'',f8.3)')parin
          endif
c              
          itrens=0
          sumi=sum(qorb)
c
c** determine fermi energy
c
          efermi = eorb(1)
          ifermi = 1
          ilumo = -1
          elumo = 0.d0
          do iorb = 1,nmos
            if ((qorb(iorb).gt.eps).and.
     + (eorb(iorb).gt.efermi)) then
              efermi = eorb(iorb)
              ifermi = iorb
            endif
            if ((abs(qorb(iorb)-occorb).gt.epsq).and.
     + (eorb(iorb).lt.elumo)) then
              elumo = eorb(iorb)
              ilumo = iorb
            endif
          enddo
c
          gap = efermi-elumo
          if ((gap.gt.epse)) then
c
c** charge transfer
c
            nq = ifermi-ilumo+1
            if (kpens.lt.0) kkk=nq
            call frcctr(nq,eorb(ilumo),qorb(ilumo),occorb,dqmax,parin,
     + epse,epsq,kkk)
            sumo=sum(qorb)
            if (abs(sumi-sumo).gt.eps) then
              write(6,'(''ERROR; number of electrons destroyed'')')
              write(6,'(''  in  : '',f6.2,/,''  out : '',f6.2)')
     + sumi,sumo
              stop
            endif
            pksks(1:norbtr)=0.d0
            do imo=1,nmos
              jmo=iorder(imo)
              occmo(jmo)=qorb(imo)
              pksks(jmo*(jmo+1)/2)=qorb(imo)
            enddo
            call tmtdag(pksks,nmos,pvksmo,norb,vksmo,c)
            call tmtdag(pvksmo,norb,pksao,norb,vmoao,c)
          endif
        endif
c
        pksmo(1:norbtr)=scfdpi*pksmo(1:norbtr)+scfdmp*pvksmo(1:norbtr)
        if (tstr.lt.tstthr) goto 110
c
c** store some intermediate results
c
        if (itr.eq.info(ippr)) then
          write(15,'('' vxc at iteration '',i2)')itr
          write(25,'('' dns at iteration '',i2)')itr
          do ip=npold,npnt
            write(15,'(f12.6)') vxc(ip)
            write(25,'(f12.6)') dns(ip)
          enddo
          if (ippr.lt.nppr) ippr=ippr+1
          write(15,'(/)')
          write(25,'(/)')
        endif
      do ipnt=1,npnt
      srho(ipnt,isp)=rho(ipnt)
      enddo       
      do i=1,norbtr
      spksmo(i,isp)=pksmo(i)
      spksks(i,isp)=pksks(i)
      enddo
      do i=1,nmos
      soccmo(i,isp)=occmo(i)
      enddo
      enddo
c
  100 continue
c
c** end of cycle
c
      write(6,'(/,'' ** no convergence **'')')
      goto 120
  110 write(6,'(/,'' ** converged **'')')
  120 continue
      do ip=npold,npnt
      write(16,'(3(1x,f12.6))')
     * (sdns(ip,1)-sdns(ip,2)),(svxc(ip,1)-svxc(ip,2)),grid(ip,3)
      write(17,'(4(1x,f12.6))')
     * (srho(ip,1)-srho(ip,2)),svxc(ip,1),svxc(ip,2),grid(ip,3)
      enddo
c
c      call clvkin(norb,npold,npnt,pnomo,pvksmo,rho,dns,ddns,grdmo
c     +,vck)
c            call grdrho(orbdnst,rho,drho,dsrho,ityp,nmos,ndvcr,norb,
c     + npnt,npntmx,pksmo,vksmo,occmo,valmo,grdmo,grid,weight,rnel)
c      allocate(am(ndvcr,ndvcr),stat=ialloc)
c      if (ialloc.ne.0) then
c      write(6,'(/'' No memory allocated for steps'')')
c      endif
c      allocate(bm(ndvcr,ndvcr),stat=ialloc)
c      if (ialloc.ne.0) then
c      write(6,'(/'' No memory allocated for steps'')')
c      endif
c      allocate(cm(ndvcr,ndvcr),stat=ialloc)
c      if (ialloc.ne.0) then
c      write(6,'(/'' No memory allocated for steps'')')
c      endif
c      do i=1,nmos
c      occtyp(i)=0.d0
c      enddo
c      do i=1,nmos
c      l=ityp(i)
c      if (l.ne.0) then
c      occtyp(l)=occtyp(l)+occmo(i)
c      endif
c      enddo
c      do i=1,ndvcr
c      anr(i)=0.d0
c      erv(i)=0.d0
c      goto 210
c      do ip=1,npnt
c      knm=(ip-1)*ndvcr
c      anr(i)=anr(i)+orbdnst(knm+i)*weight(ip)
c      vresp(ip)=vxc(ip)-vcond(ip)+vhartr(ip)-vck(ip)
c      erv(i)=erv(i)+orbdnst(knm+i)*weight(ip)
c     +*(vxc(ip)-vcond(ip)+vhartr(ip)-vck(ip))
c      enddo
c  210 continue
c      enddo
c      write(6,'(/''  anr  :'',10(1x,f11.6))')
c     +(anr(i),i=1,ndvcr)
c      write(6,'(/''  erv  :'',10(1x,f11.6))')
c     +(erv(i),i=1,ndvcr)
c      write(6,'(/''  occtyp  :'',10(1x,f11.6))')
c     +(occtyp(i),i=1,ndvcr)
c      do i=1,ndvcr
c      do j=1,ndvcr
c      am(i,j)=0.d0
c      enddo
c      enddo
c      do i=1,ndvcr
c      do j=1,i
c      do ip=1,npnt
c      knm=(ip-1)*ndvcr
c      am(i,j)=am(i,j)+orbdnst(knm+i)*orbdnst(knm+j)
c     +*weight(ip)/rho(ip)
c      enddo
c      enddo
c      enddo
c      do i=1,ndvcr
c      do j=1,i-1
c      am(j,i)=am(i,j)
c      enddo
c      enddo
c      do i=1,ndvcr
c      do j=1,ndvcr
c      bm(i,j)=am(i,j)
c      cm(i,j)=0.d0
c      enddo
c      enddo
c      tol=1.d-20
c      call minvr(am,tol,det,ier,ndvcr)
c      if(ier.eq.1) then
c      write(6,'(/,''ERROR ; UNABLE TO INNVERT M MATRIX'')')
c      endif
c      do i=1,ndvcr
c      write(6,'(/''  M  : '',10(1x,f20.6))')
c     +(bm(i,j),j=1,ndvcr)
c      enddo
c      do i=1,ndvcr
c      write(6,'(/''  M-1  : '',10(1x,f20.6))')
c     +(am(i,j),j=1,ndvcr)
c      enddo
c      do i=1,ndvcr
c      do j=1,ndvcr
c      do k=1,ndvcr
c      cm(i,j)=cm(i,j)+am(i,k)*bm(k,j)
c      enddo
c      enddo
c      enddo
c      do i=1,ndvcr
c      write(6,'(/''  (M-1)*M  : '',10(1x,f11.6))')
c     +(cm(i,j),j=1,ndvcr)
c      enddo
c      do i=1,ndvcr
c      stp(i)=0.0d0
c      do j=1,ndvcr
c      stp(i)=stp(i)+erv(j)*am(i,j)
c      enddo
c      enddo
c      write(6,'(/'' stp : '',10(1x,f11.6))')
c     +(stp(i),i=1,ndvcr)
c      goto 211
c      do i=1,nmos
c      do j=1,nroot
c      skd(i,j)=0.0d0
c      enddo
c      enddo
c      do ipnt=1,npnt
c      k=(ipnt-1)*nroot
c      j=(ipnt-1)*norb
c      do ks=1,nmos
c      tmp=0.0d0
c      m=(ks-1)*norb
c      do nor=1,norb
c      tmp=tmp+vksmo(m+nor)*valmo(j+nor)
c      enddo
c      do id=1,nroot 
c      skd(ks,id)=skd(ks,id)+tmp*valdo(k+id)*weight(ipnt)
c      enddo
c      enddo
c      enddo
c      open(67,file='skd.dat')
c      rewind(67)
c      do i=1,nmos
c      do j=1,nroot
c      write(67,*) skd(i,j)
c      enddo
c      enddo
c      close(67) 
c      allocate(pkd(ndvcr,nroot),stat=ialloc)
c      if (ialloc.ne.0) then
c      write(6,'(/" No memory allocated for P")')
c      endif
c      allocate(apkd(ndvcr,nroot),stat=ialloc)
c      if (ialloc.ne.0) then
c      write(6,'(/" No memory allocated for (M-1)*P")') 
c      endif
c      do i=1,ndvcr
c      do j=1,nroot
c      pkd(i,j)=0.0d0
c      apkd(i,j)=0.0d0
c      enddo
c      enddo
c      do ip=1,npnt
c      k=(ip-1)*nroot
c      j=(ip-1)*ndvcr
c      do ks=1,ndvcr
c      do id=1,nroot
c      pkd(ks,id)=pkd(ks,id)+orbdnst(j+ks)*valdo(k+id)
c     +*oc(id)*valdo(k+id)*2.0d0*weight(ip)/rho(ip)
c      enddo
c      enddo
c      enddo
c      do i=1,ndvcr
c      do j=1,nroot
c      do k=1,ndvcr
c      apkd(i,j)=apkd(i,j)+am(i,k)*pkd(k,j)
c      enddo
c      enddo
c      enddo
c      open(68,file='pkd.dat')
c      rewind(68)
c      do i=1,ndvcr
c      do j=1,nroot
c      write(68,*) pkd(i,j)
c      enddo
c      enddo
c      close(68)
c      open(69,file='apkd.dat')
c      rewind(69)
c      do i=1,ndvcr
c      do j=1,nroot
c      write(69,*) apkd(i,j)
c      enddo
c      enddo
c      close(69)
c      open(91,file='vipag.dat')
c      rewind(91)
c      do i=1,nroot
c      read(91,*) vip(i)
c      enddo
c      close(91)
c      do i=1,ndvcr
c      evip(i)=0.0d0
c      do j=1,nroot
c      evip(i)=evip(i)-apkd(i,j)*vip(j)
c      enddo
c      enddo
c      write(6,'(/'' evip : '',10(1x,f11.6))')
c     +(evip(i),i=1,ndvcr)
c      temp=0.0d0
c      do i=1,nroot
c      temp=temp+oc(i)
c      enddo
c      write(6,'(/'' N/2 Dyson : '',10(1x,f11.6))') temp
c  211 continue
c
c** write data to the files vxc.dat, dnst.dat
c
      write(15,'('' vxc at iteration '',i2)')itr
      write(25,'('' dns at iteration '',i2)')itr
      do ip=npold,npnt
        write(15,'(f12.6)') vxc(ip)
        write(25,'(f12.6)') dns(ip)
      enddo
      close(15)
      close(25)
c
      open(16,file='vxc.dat')
      rewind(16)
      open(26,file='dns.dat')
      rewind(26)
      do ip=1,npnt
        write(16,*) vxc(ip)
        write(26,*) grid(ip,3),vresp(ip)
      enddo
      close(16)
      close(26)
c
      call wrtorb(vksmo,occmo,ityp,ndvcr,nmos,norb,npnt,npold,valmo)
c
      if (lvhart) then
        print*,'reconstruct the Kohn-Sham potential from the orbitals'
        call clcpot(norb,nmos,npnt,vksmo,hksmat,occmo,vhartr,vnuc,fxyz,
     + valmo,grdmo(1,4),grid)
      endif
c
c** store ks-vectors and occupation on dumpfile.
c
      if (isks.gt.0) then
        call wrtkso(vksmo,vmoao,hksmat,occmo,nmos,norb,idmp,isks)
      endif
c
c** use hksmat as scratch to get the right order
c
      hksmat(1:norb*norb)=0.d0
      do imo=1,nmos
        jmo=iorder(imo)
        iimo=(imo-1)*norb
        jjmo=(jmo-1)*norb
        hksmat(iimo+1:iimo+norb)=vksmo(jjmo+1:jjmo+norb)
      enddo
      if (nvpr.gt.nmos) then
        nstrt=nmos*norb
        nfinl=nvpr*norb
        hksmat(nstrt+1:nfinl)=vksmo(nstrt+1:nfinl)
      endif
clmm      call prvec(hksmat,norb,nvpr,'Kohn-Sham orbitals in mo basis')
clmm      call wmtrx('####  ks-density matrix in mo basis  ####',
clmm     + pvksmo,norb,1.d-2)
c
c** calculate one electron energy, kinetic energy and electron nuclear 
c** contribution to one electron energy
c
      eonks=0.d0
      eonmo=0.d0
      eonno=0.d0
      etsmo=0.d0
      etsno=0.d0
      evnmo=0.d0
      evnno=0.d0
c
c** first transform density from mo basis to ao basis.
c
      call tmtdag(pmo,norb,pmoao,norb,vmoao,c)
      call tmtdag(pnomo,norb,pnoao,norb,vmoao,c)
c
      do i=1,norbtr
        eonks=eonks+hmatx(i)*pvksmo(i)
        eonmo=eonmo+hmatx(i)*pmo(i)
        eonno=eonno+hmatx(i)*pnomo(i)
        etsmo=etsmo+tsmat(i)*pmoao(i)
        etsno=etsno+tsmat(i)*pnoao(i)
        evnmo=evnmo+vnmat(i)*pmoao(i)
        evnno=evnno+vnmat(i)*pnoao(i)
      enddo
c
c** calculate the Coulomb and Exchange energy associated with the 
c** diagonal one-particle density matrix
c
      call erep(pmo,norb,iint,ehmo,exmo)
      call erep(pnomo,norb,iint,ehno,exno)
c
      write(6,'(///,'' HF kinetic energy                      :'',
     +f16.10)')etsmo
      write(6,'('' HF electron-nuclear attraction energy  :'',
     +f16.10)')evnmo
      write(6,'('' HF coulomb repulsion energy            :''
     +,f16.10)')ehmo
      write(6,'('' HF exchange energy                     :''
     +,f16.10)')exmo
      write(6,'('' HF total electronic energy             :''
     +,f16.10)')eonmo+exmo+ehmo+enuc
c
      write(6,'(/,'' CI kinetic energy                      :'',
     +f16.10)')etsno
      write(6,'('' CI electron-nuclear attraction energy  :'',
     +f16.10)')evnno
      write(6,'('' CI coulomb repulsion energy            :''
     +,f16.10)')ehno
c
c** calculate Kohn-Sham kinetic energy which is the
c** one-electron energy minus the electron nuclear
c** attraction energy
c
      write(6,'(/,'' KS kinetic energy                      :'',
     +f16.10)') etsks
      write(6,'('' KS electron-nuclear attraction energy  :'',
     +f16.10)')evnks
      write(6,'('' KS coulomb repulsion energy            :''
     +,f16.10)')ehks
      write(6,'('' KS exchange energy                     :''
     +,f16.10)')exks
      write(6,'('' KS energy expectation value            :''
     +,f16.10)')eonks+exks+ehks+enuc
      do imo=1,nmos
        if ((abs(occmo(imo)-occorb).gt.eps).and.
     + (abs(occmo(imo)).gt.eps)) nfomos=nfomos+1
      enddo
      if ((nfomos.gt.0).and.(nmos.gt.1)) then
        if (nfomos.eq.1) then
c
          exks=0.d0
          pksiks(1:norbtr)=0.d0
          do is=1,2
            do i=1,nmos
              ii=i*(i+1)/2
              if (occmo(i).gt.1.d-6) then
                if (is.eq.1) then
                  pksiks(ii)=1.d0
                else
                  pksiks(ii)=occmo(i)-1.d0     
                endif
              endif
            enddo
c
            write(6,'(/)')
            call wmtrx('####  ks-density matrix in ks basis  ####',
     + pksiks,norb,1.d-2)
c
            call tmtdag(pksiks,norb,pksimo,norb,vksmo,5.d-1)
            call erep(pksimo,norb,iint,rdum,extmp)
            exks=exks+2.d0*extmp
          enddo
c
      write(6,'(/'' KS exchange energy                     :''
     +,f16.10)')exks
      write(6,'('' KS energy expectation value            :''
     +,f16.10)')eonks+exks+ehks+enuc
c
c        else
c          call clcens(vksmo,pksks,occmo,occorb,eonks,norb,nmos,nfomos,
c     + iint)
        endif
      endif
c
c** Calculate orbital kinetic energies
c
      c=5.d-1
      pksiks(1:norbtr)=0.d0
      do imo=1,nmos
        etsi=0.d0
        iimo=imo*(imo+1)/2
        pksiks(iimo)=occmo(imo)
        call tmtdag(pksiks,nmos,pksimo,norb,vksmo,c)
        call tmtdag(pksimo,norb,pksiao,norb,vmoao,c)
        do iorb=1,norbtr
          etsi=etsi+tsmat(iorb)*pksiao(iorb)
        enddo
        tsorb(imo)=etsi
        pksiks(iimo)=0.d0
      enddo
      write(6,'(/''  ev : '',10(1x,f11.5))')
     + (eorb(i),i=1,nmopr)
      write(6,'(''  ts : '',10(1x,f11.5))')
     + (tsorb(iorder(i)),i=1,nmopr)
      write(6,'(///,'' CI energies                      :'',
     +4f16.10)') etsno,evnno,ehno,enuc
      do i=1,norb
      ii=(i-1)*i/2
      pnomoa(ii+i)=0.5d0*(pnomo(ii+i)+psnomo(ii+i))
      pnomob(ii+i)=0.5d0*(pnomo(ii+i)-psnomo(ii+i)) 
      if(i.gt.1) then
      do j=1,i-1
      pnomoa(ii+j)=0.25d0*(pnomo(ii+j)+psnomo(ii+j))
      pnomob(ii+j)=0.25d0*(pnomo(ii+j)-psnomo(ii+j))
      pnomo(ii+j)=0.5d0*pnomo(ii+j)
      psnomo(ii+j)=0.5d0*psnomo(ii+j)
      enddo
      endif
      enddo
c      pnomoa(1:norbtr)=0.5d0*(pnomo(1:norbtr)+psnomo(1:norbtr))
c      pnomob(1:norbtr)=0.5d0*(pnomo(1:norbtr)-psnomo(1:norbtr))
c      pnomob(1:norbtr)=psnomo(1:norbtr)
c      do i=1,norb
c      ii=(i-1)*i/2
c      pnomoa(ii+i)=pnomo(ii+i)
c      if(i.gt.1) then
c      do j=1,i-1
c      pnomoa(ii+j)=0.5d0*pnomo(ii+j)
c      enddo
c      endif
c      enddo
c      do i=1,norb
c      ii=(i-1)*i/2
c      pnomo(ii+i)=psnomo(ii+i)
c      if(i.gt.1) then
c      do j=1,i-1
c      pnomob(ii+j)=0.5d0*psnomo(ii+j)
c      enddo
c      endif
c      enddo
      vnomoa(1:norb*norb)=0.d0
      vnomob(1:norb*norb)=0.d0
      vnomo(1:norb*norb)=0.d0
      vnomos(1:norb*norb)=0.d0
      vnoma(1:norb*norb)=0.d0
      vnomb(1:norb*norb)=0.d0
      vnom(1:norb*norb)=0.d0
      vnoms(1:norb*norb)=0.d0
c      write(6,'(''pnomoa   : '',10(1x,f11.5))')
c     + (pnomoa(i),i=1,norbtr)
c      write(6,'(''pnomob   : '',10(1x,f11.5))')
c     + (pnomob(i),i=1,norbtr)
      call sqrmat(pnomoa,norb)
      call sqrmat(pnomob,norb)
      call sqrmat(pnomo,norb)
      call sqrmat(psnomo,norb)
      call diaglz(pnomoa,vnoma,norb)
      call diaglz(pnomob,vnomb,norb)
      call diaglz(pnomo,vnom,norb)
      call diaglz(psnomo,vnoms,norb)
c      vnomb(1:norb*norb)=vnoma(1:norb*norb) 
c      do i=1,norb
c      write(6,'(''numbers   : '',i5)') i
c      write(6,'(''vnoma   : '',14(1x,f11.5))')
c     + (vnoma((i-1)*norb+j),j=1,norb)
c      enddo
c      do i=1,norb
c      write(6,'(''numbers   : '',i5)') i
c      write(6,'(''vnomb   : '',14(1x,f11.5))')
c     + (vnomb((i-1)*norb+j),j=1,norb)
c      enddo
      do iorb=1,norb
      occna(iorb)=abs(pnomoa((iorb-1)*norb+iorb))
      occnb(iorb)=abs(pnomob((iorb-1)*norb+iorb))
      occn(iorb)=abs(pnomo((iorb-1)*norb+iorb))
      occns(iorb)=psnomo((iorb-1)*norb+iorb)
      enddo
c      write(6,'(''occna   : '',10(1x,f11.7))')
c     + (occna(i),i=1,norb)
c      write(6,'(''occnb   : '',10(1x,f11.7))')
c     + (occnb(i),i=1,norb)
c      write(6,'(''occn  : '',10(1x,f11.7))')
c     + (occn(i),i=1,norb)
      do i=1,norb
      ic=1
      do j=1,norb
      if(i.ne.j) then
      if(occns(i).lt.occns(j)) then
      ic=ic+1
      endif
      if(occns(i).eq.occns(j)) then
      if(i.gt.j) then
      ic=ic+1
      endif
      endif
      endif
      enddo
      inew(i)=ic
      enddo
      do i=1,norb
      is=inew(i)
      occnos(is)=occns(i)
      do j=1,norb
      vnomos((is-1)*norb+j)=vnoms((i-1)*norb+j)
      enddo
      enddo
      do i=1,norb
      ic=1
      do j=1,norb
      if(i.ne.j) then
      if(occna(i).lt.occna(j)) then
      ic=ic+1
      endif
      if(occna(i).eq.occna(j)) then
      if(i.gt.j) then
      ic=ic+1
      endif
      endif
      endif
      enddo
      inew(i)=ic
      enddo
      do i=1,norb
      is=inew(i)
      occnoa(is)=occna(i)
      do j=1,norb
      vnomoa((is-1)*norb+j)=vnoma((i-1)*norb+j)
      enddo
      enddo
      do i=1,norb
      ic=1
      do j=1,norb
      if(i.ne.j) then
      if(occnb(i).lt.occnb(j)) then
      ic=ic+1
      endif
      if(occnb(i).eq.occnb(j)) then
      if(i.gt.j) then
      ic=ic+1
      endif
      endif
      endif
      enddo
      inew(i)=ic
      enddo
      do i=1,norb
      is=inew(i)
      occnob(is)=occnb(i)
      do j=1,norb
      vnomob((is-1)*norb+j)=vnomb((i-1)*norb+j)
      enddo
      enddo
      do i=1,norb
      ic=1
      do j=1,norb
      if(i.ne.j) then
      if(occn(i).lt.occn(j)) then
      ic=ic+1
      endif
      if(occn(i).eq.occn(j)) then 
      if(i.gt.j) then
      ic=ic+1
      endif
      endif
      endif
      enddo 
      inew(i)=ic
      enddo
      do i=1,norb
      is=inew(i)
      occno(is)=occn(i) 
      do j=1,norb
      vnomo((is-1)*norb+j)=vnom((i-1)*norb+j)
      enddo
      enddo
      write(6,'(''occnoa   : '',8(1x,f12.10))')
     + (occnoa(i),i=1,norb)
c      write(6,'(''occnob   : '',10(1x,f11.7))')
c     + (occnob(i),i=1,norb)
c      vnoaoa(1:norb*norb)=0.d0
c      call matml(vnomo,vmoao,vnoaoa,norb) 
      do i=1,norb
      write(6,'(''numbers   : '',i5)') i
      write(6,'(''vnoaoa   : '',14(1x,f11.5))')
     + (vnoaoa((i-1)*norb+j),j=1,norb)
      enddo
c      write(6,'(''numbers   : '',i5)') i
c      write(6,'(''vnomoa   : '',14(1x,f11.5))')
c     + (vnomoa((i-1)*norb+j),j=1,norb)
c      do i=1,norb
c      write(6,'(''numbers   : '',i5)') i
c      write(6,'(''vnomob   : '',14(1x,f11.5))')
c     + (vnomob((i-1)*norb+j),j=1,norb)
c      enddo
c      write(6,'(''occnoa   : '',10(1x,f11.7))')
c     + (occnoa(i),i=1,norb)
c      write(6,'(''occnob   : '',10(1x,f11.7))')
c     + (occnob(i),i=1,norb)
      do i=1,norb
      ii=(i+1)*i/2
      pnono(ii)=sqrt(2.0d0*occno(i))
      enddo
      pnomo(1:norbtr)=0.d0
      extot=0.d0
      call tmtdag(pnono,norb,pnomo,norb,vnomo,c)
      call erep(pnomo,norb,iint,ehdum,extot)
      write(6,'(/,"BB-tot energy:",5f16.10)') aa,bb,ehdum,extot,
     +etsno+evnno+ehno+extot+enuc
      do i=1,norb
      ii=(i+1)*i/2
      if(occnos(i).ge.0.d0) then
      si=1.0d0
      else
      si=-1.0d0
      endif 
      pnons(ii)=occnos(i)
      enddo
      psnomo(1:norbtr)=0.d0
      exs=0.d0
      call tmtdag(pnons,norb,psnomo,norb,vnomos,c)
      call erep(psnomo,norb,iint,ehdum,exs)
c      extot=extot+exs
      write(6,'(/,"BB-tots energy:",5f16.10)') aa,bb,ehdum,extot,
     +etsno+evnno+ehno+extot+enuc
      do i=1,norb
      ii=(i+1)*i/2
      pnona(ii)=sqrt(2.0d0*occnoa(i))
      enddo
      pnomoa(1:norbtr)=0.d0
      exno=0.d0
      call tmtdag(pnona,norb,pnomoa,norb,vnomoa,c)
      call erep(pnomoa,norb,iint,ehdum,exno)
      do i=1,norb
      ii=(i+1)*i/2
      pnonb(ii)=sqrt(2.0d0*occnob(i))
      enddo
      pnomob(1:norbtr)=0.d0
      exi=0.d0
      call tmtdag(pnonb,norb,pnomob,norb,vnomob,c)
      call erep(pnomob,norb,iint,ehdum,exi)
      exno=exno+exi
      write(6,'(/,"BB energy:",7f16.10)') etsno,evnno,ehno,exno,
     +ehno+exno,enuc,etsno+evnno+ehno+exno+enuc
      ii1(1:8)=0
      jj1(1:8)=0
      ll1(1:8)=0
      kk1(1:8)=0
      gaa(1:norb*norb)=0.d0
      gbb(1:norb*norb)=0.d0
      gab(1:norb*norb)=0.d0
      call search(1,iint)
  508 call find(iint)
      call get(gin,nw)
      if (nw.eq.0) goto 509
      call unpack(gijkl,8,ijkl,4*nword)
      do m=1,nword
      if (abs(gin(m)).gt.1.d-12) then
      ii1(1)=ijkl(1,m)
      jj1(1)=ijkl(2,m)
      kk1(1)=ijkl(3,m)
      ll1(1)=ijkl(4,m)
      call nrperm(ii1(1),jj1(1),kk1(1),ll1(1),np)
      do iperm=1,np
      k=0
      do i=1,norb
      ii=(i-1)*norb
      if(k.eq.0) then
      iin=ii+norb
      else
      iin=ii-norb
      endif
      do j=1,norb,2
      jj=(j-1)*norb
      jjn=jj+norb
      gaa(ii+j)=gaa(ii+j)+vnomoa(ii+ii1(iperm))*vnomoa(jj+
     +jj1(iperm))*vnomoa(iin+kk1(iperm))*vnomoa(jjn
     ++ll1(iperm))*gin(m)
      enddo
      if(k.eq.0) then
      k=k+1
      else
      k=k-1
      endif
      enddo
      enddo
      endif
      enddo
      goto 508
  509 eetot=0.d0
      k=0
      do i=1,norb
      ii=(i-1)*norb
      if(k.eq.0) then
      sig=1.0d0
      else
      sig=-1.0d0
      endif
      do j=1,norb,2
      eetot=eetot+sig*sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
      enddo
      if(k.eq.0) then
      k=k+1
      else
      k=k-1
      endif
      enddo
      esum=enuc+eetot+evnno+etsno
      write(6,'(/,"DMFT energy:",5f16.10)') etsno,evnno,eetot,
     +enuc,esum
      kk=100
      do ll=1,norb
      kk=kk+1
      enddo
      write(6,'(/,"DMFT energy:",5f16.10)') etsno,evnno,eetot
c      gaa(ii+i)=gaa(ii+i)+vnomoa(ii+ii1(iperm))*vnomoa(ii+
c     +jj1(iperm))*vnomoa(ii+kk1(iperm))*vnomoa(ii+ll1(iperm))*
c     +gin(m)
c      if(i.lt.nmosa) then
c      do j=1,nmosa
c      j=nmosa
c      if(i.ne.j) then
c      jj=(j-1)*norb
c      gaa(ii+j)=gaa(ii+j)+vnomoa(ii+ii1(iperm))*vnomoa(jj+
c     +jj1(iperm))*vnomoa(ii+kk1(iperm))*vnomoa(jj+ll1(iperm))*
c     +gin(m)
c      endif
c      enddo
c      endif
c      enddo
c      do i=1,norb
c      ii=(i-1)*norb
c      gbb(ii+i)=gbb(ii+i)+vnomob(ii+ii1(iperm))*vnomob(ii+
c     +jj1(iperm))*vnomob(ii+kk1(iperm))*vnomob(ii+ll1(iperm))*
c     +gin(m)
c      gab(ii+i)=gab(ii+i)+vnomo(ii+ii1(iperm))*vnomo(ii+
c     +jj1(iperm))*vnomo(ii+kk1(iperm))*vnomo(ii+ll1(iperm))*
c     +gin(m)
c      if(i.eq.1) then
c      gabs(1)=gabs(1)+vnomos(ii1(iperm))*vnomos(jj1(iperm))*
c    +vnomos(kk1(iperm))*vnomos(ll1(iperm))*gin(m)
c      endif
c      if(i.le.nmosa) then
c      noi=1
c      nof=nmosb
c      else
c      noi=nmosa
c      nof=nmosa
c      endif
c      do j=1,nmosb
c      do j=noi,nof
c      if(i.ne.j) then
c      jj=(j-1)*norb
c      gbb(ii+j)=gbb(ii+j)+vnomob(ii+ii1(iperm))*vnomob(jj+
c     +jj1(iperm))*vnomob(ii+kk1(iperm))*vnomob(jj+ll1(iperm))*
c     +gin(m)
c      gab(ii+j)=gab(ii+j)+vnomo(ii+ii1(iperm))*vnomo(jj+
c     +jj1(iperm))*vnomo(ii+kk1(iperm))*vnomo(jj+ll1(iperm))*
c     +gin(m)
c      endif
c      enddo
c      endif
c      enddo
c      enddo
c      endif
c      enddo
c      goto 508
c  509 exph=0.d0
c      extot=extot+0.25d0*occnos(1)*occnos(1)*gabs(1)
c      write(6,'(/,"BB-tos energy:",5f16.10)') aa,bb,ehdum,extot,
c     +etsno+evnno+ehno+extot+enuc
c      do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,nmosb
c      extot=extot+(1.0d0-1.0d0/sqrt(2.0d0))*
c     +sqrt(occno(i)*occno(j))*gab(ii+j)
c      extot=extot+(sqrt(occno(i)*occno(j))-0.5d0*
c     +occno(i)*occno(j))*gab(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"BB-tcor energy:",5f16.10)') aa,bb,ehdum,extot,
c     +etsno+evnno+ehno+extot+enuc
c      egu=exno
c      do i=1,norb
c      ii=(i-1)*norb
c      egu=egu-(0.5*occnoa(i)*occnoa(i)-0.5d0*occnoa(i))*gaa(ii+i)
c     +-(0.5*occnob(i)*occnob(i)-0.5d0*occnob(i))*gbb(ii+i)
c      enddo
c      write(6,'(/,"GU energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+egu+enuc
c      nvirt=norb-nmosa
c      pnovo(1:nvirt*(nvirt+1)/2)=0.d0
c      vinomo(1:nvirt*norb)=0.d0
c      pnomvo(1:norbtr)=0.d0
c      do i=nmosa+1,norb
c      ivirt=i-nmosa
c      pnovo((ivirt+1)*ivirt/2)=2.d0*sqrt(occno(i))
c      iivir=(ivirt-1)*norb
c      iitot=(i-1)*norb
c      do j=1,norb
c      vinomo(iivir+j)=vnomo(iitot+j)
c      enddo
c      enddo
c      call tmtdag(pnovo,nvirt,pnomvo,norb,vinomo,c)
c      exph=0.d0
c      call erep(pnomvo,norb,iint,ehdum,exph)
c      extot=extot-exph
c            do i=nmosa+1,norb
c      ii=(i-1)*norb
c      extot=extot-occno(i)*gab(ii+i)
c+sqrt(occno(i)*
c     +occno(nmosa))*gab(ii+nmosa)
c      enddo
c      write(6,'(/,"BBC1-tot energy:",5f16.10)') aa,bb,ehdum,extot,
c     +etsno+evnno+ehno+extot+enuc
c      pnova(1:nvirt*(nvirt+1)/2)=0.d0
c      vinoma(1:nvirt*norb)=0.d0
c      pnomva(1:norbtr)=0.d0
c      do i=nmosa+1,norb
c      ivirt=i-nmosa
c      pnova((ivirt+1)*ivirt/2)=2.d0*sqrt(occnoa(i))
c      iivir=(ivirt-1)*norb
c      iitot=(i-1)*norb
c      do j=1,norb
c      vinoma(iivir+j)=vnomoa(iitot+j)
c      enddo
c      enddo
c      call tmtdag(pnova,nvirt,pnomva,norb,vinoma,c)
c      exph=0.d0
c      call erep(pnomva,norb,iint,ehdum,exph)
c      nvirt=norb-nmosb
c      pnovb(1:nvirt*(nvirt+1)/2)=0.d0
c      vinomb(1:nvirt*norb)=0.d0
c      pnomvb(1:norbtr)=0.d0
c      do i=nmosb+1,norb
c      ivirt=i-nmosb
c      pnovb((ivirt+1)*ivirt/2)=2.d0*sqrt(occnob(i))
c      iivir=(ivirt-1)*norb
c      iitot=(i-1)*norb
c      do j=1,norb
c      vinomb(iivir+j)=vnomob(iitot+j)
c      enddo
c      enddo
c      call tmtdag(pnovb,nvirt,pnomvb,norb,vinomb,c)
c      expi=0.d0
c      call erep(pnomvb,norb,iint,ehdum,expi)
c      exno=exno-exph-expi
c      exc=exno
c            do i=nmosa+1,norb
c      ii=(i-1)*norb
c      exc=exc-occnoa(i)*gaa(ii+i)
c      j=nmosa
c      exc=exc+0.5d0*sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
c      enddo
c      do i=nmosb+1,norb
c      ii=(i-1)*norb
c      exc=exc-occnob(i)*gbb(ii+i)
c      enddo
c      write(6,'(/,"BBC1 energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+exc+enuc
c      do i=nmosa+1,norb
c      ii=(i-1)*norb
c      do j=1,nmosa
c      j=nmosa
c      exc=exc+0.5d0*sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"BBC1m energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+exc+enuc
c      emp=exc
c      do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,nmosb
c      emp=emp-(0.5d0*occnob(i)*occnob(j)
c     +-0.5d0*sqrt(occnob(i)*occnob(j)))*gbb(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"BBC1.5 energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+emp+enuc
c      do i=1,nmosa-1
c      ii=(i-1)*norb
c      do j=i+1,nmosa
c      exc=exc+(sqrt(occnoa(i))*sqrt(occnoa(j))-occnoa(i)
c     +*occnoa(j))*gaa(ii+j)
c      enddo
c      enddo
c      if(nmosb.gt.1) then
c      do i=1,nmosb-1
c      ii=(i-1)*norb
c      do j=i+1,nmosb
c      exc=exc+(sqrt(occnob(i))*sqrt(occnob(j))-occnob(i)
c     +*occnob(j))*gbb(ii+j)
c      enddo
c      enddo
c      endif
c      write(6,'(/,"BBC2 energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+exc+enuc
c      emp=exc
c      do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,nmosb
c      emp=emp-(0.5d0*occnob(i)*occnob(j)
c     +-0.5d0*sqrt(occnob(i)*occnob(j)))*gbb(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"BBC2.5 energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+emp+enuc
c      do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,nmosb-1
c      emp=emp-(0.5d0*occnob(i)*occnob(j)
c     +-0.5d0*sqrt(occnob(i)*occnob(j)))*gbb(ii+j)
c      enddo
c      enddo
c      do i=1,norb
c      ii=(i-1)*norb
c      if(i.ne.nmosa) then
c      egu=egu-(0.5*occnoa(i)*occnoa(i)-0.5d0*occnoa(i))*gaa(ii+i)
c     +-(0.5*occnob(i)*occnob(i)-0.5d0*occnob(i))*gbb(ii+i)
c      enddo
c      write(6,'(/,"GU energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+egu+enuc
c      write(6,'(/,"BBC4 energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+emp+enuc
c      coulom=0.d0
c      do i=1,norb
c      ii=(i-1)*norb
c      do j=1,norb
c      coulom=coulom+0.5d0*occnoa(i)*occnoa(j)*gaa(ii+j)+
c     +0.5d0*occnob(i)*occnob(j)*gbb(ii+j)+occnoa(i)*
c     +occnob(j)*gab(ii+j)
c      enddo
c      enddo
c      do i=1,norb
c      do j=1,norb
c      write(6,'(/,"Coulomb energy:",2i5,3f16.10)') i,j,gaa(ii+j),
c     +gbb(ii+j),gab(ii+j)
c      enddo
c      enddo 
c      excbbc=0.d0
c      exc=0.d0
c      do i=1,norb
c      ii=(i-1)*norb
c      do j=1,norb
c      exc=exc-0.5d0*(sqrt(occnoa(i)*occnoa(j))*
c     +gaa(ii+j)+sqrt(occnob(i)*occnob(j))*gbb(ii+j))
c      enddo
c      enddo
c      write(6,'(/,"BB-spin energy:",10f16.10)') exc,
c     +etsno+evnno+ehno+exc+enuc
c      do i=1,norb
c      ii=(i-1)*norb
c      do j=1,norb
c      excbbc=excbbc-0.5d0*occno(i)*gab(ii+i)
c      excbbc=excbbc-0.5d0*(occnoa(i)*gaa(ii+i)+occnob(i)*gbb(ii+i))
c      excbbc=excbbc-0.5d0*(sqrt(occnoa(i)*occnoa(j))*
c     +gaa(ii+j)+sqrt(occnob(i)*occnob(j))*gbb(ii+j))
c      enddo
c      enddo
c      write(6,'(/,"E1 energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+excbbc+enuc
c      egu=exc
c      do i=1,norb
c      ii=(i-1)*norb
c      egu=egu-(0.5*occnoa(i)*occnoa(i)-0.5d0*occnoa(i))*gaa(ii+i)
c     +-(0.5*occnob(i)*occnob(i)-0.5d0*occnob(i))*gbb(ii+i)
c      enddo
c      write(6,'(/,"GU-spin energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+egu+enuc
c      if(nmosa.gt.1) then
c      do i=1,nmosa-1
c      ii=(i-1)*norb
c      do j=i+1,nmosa
c      excbbc=excbbc-0.5d0*sqrt(occno(i)*occno(j))*gab(ii+j)
c      excbbc=excbbc-sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
c      enddo
c      enddo
c      endif
c      write(6,'(/,"E2:",10f16.10)') excbbc
c      if(nmosb.gt.1) then
c      do i=1,nmosb-1
c      ii=(i-1)*norb
c      do j=i+1,nmosb
c      excbbc=excbbc-0.5d0*sqrt(occno(i)*occno(j))*gab(ii+j)
c      excbbc=excbbc-sqrt(occnob(i)*occnob(j))*gbb(ii+j)
c      enddo
c      enddo
c      endif
c      write(6,'(/,"E3:",10f16.10)') excbbc
c      do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=nmosa+1,norb
c      excbbc=excbbc-0.5d0*sqrt(occno(i)*occno(j))*gab(ii+j)
c      excbbc=excbbc-sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
c      enddo
c      enddo
c      do i=nmosa+1,norb
c      ii=(i-1)*norb
c      do j=nmosa+1,norb
c      if(i.ne.j) then
c      exc=exc+sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
c      endif
c
c      enddo
c      enddo
c      do i=nmosb+1,norb
c      ii=(i-1)*norb
c      do j=nmosb+1,norb
c      if(i.ne.j) then
c      exc=exc+sqrt(occnob(i)*occnob(j))*gbb(ii+j)
c      endif
c
c      enddo
c      enddo
c      write(6,'(/,"BBC1-spin energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+exc+enuc
c      ec=exc
c      ek=exc 
c      do i=1,nmosa-1
c      ii=(i-1)*norb
c      do j=i+1,nmosa
c      ec=ec+(sqrt(occnoa(i))*sqrt(occnoa(j))-occnoa(i)
c     +*occnoa(j))*gaa(ii+j)
c      enddo
c      enddo
c      if(nmosb.gt.1) then
c      do i=1,nmosb-1
c      ii=(i-1)*norb
c      do j=i+1,nmosb
c      ec=ec+(sqrt(occnob(i))*sqrt(occnob(j))-occnob(i)
c     +*occnob(j))*gbb(ii+j)
c      enddo
c      enddo
c      endif
c      write(6,'(/,"BBC2-spin energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+ec+enuc
c      do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=nmosa+1,norb
c      ek=ek-(occnoa(i)*occnoa(j)
c     +-sqrt(occnoa(i)*occnoa(j)))*gaa(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"BBC1+ac energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+ek+enuc
c      do i=1,norb
c      ii=(i-1)*norb
c       if(i.ne.1) then
c      if((i.lt.nmosa-3).or.(i.gt.nmosa)) then 
c      if((i.ne.nmosa).and.(i.ne.(nmosa+1))) then
c      write(6,'(/,"BBC1 energy:",10f16.10)') excbbc,
     +etsno+evnno+ehno+exc+enuc
      ec=ec-(0.5*occnoa(i)*occnoa(i)-0.5d0*occnoa(i))*gaa(ii+i)
c      if(i.ne.nmosb) then
c      if((i.ne.nmosb).and.(i.ne.(nmosb+1))) then
c      ec=ec-(0.5*occnob(i)*occnob(i)-0.5d0*occnob(i))*gbb(ii+i)
c      endif
c      endif
c      enddo
c      write(6,'(/,"BBC2+selGU energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+ec+enuc
c      do j=1,nmosb
c      ii=nmosb*norb
c      ec=ec-(occnob(nmosb+1)*occnob(j)
c     +-sqrt(occnob(nmosb+1)*occnob(j)))*gbb(ii+j)
c      enddo
c      write(6,'(/,"BBC3-spin energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+ec+enuc,gbb(norb+1),gbb(2) 
c      do i=1,nmosb
c      ii=(i-1)*norb
c      do j=nmosb+1,nmosa
c      jj=9
c      exc=exc-(0.5d0*occnob(i)*occnob(jj)
c     +-0.5d0*sqrt(occnob(i)*occnob(jj)))*gbb(ii+jj)
c      enddo
c      enddo
c      write(6,'(/,"BBC1+emp. energy:",10f16.10)') excbbc,
c     +etsno+evnno+ehno+exc+enuc
c      excbbc=excbbc-sqrt(occnob(1)*occnob(2))*gbb(2)
c      excbbc=0.d0
c            do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,norb
c      excbbc=excbbc-0.5d0*occno(i)*gab(ii+i)
c      excbbc=excbbc-0.5d0*(occnoa(i)*gaa(ii+i)+occnob(i)*gbb(ii+i))
c      excbbc=excbbc-0.5d0*(sqrt(occnoa(i)*occnoa(j))*
c     +gaa(ii+j)+sqrt(occnob(i)*occnob(j))*gbb(ii+j))
c      enddo
c      enddo
c            do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,nmosa
c      if(j.ne.i) then
c      excbbc=excbbc-0.5d0*occno(i)*gab(ii+i)
c      excbbc=excbbc-sqrt(occnoa(i)*occnoa(j))*gaa(ii+j)
c      excbbc=excbbc-0.5d0*(sqrt(occnoa(i)*occnoa(j))*
c     +gaa(ii+j)+sqrt(occnob(i)*occnob(j))*gbb(ii+j))
c      endif
c      enddo
c      enddo
c            do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=nmosa+1,norb
c      if(j.ne.i) then
c      excbbc=excbbc-0.5d0*occno(i)*gab(ii+i)
c      excbbc=excbbc-(occnoa(i)*occnoa(j))*gaa(ii+j)
c      excbbc=excbbc-0.5d0*(sqrt(occnoa(i)*occnoa(j))*
c     +gaa(ii+j)+sqrt(occnob(i)*occnob(j))*gbb(ii+j))
c      endif
c      enddo
c      enddo
c            do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=1,nmosb
c      excbbc=excbbc-0.5d0*occno(i)*gab(ii+i)
c      excbbc=excbbc-sqrt(occnob(i)*occnob(j))*gbb(ii+j)
c      excbbc=excbbc-0.5d0*(sqrt(occnoa(i)*occnoa(j))*
c     +gaa(ii+j)+sqrt(occnob(i)*occnob(j))*gbb(ii+j))
c      enddo
c      enddo
c            do i=nmosb+1,nmosa
c      ii=(i-1)*norb
c      do j=nmosb+1,norb
c      if(j.ne.i) then
c      excbbc=excbbc-0.5d0*
       return
       end
