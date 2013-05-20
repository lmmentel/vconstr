      subroutine dbrain(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmn,
     + crrmx,thresh,dqmax,dvdmp,scfdmp,kpens,info,nppr,nvpr,npnt,npold,
     + norb,nmos,nmomx,itrx,ibcens,iscens,
     + lsym,lintsm,lrdocc,lrfun,nele)
      use integralsModule
      use hartreePotentialModule
      use coulombAndExchangeModule
      use ioModule
      use commonsModule
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
!      parameter(epsilon = 1.0d-8)
      parameter(fparin = 1.01d0)
c
      logical lsym,lintsm,lrdocc,lrfun,lincpr
      dimension info(nppr)
      dimension occmo(nmomx)
c
      logical ltrian,lcens,lvhart
      dimension fxyz(3),dxyz(3)
      dimension iorder(nmomx),ityp(nmomx),iord(norb)
      dimension qorb(nmomx),qorbld(kpens*nmomx),chkval(nmomx),
     + tsorb(nmomx),mm(6),d1(norb),d2(norb),d3(norb)
      dimension pnomoa(norb*norb),vnomoa(norb*norb),occnoa(norb)
      dimension pksks(norb*(norb+1)/2),pmo(norb*(norb+1)/2),
     + pmoao(norb*(norb+1)/2),pnomo(norb*(norb+1)/2),
     + pnoao(norb*(norb+1)/2),pksmo(norb*(norb+1)/2),
     + pksao(norb*(norb+1)/2),vxcmat(norb*(norb+1)/2),
     + vhrmat(norb*(norb+1)/2),hmatx(norb*(norb+1)/2),
     + vnmat(norb*(norb+1)/2),tsmat(norb*(norb+1)/2),
     + pksiks(norb*(norb+1)/2),pksimo(norb*(norb+1)/2),
     + pksiao(norb*(norb+1)/2),pvksmo(norb*(norb+1)/2)
      dimension vksmo(norb*norb),vmoao(norb*norb),hksmat(norb*norb),
     + vmopao(norb*norb),vksmoo(norb*norb),vnoao(norb*norb),
     +vnomo(norb*norb),scrtc(norb*norb),pnono(norb*(norb+1)/2),
     +vksno(norb*norb),vksao(norb*norb)
      dimension rho(npnt),vxc(npnt),vnuc(npnt)
      dimension eorb(norb),iorbsm(norb),an(norb),eorbo(norb),
     +occno(norb)
      dimension grid(npnt,3),weight(npnt)
      dimension dns(npnt),ddns(npnt,3),dsdns(npnt)
      dimension valmo(norb*npnt)
      dimension occtyp(nmos),anr(nmos),stp(nmos) 
      dimension vck(npnt),vresp(npnt),erv(norb),orbdnst(npnt*nmos)
      dimension dymo(784),valdo(28*npnt),ia(28),anorm(28)
      dimension skd(nmos,28),oc(28),vip(28),evip(nmos)
      dimension psnomo(norb*(norb+1)/2),srho(npnt),fij(norb*norb)
      dimension ii1(8),jj1(8),ll1(8),kk1(8),iz(norb)
      dimension wij(norb),pnovi((norb-nmos)*(norb-nmos+1)/2),
     +vinomo((norb-nmos)*norb),pnomov(norb*(norb+1)/2)
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
      allocatable grdmo(:,:),norbsm(:,:),drho(:,:),am(:,:),
     + bm(:,:),cm(:,:),dsrho(:),vhartr(:),vcond(:),
     +pkd(:,:),apkd(:,:)
c
      write(*,'(/"now debuggin in dbrain"/)') 
      write(*,*) 'occmo = ', occmo
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
c
c** get one electron matrix (on mo basis) from dumpfile
c
clmm..call hmatmo(idmp,isoe,hmatx,norb)
clmm..added a call for reading gamess-us core hamiltonian matrix in MO basis (need transfomation)
      if (lsym) then
clmm..write an information about features that are not implemented yet
        write(*,*) 'Option lsym=.true. not implemented yet'
        stop 'bye bye from -dbrain-'
c
c** determine symmetry of the orbitals, using the orbitals itself if
c** lintsm, else using the hamiltonian matrix
c
        if (lintsm) then
        write(*,*) 'Option lintsm=.true. not implemented yet'
        stop 'bye bye from -dbrain-'
!          call orbsym(iorbsm,nsym,norb,thresh,iint)
!        else
!          call hsym(iorbsm,nsym,hmatx,norb,thresh)
        endif
        do i=1,5
        j=i
        enddo
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
clmm..call rdmat(idmp,norb,tsmat,vnmat,enuc)
clmm..write(6,'(/''  nuclear repulsion energy : '',f16.8,/)')enuc
clmm..
clmm..read kinetic energy electron-nucleus attarction integrals 
clmm..and nuclear repulsion energy from gamess-us dictionary
      call readOneEintegrals(tsmat, vnmat, enuc)
      if (printLevel >= 1) then 
        call matPrint(tsmat, norb,'Kinetic energy integrals in AO')
        call matPrint(vnmat, norb,'Potential energy integrals in AO')
        write(*,'(/"Nuclear repulsion energy = ",f14.10/)') enuc
      endif
c
c** read dumpfile : vmopao - molecular orbitals in primitive ao basis
c**                 vmoao  - molecular orbitals in ao basis
c**                 pmo    - mo density matrix in mo basis
c**                 pnomo  - ci density matrix in mo basis
c
clmm..call rddmp(vmopao,vmoao,pmo,pnomo,norb,nvpr,rnel,idmp,ismo,isao,
clmm..+ isno)
clmm..
clmm..now read the orbitals
      call getOrbitalsAndDensities(vmopao, 
     & vmoao, pmo, pnomo, norb, 
     & rnel, nele)
      if (printLevel >= 2) then 
        call matPrint(reshape(vmopao, (/norb, norb/)), 
     & 'Vmo in primitive ao, vmopao')
      endif
c
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
          if (pksks(jmo).gt.epsilon) then
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
        if (abs(rnelo-rnel).gt.epsilon) then
          write(6,'(/''ERROR; sum of Kohn-Sham occupation numbers :'',
     + f8.3)') rnelo
          stop
        endif
      endif
      pksmo(1:norbtr)=pksks(1:norbtr)
      rnel=sum(occmo)
      if (nvpr.lt.nmos) nvpr=nmos
c
c** dipole : transform KS-density matrix from mo to ao basis
c     
      call tmtdag(pksmo,norb,pksao,norb,vmoao,5.d-1)
clmm      call clcdip(dxyz,pksao,idmp,norb)
clmm      write(6,'(/'' initial dipole moment'')')
clmm      write(6,'(''  dx = '',f14.10,''    dy = '',f14.10, 
clmm     + ''    dz = '',f14.10)')dxyz(1),dxyz(2),dxyz(3)
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
      call clcvls(dns,ddns,dsdns,vnuc,valmo,grdmo,npnt,npntmx,
     + norb,pnomo,vmopao,grid,weight,rnel)
      open(76,file='vnuc.dat')
      rewind(76)
      do ip=1,npnt
        write(76,*) vnuc(ip)
      enddo
c
clmm      nroot=28
clmm      na1=10
clmm      ma1=28
clmm      i=1
clmm      open(88,file='dag.dat')
clmm      do id=1,nroot
clmm      do irow=1,na1
clmm      if(irow.eq.10) then
clmm      read(88,*) dymo(i)
clmm      i=i+1
clmm      else
clmm      read (88,*) dymo(i),dymo(i+1),dymo(i+2)
clmm      i=i+3
clmm      endif
clmm      enddo
clmm      enddo
clmm      close(88)
clmm      open(77,file='dymo.dat')
clmm      rewind(77)
clmm      do id=1,nroot
clmm      do irow=1,28
clmm      i=(id-1)*ma1+irow
clmm      write(77,*) dymo(i)
clmm      enddo
clmm      enddo
clmm      close(77)
clmm      open(66,file='iag.dat')
clmm      rewind(66)
clmm      do i=1,ma1
clmm      read(66,*) ia(i)
clmm      enddo
clmm      close(66) 
clmm      open(89,file='ocag.dat')
clmm      rewind(89)
clmm      do i=1,nroot
clmm      read(89,*) oc(i)
clmm      enddo
clmm      close(89)
clmm      do i=1,npnt
clmm      do j=1,nroot
clmm      ij=(i-1)*nroot+j
clmm      valdo(ij)=0.0d0
clmm      enddo
clmm      enddo
clmm      do ipnt=1,npnt
clmm      k=(ipnt-1)*nroot
clmm      j=(ipnt-1)*norb
clmm      do id=1,nroot
clmm      ic=(id-1)*ma1
clmm      do im=1,ma1
clmm      imo=ia(im)
clmm      valdo(k+id)=valdo(k+id)+valmo(j+imo)*dymo(ic+im)
clmm      enddo
clmm      enddo
clmm      enddo
clmm      do i=1,nroot
clmm      anorm(i)=0.0d0
clmm      enddo
clmm      do ipnt=1,npnt
clmm      k=(ipnt-1)*nroot
clmm      do id=1,nroot
clmm      anorm(id)=anorm(id)+valdo(k+id)*valdo(k+id)*weight(ipnt)
clmm      enddo
clmm      enddo
clmm      write(6,'(/"Dyson norm :",10(1x,f11.6))')
clmm     +(anorm(i),i=1,nroot) 
c
      open(11,file='rhcc.dat')
      rewind(11)
      do ipnt=npold,npnt
        write(11,*) dns(ipnt),ddns(ipnt,1),ddns(ipnt,2),ddns(ipnt,3),
     + dsdns(ipnt)
      enddo
      close(11)
c
      rho(1:npnt)=dns(1:npnt)
      open(90,file='vntial.dat',status='old',err=222)
      write(6,'(/'' initial potential read from file''/)')
      rewind(90)
      do ip=1,npnt
        read(90,*) vxc(ip)
      enddo
      close(90)
      goto 223
  222 call vntial(npnt,dns,ddns,vxc,alpha,beta,gamma)
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
c      do ip=1,npnt
c      vxc(ip)=vcond(ip)-vhartr(ip)
c      enddo
c
c** In the following cycles, pksmo will converge to the KS-density matrix
c
      itrens=itrx+1
      ippr=1
c
      do 100 itr=0,itrx
        itrens=itrens+1
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
            call vlchxc(vxc,rho,dns,dvmax,drhmx,idvmax,idrhmx,df,
     + crrmn,crrmx,npnt,intpnt)
          endif
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
clmmm..replacement for clcvhr is getHartreePot
          call getHartreepot(vhrmat, pksmo, norb) 
            if (printLevel >= 4) then 
              call matPrint(vhrmat, norb, 'Hartree potential matrix')
            endif 
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
          if (occmo(iorder(imo)).gt.epsilon) then
            homo=eorb(imo)
            exit
          endif
        enddo
        do imo=nmos+1,norb
          if ((eorb(imo)-homo).lt.epsilon) then
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
          if (qorb(imo).gt.epsilon) then
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
              write(6,'(''WARNING Orbital '',i3,'' has '',g9.2,
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
        call getCoulombAndExchange(pvksmo, norb, ehks, exks)
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
        write(6,'(/''  error :'',f16.8,'' ('',g16.8,'')'')')tester,tstr
        write(6,'(''  ts    : '',f12.6)')etsks
        write(6,'(''  en    : '',f12.6)')evnks
        write(6,'(''  eh    : '',f12.6)')ehks
        write(6,'(''  ex    : '',f12.6)')exks
        write(6,'(''  total : '',f12.6)')etsks+evnks+ehks+exks+enuc
        write(6,'(/''  elumo : '',f8.4)')elumo
c
c** dipole
c
clmm        call clcdip(dxyz,pksao,idmp,norb)
clmm        write(6,'(/''  dx = '',f14.10,''    dy = '',f14.10,
clmm     + ''    dz = '',f14.10)')dxyz(1),dxyz(2),dxyz(3)
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
            if ((qorb(iorb).gt.epsilon).and.
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
            if (abs(sumi-sumo).gt.epsilon) then
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
c
  100 continue
c
c** end of cycle
c
      write(6,'(/,'' ** no convergence **'')')
      goto 120
  110 write(6,'(/,'' ** converged **'')')
  120 continue
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
      do ip=npold,npnt
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
clmm      if (isks.gt.0) then
clmm        call wrtkso(vksmo,vmoao,hksmat,occmo,nmos,norb,idmp,isks)
clmm      endif
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
      call getCoulombAndExchange(pmo, norb, ehmo, exmo)
      call getCoulombAndExchange(pnomo, norb, ehno, exno)
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
      open(121,file='test')
      write(121,*) (pnomo(i),i=1,norbtr)
      close(121) 
      write(6,'(''pnomo   : '',10(1x,f11.5))')
     + (pnomo(i),i=1,norbtr)
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
clmm..
      write(*,*) 'nmos before -clcens- = ', nmos
      write(*,*) 'epsilon before -clcens- = ', epsilon
      do imo=1,nmos
        if ((abs(occmo(imo)-occorb).gt.epsilon).and.
     + (abs(occmo(imo)).gt.epsilon)) nfomos=nfomos+1
      enddo
clmm..print the nfomos variable
      write(*,*) 'nfomos before -clcens- = ', nfomos
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
            call getCoulombAndExchange(pksimo, norb, rdum, extmp)
            exks=exks+2.d0*extmp
          enddo
c
      write(6,'(/'' KS exchange energy                     :''
     +,f16.10)')exks
      write(6,'('' KS energy expectation value            :''
     +,f16.10)')eonks+exks+ehks+enuc
c
        else
          call clcens(vksmo,pksks,occmo,occorb,eonks,norb,nmos,nfomos)
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
     +(eorb(i),i=1,norb)
c     + (eorb(i),i=1,nmopr)
      write(6,'(''  ts : '',10(1x,f11.5))')
     + (tsorb(iorder(i)),i=1,nmopr)
c
      return
      end
