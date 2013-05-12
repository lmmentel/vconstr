      subroutine dbrain(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmn,
     + crrmx,thresh,dqmax,dvdmp,scfdmp,kpens,info,nppr,nvpr,npnt,npold,
     + norb,nmos,nmomx,itrx,ibcens,iscens,iint,
     + isks,lsym,lintsm,lrdocc,lrfun,atmol4,gdictfile,gintegfile,
     + nele)
      use integralsModule
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.0d-8)
      parameter(fparin = 1.01d0)
c
      logical lsym,lintsm,lrdocc,lrfun,lincpr,atmol4
      character(*) gdictfile, gintegfile
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
clmm..      call new(dictionary, gdictnfile)
clmm..      call readH(hmatx, dictionary)
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
      call readOneEintegrals(tsmat, vnmat, enuc, gDictFile)
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
     & rnel, nele, gDictFile)
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
     + norb,pnomo,vmopao,grid,weight,rnel,atmol4)
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
          call clcvhr(pksmo,norb,iint,vhrmat)
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
        else
          call clcens(vksmo,pksks,occmo,occorb,eonks,norb,nmos,nfomos,
     + iint)
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
      open(36,file='vmoao.dat')
      rewind(36)
      do i=1,norb
      do j=1,norb
      write (36,*) vmoao(norb*(i-1)+j)
      enddo
      enddo
      close(36)
      do i=1,norb
      iord(i)=1
      do j=1,norb
      if(i.ne.j) then
      if(eorb(i).gt.eorb(j)) then
      iord(i)=iord(i)+1
      endif
      endif
      enddo
      ij=iord(i)
      eorbo(ij)=eorb(i)
      enddo
      k=0
      do i=1,norb
      icount=0
      k=k+1
      if(k.gt.norb) goto 506
      do j=k,norb
      if(i.ne.j) then
      if(iord(i).eq.iord(j)) then
      icount=icount+1
      iord(j)=iord(j)+icount
      k=k+1
      ij=iord(j)
      eorbo(ij)=eorb(j)
      endif
      endif
      enddo
  506 continue
      enddo
      do i=1,norb
      ij=iord(i)
      ii0=(ij-1)*norb
      ii=(i-1)*norb
      eorb(i)=eorbo(i)
      do j=1,norb
      vksmoo(ii0+j)=vksmo(ii+j)
      enddo
      enddo
      do i=1,norb
      ij=(i-1)*norb
      do j=1,norb
      vksmo(ij+j)=vksmoo(ij+j)
      enddo
      enddo
      write(6,'(/''  ev : '',10(1x,f11.5))')
     +(eorb(i),i=1,norb)
      goto 555
      mn=10
      anel=2*nmos
      thr=1.d-6
      aa=0.008d0
      bb=0.045d0
      do i=1,norb
      an(i)=0.d0 
      enddo
      efl=0.5d0*(elumo+homo)
      dnell=0.d0
      call calcni(eorb,an,homo,elumo,efl,aa,bb,anel,dnell,norb,mn)
      if(abs(dnell).le.thr) goto 501
      efr=efl+0.001d0
      dnelr=0.d0
      call calcni(eorb,an,homo,elumo,efr,aa,bb,anel,dnelr,norb,mn)
      if(abs(dnelr).le.thr) goto 501
      if((dnell*dnelr).lt.0.d0) goto 502
      if(abs(dnelr).lt.abs(dnell)) then
  503 dnell=dnelr
      efl=efr
      efr=efr+0.001d0
      call calcni(eorb,an,homo,elumo,efr,aa,bb,anel,dnelr,norb,mn)
      if(abs(dnelr).le.thr) goto 501
      if((dnell*dnelr).lt.0.d0) goto 502
      goto 503
      else
  504 dnelr=dnell
      efr=efl
      efl=efl-0.001d0
      call calcni(eorb,an,homo,elumo,efl,aa,bb,anel,dnell,norb,mn)
      if(abs(dnell).le.thr) goto 501
      if((dnell*dnelr).lt.0.d0) goto 502
      goto 504
      endif
  502 dnelm=0.d0
  505 efm=0.5d0*(efl+efr)
      call calcni(eorb,an,homo,elumo,efm,aa,bb,anel,dnelm,norb,mn)
      if(abs(dnelm).le.thr) goto 501
      if((dnelm*dnell).lt.0.d0) then
      efr=efm
      else
      efl=efm
      endif
      goto 505
  501 ss=0.d0
      do i=1,mn
      ss=ss+an(i)
      enddo
      write(6,'(/"naccuracy",4(1x,f11.6))') ss,anel,(ss-anel),efm
      write(6,'(/"ni",10(1x,f11.6))')
     +(an(i),i=1,norb)
      pksks(1:norbtr)=0.d0
      pvksmo(1:norbtr)=0.d0   
      do i=1,norb
      ii=i*(i+1)/2
      pksks(ii)=sqrt(2.0d0*an(i))
      enddo
      call tmtdag(pksks,norb,pvksmo,norb,vksmo,c)
      ehdum=0.d0
      excbb=0.d0
      call erep(pvksmo,norb,iint,ehdum,excbb)
      write(6,'(/,"BB energy:",5f16.10)') aa,bb,efm,excbb,
     +etsks+evnno+ehno+excbb+enuc
  555 vnomo(1:norb*norb)=0.d0
      vnoao(1:norb*norb)=0.d0
      pnono(1:norbtr)=0.d0
      scrtc(1:norb*norb)=0.d0
      occno(1:norb)=0.d0
clmm      call getdmp(idmp,isno,vnoao,occno,norb)
      open(123,file='vnoao')
      do i=1,norb
      ii=(i-1)*norb
      write(123,*) i, occno(i)
      write(123,'(6(1x,f14.10))')
     + (vnoao(ii+j),j=1,norb)
      enddo
      close(123)
c      write(6,'(''vnoao   : '',10(1x,f11.7))')
c     + (vnoao(i),i=1,norb*norb)
      scrtc(1:norb*norb)=vmoao(1:norb*norb)
      tol=1.0d-7
      call minvr(scrtc,tol,del,ier,norb) 
      open(122,file='vmoao1')
      write(122,*) (scrtc(i),i=1,norb*norb)
      close(122)
c      write(6,'(''pnomo   : '',10(1x,f11.5))')
c     + (pnomo(i),i=1,norbtr)
      call matml(vnoao,scrtc,vnomo,norb)
      open(124,file='vnomo')
      do i=1,norb
      ii=(i-1)*norb
      write(124,*) i, occno(i)
      write(124,'(6(1x,f14.10))')
     + (vnomo(ii+j),j=1,norb)
      enddo
      close(124)
      vksao(1:norb*norb)=0.d0
      call matml(vksmo,vmoao,vksao,norb)
      open(125,file='vksao')
      do i=1,norb
      ii=(i-1)*norb
      write(125,*) i
      write(125,'(6(1x,f14.10))')
     + (vksao(ii+j),j=1,norb)
      enddo
      close(125)
      open(127,file='vmoao')
      do i=1,norb
      ii=(i-1)*norb
      write(127,*) i
      write(127,'(6(1x,f14.10))')
     + (vmoao(ii+j),j=1,norb)
      enddo
      close(127)
c      occno(1:norb)=0.d0
c      occno(1)=0.97494d0
c      occno(2)=0.97295d0
c      occno(3)=0.012d0
c      occno(4)=0.00338d0
c      occno(5)=0.0002d0
c      occno(6)=0.00436d0
c      occno(7)=0.00d0
c      occno(8)=0.00d0
c      occno(9)=0.00042d0 
c      occno(10)=0.00039d0
      do i=1,norb
      ii=(i+1)*i/2
      pnono(ii)=sqrt(2.0d0*occno(i))
      enddo
      write(6,'(''occupations  : '',10(1x,f11.5))')
     + (occno(i),i=1,norb)
      pnomo(1:norbtr)=0.d0
      exno=0.d0
      call tmtdag(pnono,norb,pnomo,norb,vnomo,c)
c      do i=1,norb
c      ii=(i-1)*i/2
c      pnomoa(ii+i)=pnomo(ii+i)
c      if(i.gt.1) then
c      do j=1,i-1
c      pnomoa(ii+j)=0.5d0*pnomo(ii+j)
c      enddo
c      endif
c      enddo 
c      vnomoa(1:norb*norb)=0.d0
c      write(6,'(''pnomoa   : '',10(1x,f11.5))')
c     + (pnomoa(i),i=1,norbtr)
c      call sqrmat(pnomoa,norb)
c      call diaglz(pnomoa,vnomoa,norb)
c      write(6,'(''vnomoa   : '',10(1x,f11.5))')
c     + (vnomoa(i),i=1,norb*norb)
c      do iorb=1,norb
c      occnoa(iorb)=pnomoa((iorb-1)*norb+iorb)
c      enddo
c      write(6,'(''occnoa   : '',10(1x,f11.7))')
c     + (occnoa(i),i=1,norb)
      call erep(pnomo,norb,iint,ehdum,exno)
      write(6,'(/,"BB energy:",5f16.10)') aa,bb,ehdum,exno,
     +etsno+evnno+ehno+exno+enuc
      exso=exno
      nmos=2
c      nvirt=norb-nmos
c      pnovi(1:nvirt*(nvirt+1)/2)=0.d0
c      vinomo(1:nvirt*norb)=0.d0
c      pnomov(1:norbtr)=0.d0
c      do i=nmos+1,norb
c      ivirt=i-nmos
c      pnovi((ivirt+1)*ivirt/2)=2.d0*sqrt(occno(i))
c      iivir=(ivirt-1)*norb
c      iitot=(i-1)*norb
c      do j=1,norb
c      vinomo(iivir+j)=vnomo(iitot+j)
c      enddo
c      enddo
c      call tmtdag(pnovi,nvirt,pnomov,norb,vinomo,c)
c      exph=0.d0
c      call erep(pnomov,norb,iint,ehdum,exph)
      ii1(1:8)=0
      jj1(1:8)=0
      ll1(1:8)=0
      kk1(1:8)=0
      fij(1:norb*norb)=0.d0
      wij(1:norb)=0.d0
      iz(1)=1
      iz(2)=2
      iz(3)=4
      iz(4)=3
      iz(5)=6
      iz(6)=8
      iz(7)=7
      iz(8)=9
      iz(9)=5
      iz(10)=10
      iz(11)=12
      iz(12)=11
      iz(13)=16
      iz(14)=13
      iz(15)=14
      iz(16)=15
      iz(17)=18
      iz(18)=17
      iz(19)=20
      iz(20)=19
      iz(21)=21
      iz(22)=24
      iz(23)=22
      iz(24)=25
      iz(25)=23
      iz(26)=26
      iz(27)=27
      iz(28)=28
      iz(29)=29
      iz(30)=33
      iz(31)=30
      iz(32)=34
      iz(33)=31
      iz(34)=32
      iz(35)=35
      iz(36)=36
      iz(37)=37
      iz(38)=38
      iz(39)=39
      iz(40)=41
      iz(41)=40
      iz(42)=42
      iz(43)=43
      iz(44)=44
      iz(45)=45
      iz(46)=47
      iz(47)=46
      iz(48)=48  
clmm      call search(1,iint)
clmm  508 call find(iint)
clmm      call get(gin,nw)
c      write(6,'(/,"nw:",i5)') nw
clmm      if (nw.eq.0) goto 509
clmm      call unpack(gijkl,8,ijkl,4*nword)
clmm      do m=1,nword
clmm      if (abs(gin(m)).gt.1.d-12) then
clmm      ii1(1)=ijkl(1,m)
clmm      jj1(1)=ijkl(2,m)
clmm      kk1(1)=ijkl(3,m)
clmm      ll1(1)=ijkl(4,m)
c      write(6,'(/,"unperm numbers:",5i5)') m,ii1(1),jj1(1)
c     +,kk1(1),ll1(1)
clmm      call nrperm(ii1(1),jj1(1),kk1(1),ll1(1),np)
clmm      do iperm=1,np
clmm      k=0
clmm      do i=1,8
clmm      ii=(iz(i)-1)*norb
c      ii=(i-1)*norb
clmm      if(k.eq.0) then
clmm      iin=(iz(i+1)-1)*norb
clmm      else
clmm      iin=(iz(i-1)-1)*norb
clmm      endif
clmm      do j=1,8,2
c      do j=1,8
clmm      jj=(iz(j)-1)*norb
c      jj=(j-1)*norb
clmm      jjn=(iz(j+1)-1)*norb
clmm      fij(ii+iz(j))=fij(ii+iz(j))+vnomo(ii+ii1(iperm))*vnomo(jj+
clmm     +jj1(iperm))*vnomo(iin+kk1(iperm))*vnomo(jjn
clmm     ++ll1(iperm))*gin(m)
c      fij(ii+j)=fij(ii+j)+vnomo(ii+ii1(iperm))*vnomo(jj+jj1(iperm))
c     +*vnomo(ii+kk1(iperm))*vnomo(jj+ll1(iperm))*gin(m) 
clmm      enddo
clmm      if(k.eq.0) then
clmm      k=k+1
clmm      else
clmm      k=k-1
clmm      endif
clmm      enddo
c      write(6,'(/,"iperm:",5i5)') iperm,ii1(iperm),jj1(iperm)
c     +,kk1(iperm),ll1(iperm)
c      do i=1,norb
c      ii=(i-1)*norb
c      do i=nmos+1,norb
c      if(i.gt.(nmos+3)) goto 511
c      do j=1,nmos
c      jj=(j-1)*norb
c      fij(ii+j)=fij(ii+j)+vnomo(ii+ii1(iperm))*vnomo(jj+jj1
c     +(iperm))*vnomo(ii+kk1(iperm))*vnomo(jj+ll1(iperm))
c     +*gin(m)
c      enddo
c  511 continue
c      jj=(j-1)*norb
c      wij(ii+j)=wij(ii+j)+vnomo(ii+ii1(iperm))*vnomo(jj+jj1
c     +(iperm))*vnomo(ii+kk1(iperm))*vnomo(jj+ll1(iperm))
c     +*gin(m)
c       wij(i)=wij(i)+vnomo(ii+ii1(iperm))*vnomo(ii+jj1
c     +(iperm))*vnomo(ii+kk1(iperm))*vnomo(ii+ll1(iperm))
c     +*gin(m)
c  510 continue
c      enddo
c      enddo
clmm      enddo
clmm      endif
clmm      enddo
clmm      goto 508
  509 eetot=0.d0
      eetot1=0.d0
      k=0
      open(126,file='intt')
      do i=1,8
      ii=(iz(i)-1)*norb
      if(k.eq.0) then
      iin=(iz(i+1)-1)*norb
      else
      iin=(iz(i-1)-1)*norb
      endif
      do j=1,8,2
      jj=(iz(j)-1)*norb
      jjn=(iz(j+1)-1)*norb
      write(126,*) ii/norb+1,jj/norb+1,iin/norb+1
     +,jjn/norb+1
      write(126,*) fij(ii+iz(j))
      enddo
      if(k.eq.0) then
      k=k+1
      else
      k=k-1
      endif
      enddo
      close(126)
      k=0
      do i=1,2
c      do i=1,20
c      ii=(i-1)*norb
c      if(i.le.2) then 
c      frt=0.5d0
c      else
c      frt=-0.5d0
c      endif
      ii=(iz(i)-1)*norb
      if(k.eq.0) then
      iin=iz(i+1)
c      if(i.eq.2) then
      sig=1.0d0
c      sig=-1.0d0
      else
      iin=iz(i-1)
      sig=-1.0d0
c      sig=1.0d0
      endif
      do j=3,8,2
c      do j=3,4
      eetot1=eetot1+sig*sqrt(sqrt(occno(iz(i))*occno(iin)*
     +occno(iz(j))*occno(iz(j+1))))*fij(ii+iz(j))
      eetot=eetot+sig*0.5d0*sqrt((occno(iz(i))+occno(iin))
     +*(occno(iz(j))+occno(iz(j+1))))*fij(ii+iz(j))
c      eetot1=eetot1+2.0d0*sig*sqrt(occno(iz(i))*occno(iz(j)))
c     +*fij(ii+iz(j))
c      eetot1=eetot1+sig*frt*(sqrt(occno(1))+sqrt(occno(2)))
c     +*sqrt(occno(j))*fij(ii+j)
c      eetot1=eetot1+sig*frt*sqrt(occno(i)*occno(j))*fij(ii+j)
      enddo
c      do j=6,8
c      eetot1=eetot1+sig*frt*sqrt(occno(i)*occno(j))*fij(ii+j)
c      enddo
      if(k.eq.0) then
      k=k+1
      else
      k=k-1
      endif
      enddo
c      eetot=occno(iz(1))*fij((iz(1)-1)*norb+iz(1))-sqrt
c     +(occno(iz(1))*occno(iz(2)))*fij((iz(2)-1)*norb+iz(1))
c      eetot=(0.5d0*occno(1)*fij(1)+0.5d0*occno(2)*fij(norb+2)
c     +-sqrt(occno(1)*occno(2))*fij(2))*0.5d0
c      eetot=0.5d0*sqrt(occno(1)*occno(2))*(fij(1)+fij(norb+2)
c     +-2.0d0*fij(2))
      esum=enuc+eetot+evnno+etsno
c      eetot1=eetot-occno(iz(1))*fij((iz(1)-1)*norb+iz(1))+sqrt
c     +(occno(iz(1))*occno(iz(2)))*fij((iz(2)-1)*norb+iz(1))
      write(6,'(/,"DMFT energy:",6f16.10)') etsno,evnno,eetot,
     +enuc,esum,eetot1 
c  509 exno=exno-exph
c  509 exro=exno
c      exno=exno-exph
c      do i=nmos+2,norb
c      ii=(i-1)*norb
c      do j=nmos+1,i-1
c      exro=exro+2.0d0*sqrt(occno(i)*occno(j))*fij(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"test energy:",6f16.10)') etsno,evnno,exno,enuc
c     +,ehno,etsno+evnno+ehno+exro+enuc
c      coulom=0.d0
c      do i=1,norb
c      ii=(i-1)*norb
c      do j=1,norb
c      coulom=coulom+0.5d0*occno(i)*occno(j)*fij(ii+j)
c      enddo
c      enddo
c      write(6,'(/,"Coulomb energy:",10f16.10)') coulom
      expo=exno-occno(nmos+1)*wij(nmos+1)
c      exno=exno-occno(nmos+2)*wij(nmos+2)
c      do i=nmos+1,norb
      do i=nmos+2,norb
      expo=expo-occno(i)*wij(i)
c       do i=nmos+2,norb
c      ii=(i-1)*norb
c      do j=nmos+1,norb
c      if(i.eq.j) goto 507
c      exno=exno-occno(i)*wij(ii+i)
c      exno=exno-0.5d0*(occno(i)+occno(i)*occno(i))*wij(i)
c       exno=exno+sqrt(occno(i)*occno(j))*wij(ii+j)
c  507 continue
c      enddo
      enddo
      write(6,'(/,"BBC1 energy:",5f16.10)') aa,bb,efm,exno,
     +etsno+evnno+ehno+expo+enuc
      do i=1,nmos-1
      ii=(i-1)*norb
      do j=i+1,nmos
      expo=expo+(0.5d0*sqrt(occno(i))*sqrt(occno(j))-0.25d0*
     +occno(i)*occno(j))*fij(ii+j)*2.0d0
      enddo
      enddo
      write(6,'(/,"BBC2 energy:",5f16.10)') aa,bb,efm,expo,
     +etsno+evnno+ehno+expo+enuc
      goto 513
      do i=1,nmos-1
      ii=(i-1)*norb
      j=7
      expo=expo+(0.5d0*sqrt(occno(i))*sqrt(occno(j))-0.25d0*
     +occno(i)*occno(j))*fij(ii+j)*2.0d0
      enddo
  513 continue
      exto=expo
      exto1=expo
      a1=20.0d0
      a2=300.0d0
      b1=1.94d0
      b2=0.03d0
      do i=1,norb
      e1=-a1*(occno(i)-b1)
      e2=-a2*(occno(i)-b2)
      if(e1.gt.40.0d0) then
      d1(i)=1.0d0
      d3(i)=0.0d0
      else
      ex1=exp(e1)
      d1(i)=ex1/(1.0d0+ex1)
      d3(i)=1.0d0/(1.0d0+ex1)
      endif
      if(e2.gt.40.0d0) then
      d2(i)=0.0d0
      else
      ex2=exp(e2)
      d2(i)=1.0d0/(1.0d0+ex2)
      endif
      c1=d1(i)*d2(i)
      c2=1.0d0-c1 
c      exto=exto-occno(i)*wij(i)
c       do i=nmos+2,norb
c      ii=(i-1)*norb
c      do j=nmos+1,norb
c      exno=exno-occno(i)*wij(ii+i)
      if((i.eq.nmos).or.(i.eq.(nmos+1))) goto 512
      exto=exto+0.5d0*(occno(i)-0.5d0*occno(i)*occno(i))*wij(i)
  512 continue      
c      exto1=exto1+0.5d0*((1.0d0-c1)*occno(i)-0.5d0*c2
c     **occno(i)*occno(i))*wij(i)
      exto1=exto1+0.5d0*((1.d0-d2(i))*occno(i)-0.5d0*(1.0d0-d2(i))
     **occno(i)*occno(i))*wij(i)
      enddo
      i=nmos+1
      ii=(i-1)*norb
      do j=1,nmos-1
      exto=exto+(sqrt(occno(i)*occno(j))-0.5d0*occno(i) 
     **occno(j))*fij(ii+j)
      enddo
      write(6,'(/,"BBC3 energy:",5f16.10)') aa,bb,efm,exto,
     +etsno+evnno+ehno+exto+enuc
c      do i=nmos+1,norb
c      ii=(i-1)*norb
c      do j=1,nmos
c      c3=d2(i)*d3(j)
c      exto1=exto1+(c3*sqrt(occno(i)*occno(j))
c     *-0.5d0*c3*occno(i)*occno(j))*fij(ii+j)
c      enddo
c      enddo
      write(6,'(/,"BBC4 energy:",5f16.10)') aa,bb,efm,exto,
     +etsno+evnno+ehno+exto1+enuc
      do i=1,norb
      exso=exso+0.5d0*(occno(i)-0.5d0*occno(i)*occno(i))*wij(i)
      enddo
      write(6,'(/,"BB-GU energy:",5f16.10)') aa,bb,efm,exso,
     +etsno+evnno+ehno+exso+enuc
      vksno(1:norb*norb)=0.d0
      scrtc(1:norb*norb)=vnomo(1:norb*norb)
      call minvr(scrtc,tol,det,ier,norb)
      call matml(vksmo,scrtc,vksno,norb)
c      call prvec(vksno,norb,norb,' ks orbitals in no basis ')
      return
      end
