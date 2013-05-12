      subroutine sbrain(occmo,nspin,nsp,nmos,nmomx,norb,npnt,npold,
     + itrx,tstthr,damp1,damp2,dqmax,kpens,nsens,idamp,nvpr,thresh,
     + iint,idmp,ismo,isao,isks,isoe,ikli,lb88,lpw91x,lprdwc,lpw91c,
     + llyp,lsym,lintsm,lrdocc,ldisk,atmol4)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-8)
c
      logical lsym,lintsm,lrdocc,lb88,lpw91x,lprdwc,lpw91c,llyp,
     + ldisk,atmol4
      dimension occmo(nmomx,nspin)
c
      logical ltrian,lcens,lcgap
c
      dimension ndvx(nspin)
      dimension dxyz(3)
      dimension wsi(nmomx),tsorb(nmomx),vnorb(nmomx)
      dimension ityp(nmomx,nspin),wxi(nmomx,nspin)
      dimension pksks(norb*(norb+1)/2,nspin),
     + pksmo(norb*(norb+1)/2,nspin),vxcmat(norb*(norb+1)/2),
     + hmatx(norb*(norb+1)/2),vnmat(norb*(norb+1)/2),
     + tsmat(norb*(norb+1)/2),vhrmat(norb*(norb+1)/2),
     + ptksmo(norb*(norb+1)/2),pmo(norb*(norb+1)/2),
     + pvksmo(norb*(norb+1)/2,nspin),pksao(norb*(norb+1)/2),
     + pksiks(norb*(norb+1)/2),pksimo(norb*(norb+1)/2),
     + pksiao(norb*(norb+1)/2)
      dimension vksmo(norb*norb,nsp),vmoao(norb*norb),hksmat(norb*norb),
     + vmopao(norb*norb),eorb(norb)
      dimension eorbin(nsens*nmomx),qorb(nsens*nmomx),chkval(nmomx)
      dimension vxc(npnt,nspin),vxcn(ikli*npnt,nspin),exc(npnt,4),
     + vresp(ikli*npnt,nspin)
      dimension rho(npnt,nspin),drho(npnt,3,nspin),d2rho(npnt,6,nspin),
     + grid(npnt,3),weight(npnt),valmo(npnt*norb)
      dimension orbdns(ikli*npnt,nmomx,nspin)
c
      dimension iorbsm(norb)
      dimension iorder(nmomx),itord(nsens*nmomx)
c
      allocatable norbsm(:,:)
c
c
c** kli
c
      parameter(lqntmx = 4)
      parameter(ncmpx = 1+3+6+10+15)
      parameter(l1l2mx = 2*lqntmx)
      dimension ncomp(0:lqntmx),ncbas(0:lqntmx),ll(ncmpx),mm(ncmpx),
     + nn(ncmpx)
      dimension xintmx(ikli*norb*(norb+1)/2),
     + pxksmo(ikli*norb*(norb+1)/2),pxksao(ikli*norb*(norb+1)/2)
c
      allocatable hxyz(:,:),prefac(:)
c
      data ncomp/1,3,6,10,15/
      data ncbas/0,1,4,10,20/
      data ll/0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,1,0,1,0,1,
     +        4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data mm/0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1,
     +        0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data nn/0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,0,1,2,2,1,
     +        0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
c
c
      ltrian = .false.
      lcens  = .false.
      lcgap  = .false.
c
      norbtr = norb*(norb+1)/2
c
      nstrt = 0
      ntmos = 0
c
      npntmx = npnt
      intpnt = npold-1
c
      occorb=2.d0
      if (nspin.eq.2) occorb=1.d0
c
      npks=0
      if ((nspin.eq.2).and.(nsp.eq.1)) npks=1
c
      kkk=1
      if (kpens.ne.0) lcens = .true.
      parin = 10.d0
      epse  = 1.d-4
      epsq  = 1.d-8
c
      ityp(1:nmomx,1:nspin) = 0
c
      do ipnt=1,npnt
        read(99,*) grid(ipnt,1),grid(ipnt,2),grid(ipnt,3),weight(ipnt)
      enddo
      close(99)
c
c** get one-electron matrix (on mo-basis) from dumpfile
c
      call hmatmo(idmp,isoe,hmatx,norb)
c
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
c** read dumpfile : pmo    - hf density matrix in mo basis
c**                 vmopao - molecular orbitals in primitive ao basis
c**                 vmoao  - molecular orbitals in ao basis
c
      call rdhfmo(pmo,vmopao,vmoao,norb,nvpr,idmp,ismo,isao)
      if (ldisk) then
        call clcmof(valmo,npnt,npntmx,norb,vmopao,grid,weight,atmol4)
      else
        call clcmod(valmo,npnt,npntmx,norb,vmopao,grid,weight,atmol4)
      endif
c
c** if not specified, prepare occupation numbers in occmo
c
      if (lrdocc) then
        do imo=norb,1,-1
          iimo=imo*(imo+1)/2
          if (abs(pmo(iimo)).gt.eps) then
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
        if (nspin.eq.1) then
          do imo=1,nmos
            iimo=imo*(imo+1)/2
            occmo(imo,1)=pmo(iimo)
          enddo
        else
          do imo=1,nmos
            iimo=imo*(imo+1)/2
            if (pmo(iimo).lt.1.d0) then
              occmo(imo,1)=pmo(iimo)
              occmo(imo,2)=0.d0
            else
              occmo(imo,1)=1.d0
              occmo(imo,2)=pmo(iimo)-1.d0
            endif
          enddo
        endif
      endif
c
c** set density matrices, only lower diagonal
c
      pksks(1:norbtr,1:nspin)=0.d0
      do is=1,nspin
        do imo = 1,nmos
          pksks(imo*(imo+1)/2,is) = occmo(imo,is)
        enddo
      enddo
      pksmo(1:norbtr,1:nspin)=pksks(1:norbtr,1:nspin)
c
      if (nsp.eq.2) then
        write(6,'(//)')
        call wmtrx('## ks-alpha-density matrix in mo basis ##',
     + pksmo(1,1),norb,1.d-2)
        write(6,'(//)')
        call wmtrx('##  ks-beta-density matrix in mo basis ##',
     + pksmo(1,2),norb,1.d-2)
        write(6,'(//)')
      else
        write(6,'(//)')
        call wmtrx('####  ks-density matrix in mo basis  ####',
     + pksmo(1,1),norb,1.d-2)
        write(6,'(//)')
      endif
c
c** store electron-nuclear attraction integrals from  dumpfile, section 192,
c** in the vector vnmat
c
      call rdmat(idmp,norb,tsmat,vnmat,enuc)
      write(6,'(/''  nuclear repulsion energy : '',f16.8,/)')enuc
c
      if (ldisk) then
        call clrhof(rho,drho,d2rho,nspin,norb,npnt,npntmx,pksmo,
     + valmo,grid,weight)
      else
        call clrhod(rho,drho,d2rho,nspin,norb,npnt,npntmx,pksmo,
     + vmopao,valmo,grid,weight,atmol4)
      endif
c
      if (ikli.ne.0) then
c
        open(90,file='vxc.dat',status='old',err=222)
        write(6,'('' initial potential read from file'')')
        rewind(90)
        do ip=1,npnt
          read(90,*) (vxc(ip,is),is=1,nspin)
        enddo
        close(90)
        goto 223
  222   call prepvs(vxc,rho,drho,nspin,npnt)
c
  223   continue
c
c** calculate some constants (kli)
c
        call fctrl
        call cnstnd
c
c*** calculate hx,hy and hz for every pair of gaussians
c
        if (atmol4) then
c
c* calculate transformation matrices for transformation from
c* cartesian to spherical gaussians.  stored in common/spheri/.
c
          call spherw
c
          call cnwhmx(nhmx,nprmx,ngrpx)
          allocate(hxyz(nhmx,3),stat=ialloc)
          if(ialloc.ne.0) then
            write(6,'(/'' No hxyz allocation '')')
            stop
          endif
          allocate(prefac(nprmx),stat=ialloc)
          if(ialloc.ne.0) then
            write(6,'(/'' No hxyz allocation '')')
            stop
          endif
c
          call calchw(hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,prefac,
     + nprmx,l1l2mx)
        else
c
c*** idddd = 0 : 6-d set is used;  idddd<>0 : 5-d set is used.
c
          call dset(idddd)
          call spherv(idddd)
c
          call cnvhmx(nhmx,nprmx,ngrpx)
          allocate(hxyz(nhmx,3),stat=ialloc)
          if (ialloc.ne.0) then
            write(6,'(/'' No hxyz allocation '')')
            stop
          endif
c
          call calchv(hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,prefac,
     + nprmx,l1l2mx)
        endif
c
      else
        call clcvxc(nspin,npntmx,npnt,rho,drho,d2rho,vxc,exc,
     + lb88,lpw91x,lprdwc,lpw91c,llyp)
      endif
c
c** In the following cycles, pksmo will converge to the KS-density matrix
c
      damp  = damp1
      dampi = 1.d0-damp
c
      do 100 itr=0,itrx
        exks=0.d0
        eonks=0.d0
        etsks=0.d0
        evnks=0.d0
        tester=0.d0
c
        if (itr.eq.idamp) then
          damp=damp2
          dampi=1.d0-damp
        endif
c
        if (npks.eq.1) call spavxc(vxc,rho,nspin,npnt,npntmx)
c
        ptksmo(1:norbtr)=0.d0
        do is=1,nspin
          ptksmo(1:norbtr)=ptksmo(1:norbtr)+pksmo(1:norbtr,is)
        enddo
        call clcvhr(ptksmo,norb,iint,vhrmat)
        write(6,'(/,'' -----  itr '',i4,'' -----'')')itr
c
c** start spin iteration
c
        do 50 is=1,nsp
c
c** calculate matrix elements vxcmat in mo basis.
c
          if (lsym) then
            call vksmts(vxcmat,vxc(1,is),valmo,weight,norbsm,nsym,
     + intpnt,norb)
          else
            call vksmat(vxcmat,intpnt,norb,vxc(1,is),valmo,weight)
          endif
c
          hksmat(1:norbtr)=hmatx(1:norbtr)+vhrmat(1:norbtr)+
     + vxcmat(1:norbtr)
c
c** expand hksmat in triangular form to full norb*norb matrix
c
          call sqrmat(hksmat,norb)
c
c** diagonalize matrix hksmat and return the eigenvalues on the main 
c** diagonal of hksmat. vksmo returns the eigenvectors in mo basis.
c
          call diaglz(hksmat,vksmo(1,is),norb)
          do iorb=1,norb
            eorb(iorb)=hksmat((iorb-1)*norb+iorb)
          enddo
c
          iorder(1:nmomx)=0
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
c
          do imo=nmos,1,-1
            if (occmo(iorder(imo),is).gt.1.d-6) then
              ehomo=eorb(imo)
              nmopr=imo
              exit
            endif
          enddo
c
          do imo=nmos+1,norb
            if (eorb(imo).lt.ehomo) then
              write(6,'(/''WARNING; Spectrum is incorrect'')')
              if (nmos.lt.nmomx) then
                nmopr=nmopr+1
                enew=eorb(imo)
                do jmo=nmos,1,-1
                  if (enew.gt.eorb(jmo)) then
                    do kmo=nmos,jmo,-1
                      eorb(kmo+1)=eorb(kmo)
                      iorder(kmo+1)=iorder(kmo)
                    enddo
                    eorb(jmo)=enew
                    iorder(jmo)=imo
                    exit
                  endif
                enddo
                nmos=nmos+1
              else
                write(6,'(/''ERROR; nmos > nmomx'')')
                write(6,'(''  nmos  = '',i3,/,''  nmomx = '',i3)')
     + nmos,nmomx
                stop
              endif
            endif
          enddo
          if (lcens) then
            eorbin(nstrt+1:nstrt+nmos)=eorb(1:nmos)
            if (is.eq.1) then
              itord(1:nmos)=iorder(1:nmos)
              nstrt=nmos
            else
              itord(nstrt+1:nstrt+nmos)=-iorder(1:nmos)
            endif
            ntmos=ntmos+nmos
          endif
c
          elumo=1.d6
          do imo=1,nmos
            if (occmo(imo,is).lt.epsq) then
              if (eorb(imo).lt.elumo) elumo=eorb(imo)
            endif
          enddo
          if (elumo.gt.1.d5) then
            do iorb=nmos+1,norb
              if (eorb(iorb).lt.elumo) elumo=eorb(iorb)
            enddo
          endif
c
          if (lsym) then
c
c** delete symmetry contaminating elements from Kohn-Sham vectors
c
            do i=1,norb
              xnorm=0.d0
              ii=(i-1)*norb
              do j=1,norb
                ij=ii+j
                if (iorbsm(i).eq.iorbsm(j)) then
                  xnorm=xnorm+vksmo(ij,is)*vksmo(ij,is)
                else
                  vksmo(ij,is)=0.d0
                endif
              enddo
            enddo
            if (abs(xnorm-1.d0).gt.1.d-6) then
              write(6,'(''WARNING Orbital '',i3,'' has '',g8.2,
     + ''% admixture of wrong symmetry'')')i,100*abs(xnorm-1.d0)
            endif
            if (xnorm.lt.1.d-1) then
              write(6,'(''ERROR in normalization. Orbital '',
     + i3,'' has norm '',f5.2)')i,xnorm
              stop
            endif
c
c** normalize KS vector
c
            xnorm=1.d0/sqrt(xnorm)
            do j=1,norb
              vksmo(ii+j,is)=vksmo(ii+j,is)*xnorm
            enddo
c
            write(6,'(/''  symm  : '',10(7x,i3,2x))')
     + (iorbsm(iorder(i)),i=1,nmopr)
c
          endif
c
c** check the KS orbital density
c
          chkval(1:nmopr)=0.d0
          do ipnt = 1,intpnt
            k=(ipnt-1)*norb
            do imo=1,nmos
              tmp=0.d0
              i=iorder(imo)
              ii=(i-1)*norb
              do j=1,norb
                tmp=tmp+vksmo(ii+j,is)*valmo(k+j)
              enddo
              chkval(imo)=chkval(imo)+tmp*tmp*weight(ipnt)
            enddo
          enddo
          chkval(1:nmopr)=chkval(1:nmopr)-1.d0
c
          write(6,'(''  integ : '',10(1x,g11.3))')
     + (chkval(i),i=1,nmopr)
          write(6,'(''  occ   : '',10(1x,f11.5))')
     + (occmo(iorder(i),is),i=1,nmopr)
          write(6,'(/''  eorb  : '',10(1x,f11.5))')
     + (eorb(i),i=1,nmopr)
c
c** transform KS-density from KS-orbital basis to HF-MO basis
c** pksks was already initialized with the occupation number input
c
          c = 0.5d0
c
c** transform matrix A to B = QAQ+ ; B(kl)=Q(ki)*A(ij)*Q(lj) sum over i,j
c**   before the matrix multiplication, the off-diagonal elements of 
c**   matrix A are multiplied by the constant c. for triangular matrix, 
c**   c should be equal to 0.5
c
          call tmtdag(pksks(1,is),nmos,pvksmo(1,is),norb,vksmo(1,is),c)
c
c** transform KS-density matrix from MO to AO basis ; pvksmo => pksao
c
          call tmtdag(pvksmo(1,is),norb,pksao,norb,vmoao,c)
c
c** calculate the energies
c
          do iorb=1,norbtr
            eonks=eonks+hmatx(iorb)*pvksmo(iorb,is)
            etsks=etsks+tsmat(iorb)*pksao(iorb)
            evnks=evnks+vnmat(iorb)*pksao(iorb)
          enddo
c
          call erep(pvksmo(1,is),norb,iint,rdum,extmp)
          if (nspin.eq.2) extmp=2*extmp
          exks=exks+extmp
c
c** Calculate orbital kinetic and electron nuclear energies
c
          pksiks(1:norbtr)=0.d0
          do imo=1,nmos
            etsi=0.d0
            evni=0.d0
            iimo=imo*(imo+1)/2
            pksiks(iimo)=occmo(imo,is)
            call tmtdag(pksiks,nmos,pksimo,norb,vksmo(1,is),c)
            call tmtdag(pksimo,norb,pksiao,norb,vmoao,c)
            do iorb=1,norbtr
              etsi=etsi+tsmat(iorb)*pksiao(iorb)
              evni=evni+vnmat(iorb)*pksiao(iorb)
            enddo
            tsorb(imo)=etsi
            vnorb(imo)=evni
            pksiks(iimo)=0.d0
          enddo
c
          write(6,'(''  ts    : '',10(1x,f11.5))')
     + (tsorb(iorder(i)),i=1,nmopr)
          write(6,'(''  vn    : '',10(1x,f11.5))')
     + (vnorb(iorder(i)),i=1,nmopr)
c
c** calculate tester
c
          do ipnt=1,intpnt
            m=(ipnt-1)*norb+1
            pks=valmat(pvksmo(1,is),norb,valmo(m),ltrian)
            tester=tester+abs(pks-rho(ipnt,is))*weight(ipnt)
          enddo
c
          if (ikli.ne.0) then
            it=0
            ityp(1:nmos,is)=0
            do imo=1,nmos
              iimo=iorder(imo)
              occimo=occmo(iimo,is)
              if (occimo.gt.epsq) then
                if (ityp(iimo,is).eq.0) then
                  it=it+1
                  ityp(iimo,is)=it
                  do jmo=imo+1,nmos
                    jjmo=iorder(jmo) 
                    occjmo=occmo(jjmo,is)
                    if (abs(occimo-occjmo).lt.epsq) then
                      if (abs(eorb(imo)-eorb(jmo)).lt.epse) then
                        if (lsym) then
                          if (iorbsm(jjmo).ne.iorbsm(iimo)) then
                            ityp(jjmo,is)=ityp(iimo,is)
                          else
                            it=it+1
                            ityp(jjmo,is)=it
                          endif
                        else
                          ityp(jjmo,is)=ityp(iimo,is)
                        endif 
                      endif
                    endif
                  enddo
                endif
              endif
            enddo
            ndvx(is)=it
            write(6,'(/''  ityp  : '',10(7x,i3,2x))')
     + (ityp(iorder(i),is),i=1,nmopr)
          endif
c
   50   continue
c
        if (kpens.ne.0) then
c
c** determine possible ensemble
c
        if ((lcens).and.(tester.lt.sqrt(tstthr))) then
          write(6,'(/'' ** change occupation numbers if necessary**'')')
          lcgap=.true.
        endif
        if (lcgap) then
c
c** calculate new occupation pattern
c
        if (nspin.eq.2) then
          do imo=1,ntmos-1
            emin=eorbin(imo)
            do jmo=imo+1,ntmos
              if (eorbin(jmo).lt.emin) then
                eorbin(imo)=eorbin(jmo)
                eorbin(jmo)=emin
                emin=eorbin(imo)
                itmp=itord(imo)
                itord(imo)=itord(jmo)
                itord(jmo)=itmp
              endif
            enddo
          enddo
          do imo=1,ntmos
            jmo=itord(imo)
          if (jmo.lt.0) then
            is=2
            jmo=abs(jmo)
          else
            is=1
            endif
            qorb(imo)=occmo(jmo,is)
          enddo
        else
          qorb(1:nmos)=occmo(iorder(1:nmos),1)
        endif
        sumi=sum(qorb)
c
c** determine fermi energy
c
        efermi = eorbin(1)
        ifermi = 1
        ilumo = -1
        elumo = 0.d0
        do iorb = 1,ntmos
          if ((qorb(iorb).gt.eps).and.(eorbin(iorb).gt.efermi)) then
            efermi = eorbin(iorb)
            ifermi = iorb
          endif
          if ((abs(qorb(iorb)-occorb).gt.eps).and.
     + (eorbin(iorb).lt.elumo)) then
            elumo = eorbin(iorb)
            ilumo = iorb
          endif
        enddo
c
        gap = efermi-elumo
        if (gap.lt.epse) then
          lcgap = .false.
        else
c
c** charge transfer
c
          nq = ifermi-ilumo+1
          if (kpens.lt.0) kkk=nq
          call frcctr(nq,eorbin(ilumo),qorb(ilumo),occorb,dqmax,parin,
     + epse,epsq,kkk)
          sumo=sum(qorb) 
          if (abs(sumi-sumo).gt.1.d-8) then
            write(6,'(''ERROR; number of electrons destroyed'')')
            write(6,'(''  in  : '',f6.2,/,''  out : '',f6.2)')sumi,sumo
            stop
          endif
          pksks(1:norbtr,1:nspin)=0.d0
          if (nspin.eq.2) then
            do imo=1,ntmos
              jmo=itord(imo)
              if (jmo.lt.0) then
                is=2
                jmo=abs(jmo)
              else
                is=1
              endif
              occmo(jmo,is)=qorb(imo)
              pksks(jmo*(jmo+1)/2,is) = qorb(imo)
            enddo
          else
            do imo=1,nmos
              jmo=iorder(imo)
              occmo(jmo,1) = qorb(imo)
              pksks(jmo*(jmo+1)/2,1) = qorb(imo)
            enddo
          endif
          do is=1,nspin
            call tmtdag(pksks(1,is),nmos,pvksmo(1,is),norb,
     + vksmo(1,is),c)
          enddo
        endif
c
          endif
        endif
        nstrt = 0
        ntmos = 0
c               
c** calculate new input density
c           
        if (itr.ne.0) then
          pksmo(1:norbtr,1:nspin)=dampi*pksmo(1:norbtr,1:nspin)+
     + damp*pvksmo(1:norbtr,1:nspin)
        else
          pksmo(1:norbtr,1:nspin)=pvksmo(1:norbtr,1:nspin)
        endif
c
        if (ldisk) then
          call clrhof(rho,drho,d2rho,nspin,norb,npnt,npntmx,
     + pksmo,valmo,grid,weight)
        else
          call clrhod(rho,drho,d2rho,nspin,norb,npnt,npntmx,
     + pksmo,vmopao,valmo,grid,weight,atmol4)
        endif
c
        if (ikli.ne.0) then
c
          do ip=1,npnt
            m=(ip-1)*norb+1
c
c*** calculate integrals on ao-basis.  store in array xintmx
c
            xpnt=grid(ip,1)
            ypnt=grid(ip,2)
            zpnt=grid(ip,3)
c
            if (atmol4) then
              call ccintw(xintmx,norb,ncomp(0),ncbas(0),lqntmx,ll,mm,nn,
     + ncmpx,hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,xpnt,ypnt,zpnt,prefac,
     + nprmx,l1l2mx,ngrpx)
            else
              call ccintv(xintmx,norb,ncomp(0),ncbas(0),lqntmx,ll,mm,nn,
     + ncmpx,hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,xpnt,ypnt,zpnt,prefac,
     + nprmx,l1l2mx,ngrpx)
            endif
c
            do is=1,nspin
              vs=0.d0
c
              call exchan(pvksmo(1,is),norb,pxksmo,valmo(m))
              call tmtdag(pxksmo,norb,pxksao,norb,vmoao,5.d-1)
c
              do i=1,norb
                do j=1,i
                  ij=i*(i-1)/2+j
                  vs=vs+pxksao(ij)*xintmx(ij)
                enddo
              enddo
              vxcn(ip,is)=vs
            enddo
          enddo
c
          call clcwxi(wxi,occmo,vksmo,nspin,nmos,nmomx,norb,iint)
c
          orbdns(1:npnt,1:nmos,1:nspin)=0.d0
          do is=1,nspin
            do imo=1,nmos
              wsitmp=0.d0
              l=ityp(imo,is)
              if (l.ne.0) then
                iimo=(imo-1)*norb
                do ip=1,npnt
                  k=(ip-1)*norb
                  tmp=0.d0
                  do j=1,norb
                    tmp=tmp+valmo(k+j)*vksmo(iimo+j,is)
                  enddo
                  wsitmp=wsitmp+tmp*tmp*vxcn(ip,is)*weight(ip)
                  orbdns(ip,l,is)=orbdns(ip,l,is)+occmo(imo,is)*tmp*tmp
                enddo
                wsi(imo)=wsitmp
              endif
            enddo
c
            write(6,'(/''  ws    : '',10(1x,f11.5))')
     + (wsi(iorder(i)),i=1,nmopr)
            write(6,'(''  wx    : '',10(1x,f11.5))')
     + (wxi(iorder(i),is),i=1,nmopr)
c
            if (ndvx(is).gt.1) then
              call clvkli(vxcn(1,is),vresp(1,is),wxi(1,is),wsi,
     + orbdns(1,1,is),weight,occmo(1,is),ityp(1,is),ndvx(is),nmos,
     + npnt,npntmx)
            endif
c
            if (itr.ne.0) then
              vxc(1:npnt,is)=dampi*vxc(1:npnt,is)+damp*vxcn(1:npnt,is)
            else
              vxc(1:npnt,is)=vxcn(1:npnt,is)
            endif
          enddo
        else
          call clcvxc(nspin,npntmx,npnt,rho,drho,d2rho,vxc,exc,
     + lb88,lpw91x,lprdwc,lpw91c,llyp)
        endif
c
c** dipole and Hartree
c
        write(6,'(/)')
        ptksmo(1:norbtr)=0.d0
        do is=1,nspin
          ptksmo(1:norbtr)=ptksmo(1:norbtr)+pksmo(1:norbtr,is)
        enddo
        call erep(ptksmo,norb,iint,ehks,rdum)
c       
c** transform KS-density matrix from MO to AO basis ; pvksmo => pksao
c
        call tmtdag(ptksmo,norb,pksao,norb,vmoao,c)  
        call clcdip(dxyz,pksao,idmp,norb)
        write(6,'(''  dx = '',f14.10,''    dy = '',f14.10,
     + ''    dz = '',f14.10)')dxyz(1),dxyz(2),dxyz(3)
c    
        write(6,'(/''  error : '',f16.10)')tester
        write(6,'(''  ts    : '',f12.6)')etsks
        write(6,'(''  en    : '',f12.6)')evnks
        write(6,'(''  eh    : '',f12.6)')ehks
        write(6,'(''  ex    : '',f12.6)')exks
        write(6,'(''  ef    : '',f12.6)')eonks-etsks-evnks
        write(6,'(''  total : '',f12.6)')eonks+ehks+exks+enuc
        write(6,'(/''  elumo : '',f8.4)')elumo
c
        if (tester.lt.tstthr) goto 110
c
  100 continue
c
c** end of cycle
c
      write(6,'(/,'' +++++ NO CONVERGENCE +++++'')')
      goto 120
  110 write(6,'(/,'' +++++ CONVERGED +++++'')')
c
c** store density on dumpfile if requested. Set section-type to 50
c
  120 continue
c
      if (nsp.eq.2) then
        write(6,'(/,''             alpha '')')
      call prvec(vksmo(1,1),norb,nvpr,'Kohn-Sham orbitals in mo basis')
        write(6,'(/,''             beta '')')
      call prvec(vksmo(1,2),norb,nvpr,'Kohn-Sham orbitals in mo basis')
        write(6,'(//)')
        call wmtrx('## ks-alpha-density matrix in mo basis ##',
     + pksmo(1,1),norb,1.d-2)
        write(6,'(//)')
        call wmtrx('##  ks-beta-density matrix in mo basis ##',
     + pksmo(1,2),norb,1.d-2)
        write(6,'(//)')
      else
        call prvec(vksmo,norb,nvpr,'Kohn-Sham orbitals in mo basis')
        write(6,'(//)')
        call wmtrx('####  ks-density matrix in mo basis  ####',
     + pksmo(1,1),norb,1.d-2)
        write(6,'(//)')
        if (isks.gt.0) then
          call wrtkso(vksmo,vmoao,hksmat,occmo,nmos,norb,idmp,isks)
        endif
        call prvec(vksmo,norb,nvpr,'Kohn-Sham orbitals in mo basis')
      endif
c
      if (ldisk) then
        call clrhof(rho,drho,d2rho,nspin,norb,npnt,npntmx,
     + pksmo,valmo,grid,weight)
      else
        call clrhod(rho,drho,d2rho,nspin,norb,npnt,npntmx,
     + pksmo,vmopao,valmo,grid,weight,atmol4)
      endif
      if (ikli.eq.0) then
        call clcvxc(nspin,npntmx,npnt,rho,drho,d2rho,vxc,exc,lb88,
     + lpw91x,lprdwc,lpw91c,llyp)
        if (npks.eq.1) call spavxc(vxc,rho,nspin,npnt,npntmx)
        call clcexc(exld,exnl,ecld,ecnl,intpnt,npntmx,exc,weight)
      endif
c
      exks=0.d0
      eonks=0.d0
      etsks=0.d0
      evnks=0.d0
c
      ptksmo(1:norbtr)=0.d0
      do is=1,nspin
        ptksmo(1:norbtr)=ptksmo(1:norbtr)+pksmo(1:norbtr,is)
      enddo
      call erep(ptksmo,norb,iint,ehks,rdum)
c
      do is=1,nspin
c
c** calculate the Coulomb and Exchange energy associated
c** with the diagonal one-particle density matrix
c
        call erep(pksmo(1,is),norb,iint,rdum,extmp)
        if (nspin.eq.2) extmp=2*extmp
        exks=exks+extmp

        call tmtdag(pksmo(1,is),norb,pksao,norb,vmoao,5.d-1)
        write(6,'(//)') 
        call wmtrx('####  ks-density matrix in ao basis  ####',
     + pksao,norb,1.d-2)
        write(6,'(//)')

        do iorb=1,norbtr
          eonks=eonks+hmatx(iorb)*pksmo(iorb,is)
          etsks=etsks+tsmat(iorb)*pksao(iorb)
          evnks=evnks+vnmat(iorb)*pksao(iorb)
        enddo
      enddo
c
c** calculate Kohn-Sham kinetic energy which is the
c** one-electron energy minus the electron nuclear
c** attraction energy
c
      write(6,'(//,'' KS total energy                        :'',
     +f16.10)') eonks+ehks+exld+exnl+ecld+ecnl+enuc
      write(6,'('' KS kinetic energy                      :'',
     +f16.10)') etsks
      write(6,'('' KS electron-nuclear attraction energy  :'',
     +f16.10)') evnks
      write(6,'('' KS coulomb repulsion energy            :'',
     +f16.10,/)') ehks
c
      if (ikli.eq.0) then
        write(6,'('' KS local density exchange energy       :'',
     +f16.10)')exld
        write(6,'('' KS non-local exchange energy           :'',
     +f16.10)')exnl
        write(6,'('' KS total exchange energy               :'',
     +f16.10)')exld+exnl
      endif
      write(6,'('' KS orbital-exchange energy             :'',
     +f16.10)')exks
c
      if (ikli.eq.0) then
        write(6,'(/'' KS local density correlation energy    :'',
     +f16.10)')ecld
        write(6,'('' KS non-local correlation energy        :'',
     +f16.10)')ecnl
        write(6,'('' KS total correlation energy            :'',
     +f16.10)')ecld+ecnl
        write(6,'(/'' KS total xc-energy                     :'',
     +f16.10)')exld+exnl+ecld+ecnl
      endif
c
c** write data to the files vxc.dat, dnst.dat
c
      if (ikli.ne.0) then
        open(70,file='dns.dat')
        rewind(70)
        open(80,file='vxc.dat')
        rewind(80)
        do ip=1,npnt
          write(70,*) (rho(ip,is),
     + (orbdns(ip,idvx,is),idvx=1,ndvx(is)),is=1,nspin)
          write(80,*) (vxc(ip,is),vresp(ip,is),rho(ip,is),
     + (orbdns(ip,idvx,is),idvx=1,ndvx(is)),is=1,nspin)
        enddo
        open(71,file='dns.plt')
        rewind(71)
        write(71,'('' dns        orbdns'')')
        open(81,file='vxc.plt')
        rewind(81)
        write(81,'('' vxckli       vrespkli'')')
        do ip=npold,npnt
          write(81,'(4(F14.8))') (rho(ip,is),
     + (orbdns(ip,idvx,is),idvx=1,ndvx(is)),is=1,nspin)
          write(81,'(4(F14.8))') (vxc(ip,is),vresp(ip,is),is=1,nspin)
        enddo
      else
        call wrtpot(rho,vxc,exc,nsp,intpnt,npntmx)
      endif
c
      return
      end
