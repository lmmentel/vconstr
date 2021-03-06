cam      subroutine vbrain(npnt,npold,norb,idmp,ismo,isno,isao,isks,atmol4)
      subroutine vbrain(npnt,npold,norb,nmomx,nele,inpksoFile)
      use commonsModule
      use integralsModule
      use ioModule
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.0d-8)
c     
      character*100 inpksoFile
c
      logical ltrian
      dimension pksks(norb*(norb+1)/2),pmo(norb*(norb+1)/2),
     + pnomo(norb*(norb+1)/2),pksmo(norb*(norb+1)/2)
      dimension vksao(norb*norb),vksmo(norb*norb),vmoao(norb*norb),
     + vmopao(norb*norb),scrtc(norb*norb),occmo(norb),orbdns(norb)

      dimension rho(npnt),vnuc(npnt)
      dimension grid(npnt,4)
      dimension dns(npnt),ddns(npnt,3),dsdns(npnt)
      dimension valmo(norb*npnt)
      allocatable grdmo(:,:)
c
      ltrian = .false.
c
      allocate(grdmo(norb*npnt,4),stat=ialloc)
      if (ialloc.ne.0) then
        write(6,'(/'' No memory allocated for derivativess'')')
        stop
      endif
c
      rewind(99)
      do ipnt = 1, npnt
        read(99,'(4e25.14)',iostat=ios) grid(ipnt,1),grid(ipnt,2),
     &                                  grid(ipnt,3),grid(ipnt,4)
      enddo
      close(99)
*c
c
c** read dumpfile : vmopao - molecular orbitals in primitive ao basis
c**                 vmoao  - molecular orbitals in ao basis
c**                 pmo    - mo density matrix in mo basis
c**                 pnomo  - ci density matrix in mo basis
c
cam      call rddmp(vmopao,vmoao,pmo,pnomo,norb,nvpr,rnel,idmp,ismo,isao,
cam     + isno)
cam..replacement for the rddmp routine now read the orbitals
      call getOrbitalsAndDensities(vmopao, 
     & vmoao, pmo, pnomo, norb, 
     & rnel, nele)
      if (printLevel >= 2) then 
        call matPrint(reshape(vmopao, (/norb, norb/)), 
     & 'Vmo in primitive ao, vmopao')
      endif
c     
c** get Kohn-Sham Orbitals in ao basis
c
cam      call getdmp(idmp,isks,vksao,occmo,norb)
      open(26,file=trim(inpksoFile))
      occmo = 0.d0
      do i=1,nmomx
         read(26,*) occmo(i)
      enddo
      do i=1,norb*norb
         read(26,*) vksao(i)
      enddo
      close(26)
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
      call clcvls(dns,ddns,dsdns,vnuc,valmo,grdmo,npnt,
     + norb,pnomo,vmopao,grid,rnel)
c
      do ip=1,npnt
        k=(ip-1)*norb
        m=k+1
        rho(ip)=valmat(pksmo,norb,valmo(m),ltrian)
        do i=1,nmos
          ii=(i-1)*norb
          tmp=0.d0
          do j=1,norb
            tmp=tmp+vksmo(ii+j)*valmo(k+j)
          enddo
          orbdns(i)=occmo(i)*tmp*tmp
        enddo
        if (ip.ge.npold) then
          z=grid(ip,3)
          print*,z,rho(ip),orbdns(1:nmos)
        endif
      enddo
c
      call clvkin(norb,npold,npnt,pnomo,pksmo,rho,dns,ddns,grdmo,grid)
c
      return
      end
