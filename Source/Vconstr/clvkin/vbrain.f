      subroutine vbrain(npnt,npold,norb,idmp,ismo,isno,isao,isks,atmol4)                                 
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.0d-8)
c
      logical atmol4,ltrian
      dimension pksks(norb*(norb+1)/2),pmo(norb*(norb+1)/2),
     + pnomo(norb*(norb+1)/2),pksmo(norb*(norb+1)/2)
      dimension vksao(norb*norb),vksmo(norb*norb),vmoao(norb*norb),
     + vmopao(norb*norb),scrtc(norb*norb),occmo(norb),orbdns(norb)

      dimension rho(npnt),vnuc(npnt)
      dimension grid(npnt,3),weight(npnt)
      dimension dns(npnt),ddns(npnt,3),dsdns(npnt)
      dimension valmo(norb*npnt)
      allocatable grdmo(:,:)
c
      ltrian = .false.
c
      npntmx=npnt
c
      allocate(grdmo(norb*npnt,4),stat=ialloc)
      if (ialloc.ne.0) then
        write(6,'(/'' No memory allocated for derivativess'')')
        stop
      endif
c
      do ip=1,npnt
        read(99,*) grid(ip,1),grid(ip,2),grid(ip,3),weight(ip)
      enddo
      close(99)
c
c
c** read dumpfile : vmopao - molecular orbitals in primitive ao basis
c**                 vmoao  - molecular orbitals in ao basis
c**                 pmo    - mo density matrix in mo basis
c**                 pnomo  - ci density matrix in mo basis
c
      call rddmp(vmopao,vmoao,pmo,pnomo,norb,nvpr,rnel,idmp,ismo,isao,
     + isno)
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
      call clcvls(dns,ddns,dsdns,vnuc,valmo,grdmo,npnt,npntmx,
     + norb,pnomo,vmopao,grid,weight,rnel,atmol4)
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
      call clvkin(norb,npold,npnt,pnomo,pksmo,rho,dns,ddns,grdmo)
c
      return
      end
