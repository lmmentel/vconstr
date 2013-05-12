      subroutine ebrain(norb,npnt,npold,nspin,idmp,ismo,isks,isao,
     + atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical atmol4
      dimension occmo(norb)
      dimension pdnmo(norb*(norb+1)/2,nspin),
     + pdndn(norb*(norb+1)/2,nspin),ptdnmo(norb*(norb+1)/2)
      dimension vmoao(norb*norb),vmopao(norb*norb),vaopao(norb*norb),
     + vdnmo(norb*norb),vdnao(norb*norb),scrtc(norb*norb)
      dimension vxc(npnt,nspin,7),exc(npnt,7),grid(npnt,3),weight(npnt)
      dimension rho(npnt,nspin),drho(npnt,3,nspin),d2rho(npnt,6,nspin)
c
      vxc(1:npnt,1:nspin,1:7)=0.d0
      exc(1:npnt,1:7)=0.d0
c
      do ipnt=1,npnt
        read(99,*) grid(ipnt,1),grid(ipnt,2),grid(ipnt,3),weight(ipnt)
      enddo
      close(99)
      npntmx=npnt
      intpnt=npold-1
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
      call getdmp(idmp,isks,vdnao,occmo,norb)
      rksnel=sum(occmo)
      write(6,'(/,'' rnel '',f16.10,'' rksnel '',f16.10)')
     + rnel,rksnel
c     
c** transform orbitals from ao to mo basis
c
      scrtc(1:norb*norb)=vmoao(1:norb*norb)
c     
c** invert the matrix scrtc => scrtc
c
      tol = 1.d-7
      call minvr(scrtc,tol,det,ier,norb)
      if (ier.ne.0) then
        write(6,'(/,'' ERROR : transformation from ao to '',
     +       ''mo basis is singular'')')
        stop 
      endif  
c     
c** transform orbitals from ao to mo basis
c**   vdnmo = vdnao*scrtc
c
      call matml(vdnao,scrtc,vdnmo,norb)
c
      mvpr = 0
      pdndn(1:norb*(norb+1)/2,1:nspin)=0.d0
c
      do i=1,norb
        ii=i*(i+1)/2
        if (occmo(i).gt.1.d-6) then
          mvpr=mvpr+1
          if (nspin.eq.2) then
            if (abs(occmo(i)-1.d0).gt.1.d-6) then
              pdndn(ii,1)=1.d0
              pdndn(ii,2)=1.d0
            else
              pdndn(ii,1)=1.d0
            endif
          else
            pdndn(ii,1)=occmo(i)
          endif
        endif
      enddo
c
      write(6,'(//)') 
      if (nspin.eq.1) then
        call wmtrx('#####  density matrix in ks basis  #####',
     + pdndn(1,1),norb,1.d-2)
      else
        call wmtrx('### alpha-density matrix in ks basis ###',
     + pdndn(1,1),norb,1.d-2)
        write(6,'(/)')
        call wmtrx('### beta-density matrix in ks basis ####',
     + pdndn(1,2),norb,1.d-2)
      endif
c     
      do imo=norb,1,-1
        if (occmo(imo).gt.1.d-6) then
          nmos=imo
          exit
        endif 
      enddo   
      call tmtdag(pdndn(1,1),norb,pdnmo(1,1),norb,vdnmo,5.d-1)
      if (nspin.eq.2) then
        call tmtdag(pdndn(1,2),norb,pdnmo(1,2),norb,vdnmo,5.d-1)
      endif
c
      write(6,'(//)')
      if (nspin.eq.1) then
        call wmtrx('#####  density matrix in mo basis  #####',
     + pdnmo(1,1),norb,1.d-2)
      else
        call wmtrx('### alpha-density matrix in mo basis ###',      
     + pdnmo(1,1),norb,1.d-2)
        write(6,'(/)')
        call wmtrx('### beta-density matrix in mo basis ####',           
     + pdnmo(1,2),norb,1.d-2)
      endif
c
      ecdn=0.d0
      exdn=0.d0
      factor=1.d0
      if (nspin.eq.2) factor=2.d0
      ptdnmo(1:norb*(norb+1)/2)=0.d0
      do is=1,nspin
c        call erep(pdnmo(1,is),norb,7,ectmp,extmp)
        ecdn=ecdn+factor*ectmp
        exdn=exdn+factor*extmp
        print*,factor*ectmp,factor*extmp
        ptdnmo(1:norb*(norb+1)/2)=ptdnmo(1:norb*(norb+1)/2)+
     + pdnmo(1:norb*(norb+1)/2,is)
      enddo
c
      write(6,'(/'' exchange energy                     :''
     +,f16.10,/)')exdn
c
      write(6,'(//)') 
      call wmtrx('#####  density matrix in mo basis  #####',
     + ptdnmo,norb,1.d-2) 
c    
      call clcrho(rho,drho,d2rho,norb,npnt,npntmx,nspin,rnel,
     + pdnmo,vmopao,grid,weight,atmol4)
c     
      call clcvxc(nspin,npntmx,npnt,rho,drho,d2rho,vxc,exc)
      call clcexc(intpnt,npnt,exc,weight)
      call wrtpot(rho,drho,exc,vxc,weight,intpnt,npnt,npntmx,
     + nspin)
c
      return
      end
