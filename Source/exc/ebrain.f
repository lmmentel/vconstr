      subroutine ebrain(norb,npnt,npold,rnel,idmp,ismo,isno,isks,isao,
     + atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-6)
c
      logical atmol4
      dimension occmo(norb)
      dimension pdnmo(norb*(norb+1)/2),pdndn(norb*(norb+1)/2)
      dimension vmoao(norb*norb),vmopao(norb*norb),vaopao(norb*norb),
     + vdnmo(norb*norb),vdnao(norb*norb),scrtc(norb*norb)
      dimension vxc(npnt,7),exc(npnt,7)
      dimension rho(npnt,1),drho(npnt,3,1),d2rho(npnt,6,1)
      dimension grid(npnt,3),weight(npnt)
c
      vxc(1:npnt,1:7)=0.d0
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
c** ensemble representation
c
      if (isks.gt.0) then
c
c** get Kohn-Sham Orbitals in ao basis
c
      call getdmp(idmp,isks,vdnao,occmo,norb)
c
      occorb=2.d0
      do imo=1,norb
        if ((abs(occmo(imo)-occorb).gt.eps).and.
     + (abs(occmo(imo)).gt.eps)) ndets=ndets+1
      enddo
c
      else
c
c** get Natural Orbitals in ao basis
c
        call getdmp(idmp,isno,vdnao,occmo,norb)
c
      endif
c
c** transform orbitals from ao to mo basis
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
c** transform orbitals from ao to mo basis
c**   vdnmo = vdnao*scrtc
c
      call matml(vdnao,scrtc,vdnmo,norb)
c
      mvpr = 0
      sno = 0.d0
      pdndn(1:norb*(norb+1)/2)=0.d0
      do i=1,norb
        ii=i*(i+1)/2
        sno = sno+occmo(i)
        if (occmo(i).gt.1.d-6) mvpr=mvpr+1
        pdndn(ii)=occmo(i)
      enddo
      if (abs(sno-rnel).gt.eps) then
        write(6,'(/''ERROR; sum of occupation numbers :'',
     + f8.3)') sno
        write(6,'(''       number of electrons       :'',
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
      call tmtdag(pdndn,norb,pdnmo,norb,vdnmo,5.d-1)
c
      call wmtrx('#####  density matrix in mo basis  #####',
     + pdnmo,norb,1.d-2)
c
      call clcrho(rho,drho,d2rho,norb,npnt,npntmx,rnel,pdnmo,
     + vmopao,grid,weight,atmol4)
      call clcvxc(npntmx,npnt,rho,drho,d2rho,vxc,exc)
      call clcexc(npnt,npntmx,exc,weight)
c
      if (isks.gt.0) then
c
        write(6,'(//''  Now for ensemble''/)')
        call cledet(occmo,pdndn,vdnmo,vmopao,rho,grid,weight,
     + rnel,ndets,nmos,norb,npnt,intpnt,vxc,exc,atmol4)
c
      endif
      call wrtpot(rho,exc,vxc,weight,intpnt,npnt,npntmx)
c
      write(6,'(//)')
c
      return
      end
