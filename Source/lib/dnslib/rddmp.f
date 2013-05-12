      subroutine rddmp(vmopao,vmoao,pmomo,pnomo,norb,nvpr,rnel,idmp,
     + ismo,isao,isno)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-8)
c
      dimension vmopao(norb*norb),vmoao(norb*norb),
     + pmomo(norb*(norb+1)/2),
     + pnomo(norb*(norb+1)/2),
     + pnono(norb*(norb+1)/2)
      dimension occ(norb),scrtc(norb*norb),vaopao(norb*norb),
     + vnomo(norb*norb),vnoao(norb*norb)
c
c** get MO basis
c
      call getdmp(idmp,ismo,vmoao,occ,norb)
      rhfnel=sum(occ)
c
c** vmoao contains the orbitals in (symmetry adapted or primitive) ao basis
c
      pmomo(1:norb*(norb+1)/2)=0.d0
      do i=1,norb
        imo=i*(i+1)/2
        pmomo(imo)=occ(i)
      enddo
c
      if (isao.gt.0) then
        call getdmp(idmp,isao,vaopao,occ,norb)
c
c** transform MO orbitals from symmetry adapted to primitive ao basis 
c**  vmopao = vmoao*vaopao
c
        call matml(vmoao,vaopao,vmopao,norb)
      else
        vmopao(1:norb*norb)=vmoao(1:norb*norb)
      endif
c
c** get Natural Orbitals in ao basis
c
      call getdmp(idmp,isno,vnoao,occ,norb)
c
      mvpr = 0
      rnel = 0.d0
      pnono(1:norb*(norb+1)/2)=0.d0
      do i=1,norb
        ii=i*(i+1)/2
        pnono(ii)=occ(i)
        rnel = rnel+occ(i)
        if (occ(i).gt.1.d-2) mvpr=mvpr+1
      enddo
      if (abs(rnel-rhfnel).gt.eps) then
        write(6,'(/''WARNING; sum of natorb occupation numbers :'',
     + f8.3)') rnel
        write(6,'(''         sum of mo occupation numbers     :'',
     + f8.3)') rhfnel
      endif
      if (nvpr.le.mvpr) nvpr=mvpr
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
c** transform natural orbitals from ao to mo basis
c**   vnomo = vnoao*scrtc
c
      call matml(vnoao,scrtc,vnomo,norb)
      call prvec(vnomo,norb,nvpr,' natural orbitals in mo-basis ')
c
c** transform the occupation matrix of the natural orbitals from no to mo basis
c**   pnono = 1/vnomo*pnono*vnomo
c
      call tmtdag(pnono,norb,pnomo,norb,vnomo,5.d-1)
      write(6,'(//)')
c
c** print density matrix in HF-MO basis
c
      write(6,'(//)')
      call wmtrx('####  ci-density matrix in mo basis  ####',
     + pnomo,norb,1.d-2)
      write(6,'(//)')
      call wmtrx('####  mo-density matrix in mo basis  ####',
     + pmomo,norb,1.d-20)
c
      return
      end
