      subroutine rdhfmo(pmo,vmopao,vmoao,norb,nvpr,idmp,ismo,isao)
c
c** Read HF-MOs in AO basis from the dumpfile and store in vmoao
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension pmo(norb*(norb+1)/2),vmoao(norb*norb),vmopao(norb*norb)
      dimension occ(norb),vaopao(norb*norb)
c
      call getdmp(idmp,ismo,vmoao,occ,norb)
      rtmp=sum(occ)
      write(6,'(/'' number of electrons found on dumpfile : '',
     + f8.4)')rtmp
c
      pmo(1:norb*(norb+1)/2)=0.d0
      do imo = 1,norb
        pmo(imo*(imo+1)/2) = occ(imo)
      enddo
c
      if (isao.gt.0) then
        call getdmp(idmp,isao,vaopao,occ,norb)
c
c** transform HF orbitals from symmetry adapted ao to primitive ao basis
c
        call matml(vmoao,vaopao,vmopao,norb)
      else
        vmopao(1:norb*norb)=vmoao(1:norb*norb)
      endif
c
c** vmoao contains the HF orbitals in AO basis
c
      call prvec(vmoao,norb,nvpr,'mo basis expressed in ao basis')
c
      return
      end
