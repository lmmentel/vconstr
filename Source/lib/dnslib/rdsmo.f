      subroutine rdsmo(vmopao,vmoao,norb,nvpr,idmp,ismo,isao,rnel)
c
c** Read HF-MOs in AO basis from the dumpfile and store in vmoao
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension vmoao(norb*norb),vmopao(norb*norb)
      dimension occ(norb),vaopao(norb*norb)
c
      call getdmp(idmp,ismo,vmoao,occ,norb)
      rtmp=sum(occ)
      if (abs(rnel-rtmp).gt.1.d-6) then
        write(6,'(/''ERROR; number of electrons found od dumpfile '',
     + ''not equal to the specified number.'',/,'' given : '',f8.4,/,
     + '' found : '',f8.4)')rnel,rtmp
        stop
      endif
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
