      subroutine wrtkso(vksmo,vmoao,hksmat,occmo,nmos,norb,idmp,isks)
c***********************************************************************
c Store ks-vectors on dumpfile. first transform from MO to AO basis.
c   The denstity matrix is already given in AO-basis
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension vksmo(norb*norb),vmoao(norb*norb),
     + hksmat(norb*norb),occmo(nmos)
      dimension scrtch(norb*norb),eorb(norb),occ(norb)
c
      occ(1:norb)=0.d0
      occ(1:nmos)=occmo(1:nmos)
      scrtch(1:norb*norb)=0.d0
      do i=1,norb
        eorb(i)=hksmat((i-1)*norb+i)
        iimo=(i-1)*norb
        iiao=(i-1)*norb
        do j=1,norb
          jjao=(j-1)*norb
          do k=1,norb
            scrtch(iiao+k) = scrtch(iiao+k)+
     + vksmo(iimo+j)*vmoao(jjao+k)
          enddo
        enddo
      enddo
      call dmpvec(idmp,isks,scrtch,eorb,occ,norb)
c
      return
      end
