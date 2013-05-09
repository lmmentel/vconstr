      subroutine calcni(eorb,an,homo,elumo,ef,aa,bb,anel,
     +dnel,norb,mn)
      implicit real*8(a-h,o-z),integer(i-n)
      dimension eorb(norb),an(norb)
      dnel=0.d0
      ee=elumo-homo
      do i=1,mn
      an(i)=2.0d0/(1.0d0+exp((eorb(i)-ef)/sqrt(aa*ee
     ++bb*bb*ee*ee)))
      dnel=dnel+an(i)
      enddo
      dnel=dnel-anel
      return
      end
