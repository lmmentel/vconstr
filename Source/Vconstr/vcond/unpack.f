      subroutine upak8v(rijkl,i205)
c***********************************************************************
c     unpacks a 64-bit real by equivalencing it to two
c     32-bit integers by use of the commonblock cmmn in
c     subroutine upak8a. this is done in order to use the
c     logical operators ishift and iand which do only work
c     for integers and not for reals.
c***********************************************************************
      implicit real*8  (a-h,p-w),integer  (i-n)
      dimension rijkl(170),i205(1360)
c
c...       nbits 8
c
         common/cmmn/r
c
         iu = 1
         do 20 i=1,170
         r=rijkl(i)
         call upak8a(i1,i2,i3,i4,i5,i6,i7,i8)
         i205(iu)=i1
         i205(iu+340)=i2
         i205(iu+680)=i3
         i205(iu+1020)=i4
         i205(iu+1)=i5
         i205(iu+341)=i6
         i205(iu+681)=i7
         i205(iu+1021)=i8
20       iu = iu + 2
c...
         return
         end
c
         subroutine upak8a(i1,i2,i3,i4,i5,i6,i7,i8)
         implicit real*8  (a-h,p-w),integer  (i-n)
c
c        unpacks two 32-bit integers jhigh and jlow which
c        by use of the commonblock cmmn are equivalenced to
c        a 64-bit real in subroutine upak8.
c
			integer*4 jhi,jlo
         common/cmmn/jhi,jlo
c
         i1=ishft(iand(jhi,ishft(255,24)),-24)
         i2=ishft(iand(jhi,ishft(255,16)),-16)
         i3=ishft(iand(jhi,ishft(255,8)),-8)
         i4=ishft(iand(jhi,ishft(255,0)),-0)
         i5=ishft(iand(jlo,ishft(255,24)),-24)
         i6=ishft(iand(jlo,ishft(255,16)),-16)
         i7=ishft(iand(jlo,ishft(255,8)),-8)
         i8=ishft(iand(jlo,ishft(255,0)),-0)
c
         return
         end
