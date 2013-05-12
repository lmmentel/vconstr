      subroutine hmatmo(iunit,isec,hmat,norb)
c***********************************************************************
c
c get h-matrix on MO basis from section isec on dumpfile
c unit iunit, starting in block 1
c
c it is assumed that the matrix elements where stored on
c dumpfile by the atmol TRANS program (SPLICE directive).
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension hmat(norb*(norb+1)/2)
c
      nn = norb*(norb+1)/2
      call secini(1,iunit)
      call secget(isec,1004,iblock)
c
c** skip first 4 blocks **
c
      iblock = iblock+3
      call search(iblock,iunit)
      call reads(hmat,nn,iunit)
c
      return
      end
