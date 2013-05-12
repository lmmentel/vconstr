      subroutine clcpci(pcimo,occvno,norb)
c***********************************************************************
c
c Calculate CI-density matrix pcimo in HF-MO basis, first convert density 
c   matrix occvno to triangular form
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension pcimo(norb*(norb+1)/2),occvno(norb*norb)
c
      do 10 i=1,norb
        do 20 j=1,i-1
          occvno((i-1)*norb+j)=occvno((i-1)*norb+j)+
     + occvno((j-1)*norb+i)
 20     continue
 10   continue
      do 30 i=1,norb
        do 40 j=1,i
          pcimo(i*(i-1)/2+j)=occvno((i-1)*norb+j)
 40     continue
 30   continue
c
      return
      end
