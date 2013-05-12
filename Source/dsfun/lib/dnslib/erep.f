      subroutine erep(pmo,norb,iunit,ecoul,eexch)
c***********************************************************************
c
c Calculate the Coulomb and Exchange energy associated with the
c one-particle density matrix stored in pmo.
c
c arguments
c  pmo   : density matrix
c  norb  : dimension of density matrix
c  iunit : unit number of integral file
c  ecoul : Coulomb energy gin (returned)
c  eexch : exchange energy gin (returned)
c
c***********************************************************************
      implicit real*8(a-h,o-z),integer(i-n)
      parameter(tol  = 1.d-12)
c
      dimension pmo(norb*(norb+1)/2)
c      dimension vksmo(norb*norb)
      dimension ij(norb,norb),ii(8),jj(8),ll(8),kk(8)
c
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
c
c** array ij is used for storage of i*(i-1)/2+j 
c
      w2244=0.d0
      w2424=0.d0 
      do i = 1,norb
        do j = 1,i
          ij(i,j) = i*(i-1)/2+j
          ij(j,i) = ij(i,j)
        enddo
      enddo
c
c** multiply off-diagonal elements of pmo with 0.5
c** (multiply again with 2.0 before leaving routine)
c
      do i = 2,norb
        do j = 1,i-1
          mn = i*(i-1)/2+j
          pmo(mn) = 0.5d0*pmo(mn)
        enddo
      enddo
c
      ecoul = 0.d0
      eexch = 0.d0
c
clmm      call search(1,iunit)
clmm  100 call find(iunit)
clmm      call get(gin,nw)
      if (nw.eq.0) goto 200
clmm      call unpack(gijkl,8,ijkl,4*nword)
      do m = 1,nword
        if (abs(gin(m)).gt.tol) then
          ii(1) = ijkl(1,m)
          jj(1) = ijkl(2,m)
          kk(1) = ijkl(3,m)
          ll(1) = ijkl(4,m)
          call nrperm(ii(1),jj(1),kk(1),ll(1),np)
          do iperm = 1,np
            ecoul = ecoul+pmo(ij(ii(iperm),jj(iperm)))*
     + pmo(ij(kk(iperm),ll(iperm)))*gin(m)
            eexch = eexch+pmo(ij(ii(iperm),ll(iperm)))*
     + pmo(ij(kk(iperm),jj(iperm)))*gin(m)
c      w2244=w2244+vksmo(norb+ii(iperm))*vksmo(norb+jj(iperm))
c     +*vksmo(3*norb+kk(iperm))*vksmo(3*norb+ll(iperm))*gin(m)
c      w2424=w2424+vksmo(norb+ii(iperm))*vksmo(3*norb+jj(iperm))
c     +*vksmo(norb+kk(iperm))*vksmo(3*norb+ll(iperm))*gin(m)
          enddo
        endif
      enddo
clmm      goto 100
c
  200 ecoul = 0.5d0*ecoul
      eexch = -0.25d0*eexch
c
c** multiply off-diagonal elements of pmo with 2.0
c
      do i = 2,norb
        do j = 1,i-1
          mn = i*(i-1)/2+j
          pmo(mn) = 2.0d0*pmo(mn)
        enddo
      enddo
c
c      write(6,'(/"KS integrals: ", 2(1x,f30.20))') w2244,w2424 
      return
      end
