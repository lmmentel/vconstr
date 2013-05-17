      subroutine clcvhr(pmo, norb, vhrmat, gintegfile)
       use integralsModule
       use varModule
       use commonsModule
c***********************************************************************
c
c Calculate the Hartree potential associated with the one-particle 
c  density matrix stored in pmo.
c
c arguments
c  pmo    : density matrix
c  norb   : dimension of density matrix
c  iunit  : unit number of integral file
c  vhrmat : Hartree potential in MO basis (returned)
c
c***********************************************************************
      implicit real*8(a-h,o-z),integer(i-n)
c
      dimension pmo(norb*(norb+1)/2),vhrmat(norb*(norb+1)/2)
      dimension ij(norb,norb),ii(8),jj(8),ll(8),kk(8)
c
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
clmm..matrix for storing two electron integrals
      real(dp), allocatable :: twoEint(:)
clmm..name of the gamess-us tow electron integral file
      character(*) gintegfile 
c
      vhrmat(1:norb*(norb+1)/2)=0.d0
c
c*** array ij is used for storage of i*(i-1)/2+j ***
c
      do i = 1,norb
        do j = 1,i
          ij(i,j) = i*(i-1)/2+j
          ij(j,i) = ij(i,j)
        enddo
      enddo
c
c* multiply off-diagonal elements of pmo with 0.5
c* (multiply again with 2.0 before leaving routine)
c
      do i = 2,norb
        do j = 1,i-1
          mn = i*(i-1)/2+j
          pmo(mn) = pmo(mn)*0.5d0
        enddo
      enddo
c
clmm      call search(1,iunit)
clmm  100 call find(iunit)
clmm      call get(gin,nw)
clmm      if (nw.eq.0) goto 200
clmm      call unpack(gijkl,8,ijkl,4*nword)
clmm      do m = 1,nword
clmm        if (abs(gin(m)).gt.1.d-10) then
clmm          ii(1) = ijkl(1,m)
clmm          jj(1) = ijkl(2,m)
clmm          kk(1) = ijkl(3,m)
clmm          ll(1) = ijkl(4,m)
clmm          call nrperm(ii(1),jj(1),kk(1),ll(1),np)
clmm          do iperm = 1,np
clmm            vhrmat(ij(kk(iperm),ll(iperm)))=
clmm     + vhrmat(ij(kk(iperm),ll(iperm)))+
clmm     + pmo(ij(ii(iperm),jj(iperm)))*gin(m)
clmm          enddo
clmm        endif
clmm      enddo
clmm      goto 100
clmcc
clmm  200 continue
clmm
clmm..calculate the dimnsion for the two electron integral vector
clmm..and allocate the vector 

      n1 = norb*(norb+1)/2
      n2 = n1*(n1+1)/2
      allocate(twoEint(n2))
      twoEint = 0.0d0

clmm..read the two electron integrals

      call readTwoIntMO(twoEint, trim(gintegfile))

clmm..calculate the  hartree potential and store it in 'vhrmat'

      ijindx = 0
      do i = 1, norb
        do j = 1, i
          ijindx = ijindx + 1
          kl = 0
          cumul = 0.0d0
          do k = 1, norb
            do l = 1, k
              kl = kl + 1
            if ((ijindx >= kl).and.
     &(abs(twoEint(addr(i,j,k,l))) > 1.0d-10)) then
          call nrperm(ii(1),jj(1),kk(1),ll(1),np)
          do iperm = 1,np
            vhrmat(ij(kk(iperm),ll(iperm)))=
     + vhrmat(ij(kk(iperm),ll(iperm)))+
     + pmo(ij(ii(iperm),jj(iperm)))*twoEint(addr(i,j,k,l))
          enddo
              endif
            enddo
          enddo
        enddo
      enddo   
      if (printLevel > 2) call print4i(twoEint, norb)
clmm..two electron integrals are no longer needed
     
      deallocate(twoEint)
c
c** multiply off-diagonal elements of pmo with 2.0
c
      do i = 2,norb
        do j = 1,i-1
          mn = i*(i-1)/2+j
          pmo(mn) = pmo(mn)*2.0d0
          vhrmat(mn)=0.5d0*vhrmat(mn)
        enddo
      enddo
c
      return
      end

      subroutine nrperm(ii,jj,kk,ll,n)
c***********************************************************************
c
c given an integral (or p2) index (ii jj kk ll) this routine
c calculates all possible permutations. the number of possible
c permutations is returned in n.
c
c***********************************************************************
      implicit real*8  (a-h,o-z),integer(i-n)
c
      dimension ii(8),jj(8),kk(8),ll(8)
c
      i = ii(1)
      j = jj(1)
      k = kk(1)
      l = ll(1)
      n = 1
c
      if (i.ne.j) then
        n = n+1
        ii(n) = j
        jj(n) = i
        kk(n) = k
        ll(n) = l
      endif
c
      if (k.ne.l) then
        n = n+1
        ii(n) = i
        jj(n) = j
        kk(n) = l
        ll(n) = k
      endif
c
      if ((i.ne.j).and.(k.ne.l)) then
        n = n+1
        ii(n) = j
        jj(n) = i
        kk(n) = l
        ll(n) = k
      endif
c
      if ((i.ne.k).or.(j.ne.l)) then
        n = n+1
        ii(n) = k
        jj(n) = l
        kk(n) = i
        ll(n) = j
        if (i.ne.j) then
          n = n+1
          ii(n) = k
          jj(n) = l
          kk(n) = j
          ll(n) = i
        endif
c
        if (k.ne.l) then
          n = n+1
          ii(n) = l
          jj(n) = k
          kk(n) = i
          ll(n) = j
        endif
c
        if ((i.ne.j).and.(k.ne.l)) then
          n = n+1
          ii(n) = l
          jj(n) = k
          kk(n) = j
          ll(n) = i
        endif
c
      endif
c
      return
      end

