module hartreePotentialModule
  use varModule
  use commonsModule
 
  implicit none

  private

  interface getHartreePot
    module procedure getHartreePotential
  end interface

  public :: getHartreePot 

contains 

  subroutine getHartreePotential(hartreePotential, pmo, norb) 
    real(DP),         intent(out)   :: hartreePotential(:)
    real(DP),         intent(inout) :: pmo(:)
    integer,          intent(in)    :: norb 

    real(DP),          allocatable  :: buffer(:)
    integer(IPgamess), allocatable  :: indexBuffer(:)

    integer(IK)       :: readStatus
    integer(IK)       :: m
    integer(IK)       :: i,j,k,l
    integer(IK)       :: bufLength, twoIntIndexBufSize, twoIntBufferSize
    integer(IPgamess) :: label, label1, label2
    integer(IPgamess) :: length, nintmx, labsiz, twoemo, iw, ipu, is, transInt
    integer(IPgamess) :: nft11,nft12,nft13,nft14,nft15,nft16
    logical           :: largeLabels    

    integer, allocatable  :: ij(:,:)
    integer, dimension(8) :: ii, jj, kk, ll
    integer :: mn, iperm, np 

    allocate(ij(norb, norb))

! array ij is used for storage of i*(i-1)/2+j 

    do i = 1,norb
        do j = 1,i 
            ij(i,j) = i*(i-1)/2+j
            ij(j,i) = ij(i,j)
        enddo
    enddo

! multiply off-diagonal elements of pmo with 0.5
! (multiply again with 2.0 before leaving routine)

    do i = 2,norb
        do j = 1,i-1
            mn = i*(i-1)/2+j
            pmo(mn) = pmo(mn)*0.5d0
        enddo
    enddo

    nintmx = 15000
    labsiz = 1 
    twoemo = twoein

    if (labsiz /= 1_IPgamess .and. labsiz /= 2_IPgamess) then
      write(*,*) 'RdTwoIntAO:  CONFUSION IN LABSIZ! '
      stop
    endif

    largeLabels = (labsiz == 2_IPgamess)

    twoIntBufferSize = int(nintmx, kind=IK)

    twoIntIndexBufSize = twoIntBufferSize
    if (largeLabels) then
      if (IPgamess == 4) twoIntIndexBufSize = 2*twoIntIndexBufSize
    else
      if (IPgamess == 8) twoIntIndexBufSize = (twoIntIndexBufSize + 1) / 2
    endif

    open(unit=twoemo, file=trim(integralsFile), status='old', form='unformatted')
    rewind(twoemo)
    read(twoemo)

    allocate(buffer(twoIntBufferSize))
    allocate(indexBuffer(twoIntIndexBufSize))

    length = 1_IPgamess
    do while (length > 0_IPgamess)
      Read(twoemo,iostat=readStatus) length,indexBuffer,buffer 
      if (readStatus /= 0) then
        if (readStatus == 2) then
          write(*,*) 'RdTwoIntMO: ENCOUNTERED UNEXPECTED END WHILE READING TWO-ELECTRON FILE'
        else
          write(*,*) 'RdTwoIntMO: ENCOUNTERED I/O ERROR WHILE READING TWO-ELECTRON FILE'
        endif
        stop
      endif

      bufLength = abs(int(length, kind=IK))
      if (bufLength > twoIntBufferSize) stop

      do m=1, bufLength
        if (IPgamess == 4) then                  ! 32-bit Gamess integers
          if (largeLabels) then       
            label1 = indexbuffer(2*m-1) 
            label2 = indexBuffer(2*m)
            i = ishft(label1, -16_IPgamess)
            j = iand( label1,  65535_IPgamess)
            k = ishft(label2, -16_IPgamess)
            l = iand( label2,  65535_IPgamess)
          else
            label = indexBuffer(m)
            i =      ishft(label, -24_IPgamess)
            j = iand(ishft(label, -16_IPgamess), 255_IPgamess)
            k = iand(ishft(label,  -8_IPgamess), 255_IPgamess)
            l = iand(      label,                255_IPgamess)
          endif  
        else                                     ! 64-bit Gamess integers
          if (largeLabels) then
            label = indexBuffer(m)
            i = int(     ishft(label,   -48_IPgamess),                  kind=IK)
            j = int(iand(ishft(label,   -32_IPgamess), 65535_IPgamess), kind=IK)
            k = int(iand(ishft(label,   -16_IPgamess), 65535_IPgamess), kind=IK)
            l = int(iand(      label,                  65535_IPgamess), kind=IK)  
          else
            if (mod(m,2) == 0) then
              label = indexBuffer(m/2)
              i = int(iand(ishft(label, -24_IPgamess ), 255_IPgamess), kind=IK)
              j = int(iand(ishft(label, -16_IPgamess ), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label,  -8_IPgamess ), 255_IPgamess), kind=IK)
              l = int(iand(      label,                 255_IPgamess), kind=IK)
            else
              label = indexBuffer(m/2+1)
              i = int(     ishft(label, -56_IPgamess),                kind=IK)
              j = int(iand(ishft(label, -48_IPgamess), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label, -40_IPgamess), 255_IPgamess), kind=IK)
              l = int(iand(ishft(label, -32_IPgamess), 255_IPgamess), kind=IK)
            endif      
          endif
        endif

        ii(1) = i
        jj(1) = j
        kk(1) = k
        ll(1) = l
        call nrperm(ii, jj, kk, ll, np)
        do iperm = 1, np
            hartreePotential(ij(kk(iperm), ll(iperm))) = hartreePotential(ij(kk(iperm), ll(iperm))) &
                + pmo(ij(ii(iperm), jj(iperm))) * buffer(m)
        enddo
 
      enddo
    enddo
   
! multiply off-diagonal elements of pmo with 2.0

    do i = 2, norb
        do j = 1, i-1
            mn = i*(i-1)/2+j
            pmo(mn) = pmo(mn)*2.0_dp
            hartreePotential(mn)=0.5_dp*hartreePotential(mn)
        enddo
    enddo

    close(twoemo)
    deallocate(buffer)
    deallocate(indexBuffer)
    deallocate(ij)
  end subroutine getHartreePotential 


  subroutine nrperm(ii,jj,kk,ll,n)
!***********************************************************************
!
! given an integral (or p2) index (ii jj kk ll) this routine
! calculates all possible permutations. the number of possible
! permutations is returned in n.
!
!***********************************************************************
    integer, intent(out) :: n
    integer, dimension(8), intent(inout) ::  ii, jj, kk, ll
    
    integer :: i, j, k, l

      i = ii(1)
      j = jj(1)
      k = kk(1)
      l = ll(1)
      n = 1

      if (i.ne.j) then
        n = n+1
        ii(n) = j
        jj(n) = i
        kk(n) = k
        ll(n) = l
      endif

      if (k.ne.l) then
        n = n+1
        ii(n) = i
        jj(n) = j
        kk(n) = l
        ll(n) = k
      endif

      if ((i.ne.j).and.(k.ne.l)) then
        n = n+1
        ii(n) = j
        jj(n) = i
        kk(n) = l
        ll(n) = k
      endif

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

        if (k.ne.l) then
          n = n+1
          ii(n) = l
          jj(n) = k
          kk(n) = i
          ll(n) = j
        endif

        if ((i.ne.j).and.(k.ne.l)) then
          n = n+1
          ii(n) = l
          jj(n) = k
          kk(n) = j
          ll(n) = i
        endif

      endif

  end subroutine nrperm

end module hartreePotentialModule
