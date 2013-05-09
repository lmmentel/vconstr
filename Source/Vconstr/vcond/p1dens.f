      subroutine p1dens(iunit,pmo,norb,rnel)
c*******************************************************************
c  calculation of 1-density from 2-density
c
c  density stored in array pmo in triangular form (that is: the sum of
c  of the off-diagonal elements (i,j) and (j,i) is stored in position
c  pmo(i*(i-1)/2+j),  assuming i>j)
c
c arguments:
c iunit : unit-number of two-matrix file
c pmo   : array for storage of density
c norb  : number of orbitals on two-matrix file 
c rnel  : number of electrons
c*******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension pmo(norb*(norb+1)/2)
      common/atbuf/coef(340),rijkl(170),nword,ndum
      common/scra/ijkl(340,4)
c
c*** initialize elements of array pmo
c
      imb=0
      nmb=0
      pmo(1:norb*(norb+1)/2)=0.d0 
      call search(1,iunit)
  100 call fget(coef,nw,iunit)
c
c* end-file block is reached if nw=0. (note: nw is the total number of
c* elements in the block, maximum 511. the number of integrals or coef-
c* elements is given by nword, maximum 340)
c
      if (nw.ne.0) then
        imb=imb+1
c
c*** unpack ijkl-indices
c
        call upak8v(rijkl(1),ijkl(1,1))
c
c*** loop over coef elements in this block
c
        do m = 1,nword
          nmb=nmb+1
          i = ijkl(m,1)
          j = ijkl(m,2)
          k = ijkl(m,3)
          l = ijkl(m,4)
          if ((i.lt.1).or.(j.lt.1).or.(k.lt.1).or.(l.lt.1).or.
     + (i.gt.norb).or.(j.gt.norb).or.(k.gt.norb).or.
     + (l.gt.norb)) then
            write(6,'(/,
     + '' ERROR; nonexistent determinant '',i3,1x,i3,1x,i3,1x,i3,
     + '' encountered'')')i,j,k,l
            stop
          else
            ij = i*(i-1)/2+j
            kl = k*(k-1)/2+l
c
c** add contribution to density matrix
c
            if (k.eq.l) pmo(ij) = pmo(ij)+coef(m)
            if (i.eq.j) pmo(kl) = pmo(kl)+coef(m)
          endif
        enddo
c
c*** proceed reading and processing blocks until end-file block is reached
c
        goto 100
      endif
c
c*** normalize density matrix
c
      fctr = 1.0d0/(rnel-1.d0)
      do m = 1,norb*(norb+1)/2
        pmo(m) = pmo(m)*fctr
      enddo
c
      rclnel = 0.d0
      do imo=1,norb
        iimo=imo*(imo+1)/2
        rclnel=rclnel+pmo(iimo)
      enddo
      write(6,'(/''  number of p2(ijkl)-elements restored from file :'',
     + i16,2x,i16)')nmb,imb
      call wmtrx(' one-electron-density matrix in mo basis ',
     + pmo,norb,1.d-2)
      if (abs(rclnel-rnel).gt.1.d-6) then
        write(6,'(''WARNING; rnel : '',f16.10,2x,g16.10)')rnel,rclnel
        write(6,'(''  diff : '',g16.6)')abs(rnel-rclnel)
      endif
c
      return
      end
