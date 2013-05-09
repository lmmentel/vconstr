      subroutine orbsym(iorbsm,ngrp,norb,thresh,iunit)
c*********************************************************************
c
c  use two-electron integrals to divide the orbitals into
c  symmetry equivalent groups and construct the multiplication
c  table. The first group is assumed to be the fully symmetric group
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(ngrpmx = 25)
c
      logical lchk 
      dimension iorbsm(norb)
      dimension ismgrp(norb),neworb(norb),ndim(ngrpmx),
     + mult(ngrpmx,ngrpmx)
c
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
c
      ngrp = 0
      nhigh = 0
      ismgrp(1:norb)=0
      mult(1:ngrpmx,1:ngrpmx)=0
c
      write(6,'(//,'' orbital symmetry and multiplication table'',
     +'' determined from two electron integrals'',
     +/,'' threshold = '',g8.2,/)')thresh
c
      call search(1,iunit)
c
c** read block with integrals ***
c
   10 call find(iunit)
      call get(gin,nw)
c
c** End-file block reached when nw=0 ***
c
      if (nw.eq.0) goto 70
c
c** unpack ijkl indices ***
c
      call unpack(gijkl,8,ijkl,4*nword)
c
c** loop over integrals
c
      do 30 m = 1,nword
        if (abs(gin(m)).gt.thresh) then
          i = ijkl(1,m)
          j = ijkl(2,m)
          k = ijkl(3,m)
          l = ijkl(4,m)
          if ((i.lt.1).or.(j.lt.1).or.(k.lt.1).or.(l.lt.1).or.
     + (i.gt.norb).or.(j.gt.norb).or.(k.gt.norb).or.
     + (l.gt.norb)) then
            write(6,'(/,
     + '' ERROR; nonexistent determinant '',i3,1x,i3,1x,i3,1x,i3,
     + '' encountered'')')i,j,k,l
            stop
          endif
c
c* search for two identical orbitals. the two remaining orbitals (ia
c* and ib) then belong to the same symmetry group
c
          if (i.eq.j) then
            ia = k
            ib = l
            goto 50
          elseif (i.eq.k) then
            ia = j
            ib = l
            goto 50
          elseif (i.eq.l) then
            ia = j
            ib = k
            goto 50
          elseif (j.eq.k) then
            ia = i
            ib = l
            goto 50
          elseif (j.eq.l) then
            ia = i
            ib = k
            goto 50
          elseif (k.eq.l) then
            ia = i
            ib = j
            goto 50
          endif
          goto 30
c
c** store highest orbital used for symmetry determination in nhigh
c
   50     if (ia.gt.nhigh) nhigh = ia
          if (ib.gt.nhigh) nhigh = ib
c
c** check if ia and ib have same symmetry ***
c
          if (ismgrp(ia).ne.0) then
c
c** orbital ia already assigned to a symmetry group ***
c
            if (ismgrp(ib).ne.0) then
              if (ismgrp(ia).ne.ismgrp(ib)) then
c
c* symmetry of orbital ia must be equal to symmetry of orbital ib.
c* this is stored in array mult
c
                mult(ismgrp(ia),ismgrp(ib)) = 1
                mult(ismgrp(ib),ismgrp(ia)) = 1
              endif
            else
c
c* No symmetry assigned to orbital ib.
c* Put orbital ib in same symmetry group as orbital ia.
c
              ismgrp(ib) = ismgrp(ia)
            endif
          else
            if (ismgrp(ib).ne.0) then
c
c* No symmetry assigned to orbital ia. 
c* Put orbital ia in same symmetry group as orbital ib
c
              ismgrp(ia) = ismgrp(ib)
            else
c
c* No symmetry assigned to both orbital ia and ib.
c* Create new symmetry group.
c
              ngrp = ngrp+1
              if (ngrp.gt.ngrpmx) then
                write(6,'(/,'' ERROR in orbsym. ngrpmx too small'')')
                stop
              endif
              ismgrp(ia) = ngrp
              ismgrp(ib) = ngrp
            endif
          endif
        endif
   30 continue
      goto 10
c
   70 if (nhigh.ne.norb) then
        write(6,'(/,'' ERROR; number of orbitals assigned '',i4,
     + '' is not equal to the total number of orbitals '',i4)')
     + nhigh,norb
      else
        write(6,'(/,'' number of orbitals on integral file :'',i4)')
     + norb
      endif
c
c** check if a symmetry is assigned to all orbitals ***
c
      do i=1,norb
        if (ismgrp(i).eq.0) then
          write(6,'(/,'' ** unable to assign symmetry to orbital :'',
     +       i3,/,''    try again with smaller threshhold'')')i
          stop
        endif
      enddo
c
c** correct for orbitals in different groups with same symmetry ***
c
      do igrp = ngrp,1,-1
        ndim(igrp) = 0
c
c** find group with same symmetry as group igrp ***
c
        do jgrp = 1,igrp-1
          if (mult(igrp,jgrp).eq.1) then
c
c** assign symmetry label jgrp to orbitals in group igrp
c
            do i = 1,norb
              if (ismgrp(i).eq.igrp) ismgrp(i) = jgrp
            enddo
c
c** no orbitals in group igrp. flag this by setting ndim(igrp) to -1
c
            ndim(igrp) = -1
            exit
          endif
        enddo
      enddo
c
c** renumber groups and delete empty groups ***
c
      kgrp = 0
      do igrp = 1,ngrp
        if (ndim(igrp).ge.0) then
          kgrp = kgrp+1
c
c** kgrp is new groupnumber for orbitals in group igrp ***
c
          do i = 1,norb
            if (ismgrp(i).eq.igrp) ismgrp(i) = kgrp
          enddo
        endif
      enddo
      ngrp = kgrp
c
      write(6,'(//,''*** orbitals divided into symmetry groups ***''/)')
c
      k = 0
      do i = 1,ngrp
        ndim(i) = 0
        k0 = k
        do j = 1,norb
          if (ismgrp(j).eq.i) then
            k = k+1
            neworb(k) = j
            ndim(i) = ndim(i)+1
          endif
        enddo
        write(6,'('' irrep : '',i2)')i
        write(6,'(12(1x,i3))')(neworb(k0+j),j=1,ndim(i))
      enddo
c
c* calculate multiplication table from two elec integrals
c* assume symmetry group 1 is of A1 symmetry
c
      mult(1:ngrpmx,1:ngrpmx)=0
      do i = 1,ngrp
        mult(i,i) = 1
        mult(i,1) = i
        mult(1,i) = i
      enddo
c
      x = 0.d0
      lchk = .false.
      call search(1,iunit)
  100 call find(iunit)
      call get(gin,nw)
      if (nw.eq.0) goto 200
      call unpack(gijkl,8,ijkl,4*nword)
      do m = 1,nword
        k = 0
        do ia = 1,4
          if (ismgrp(ijkl(ia,m)).eq.1) then
c
c** switch element ijkl(ia,m) and element ijkl(k,m)
c
            k = k+1
            ii = ijkl(ia,m)
            ijkl(ia,m) = ijkl(k,m)
            ijkl(k,m)  = ii
          endif
        enddo
c
        i = ijkl(1,m)
        j = ijkl(2,m)
        k = ijkl(3,m)
        l = ijkl(4,m)
c
        if (ismgrp(i).eq.1) then
          if (ismgrp(j).eq.1) then
c
c*** check if k and l belong to same symmetry ***
c
            if (ismgrp(k).ne.ismgrp(l)) then
              if (abs(gin(m)).gt.abs(x)) then
                x = gin(m)
                lchk = .true.
              endif
              if (abs(gin(m)).gt.thresh) then
                write(6,'(/'' forbidden integral ('',i3,2(1x,i3),1x,i3,
     + '') = '',g12.4,'' larger than '',g12.4)')i,j,k,l,gin(m),thresh
                write(6,'('' symmetry :         '',4(1x,i3),/)')
     + ismgrp(i),ismgrp(j),ismgrp(k),ismgrp(l)
                stop
              endif
            endif
          else
            isymkl = mult(ismgrp(k),ismgrp(l))
            if (abs(gin(m)).gt.thresh) then
              if ((isymkl.eq.0)) then
                mult(ismgrp(k),ismgrp(l)) = ismgrp(j)
                mult(ismgrp(l),ismgrp(k)) = ismgrp(j)
              else
                if (isymkl.ne.ismgrp(j)) then
                  write(6,'(/'' multiplication table error'')')
                  write(6,'('' forbidden integral ('',
     + i3,2(1x,i3),1x,i3,'') = '',g12.4,'' larger than '',g12.4)')
     + i,j,k,l,gin(m),thresh
                  write(6,'('' symmetry          :'',4(1x,i3),/)')
     + ismgrp(i),ismgrp(j),ismgrp(k),ismgrp(l)
                  stop
                endif
              endif
            else
c
c*** check for forbidden element ***
c
              if ((isymkl.ne.0).and.(isymkl.ne.ismgrp(j))) then
                if (abs(gin(m)).gt.abs(x)) then
                  x = gin(m)
                  lchk = .true.
                endif
              endif
            endif
          endif
        endif
c
      enddo
      goto 100
c
  200 do i = 1,ngrp
        do j = 1,ngrp
          if (mult(i,j).eq.0) then
            write(6,'('' *** cannot determine multiplication table'')')
            stop
          endif
        enddo
      enddo
c
      write(6,'(/,'' multiplication table : '',/)')
      write(6,'(''      '',20i3)')(i,i=1,ngrp)
      write(6,'(/)')
      do i = 1,ngrp
        write(6,'(i4,''  '',20i3)')i,(mult(i,j),j=1,ngrp)
      enddo
      if (lchk) then
        write(6,'(/,'' largest forbidden integral : '',g8.2)')x
      else
        write(6,'(/,'' No forbidden integrals'')')
      endif
c
c** symmetry of individual orbitals is stored in array iorbsm
c
        k = 0
        do isym = 1,ngrp
          do j = 1,ndim(isym)
            k = k+1
            iorbsm(neworb(k)) = isym
          enddo
        enddo
c
      return
      end
c
      subroutine hsym(iorbsm,nsym,hmat,norb,thresh)
c********************************************************************
c
c Determine orbital symmetry from h-integrals.
c Array iorbsm returns with a symmetry-number for every orbital.
c
c Method used: If the integral value <i|h|j> is smaller than a given
c threshhold then orbitals i and j are probably of different symmetry.
c Now loop over all other orbitals to check if an orbital k exists for
c which <i|h|k> and <j|h|k> are bigger than the threshhold value. If
c such an orbital exists then i and j (and k) belong to the same
c symmetry, otherwise i and j are assigned different symmetry numbers.
c
c*********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical lbadsm
c
      dimension hmat(norb*(norb+1)/2),iorbsm(norb),lll(200)
c
      nsym = 0
      do 310 i = 1,norb
c
c*** initially set the symmetry number of orbital i at 0.d0 ***
c
        iorbsm(i) = 0
        ii = i*(i-1)/2
        do 320 j = 1,i-1
          ij = ii + j
          if (abs(hmat(ij)).lt.thresh) then
            lbadsm = .true.
c
c orbitals i and j have probably different symmetries.
c check if orbital k exists for which h(i,k) and h(j,k) <> 0
c
            do 330 k = 1,norb
              if ( k .gt. i ) then
                ki = k*(k-1)/2+i
              else
                ki = i*(i-1)/2+k
              end if
              if (abs(hmat(ki)).gt.thresh) then
                if ( k .gt. j ) then
                  kj = k*(k-1)/2+j
                else
                  kj = j*(j-1)/2+k
                end if
                if (abs(hmat(kj)).gt.thresh) then
c
c orbital k has matrix elements with both orbitals i and j so
c i and j have same symmetry.
c
                  lbadsm = .false.
                  goto 331
                end if
              end if
  330       continue
          else
c
c*** i and j have a nonvanishing matrix element and belong to same symmetry
c
            lbadsm = .false.
          end if
c
c*** if i and j belong to same symmetry then give i same sym. number as j
c
  331     if (.not.lbadsm) then
            iorbsm(i) = iorbsm(j)
          end if
  320   continue
        if (iorbsm(i).eq.0) then
          nsym = nsym+1
          iorbsm(i) = nsym
        end if
  310 continue
c
      write(6,'(//'' orbital symmetry determined from H-matrix'',/)')
      do 350 i = 1,nsym
        k = 0
        do 360 j = 1,norb
          if (iorbsm(j).eq.i) then
            k = k + 1
            lll(k) = j
          end if
  360   continue
        write(6,'(i3,10(''   '',20i4))')i,(lll(j),j=1,k)
  350 continue
c
      return
      end
