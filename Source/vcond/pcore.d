      subroutine p2core(iunp2,iblp2,p,norb,iunnew,iblnew,ncore)
c***********************************************************************
c expand the p2-matrix with ncore core functions.
c see comment in routines core and corval
c
c arguments:
c iunp2, iblp2  : unit-number and startblock of p2 file
c p, norb       : density matrix and dimension
c iunnew, iblnew: unit-number and startblock of expanded p2-file
c ncore         : number of core orbitals
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension p(norb*(norb+1)/2)
c
      call search(iblnew,iunnew)
c
c** core part
c
      call core(iunnew,ncore)
c
c** core valance part
c
      call corval(iunp2,iblp2,iunnew,p(1),norb,ncore)
c
      return
      end
c
      subroutine core(iunit2,ncore)
c********************************************************************
c This routine and routine corval are both called from p2core to expand
c the p2-matrix with ncore core orbitals.
c
c In this routine the hf-like core part is calculated while the core-
c valence part (which makes use of the one-density) is calculated in
c corval.
c
c the p2-elements are stored on atmol file with unit number iunit2.
c********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(one  = 1.0d0)
      parameter(two  = 2.0d0)
      parameter(four = 4.0d0)
c
c*** common/atbuf/ used to store blocks of p2-file.
c
      common/atbuf2/coef2(340),rijkl2(170),icoef,ndum
      common/scra/ijkl(340,4)
c
      icoef = 0
c
c*** loop over core orbitals
c
      do 100 i = 1,ncore
c
c*** first exchange elements [i 1,i 1], [i 2,i 2], ...[i i-1,i i-1]
c
        do 10 j = 1,i-1
          icoef = icoef+1
c
c*** pack indices
c
          ijkl(icoef,1) = i
          ijkl(icoef,2) = j
          ijkl(icoef,3) = i
          ijkl(icoef,4) = j
          coef2(icoef) = -two
          if (icoef.eq.340) then
            call pak8x(ijkl(1,1),rijkl(1))
            call put(coef2(1),511,iunit2)
            icoef = 0
          endif
   10   continue
c
c*** now coulomb elements [i i,1 1], [i i,2 2] ... [i i,i i]
c
        do 20 j = 1,i
          icoef = icoef+1
c
c*** pack indices
c
          ijkl(icoef,1) = i
          ijkl(icoef,2) = i 
          ijkl(icoef,3) = j
          ijkl(icoef,4) = j
          if (i.eq.j) then
            coef2(icoef) = one
          else
            coef2(icoef) = four
          endif
          if (icoef.eq.340) then
            call put(coef2(1),511,iunit2)
            icoef = 0
          endif
   20   continue
  100 continue
c
      return
      end
c
      subroutine corval(iunit1,iblk1,iunit2,p,norb,ncore)
      implicit real*8 (a-h,o-z),integer(i-n)
      logical lfp2,odd
      dimension p(*)
c
      common/atbuf/coef1(340),rijkl1(170),nword1,ndum1
      common/atbuf2/coef2(340),rijkl2(170),nword2,ndum2
      common/scra/ijkl(340,4)
c
      parameter(one = 1.0d0) 
      parameter(two = 2.0d0)
c
      call search(iblk1,iunit1)
      call find(iunit1)
      call get(coef1(1),nw)
      if (nw.ne.0) then
c
c*** unpack ijkl-indices
c
        call unpack(rijkl1(1),8,ijkl(1,1),nword1*4)
        m1 = 1
        ll1 = ijkl(1,4)+ncore
        kk1 = ijkl(1,3)+ncore
        jj1 = ijkl(1,2)+ncore
        ii1 = ijkl(1,1)+ncore
        ijkl1 = ior(ior(ishft(ii1,24),ishft(jj1,16)),
     1 ior(ishft(kk1,8),ll1))
        lfp2 = .false.
c
c*** loop
c
        do 500 ip = 1,norb
          ipc = ip+ncore
c
c*** first exchange elements
c
      do 40 ii = 1,ncore
        do 50 iq = 1,ip
          iqc = iq+ncore
          val = -one*p(ip*(ip-1)/2+iq)
          if (abs(val).gt.1.0d-12) then
            ijklx = ior(ior(ishft(ipc,24),ishft(ii,16)),
     1 ior(ishft(iqc,8),ii))
c
            nword2 = nword2+1
            if ((ijklx.gt.ijkl1)) then
              coef2(nword2) = coef1(m1)
              rijkl2(nword2) = rijkl1(m1)
              idx = (nword2+1)/2
              if (m1.eq.nword1) then
                call find(iunit1)
                call get(coef1(1),nw)
                if (nw.eq.0) then
                  lfp2 = .true.
                else
                  call unpack(rijkl1(1),8,ijkl(1,1),nword1*4)
                endif
                m1 = 0
              endif
         if (lfp2) goto 170
c
         m1 = m1 + 1
         ll1 = ijkl(m1,4) + ncore
         kk1 = ijkl(m1,3) + ncore
         jj1 = ijkl(m1,2) + ncore
         ii1 = ijkl(m1,1) + ncore
         ijkl1 = ior(ior(ishft(ii1,24),ishft(jj1,16)),
     1               ior(ishft(kk1,8),ll1))
      goto 150
c
  150 if (nword2.eq.340) then
        call put(coef2(1),511,iunit2)
        nword2 = 0
      end if
c

c
  170    coef2(nword2) = val
         idx = (nword2+1)/2
         if (odd(nword2)) then
           rijkl2(idx) = pack2(ijklx,0)
         else
           rijkl2(idx) = pack2(rijkl2(idx),ijklx)
         end if
   50 continue
c
   40 continue
c
c now coulomb elements
c
      do 90 iq = 1,ip
      iqc = iq+ncore
      val =  two*p(ip*(ip-1)/2+iq)
      if (abs(val).lt.1.0e-10) goto 90
c
      do 80 ii = 1,ncore
      ijklx = ior(ior(ishft(ipc,24),ishft(iqc,16)),
     1 ior(ishft(ii,8),ii))
c
  200 if (nword2.eq.340) then
        call put(coef2(1),511,iunit2)
        nword2 = 0
      end if
c
      nword2 = nword2 + 1
      if ((ijklx.lt.ijkl1).or.lfp2) goto 250
         coef2(nword2) = coef1(m1)
         idx = (nword2+1)/2
         if (odd(nword2)) then
           rijkl2(idx) = pack2(ijkl1,0)
         else
           rijkl2(idx) = pack2(rijkl2(idx),ijklx)
         end if
c
         if (m1.eq.nword1) then
            call find(iunit1)
            call get(coef1(1),nw)
            if (nw.eq.0) then
              lfp2 = .true.
            else
              call upak8v(rijkl1(1),ijkl(1,1),nword1)
            end if
            m1 = 0
         end if
         if (lfp2) goto 250
c
         m1 = m1 + 1
         ll1 = ijkl(m1,4) + ncore
         kk1 = ijkl(m1,3) + ncore
         jj1 = ijkl(m1,2) + ncore
         ii1 = ijkl(m1,1) + ncore
         ijkl1 = ior(ior(ishft(ii1,24),ishft(jj1,16)),
     1               ior(ishft(kk1,8),ll1))
      goto 200
c
  250    coef2(nword2) = val
         idx = (nword2+1)/2
         if (odd(nword2)) then
           rijkl2(idx) = pack2(ijklx,0)
         else
           rijkl2(idx) = pack2(rijkl2(idx),ijklx)
         end if
   80 continue
   90 continue
c
  100 continue
c
      if (lfp2) goto 400
  300 if (nword2.eq.340) then
        call put(coef2(1),511,iunit2)
        nword2 = 0
      end if
c
      nword2 = nword2 + 1
         coef2(nword2) = coef1(m1)
         idx = (nword2+1)/2
         if (odd(nword2)) then
           rijkl2(idx) = pack2(ijkl1,0)
         else
           rijkl2(idx) = pack2(rijkl2(idx),ijkl1)
         end if
c
         if (m1.eq.nword1) then
            call find(iunit1)
            call get(coef1(1),nw)
            if (nw.eq.0) then
              lfp2 = .true.
            else
              call upak8v(rijkl1(1),ijkl(1,1),nword1)
            end if
            m1 = 0
         end if
         if (lfp2) goto 400
c
         m1 = m1 + 1
         ll1 = ijkl(m1,4) + ncore
         kk1 = ijkl(m1,3) + ncore
         jj1 = ijkl(m1,2) + ncore
         ii1 = ijkl(m1,1) + ncore
         ijkl1 = ior(ior(ishft(ii1,24),ishft(jj1,16)),
     1               ior(ishft(kk1,8),ll1))
      goto 300
c
  400 call put(coef2(1),511,iunit2)
      call put(coef2(1),0,iunit2)
      return
      end
