      subroutine cnwhmx(nhmx,nprmx,ngrpx)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq
c
      nhmx  = 0
      nprmx = 0
      ngrpx = ngrp
      do igrp = 1,ngrp
        lqnti = lquant(igrp)
        nitrm = ntrm(igrp)
        do jgrp = 1,igrp
          lqntj = lquant(jgrp)
          njtrm = ntrm(jgrp)
          nlm = 0
          do la = 0,lqnti
            do lb = 0,lqntj
              lalb = la+lb
              nlm=nlm+(lalb+1)*(lalb+1)
            enddo
          enddo
          nhmx = nhmx+nitrm*njtrm*nlm
          nprmx=nprmx+nitrm*njtrm
        enddo
      enddo
c
      return
      end
c
      subroutine cnvhmx(nhmx,nprmx,ngrpx)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/gaus3/zeta(10,160),coef(10,160),cx(100),cy(100),cz(100),
     1 zz(100),itype(160),kad(160),icentr(160),mbase(160),ntrm(160),
     2 title(10),acc1,acc2,imin(4),imax(4),nbasis,natoms,iblock,maxblk,
     3 irun,ngroup,ncentg,nosym,cgx(12),cgy(12),cgz(12),idddd,norbd,
     4 nnbase(160)      
c
      nhmx  = 0
      nprmx = 0
      ngrpx = ngroup
      do igrp = 1,ngroup
        lqnti = itype(igrp)-1
        nitrm = ntrm(igrp)
        do jgrp = 1,igrp
          lqntj = itype(jgrp)-1
          njtrm = ntrm(jgrp)
          nlm = 0
          do la = 0,lqnti
            do lb = 0,lqntj
              lalb = la+lb
              nlm=nlm+(lalb+1)*(lalb+1)
            enddo
          enddo
          nhmx = nhmx+nitrm*njtrm*nlm
          nprmx = nprmx+nitrm*njtrm
        enddo
      enddo
c
      return
      end

