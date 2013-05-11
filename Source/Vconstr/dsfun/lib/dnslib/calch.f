      subroutine calchw(hx,hy,hz,nhmx,prefac,nprmx,l1l2mx)
c***********************************************************************
c calculate hx, hy and hz for every gaussian-pair
c
c gaussians in common/gauss/ with integw-dumpfile layout
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(twopi = 6.28318530717959647692d0)
c
      dimension hx(nhmx),hy(nhmx),hz(nhmx),prefac(nprmx)
      dimension gl(0:l1l2mx)
c
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq
c
      nhstrt = 1
c
c** loop over groups
c
      iprf = 0
      do 100 igrp = 1,ngrp
        nitrm = ntrm(igrp)
        ax = cart(icen3(igrp)+1)
        ay = cart(icen3(igrp)+2)
        az = cart(icen3(igrp)+3)
c
        do 90 jgrp = 1,igrp
          njtrm = ntrm(jgrp)
          bx = cart(icen3(jgrp)+1)
          by = cart(icen3(jgrp)+2)
          bz = cart(icen3(jgrp)+3)
c
c** loop over gaussians in contractie
c
      do 80 it = 1,nitrm
        alpha = zeta(ilifc(igrp)+it)
c
        do 70 jt = 1,njtrm
          alphb = zeta(ilifc(jgrp)+jt)
          cff = ctran(ilifc(igrp)+it)*ctran(ilifc(jgrp)+jt)
          somexp = alpha+alphb
          smexpi = 1.d0/somexp
          eps = 1.d0/(4.d0*somexp)
          px = (alpha*ax+alphb*bx)/somexp
          py = (alpha*ay+alphb*by)/somexp
          pz = (alpha*az+alphb*bz)/somexp
          pax = px-ax
          pay = py-ay
          paz = pz-az
          pbx = px-bx
          pby = py-by
          pbz = pz-bz
          ab2 = (ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz)
          iprf = iprf+1
          prefac(iprf) = twopi*exp(-alpha*alphb*ab2*smexpi)*
     + cff*smexpi
c
      do la = 0,lquant(igrp)
        do lb = 0,lquant(jgrp)
          lalb = la+lb
          call g(la,lb,pax,pbx,eps,gl(0))
          call h(gl(0),eps,hx(nhstrt),lalb)
          call g(la,lb,pay,pby,eps,gl(0))
          call h(gl(0),eps,hy(nhstrt),lalb)
          call g(la,lb,paz,pbz,eps,gl(0))
          call h(gl(0),eps,hz(nhstrt),lalb)
          nhstrt = nhstrt+(lalb+1)*(lalb+1)
        enddo
      enddo
c
   70   continue
   80 continue
c
   90   continue
  100 continue 
c
      return
      end
c
      subroutine calchv(hx,hy,hz,nhmx,prefac,nprmx,l1l2mx)
c***********************************************************************
c calculate hx, hy and hz for every gaussian-pair
c
c gaussians in common/gauss/ with integv-dumpfile layout
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(twopi = 6.28318530717959647692d0)
c
      dimension hx(nhmx),hy(nhmx),hz(nhmx),prefac(nprmx)
      dimension gl(0:l1l2mx)
c
      common/gaus3/zeta(10,160),coef(10,160),cx(100),cy(100),cz(100),
     1 zz(100),itype(160),kad(160),icentr(160),mbase(160),ntrm(160),
     2 title(10),acc1,acc2,imin(4),imax(4),nbasis,natoms,iblock,maxblk,
     3 irun,ngroup,ncentg,nosym,cgx(12),cgy(12),cgz(12),idddd,norbd,
     4 nnbase(160)      
c
      nhstrt = 1
c
c** loop over groups
c
      iprf = 0
      do 100 igrp = 1,ngroup
        lqnti = itype(igrp)-1
        nitrm = ntrm(igrp)
        ax = cx(icentr(igrp))
        ay = cy(icentr(igrp))
        az = cz(icentr(igrp))
c
        do 90 jgrp = 1,igrp
          lqntj = itype(jgrp)-1
          njtrm = ntrm(jgrp)
          bx = cx(icentr(jgrp))
          by = cy(icentr(jgrp))
          bz = cz(icentr(jgrp))
c
c** loop over gaussians in contraction
c
      do 80 it = 1,nitrm
        alpha = zeta(it,igrp)
c
        do 70 jt = 1,njtrm
          cff = coef(jt,jgrp)*coef(it,igrp)
          alphb = zeta(jt,jgrp)
          somexp = alpha+alphb
          smexpi = 1.d0/somexp
          eps = 1.d0/(4.d0*somexp)
          px = (alpha*ax+alphb*bx)/somexp
          py = (alpha*ay+alphb*by)/somexp
          pz = (alpha*az+alphb*bz)/somexp
          pax = px-ax
          pay = py-ay
          paz = pz-az
          pbx = px-bx
          pby = py-by
          pbz = pz-bz
          ab2 = (ax-bx)**2+(ay-by)**2+(az-bz)**2
          iprf = iprf+1
          prefac(iprf) = twopi*exp(-alpha*alphb*ab2*smexpi)*
     + cff*smexpi
c
      do la = 0,lqnti
        do lb = 0,lqntj
          lalb = la + lb
          call g(la,lb,pax,pbx,eps,gl(0))
          call h(gl(0),eps,hx(nhstrt),lalb)
          call g(la,lb,pay,pby,eps,gl(0))
          call h(gl(0),eps,hy(nhstrt),lalb)
          call g(la,lb,paz,pbz,eps,gl(0))
          call h(gl(0),eps,hz(nhstrt),lalb)
          nhstrt = nhstrt+(lalb+1)*(lalb+1)
        enddo
      enddo
c
   70   continue
   80 continue
c
   90   continue
  100 continue
c
      return
      end
