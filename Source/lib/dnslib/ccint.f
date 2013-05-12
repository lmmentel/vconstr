      subroutine ccintw(valint,norbmx,ncomp,ncbas,lqntmx,ll,mm,nn,
     + ncmpx,hx,hy,hz,nhmx,xpnt,ypnt,zpnt,prefac,nprmx,l1l2mx,ngrpx)
c***********************************************************************
c calculate one-electron integrals of the electron-nuclear attraction
c type in a basis of cgtf's.  the 'nucleus' is replaced by a negative
c unit point charge at position (xpnt,ypnt,zpnt)
c
c use this routine for gaussians that are restored from an 'integw'
c dumpfile.
c
c routine calchw must have been executed (filling arrays hx,hy and hz)
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension valint(norbmx*(norbmx+1)/2),ncomp(0:lqntmx),
     + ncbas(0:lqntmx),ll(ncmpx),mm(ncmpx),nn(ncmpx),
     + hx(nhmx),hy(nhmx),hz(nhmx),prefac(nprmx)
      dimension nbsgrp(ngrpx),xint(15,15),fnutau(0:l1l2mx),
     + xi(0:l1l2mx,0:lqntmx,0:lqntmx),yj(0:l1l2mx,0:lqntmx,0:lqntmx),
     + zk(0:l1l2mx,0:lqntmx,0:lqntmx)
c
      common/gaus4/ctran(1300),zeta(1300),ilifc(340),ntrm(340),
     1 lquant(340),nquant(340),icen3(340),iext(340),iflaga(340),
     2 iflagb(340),icont(340),norb(340),ngrp,nbasis,title(10),
     3 cart(600),zz(200),tag(200),ncen,ndum,acc1,acc2,acc3,
     4 radius,radsq
c
      nhstrt = 1
      nbsgrp(1) = 0
c
      do i = 2,ngrp
        nbsgrp(i) = nbsgrp(i-1)+ncomp(lquant(i-1))
      enddo
c
c** loop over groups
c
      iprf = 0
      do 100 igrp = 1,ngrp
        lqnti = lquant(igrp)
        nicomp = ncomp(lqnti)
        nitrm = ntrm(igrp)
        ax = cart(icen3(igrp)+1)
        ay = cart(icen3(igrp)+2)
        az = cart(icen3(igrp)+3)
c
        do 90 jgrp = 1,igrp
          lqntj = lquant(jgrp)
          njcomp = ncomp(lqntj)
          njtrm = ntrm(jgrp)
          bx = cart(icen3(jgrp)+1)
          by = cart(icen3(jgrp)+2)
          bz = cart(icen3(jgrp)+3)
c
          lmnsum = lqnti+lqntj
          xint(1:nicomp,1:njcomp)=0.d0
c
c** loop over gaussians in contraction
c
      do 80 it = 1,nitrm
        alpha = zeta(ilifc(igrp)+it)
c
        do 70 jt = 1,njtrm
          alphb = zeta(ilifc(jgrp)+jt)
          somexp = alpha+alphb
          smexpi = 1.d0/somexp
          px = (alpha*ax+alphb*bx)*smexpi
          py = (alpha*ay+alphb*by)*smexpi
          pz = (alpha*az+alphb*bz)*smexpi
c
          iprf = iprf+1
          prf = prefac(iprf)
c
          cpx = xpnt-px
          cpy = ypnt-py
          cpz = zpnt-pz
          cp2 = cpx**2+cpy**2+cpz**2
          tau = somexp*cp2
c
c** calculate possible values for f(nu,tau)
c
      do nu = 0,lmnsum
        fnutau(nu) = fnusp(nu,tau)
      enddo
c
      do la = 0,lqnti
        do lb = 0,lqntj
          lalb = la+lb
          call xyzfac(hx(nhstrt),lalb,cpx,xi(0,la,lb))
          call xyzfac(hy(nhstrt),lalb,cpy,yj(0,la,lb))
          call xyzfac(hz(nhstrt),lalb,cpz,zk(0,la,lb))
          nhstrt = nhstrt+(lalb+1)*(lalb+1)
        enddo
      enddo
c
c** loop over components in group
c
      do 60 ic = 1,nicomp
        numi = nbsgrp(igrp)+ic
        li = ll(ncbas(lqnti)+ic)
        mi = mm(ncbas(lqnti)+ic)
        ni = nn(ncbas(lqnti)+ic)
c
        do 50 jc = 1,njcomp
          numj = nbsgrp(jgrp)+jc
c
          if (numj.le.numi) then
            lj = ll(ncbas(lqntj)+jc)
            mj = mm(ncbas(lqntj)+jc)
            nj = nn(ncbas(lqntj)+jc)
            xxx = 0.d0
            do i = 0,li+lj
              do j = 0,mi+mj
                do k = 0,ni+nj
                  xxx = xxx+xi(i,li,lj)*yj(j,mi,mj)*
     + zk(k,ni,nj)*fnutau(i+j+k)
                enddo
              enddo
            enddo
            xint(ic,jc) = xint(ic,jc)+xxx*prf
            if (igrp.eq.jgrp) xint(jc,ic) = xint(ic,jc)
          endif
c
   50   continue
   60 continue
c
   70   continue
   80 continue
c
c** transform integrals in xint from cartesian to sperical harmonic
c** basis and store transformed integrals in array valint
c
      call carsph(xint(1,1),lqnti,lqntj,ncomp(0),lqntmx)
c
      do i = 1,2*lqnti+1
        do j = 1,2*lqntj+1
          iorb = norb(igrp)+i
          jorb = norb(jgrp)+j
          if (jorb.le.iorb) then
            ij = iorb*(iorb-1)/2+jorb
            valint(ij) = xint(i,j)
          endif
        enddo
      enddo
c
   90   continue
  100 continue
c
      return
      end
c
      subroutine ccintv(valint,norbmx,ncomp,ncbas,lqntmx,ll,mm,nn,
     + ncmpx,hx,hy,hz,nhmx,xpnt,ypnt,zpnt,prefac,nprmx,l1l2mx,ngrpx)
c***********************************************************************
c calculate one-electron integrals of the electron-nuclear attraction
c type in a basis of cgtf's.  the 'nucleus' is replaced by a negative
c unit point charge at position (xpnt,ypnt,zpnt)
c
c use this routine for gaussians that are restored from an 'integv'
c dumpfile.
c
c routine calchv must have been executed (filling arrays hx,hy and hz)
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension valint(norbmx*(norbmx+1)/2),ncomp(0:lqntmx),
     + ncbas(0:lqntmx),ll(ncmpx),mm(ncmpx),nn(ncmpx),
     + hx(nhmx),hy(nhmx),hz(nhmx),prefac(nprmx)
      dimension nbsgrp(ngrpx),xint(15,15),fnutau(0:l1l2mx),
     + xi(0:l1l2mx,0:lqntmx,0:lqntmx),yj(0:l1l2mx,0:lqntmx,0:lqntmx),
     + zk(0:l1l2mx,0:lqntmx,0:lqntmx)
c
      common/gaus3/zeta(10,160),coef(10,160),cx(100),cy(100),cz(100),
     1 zz(100),itype(160),kad(160),icentr(160),mbase(160),ntrm(160),
     2 title(10),acc1,acc2,imin(4),imax(4),nbasis,natoms,iblock,maxblk,
     3 irun,ngroup,ncentg,nosym,cgx(12),cgy(12),cgz(12),idddd,norbd,
     4 nnbase(160)      
c
      nhstrt = 1
      nbsgrp(1) = 0
      do i = 2,ngroup
        nbsgrp(i) = nbsgrp(i-1)+ncomp(itype(i-1)-1)
      enddo
c
c** loop over groups
c
      iprf = 0
      do 100 igrp = 1,ngroup
        lqnti = itype(igrp)-1
        nicomp = ncomp(lqnti)
        nitrm = ntrm(igrp)
        ax = cx(icentr(igrp))
        ay = cy(icentr(igrp))
        az = cz(icentr(igrp))
c
        do 90 jgrp = 1,igrp
          lqntj = itype(jgrp)-1
          njcomp = ncomp(lqntj)
          njtrm = ntrm(jgrp)
          bx = cx(icentr(jgrp))
          by = cy(icentr(jgrp))
          bz = cz(icentr(jgrp))
c
          lmnsum = lqnti+lqntj
          xint(1:nicomp,1:njcomp)=0.d0
c
c** loop over gaussians in contraction
c
      do 80 it = 1,nitrm
        alpha = zeta(it,igrp)
c
        do 70 jt = 1,njtrm
          alphb = zeta(jt,jgrp)
          somexp = alpha+alphb
          smexpi = 1.d0/somexp
          px = (alpha*ax+alphb*bx)*smexpi
          py = (alpha*ay+alphb*by)*smexpi
          pz = (alpha*az+alphb*bz)*smexpi
c
          iprf = iprf+1
          prf = prefac(iprf)
c
          cpx = xpnt-px
          cpy = ypnt-py
          cpz = zpnt-pz
          cp2 = cpx*cpx+cpy*cpy+cpz*cpz
          tau = somexp*cp2
c
c** calculate possible values for f(nu,tau)
c
      do nu = 0,lmnsum
        fnutau(nu) = fnusp(nu,tau)
      enddo
c
      do la = 0,lqnti
        do lb = 0,lqntj
          lalb = la + lb
          call xyzfac(hx(nhstrt),lalb,cpx,xi(0,la,lb))
          call xyzfac(hy(nhstrt),lalb,cpy,yj(0,la,lb))
          call xyzfac(hz(nhstrt),lalb,cpz,zk(0,la,lb))
          nhstrt = nhstrt+(lalb+1)*(lalb+1)
        enddo
      enddo
c
c** loop over components in group
c
      do 60 ic = 1,nicomp
        numi = nbsgrp(igrp)+ic
        li = ll(ncbas(lqnti)+ic)
        mi = mm(ncbas(lqnti)+ic)
        ni = nn(ncbas(lqnti)+ic)
c
        do 50 jc = 1,njcomp
          numj = nbsgrp(jgrp) + jc
c
          if (numj.le.numi) then
            lj = ll(ncbas(lqntj)+jc)
            mj = mm(ncbas(lqntj)+jc)
            nj = nn(ncbas(lqntj)+jc)
            xxx = 0.d0
            do i = 0,li+lj
              do j = 0,mi+mj
                do k = 0,ni+nj
                  xxx = xxx+xi(i,li,lj)*yj(j,mi,mj)*
     + zk(k,ni,nj)*fnutau(i+j+k)
                enddo
              enddo
            enddo
            xint(ic,jc) = xint(ic,jc)+xxx*prf
            if (igrp.eq.jgrp) xint(jc,ic) = xint(ic,jc)
          endif
c
   50   continue
   60 continue
c
   70   continue
   80 continue
c
c** in case of d-orbitals, normalize and (if necessary) transform
c** from 6-d set to 5-d set.
c
      if ( (lqnti.eq.2).or.(lqntj.eq.2)) then
        if (idddd.eq.0) then
          call cart6d(xint(1,1),ncomp(0),lqi,lqj,lqntmx)
        else
          call carsph(xint(1,1),lqnti,lqntj,ncomp(0),lqntmx)
        endif
      endif
c
      if ((idddd.eq.0).and.(lqnti.eq.2))  then
c** 6-d set is used
        iimax = 6
      else
        iimax = 2*lqnti+1
      endif
      if ((idddd.eq.0).and.(lqntj.eq.2))  then
c** 6-d set is used
        jmax = 6
      else
        jmax = 2*lqntj+1
      endif
c
      if (idddd.eq.0) then
        nbasi = mbase(igrp)
        nbasj = mbase(jgrp)
      else
        nbasi = nnbase(igrp)
        nbasj = nnbase(jgrp)
      endif
c
      do i = 1,iimax
        do j = 1,jmax
          iorb = nbasi+i
          jorb = nbasj+j
          if (jorb.le.iorb) then
            ij = iorb*(iorb-1)/2+jorb
            valint(ij) = xint(i,j)
          endif
        enddo
      enddo
c
   90 continue
  100 continue
c
      return
      end
