      subroutine cond(iunit,pcond,norb,rnel,valmo)
c *******************************************************************
c  calculate hole density matrix from coef-matrix :
c
c        pcond(x1|x2) = coef(x1,x2) / p1(x1)
c
c  the density p1(x1) is calculated from coef(x1,x2)
c
c  reference electron in point 1, coordinates (x,y,z)
c
c  coef-matrix on unit "iunit" (ed0=1,ed1=2 etc), pcond matrix in array 
c  pcond (in triangular form, pcond(i*(i-1)/2+j) = sum of off-diagonal 
c  elements (i,j) and (j,i) ). pcond matrix is calculated in same basis 
c  as coef-matrix (usually mo's).
c
c  normalization :  coef normalized to n(n-1)
c                   p1 normalized to n
c                   pcond normalized to n-1
c
c  further arguments :  norb  : highest orbital number on coef-file
c                       rnel  : number of electrons
c                       valmo : array to hold mo-values in ref. pos.
c *******************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(explim = 300.0d0)
      dimension pcond(norb*(norb+1)/2),valmo(norb)
c
c*** common/atbuf/ used to store blocks of coef-file.
c
      common/atbuf/coef(340),rijkl(170),nword,ndum
      common/scra/ijkl(340,4)
c
      tol = exp(-explim)
c
      ijmx = norb*(norb+1)/2
c
      p1ex = 0.d0
      pcond(1:norb*(norb+1)/2)=0.d0 
c
c*** find startblock of coef-file, and read blocks
c
      call search(1,iunit)
  100 call fget(coef(1),nw,iunit)
      if (nw.ne.0) then
c
c*** unpack ijkl-indices
c
        call upak8v(rijkl(1),ijkl(1,1))
        do m = 1,nword
          i = ijkl(m,1)
          j = ijkl(m,2)
          k = ijkl(m,3)
          l = ijkl(m,4)
          ij = i*(i-1)/2+j
          if ((ij.le.0).or.(ij.gt.ijmx)) then
            write(6,'('' ERROR; ij too big:'',i16,
     + '' ijmx '',i16,2x,i12,2x,i12)')ij,ijmx,i,j
            stop
          endif
          kl = k*(k-1)/2+l
          if ((kl.le.0).or.(kl.gt.ijmx)) then
            write(6,'('' ERROR; kl too big:'',i16,
     + '' ijmx '',i16,2x,i12,2x,i12)')kl,ijmx,k,l
            stop
          endif 
          valij = coef(m)*valmo(i)*valmo(j)
          valkl = coef(m)*valmo(k)*valmo(l)
          if (k.eq.l) p1ex = p1ex+valij
          if (i.eq.j) p1ex = p1ex+valkl
          pcond(ij) = pcond(ij)+valkl
          pcond(kl) = pcond(kl)+valij
        enddo
        goto 100
      endif
      p1ex = p1ex/(rnel-1.d0)
c
c*** normalizing and dividing by exact density in reference position 1
c
      if (p1ex.lt.tol) then
        write(6,'(''WARNING; small density in reference position'')')
        write(6,'(''  rho : '',g16.8)')p1ex
      endif
      do iorb = 1,norb*(norb+1)/2
        pcond(iorb) = pcond(iorb)/p1ex
      enddo
c
      cnel=0.d0
      do iorb=1,norb
        cnel=cnel+pcond(iorb*(iorb+1)/2)
      enddo
      if (abs(cnel+1.d0-rnel).gt.1.d-6) then
        write(6,'(''ERROR; cnel = '',f16.10,'' rnel = '',f16.10,
     + '' diff = '',g16.8)') cnel,rnel,cnel+1.d0-rnel
        stop
      endif
c
      return
      end
