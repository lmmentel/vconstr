      program dsfun
c***********************************************************************
c Calculate Kohn-Sham density functional
c energy and wavefunction in basis of HF-orbitals.
c
c The one-electron h-matrix is restored from dumpfile (from splice-
c section). The Coulomb repulsion matrix is calculated from the two-
c electron integrals and the Exchange (rho**1/3) matrix is calculated
c numerically.
c The points and weights that are needed in the numerical intgration are
c read from file POINTS
c
c The orbital occupations are given as input and can have fractional values.
c
c It is possible to utilize the symmetry of the HF (basis) orbitals. The
c symmetry of the HF orbs is then determined, using either the two-electron
c integrals or the h-matrix, and after each cycle the HFS orbitals are
c "cleaned", that is symmetry contaminating contributions are removed.
c The final symmetry of the KS orbs is then equal to the symmetry of the
c HF orbs.
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(nmomx  = 60)
      parameter(infmx  = 10)
      parameter(eps = 1.d-6)
c
      character*4 ied,iednum 
      character*44 fname
      common/discc/ied(16),fname(16)
c
      character*8 smtype,hmat,int2e
      logical lsym,lintsm,lrdocc,lrfun,lfield,atmol4
      dimension occmo(nmomx),info(infmx)
      dimension fxyz(3)
c
      common /titel/title(10)
c
      data hmat,int2e/'h-matrix','2-elec'/
c
      call prep99
      call tidajt(date,time,accno,anam,idum)
      write(6,666)anam,date,time,accno
  666 format(1h1//
     *24x,18('*'),' dsfun ',18('*') //
     *24x,'job name   ',a8//
     *24x,'date       ',a8//
     *24x,'time       ',a8//
     *24x,'acct       ',a8//
     *24x,44('*')///)
c
      nmos = -1
      nppr = 0
      nvpr = 0
      itrx  = 75
      kpens = 0
c
      idmp = 2
      ismo = 1
      isno = 2
      isks = 3
      isao = -1
      isoe = 198
      iint = 7
c
      lsym   = .false.
      lintsm = .false.
      lrdocc = .true.
      lrfun  = .false.
      lfield = .false.
c
      df     = 0.25d0
c      df=0.5d0
      dvdmp  = 0.7d0
      scfdmp = 3.d-1
      tstthr = eps
      thresh = eps
      crrmn  = 0.95d0
      crrmx  = 1.05d0
      rnel   = 0.d0
c
      info(1:infmx)=-1
      occmo(1:nmomx)=0.d0
      fxyz(1:3)=0.d0
c
      write(6,'(''  Selfconsistent Kohn-Sham calculation '')')
      call givtim(top,bot)
      write(6,'(/''****input read at'',f13.3,'' wall'',f13.3,
     + '' secs''/)')top,bot
      goto 55
c
   50 call input
   55 call passdf(ii)
      if (ii.le.0) then
         call cinput(jrec,jump,0)
         call erroutf(jrec)
         call caserr(' password not recognised ')
      endif
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)ii
c
c*** title
    1 read(5,'(10a8)') title
      goto 50
c
c*** nmos
    2 call inpi(nmos)
      if (nmos.gt.nmomx) then
        write(6,'(''ERROR; Too many Kohn-Sham orbitals specified; '',
     + i4,'' > '',i4)')nmos,nmomx
        stop
      endif
      lrdocc=.false.
      do imo=1,nmos
        occmo(imo)=2.d0
      enddo
      goto 50
c
c*** occ
    3 if (nmos.le.0) then
        write(6,'(''ERROR; nmos not yet defined'')')
        stop
      endif
      do i=1,nmos
        call inpf(occmo(i))
        rnel=rnel+occmo(i)
      enddo
      if (abs(modulo(rnel,1.d0)).gt.eps) then
        write(6,'(''ERROR; Number of electrons incorrect'')')
        write(6,'('' rnel = '',f6.2)')rnel
        stop
      endif
      goto 50
c
c*** conv
    4 call inpi(itrx)
      call inpf(tstthr)
      goto 50
c
c*** shift
    5 call inpf(df)
      goto 50
c
c*** vcrrmx
    6 call inpf(crrmn)
      crrmn=abs(crrmn)
      call inpf(crrmx)
      crrmx=abs(crrmx)
      if (crrmn.ge.1.d0) then
        write(6,'(/,''ERROR; crrmn must be less than 1.0'')')
        crrmn = 0.95d0
      endif
      if (crrmx.le.1.d0) then
        write(6,'(/,''ERROR; crrmx must be greater than 1.0'')')
        crrmx = 1.05d0
      endif
      goto 50
c
c*** symdet
    7 lsym = .true.
      call inpa(smtype)
      if (smtype.eq.int2e) then
        lintsm = .true.
      elseif (smtype.ne.hmat) then
        write(6,'(/,''WARNING; symmetry determination not '',
     + ''properly specified'')')
      endif
      call inpf(thresh)
      goto 50
c
c*** vxalph
    8 call inpf(alpha)
      call inpf(beta)
      call inpf(gamma)
      goto 50
c
c*** pprint
    9 call inpi(nppr)
      if (nppr.gt.infmx) then
        write(6,'(''ERROR; Requested too many intermediate storages ''
     +,''of the potential ; '',i4 ,'' > '',i4)') nppr,infmx
        stop
      endif
      do i=1,nppr
        call inpi(info(i))
      enddo
      goto 50
c
c*** vprint 
   10 call inpi(nvpr)
      goto 50
c
c*** dmpfile
   11 call inpi(idmp)
      call inpi(ismo)
      if ((ismo.lt.1).and.(ismo.gt.190)) then
        write(6,'(''ERROR; mo section not properly specified'')')
        stop
      endif
      call inpi(isno)
      if ((isno.lt.1).and.(isno.gt.190)) then
        write(6,'(''ERROR; no section not properly specified'')')
        stop
      endif
      call inpi(isoe)
      if ((isoe.lt.1).and.(isoe.gt.190)) then
        write(6,'(''ERROR; splice section not properly specified'')')
        stop
      endif
      goto 50
c
c*** intfile
   12 call inpi(iint)
      goto 50
c
c*** ksdmp
   13 call inpi(iks)
      if ((isks.lt.1).and.(isks.gt.190)) then
        write(6,'(''ERROR; ks section not properly specified'')')
        stop
      endif
      goto 50
c
c** adapt
   14 call inpi(isao)
      if ((isao.lt.1).and.(isao.gt.190)) then
        write(6,'(''ERROR; isao not properly specified'')')
        stop
      endif
      goto 50
c
c** ensdet
   15 call inpi(ibcens)
      call inpi(iscens)
      call inpi(kpens)
      if (kpens.lt.0) then
        kpens=-1
      else
        kpens=1
      endif
      call inpf(dqmax)
      if ((dqmax.le.0.d0).or.(dqmax.ge.2.d0)) then
        write(6,'(/''WARNING; dqmax > 0.3'')')
        write(6,'(''  default value taken instead'')')
        dqmax=1.d-1
      endif
      goto 50
c
c*** dipole
   16 lfield=.true.
      do i=1,3
        call inpf(fxyz(i))
      enddo
      goto 50
c
c** lrfun
   17 lrfun=.true.
      call inpf(dvdmp)
      goto 50
c
c** damping
   18 call inpf(scfdmp)
      goto 50
c
c*** enter
   19 call getgss(idmp,norb,atmol4)
      write(6,'(//''***** title : '',10a8,''*****'')')title
      write(6,'(/''  norb = '',i3)')norb
      if (nmos.gt.0) then
        if (nmos.gt.nmomx) then
          write(6,'(''ERROR; nmos > nmomx'')')
          write(6,'(''  nmos  = '',i3,/,''  nmomx = '',i3)')
     + nmos,nmomx
          stop
        endif
        if (norb.lt.nmos) then
          write(6,'(''ERROR; norb < nmos'')')
          write(6,'(''  nmos  = '',i3,/,''  norb = '',i3)')
     + nmos,norb
          stop
        endif
        write(6,'(''  nmos = '',i2)')nmos
        write(6,'(''  orbitals   : '',3(15(i2,''   ''),/,''     ''))')
     + (i,i=1,nmos)
        write(6,'(''  occupation : '',3(15f5.2,/,''     ''))')
     +(occmo(i),i=1,nmos)
      endif
      if (nvpr.lt.nmos) nvpr = nmos
      ispi=0
      do i=1,nmos
      if(occmo(i).lt.1.5d0) ispi=1
      enddo
c
      open(90,file='vntial.dat',status='old',err=222)
      write(6,'(/'' Initial potential Vxc read from file'')')
      close(90)
      goto 223
  222 write(6,'(/'' Initial potential Vxc = '',f6.3,
     + '' X-alpha + '',f6.3,'' Vosko-Wilk-Nusair + '',f6.3,
     + '' Becke '')')alpha,beta,gamma
c
  223 if (lsym) then
        if (lintsm) then
          write(6,'(/'' Symmetry determined from orbitals'')')
        else
          write(6,'(/'' Symmetry determined from one-electron '',
     + ''h-matrix'')')
        endif
      endif
      write(6,'(/'' Convergence'',/,''  itrx = '',i4,/,
     + ''  threshold = '',g8.2)')itrx,tstthr
      if (lrfun) then
        write(6,'(/'' Linear response procedure''/)')
        write(6,'(''  damping = '',f5.2)')dvdmp
      else
        write(6,'(''  shift   = '',f5.2)')df
        write(6,'(''  bounds  = '',f5.2,'','',f5.2)')crrmn,crrmx
      endif
      if (scfdmp.lt.1.d0) write(6,'(''  density damping : '',f5.2)')
     + scfdmp
      if (kpens.ne.0) then
        write(6,'(/'' Calculate ensemble from itr = '',i4,
     + ''   once every  '',i4,'' iterations'')')ibcens,iscens
        write(6,'(''   kpens = '',i4)')kpens
        write(6,'(''   dqmax = '',f4.2)')dqmax
      endif
      if (nppr.ne.0) write(6,'(/''  Intermediate storage of Vxc '',
     + ''required at itr = '',10(i4,'' ''))')(info(i),i=1,nppr)
      iednum=ied(idmp)
      write(6,'(/'' Summary of dumpfile  '',a4,'' sections''/,
     + ''  Molecular Orbital basis :'',i3,/,
     + ''  Natural Orbitals        :'',i3,/,
     + ''  one-electron h-matrix   :'',i3)')iednum,ismo,isno,isoe
      if (isao.gt.0) then
        write(6,'(''  Symmetry adaptation     :'',i3)')isao
      endif
      if (isks.gt.0) then
        write(6,'(''  Storage of KS orbitals  :'',i3)')isks
      endif
      iednum=ied(iint)
      write(6,'(/'' Integrals on file  '',a4)')iednum
c
      open(99,file='points',status='old',err=444)
      rewind(99)
      read(99,*)npnt,npold
      write(6,'(/,'' Number of gridpoints in numerical'',
     + '' integration :'',i5)')npold-1
      write(6,'('' Number of dummy points :'',i5)')npnt+1-npold
c
      if (lfield) write(6,'(/,'' Electric field : '',f6.4,
     + ''x + '',f6.4,''y + '',f6.4,''z'')')(fxyz(i),i=1,3)

c      if (ispi.gt.0) then
c      nmosa=2
c      nmosb=0
c     call brains(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmn,crrmx,
c     + thresh,dqmax,dvdmp,scfdmp,kpens,info,nppr,nvpr,npnt,npold,
c     + norb,nmos,nmomx,itrx,ibcens,iscens,iint,idmp,ismo,isno,isao,
c     + isoe,isks,lsym,lintsm,lrdocc,lrfun,nmosa,nmosb,atmol4)
c      else
      nmos=2
      call dbrain(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmn,crrmx,
     + thresh,dqmax,dvdmp,scfdmp,kpens,info,nppr,nvpr,npnt,npold,
     + norb,nmos,nmomx,itrx,ibcens,iscens,iint,idmp,ismo,isno,isao,
     + isoe,isks,lsym,lintsm,lrdocc,lrfun,atmol4) 
c      endif

      call revind
      call secsum
      call whtps
      call givtim(top,bot)
      write(6,10101)top,bot
10101 format(//' end of dnstyfun at',f13.3,' wall',f13.3,' secs')
      stop
c
  444 write(6,'(''ERROR; file "points" does not exist'')')
      stop
c
      end
c
c***********************************************************************
c
      subroutine passdf(ip)
c     get data field of a card and see what it is
c     ip =  0  : not recognised
c     ip >  0  : it was item ip of the list
      implicit real*8(a-h,o-z)
      parameter (ncnt=19)
      character*8 cnt(ncnt),word
      data cnt/ 'title','nmos','occ','conv','shift','vcrrmx','symdet',
     * 'vxalph','pprint','vprint','dmpfile','intfile','ksdmp','adapt',
     * 'ensdet','dipole','lrfun','damping','enter' /
c
      call reada(word)
      ip=locatc(cnt,ncnt,word)
      return
      end
c
c***********************************************************************
c
      block data prg
      implicit real *8   (a-h,p-w), integer (i-n), logical (o)
      common /prgnam/prgtit
      data prgtit/'dsfun'/
      end
