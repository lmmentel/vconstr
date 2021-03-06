      program dsfun
       use basisModule 
       use commonsModule
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
c
      character*8 smtype,hmat,int2e
      logical lsym,lintsm,lrdocc,lrfun,lfield
      dimension occmo(nmomx),info(infmx)
      dimension fxyz(3)
c
clmm..stuff for input processing 
      character*100 buffer, inptf, gbasisfile, gdictnfile, gintegfile, 
     +     outksoFile
      character*100 gridfile
      character*50  title
      character*80  gamtitle
      namelist /input/ 
     & alpha, 
     & beta, 
     & fxyz,
     & gamma,
     & gbasisfile, 
     & gdictnfile, 
     & gintegfile, 
     & gridfile,
     & outksoFile,
     & iprint, 
     & lfield, 
     & lintsm, 
     & lrfun, 
     & lsym,
     & maxiters, 
     & nmos,
     & nppr, 
     & scfdamp, 
     & smtype, 
     & thresh, 
     & title,
     & tstthr 
      namelist /occupations/ occmo
      namelist /linresp/ df, dvdamp, crrmin, crrmax 
c
      gridfile='points'
      title = 'default '
      outksoFile = ''
      nmos = -1
      nppr = 0
      nvpr = 0
      maxiters  = 75
      kpens = 0
c
      iprint = 0
c
      lsym   = .false.
      lintsm = .false.
      lrdocc = .true.
      lrfun  = .false.
      lfield = .false.
c
      df     = 0.25d0
      dvdamp  = 0.7d0
      scfdamp = 3.d-1
      tstthr = epsilon
      thresh = epsilon
      crrmin  = 0.95d0
      crrmax  = 1.05d0
      rnel   = 0.d0
c
      info(1:infmx)=-1
      occmo(1:nmomx)=0.d0
      fxyz(1:3)=0.d0
c
      write(6,'(''  Selfconsistent Kohn-Sham calculation '')')
c
c.....lmm start.....commented this horrible input processing 
c
c*** vcrrmx
c    6 call inpf(crrmn)
c      crrmn=abs(crrmn)
c      call inpf(crrmx)
c      crrmx=abs(crrmx)
c      if (crrmn.ge.1.d0) then
c        write(6,'(/,''ERROR; crrmn must be less than 1.0'')')
c        crrmn = 0.95d0
c      endif
c      if (crrmx.le.1.d0) then
c        write(6,'(/,''ERROR; crrmx must be greater than 1.0'')')
c        crrmx = 1.05d0
c      endif
c      goto 50
c
c
c*** pprint
c    9 call inpi(nppr)
c      if (nppr.gt.infmx) then
c        write(6,'(''ERROR; Requested too many intermediate storages ''
c     +,''of the potential ; '',i4 ,'' > '',i4)') nppr,infmx
c        stop
c      endif
c      do i=1,nppr
c        call inpi(info(i))
c      enddo
c      goto 50
c
c** ensdet
c   15 call inpi(ibcens)
c      call inpi(iscens)
c      call inpi(kpens)
c      if (kpens.lt.0) then
c        kpens=-1
c      else
c        kpens=1
c      endif
c      call inpf(dqmax)
c      if ((dqmax.le.0.d0).or.(dqmax.ge.2.d0)) then
c        write(6,'(/''WARNING; dqmax > 0.3'')')
c        write(6,'(''  default value taken instead'')')
c        dqmax=1.d-1
c      endif
c      goto 50
c
clmm..get input file name form command line open it, read the 
clmm..contents and close
      call getarg(1, buffer)
      read(buffer, *) inptf
      open(unit=11, file=trim(inptf), status='old',form='formatted',
     &     delim='apostrophe', iostat=ios)
      read(11, nml=input)
      read(11, nml=occupations)
clmm..get linear response input if linear reponse is being done
      if (lrfun) read(11, nml=linresp)
      close(11)
clmm..rewrite some information to commons
      basisInfoFile  = gbasisfile
      dictionaryFile = gdictnfile
      integralsFile  = gintegfile
      printLevel     = iprint
      write(6,'(/,''grid file from nml: '',a)') gridfile
c.....lmm end of input reading 
clmm..
clmm..get the gaussian basis information from the basis file
clmm..
      call read_system_info(trim(gbasisfile), gamtitle, natoms, icharge, 
     & mult, nbf, nx, nele, na, nb, nshell, nprimi) 
clmm..not sure if norb should be set to the number of cartesian gaussians
clmm..or number of orbitals used in the calculation
      norb = nbf 
c.....nmos
      if (nmos.gt.-1) then
        if (nmos.gt.nmomx) then
          write(6,'(''ERROR; Too many Kohn-Sham orbitals specified; '',
     + i4,'' > '',i4)')nmos,nmomx
          stop
        endif
        lrdocc=.false.
        do imo=1,nmos
          occmo(imo)=2.d0
        enddo
      endif
c.....symdet
      if (smtype.eq.int2e) then
        lintsm = .true.
      elseif (smtype.ne.hmat) then
        write(6,'(/,''WARNING; symmetry determination not '',
     + ''properly specified'')')
      endif
c.....lmm end of check for input consistency 
      write(6,'(//''***** title : '',a50,''*****'')')title
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
      write(6,'(/'' Convergence'',/,''  maxiters = '',i4,/,
     + ''  threshold = '',g8.2)')maxiters,tstthr
      if (lrfun) then
        write(6,'(/'' Linear response procedure''/)')
        write(6,'(''  damping = '',f5.2)')dvdamp
      else
        write(6,'(''  shift   = '',f5.2)')df
        write(6,'(''  bounds  = '',f5.2,'','',f5.2)')crrmin,crrmax
      endif
      if (scfdamp.lt.1.d0) write(6,'(''  density damping : '',f5.2)')
     + scfdamp
      if (kpens.ne.0) then
        write(6,'(/'' Calculate ensemble from itr = '',i4,
     + ''   once every  '',i4,'' iterations'')')ibcens,iscens
        write(6,'(''   kpens = '',i4)')kpens
        write(6,'(''   dqmax = '',f4.2)')dqmax
      endif
clmm      if (isks.gt.0) then
clmm        write(6,'(''  Storage of KS orbitals  :'',i3)')isks
clmm      endif
c
      ios = 0
      icounter = 0
      open(99,file=trim(gridfile),status='old',err=444)
      rewind(99)
      do while (ios == 0)
        read(99,'(4e25.14)',iostat=ios) x,y,z,w
        icounter = icounter + 1
      enddo
      npnt = icounter - 1  
      write(6,'(/,''Grid points read from: '',a)') gridfile
      write(6,'(/,'' Number of gridpoints in numerical'',
     + '' integration :'',i9)') npnt
c
      if (lfield) write(6,'(/,'' Electric field : '',f6.4,
     + ''x + '',f6.4,''y + '',f6.4,''z'')')(fxyz(i),i=1,3)

c      if (ispi.gt.0) then
c      nmosa=2
c      nmosb=0
c     call brains(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmn,crrmx,
c     + thresh,dqmax,dvdamp,scfdamp,kpens,info,nppr,nvpr,npnt,npold,
c     + norb,nmos,nmomx,maxiters,ibcens,iscens,iint,idmp,ismo,isno,isao,
c     + isoe,isks,lsym,lintsm,lrdocc,lrfun,nmosa,nmosb,atmol4)
c      else
c      nmos=2
clmm..initialized kpens variable to avoid a memory leak later  
      kpens = 1
      call dbrain(occmo,fxyz,tstthr,alpha,beta,gamma,df,crrmin,crrmax,
     + thresh,dqmax,dvdamp,scfdamp,kpens,info,nppr,nvpr,npnt,
     + norb,nmos,nmomx,maxiters,ibcens,iscens,
     + lsym,lintsm,lrdocc,lrfun,nele,outksoFile) 
c      endif

clmm      call revind
clmm      call secsum
clmm      call whtps
clmm      call givtim(top,bot)
      write(6,10101)top,bot
10101 format(//' end of dnstyfun at',f13.3,' wall',f13.3,' secs')
      stop
c
  444 write(6,'(''ERROR; file '',a,'' does not exist'')') trim(gridfile)
      stop
c
      end
