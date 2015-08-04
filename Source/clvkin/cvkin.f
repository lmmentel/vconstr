cam... program that computes kinetic potential from a GAMESS-US
cam... or KS probability amplitude \Phi(2...N|N)
cam... 
cam... v_kin(1) = 1/2 \int |\nabla_1 \Phi(2...N|1)|^2 d2...dN
cam...
cam... v_s,kin(1) = 1/2 \int |\nabla_1 \Phi_s(2...N|1)|^2 d2...dN
cam...
cam... also the kinetic correlation potential is computed
cam... 
cam... v_c,kin(1) = v_kin(1) - v_c,kin(1)
cam...
cam... the potentials are computed for the points provided in the
cam... gridfile *.xyzw
cam...
cam... further a *.basinfo and *.F10 file from the GAMESS-US run 
cam... and a *.kso file from a dsfun run need to be present
      program cvkin
      use basisModule
      use commonsModule
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(nmomx  = 60)
c
      character*4 ied,iednum 
      character*44 fname
      common/discc/ied(16),fname(16)
c
cam      logical atmol4
c
cam      common /titel/title(10)
cam..input processing variables
      character*100 buffer, inptf, gbasisfile, gdictnfile,
     + gridfile, inpksoFile
      character*50  title
      character*80  gamtitle
      namelist /input/ 
     & gbasisfile, 
     & gdictnfile,
     & gridfile,
     & inpksoFile,
     & iprint, 
     & nmos,
     & nppr, 
     & title
cam
      gridfile='points'
      title = 'default'
      nmos = -1
      nppr = 0
      nvpr = 0
cam
      iprint = 0
cam
      rnel   = 0.d0
cam
cam      info(1:infmx)=-1
c
      write(6,'(''  Calculation of kinetic potential '')')
c
cam......input processing
c$$$c
c$$$      idmp = 2
c$$$      ismo = 1
c$$$      isno = 2
c$$$      isks = 3
c$$$      isao = -1
c$$$c*** dmpfile
c$$$    2 call inpi(idmp)
c$$$      call inpi(ismo)
c$$$      if ((ismo.lt.1).and.(ismo.gt.190)) then
c$$$        write(6,'(''ERROR; mo section not properly specified'')')
c$$$        stop
c$$$      endif
c$$$      call inpi(isno)
c$$$      if ((isno.lt.1).and.(isno.gt.190)) then
c$$$        write(6,'(''ERROR; no section not properly specified'')')
c$$$        stop
c$$$      endif
c$$$      call inpi(isks)
c$$$      if ((isks.lt.1).and.(isks.gt.190)) then
c$$$        write(6,'(''ERROR; ks section not properly specified'')')
c$$$        stop
c$$$      endif
c$$$      goto 50
c$$$c
c$$$c** adapt
c$$$    3 call inpi(isao)
c$$$      if ((isao.lt.1).and.(isao.gt.190)) then
c$$$        write(6,'(''ERROR; isao not properly specified'')')
c$$$        stop
c$$$      endif
c$$$      goto 50
c$$$c
c$$$c*** enter
c$$$    4 call getgss(idmp,norb,atmol4)
c$$$      iednum=ied(idmp)
cam..get input file name form command line open it, read the 
cam..contents and close
      call getarg(1, buffer)
      read(buffer, *) inptf
      open(unit=11, file=trim(inptf), status='old',form='formatted',
     &     delim='apostrophe', iostat=ios)
      read(11, nml=input)
      close(11)
cam..rewrite some information to commons
      basisInfoFile  = gbasisfile
      dictionaryFile = gdictnfile
      printLevel     = iprint
      write(6,'(/,''grid file from nml: '',a)') gridfile
c.....am end of input reading 
cam..
cam..get the gaussian basis information from the basis file
cam..
      call read_system_info(trim(gbasisfile), gamtitle, natoms, icharge, 
     & mult, nbf, nx, nele, na, nb, nshell, nprimi) 
cam..not sure if norb should be set to the number of cartesian gaussians
cam..or number of orbitals used in the calculation
      norb = nbf 
c.....nmos
      if (nmos.gt.-1) then
        if (nmos.gt.nmomx) then
          write(6,'(''ERROR; Too many Kohn-Sham orbitals specified; '',
     + i4,'' > '',i4)')nmos,nmomx
          stop
        endif
      endif
c$$$c.....symdet
c$$$      if (smtype.eq.int2e) then
c$$$        lintsm = .true.
c$$$      elseif (smtype.ne.hmat) then
c$$$        write(6,'(/,''WARNING; symmetry determination not '',
c$$$     + ''properly specified'')')
c$$$      endif
c.....amm end of check for input consistency 
      write(6,'(//''***** title : '', a50,''*****'')')title
      write(6,'(/''  norb = '',i3)')norb
c$$$      write(6,'(/'' Summary of dumpfile  '',a4,'' sections''/,
c$$$     + ''  Molecular Orbital basis :'',i3,/,
c$$$     + ''  Natural Orbitals        :'',i3,/,
c$$$     + ''  Kohn-Sham orbitals      :'',i3)')iednum,ismo,isno,isks
c$$$      if (isao.gt.0) then
c$$$        write(6,'(''  Symmetry adaptation     :'',i3)')isao
c$$$      endif
c
cam      open(99,file='points',status='old',err=222)
cam      rewind(99)
cam      read(99,*)npnt,npold
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
      write(6,'(/,'' computing vkin at '', i7, '' grid points'')')  npnt
cam      write(6,'(/,'' Number of gridpoints in numerical'',
cam     + '' integration :'',i5)')npold-1
cam      write(6,'('' Number of dummy points :'',i5)')npnt+1-npold
c
cam      call vbrain(npnt,npold,norb,idmp,ismo,isno,isao,isks,atmol4) 
      call vbrain(npnt,npold,norb,nmomx,nele,inpksoFile) 
c
cam      call revind
cam      call secsum
cam      call whtps
      stop
c
cam  222 write(6,'(''ERROR; file "points" does not exist'')')
cam      stop
c
  444 write(6,'(''ERROR; grid file '',a,'' does not exist'')') 
     + trim(gridfile)
      stop
c
      end
c
c
c***********************************************************************
c
      block data prg
      implicit real *8   (a-h,p-w), integer (i-n), logical (o)
      common /prgnam/prgtit
      data prgtit/'cvkin'/
      end
