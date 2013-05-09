      program scf
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
c
      character*4 yname,ydnam
      character*44 fname
      common/discc/yname(16),fname(16)
c
      character*8 smtype,int2e,scftyp,restr,unrestr,avrgvxc
      logical ldisk,lsym,lintsm,lrdocc,lb88,lpw91x,lprdwc,lpw91c,llyp,
     + atmol4
c
      dimension occmo(nmomx,2)
c
      common /titel/title(10)
c
      data int2e/'2-elec'/
      data restr,unrestr,avrgvxc/'rest','unrest','avrgxc'/
c
      call prep99
      call tidajt(date,time,accno,anam,idum)
      write(6,666)anam,date,time,accno
 666  format(1h1//
     *24x,18('*'),' scf ',18('*') //
     *24x,'job name   ',a8//
     *24x,'date       ',a8//
     *24x,'time       ',a8//
     *24x,'acct       ',a8//
     *24x,44('*')///)
      call givtim(top,bot)
      write(6,444)top,bot
 444  format(/
     *' ****input read at',f13.3,' wall',f13.3,' secs')
c
      nmos = -1
      nvpr = -1
      idamp = -1
      nsp   = 1
      nspin = 1
      itrx  = 50
      ispin = 0
c
      ikli = 0
c
      idmp = 2
      ismo = 1
      isoe = 198
      isao = -1
      isks = -1
c
      iint = 7
c
      kpens = 0
      nsens = 0
c
      lsym   = .false.
      lintsm = .false.
      lrdocc = .true.
c
      ldisk  = .true.
      lb88   = .false.
      lpw91x = .false.
      lpw91c = .false.
      llyp   = .false.
c
      damp1 = 0.3d0
      damp2 = 0.3d0
      tstthr = 1.0d-6
      thresh = 1.0d-6
c
      occmo(1:nmomx,1:2) = 0.d0
c
      write(6,'(''  Selfconsistent Kohn-Sham calculation '')')
      goto 55
c
   50 call input
   55 call passdf(ii)
      if (ii.le.0) then
         call cinput(jrec,jump,0)
         call erroutf(jrec)
         call caserr(' password not recognised ')
      endif
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)ii
c
c*** title
    1 read(5,'(10a8)') title
      goto 50
c
c*** spin
    2 call inpa(scftyp)
      if (scftyp.eq.unrestr) then
        nsp  =2
        nspin=2
      elseif(scftyp.eq.avrgvxc) then
        nsp  =1
        nspin=2
      elseif (scftyp.ne.restr) then
        write(6,'(/,''WARNING; type of calculation not '',
     + ''properly specified.'')')
      endif
      goto 50
c
c*** ensdet
    3 call inpi(kpens)
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
c*** nmos
    4 call inpi(nmos)
      if (nmos.gt.nmomx) then
        write(6,'(''ERROR; Too many Kohn-Sham orbitals specified; '',
     + i4,'' > '',i4)')nmos,nmomx
        stop
      endif
      goto 50
c
c*** occ
    5 if (nmos.le.0) then
        write(6,'(''ERROR; nmos not yet defined'')')
        stop
      endif
      if (ispin.le.nspin) then
        ispin=ispin+1
        do i=1,nmos
          call inpf(occmo(i,ispin))
        enddo
      else
        write(6,'(''ERROR; Too many occupation patterns '',
     + ''specified'')')
        stop
      endif
      lrdocc=.false.
      goto 50
c
c*** conv
    6 call inpi(itrx)
      call inpf(tstthr)
      goto 50
c
c*** damping
    7 call inpf(damp1)
      call inpf(damp2)
      call inpi(idamp)
      if (idamp.le.0) damp2 = damp1
      goto 50
c
c*** symdet
    8 lsym = .true.
      call inpa(smtype)
      if(smtype.eq.int2e) then
        lintsm = .true.
      endif
      call inpf(thresh)
      goto 50
c
c*** becke
    9 lb88=.true.
      goto 50
c
c*** pw91x
   10 lpw91x=.true.
      goto 50
c
c*** perdew
   11 lprdwc=.true.
      goto 50
c
c*** pw91c 
   12 lpw91c=.true.
      goto 50 
c
c*** lyp
   13 llyp=.true.
      goto 50
c
c*** vprint 
   14 call inpi(nvpr)
      goto 50
c
c*** dmpfile
   15 call inpi(idmp)
      call inpi(ismo)
      call inpi(isoe)
      goto 50
c
c*** intfile
   16 call inpi(iint)
      goto 50
c
c*** ksdmp
   17 call inpi(isks)
      goto 50
c
c** adapt
   18 call inpi(isao)
      goto 50
c
c** kli
   19 ikli=1
      goto 50
c
c** direct
   20 ldisk=.false.
      goto 50
c
c*** enter
   21 call getgss(idmp,norb,atmol4)
      write(6,'(//''***** title : '',10a8,''*****'')')title
      write(6,'(/''  norb = '',i3)')norb
      if (nmos.le.0) then
        nmos=nmomx
      else
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
        do is=1,nspin
          write(6,'(''  occupation : '',3(15f5.2,/,''     ''))')
     +(occmo(i,is),i=1,nmos)
        enddo
      endif
      if (nvpr.le.0) nvpr = nmos
c
      if (scftyp.eq.unrestr) then
        write(6,'(/,''  Spin unrestricted calculation'',/)')
      elseif(scftyp.eq.restr) then
        write(6,'(/,''  Spin restricted calculation'',/)')
      elseif(scftyp.eq.avrgvxc) then
        write(6,'(/,'' Exchange-correlation potential is '',
     + ''spin-averaged (sort of restricted calculation)'',/)')
      else
        write(6,'(/,''  Spin restricted calculation is '',
     + ''assumed'',/)')
      endif
      if (ikli.ne.0) then
        write(6,'('' KLI'')')
      elseif (lb88.and.lpw91x) then 
        write(6,'('' ERROR; Two exchange potentials specified'')')
        stop
      elseif ((lprdwc.and.lpw91c.and.llyp).or.(lprdwc.and.lpw91c).or.
     + (lprdwc.and.llyp).or.(lpw91c.and.llyp)) then
        write(6,'('' ERROR; More than one correlation potential'',
     + '' specified'')')
        stop
      elseif (lb88) then
        if (lprdwc) then
          write(6,'(/,'' LDA + Becke exchange + Perdew '',
     + ''correlation'')')
        elseif (lpw91c) then
          write(6,'(/,'' LDA + Becke exchange + Perdew-Wang '',
     + ''correlation'')')
        elseif (llyp) then
          write(6,'(/,'' LDA + Becke exchange + Lee-Yang-Parr '',
     + ''correlation'')')
        else
          write(6,'(/,'' LDA + Becke exchange'')')
        endif
      elseif (lpw91x) then
        if (lprdwc) then
          write(6,'(/,'' LDA + Perdew-Wang exchange and '',
     + ''Perdew correlation'')')
        elseif (lpw91c) then
          write(6,'(/,'' LDA + Perdew-Wang exchange and '',
     + ''correlation'')')
        elseif (llyp) then
          write(6,'(/,'' LDA + Perdew-Wang exchange + '',
     + ''Lee-Yang-Parr correlation'')')
        else
          write(6,'(/,'' LDA + Perdew-Wang exchange'')')
        endif
      else
        write(6,'(/,'' LDA '')')
      endif
      write(6,'(/'' Convergence reached when the absolute integrated'',
     + '' density difference is less than '',g8.2)') tstthr
      write(6,'(''  maximum # iterations = '',i4)')itrx
      if (idamp.lt.0) then
        write(6,'(/''  rho(n+1) = '',f5.2,'' rho(n) + '',f5.2,
     + '' rho(n-1)'')')damp1,1.d0-damp1
      else
        write(6,'(/''  rho(n+1) = '',f5.2,'' rho(n) + '',f5.2,
     + '' rho(n-1) for itr < '',i3)')damp1,1.d0-damp1,idamp
        write(6,'(''  rho(n+1) = '',f5.2,'' rho(n) + '',f5.2,
     + '' rho(n-1) for itr > '',i3)')damp2,1.d0-damp2,idamp
      endif
      if (kpens.ne.0) then
        write(6,'(/''  calculate possible ensemble'')')
        write(6,'(''    kpens = '',i4)')kpens
        write(6,'(''    dqmax = '',f4.2)')dqmax
        nsens = nspin
      endif
      if (lsym) then
        if (lintsm) then
          write(6,'(/'' Symmetry determined from orbitals'')')
        else
          write(6,'(/'' Symmetry determined from one-electron '',
     + ''h-matrix'')')
        endif
      endif
      ydnam=yname(idmp)
      write(6,'(/'' Summary of dumpfile  '',a4,'' sections''/,
     + ''  Molecular Orbital basis :'',i3,/,
     + ''  one-electron h-matrix   :'',i3)')ydnam,ismo,isoe
      if (isao.gt.0) then
        write(6,'(''  Symmetry adaptation     :'',i3)')isao
      endif
      if (isks.gt.0) then
        write(6,'(''  Storage of KS orbitals  :'',i3)')isks
      endif
      ydnam=yname(iint)
      write(6,'(/'' Integrals on file  '',a4)')ydnam
      if (.not.ldisk) then
        write(6,'(/'' MO-basis calculated each iteration'')')
      else
        write(6,'(/'' MO-basis stored in file hfvls'')')
      endif
c
      open(99,file='points',status='old',err=222)
      rewind(99)
      read(99,*)npnt,npold
      write(6,'(/,'' number of gridpoints in numerical'',
     + '' integration :'',i5)')npold-1
      write(6,'('' number of dummy points :'',i5)')npnt+1-npold
c
      call sbrain(occmo,nspin,nsp,nmos,nmomx,norb,npnt,npold,
     + itrx,tstthr,damp1,damp2,dqmax,kpens,nsens,idamp,nvpr,thresh,
     + iint,idmp,ismo,isao,isks,isoe,ikli,lb88,lpw91x,lprdwc,
     + lpw91c,llyp,lsym,lintsm,lrdocc,ldisk,atmol4)
c
      call revind
      call secsum
      call whtps
      call givtim(top,bot)
      write(6,'(//'' end of scf at'',f13.3,'' wall'',f13.3,
     + '' secs'')')top,bot
      stop
c
  222 write(6,'(''ERROR; file "points" does not exist'')')
      stop
c
      end
c
      subroutine passdf(ip)
c     get data field of a card and see what it is
c     ip =  0  : not recognised
c     ip >  0  : it was item ip of the list
      implicit real*8(a-h,o-z)
      parameter (ncnt=21)
      character*8 cnt(ncnt),word
      data cnt/ 'title','spin','ensdet','nmos','occ','conv','damping',
     * 'symdet','becke','pw91x','perdew','pw91c','lyp','vprint',
     * 'dmpfile','intfile','ksdmp','adapt','kli','direct','enter' /
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
      data prgtit/'scf'/
      end
