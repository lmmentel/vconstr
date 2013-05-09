      program vcndex
c***********************************************************************
c Program to calculate values of conditional potential in a number
c of points. The points usually represent a numerical integration
c grid and are read from file "points". The calculation of vcond in all
c points (x,y,z) involves integrals of the type :  phi(i)*phi(j)/r12,
c integrate over 2.  These integrals are of the same type as the
c electron-nuclear attraction integrals and are evaluated exactly.
c
c the program is capable of handling atmol3 (integv) dumpfiles and
c atmol4 (atmol4) dumpfiles.  in the last case s,p,d,f and g functions
c can be used.
c !! g-functions not implemented at the moment.  see routine spherw !!
c
c parameters :
c   nhdimx : dimension of array's hx,hy and hz. these array's are
c            used for storage of point-independent constants.
c   lqntmx : maximum l-quantum number (l(s)=0, l(p)=1, ... etc )
c   ncmpmx : total number of possible gtf-type components (1 for s-type
c            gtf, 3 for p-type gtf (x,y,z), 6 for d-type gtf (xx,yy,zz,
c            xy,xz,yz), ... etc )
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      character*4 ied,iednum
      character*44 fname
      common/discc/ied(16),fname(16)
c
      common /titel/title(10)
c
      logical lvhart,lvcond,lvexch,atmol4
c
      call prep99
      call tidajt(date,time,accno,anam,idum)
      write(6,666)anam,date,time,accno
  666 format(1h1//
     *24x,18('*'),' vcond ',18('*') //
     *24x,'job name   ',a8//
     *24x,'date       ',a8//
     *24x,'time       ',a8//
     *24x,'acct       ',a8//
     *24x,44('*')///)
c
      ncore = 0
      nvpr = 100
c
      idmp = 2
      ismo = 1
      isao = -1
c
      ip2un = 6
c
      lvhart = .false.
      lvcond = .false.
      lvexch = .false.
c
      rnel = 0.d0
c
      call givtim(top,bot)
      write(6,444)top,bot
 444  format(/
     *' ****input read at',f13.3,' wall',f13.3,' secs')
      goto 55
c
   50 call input
      write(6,'(''I passed input'')')
   55 call passdf(ii)
      if (ii.le.0) then
         call cinput(jrec,jump,0)
         call erroutf(jrec)
         call caserr(' password not recognised ')
      endif
      go to (1,2,3,4,5,6,7,8,9,10,11)ii
c
c*** title
    1 read(5,599) title
599   format(10a8)
      write(6,'(/,''  title= '',10a8,/)')title
      goto 50
c
c*** ncore
    2 call inpi(ncore)
      goto 50
c
c*** iprint
    3 call inpi(nvpr)
      goto 50
c
c*** dmpfile
    4 call inpi(idmp)
      call inpi(ismo)
      goto 50
c
c*** p2file
    5 call inpi(ip2un)
      goto 50
c
c*** vhartr
    6 lvhart=.true.
      goto 50
c
c*** vcond
    7 lvcond=.true.
      goto 50
c
c** vexch
    8 lvexch=.true.
      goto 50
c
c** adapt
    9 call inpi(isao)
      goto 50
c
c** rnel
   10 call inpf(rnel)
      goto 50
c
c*** enter
   11 write(6,'(''I approach getgss'')')
      call getgss(idmp,norb,atmol4)
      write(6,'(''I passed getgss'')')
      if (norb.le.0) then
        write(6,'(''ERROR; norb not properly specified'')')
        stop
      endif
      if (ncore.lt.0) then
        write(6,'(''ERROR; ncore not properly specified'')')
        stop
      endif
      write(6,'(//''***** title : '',10a8,''*****'')')title
      write(6,'(/''  ncore = '',i3)')ncore
      write(6,'(''  norb  = '',i3)')norb
      if (rnel.le.0.d0) then
        write(6,'(''ERROR; rnel not properly specified'')')
        stop
      else
        write(6,'(''  # of electrons : '',f8.2)')rnel
      endif
      iednum=ied(idmp)
      write(6,'(/'' Summary of dumpfile  '',a4,'' sections''/,
     + ''  Molecular Orbital basis :'',i3)')iednum,ismo
      if (isao.gt.0) then
        write(6,'(''  Symmetry adapted basis :'',i3)')isao
      endif
      iednum=ied(ip2un)
      write(6,'(/'' two-electron density matrix on on file  '',
     + a4)')iednum
c
      if ((.not.lvhart).and.(.not.lvcond).and.(.not.lvexch)) then
        write(6,'(''ERROR; no potential specified'')')
        stop
      elseif (lvhart) then
         write(6,'('' calculate Hartree potential'')')
      elseif (lvcond) then
         write(6,'('' calculate Conditional potential'')')
      elseif (lvexch) then
         write(6,'('' calculate Exchange potential'')')
      endif
c
      if (atmol4) then
        call cnwhmx(nhmx,nprmx,ngrpx)
      else
        call cnvhmx(nhmx,nprmx,ngrpx)
      endif
      write(6,'(//,'' # of H-matrix elements : '',i6,/,
     + '' # of prefactors        : '',i6)')nhmx,nprmx
c
      open(10,file='points',status='old',err=222)
      rewind(10)
c
      call vbrain(nhmx,nprmx,ngrpx,norb,ncore,nvpr,ip2un,
     + idmp,ismo,isao,rnel,lvhart,lvcond,lvexch,atmol4)
c
      call givtim(top,bot)
      write(6,'(//'' end of vcond at'',f13.3,'' wall'',
     + f13.3,'' secs'')')top,bot
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
      parameter (ncnt=11)
      character*8 cnt(ncnt),word
      data cnt/ 'title','ncore','iprint','dmpfile','p2file',
     * 'vhartr','vcond','vexch','adapt','rnel','enter' /
      call reada(word)
      ip=locatc(cnt,ncnt,word)
c
      return
      end
