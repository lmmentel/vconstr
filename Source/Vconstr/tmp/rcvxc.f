      program rcvxc
      implicit real*8 (a-h,o-z),integer(i-n)
c
      character*4 ied,iednum 
      character*44 fname
      common/discc/ied(16),fname(16)
c
      logical lfield,atmol4
      dimension fxyz(3)
c
      common /titel/title(10)
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
      idmp = 2
      ismo = 1
      isno = 2
      isao = -1
c
      lfield = .false.
c
      fxyz(1:3)=0.d0
c
      write(6,'(''  Reconstruction of Kohn-Sham potential'')')
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
      go to (1,2,3,4,5)ii
c
c*** title
    1 read(5,'(10a8)') title
      goto 50
c
c*** dmpfile
    2 call inpi(idmp)
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
      goto 50
c
c** adapt
    3 call inpi(isao)
      if ((isao.lt.1).and.(isao.gt.190)) then
        write(6,'(''ERROR; isao not properly specified'')')
        stop
      endif
      goto 50
c
c*** dipole
    4 lfield=.true.
      do i=1,3
        call inpf(fxyz(i))
      enddo
      goto 50
c
c*** enter
    5 call getgss(idmp,norb,atmol4)
      write(6,'(//''***** title : '',10a8,''*****'')')title
      write(6,'(/''  norb = '',i3)')norb
c
      iednum=ied(idmp)
      write(6,'(/'' Summary of dumpfile  '',a4,'' sections''/,
     + ''  Molecular Orbital basis :'',i3,/,
     + ''  Natural Orbitals        :'',i3)')iednum,ismo,isno
      if (isao.gt.0) then
        write(6,'(''  Symmetry adaptation     :'',i3)')isao
      endif
c
      open(99,file='points',status='old',err=444)
      rewind(99)
      read(99,*)npnt,npold
      write(6,'('' Number of dummy points :'',i5)')npnt+1-npold
c
      if (lfield) write(6,'(/,'' Electric field : '',f6.4,
     + ''x + '',f6.4,''y + '',f6.4,''z'')')(fxyz(i),i=1,3)
c
      call cbrain(fxyz,npnt,norb,idmp,ismo,isno,isao,atmol4) 
c
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
      parameter (ncnt=5)
      character*8 cnt(ncnt),word
      data cnt/ 'title','dmpfile','adapt','dipole','enter' /
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
