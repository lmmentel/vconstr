      program ensnergy
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical atmol4
      character*4 ied,iednum
      character*44 fname
      common/discc/ied(16),fname(16)
      common /titel/title(10)
c
      call prep99
      call tidajt(date,time,accno,anam,idum)
      write(6,666)anam,date,time,accno
  666 format(1h1//
     *24x,18('*'),' exc ',18('*') //
     *24x,'job name   ',a8//
     *24x,'date       ',a8//
     *24x,'time       ',a8//
     *24x,'acct       ',a8//
     *24x,44('*')///)
c
      idmp = 2
      ismo = 1
      isao = -1
      isks = 3
c
      call givtim(top,bot)
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
      go to (1,2,3,4)ii
c
c*** title
    1 read(5,599) title
599   format(10a8)
      write(6,'(/,''  title= '',10a8)')title
      goto 50
c
c*** dmpfile
    2 call inpi(idmp)
      call inpi(ismo)
      if ((ismo.lt.1).and.(ismo.gt.190)) then
        write(6,'(''ERROR; mo section not properly specified'')')
        stop
      endif
      call inpi(isks)
      if ((isks.lt.1).and.(isks.gt.190)) then
        write(6,'(''ERROR; no section not properly specified'')')
        stop
      endif
      goto 50
c
c** adapt
    3 call inpi(isao)
      write(6,'(/,''  AO orbitals are symmetry adapted '',/)')
      goto 50
c
c*** enter
    4 call getgss(idmp,norb,atmol4)
      write(6,'(//''***** title : '',10a8,''*****'')')title
      write(6,'(/''  norb = '',i3)')norb
      iednum=ied(idmp)
      write(6,'(/'' Summary of dumpfile  '',a4,'' sections''/,
     + ''  Molecular Orbital basis :'',i3,/,
     + ''  Kohn-Sham Orbitals      :'',i3)')iednum,ismo,isks
      if (isao.gt.0) then
        write(6,'(''  Symmetry adaptation     :'',i3)')isao
      endif
c
      call ebrain(norb,rnel,idmp,ismo,isks,isao,atmol4)
c
      call revind
      call secsum
      call whtps
      call givtim(top,bot)
      write(6,10101)top,bot
10101 format(//' end of exc at',f13.3,' wall',f13.3,' secs')
c
      stop
      end
c
      subroutine passdf(ip)
c     get data field of a card and see what it is
c     ip =  0  : not recognised
c     ip >  0  : it was item ip of the list
      implicit real*8(a-h,o-z)
      parameter (ncnt=4)
      character*8 cnt(ncnt),word
      data cnt/ 'title','dmpfile','adapt','enter' /
      call reada(word)
      ip=locatc(cnt,ncnt,word)
      return
      end
