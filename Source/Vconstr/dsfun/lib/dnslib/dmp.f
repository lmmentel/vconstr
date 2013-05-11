      subroutine dmpvec(iunit,isec,vksmo,eig,docc,norb)
c
c-----------------------------------------------------------------------
c
      implicit real*8  (a-h,o-z),integer(i-n)
c
      logical iftran
      dimension vksmo(norb*norb),docc(norb),eig(norb)
c
      character*4 ied,iednum
      character*44 fname
      common/discc/ied(16),fname(16)
c
      common/scratc/ilifc(256),ntran(256),itran(600),ctran(600),
     + iftran,idum
      common/scrtc/com(19),tit(10),value(255),occ(255),
     + nbasis,newbas,ncol,ivalue,iocc,idumm
      common/titel/title(10)
      common/prgnam/prgtit
      data akso/'  kso'/
c
      iftran = .true.
c
      nbasis = norb
      newbas = norb
      nbsq  = norb**2
      n544  = 2*255+29+(5-1)/nipw()+1
      n1713 = (1113-1)/nipw()+1+600
      nblkk = (nbsq-1)/511+1+(n544-1)/511+1+(n1713-1)/511+1
      call fmove(title,tit,10)
      call fmove(eig,value,norb)
      call fmove(docc,occ,norb)
c
c** determine date and time 
c
      call tidajt(date,time,accno,ajnam,idum)
      com(1) = ajnam
      com(2) = date
      com(3) = time
      com(4) = prgtit
      com(5) = akso
      com(6) = accno
c
c** write to dumpfile 
c
      call secini(1,iunit)
      call secput(isec,3,nblkk,iblkk)
      call wrt3(com,n544,iblkk,iunit)
      call wrt3s(ilifc,n1713,iunit)
      call wrt3s(vksmo,nbsq,iunit)
      iednum=ied(iunit)
      write(6,'(/,i6,''  vectors stored in section'',i4,
     *'' of dumpfile starting at block'',i6,'' of '',a4)')
     *norb,isec,iblkk,iednum 
      write(6,'(//'' header block information:''/
     *a7,''vectors created by '',a8,'' program at '',a8,
     *'' on '',a8,'' in the job '',a8,'' under acct. '',a8/
     *'' with the title '',10a8)')com(5),com(4),com(3),com(2),
     *com(1),com(6),tit
c
      call revind
c
      return
      end
c
      subroutine prvec(vector,norb,nmos,title)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      character*30 title
      dimension vector(norb*nmos)
c
      write(6,'(//''  ******* '',a30,'' ********'',/,
     +''    basis functions :'',i4,''  orbitals :'',i4)')title,norb,nmos
      n1 = nmos
      n2 = (n1-1) / 10  + 1
      do i=1,n2
        k=10
        if (i.eq.n2) k=n1-(i-1)*10
        ibas=(i-1)*10*norb
        write(6,'(//10i10)')((i-1)*10+l,l=1,k)
        write(6,'(/)')
        do j=1,norb
          write(6,'(i3,10(2x,f8.4))')j,(vector(ibas+j+(l-1)*norb),l=1,k)
        enddo
      enddo
c
      return
      end
