      subroutine getdmp(iunit,isec,vec,occ,norb)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical iftran
      character*4 ied,iednum
      character*44 fname
      dimension vec(norb*norb),occ(norb)
      dimension tmpvec(norb)
c
      common/discc/ied(16),fname(16)
      common/tran/ilifc(256),ntran(256),itran(600),ctran(600),
     +iftran,ndum
      common/scrtc/com(19),tit(10),value(255),docc(255),
     * nbasis,newbas,idum(4)
c
      call secini(1,iunit)
      call secget(isec,3,iblock)
      iednum=ied(iunit)
      n544 = 29+2*255+(5-1)/nipw()+1
      call rdedx(com,n544,iblock,iunit)
      if ((norb.ne.nbasis).or.(nbasis.ne.newbas)) then
        write(6,'(/,'' ERROR; norb = '',i4,/
     +              ''        nbas = '',i4,/
     +              ''        nwbas= '',i4)')
     + norb,nbasis,newbas
        stop
      endif
      occ(1:norb)=docc(1:norb)
c
      n1713 = (1113-1)/nipw()+600+1
      call reads(ilifc,n1713,iunit)
c
c** total number of vector elements : norb*norb
c
      ndim = norb**2
      call reads(vec,ndim,iunit)
      write(6,'(/''  '',i4,
     *'' vectors and occupation numbers restored from section'',
     *i4,'' of dumpfile starting at block'',i6,'' of '',a4//
     *'' header block information:''/
     *a7,''vectors created by '',a8,'' atmol program at '',a8,
     *'' on '',a8,'' in the job '',a8,'' under acct. '',a8/
     *'' with the title '',10a8)')norb,isec,iblock,iednum,com(5),
     *com(4),com(3),com(2),com(1),com(6),tit
c
      if (iftran) return
c
**  ctrans option was used.  Now do the same transformation 
**  backwards. If newbas.lt.nbas then first adapt the vector-lenght.

      write(6,'(/////14x,''list of lcbf''//
     *'' coefficient old orbital nterm new orbital'')')
      do i=1,newbas
        n=ntran(i)
        j=ilifc(i)
        write(6,'(f12.7,i12,i6,i12)')ctran(j+1),itran(j+1),n,i
        do k=2,n
          write(6,'(f12.7,i12)')ctran(j+k),itran(j+k)
        enddo
      enddo
      write(6,'(/'' no. of lcbf='',i4)')newbas
c
      do i=1,norb
      tmpvec(1:norb)=0.d0
        ibas=(i-1)*nbasis
        do j=1,newbas
          veccf = vec(ibas+j)
           do k=1,ntran(j)
             iao = itran(ilifc(j)+k)
             cao = ctran(ilifc(j)+k)
             tmpvec(iao) = tmpvec(iao)+cao*veccf
       if(i.eq.3) then
       if((iao.eq.1).or.(iao.eq.24)) then
       write(6,'(2i2,3f20.10)')iao,k,cao,veccf,tmpvec(iao)
       endif
       endif
           enddo
         enddo
c
        do j = 1,nbasis
          vec(ibas+j) = tmpvec(j)
        enddo
      enddo
c
      return
      end
