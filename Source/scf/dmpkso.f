      subroutine dmpkso(iunit,iblk,isec,vecks,eig,docc,nbas,norb)
c***********************************************************************
c
c Dump KS vectors in AO basis to dumpfile iunit, section isec, 
c   starting at block iblk.
c
c***********************************************************************
      implicit real*8  (a-h,o-z),integer(i-n)
c
      logical iftran
c
      dimension vecks(norb*norb),docc(norb),eig(norb)
c
      common/scratc/ilifc(256),ntran(256),itran(600),ctran(600),
     + iftran,idum
      common/scrtc/com(19),tit(10),value(255),occ(255),
     + nbasis,newbas,ncol,ivalue,iocc,idumm
      common/titel/title(10)
      data scf,akso/' scf',' kso'/
c
      nbasis = nbas
      newbas = nbas
      ncol   = norb
      ivalue = -1
      iocc   = 1
      nbsq  = norb*nbas
      n544  = 2*255+29+(5-1)/nipw()+1
      n1713 = (1113-1)/nipw()+1+600
      nblkk = (nbsq-1)/511+1+(n544-1)/511+1+(n1713-1)/511+1
      call fmove(title,tit,10)
      call fmove(eig,value,norb)
      call fmove(docc,occ,norb)
c
c*** determine date and time ***
c
      call tidajt(date,time,accno,ajnam,idum)
      com(1) = ajnam
      com(2) = date
      com(3) = time
      com(4) = scf
      com(5) = akso
      com(6) = accno
c
c*** write to dumpfile ***
c
      call secini(iblk,iunit)
      call secput(isec,3,nblkk,iblkk)
      call wrt3(com,n544,iblkk,iunit)
      call wrt3s(ilifc,n1713,iunit)
      call wrt3s(vecks,nbsq,iunit)
      write(6,987)norb,iunit,isec
  987 format(/,i6,' orbitals written to unit',i2,
     + ', section ',i4)
c
      call revind
c
      return
      end
