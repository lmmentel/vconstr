      subroutine rdmat(iunit,norb,tsmat,vnmat,enuc)
c***********************************************************************
c
c read matrix of kinetic and electron nuclear attraction integrals from
c dumpfile, section 192
c
c  see common/atbuf/ for layout of each block in section 192.
c  v(ij) = fm(ij) - tm(ij).   the value of i*(i+1)/2+j
c  is stored in array im.
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(nints = 77)
c
      common/atbuf/potnuc,dx,dy,dz,sm(nints),tm(nints),fm(nints),
     * xm(nints),ym(nints),zm(nints),im(nints),idum,rdum(6)
c
      dimension vnmat(norb*(norb+1)/2),tsmat(norb*(norb+1)/2)
c
c** find integrals
c
      call secini(1,iunit)
      call secget(192,2,iblock)
      call search(iblock,iunit)
      call reads(potnuc,1,iunit)
      enuc=potnuc
      ism=0
      do iorb=1,norb*(norb+1)/2
        ism=ism+1
        if (ism.gt.nints) then
          call reads(potnuc,1,iunit)
          ism=1
        endif
c
c** fm and tm are respectively the the one-electron fock matrix and the 
c** kinetic energy in ao-basis
c
        vnmat(im(ism))=fm(ism)-tm(ism)
        tsmat(im(ism))=tm(ism)
      enddo
c
      return
      end
