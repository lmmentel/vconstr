      subroutine clcdip(dxyz,pmo,iunit,norb)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(nints = 77)
c
      common/atbuf/potnuc,dx,dy,dz,sm(nints),tm(nints),fm(nints),
     * xm(nints),ym(nints),zm(nints),im(nints),idum,rdum(6)
c
      dimension dxyz(3),pmo(norb*(norb+1)/2)
c
      dxyz(1:3)=0.d0
c
c** find integrals
c
      call secini(1,iunit)
      call secget(192,2,iblock)
      call search(iblock,iunit)
      call reads(potnuc,1,iunit)
      ism=0
      do iorb=1,norb*(norb+1)/2
        ism=ism+1
        if (ism.gt.nints) then
          call reads(potnuc,1,iunit)
          ism=1
        endif
        dxyz(1)=dxyz(1)+pmo(im(ism))*xm(ism)
        dxyz(2)=dxyz(2)+pmo(im(ism))*ym(ism)
        dxyz(3)=dxyz(3)+pmo(im(ism))*zm(ism)
      enddo
c
      return
      end
