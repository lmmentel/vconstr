      subroutine clcwxi(wxi,occmo,vksmo,nspin,nmos,nmomx,norb,iunit)
      implicit real*8(a-h,o-z),integer(i-n)
      parameter(tol  = 1.d-12)
c
      dimension wxi(nmomx,nspin),occmo(nmomx,nspin),
     + vksmo(norb*norb,nspin)
      dimension ii(8),jj(8),ll(8),kk(8)
c
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
c
      wxi(1:nmos,1:nspin) = 0.d0
c
      call search(1,iunit)
  100 call find(iunit)
      call get(gin,nw)
      if (nw.eq.0) goto 200
      call unpack(gijkl,8,ijkl,4*nword)
      do m = 1,nword
        if (abs(gin(m)).gt.tol) then
          ii(1) = ijkl(1,m)
          jj(1) = ijkl(2,m)
          kk(1) = ijkl(3,m)
          ll(1) = ijkl(4,m)
          vijkl = gin(m)
          call nrperm(ii(1),jj(1),kk(1),ll(1),np)
          do is=1,nspin
            do imo = 1,nmos
              wxii = 0.d0
              iimo=(imo-1)*norb
              do jmo=1,nmos
                jjmo=(jmo-1)*norb
                do iperm = 1,np
                  wxii = wxii+
     + occmo(jmo,is)*vksmo(iimo+ii(iperm),is)*
     + vksmo(jjmo+jj(iperm),is)*vksmo(jjmo+kk(iperm),is)*
     + vksmo(iimo+ll(iperm),is)*vijkl
                enddo
              enddo
              wxi(imo,is) = wxi(imo,is)-wxii
            enddo
          enddo
        endif
      enddo
      goto 100
c
  200 wxi(1:nmos,1:nspin)=wxi(1:nmos,1:nspin)/2.d0
c
      return
      end
