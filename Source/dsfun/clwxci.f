      subroutine clwxci(occmo,vksmo,nmos,norb,iunit)
      implicit real*8(a-h,o-z),integer(i-n)
      parameter(tol  = 1.d-12)
c
      dimension wxci(nmos),occmo(nmos),vksmo(nmos*norb)
      dimension ii(8),jj(8),ll(8),kk(8)
c
      common/atbuf/gin(340),gijkl(170),nword,ndum
      common/scijkl/ijkl(4,340)
c
      wxci(1:nmos) = 0.d0
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
          do imo = 1,nmos
            wxcii = 0.d0
            iimo=(imo-1)*norb
            do jmo=1,nmos
              jjmo=(jmo-1)*norb
              occj=occmo(jmo)
              do iperm = 1,np
                wxcii = wxcii+
     + occj*vksmo(iimo+ii(iperm))*vksmo(jjmo+jj(iperm))*
     + vksmo(jjmo+kk(iperm))*vksmo(iimo+ll(iperm))*vijkl
              enddo
            enddo
            wxci(imo) = wxci(imo)-wxcii
          enddo
        endif
      enddo
      goto 100
c
  200 wxci(1:nmos)=wxci(1:nmos)/4.d0
      write(6,'(//'' The weights of the orbital dependent Hartree '',
     + ''Fock exchange potentials are :'')')
      write(6,'(10(2x,g14.6))')(wxci(imo),imo=1,nmos)
c
      return
      end
