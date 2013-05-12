      subroutine grdrho(orbdnst,rho,drho,dsrho,ityp,nmos,ndvcr,norb,
     + npnt,npntmx,pmo,vksmo,occmo,valmo,grdmo,grid,weight,rnel)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.0d-4)
c
      logical ltrian
      dimension orbdnst(ndvcr*npnt),rho(npnt),drho(npntmx,3),
     + dsrho(npnt)
      dimension pmo(norb*(norb+1)/2),vksmo(norb*norb),occmo(nmos),
     + ityp(nmos)
      dimension valmo(npnt*norb),grdmo(npntmx*norb,4),
     + grid(npntmx,3),weight(npnt)
      dimension chkval(nmos)
c
      ltrian = .false.
c
      chkval(1:nmos)=0.d0
      orbdnst(1:ndvcr*npnt)=0.d0
c
      tmpd=0.d0
      tmpr=0.d0
      tmps=0.d0
      tmpf=0.d0
      do 100 ipnt=1,npnt
        k=(ipnt-1)*norb
        knm=(ipnt-1)*ndvcr
        m=k+1
        rho(ipnt)=valmat(pmo(1),norb,valmo(m),ltrian)
        call valgrd(gx,gy,gz,pmo(1),norb,valmo(m),grdmo(m,1),
     + grdmo(m,2),grdmo(m,3))
        drho(ipnt,1)=gx
        drho(ipnt,2)=gy
        drho(ipnt,3)=gz
        dsrho(ipnt)=valsgrd(pmo(1),valmo(m),grdmo(m,1),grdmo(m,2),
     + grdmo(m,3),grdmo(m,4),norb)
        tmpd=tmpd+rho(ipnt)*weight(ipnt)
        tmpr=tmpr+(grid(ipnt,1)*drho(ipnt,1)+grid(ipnt,2)*drho(ipnt,2)+
     + grid(ipnt,3)*drho(ipnt,3))*weight(ipnt)
        tmps=tmps+dsrho(ipnt)*
     + (grid(ipnt,1)**2+grid(ipnt,2)**2+grid(ipnt,3)**2)*weight(ipnt)
        tmpf=tmpf+rho(ipnt)*(drho(ipnt,1)+drho(ipnt,2)+
     + drho(ipnt,3))*weight(ipnt)
        do i=1,nmos
          tmp=0.d0
          ii=(i-1)*norb
          l=ityp(i)
          if (l.ne.0) then
            do j=1,norb
              tmp=tmp+vksmo(ii+j)*valmo(k+j)
            enddo
            orbdnst(knm+l)=orbdnst(knm+l)+occmo(i)*tmp*tmp
            chkval(i)=chkval(i)+tmp*tmp*weight(ipnt)
          endif
        enddo
  100 continue
      write(6,'(/''  Integration accuracy : '',4(g8.2))')
     + tmpd-rnel,tmpr/3+rnel,tmps/6-rnel,tmpf
c
      do i=1,nmos
        if (ityp(i).ne.0) then
          tmp=chkval(i)-1.d0
          if (tmp.gt.eps) write(6,'('' WARNING: KS orbital '',
     + i4,'' yields an error upon integration of '',g8.2,
     + ''%'')')i,tmp*100
        endif
      enddo
c
      return
      end
