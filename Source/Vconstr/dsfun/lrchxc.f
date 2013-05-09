      subroutine lrchxc(vxc,rho,orbdnst,drho,dsrho,dns,ddns,dsdns,
     + weight,occmo,damp,dvmax,drhmx,idvmax,idrhmx,df,ityp,intpnt,
     + npnt,npntmx,nmos,ndvcr)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(eps    = 1.0d-25)
      parameter(fourth = 1.0d0/4.0d0)
      parameter(zero   = 0.0d0)
c
      dimension vxc(npnt)
      dimension rho(npnt),orbdnst(ndvcr*npnt),drho(npntmx,3),
     + dsrho(npnt),dns(npnt),ddns(npntmx,3),dsdns(npnt),weight(npnt)
      dimension occmo(nmos),ityp(nmos)
      dimension dvxc(ndvcr)
      dimension Avec(npnt)
c
      idrhmx = 0
      idvmax = 0
      dvmax  = 0.d0
      drhmx  = 0.d0
      Avec(1:npnt)  = 0.d0
      dvxc(1:ndvcr) = 0.d0
c
      do ipnt=1,npnt
        crr=(dns(ipnt)+df)/(rho(ipnt)+df)
        if (ipnt.le.intpnt) then
          if(crr.gt.drhmx) then
            drhmx=crr
            idrhmx=ipnt
          endif
        endif
        if (rho(ipnt).gt.eps) then
          factor=dns(ipnt)/rho(ipnt)
          Avec(ipnt)=dsdns(ipnt)-factor*dsrho(ipnt)-
     + (ddns(ipnt,1)*drho(ipnt,1)+ddns(ipnt,2)*drho(ipnt,2)+
     + ddns(ipnt,3)*drho(ipnt,3)-factor*(drho(ipnt,1)**2+
     + drho(ipnt,2)**2+drho(ipnt,3)**2))/rho(ipnt)
        endif
      enddo
c
      if (ndvcr.ne.1) then
        call dvksmt(dvxc,intpnt,nmos,ndvcr,rho,orbdnst,Avec,weight,
     + occmo,ityp)
        do i=1,ndvcr
          write(6,'(''  orbital'',i2,'' correction'',g10.2)')i,dvxc(i)
        enddo
      endif
c
c new potential
c
      dmxv=zero
      do ipnt=1,npnt
        if (rho(ipnt).gt.eps) then
          knmo=(ipnt-1)*ndvcr
          dvcr=fourth*Avec(ipnt)/rho(ipnt)
          do i=1,ndvcr
            dvcr=dvcr+orbdnst(knmo+i)*dvxc(i)/rho(ipnt)
          enddo
          dvcr=dvcr-dvxc(ndvcr)
          if (abs(dvcr).gt.dmxv) then
            dmxv=abs(dvcr)
          endif
          if (ipnt.le.intpnt) then
            err=abs(damp*dvcr/vxc(ipnt))
            if(err.gt.dvmax) then
              dvmax=err
              idvmax=ipnt
            endif
          endif
          vxc(ipnt)=vxc(ipnt)+damp*dvcr
        endif
      enddo
c
      return
      end
