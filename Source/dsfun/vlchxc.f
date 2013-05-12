      subroutine vlchxc(vxc,rho,dns,dvmax,drhmx,idvmax,idrhmx,df,
     + crrmn,crrmx,npnt,intpnt)
c***********************************************************************
c
c Calculate the corrections to the KS-potential
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension rho(npnt),dns(npnt),vxc(npnt)
c
      idrhmx = 0
      idvmax = 0
      dvmax  = 0.d0
      drhmx  = 0.d0
c
      do ipnt=1,npnt
        cr=(dns(ipnt)+df)/(rho(ipnt)+df)
        if (ipnt.le.intpnt) then
          crr=max(cr,1.0d0/cr)
          if(crr.gt.drhmx) then
            drhmx=crr
            idrhmx=ipnt
          endif
        endif
        if(cr.gt.1.d0) then
          cr=min(cr,crrmx)
        else if(cr.lt.1.d0) then
          cr=max(cr,crrmn)
        endif
        vxcn=vxc(ipnt)*cr
        if (ipnt.le.intpnt) then
          err=abs((vxcn-vxc(ipnt))/vxc(ipnt))
          if(err.gt.dvmax) then
            dvmax=err
            idvmax=ipnt
          endif
        endif
        vxc(ipnt)=vxcn
      enddo
c
      return
      end
