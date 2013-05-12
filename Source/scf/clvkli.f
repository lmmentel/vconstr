      subroutine clvkli(vxc,vresp,wxi,wsi,orbdns,weight,occmo,ityp,
     + ndvx,nmos,npnt,npntmx)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(eps = 1.d-12)
c
      dimension vxc(npnt),vresp(npnt)
      dimension ityp(nmos)
      dimension wxi(nmos),wsi(nmos),occmo(nmos)
      dimension orbdns(npntmx,nmos)
      dimension weight(npnt)
      dimension wkli(ndvx),wovs(ndvx),occtyp(ndvx)
c
      dimension rho(npnt)
c
      wovs(1:ndvx)=0.d0
      occtyp(1:ndvx)=0.d0
      exhf=0.d0
      exks=0.d0
      do i=1,nmos
        exhf=exhf+wxi(i)
        exks=exks+wsi(i)
        l=ityp(i)
        if (l.ne.0) then
          wovs(l)=wovs(l)+occmo(i)*(wsi(i)-wxi(i))
          occtyp(l)=occtyp(l)+occmo(i)
        endif
      enddo
      print*,'exhf',exhf,'exks',exks,'diff',exhf-exks
c
      rho(1:npnt)=0.d0
      do i=1,ndvx
        rho(1:npnt)=rho(1:npnt)+orbdns(1:npnt,i)
      enddo
c
      call clcdvx(wkli,wovs,rho,orbdns,weight,occtyp,npnt,npntmx,ndvx)
      write(6,'(/''  wvs   : '',10(1x,f11.6))')
     + (wovs(i)/occtyp(i),i=1,ndvx)
      write(6,'(''  wkli  : '',10(1x,f11.6))')
     + (wkli(i),i=1,ndvx)
      write(6,'(''  wxi  : '',10(1x,f11.6))')
     + (wxi(i),i=1,ndvx)
c
c** new potential
c
      do ip=1,npnt
        vkli=0.d0
        if (rho(ip).gt.eps) then
          do i=1,ndvx
            vkli=vkli+orbdns(ip,i)*wkli(i)/rho(ip)
          enddo
          vkli=vkli-wkli(ndvx)
        endif
        vresp(ip)=vkli
        vxc(ip)=vxc(ip)+vkli
      enddo
c
      return
      end
