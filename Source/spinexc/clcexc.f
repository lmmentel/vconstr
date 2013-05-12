      subroutine clcexc(npnt,npntmx,exc,weight)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension exc(npntmx,7),weight(npnt)
c
      exld=0.d0
      ecld=0.d0
      exb=0.d0
      expw=0.d0
      ecp=0.d0
      ecpw=0.d0
      eclyp=0.d0
      do ip=1,npnt
        exld=exld+exc(ip,1)*weight(ip)
        ecld=ecld+exc(ip,2)*weight(ip)
        exb=exb+exc(ip,3)*weight(ip)
        expw=expw+exc(ip,4)*weight(ip)
        ecp=ecp+exc(ip,5)*weight(ip)
        ecpw=ecpw+exc(ip,6)*weight(ip)
        eclyp=eclyp+exc(ip,7)*weight(ip)
      enddo
c
      write(6,'(/''** Exchange energies **'')')
      write(6,'(6x,''LDA    :'',12x,f10.6)')exld
      write(6,'(6x,''Becke  :'',f10.6,2x,f10.6)')exb,exld+exb
      write(6,'(6x,''PW91   :'',f10.6,2x,f10.6)')expw,exld+expw
c
      write(6,'(/''** Correlation energies **'')')
      write(6,'(6x,''LDA    :'',12x,f10.6)')ecld
      write(6,'(6x,''LYP    :'',12x,f10.6)')eclyp
      write(6,'(6x,''Perdew :'',f10.6,2x,f10.6)')ecp,ecld+ecp
      write(6,'(6x,''PW91   :'',f10.6,2x,f10.6)')ecpw,ecld+ecpw
c
      return
      end
