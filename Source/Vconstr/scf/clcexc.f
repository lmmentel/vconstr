      subroutine clcexc(exld,exnl,ecld,ecnl,intpnt,npntmx,exc,weight)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(zero=0.0d0)
c
      dimension exc(npntmx,4),weight(npntmx)
c
      exld=zero
      ecld=zero
      exnl=zero
      ecnl=zero
      do ip=1,intpnt
        exld=exld+exc(ip,1)*weight(ip)
        ecld=ecld+exc(ip,2)*weight(ip)
        exnl=exnl+exc(ip,3)*weight(ip)
        ecnl=ecnl+exc(ip,4)*weight(ip)
      enddo
c
      return
      end
