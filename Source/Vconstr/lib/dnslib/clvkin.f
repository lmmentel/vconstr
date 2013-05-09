      subroutine clvkin(norb,npold,npnt,pnomo,pksmo,rho,dns,ddns,grdmo
     +,vck)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero=0.0d0)
c
      dimension pnomo(norb*(norb+1)/2),pksmo(norb*(norb+1)/2)
      dimension rho(npnt),dns(npnt),ddns(npnt,3),grdmo(norb*npnt,3)
      dimension vck(npnt)
c
      open(25,file='vkin.plt')
      open(26,file='vkin.dat')
      rewind(25)
      rewind(26)
      write(25,'(''  vckin        vkin        vskin'')')
      do 100 ip=1,npnt
        vci=zero
        vks=zero
        k=(ip-1)*norb
        do i=1,norb
          ii=i*(i-1)/2
          tmp=grdmo(k+i,1)**2+grdmo(k+i,2)**2+grdmo(k+i,3)**2
          vci=vci+pnomo(ii+i)*tmp
          vks=vks+pksmo(ii+i)*tmp
          do j=1,i-1
            tmp=grdmo(k+i,1)*grdmo(k+j,1)+
     + grdmo(k+i,2)*grdmo(k+j,2)+grdmo(k+i,3)*grdmo(k+j,3)
            vci=vci+pnomo(ii+j)*tmp
            vks=vks+pksmo(ii+j)*tmp
          enddo
        enddo
        twrh=2.d0*rho(ip)
        twdn=2.d0*dns(ip)
        dsdn=ddns(ip,1)**2+ddns(ip,2)**2+ddns(ip,3)**2
        vwke=dsdn/(8.d0*dns(ip)**2)
        vskin=vks/twrh
        vkin=vci/twdn
        vckin=vkin-vskin
        vskin=vskin-vwke
        vck(ip)=vckin
        vkin=vkin-vwke
        if (ip.ge.npold) write(25,'(3(f12.6,2x))') vckin,vkin,vskin
        write(26,*) vckin,vkin,vskin
  100 continue
      close(25)
      close(26)
c
      return
      end
