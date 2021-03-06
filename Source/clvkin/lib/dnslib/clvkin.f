 1    subroutine clvkin(norb,npold,npnt,pnomo,pksmo,rho,dns,ddns,grdmo,
     + grid)
cam     + vck)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero=0.0d0)
c
      dimension pnomo(norb*(norb+1)/2),pksmo(norb*(norb+1)/2)
      dimension rho(npnt),dns(npnt),ddns(npnt,3),grdmo(norb*npnt,4)
      dimension grid(npnt,3)
cam      dimension vck(npnt)
c
cam      open(25,file='vkin.plt')
      open(26,file='vkin.dat')
      open(27,file='dns.dat')
cam      rewind(25)
      rewind(26)
      rewind(27)
cam      write(25,'(''  vckin        vkin        vskin'')')
      write(26,'(''#      x         y          z        
     + vckin        vkin        vskin'')')
      write(27,'(''#      x         y          z        
     + density'')')
      ekin = 0.0d0
      eskin = 0.0d0
      eckin = 0.0d0
      evw = 0.0d0
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
cam        vck(ip)=vckin
        vkin=vkin-vwke
cam        if (ip.ge.npold) write(25,'(3(f12.6,2x))') vckin,vkin,vskin
           write(26,'(6es18.10)') grid(ip,1), grid(ip,2), grid(ip,3), 
     + vckin, vkin,vskin
           write(27,'(6es18.10)') grid(ip,1), grid(ip,2), grid(ip,3), 
     + dns(ip)
cam... check for accuracy of grid by computing the kinetic energy
cam...
cam... T = \int dr \rho(r) ( v_x(r) + vwke(r))
cam...
cam... where x... kin or s,kin from the CI or KS wavefunction resp.
cam... to be compared with dsfun and GAMESS
cam... makes only sense if large grid is used, 
cam... usually larger than ploting grid
c$$$        ekin = ekin + vci * 0.5d0 * grid(ip,4)
c$$$        eskin = eskin + vks * grid(ip,4) * 0.5d0
c$$$        eckin = eckin + ( vci - vks) * grid(ip,4) * 0.5d0
c$$$        evw = evw + dns(ip) * vwke * grid(ip,4)
  100 continue
cam      close(25)
      close(26)
      close(27)
c$$$      write(6,'(/''CI kinetic energy            :'',f16.10)') ekin
c$$$      write(6,'(''KS kinetic energy            :'',f16.10)') eskin
c$$$      write(6,'(''Correlation kinetic energy   :'',f16.10)') eckin
c$$$      write(6,'(''Van Weizacker kinetic energy :'',f16.10)') evw
c
      return
      end
