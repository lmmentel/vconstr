      subroutine spavxc(vxc,rho,nspin,npnt,npntmx)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension rho(npntmx,nspin),vxc(npntmx,nspin)
c
      do 100 i=1,npnt
        dnst=rho(i,1)+rho(i,2)
        tmp=vxc(i,1)*rho(i,1)+vxc(i,2)*rho(i,2)
        vxc(i,1)=tmp/dnst
  100 continue
c
      return
      end
