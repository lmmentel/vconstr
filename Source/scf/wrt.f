      subroutine wrtpot(rho,vxc,exc,nspin,intpnt,npntmx)
c***********************************************************************
c
c Write data to the files vxc.dat, dnst.dat
c
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension rho(npntmx,nspin),vxc(npntmx,nspin),exc(npntmx,4)
c
      open(93,file='dnst.dat')
      open(94,file='exc.dat')
      open(95,file='vxc.dat')
      rewind(93)
      rewind(94)
      rewind(95)
      if (nspin.eq.2) then
        do 100 i=intpnt+1,npntmx
          dns=rho(i,1)+rho(i,2)
          vspint=(vxc(i,1)*rho(i,1)+vxc(i,2)*rho(i,2))/dns
          write(93,*) dns,(rho(i,is),is=1,nspin)
          write(94,*) (exc(i,1)+exc(i,2)+exc(i,3)+exc(i,4))/dns,
     + exc(i,1)/dns,exc(i,2)/dns,exc(i,3)/dns,exc(i,4)/dns
          write(95,*) vspint,(vxc(i,is),is=1,nspin)
  100   continue
      else
        do 200 i=intpnt+1,npntmx
          write(93,*) rho(i,1)
          write(94,*) (exc(i,1)+exc(i,2)+exc(i,3)+exc(i,4))/rho(i,1),
     + exc(i,1)/rho(i,1),exc(i,2)/rho(i,1),exc(i,3)/rho(i,1),
     + exc(i,4)/rho(i,1)
          write(95,*) vxc(i,1)
  200   continue
      endif
c
      close(93)
      close(94)
      close(95)
c
      return
      end
