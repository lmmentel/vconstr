      subroutine wrtpot(rho,exc,vxc,weight,intpnt,npnt,npntmx)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension rho(npntmx,1),exc(npntmx,7),vxc(npntmx,7),
     + weight(npnt)
c
      ek=0.d0
      eh=0.d0
      ec=0.d0
      ehks=0.d0
      ecks=0.d0
      open(90,file='vkin.dat')
      open(91,file='vhartr.dat')
      open(92,file='vcond.dat')
      open(93,file='vhartr.ks')
      open(94,file='vcond.ks')
      rewind(90)
      rewind(91)
      rewind(92)
      rewind(93)
      rewind(94)
      do ip=1,intpnt
        read(90,*) vckin
        read(91,*) vH
        read(92,*) vcnd
        read(93,*) vHks
        read(94,*) vcndks
        ek=ek+rho(ip,1)*vckin*weight(ip)
        eh=eh+rho(ip,1)*vH*weight(ip)
        ec=ec+rho(ip,1)*vcnd*weight(ip)
        ehks=ehks+rho(ip,1)*vHks*weight(ip)
        ecks=ecks+rho(ip,1)*vcndks*weight(ip)
      enddo
      write(6,'(/'' Kinetic correlation energy :'',f10.6)')ek
      write(6,'('' Coulomb repulsion energy   :'',f10.6)')5.d-1*eh
      write(6,'('' KS                         :'',f10.6)')5.d-1*ehks
      write(6,'('' XC hole energy             :'',f10.6)')
     + 5.d-1*(ec-eh)
      write(6,'('' X hole energy              :'',f10.6)')
     + 5.d-1*(ecks-ehks)
c
      exc(1:npnt,3)=exc(1:npnt,1)+exc(1:npnt,3)
      exc(1:npnt,4)=exc(1:npnt,1)+exc(1:npnt,4)
      exc(1:npnt,5)=exc(1:npnt,2)+exc(1:npnt,5)
      exc(1:npnt,6)=exc(1:npnt,2)+exc(1:npnt,6)
c
      open(83,file='dnst.dat')
      open(84,file='fexc.dat')
      open(85,file='fvxc.dat')
      rewind(83)
      rewind(84)
      rewind(85)
      do i=intpnt+1,npnt
        read(90,*) vckin
        read(91,*) vH
        read(92,*) vcnd
        write(83,*) rho(i,1)
        write(84,*) 5.d-1*(vcnd-vH)+vckin,exc(i,1)/rho(i,1),
     + exc(i,2)/rho(i,1),exc(i,3)/rho(i,1),exc(i,4)/rho(i,1),
     + exc(i,5)/rho(i,1),exc(i,6)/rho(i,1),exc(i,7)/rho(i,1)
        write(85,*) vxc(i,1),vxc(i,2),
     + vxc(i,3),vxc(i,4),vxc(i,5),vxc(i,6),vxc(1,7)
      enddo
      close(83)
      close(84)
      close(85)
      close(90)
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
c
      return
      end
