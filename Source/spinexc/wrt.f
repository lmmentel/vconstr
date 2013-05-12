      subroutine wrtpot(rho,drho,exc,vxc,weight,intpnt,npnt,npntmx,
     + nspin)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      dimension vxc(npntmx,nspin,7),exc(npntmx,7),rho(npntmx,nspin),
     + drho(npntmx,3,nspin),weight(npnt)
c
      ek=0.d0
      eh=0.d0
      ec=0.d0
      ehks=0.d0
      ecks=0.d0
      ehcth=0.d0
c
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
        dnst=rho(ip,1)
        if (nspin.eq.2) dnst=dnst+rho(ip,2)
        ek=ek+dnst*vckin*weight(ip)
        eh=eh+dnst*vH*weight(ip)
        ec=ec+dnst*vcnd*weight(ip)
        ehks=ehks+dnst*vHks*weight(ip)
        ecks=ecks+dnst*vcndks*weight(ip)
        if (nspin.eq.2) then
          rha=rho(ip,1)
          drxa=drho(ip,1,1)
          drya=drho(ip,2,1)
          drza=drho(ip,3,1)
          rhb=rho(ip,2)
          drxb=drho(ip,1,2)
          dryb=drho(ip,2,2)
          drzb=drho(ip,3,2)
        else
          rha=5.d-1*rho(ip,1) 
          drxa=5.d-1*drho(ip,1,1)
          drya=5.d-1*drho(ip,2,1)
          drza=5.d-1*drho(ip,3,1)
          rhb=5.d-1*rho(ip,1)
          drxb=5.d-1*drho(ip,1,1)
          dryb=5.d-1*drho(ip,2,1)
          drzb=5.d-1*drho(ip,3,1)
        endif
        za=sqrt(drxa**2+drya**2+drza**2)
        zb=sqrt(drxb**2+dryb**2+drzb**2)
        zab=sqrt(drxa*drxb+drya*dryb+drza*drzb)
        call hcth(dfdra,dfdza,dfdrb,dfdzb,dfdzab,rha,rhb,
     + za,zb,zab,.true.,totalF_xc)
        ehcth=ehcth+totalF_xc*weight(ip)
      enddo
c     
      write(6,'(/6x,''hcthxc :'',12x,f10.6)')ehcth
c     
      write(6,'(/'' Kinetic correlation energy :'',f10.6)')ek
      write(6,'('' Coulomb repulsion energy   :'',f10.6)')5.d-1*eh
      write(6,'('' KS                         :'',f10.6)')5.d-1*ehks
      write(6,'('' XC hole energy             :'',f10.6)')
     + 5.d-1*(ec-eh)
      write(6,'('' X hole energy              :'',f10.6)')   
     + 5.d-1*(ecks-ehks)
      write(6,'('' 2-electron energy          :'',f10.6)')5.d-1*ec
c
      exc(1:npnt,3)=exc(1:npnt,1)+exc(1:npnt,3)
      exc(1:npnt,4)=exc(1:npnt,1)+exc(1:npnt,4)
      exc(1:npnt,5)=exc(1:npnt,2)+exc(1:npnt,5)
      exc(1:npnt,6)=exc(1:npnt,2)+exc(1:npnt,6)
c
      open(84,file='fexc.dat')
      open(85,file='fvxc.dat')
      rewind(84)
      rewind(85)
      do ip=intpnt+1,npnt
        read(90,*) vckin
        read(91,*) vH   
        read(92,*) vcnd
        read(93,*) vHks
        read(94,*) vcndks
        dnst=rho(ip,1)
        if (nspin.eq.2) dnst=dnst+rho(ip,2)
        vldax=0.d0
        vldac=0.d0
        vbx=0.d0
        vpwx=0.d0
        vpc=0.d0
        vpwc=0.d0
        vlyp=0.d0
        write(84,*) 5.d-1*(vcnd-vcndks-vH+vHks)+vckin,
     + 5.d-1*(vcndks-vHks),exc(ip,1)/dnst,exc(ip,2)/dnst,
     + exc(ip,3)/dnst,exc(ip,4)/dnst,exc(ip,5)/dnst,exc(ip,6)/dnst,
     + exc(ip,7)/dnst
        if (nspin.eq.2) then
          rha=rho(ip,1)
          drxa=drho(ip,1,1) 
          drya=drho(ip,2,1)
          drza=drho(ip,3,1)
          rhb=rho(ip,2)
          drxb=drho(ip,1,2)
          dryb=drho(ip,2,2)
          drzb=drho(ip,3,2)
        else
          rha=5.d-1*rho(ip,1)
          drxa=5.d-1*drho(ip,1,1)
          drya=5.d-1*drho(ip,2,1)
          drza=5.d-1*drho(ip,3,1)
          rhb=5.d-1*rho(ip,1)
          drxb=5.d-1*drho(ip,1,1)
          dryb=5.d-1*drho(ip,2,1)
          drzb=5.d-1*drho(ip,3,1)
        endif
        za=sqrt(drxa**2+drya**2+drza**2)
        zb=sqrt(drxb**2+dryb**2+drzb**2)
        zab=sqrt(drxa*drxb+drya*dryb+drza*drzb)
        call hcth(dfdra,dfdza,dfdrb,dfdzb,dfdzab,rha,rhb,
     + za,zb,zab,.true.,totalF_xc)
        do is=1,nspin
          vldax=vldax+vxc(ip,is,1)/2
          vldac=vldac+vxc(ip,is,2)/2
          vbx=vbx+vxc(ip,is,3)/2
          vpwx=vpwx+vxc(ip,is,4)/2
          vpc=vpc+vxc(ip,is,5)/2
          vpwc=vpwc+vxc(ip,is,6)/2
          vlyp=vlyp+vxc(ip,is,7)/2
        enddo
        write(85,*) vldax,vldac,vbx,vpwx,vpc,vpwc,vlyp
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
