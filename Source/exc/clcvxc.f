      subroutine clcvxc(npntmx,npnt,rho,drho,d2rho,vxc,exc)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(nspin  = 1)
c
      dimension rho(npntmx,nspin),drho(npntmx,3,nspin),
     + d2rho(npntmx,6,nspin),vxc(npntmx,7),exc(npntmx,7)
c
      logical lcalcv,nilrho
      dimension nilrho(npnt,3)
      dimension frho(npnt,nspin),fdrho(npnt,3,nspin),
     + fd2rho(npnt,6,nspin),rhotrd(npnt,nspin),relrho(npnt,nspin),
     + sum(npnt),dsum(npnt,3),d2sum(npnt,6),sumtrd(npnt),
     + sumsxt(npnt)
c
      exc(1:npnt,1:7)=0.d0
      vxc(1:npnt,1:7)=0.d0
c
c     ===============================
c     PREPARATION FOR THE CORRECTIONS
c     ===============================
c
      lcalcv = .true.
c
c     ================================
c     COMPUTATION OF ENERGY, POTENTIAL
c     ================================
c
      dum=1.d0
      call xcpot(nspin,npnt,npntmx,rho,1,dum,vxc(1,1),vxc(1,2))
      call xcener(nspin,npnt,npntmx,rho,1,dum,exc(1,1),exc(1,2))
c
      call xcprep(lcalcv,nspin,npntmx,npnt,rho,drho,d2rho,nilrho,
     + frho,fdrho,fd2rho,rhotrd,relrho,sum,dsum,d2sum,sumtrd,sumsxt)
c
      call xcbeck(lcalcv,nspin,npntmx,npnt,frho,fdrho,fd2rho,rhotrd,
     + nilrho,exc(1,3),vxc(1,3))
c
      call pw91x(lcalcv,nspin,npntmx,npnt,rho,drho,d2rho,nilrho,
     + exc(1,4),vxc(1,4))
c
      call xcprdw(lcalcv,nspin,npntmx,npnt,frho,fdrho,rhotrd,
     + relrho,sum,dsum,d2sum,sumtrd,sumsxt,nilrho,exc(1,5),vxc(1,5))
c
      call pw91c(lcalcv,nspin,npntmx,npnt,rho,drho,d2rho,nilrho,
     + exc(1,6),vxc(1,6))
c
      call lyp(lcalcv,nspin,npntmx,npnt,frho,fdrho,fd2rho,relrho,
     + sum,dsum,d2sum,sumtrd,nilrho,exc(1,7),vxc(1,7))
c
      return
      end
