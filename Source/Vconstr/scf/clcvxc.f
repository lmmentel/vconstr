      subroutine clcvxc(nspin,npntmx,npnt,rho,drho,d2rho,vxc,exc,
     + lb88,lpw91x,lprdwc,lpw91c,llyp)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical lcalcv,lb88,lpw91x,lprdwc,lpw91c,llyp
      dimension rho(npntmx,nspin),drho(npntmx,3,nspin),
     + d2rho(npntmx,6,nspin),vxc(npntmx,nspin),exc(npntmx,4)
c
      logical nilrho
      dimension nilrho(npnt,3)
      dimension frho(npnt,nspin),fdrho(npnt,3,nspin),
     + fd2rho(npnt,6,nspin),rhotrd(npnt,nspin),relrho(npnt,nspin),
     + sum(npnt),dsum(npnt,3),d2sum(npnt,6),sumtrd(npnt),
     + sumsxt(npnt)

c
      exc(1:npnt,1:4)=0.d0
      vxc(1:npnt,1:nspin)=0.d0
c
c     ================================
c     COMPUTATION OF ENERGY, POTENTIAL
c     ================================
c
      if (llyp) then
        xcpar=2d0/3d0
        call xcvpot(nspin,npnt,npntmx,rho,0,xcpar,vxc)
        call xcener(nspin,npnt,npntmx,rho,0,xcpar,exc(1,1),
     + exc(1,2))
      else
        dum=0.d0
        call xcvpot(nspin,npnt,npntmx,rho,1,dum,vxc)
        call xcener(nspin,npnt,npntmx,rho,1,dum,exc(1,1),
     + exc(1,2))
      endif
c
c     ===============================
c     PREPARATION FOR THE CORRECTIONS
c     ===============================
c
      lcalcv = .true.
c
      call xcprep(lcalcv,nspin,npntmx,npnt,rho,drho,d2rho,nilrho,
     + frho,fdrho,fd2rho,rhotrd,relrho,sum,dsum,d2sum,sumtrd,sumsxt)
c
c     -------------------
c     EXCHANGE CORRECTION
c     -------------------
c
c
      if (lb88) then
c
        call xcbeck(lcalcv,nspin,npntmx,npnt,frho,fdrho,fd2rho,rhotrd,
     + nilrho,exc(1,3),vxc)
c
      elseif (lpw91x) then
c
        call pw91x(lcalcv,nspin,npntmx,npnt,rho,drho,d2rho,nilrho,
     + exc(1,3),vxc)
c
      endif
c
c     ----------------------
c     CORRELATION CORRECTION
c     ----------------------
c
      if (lprdwc) then
c
        call xcprdw(lcalcv,nspin,npntmx,npnt,frho,fdrho,rhotrd,
     + relrho,sum,dsum,d2sum,sumtrd,sumsxt,nilrho,exc(1,4),vxc)
c
      elseif (lpw91c) then
c
        call pw91c(lcalcv,nspin,npntmx,npnt,rho,drho,d2rho,nilrho,
     + exc(1,4),vxc)
c
      elseif (llyp) then
c
        call lyp(lcalcv,nspin,npntmx,npnt,frho,fdrho,fd2rho,relrho,
     + sum,dsum,d2sum,sumtrd,nilrho,exc(1,4),vxc)
c
      endif
c
      return
      end
