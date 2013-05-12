      subroutine vbrain(nhmx,nprmx,ngrpx,norb,ncore,nvpr,ip2un,
     + idmp,ismo,isao,rnel,lvhart,lvcond,lvexch,atmol4)
c
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(half = 5.d-1)
c
c** up to g-functions
c
      parameter(lqntmx = 4)
      parameter(ncmpx = 1+3+6+10+15)
      parameter(l1l2mx = 2*lqntmx)
c
      logical atmol4,lvhart,lvcond,lvexch
      dimension ncomp(0:lqntmx),ncbas(0:lqntmx),ll(ncmpx),mm(ncmpx),
     1 nn(ncmpx),pao(norb*(norb+1)/2),pmo(norb*(norb+1)/2),
     2 pcndao(norb*(norb+1)/2),pxchao(norb*(norb+1)/2),valao(norb),
     3 valmo(norb),vmoao(norb*norb),scrtch(norb*(norb+1)/2),
     4 xintmx(norb*(norb+1)/2),prefac(nprmx)
c
      allocatable hxyz(:,:)
c
      data ncomp/1,3,6,10,15/
      data ncbas/0,1,4,10,20/
      data ll/0, 1,0,0, 2,0,0,1,1,0, 3,0,0,2,2,1,0,1,0,1,
     +        4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data mm/0, 0,1,0, 0,2,0,1,0,1, 0,3,0,1,0,2,2,0,1,1,
     +        0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data nn/0, 0,0,1, 0,0,2,0,1,1, 0,0,3,0,1,0,1,2,2,1,
     +        0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
c
      allocate(hxyz(nhmx,3),stat=ialloc)
      if (ialloc.eq.0) then 
        write(6,'(/'' Memory allocation successfull '')')
      else
        write(6,'(/'' No memory allocation '')')
        stop
      endif
c
c** initialise files, and read number of points
c
      read(10,*)npnt
      write(6,'(/,'' # of grid points       : '',i6,//)')npnt
      open(11,file='vhartr.dat')
      open(12,file='vcond.dat')
      open(13,file='vexch.dat')
      rewind(11)
      rewind(12)
      rewind(13)
c
      if (ncore.eq.0) then
        ip2unw = ip2un
        call p1dens(ip2un,pmo,norb,rnel)
      else
        ip2unw = 1
        print*,'check implementation of core-routines'
        stop
c        norb = norb+ncore
c        rnel = rnel+2.d0*ncore
c
c*** calculate density from core-expanded two-matrix
c
c        call p1dens(ip2unw,pmo,norb,rnel)
      endif
c
c** get gaussians and HF-vmoaos from dumpfile
c
      call rdsamo(idmp,ismo,isao,vmoao,norb,rnel)
c
c** transform density from mo to ao basis
c
      call tmtdag(pmo,norb,pao,norb,vmoao,half)
c
c** calculate some constants
c
      call fctrl
      call cnstnd
c
c*** calculate hx,hy and hz for every pair of gaussians
c
      if (atmol4) then
c
c* calculate transformation matrices for transformation from
c* cartesian to spherical gaussians.  stored in common/spheri/.
c
        call spherw
        call calchw(hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,prefac,
     + nprmx,l1l2mx)
      else
c
c*** idddd = 0 : 6-d set is used;  idddd<>0 : 5-d set is used.
c
        call dset(idddd)
        call spherv(idddd)
        call calchv(hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,prefac,
     + nprmx,l1l2mx)
      endif
c
c*** loop over points
c
      do 100 ipnt = 1,npnt
c
        vhartr = 0.d0
        vcond = 0.d0
        vexch = 0.d0
c
        read(10,*) xpnt,ypnt,zpnt
c
c*** calculate integrals on ao-basis.  store in array xintmx
c
        if (atmol4) then
          call aovlw(xpnt,ypnt,zpnt,valao)
          call ccintw(xintmx,norb,ncomp(0),ncbas(0),lqntmx,ll,mm,nn,
     + ncmpx,hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,xpnt,ypnt,zpnt,prefac,
     + nprmx,l1l2mx,ngrpx)
        else
          call aovlv(xpnt,ypnt,zpnt,valao)
          call ccintv(xintmx,norb,ncomp(0),ncbas(0),lqntmx,ll,mm,nn,
     + ncmpx,hxyz(:,1),hxyz(:,2),hxyz(:,3),nhmx,xpnt,ypnt,zpnt,prefac,
     + nprmx,l1l2mx,ngrpx)
        endif
        call vecmat(valao,vmoao,norb,valmo)
c
        if (lvhart.and.lvcond.and.lvexch) then 
          call cond(ip2unw,scrtch,norb,rnel,valmo)
          call tmtdag(scrtch,norb,pcndao,norb,vmoao,half)
          call exchan(pmo,norb,scrtch,valmo)
          call tmtdag(scrtch,norb,pxchao,norb,vmoao,half)
          do i = 1,norb
            do j = 1,i
              ij = i*(i-1)/2+j
              vhartr = vhartr+pao(ij)*xintmx(ij)
              vcond = vcond+pcndao(ij)*xintmx(ij)
              vexch = vexch+pxchao(ij)*xintmx(ij)
            enddo
          enddo
          write(11,*) vhartr
          write(12,*) vcond
          write(13,*) vexch
        elseif (lvhart.and.lvcond) then
          call cond(ip2unw,scrtch,norb,rnel,valmo)
          call tmtdag(scrtch,norb,pcndao,norb,vmoao,half)
          do i = 1,norb
            do j = 1,i
              ij = i*(i-1)/2+j
              vhartr = vhartr+pao(ij)*xintmx(ij)
              vcond = vcond+pcndao(ij)*xintmx(ij)
            enddo
          enddo
          write(11,*) vhartr
          write(12,*) vcond
        elseif (lvhart) then
          do i = 1,norb
            do j = 1,i
              ij = i*(i-1)/2+j
              vhartr = vhartr+pao(ij)*xintmx(ij)
            enddo
          enddo
          write(11,*) vhartr
        elseif (lvcond) then
          call cond(ip2unw,scrtch,norb,rnel,valmo)
          call tmtdag(scrtch,norb,pcndao,norb,vmoao,half)
          do i = 1,norb
            do j = 1,i
              ij = i*(i-1)/2+j
              vcond = vcond+pcndao(ij)*xintmx(ij)
            enddo
          enddo
          write(12,*) vcond
        elseif (lvexch) then
c
c*** (this only works well for densities with equal alpha and beta part)
c

          call exchan(pmo,norb,scrtch,valmo)
          call tmtdag(scrtch,norb,pxchao,norb,vmoao,half)
          do i = 1,norb
            do j = 1,i
              ij = i*(i-1)/2+j
              vexch = vexch+pxchao(ij)*xintmx(ij)
            enddo
          enddo
          write(13,*) vexch
        endif
c
        nmbr = mod(ipnt,nvpr)
        if (nmbr.eq.0) then
          call givtim(top,bot)
          write(6,'(''   '',i6,'' points calculated at '',f13.3,
     + '' wall '',f13.3,'' secs'')')ipnt,top,bot
          write(6,'(''vhartr '',g8.2,''vcond '',g8.2,
     + '' vexch '',g8.2)')vhartr,vcond,vexch
        endif
c
  100 continue
c
      close(11)
      close(12)
      close(13)
c
      return
      end
