      subroutine cbrain(fxyz,npnt,norb,idmp,ismo,isno,isao,atmol4)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z),integer(i-n)
c
      logical atmol4
c
      dimension fxyz(3)
      dimension pmo(norb*(norb+1)/2),pnomo(norb*(norb+1)/2)
      dimension vmoao(norb*norb), vmopao(norb*norb)
      dimension vnuc(npnt)
      dimension grid(npnt,3),weight(npnt)
      dimension dns(npnt),ddns(npnt,3),dsdns(npnt)
      dimension valmo(norb*npnt)
      allocatable grdmo(:,:)
c
      npntmx=npnt
c
      allocate(grdmo(norb*npnt,4),stat=ialloc)
      if (ialloc.ne.0) then
        write(6,'(/'' No memory allocated for derivativess'')')
        stop
      endif
c
      do ipnt=1,npnt
        read(99,*) grid(ipnt,1),grid(ipnt,2),grid(ipnt,3),weight(ipnt)
      enddo
      close(99)
c
c** read dumpfile : vmopao - molecular orbitals in primitive ao basis
c**                 vmoao  - molecular orbitals in ao basis
c**                 pmo    - mo density matrix in mo basis
c**                 pnomo  - ci density matrix in mo basis
c
      call rddmp(vmopao,vmoao,pmo,pnomo,norb,nvpr,rnel,idmp,ismo,isao,
     + isno)
c
      call clcvls(dns,ddns,dsdns,vnuc,valmo,grdmo,npnt,npntmx,
     + norb,pnomo,vmopao,grid,weight,rnel,atmol4)
c
      open(45,file='vhartr.dat')
      rewind(45)
      do ip=1,npnt
        read(45,*) vhartr
        print*,dsdns(ip)/(4*dns(ip))-
     + (ddns(ip,1)**2+ddns(ip,2)**2+ddns(ip,3)**2)/(8*dns(ip)**2)+
     + vnuc(ip)-vhartr+fxyz(1)*grid(ip,1)+fxyz(2)*grid(ip,2)+
     + fxyz(3)*grid(ip,3)
      enddo
c
      return
      end
