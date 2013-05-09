      subroutine clcmod(valmo,npnt,npntmx,norb,vmopao,grid,weight,
     + atmol4)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      parameter(eps  = 1.d-6)
c
      logical lorb,atmol4
      dimension valmo(npnt*norb)
      dimension vmopao(norb*norb)
      dimension grid(npntmx,3),weight(npnt)
      dimension chkvl(norb)
      dimension valao(norb)
c
      chkvl(1:norb)  = 0.d0
c
      do 100 ipnt=1,npnt
        k=(ipnt-1)*norb
        m=k+1
        x=grid(ipnt,1)
        y=grid(ipnt,2)
        z=grid(ipnt,3)
        w=weight(ipnt)
        if (atmol4) then
          call aovlw(x,y,z,valao)
        else
          call aovlv(x,y,z,valao)
        endif
        chkvl(1:norb)=chkvl(1:norb)+w*valao(1:norb)**2
        call vecmat(valao,vmopao,norb,valmo(m))
  100 continue
c
      lorb=.true.
      write(6,'(/'' Orbital errors upon integration'')')
      do i=1,norb
        t1=abs(chkvl(i)-1.d0)
        if (t1.gt.eps) then
          write(6,'(i3,1x,g8.2)')i,t1
          lorb=.false.
        endif
      enddo
      if (lorb) write(6,'('' No errors greater than '',g12.2)')eps
c
      return
      end
