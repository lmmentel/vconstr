      program cdat
      implicit double precision (a-h,o-z)
c
      dz=2.d-3
c
      open(85,file='pnts')
      rewind(85)
      z=-.001d0
      do ip=1,500
        z=z+dz
        write(85,'(4(2x,f12.6))')0.d0,0.d0,z,0.d0
      enddo
c
      end

