      subroutine copy(a,b,n)
c
c copy n elements of array a to array b 
c
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension a(n),b(n)
      do 10 i = 1,n
        b(i) = a(i)
   10 continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dcopy(n,dx,incx,dy,incy)
      implicit real*8  (a-h,o-z)
c
c     copies a vector, dx, to a vector, dy.
c
      dimension dx(*),dy(*)
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 20 i = 1,n
        dy(iy) = dx(ix)
        ix = ix+incx
        iy = iy+incy
   20 continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine clearr(n,a)
      implicit real*8 (a-h,o-z),integer(i-n)
      parameter(zero = 0.0d0)
      dimension a(*)
      do 110 i=1,n
        a(i)=zero
 110  continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine cleari(n,m)
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension m(*)
      do 120 i=1,n
        m(i)=0
 120  continue
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine clearl(n,l)
      implicit real*8 (a-h,o-z),integer(i-n)
      logical l
      dimension l(*)
      do 110 i=1,n
        l(i)=.false.
 110  continue
      return
      end
