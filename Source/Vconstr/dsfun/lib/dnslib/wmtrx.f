      subroutine wmtrx(a,p,ndim,eps)
c *****************************************************************
c
c   print matrix p, dimension ndim
c
c *****************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      character*40 a
      dimension p(ndim*(ndim+1)/2)
      dimension i(4),j(4),val(4)
      write(6,80)a
   80 format(/'  --------------  ',a40,' -------------',/)
      n = 1
      do 10 k = 1,ndim
        do 5 l = 1,k
          i(n) = k
          j(n) = l
          val(n) = p(k*(k-1)/2+l)
          if (abs(val(n)).gt.eps)  n = n + 1
          if (n.eq.5) then
            write(6,102)(i(iz),j(iz),val(iz),iz=1,4)
            n=1
          endif
    5   continue
   10 continue
  102 format(4(i5,i4,g14.6))
      write(6,102)(i(iz),j(iz),val(iz),iz=1,n-1)
      write(6,103)
  103 format(/,2x,70('-'))
c
      return
      end
