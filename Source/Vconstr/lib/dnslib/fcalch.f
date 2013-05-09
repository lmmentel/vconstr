      subroutine g(l1,l2,a,b,eps,gl)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/factor/fac(0:40),faci(0:40)
      dimension gl(0:l1+l2)
c
      gl(0:l1+l2)=0.d0
      do l = 0,l1+l2
        iqhigh = l/2
        aux1 = aux(l,l1,l2,a,b)*fac(l)
        xeps = 1.d0
        do iq = 0,iqhigh
          ldex = l-2*iq
          gl(ldex) = gl(ldex)+faci(iq)*xeps*aux1
          xeps = xeps*eps
        enddo
      enddo
c
      return
      end
c
      function aux(k,l1,l2,a,b)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/factor/fac(0:40),faci(0:40)
c
      sum = 0.d0
      a1 = a
      if (abs(a1).lt.1.d-14) a1 = 1.d-14
      b1 = b
      if (abs(b1).lt.1.d-14) b1 = 1.d-14
      do i = 0,l1
        x = a1**(l1-i)*fac(l1)*faci(i)*faci(l1-i)
        do j = 0,l2
          if ((i+j).eq.k) sum = sum+
     + x*b1**(l2-j)*fac(l2)*faci(j)*faci(l2-j)
        enddo
      enddo
      aux = sum
c
      return
      end
c
      subroutine h(gl,eps,hh,l1l2)
      implicit real*8 (a-h,o-z),integer(i-n)
c
      common/factor/fac(0:40),faci(0:40)
      dimension gl(0:l1l2),hh(0:l1l2,0:l1l2)
c
      lll = l1l2+1
      goto(1,2,3,4,5,6,6,6,6)lll
c
    6 hh(0:l1l2,0:l1l2)=0.d0
      do l = 0,l1l2
        ithigh = l/2
        do it = 0,ithigh
          i = l-it
          labda = i-it
          hh(i,labda) = hh(i,labda)+gl(l)*(-eps)**it*
     + faci(it)*faci(labda)
        enddo
      enddo
      return
c
    5 hh(4,4) = gl(4)*0.0416666666667d0
      hh(3,2) = 5.d-1*gl(4)*(-eps)
      hh(2,0) = 5.d-1*gl(4)*eps*eps
    4 hh(3,3) = gl(3)*0.1666666666667d0
      hh(2,1) = gl(3)*(-eps)
    3 hh(2,2) = 5.d-1*gl(2)
      hh(1,0) = gl(2)*(-eps)
    2 hh(1,1) = gl(1)
    1 hh(0,0) = gl(0)
      return
c
      end
