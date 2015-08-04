      function fnusp(m,z)
c***********************************************************************
c***** see j. chem. phys. (1971) p. 1965
c***** for z.le.zmax, fmch(m,z) is given by ...
c***** ((a(1)+a(2)*z+ ... +a(n)*z**(n-1))/(1+b(1)*z+ ... +b(n)   cont ..
c*****                      *z**n))**(m+0.5)
c***** for z.gt.zmax, fmch(m,z) is given by ...
c***** df(2m-1)*(pi/2)**0.5/(2z)**(m+0.5), where df is double factorial
c***** n,a(2)-a(n),b(1)-b(n),zmax tabulated in reference
c **** a(1) is accurate to 12 figures ... (2m+1)**(-2/(2m+1))
c **** 1.2533... is (pi/2)**0.5 to 12 figures
c***********************************************************************
      implicit real*8 (a-h,o-z),integer(i-n)
      common/fconst/a(0:9,6),b(0:9,6),df(0:9),n(0:9),zmax(0:9)
c
      if (z.gt.zmax(m)) then
        z2 = 2.d0*z
        x = sqrt(z2)
      else
        t = a(m,1)
        s = 1.d0+b(m,1)*z
        x = z
        do i = 2,n(m)
          t = t+a(m,i)*x
          x = x*z
          s = s+b(m,i)*x
        enddo
        z2 = t/s
        x = sqrt(z2)
      endif
      m1 = m+1
      goto(1,2,3,4,5,6,7,8,9,10)m1
   10 x =  x*z2
    9 x =  x*z2
    8 x =  x*z2
    7 x =  x*z2
    6 x =  x*z2
    5 x =  x*z2
    4 x =  x*z2
    3 x =  x*z2
    2 x =  x*z2
    1 if (z.gt.zmax(m)) then
        fnusp = df(m)*1.253314137316d0/x
      else
        fnusp = x
      endif
c
      return
      end
c
      subroutine xyzfac(h,l1l2,a,x)
      implicit real*8 (a-h,o-z),integer(i-n)
      dimension x(0:l1l2),h(0:l1l2,0:l1l2)
      dimension aa(0:l1l2)
c
      aa(0) = 1.d0
      do l = 1,l1l2
        aa(l) = aa(l-1)*a
      enddo
c
      ll = l1l2+1
      goto (1,2,3,4,5,6,6,6,6)ll
c
    1 x(0) = h(0,0)
      return
    2 x(0) = h(0,0)
      x(1) = h(1,1)*aa(1)
      return
    3 x(0) = h(0,0)
      x(1) = h(1,0)+h(1,1)*aa(1)
      x(2) = h(2,2)*aa(2)
      return
    4 x(0) = h(0,0)
      x(1) = h(1,0)+h(1,1)*aa(1)
      x(2) = h(2,1)*aa(1)+h(2,2)*aa(2)
      x(3) = h(3,3)*aa(3)
      return
    5 x(0) = h(0,0)
      x(1) = h(1,0)+h(1,1)*aa(1)
      x(2) = h(2,0)+h(2,1)*aa(1)+h(2,2)*aa(2)
      x(3) = h(3,2)*aa(2)+h(3,3)*aa(3)
      x(4) = h(4,4)*aa(4)
      return
    6 do l = 0,l1l2
        x(l) = 0.d0
        ithigh = l/2
        do it = 0,ithigh
          i = l-it
          labda = i-it
          x(i) = x(i)+h(i,labda)*aa(labda)
        enddo
      enddo
c
      return
      end
