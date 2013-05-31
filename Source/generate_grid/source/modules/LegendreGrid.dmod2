module LegendreGridModule

   use vartypes

   implicit none 

   private
   public :: legpnt, legend, aslgd, legval

   real(KREAL), parameter :: zero=0.0_KREAL, one=1.0_KREAL, eps=1E-12_KREAL, half=0.5_KREAL
   real(KREAL), parameter :: conv=1E-12_KREAL, p125=0.125_KREAL, two=2.0_KREAL, three=3.0_KREAL
   real(KREAL), parameter :: pi= 3.141592653589793238462643383_KREAL

contains 


   subroutine legpnt (n, x, w)

      integer(KINT), intent(in)  :: n
      real(KREAL),   intent(out) :: x(n)
      real(KREAL),   intent(out) :: w(n)

   !  ======================================================================
   !  purpose:  determine all zeros of the legendre-polynomial of order n
   !            on the standard interval (-1, 1), with the associated
   !            weights (for gauss-legendre integration).
   !
   !  input  :  n      - order of the polynomial (= nr. of points/weights)
   !
   !  output :  x      - legendre points
   !            w      - weights of legendre points
   !
   !
   !  let the legendre polynomial be denoted by p(x).
   !
   !  for each zero a first guess xx is made with the help of a
   !  distribution-function for the zeros, to be given below.
   !  the correctness of the final results depends on this first
   !  guess. no guarantee can be given, but tests for n-values up to 10000
   !  produced correct results, and they indicate, that this first
   !  guess is better for larger n, which suggests, that the algorithm
   !  will also converge towards the correct x-value for larger n.
   !
   !  the weights, which depend on dp/dx, tend to be more inprecise for
   !  larger n, as the function fluctuates rapidly (near x=1), for large n.
   !  e.g. the error in the sum of (unnormalized) weights is 5e-13 for
   !  n=100, 7e-12 for n=1000, 7e-11 for n=10000.
   !  the routine normalizes the sum of weights to 2.
   !
   !  after the first guess, the value is updated iteratively by computing
   !  pn(x) and d/dx (pn(x)) at x=xx, and using the newton-raphson
   !  approximation (first-order taylor-expansion), until convergence is
   !  reached.
   !  the newton-raphson formula is
   !       x(new)    =    x(old)  -  f / (df/dx)
   !  for the problem:  find x : f(x) = 0
   !  where f(x) and df(x)/dx are to be evaluated at the current
   !  approximation of the zero (x-old)
   !
   !  accuracy: the convergence test-value is a reasonable indication
   !  of the accuracy of the weights. the x-values are much better.
   !  for large n (greater than 10,000 say) the routine should be made
   !  double precision, with an adapted convergence criterium,
   !  to suppress the errors in the weights.
   !
   !  for computing pn(x) and dp/dx the recursion formulas for the
   !  legendre-polynomials are used; they can be found in
   !  mathematical methods for physicists arfken.
   !     p(n)     = 1/n  *  ( (2n-1)*x*p(n-1) - (n-1)*p(n-2) )
   !    dp(n)/dx  = n* (p(n-1) - x*p(n)) / (1-x*x)
   !  ( p0(x)=1;  p1(x)=x;  p2(x)=0.5 * (3*x*x-1) )
   !
   !  finally (krylov approximate calculation of integrals), the weights
   !  are given by
   !     w       = 2 / ( (1-x**2) * (dp/dx)**2 )
   !  this can be rewritten in a lot of different forms, using recurrency-
   !  and other relations for the legendre-polynomials.
   !  as x is only approximately the sought-for zero, w is also only an
   !  approximation, and the accuracy depends on the formula choosen
   !  trial and error taught us to adopt
   !      w      = (1-x**2) / ( n**2 * (p1-x*p2)**2 )
   !  where p2=p(x), and p1=the value of the (n-1)-th order polynomial
   !  at x.
   !
   !  distribution of the zeros:
   !  the set of zeros (and weights) is symmetrical around x=0.
   !  the distribution-function of the zeros x(i),i=1,n/2 in the interval
   !  (0,1) is purely empirical, and given by
   !
   !  x(i)= sin( pi/2 * (2*i-1)/( n+.5+1/(8*f(n,i)) ) )
   !  f(n,i)=0.5 + n * cos( pi/2 * (2*i-1)/(n+.5) )    i=1,n/2   n is even
   !
   !  x(i)= sin( pi/2 * (2*i  )/( n+.5+1/(8*f(n,i)) ) )
   !  f(n,i)=0.5 + n * cos( pi/2 * (2*i  )/(n+.5) )    i=1,n/2   n is odd
   !
   !  ======================================================================

      integer(KINT) :: istart, jstart, i, j
      real(KREAL) :: phi, wtotal, p1, p2, hn, dphi, xx, correc, s

      if (n<=0) return

   !  ------------------------------------------------------------
   !  index i will run over the zeros, beginning with istart, the
   !  first zero in the interval (0,1).
   !  index j will run over the orders of the polynomials, in the
   !  recursion-formulas. the recursion will take two new terms at
   !  a time, so j must start with 1 or 2, depending on n being
   !  odd or even, respectively.
   !  phi and dphi are used for the first guess of the zeros.
   !  wtotal is the sum-of weights.
   !
   !
   !  preparation-part for odd/even n
   !  for odd n one of the zeros is x=0
   !  ------------------------------------------------------------

      if (mod(n,2)==0) then
         istart = n/2 + 1
         jstart = 4
         phi    = pi/2.0_KREAL
         wtotal = zero
      else
         istart      = (n+3)/2
         jstart      = 3
         phi         = pi
         x(istart-1) = zero
         p1          = - half
         do i = 4, n - 1, 2
            p1 = p1*(i-one)/i
         end do 
         w(istart-1) = two/(n*n*p1*p1)
         wtotal      = w(istart-1)
      end if
   
      hn   = n + half
      dphi = pi

   !  -------------------
   !  loop over the zeros
   !  -------------------

      do i = istart, n

   !     --------------------------------------
   !     the first guess of xx as the i-th zero
   !     --------------------------------------

         xx  = sin(phi/(hn+p125/(half+n*cos(phi/hn))))
         phi = phi + dphi

   !     ---------------------------------------------------------
   !     newton-raphson iteration.
   !     calculate p0 and p1 / p1 and p2, for n odd / even , resp.
   !     ---------------------------------------------------------

20       if (mod(n,2)==0) then
            p1 = xx
            p2 = half*(three*xx**2-one)
         else
            p1 = one
            p2 = xx
         end if

   !     -----------------------------------------------------
   !     calculate p(j-1) and p(j), until j=n, with recursion.
   !     -----------------------------------------------------

         do j = jstart, n, 2
            p1 = ((j+j-3)*xx*p2-(j-2)*p1)/(j-1)
            p2 = ((j+j-1)*xx*p1-(j-1)*p2)/j
         end do 

   !     ------------------------------------
   !     the newton-raphson correction on xx.
   !     test on convergence
   !     ------------------------------------

         correc = (xx**2-one)*p2/(n*(p1-xx*p2))
         xx     = xx + correc
         if (abs(correc)>conv) goto 20

   !     ----------------------
   !     final zero and weight.
   !     ----------------------

         x(i)   = xx
         s      = n*(p1-xx*p2)
         w(i)   = (one-xx*xx)*two/(s*s)
         wtotal = wtotal + w(i) + w(i)
   
      end do 

   !  --------------------------------------------
   !  deal with symmetry of the points and weights
   !  normalize the sum of weights to 2
   !  --------------------------------------------

      correc = two/wtotal
   
      do i = 1, (n+1)/2
         x(i)     = - x(n+1-i)
         w(n+1-i) = correc*w(n+1-i)
         w(i)     = w(n+1-i)
      end do 

   end subroutine 




   subroutine legend (n, point, weight, xmin, xmax)

      integer(KINT), intent(in)  :: n
      real(KREAL),   intent(out) :: point(n)
      real(KREAL),   intent(out) :: weight(n)
      real(KREAL),   intent(in)  :: xmin
      real(KREAL),   intent(in)  :: xmax

   !  ======================================================================
   !  purpose: calculate legendre points and weights on a given interval.
   !
   !  input  : n       nr. of points (and weights)
   !           xmin    lower boundary of interval
   !           xmax    upper boundary of interval
   !
   !  output : point   legendre points
   !           weight  weights of legendre points
   !
   !  method : the points and weights generated by routine legpnt are
   !           transformed from the standard interval (-1, 1) to the
   !           interval (xmin, xmax).
   !
   !  remark * if xmax<xmin the ordering of the points is such that
   !           the first points are still close to xmin. All weights are
   !           positive still
   !
   !  history: 1996.06.19, GtV: allow xmin>xmax: points in reverse order,
   !           but 'normal' (i.e.: positive) weights
   !
   !  ======================================================================

      integer(KINT) :: i
      real(KREAL) :: scale

      if (n>0) then
   
         call legpnt (n, point, weight)

   !     --------------------------------------------------------
   !     the scale-factor is the length of the interval, relative
   !     to that of the standard interval, which is 2.0
   !     --------------------------------------------------------

         scale = half*(xmax-xmin)
   
         do i = 1, n
            weight(i) = weight(i)*abs(scale)
            point(i)  = xmin + scale*(point(i)+one)
         end do 
      end if
   
   end subroutine




   subroutine legval (n, x, p)

      integer(KINT), intent(in) :: n
      real(KREAL),   intent(in) :: x
      real(KREAL),   intent(out) :: p(0:n)

   !  ==============================================================
   !  purpose:  calculate values of the legendre polynomials up to a
   !            given order.
   !
   !  input  :  n - maximum order of legendre polynomials
   !            x - argument of the polynomials
   !
   !  output :  p - polynomial values (0..n)
   !  ==============================================================

      integer(KINT) :: iorder 

      if (n>=0) then
         p(0) = one
         if (n>=1) then
            p(1) = x
            do iorder = 2, n
               p(iorder) = ((2*iorder-1)*x*p(iorder-1)-(iorder-1)*p(iorder-2))/iorder
            end do
         end if
      end if

   end subroutine




   subroutine aslgd (nmax, x, asl)

      integer(KINT), intent(in)  :: nmax
      real(KREAL),   intent(in)  :: x
      real(KREAL),   intent(out) :: asl((nmax+1)*(nmax+2)/2)

   !  ======================================================================
   !  purpose: calculate value of associated legendre function for x
   !  ======================================================================

      integer(KINT) :: n, ixlg, m, nua, nub, ifn, ifnn, ixlgs

      asl(1) = one
      asl(2) = x
   
      if (one-abs(x)<=eps) then
   
         asl(3) = zero
         ixlg   = 3
         do n = 2, nmax
            ixlg      = ixlg + 1
            asl(ixlg) = one
            do m = 1, n
               ixlg      = ixlg + 1
               asl(ixlg) = zero
            end do 
         end do 
   
      else
   
         asl(3) = - sqrt(one-x*x)
         ifn    = 2
         ifnn   = 3
         ixlgs  = 4
   
         do n = 2, nmax
            ixlg      = ixlgs
            nua       = ixlgs - n
            nub       = nua - n
            asl(ixlg) = (ifnn*x*asl(nua)-(ifn-1)*asl(nub+1))/ifn
            ixlg      = ixlg + 1
            nua       = nua + 1
            nub       = nub + 1
            asl(ixlg) = ifnn*asl(3)*asl(nua-1)
            if (n/=2) asl(ixlg) = asl(nub+1) + asl(ixlg)
   
            do m = 2, n
               ixlg      = ixlg + 1 
               nua       = nua + 1
               nub       = nub + 1
               asl(ixlg) = - ((n-m+1)*(n-m+2))*asl(ixlg-2) + ((n+m-2)*(n+m-3))*asl(nub-1)
               if (n-2>=m) asl(ixlg) = asl(nub+1) + asl(ixlg)
            end do 
   
            ifn   = ifn + 1
            ifnn  = ifnn + 2
            ixlgs = ixlgs + n + 1
         end do 

      end if

   end subroutine

end module


