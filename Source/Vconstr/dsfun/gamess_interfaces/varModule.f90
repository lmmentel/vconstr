module VarModule
  use ISO_FORTRAN_ENV
  implicit none
!	----------------------
!	Define KIND parameters
! ----------------------
  integer,     parameter :: IK = kind(1)
  integer(IK), parameter :: LK2 = kind(2)
  integer(IK), parameter :: LK8 = kind(8)
  integer(IK), parameter :: SP = kind(1.0)
  integer(IK), parameter :: DP = kind(1.0D0)
! ----------------------------
! Specify kinds used in Gamess
! ----------------------------
  integer(IK), parameter :: IPgamess = 8  !4 = 32 bit integer, 8 = 64 bit integer
  integer(IK), parameter :: DPgamess = DP
  integer(IK), parameter :: LKgamess = 8
! -------------------------------
! Define length of default string
! -------------------------------
  integer(IK), parameter :: LCHARS = 60
!	---------------------------------------------
!	Set std unit numbers to system specified ones
!	---------------------------------------------
  integer(IK), parameter :: stdin  = INPUT_UNIT
  integer(IK), parameter :: stdout = OUTPUT_UNIT 
  integer(IK), parameter :: stderr = ERROR_UNIT
  integer(IK), parameter :: twoein = 9
!	---------
!	Constants
!	---------
  real(DP), parameter :: PI    = 3.1415926535897932384626433_DP
  real(DP), parameter :: TWOPI = PI+PI
  real(DP), parameter :: GOLD  = 1.6180339887498948482045868_DP
  real(DP), parameter :: CGOLD = 0.3819660112501051517954132_DP
! ------------------------------  
! some common format definitions
! ------------------------------  
 character(len=20), parameter :: plain   = '(a100)' 
 character(len=20), parameter :: mofmt   = '(i6,es22.12)' 
 character(len=20), parameter :: orbefmt = '(a4,f10.6)' 
 character(len=20), parameter :: gmsvfmt = '(i2,i3,5(es15.8))' 
! ------------------
! Define vector norm
! ------------------
  interface Norm2
    module procedure Norm2Real
    module procedure Norm2Complex
  end interface

  private Norm2Real, Norm2Complex

! -------------------------------------
! Interface to Gamess' stopping routine
! -------------------------------------
  interface
    subroutine abrt()
    end subroutine
  end interface

  interface matTrans
!      module procedure functionMatrixTransform
      module procedure subroutineMatrixTransform
      module procedure subroutineTriangularMatrixTransform
  end interface

  interface inverseMatrix
      module procedure LU_inv
  end interface

  private :: functionMatrixTransform, subroutineMatrixTransform, LU_inv

  public :: inverseMatrix

contains

  real(DP) function Norm2Real(vec)
    real(DP), intent(in) :: vec(:)
    Norm2Real = dot_product(vec,vec)
  end function
  real(DP) function Norm2Complex(vec)
    complex(DP), intent(in) :: vec(:)
    Norm2Complex = real(dot_product(vec,vec), kind=DP) !dot_product already takes conjugate
  end function

  function OuterProduct(vec1, vec2) result(matrix)
    real(DP), intent(in) :: vec1(:), vec2(:)
    real(DP)             :: matrix(size(vec1),size(vec2))
    matrix = spread(vec1, 2, size(vec2)) * spread(vec2, 1, size(vec1))
  end function

  elemental complex(DP) function ImExp(x)
    real(DP), intent(in) :: x

    ImExp = exp(cmplx(0.0_DP, x, kind=DP))
  end function

  elemental real(DP) function phase(re, im)
    real(DP), intent(in) :: re, im

    if (re > epsilon(1.0_DP)) then
      phase = atan(im/re)
    elseif (re < -epsilon(1.0_DP)) then
      if (im < 0.0_DP) then
        phase = atan(im/re) - PI
      else
        phase = atan(im/re) + PI
      endif
    else
      if (im > 0.0_DP) then
        phase =  0.5_DP*PI
      elseif (im < 0.0_DP) then
        phase = -0.5_DP*PI
      else
        phase =  0.0_DP
      endif
    endif
  end function

   function functionMatrixTransform(X,A,Y) result(R)
! R = X^T * A * Y
    real(dp), dimension(:,:), intent(in)  :: X,Y,A
    real(dp)                              :: R(size(A,1),size(A,2)) 
    real(dp), dimension(:,:), allocatable :: tmp 

  allocate(tmp(size(A,1),size(A,2)))

  tmp = matmul(A,Y)
  R   = matmul(transpose(X),tmp)
 
  deallocate(tmp)
  end function FunctionMatrixTransform
  
  subroutine SubroutineMatrixTransform(X,A,Y) 
! A = X * A * Y^T
    real(dp), dimension(:,:), intent(in)    :: X,Y 
    real(dp), dimension(:,:), intent(inout) :: A
    real(dp), dimension(:,:), allocatable   :: tmp 

  allocate(tmp(size(A,1),size(A,2)))

  tmp = matmul(A,transpose(Y))
  A   = matmul(X,tmp)

  deallocate(tmp)
  end subroutine SubroutineMatrixTransform

  subroutine SubroutineTriangularMatrixTransform(A, X) 
! A = X^T * A * Y
  real(dp), dimension(:,:), intent(in)    :: X 
  real(dp), dimension(:),   intent(inout) :: A
  real(dp), dimension(:),   allocatable   :: Atmp 
  integer(ik) :: i, j, k, l, ij, kl

  allocate(Atmp(size(A)))
  
  Atmp = A
  ij = 0
  do i = 1, size(X, 1)
    do j = 1, i
      ij = ij + 1
      A(ij) = 0.0_dp
      do k = 1, size(X, 1)
        do l = 1, size(X, 1)
          kl = (max(l,k)*(max(l,k)-1))/2 + min(l,k)
          A(ij) = A(ij) + X(k,i)*Atmp(kl)*X(l,j)
        enddo
      enddo
    enddo
  enddo

  deallocate(Atmp)
  end subroutine SubroutineTriangularMatrixTransform

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
    IMPLICIT NONE
    !Declarations
    INTEGER(ik), INTENT(IN) :: n
    INTEGER(ik), INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
    REAL(dp), INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
    REAL(dp), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
     
    LOGICAL :: FLAG = .TRUE.
    INTEGER(ik) :: i, j, k, l
    REAL(dp) :: m
    REAL(dp), DIMENSION(n,2*n) :: augmatrix !augmented matrix
     
    !Augment input matrix with an identity matrix
    DO i = 1, n
        DO j = 1, 2*n
            IF (j <= n ) THEN
                augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
                augmatrix(i,j) = 1
            Else
                augmatrix(i,j) = 0
            ENDIF
        END DO
    END DO
     
    !Reduce augmented matrix to upper traingular form
    DO k =1, n-1
        IF (augmatrix(k,k) == 0) THEN
            FLAG = .FALSE.
            DO i = k+1, n
                IF (augmatrix(i,k) /= 0) THEN
                    DO j = 1,2*n
                        augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                    END DO
                    FLAG = .TRUE.
                    EXIT
                ENDIF
                IF (FLAG .EQV. .FALSE.) THEN
                    PRINT*, "Matrix is non - invertible"
                    inverse = 0
                    errorflag = -1
                    return
                ENDIF
            END DO
        ENDIF
        DO j = k+1, n           
            m = augmatrix(j,k)/augmatrix(k,k)
            DO i = k, 2*n
                augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
        END DO
    END DO
     
    !Test for invertibility
    DO i = 1, n
        IF (augmatrix(i,i) == 0) THEN
            PRINT*, "Matrix is non - invertible"
            inverse = 0
            errorflag = -1
            return
        ENDIF
    END DO
     
    !Make diagonal elements as 1
    DO i = 1 , n
        m = augmatrix(i,i)
        DO j = i , (2 * n)              
               augmatrix(i,j) = (augmatrix(i,j) / m)
        END DO
    END DO
     
    !Reduced right side half of augmented matrix to identity matrix
    DO k = n-1, 1, -1
        DO i =1, k
        m = augmatrix(i,k+1)
            DO j = k, (2*n)
                augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
        END DO
    END DO              
     
    !store answer
    DO i =1, n
        DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
        END DO
    END DO
    errorflag = 0
END SUBROUTINE FINDinv


  subroutine LU_inv(A)
  real(dp), dimension(:,:), intent(inout) :: A

  integer(ik),dimension(size(A,1))           :: indx
  real(dp),   dimension(size(A,1))           :: b
  real(dp),   dimension(size(A,1),size(A,1)) :: Ainv
  integer(ik) ::i,n
  real(dp)    :: d

  n=size(A,1)
  Ainv=0.0_dp
  call LU_roz (A,indx,d)
  do i=1,n
    b=0.0_dp
    b(i)=1.0_dp
    call LU_pod(A,indx,b)
    Ainv(:,i)=b
  enddo
  A=Ainv
  end subroutine LU_inv
  

  subroutine LU_roz(A,indx,d)
  real(dp),   dimension(:,:), intent(inout) :: A
  integer(ik),dimension(:),   intent(out)   :: indx
  real(dp),                   intent(out)   :: d

  real(dp), dimension(size(A,1)) ::vv,tmp
  real(dp), parameter            ::tiny=10.0_dp**(-20)
  integer(ik)                    :: j,n,imax

  n=size(A,1)
  d=1.0_dp
  vv=maxval(abs(A),dim=2)
  vv=1.0_dp/vv

  do j=1,n
    imax=(j-1)+sum(maxloc(vv(j:n)*abs(A(j:n,j))))
!   imax=j
    if (j/= imax) then
      tmp(:)=a(imax,:)
      a(imax,:)=a(j,:)
      a(j,:)=tmp(:)
      d=-d
      vv(imax)=vv(j)
    endif

    indx(j)=imax
    if (A(j,j)==0.0_dp) A(j,j)=Tiny
    A(j+1:n,j)=A(j+1:n,j)/A(j,j)
    A(j+1:n,j+1:n)=A(j+1:n,j+1:n)-spread(A(j+1:n,j),dim=2,ncopies=n-j)*&
                                  spread(A(j,j+1:n),dim=1,ncopies=n-j)
  enddo
  end subroutine LU_roz


  subroutine LU_pod (A,indx,b)
  real(dp),   dimension(:,:), intent(in)    ::A
  integer(ik),dimension(:),   intent(in)    ::indx
  real(dp),   dimension(:),   intent(inout) ::b

  integer(ik) :: i,n,ii,ll
  real(dp)    :: summ

  n=size(A,1)
  ii=0
  do i=1,n
    ll=indx(i)
    summ=b(ll)
    b(ll)=b(i)
    if (ii/=0) then
      summ=summ-dot_product(A(i,ii:i-1),b(ii:i-1))
    elseif (summ /=0.0_dp) then
      ii=i
    endif
    b(i)=summ
  enddo
  do i=n,1,-1
   b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
  enddo
  end subroutine LU_pod

end module varModule
