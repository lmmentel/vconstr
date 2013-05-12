module matrixModule
use varModule
implicit none

contains 

  subroutine scaler(n, scalar, r, b)
    real(DP), intent(in)  :: b(:)
    real(DP), intent(out) :: r(:)
    real(DP), intent(in)  :: scalar
    integer,  intent(in)  :: n

    r = scalar*b
  end subroutine scaler

  subroutine fmove(a, b, n)
    real(DP), intent(in)  :: a(:)
    real(DP), intent(out) :: b(:)
    integer,  intent(in)  :: n

    b = a
  end subroutine fmove

  subroutine gtriad(n, scalar, r, a, b)
    real(DP), intent(in)  :: a(:), b(:)
    real(DP), intent(out) :: r(:)
    real(DP), intent(in)  :: scalar
    integer,  intent(in)  :: n

    r = a + scalar*b
  end subroutine gtriad

end module matrixModule
