module commonsModule
  use varModule
!===============================================================================
! a module for storing some common variables used throughout the program at 
! various places
!===============================================================================
  implicit none

  character(len=100) :: basisInfoFile
  character(len=100) :: dictionaryFile
  character(len=100) :: integralsFile

  integer :: printLevel

! maximal exponent handled in the Grid calcalation in grid module
  real(dp), parameter :: explim = 50000.0_dp
! epsilon 
  real(dp), parameter :: epsilon = 1.0e-6_dp

end module commonsModule
