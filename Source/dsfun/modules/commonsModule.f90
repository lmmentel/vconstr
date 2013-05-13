module commonsModule
!===============================================================================
! a module for storing some common variables used throughout the program at 
! various places
!===============================================================================
  implicit none

  character(len=100) :: basisInfoFile
  character(len=100) :: dictionaryFile
  character(len=100) :: integralsFile

  integer :: print_level

end module commonsModule
