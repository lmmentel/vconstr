module gridAOvaluesModule
  use varModule
  use basisModule
  use commonsModule
  implicit none

  private
 
  interface AOvalueAtPoint
    module procedure getAOvaluesAtPoint
  end interface

  public :: AOvalueAtPoint


contains

  subroutine getAOvaluesAtPoint(system, basis, x, y, z, AOvalues, gradX, gradY, gradZ, gradS, nuclearPotential)
!===============================================================================
! subroutine to calculate values, gradients and Laplacian of all AO's at grid 
! point (x, y, z)
!
! input  :
!   system           : system information read from gamess-us $JOB.basinfo file
!                      see basisModule for details
!   basis            : basis information read from gamess-us $JOB.basinfo file
!                      see basisModule for details
!   x, y, z          : real, coordianted of the point
!
! output :
!   AOvalues         :
!   gradX            :
!   gradY            :
!   gradZ            :
!   gradS            :
!   nuclearPotential :
!===============================================================================
    type(systemType), intent(in) :: system
    type(basisType),  intent(in) :: basis    

    real(dp), intent(in)  :: x, y, z
    real(dp), intent(out) :: nuclearPotential
    real(dp), intent(out) :: AOvalues(:), gradX(:), gradY(:), gradZ(:), gradS(:)

    integer(ik) :: i

    AOvalues = 0.0_dp
    gradX    = 0.0_dp
    gradY    = 0.0_dp
    gradZ    = 0.0_dp
    gradS    = 0.0_dp

    nuclearPotential = 0.0_dp

! calculate nuclear potential

!    nuclearPotential = getNuclearPotential(x, y, z, basis%coords, basis%znuc)

  end subroutine getAOvaluesAtPoint

  real(dp) function getNuclearPotential(x, y, z, coords, znuc)
!==============================================================================
! calculate nuclear potential at point (x,y,z)
! 
! input  :
!   coords : array storing atomic coordinates, the format is 
!            coords(1, natom) : x coordinate of atom 'natom'
!            coords(2, natom) : y coordinate of atom 'natom'
!            coords(3, natom) : z coordinate of atom 'natom'
!   znuc   : nuclear charges, znuc(natom) is the nuclear charge of 'natom'
!==============================================================================
    real(dp), intent(in) :: x, y, z
    real(dp), intent(in) :: coords(:,:), znuc(:)

    integer(ik) :: i
    real(dp)    :: x2, y2, z2, r, nuclearPotential

    if (size(coords, 2) /= size(znuc)) then
        write(*,'(/"in -getNuclearPotential- different number of atoms: ")') &
                size(coords, 2), size(znuc) 
        stop 'exiting program...'
    else
        nuclearPotential = 0.0_dp
        do i = 1, size(coords, 2)
            x2 = (x -coords(1, i))**2  
            y2 = (y -coords(2, i))**2  
            z2 = (z -coords(3, i))**2
            r  = dsqrt(x2 + y2 + z2)
            if (r > tiny(1.0e0)) then
                nuclearPOtential = nuclearPotential + znuc(i)/r
            else
                nuclearPotential = nuclearPotential + huge(1.0e0)
            endif
        enddo
        getNuclearPotential = nuclearPotential
    endif
  end function getNuclearPotential
          
end module gridAOvaluesMOdule
