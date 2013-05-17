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

  subroutine getAOvaluesAtPoint(basis, xp, yp, zp, AOvalues, gradX, gradY, gradZ, gradS, nuclearPotential)
!===============================================================================
! subroutine to calculate values, gradients and Laplacian of all AO's at grid 
! point (x, y, z)
!
! input  :
!   basis            : basis set and system information read from gamess-us 
!                      $JOB.basinfo file see basisModule for details
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
    type(basisType),  intent(in) :: basis    

    real(dp), intent(in)  :: xp, yp, zp
    real(dp), intent(out) :: nuclearPotential
    real(dp), intent(out) :: AOvalues(:), gradX(:), gradY(:), gradZ(:), gradS(:)

    integer(ik) :: i, j, k, ii, ioshell, AOindex
    real(dp)    :: x, y, z, r2, norm, gaus, exponent

    AOvalues = 0.0_dp
    gradX    = 0.0_dp
    gradY    = 0.0_dp
    gradZ    = 0.0_dp
    gradS    = 0.0_dp

    nuclearPotential = 0.0_dp

! calculate nuclear potential

    nuclearPotential = getNuclearPotential(x, y, z, basis%coords, basis%znuc)

! calculate value gradients and laplacians

    ii = 0
    do i = 1, basis%nshell

! tranform coordinate (xp, yp, zp) from main coordinate frame to coordinate (x, y, z)
! in a frame with the atomic orbital center as an origin

    x = xp - basis%coords(1, basis%katom(i))
    y = yp - basis%coords(2, basis%katom(i))
    z = zp - basis%coords(3, basis%katom(i))
    r2 = x**2 + y**2 + z**2 

    do j = 1, basis%kng(i)
        ioshell = 0
        do k = kmin(basis%intyp(i)), kmax(basis%intyp(i))
            ioshell = ioshell + 1
            AOindex = basis%AOlocation(i) + ioshell - 1
            ii = ii + 1
            !write(*,*) ii, AOindex, i, j, k
            
            exponent = basis%exponent(ii)*r2
            if (exponent < explim) then
                norm = normalizeGaussian(basis%intyp(i), k, basis%exponent(ii)) 
                gaus = basis%coefficient(ii)*norm*exp(-exponent) 

        orbitalType : select case(basis%intyp(i))
            case(1)       ! S orbitals
                write(*,'(5i4,3x,a7)') ii, AOindex, i, j, k, 'S shell'
!                !valao(AOindex) = valao(AOindex) + 
!            case(2)       ! P orbitals
!                write(*,'(2i4,3x,"Orbital: ",i4,3x,a7)') ij, k, AOindex, 'P shell'
!            case(3)       ! D orbitals
!                write(*,'(2i4,3x,"Orbital: ",i4,3x,a7)') ij, k, AOindex, 'D shell'
!            case(4)       ! F orbitals
!                write(*,'(2i4,3x,"Orbital: ",i4,3x,a7)') ij, k, AOindex, 'F shell'
!            case(5)       ! G orbitals
!                write(*,'(2i4,3x,"Orbital: ",i4,3x,a7)') ij, k, AOindex, 'G shell'
            case default  ! higher angular momenta are not handled by gamess-us
                write(*,*) 'maximum angular momentum for a basis function is 4 (g functions)'
        end select orbitalType   
            endif
        enddo 
    enddo
  enddo
  stop 'in grid'
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
