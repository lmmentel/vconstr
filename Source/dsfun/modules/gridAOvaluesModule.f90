module gridAOvaluesModule
  use varModule
  use basisModule
  use commonsModule
  implicit none

  private

contains

  subroutine getAOvaluesAtPoint(x, y, z, AOvalues, gradX, gradY, gradZ, gradS, nuclearPotential)
!===============================================================================
! subroutine to calculate values, gradients and Laplacian of all AO's at grid 
! point (x, y, z)
!
! input  :
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

    real(dp), intent(in)  :: x, y, z
    real(dp), intent(out) :: nuclearPotential
    real(dp), intent(out) :: AOvalues(:), gradX(:), gradY(:), gradZ(:), gradS(:)

! variables holding basis set information form gamess-us $JOB.basinfo file
   
    character(len=80) :: title
    integer :: natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi

    real(dp),          allocatable :: znuc(:), coords(:,:), evec(:)
    real(dp),          allocatable :: expon(:), contrc1(:), contrc2(:) 
    integer, allocatable :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:)

! local variables 

    integer(ik) :: i

    write(*,*) 'entering -getAOvaluesAtPoint-', print_level


! get scalar variables from gamess-us $JOB.basinfo file

    call read_job_info(trim(basisInfoFile), title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi) 

! allocate the arrays storing system and absis set information

    allocate(znuc(natoms), coords(3, natoms), imin(natoms), imax(natoms), evec(3))
    allocate(katom(nshell), intyp(nshell)) 
    allocate(ish(nprimi), ityp(nprimi), expon(nprimi), contrc1(nprimi), contrc2(nprimi)) 

! read the remaining basis set information

    call read_basis_info(trim(basisInfoFile), natoms, nshell, nprimi, znuc, coords, evec,          &
                         expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)

! print the data read from $JOB.basinfo file
!    if (print_level > 1) then 
        call write_basis_info(title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi,    &
             znuc, coords, evec, expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
!    endif 

    AOvalues = 0.0_dp
    gradX    = 0.0_dp
    gradY    = 0.0_dp
    gradZ    = 0.0_dp
    gradS    = 0.0_dp

    nuclearPotential = 0.0_dp

! calculate nuclear potential

    nuclearPotential = getNuclearPotential(x, y, z, coords, znuc)

  deallocate(ish, ityp, expon, contrc1, contrc2)
  deallocate(katom, intyp)
  deallocate(znuc, coords, imin, imax, evec)

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
