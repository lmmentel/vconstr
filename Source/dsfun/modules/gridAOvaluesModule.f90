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

  subroutine getAOvaluesAtPoint(basis, xp1, yp1, zp1, AOvalue, gradX, gradY, gradZ, gradS, nuclearPotential)
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
!   AOvalues         : value of the AO at the point (x, y, z)
!   gradX            : x gradient component of AO at the point (x, y, z) 
!   gradY            : y gradient component of AO at the point (x, y, z)
!   gradZ            : x gradient component of AO at the point (x, y, z)
!   gradS            : laplacian of the AO at the point (x, y, z)
!   nuclearPotential : nucelearPotential at th point (x, y, z)
!===============================================================================
    type(basisType),  intent(in) :: basis    

    real(dp), intent(in)  :: xp1, yp1, zp1
    real(dp), intent(out) :: nuclearPotential
    real(dp), intent(out) :: AOvalue(:), gradX(:), gradY(:), gradZ(:), gradS(:)

    integer(ik) :: i, j, k, ii, ioshell, AOindex
    real(dp)    :: x, y, z, r2, norm, gaus, exponent, zeta, drg, dxg, dyg, dzg
    real(dp)    :: mtx, mty, mtz, xp, yp, zp

    AOvalue = 0.0_dp
    gradX   = 0.0_dp
    gradY   = 0.0_dp
    gradZ   = 0.0_dp
    gradS   = 0.0_dp
! just as a test set point to 0.0, 0.0, 0.0
    xp = 0.0_dp
    yp = 0.0_dp
    zp = 0.0_dp


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
            
            zeta = basis%exponent(ii)
            exponent = zeta*r2
            if (exponent < explim) then
                norm = normalizeGaussian(basis%intyp(i), k, zeta) 
                gaus = basis%coefficient(ii)*norm*exp(-exponent)
                drg  =  2.0_dp*zeta*r2
                dxg  = -2.0_dp*zeta*x 
                dyg  = -2.0_dp*zeta*y 
                dzg  = -2.0_dp*zeta*z 

        orbitalType : select case(basis%intyp(i))
            case(1)       ! S orbitals
                AOvalue(AOindex) = AOvalue(AOindex) + gaus
                gradx(AOindex)   = gradx(AOindex)   + gaus*dxg 
                grady(AOindex)   = grady(AOindex)   + gaus*dyg 
                gradz(AOindex)   = gradz(AOindex)   + gaus*dzg 
                grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*zeta*(drg-3.0_dp) 
            case(2)       ! P orbitals
                select case (k)     ! P 
                    case (2)    ! P X
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*x
                        gradx(AOindex)   = gradx(AOindex)   + gaus*dxg*x + gaus 
                        grady(AOindex)   = grady(AOindex)   + gaus*dyg*x 
                        gradz(AOindex)   = gradz(AOindex)   + gaus*dzg*x 
                        grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*x*zeta*(drg-5.0_dp) 
                    case (3)    ! P Y
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*y
                        gradx(AOindex)   = gradx(AOindex)   + gaus+dxg*y 
                        grady(AOindex)   = grady(AOindex)   + gaus*dyg*y + gaus 
                        gradz(AOindex)   = gradz(AOindex)   + gaus*dzg*y 
                        grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*y*zeta*(drg-5.0_dp) 
                    case (4)    ! P Z
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*z
                        gradx(AOindex)   = gradx(AOindex)   + gaus+dxg*z 
                        grady(AOindex)   = grady(AOindex)   + gaus*dyg*z + gaus 
                        gradz(AOindex)   = gradz(AOindex)   + gaus*dzg*z 
                        grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*z*zeta*(drg-5.0_dp) 
                    case default 
                        write(*,*) 'wrong subshell for P shell in -grid-: ', k
                end select
            case(3)       ! D orbitals
                mtx = -2.0_dp*x 
                mty = -2.0_dp*y 
                mtz = -2.0_dp*z 
                select case (k)
                    case ( 5)    ! D XX
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*x*x
                        gradx(AOindex)   = gradx(AOindex)   + gaus*mtx*(zeta*x*x-1.0_dp) 
                        grady(AOindex)   = grady(AOindex)   + gaus*dyg*x*x 
                        gradz(AOindex)   = gradz(AOindex)   + gaus*dzg*x*x 
                        grads(AOindex)   = grads(AOindex)   + gaus*(2.0_dp+2.0_dp*x*x*zeta*(drg-7.0_dp)) 
                    case ( 6)    ! D YY
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*y*y
                        gradx(AOindex)   = gradx(AOindex)   + gaus*dxg*y*y
                        grady(AOindex)   = grady(AOindex)   + gaus*mty*(zeta*y*y-1.0_dp) 
                        gradz(AOindex)   = gradz(AOindex)   + gaus*dzg*y*y 
                        grads(AOindex)   = grads(AOindex)   + gaus*(2.0_dp+2.0_dp*y*y*zeta*(drg-7.0_dp)) 
                    case ( 7)    ! D ZZ
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*z*z
                        gradx(AOindex)   = gradx(AOindex)   + gaus*dxg*z*z
                        grady(AOindex)   = grady(AOindex)   + gaus*dyg*z*z 
                        gradz(AOindex)   = gradz(AOindex)   + gaus*mtz*(zeta*z*z-1.0_dp) 
                        grads(AOindex)   = grads(AOindex)   + gaus*(2.0_dp+2.0_dp*z*z*zeta*(drg-7.0_dp)) 
                    case ( 8)    ! D XY
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*x*y
                        gradx(AOindex)   = gradx(AOindex)   + gaus*(y-dxg*x*y) 
                        grady(AOindex)   = grady(AOindex)   + gaus*(x-dyg*x*y)
                        gradz(AOindex)   = gradz(AOindex)   + gaus*dzg*x*y 
                        grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*x*y*zeta*(drg-7.0_dp) 
                    case ( 9)    ! D XZ
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*x*z
                        gradx(AOindex)   = gradx(AOindex)   + gaus*(z-dxg*x*z) 
                        grady(AOindex)   = grady(AOindex)   + gaus*dyg*x*z
                        gradz(AOindex)   = gradz(AOindex)   + gaus*(x-dzg*x*z) 
                        grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*x*z*zeta*(drg-7.0_dp) 
                    case (10)    ! D YZ
                        AOvalue(AOindex) = AOvalue(AOindex) + gaus*y*z
                        gradx(AOindex)   = gradx(AOindex)   + gaus*dxg*y*z 
                        grady(AOindex)   = grady(AOindex)   + gaus*(z-dyg*y*z)
                        gradz(AOindex)   = gradz(AOindex)   + gaus*(y-dzg*y*z) 
                        grads(AOindex)   = grads(AOindex)   + gaus*2.0_dp*y*z*zeta*(drg-7.0_dp) 
                    case default
                        write(*,*) 'wrong subshell for D shell in -grid-: ', k
                end select
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

  call printAOvaluesAtPoint(basis, xp, yp, zp, AOvalue, gradX, gradY, gradZ, gradS, nuclearPotential)
  stop 'in grid'
  end subroutine getAOvaluesAtPoint

  subroutine printAOvaluesAtPoint(basis, x, y, z, AOvalue, gradX, gradY, gradZ, gradS, nuclearPotential)
    type(basisType),  intent(in) :: basis    

    real(dp), intent(in)  :: x, y, z
    real(dp), intent(in) :: nuclearPotential
    real(dp), intent(in) :: AOvalue(:), gradX(:), gradY(:), gradZ(:), gradS(:)

    integer(ik) :: i,  k, styp, ioshell, AOi, atom
    character(len=4) :: alab
    character(len=2) :: otyp


    call print_header('AO values at point (x, y, z)')
    write(*,'(/"Coordinates of the point (x, y, z)",3f15.10/)') x, y, z
    write(*,'("Nuclear Potential at (",3f15.10,") = ",es24.10/)') x, y, z, nuclearPotential

    write(*,'("AO ",2x,"atom",2x,"cont",2x,"shell",3x,"shell",2x,"subshell")')
    write(*,'("No.",2x,"name",2x,"gaus",2x,"No.  ",3x,"type ",2x,"  type ",5a24)') & 
        "AO value", "Gradient in X", "Gradient in Y", "Gradient in Z", "Laplacian"
    write(*,'(a)') repeat('-', 159)

    do i = 1, basis%nshell
        atom = basis%katom(i)
        ioshell = 0
        do k = kmin(basis%intyp(i)), kmax(basis%intyp(i))
            ioshell = ioshell + 1
            AOi = basis%AOlocation(i) + ioshell - 1
            alab = atomlabel(int(basis%znuc(basis%katom(i))))
            styp = basis%intyp(i)
            otyp = label(basis%intyp(i))
            write(*,'(i3,2x,a4,2x,i4,3x,i4,5x,a2,3x,a4,3x,5f24.10)') &
                AOi, alab, basis%kng(i), i, otyp, bfnam1(k), AOvalue(AOi), gradX(AOi), gradY(AOi), gradZ(AOi), gradS(AOi)
        enddo
    enddo
  end subroutine printAOvaluesAtPoint



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
