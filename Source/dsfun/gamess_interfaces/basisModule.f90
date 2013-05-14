module basisModule
  use varModule
  use commonsModule
  implicit none

  type systemType
    character(len=80) :: title
    integer :: natoms
    integer :: charge
    integer :: mult
    integer :: nbf
    integer :: nx
    integer :: ne
    integer :: na 
    integer :: nb
    integer :: nshell
    integer :: nprimi
  end type systemType

  type basisType
    real(dp), allocatable :: znuc(:) 
    real(dp), allocatable :: coords(:,:)
    real(dp), allocatable :: evec(:)
    real(dp), allocatable :: expon(:)
    real(dp), allocatable :: contrc1(:)
    real(dp), allocatable :: contrc2(:) 
    integer, allocatable  :: imin(:)
    integer, allocatable  :: imax(:)
    integer, allocatable  :: katom(:)
    integer, allocatable  :: intyp(:)
    integer, allocatable  :: ish(:)
    integer, allocatable  :: ityp(:) 
  end type basisType


contains

 subroutine print_header(header)
  character(len=*), intent(in) :: header

  write(*,'(/5x,a)') repeat('=',len(header))
  write(*,'(5x,a)') header
  write(*,'(5x,a/)') repeat('=',len(header))
 end subroutine print_header

  subroutine newSystem(self)
    use commonsModule
    type(systemType), intent(out) :: self

    call read_system_info(trim(basisInfoFile), self%title, self%natoms, self%charge, self%mult,       &
                       self%nbf, self%nx, self%ne, self%na, self%nb, self%nshell, self%nprimi)
! if debugging on print the information
    if (printLevel > 1) call print_system_info(self)
  end subroutine newSystem

  subroutine newBasis(self, sys)
    type(basisType),  intent(out) :: self
    type(systemType), intent(out) :: sys

! first get information about the system
    call newSystem(sys)

! allocate the arrays
    
    allocate(self%znuc(sys%natoms))
    allocate(self%coords(3, sys%natoms))
    allocate(self%imin(sys%natoms))
    allocate(self%imax(sys%natoms))
    allocate(self%evec(3))
    allocate(self%katom(sys%nshell))
    allocate(self%intyp(sys%nshell))
    allocate(self%ish(sys%nprimi))
    allocate(self%ityp(sys%nprimi))
    allocate(self%expon(sys%nprimi))
    allocate(self%contrc1(sys%nprimi))
    allocate(self%contrc2(sys%nprimi))

    call read_basis_info(trim(basisInfoFile), sys%natoms, sys%nshell, sys%nprimi,                  &
                         self%znuc, self%coords, self%evec, self%expon, self%contrc1,         &
                         self%contrc2, self%imin, self%imax, self%katom, self%intyp,          &
                         self%ish, self%ityp)

    if (printLevel > 1) call print_basis_info(self, sys) 
  end subroutine newBasis

  subroutine deleteBasis(self)
    type(basisType) :: self

    
    if (allocated(self%ish))     deallocate(self%ish)
    if (allocated(self%ityp))    deallocate(self%ityp) 
    if (allocated(self%expon))   deallocate(self%expon) 
    if (allocated(self%contrc1)) deallocate(self%contrc1)
    if (allocated(self%contrc2)) deallocate(self%contrc2)
    if (allocated(self%katom))   deallocate(self%katom)
    if (allocated(self%intyp))   deallocate(self%intyp)
    if (allocated(self%znuc))    deallocate(self%znuc)
    if (allocated(self%coords))  deallocate(self%coords)
    if (allocated(self%imin))    deallocate(self%imin)
    if (allocated(self%imax))    deallocate(self%imax)
    if (allocated(self%evec))    deallocate(self%evec)
  end subroutine deleteBasis

subroutine read_system_info(filename, title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi)
 implicit none 
  character(len=*),  intent(in)  :: filename
  character(len=80), intent(out) :: title
  integer, intent(out) :: natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi

  open(unit=300, file=filename, form='formatted', status='old')
  read(300,'(a80)') title
  read(300,*) natoms
  read(300,*) charge
  read(300,*) mult
  read(300,*) nbf
  read(300,*) nx
  read(300,*) ne
  read(300,*) na
  read(300,*) nb
  read(300,*) nshell
  read(300,*) nprimi
  close(300)
  return 
end subroutine read_system_info

subroutine read_basis_info(filename, natoms, nshell, nprimi, znuc, coords, evec,                   &
                           expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
 implicit none 
  character(len=*),  intent(in)  :: filename
  integer, intent(in)  :: natoms, nshell, nprimi
  real(dp),          intent(out) :: znuc(:), coords(:,:), evec(:)
  real(dp),          intent(out) :: expon(:), contrc1(:), contrc2(:) 
  integer, intent(out) :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:)

  integer :: i, j 

  open(unit=300, file=filename, form='formatted', status='old')
  rewind(300)
  do i = 1, 11
    read(300,*)
  enddo
      read(300,*) (znuc(i),i=1,natoms)
      read(300,*) ((coords(j,i),j=1,3),i=1,natoms)
      read(300,*) (imin(i),i=1,natoms)
      read(300,*) (imax(i),i=1,natoms)
      read(300,*) (evec(i),i=1,3)
      read(300,*) (katom(i),i=1,nshell)
      read(300,*) (intyp(i),i=1,nshell)
      read(300,*) (ish(i),i=1,nprimi)
      read(300,*) (ityp(i),i=1,nprimi)
      read(300,*) (expon(i),i=1,nprimi)
      read(300,*) (contrc1(i),i=1,nprimi) 
      read(300,*) (contrc2(i),i=1,nprimi)
  close(300)
  return 
end subroutine read_basis_info

  subroutine print_system_info(self)
    type(systemType), intent(in) :: self

    call print_header('Gamess-US job parameters')
    write(*,'(a47,a80)') 'Gamess-US job title                          = ', self%title
    write(*,'(a47,i5)') 'Number of atoms                              = ', self%natoms
    write(*,'(a47,i5)') 'Total charge                                 = ', self%charge
    write(*,'(a47,i5)') 'Multiplicity                                 = ', self%mult
    write(*,'(a47,i5)') 'Number of basis functions                    = ', self%nbf
    write(*,'(a47,i5)') 'Number of orbitals (spherical)               = ', self%nx
    write(*,'(a47,i5)') 'Total number of electrons                    = ', self%ne
    write(*,'(a47,i5)') 'Number of alpha electrons                    = ', self%na
    write(*,'(a47,i5)') 'Number of beta electrons                     = ', self%nb
    write(*,'(a47,i5)') 'Number of basis set shells                   = ', self%nshell
    write(*,'(a47,i5)') 'Number of cartesian gaussian basis functions = ', self%nprimi
  end subroutine print_system_info

  subroutine print_basis_info(self, system)
    type(basisType),  intent(in) :: self
    type(systemType), intent(in) :: system

    character(len=1), dimension(8) :: label=(/'S','P','D','F','G','H','I','L'/)
    integer :: i, j
    call print_header('Geometry')
    write(*,'(5x,a10,5x,3a10)') 'Charge', 'X', 'Y', 'Z' 
    write(*,'(5x,a10,5x,a30)') repeat('-', 10), repeat('-', 30) 
    do i = 1, system%natoms
        write(*,'(5x,f10.2,5x,3f10.5)') self%znuc(i), self%coords(1, i), self%coords(2, i), self%coords(3, i)
    enddo
    call print_header('Basis set')
    write(*,'(2x,"Shell",3x,"Type",3x,"Primitive",7x,"Exponent",5x,"Contraction Coefficient(s)")')
    do i = 1, system%natoms
        write(*,'(/5x,a4,i5,5x,a9,f5.2/)') 'Atom', i, 'Charge = ', self%znuc(i)
            do j = self%imin(i), self%imax(i)
                if (self%ityp(j) < 8) then
                    write(*,'(1x,i6,3x,a4,i7,f22.7,2f18.12)') & 
                        self%ish(j), label(self%ityp(j)), j, self%expon(j), self%contrc1(j)
                else
                    write(*,'(1x,i6,3x,a4,i7,f22.7,2f18.12)') & 
                        self%ish(j), label(self%ityp(j)), j, self%expon(j), self%contrc1(j), self%contrc2(j)
                endif
            enddo
    enddo
  end subroutine print_basis_info

subroutine write_basis_info(title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi,      &
       znuc, coords, evec, expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
 implicit none 
  character(len=80), intent(in) :: title
  integer, intent(in) :: natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi
  real(dp),          intent(in) :: znuc(:), coords(:,:), evec(:)
  real(dp),          intent(in) :: expon(:), contrc1(:), contrc2(:) 
  integer, intent(in) :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:) 

  character(len=1), dimension(8) :: label=(/'S','P','D','F','G','H','I','L'/)
  integer :: i, j

  call print_header('Gamess-US job parameters')
  write(*,'(a47,a80)') 'Gamess-US job title                          = ', title
  write(*,'(a47,i5)') 'Number of atoms                              = ', natoms
  write(*,'(a47,i5)') 'Total charge                                 = ', charge
  write(*,'(a47,i5)') 'Multiplicity                                 = ', mult
  write(*,'(a47,i5)') 'Number of basis functions                    = ', nbf
  write(*,'(a47,i5)') 'Number of orbitals (spherical)               = ', nx
  write(*,'(a47,i5)') 'Total number of electrons                    = ', ne
  write(*,'(a47,i5)') 'Number of alpha electrons                    = ', na
  write(*,'(a47,i5)') 'Number of beta electrons                     = ', nb
  write(*,'(a47,i5)') 'Number of basis set shells                   = ', nshell
  write(*,'(a47,i5)') 'Number of cartesian gaussian basis functions = ', nprimi
  call print_header('Geometry')
  write(*,'(5x,a10,5x,3a10)') 'Charge', 'X', 'Y', 'Z' 
  write(*,'(5x,a10,5x,a30)') repeat('-', 10), repeat('-', 30) 
  do i = 1, natoms
    write(*,'(5x,f10.2,5x,3f10.5)') znuc(i), coords(1, i), coords(2, i), coords(3, i)
  enddo
  call print_header('Basis set')
  write(*,'(2x,"Shell",3x,"Type",3x,"Primitive",7x,"Exponent",5x,"Contraction Coefficient(s)")')
  do i = 1, natoms
    write(*,'(/5x,a4,i5,5x,a9,f5.2/)') 'Atom', i, 'Charge = ', znuc(i)
    do j = imin(i), imax(i)
      if (ityp(j) < 8) then
        write(*,'(1x,i6,3x,a4,i7,f22.7,2f18.12)') ish(j), label(ityp(j)), j, expon(j), contrc1(j)
      else
        write(*,'(1x,i6,3x,a4,i7,f22.7,2f18.12)') ish(j), label(ityp(j)), j, expon(j), contrc1(j), contrc2(j)
      endif
    enddo
  enddo
end subroutine write_basis_info

end module basisModule

!program test_basis_info
! use basisModule
! implicit none
!  character(len=100) :: filename
!  character(len=80)  :: title
!  integer            :: natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi
!  real(kind(1.0d0)), allocatable :: znuc(:), coords(:,:), evec(:)
!  real(kind(1.0d0)), allocatable :: expon(:), contrc1(:), contrc2(:) 
!  integer,           allocatable :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:) 
! 
!  filename = 'basis.info'
!
!  call read_job_info(filename, title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi)
!  allocate(znuc(natoms), coords(3,natoms), imin(natoms), imax(natoms), evec(3))
!  allocate(katom(nshell), intyp(nshell))
!  allocate(ish(nprimi), ityp(nprimi), expon(nprimi), contrc1(nprimi), contrc2(nprimi))
!
! call read_basis_info(filename, natoms, nshell, nprimi, znuc, coords, evec,                        &
!                           expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
!  call write_basis_info(title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi,          &
!       znuc, coords, evec, expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
!
!  deallocate(ish, ityp, expon, contrc1, contrc2)
!  deallocate(katom, intyp)
!  deallocate(znuc, coords, imin, imax, evec)
!end program test_basis_info
