module basisModule

contains

 subroutine print_header(header)
  character(len=*), intent(in) :: header

  write(*,'(/5x,a)') repeat('=',len(header))
  write(*,'(5x,a)') header
  write(*,'(5x,a/)') repeat('=',len(header))
 end subroutine print_header

subroutine read_job_info(filename, title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi)
 implicit none 
  character(len=*),  intent(in)  :: filename
  character(len=80), intent(out) :: title
  integer,           intent(out) :: natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi

  open(unit=300, file=filename, form='formatted', status='old')
  read(300,*) title
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
end subroutine read_job_info

subroutine read_basis_info(filename, natoms, nshell, nprimi, znuc, coords, evec,                   &
                           expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
 implicit none 
  character(len=*),  intent(in)  :: filename
  integer,           intent(in)  :: natoms, nshell, nprimi
  real(kind(1.0d0)), intent(out) :: znuc(:), coords(:,:), evec(:)
  real(kind(1.0d0)), intent(out) :: expon(:), contrc1(:), contrc2(:) 
  integer,           intent(out) :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:)

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

subroutine write_basis_info(title, natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi,      &
       znuc, coords, evec, expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp)
 implicit none 
  character(len=80), intent(in) :: title
  integer           :: natoms, charge, mult, nbf, nx, ne, na, nb, nshell, nprimi
  real(kind(1.0d0)), intent(in) :: znuc(:), coords(:,:), evec(:)
  real(kind(1.0d0)), intent(in) :: expon(:), contrc1(:), contrc2(:) 
  integer,           intent(in) :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:) 

  character(len=1), dimension(8) :: label=(/'S','P','D','F','G','H','I','L'/)
  integer :: i, j

  call print_header('Gamess-US job parameters')
  write(*,'(a47,a)') 'Gamess-US job title                          = ', trim(title)
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
