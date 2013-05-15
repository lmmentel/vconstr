module basisModule
  use varModule
  use commonsModule
  implicit none

  character(len=4), dimension(106) :: atomLabel = (/'H   ','HE  ','LI  ','BE  ','B   ','C   ',     &
                                                    'N   ','O   ','F   ','NE  ','NA  ','MG  ',     &
                                                    'AL  ','SI  ','P   ','S   ','CL  ','AR  ',     &
                                                    'K   ','CA  ','SC  ','TI  ','V   ','CR  ',     &
                                                    'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',     &
                                                    'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',     &
                                                    'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ',     &
                                                    'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',     &
                                                    'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',     &
                                                    'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',     &
                                                    'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',     &
                                                    'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',     &
                                                    'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',     &
                                                    'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',     &
                                                    'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',     &
                                                    'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',     &
                                                    'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',     &
                                                    'LR  ','RF  ','X   ','BQ  '/)

  character(len=2), dimension(8)  :: LABEL = (/'S ','P ','D ','F ','G ','H ','I ','L '/)

  character(len=4), dimension(35) :: BFNAM1 = (/'  S ','  X ','  Y ','  Z ',                       &
                                                ' XX ',' YY ',' ZZ ',' XY ',' XZ ',' YZ ',         &
                                                ' XXX',' YYY',' ZZZ',' XXY',' XXZ',                &
                                                ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',                &
                                                'XXXX','YYYY','ZZZZ','XXXY','XXXZ',                &
                                                'YYYX','YYYZ','ZZZX','ZZZY','XXYY',                &
                                                'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/)

  character(len=6), dimension(49) :: BFNAM2 = (/' XXXXX',' YYYYY',' ZZZZZ',' XXXXY',' XXXXZ',      &
                                                ' YYYYX',' YYYYZ',' ZZZZX',' ZZZZY',' XXXYY',      &
                                                ' XXXZZ',' YYYXX',' YYYZZ',' ZZZXX',' ZZZYY',      &
                                                ' XXXYZ',' YYYXZ',' ZZZXY',' XXYYZ',' XXZZY',      &
                                                ' YYZZX',                                          & 
                                                '    X6','    Y6','    Z6','   X5Y','   X5Z',      &
                                                '   Y5X','   Y5Z','   Z5X','   Z5Y','  X4Y2',      &
                                                '  X4Z2','  Y4X2','  Y4Z2','  Z4X2','  Z4Y2',      &
                                                '  X4YZ','  Y4XZ','  Z4XY','  X3Y3','  X3Z3',      &
                                                '  Y3Z3',' X3Y2Z',' X3Z2Y',' Y3X2Z',' Y3Z2X',      &
                                                ' Z3X2Y',' Z3Y2X','X2Y2Z2'/)

  integer, dimension(8) :: KMIN = (/1,2, 5,11,21,34,57,1/)
  integer, dimension(8) :: KMAX = (/1,4,10,20,35,56,84,4/)


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
    integer, allocatable  :: imin(:)
    integer, allocatable  :: imax(:)
    integer, allocatable  :: katom(:)
    integer, allocatable  :: intyp(:)
    integer, allocatable  :: ityp(:) 
    integer, allocatable  :: kng(:) 
    real(dp), allocatable :: exponent(:)
    real(dp), allocatable :: coefficient(:)
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

    call read_system_info(trim(basisInfoFile), self%title, self%natoms, self%charge, self%mult,    &
                       self%nbf, self%nx, self%ne, self%na, self%nb, self%nshell, self%nprimi)
  end subroutine newSystem

  subroutine newBasis(self, sys)
    type(basisType),  intent(out) :: self
    type(systemType), intent(out) :: sys

    integer :: i, j, AOindex, nprimitives

    integer, allocatable :: ish(:)
    real(dp), allocatable :: expon(:), contrc1(:), contrc2(:)

! first get information about the system
    call newSystem(sys)

! allocate first batch of arrays 
    
    allocate(self%znuc(sys%natoms))
    allocate(self%coords(3, sys%natoms))
    allocate(self%imin(sys%natoms))
    allocate(self%imax(sys%natoms))
    allocate(self%evec(3))
    allocate(self%katom(sys%nshell))
    allocate(self%intyp(sys%nshell))
    allocate(self%ityp(sys%nprimi))
    allocate(self%kng(sys%nshell))
! allocate temporary arrays 
    allocate(ish(sys%nprimi))
    allocate(expon(sys%nprimi))
    allocate(contrc1(sys%nprimi))
    allocate(contrc2(sys%nprimi))

    call read_basis_info(trim(basisInfoFile), sys%natoms, sys%nshell, sys%nprimi,                  &
                         self%znuc, self%coords, self%evec, expon, contrc1,              &
                         contrc2, self%imin, self%imax, self%katom, & 
                         self%intyp, ish, self%ityp, self%kng)
! calculate the total number of primitives (non-unique)
    nprimitives = 0
    do i = 1, sys%nshell
            nPrimitives = nPrimitives + self%kng(i)
    enddo
    
    allocate(self%exponent(nPrimitives), self%coefficient(nPrimitives))
    AOindex = 0
    do i = 1, sys%natoms
        do j = self%imin(i), self%imax(i)
            AOindex = AOindex + 1 
            self%exponent(AOindex)    = expon(j)
            self%coefficient(AOindex) = contrc1(j)
        enddo 
    enddo

    if (printLevel > 1) call print_basis_info(self, sys)
! deallocate temporary arrays
    deallocate(ish, expon, contrc1, contrc2) 
  end subroutine newBasis

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
                           expon, contrc1, contrc2, imin, imax, katom, intyp, ish, ityp, kng)
 implicit none 
  character(len=*),  intent(in)  :: filename
  integer, intent(in)  :: natoms, nshell, nprimi
  real(dp),          intent(out) :: znuc(:), coords(:,:), evec(:)
  real(dp),          intent(out) :: expon(:), contrc1(:), contrc2(:) 
  integer, intent(out) :: imin(:), imax(:), katom(:), intyp(:), ish(:), ityp(:), kng(:)

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
      read(300,*) (kng(i),i=1,nshell)
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

    integer :: i, j, ij, k, atom, styp, jmin, jmax
    character(len=4) :: alab
    character(len=2) :: otyp

    call print_header('Geometry in a.u.')
    write(*,'(5x,a4,3x,a10,5x,3a15)') 'Atom','Charge  ', 'X     ', 'Y     ', 'Z     ' 
    write(*,'(5x,a4,3x,a10,5x,a45)') repeat('-',4), repeat('-', 10), repeat('-', 45) 
    do i = 1, system%natoms
        write(*,'(5x,a4,3x,f10.2,5x,3f15.8)') atomLabel(int(self%znuc(i))), self%znuc(i),         &
                                        self%coords(1, i), self%coords(2, i), self%coords(3, i)
    enddo
    call print_header('Basis set')
    write(*,'(2x,"Shell",3x,"Type",3x,"Primitive",7x,"Exponent",5x,"Contraction Coefficient(s)")')
    ij = 0
    do i = 1, system%natoms
        write(*,'(/a4,a10,i5,5x,a9,f5.2/)') atomLabel(int(self%znuc(i))), 'Atom No.', i, 'Charge = ', self%znuc(i)
        do j = self%imin(i), self%imax(i)
            ij = ij + 1
            if (self%ityp(j) < 8) then
                write(*,'(1x,i6,3x,a4,i7,f22.7,f18.12,3x,i3)') & 
                    i, label(self%ityp(j)), j, self%exponent(j), self%coefficient(j), ij
            else
                write(*,'(1x,i6,3x,a4,i7,f22.7,2f18.12)') & 
                    i, label(self%ityp(j)), j, self%exponent(j), self%coefficient(j)
            endif
        enddo
    enddo

    write(*,'("atom",3x,"No. ",2x,"shell",3x,"type")') 
    write(*,'(a)') repeat('-', 30)
    ij = 0
    do i = 1, system%nshell
                atom = self%katom(i)
        do j = 1, self%kng(i)
               ij  = ij + 1 
            do k = kmin(self%intyp(i)), kmax(self%intyp(i))
                alab = atomlabel(int(self%znuc(self%katom(i))))
                styp = self%intyp(i)
                otyp = label(self%intyp(i)) 
                write(*,'(a4,3x,i4,3x,i4,5x,a2,3x,a4,3x,2f24.10)') &
                    alab, atom, i, otyp, bfnam1(k), self%exponent(ij), self%coefficient(ij)
            enddo
        enddo
    enddo

  end subroutine print_basis_info

  subroutine deleteBasis(self)
    type(basisType) :: self

    if (allocated(self%ityp))        deallocate(self%ityp) 
    if (allocated(self%exponent))    deallocate(self%exponent) 
    if (allocated(self%coefficient)) deallocate(self%coefficient)
    if (allocated(self%katom))       deallocate(self%katom)
    if (allocated(self%intyp))       deallocate(self%intyp)
    if (allocated(self%znuc))        deallocate(self%znuc)
    if (allocated(self%coords))      deallocate(self%coords)
    if (allocated(self%imin))        deallocate(self%imin)
    if (allocated(self%imax))        deallocate(self%imax)
    if (allocated(self%evec))        deallocate(self%evec)
  end subroutine deleteBasis

end module basisModule
