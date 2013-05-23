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
    
  type basisType
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
    integer :: nPrimitives
    real(dp), allocatable :: znuc(:) 
    real(dp), allocatable :: coords(:,:)
    real(dp), allocatable :: evec(:)
    integer, allocatable  :: imin(:)
    integer, allocatable  :: imax(:)
    integer, allocatable  :: katom(:)
    integer, allocatable  :: intyp(:)
    integer, allocatable  :: shellType(:) 
    integer, allocatable  :: kng(:) 
    integer, allocatable  :: AOlocation(:) 
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

  subroutine newBasis(self)
    type(basisType),  intent(out) :: self

    integer :: i, j, AOindex
! local arrays
    integer,  allocatable :: ityp(:)
    real(dp), allocatable :: expon(:), contrc1(:), contrc2(:)

! first get information about the system

    call read_system_info(trim(basisInfoFile), self%title, self%natoms, self%charge, self%mult,    &
                       self%nbf, self%nx, self%ne, self%na, self%nb, self%nshell, self%nprimi)

! allocate first batch of arrays 
    
    allocate(self%znuc(self%natoms))
    allocate(self%coords(3, self%natoms))
    allocate(self%imin(self%natoms))
    allocate(self%imax(self%natoms))
    allocate(self%evec(3))
    allocate(self%katom(self%nshell))
    allocate(self%intyp(self%nshell))
    allocate(self%kng(self%nshell))
    allocate(self%AOlocation(self%nshell))
! allocate temporary arrays 
    allocate(ityp(self%nprimi))
!    allocate(ish(self%nprimi))
    allocate(expon(self%nprimi))
    allocate(contrc1(self%nprimi))
    allocate(contrc2(self%nprimi))

    call read_basis_info(trim(basisInfoFile), self%natoms, self%nshell, self%nprimi,                  &
                         self%znuc, self%coords, self%evec, expon, contrc1,              &
                         contrc2, self%imin, self%imax, self%katom, & 
                         self%intyp, ityp, self%kng, self%AOlocation)

! calculate the total number of primitives (non-unique)
    self%nprimitives = 0
    do i = 1, self%nshell
            self%nPrimitives = self%nPrimitives + self%kng(i)
    enddo
    
    allocate(self%exponent(self%nPrimitives))
    allocate(self%coefficient(self%nPrimitives)) 
    allocate(self%shellType(self%nPrimitives))
    AOindex = 0
    do i = 1, self%natoms
        do j = self%imin(i), self%imax(i)
            AOindex = AOindex + 1 
            self%exponent(AOindex)    = expon(j)
            self%coefficient(AOindex) = contrc1(j)
            self%shellType(AOindex)   = ityp(j)
        enddo 
    enddo

    if (printLevel >= 3) call print_basis_info(self)
! deallocate temporary arrays
    deallocate(ityp, expon, contrc1, contrc2) 
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
                           expon, contrc1, contrc2, imin, imax, katom, intyp, ityp, kng, kloc)
 implicit none 
  character(len=*),  intent(in)  :: filename
  integer, intent(in)  :: natoms, nshell, nprimi
  real(dp),          intent(out) :: znuc(:), coords(:,:), evec(:)
  real(dp),          intent(out) :: expon(:), contrc1(:), contrc2(:) 
  integer, intent(out) :: imin(:), imax(:), katom(:), intyp(:), ityp(:), kng(:), kloc(:)

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
!      read(300,*) (ish(i),i=1,nprimi)
      read(300,*) (ityp(i),i=1,nprimi)
      read(300,*) (expon(i),i=1,nprimi)
      read(300,*) (contrc1(i),i=1,nprimi) 
      read(300,*) (contrc2(i),i=1,nprimi)
      read(300,*) (kng(i),i=1,nshell)
      read(300,*) (kloc(i),i=1,nshell)
  close(300)
  return 
end subroutine read_basis_info

  subroutine print_basis_info(self)
    type(basisType),  intent(in) :: self

    integer :: i, j, ij, k, ioshell, atom, styp, AOindex
    character(len=4) :: alab
    character(len=2) :: otyp

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
    write(*,'(a47,i5)') 'Number of cartesian gaussian basis functions = ', self%nPrimitives

    call print_header('Geometry in a.u.')
    write(*,'(5x,a4,3x,a10,5x,3a15)') 'Atom','Charge  ', 'X     ', 'Y     ', 'Z     ' 
    write(*,'(5x,a4,3x,a10,5x,a45)') repeat('-',4), repeat('-', 10), repeat('-', 45) 
    do i = 1, self%natoms
        write(*,'(5x,a4,3x,f10.2,5x,3f15.8)') atomLabel(int(self%znuc(i))), self%znuc(i),         &
                                        self%coords(1, i), self%coords(2, i), self%coords(3, i)
    enddo
    call print_header('Basis set')
    write(*,'(2x,"Shell",3x,"Type",3x,"Primitive",7x,"Exponent",5x,"Contraction Coefficient(s)")')
    ij = 0
    do i = 1, self%natoms
        write(*,'(/a4,a10,i5,5x,a9,f5.2/)') atomLabel(int(self%znuc(i))), 'Atom No.', i, 'Charge = ', self%znuc(i)
        do j = self%imin(i), self%imax(i)
            ij = ij + 1
            if (self%shellType(j) < 8) then
                write(*,'(1x,i6,3x,a4,i7,f22.7,f18.12,3x,i3)') & 
                    i, label(self%shellType(j)), j, self%exponent(j), self%coefficient(j), ij
            else
                write(*,'(1x,i6,3x,a4,i7,f22.7,2f18.12)') & 
                    i, label(self%shellType(j)), j, self%exponent(j), self%coefficient(j)
            endif
        enddo
    enddo

    call print_header('Basis set in detail')
    write(*,'(" AO ",2x,"atom",3x,"atom",2x,"shell",3x,"shell",2x,"subshell",12x,"         ",14x,"contraction")') 
    write(*,'("indx",2x,"name",3x,"No. ",2x,"No.  ",3x,"type ",2x,"  type  ",12x,"exponents",14x,"coefficients")') 
    write(*,'(a)') repeat('-', 89)
    ij = 0
    do i = 1, self%nshell
        atom = self%katom(i)
        do j = 1, self%kng(i)
                ij  = ij + 1
                ioshell = 0 
            do k = kmin(self%intyp(i)), kmax(self%intyp(i))
                ioshell = ioshell + 1
                alab = atomlabel(int(self%znuc(self%katom(i))))
                styp = self%intyp(i)
                otyp = label(self%intyp(i)) 
                AOindex = self%AOlocation(i) + ioshell - 1
                write(*,'(i4,2x,a4,3x,i4,3x,i4,5x,a2,3x,a4,3x,2f24.10,3x,i3)') &
                    AOindex, alab, atom, i, otyp, bfnam1(k), self%exponent(ij), self%coefficient(ij)
            enddo
        enddo
    enddo

  end subroutine print_basis_info

  real(dp) function normalizeGaussian(ltype, mtype, zeta)
!===============================================================================
! calculate normalization constant fo a given cartesian primitive AO
!
!   input  :
!       ltype : l quantum number + 1, meaning 1, 2, 3, 4, 5, 6, 7, 8
!                                         for S, P, D, F, G, H, I, L
!               where L is the SP shell, ceses with ltype > 5 are not handled
!       mtype : number for the component of the primitive as stored in bfnam1
!               (see preamble for this module 
!   output :
!       normalization cosntant for a given AO
!===============================================================================
    integer,  intent(in) :: ltype, mtype
    real(dp), intent(in) :: zeta
    
    real(dp) :: norm = 0.0_dp

  select case (ltype)
    case (1)       ! S
        norm = (2.0_dp/Pi)**0.75_dp/dSqrt(zeta**(-1.5_dp))
    case (2)       ! P
        norm = (2.0_dp*(2.0_dp/Pi)**0.75_dp)/dSqrt(zeta**(-2.5_dp))
    case (3)       ! D
        select case (mtype)
            case (5:7)         ! D xx, yy, zz
                norm = (4.0_dp*(2.0_dp/Pi)**0.75_dp)/(Sqrt(3.0_dp)*dSqrt(zeta**(-3.5_dp)))
            case (8:10)        ! D xy, yz, zx
                norm = (4.0_dp*(2.0_dp/Pi)**0.75_dp)/dSqrt(zeta**(-3.5_dp))
            case default 
                write(*,*) "something wrong with mtype for D shell: ", mtype
        end select
    case (4)       ! F
        select case (mtype)
            case (11:13)       ! F xxx, yyy, zzz
                norm = (8.0_dp*(2.0_dp/Pi)**0.75_dp)/(dSqrt(15.0_dp)*dSqrt(zeta**(-4.5_dp)))
            case (14:19)       ! F xxy, xxz, yyx, yyz, zzx, zzy
                norm = (8.0_dp*(2.0_dp/Pi)**0.75_dp)/(dSqrt(3.0_dp)*dSqrt(zeta**(-4.5_dp)))
            case (20)          ! F xyz
                norm = (8.0_dp*(2.0_dp/Pi)**0.75_dp)/dSqrt(zeta**(-4.5_dp))
            case default 
                write(*,*) "something wrong with mtype for F shell: ", mtype
        end select
    case (5)       ! G
        select case (mtype)
            case (21:23)       ! G xxxx, yyyy, zzzz
                norm = (16.0_dp*(2.0_dp/Pi)**0.75_dp)/(dSqrt(105.0_dp)*dSqrt(zeta**(-5.5_dp)))
            case (24:29)       ! G xxxy, xxxz, yyyx, yyyz, zzzx, zzzy
                norm = (16.0_dp*(2.0_dp/Pi)**0.75_dp)/(dSqrt(15.0_dp)*dSqrt(zeta**(-5.5_dp)))
            case (30:32)       ! G xxyy, xxzz, yyzz
                norm = (16.0_dp*(2.0_dp/Pi)**0.75_dp)/(3.0_dp*dSqrt(zeta**(-5.5_dp)))
            case (33:35)       ! G xxyz, yyxz, zzxy
                norm = (16.0_dp*(2.0_dp/Pi)**0.75_dp)/(dSqrt(3.0_dp)*dSqrt(zeta**(-5.5_dp)))
            case default 
                write(*,*) "something wrong with mtype for G shell: ", mtype
        end select
    case default 
        write(*,*) 'Something wrong with the ltype in -normalizeGaussian-: ', ltype
    end select

    normalizeGaussian = norm

  end function normalizeGaussian

  subroutine deleteBasis(self)
    type(basisType) :: self

    if (allocated(self%exponent))    deallocate(self%exponent) 
    if (allocated(self%coefficient)) deallocate(self%coefficient)
    if (allocated(self%shellType))   deallocate(self%shellType)
    if (allocated(self%AOlocation))  deallocate(self%AOlocation)
    if (allocated(self%katom))       deallocate(self%katom)
    if (allocated(self%intyp))       deallocate(self%intyp)
    if (allocated(self%znuc))        deallocate(self%znuc)
    if (allocated(self%coords))      deallocate(self%coords)
    if (allocated(self%imin))        deallocate(self%imin)
    if (allocated(self%imax))        deallocate(self%imax)
    if (allocated(self%evec))        deallocate(self%evec)
  end subroutine deleteBasis

end module basisModule
