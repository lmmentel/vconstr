module dictionaryModule
  use varModule

  implicit none 

  private

  type dictionaryType
      integer(IPgamess) :: idaf
      integer(IPgamess) :: irecln
      integer(IPgamess) :: irecst
      integer(IPgamess) :: ioda(950)
      integer(IPgamess) :: ifilen(950)
      integer(IPgamess) :: is
      integer(IPgamess) :: ipk
  end type
      
  public :: dictionaryType

  interface new
      module procedure newPrivate
  end interface
  interface delete
      module procedure deletePrivate
  end interface

  public :: new, delete

  interface readFrom
      module procedure ReadReals
      module procedure ReadRealMatrix
      module procedure ReadIntegers
  end interface

  public :: readFrom
contains

 subroutine newPrivate(self)
!===============================================================================
! create a dictionary type object that will hold the basic information about
! gamess-us dictionary file, open the dictionary file and read the first 
! record containing irect, ioda, ifilen, is , ipk varibles and set all the 
! remaining required variables, structure of the dictionary type is given in 
! the preamble of this module
!===============================================================================
  use commonsModule
  type(dictionaryType) :: self
  integer(IPgamess) IRECLN,IRECST,IFILEN(950)
  integer(IPgamess) IR,IW,IP,IIS,IPK,IDAFX,NAV,IODAX(950)
  integer(IPgamess) ME, MASTER, NPROC, IBTYP, IPTIM 
  logical(LKgamess) GOPARR, DSKWRK, MASWRK
  COMMON /DAIOLN/ IRECLN,IRECST,IFILEN
  COMMON /IOFILE/ IR,IW,IP,IIS,IPK,IDAFX,NAV,IODAX
  COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

  self%idaf   = 10_IPgamess
  self%irecln = 4090_IPgamess

  open(unit=self%idaf,file=trim(dictionaryFile),access='direct',recl=8*self%irecln,form='unformatted')
  read(10,rec=1) self%irecst, self%ioda, self%ifilen, self%is, self%ipk
  ! set the required defaults to pass to the gamess-us common blocks
  IRECLN = self%irecln
  IRECST = self%irecst
  IFILEN = self%ifilen
  IR     = 5_IPgamess
  IW     = 6_IPgamess
  IP     = 7_IPgamess
  ME     = 0_IPgamess
  MASTER = 0_IPgamess
  NPROC  = 1_IPgamess
  GOPARR = .false.
  DSKWRK = .true.
  MASWRK = .true.
 end subroutine newPrivate 

 subroutine deletePrivate(self)
  type(dictionaryType) :: self

  close(self%idaf)
 end subroutine deletePrivate

 subroutine ReadReals(self, section, array)
    type(dictionaryType), intent(in)  :: self
    character(*),         intent(in)  :: section
    real(DP),             intent(out) :: array(:)
    integer(IPgamess) ME, MASTER, NPROC, IBTYP, IPTIM 
    logical(LKgamess) GOPARR, DSKWRK, MASWRK
    COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

    interface
      subroutine DAREAD(idaf, ioda, firstElement, length, section, realOrInt)
        use VarModule

        integer(IPgamess), intent(in)  :: idaf, ioda(950)
        real(DPgamess),    intent(out) :: firstElement
        integer(IPgamess), intent(in)  :: length, section, realOrInt
      end subroutine
    end interface

!   Ordered alphabetically
    select case(trim(section))
      case ('alpha energies')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 17_IPgamess, 0_IPgamess)
      case ('core hamiltonian')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 11_IPgamess, 0_IPgamess)
      case ('energies')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess),  2_IPgamess, 0_IPgamess)
      case ('kinetic energy')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 13_IPgamess, 0_IPgamess)
      case ('occupation numbers')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 21_IPgamess, 0_IPgamess)
      case ('overlap')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 12_IPgamess, 0_IPgamess)
      case ('X dipole')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 95_IPgamess, 0_IPgamess)
      case ('Y dipole')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 96_IPgamess, 0_IPgamess)
      case ('Z dipole')
        call DAREAD(self%idaf, self%ioda, array(1), int(size(array), kind=IPgamess), 97_IPgamess, 0_IPgamess)
      case default
        write(stderr,*) 'No valid section for reading the dictionary file selected'
        call abrt()
    end select
 end subroutine

 subroutine ReadRealMatrix(self, section, matrix)
    type(DictionaryType), intent(in)  :: self
    character(*),         intent(in)  :: section
    real(DP),             intent(out) :: matrix(:,:)
!    integer(IPgamess) ME, MASTER, NPROC, IBTYP, IPTIM 
!    logical(LKgamess) GOPARR, DSKWRK, MASWRK
!    COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

    interface
      subroutine DAREAD(idaf, ioda, firstElement, length, section, realOrInt)
        use VarModule

        integer(IPgamess),  intent(in)  :: idaf, ioda(950)
        real(DPgamess),     intent(out) :: firstElement
        integer(IPgamess),  intent(in)  :: length, section, realOrInt
      end subroutine
    end interface

    select case(trim(section))
      case ('guess')
        call DAREAD(self%idaf, self%ioda, matrix(1,1), int(size(matrix), kind=IPgamess), 265_IPgamess, 0_IPgamess)
      case ('HF MOs')
        call DAREAD(self%idaf, self%ioda, matrix(1,1), int(size(matrix), kind=IPgamess),  15_IPgamess, 0_IPgamess)
      case ('natural orbitals')
        call DAREAD(self%idaf, self%ioda, matrix(1,1), int(size(matrix), kind=IPgamess),  19_IPgamess, 0_IPgamess)
      case ('salc')
        call DAREAD(self%idaf, self%ioda, matrix(1,1), int(size(matrix), kind=IPgamess),  44_IPgamess, 0_IPgamess)
      case ('orthonormal symmetrized AOs')
        call DAREAD(self%idaf, self%ioda, matrix(1,1), int(size(matrix), kind=IPgamess),  45_IPgamess, 0_IPgamess)
      case default
        write(stderr,*) 'No valid section for reading the dictionary file selected in -ReadRealMatrix-'
        call abrt()
    end select
 end subroutine
 

 subroutine ReadIntegers(self, section, array)
!   The reading of integers is more tricky, since they are stored as reals
    type(DictionaryType), intent(in)  :: self
    character(*),         intent(in)  :: section
    integer(IK),          intent(out) :: array(:)

    integer(IPgamess), allocatable :: gamessInts(:)

    interface
      subroutine DAREAD(idaf, ioda, firstElement, length, section, realOrInt)
        use VarModule

        integer(IPgamess), intent(in)  :: idaf, ioda(950)
        integer(IPgamess), intent(out) :: firstElement
        integer(IPgamess), intent(in)  :: length, section, realOrInt
      end subroutine
    end interface

    allocate(gamessInts(size(array)))

    select case(trim(section))
      case ('hoi')
        call DAREAD(self%idaf, self%ioda, gamessInts(1), int(size(array), kind=IPgamess), -1_IPgamess, 1_IPgamess)
       case default
        write(stderr,*) 'No valid section for reading the dictionary file selected'
          call abrt()
    end select

    array = int(gamessInts, kind=IK)

    deallocate(gamessInts)
 end subroutine
end module dictionaryModule
