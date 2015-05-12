program MyGrid

   ! =================
   ! DO NOT DISTRIBUTE
   ! =================

  use Vartypes
  use BeckeGridConfigModule
  use BeckeGridModule
  use BasisModule

    implicit none

    type(BeckeConfigType) :: config
    type(basisType)       :: basis

    real(KREAL), allocatable  :: xyzAtoms(:,:)   ! size(3,nAtoms)
    real(KREAL), allocatable  :: qAtoms(:)       ! size(nAtoms)
    real(KREAL), allocatable :: xyzw(:,:), r(:), fun(:)
    integer(KINT)     :: nAtoms, totNumPoints, iounit, angularLOrder, nRadialPoints
    real(KREAL)       :: sumWeights, ipnt, integral
    logical           :: pruning, debug

    real(KREAL),   parameter :: pi=3.1415926535897931_KREAL
    integer(KINT), parameter :: PossibleLebedevGridOrders(15) = (/5,7,11,13,17,19,21,23,29,31,35,41,47,53,59/)

    integer(KINT),  parameter :: nlunit = 10
    integer(KINT)             :: ios
    character(len=100)        :: buffer, inputfile, basinfo, gridfile
    namelist /input/ basinfo, gridfile, angularLOrder, nRadialPoints, pruning, debug

    ! set the default values

    gridfile = 'IntegrationPoints.xyzw'
    debug = .false.

    ! =====================================================================
    ! Grid parameters:
    ! - AngularLOrder is the order of the ochtehedral lebedev grid,
    ! and it shoud be one of the number listed in PossibleLebedevGridOrders
    ! - nRadialPoints is the number of radial points/atom
    ! - With the "pruning" option on it will use a coarser angular grid
    ! for points near the nuclei.
    ! =====================================================================

    angularLOrder  = 29_KINT
    nRadialPoints  = 100_KINT
    pruning        = .true.

    ! read the input file with the details

    call getarg(1, buffer)
    read(buffer, *) inputfile
    open(unit=nlunit, file=trim(inputfile), status='old', form='formatted',   &
&        delim='apostrophe', iostat=ios, action='read')
    read(nlunit, nml=input)
    close(nlunit)

    if (debug) then
        write(*,*) 'basinfo  : ', basinfo
        write(*,*) 'gridfile : ', gridfile
        write(*,*) 'angularLOrder : ', angularLOrder
        write(*,*) 'nRadialPoints : ', nRadialPoints
        write(*,*) 'pruning : ', pruning
    endif

    ! the system information will be read from the gamess-us $JOB.basinfo file
    ! specified in the input file

    call newBasis(basis, basinfo, debug)

   iounit = 666
   open (unit=iounit, file=trim(gridfile), action='write' )

   call New (config, basis%coords, basis%znuc, angularLOrder, nRadialPoints, pruning)

   if (debug) call PrintInfo (config)

   call CreateBeckeGrid (config, iounit, totNumPoints)

    call deleteBasis(basis)
   call Delete (config)
   close (iounit)

end program

