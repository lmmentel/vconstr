program MyGrid

   ! =================
   ! DO NOT DISTRIBUTE
   ! =================

   use Vartypes
   use BeckeGridConfigModule
   use BeckeGridModule

   implicit none

   type(BeckeConfigType) :: config

   real(KREAL), allocatable  :: xyzAtoms(:,:)   ! size(3,nAtoms)
   real(KREAL), allocatable  :: qAtoms(:)       ! size(nAtoms)
   real(KREAL), allocatable :: xyzw(:,:), r(:), fun(:)
   integer(KINT)     :: nAtoms, totNumPoints, iounit, angularLOrder, nRadialPoints
   real(KREAL)       :: sumWeights, ipnt, integral
   logical           :: pruning, debug 
   character(LCHARS) :: filename

   real(KREAL),   parameter :: pi=3.1415926535897931_KREAL
   integer(KINT), parameter :: PossibleLebedevGridOrders(15) = (/5,7,11,13,17,19,21,23,29,31,35,41,47,53,59/)


   filename = 'IntegrationPoints.xyzw'
   debug = .false.

   ! ========================================
   ! Here you should define your system: 
   ! (in this example is just a simple dimer)
   ! ========================================

   nAtoms = 2
   allocate(xyzAtoms(3,nAtoms), qAtoms(nAtoms))
   qAtoms         = 1.0_KREAL
   xyzAtoms       = 0.0_KREAL
   xyzAtoms(3,2)  = 1.0_KREAL

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

   ! ================================================
   ! Create (and write to file) the integration grid.
   ! (Format: x y z integrationWeight)
   ! ================================================

   iounit = 666
   open (unit=iounit, file=filename, action='write' )

   call New (config, xyzAtoms, qAtoms, angularLOrder, nRadialPoints, pruning)

   if (debug) call PrintInfo (config)

   call CreateBeckeGrid (config, iounit, totNumPoints)

   write(*,"(A,A25)") " Grid written to file: ", filename 
   write(*,*) "Total number of integration points: ", totNumPoints

   call Delete (config)
   close (iounit)
   deallocate(xyzAtoms, qAtoms)

   ! ==============
   ! Test integral: 
   ! ==============

   open (unit=iounit, file=filename, action='read' )
   allocate(xyzw(totNumPoints,4), r(totNumPoints), fun(totNumPoints))

   do ipnt=1, totNumPoints
      read(iounit,*) xyzw(ipnt,:)
   end do

   r(:) = sqrt(xyzw(:,1)**2 + xyzw(:,2)**2 + xyzw(:,3)**2)
   fun = exp(-r)
   integral = sum( fun*xyzw(:,4) ) / (4.0_KREAL*pi)

   write(*,*) "Test integral, relative error:", (2.0_KREAL-integral)/2.0_KREAL

   deallocate(xyzw,r,fun)
   close (iounit)



end program 

