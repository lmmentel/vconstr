module BeckeGridModule

   use Vartypes
   use AtomicPartitionFun
   use LebedevAngGridModule
   use M4RadialGridModule
   use BeckeGridConfigModule

   implicit none

   private
   save

   public :: CreateBeckeGrid

   ! ========================================================================
   ! Purpose: generation of points and weights for numerical integration 
   !          in polyatomic systems according to the becke grid (fuzzy cells) 
   !          scheme. 
   !          Reference papers:
   !          - [1] Becke, J. Chem. Phys. 88, 2547 (1988)
   !          - [2] Treutler O., Ahlrichs R., J. Chem. Phys. 102, 346 (1995)
   !          - [3] V. Lebedev, D. Laikov, Doklady Math. 59 (1999) 477–481
   !
   ! Implementation details:
   !          - Radial mapping: M4 mapping (ref.[2])
   !          - Angular Grid: Lebedev grid (ref.[3])
   !
   ! ========================================================================
   
   integer(KINT), parameter             :: maxRadPoints = 50000_KINT    ! 1.3^6 * 1000 \approx 5000 (just to be sure *10)
   integer(KINT), parameter             :: maxAngPoints = 50000_KINT    ! it has to be larger then 48*36*MaxNrOfOper

   !   (\__/)  .~    ~. )) ! 
   !   /O O  ./      .'    ! 
   !  {O__,   \    {       ! 
   !    / .  . )    \      ! 
   !    |-| '-' \    } ))  !
   !   .(   _(   )_.'      !
   !  '---.~_ _ _          !

contains

   subroutine CreateBeckeGrid(beckeConfig, iupto1, totNumPoints, minimumRadius)

      ! ============================================================
      ! purpose : creates and writes to file the points and weights.
      ! ============================================================

      type(BeckeConfigType),    intent(in)  :: beckeConfig   ! object containing the configuration info
      integer(KINT),            intent(in)  :: iupto1        ! the points and weights will be written to this tape
      integer(KINT),            intent(out) :: totNumPoints 
      real(KREAL)  , optional,  intent(in)  :: minimumRadius ! minimum allowed distance of a point from an atom

      real(KREAL), allocatable   :: radPoints(:), radWeights(:), angPoints(:,:), angWeights(:) 
      real(KREAL), allocatable   :: angPointsUnitSph(:,:), angWeightsUnitSph(:) 
      real(KREAL), allocatable   :: pntsAndWeights(:,:), partitionFunction(:), localOper(:,:,:)
      real(KREAL)                :: octaOper(3,3,48)
      integer(KINT)              :: numAngPoints, numLocalOper, numOctaOper, actualNumPoints
      integer(KINT)              :: iRadShell, iAtom, ipnt, istat, previousAngOrder

   !  -------------------------------------------------------------------------------------
   !  parallel execution: I use the same parallel implementation used by the voronoi scheme
   !  -------------------------------------------------------------------------------------
  
      allocate(radPoints(maxRadPoints),  radWeights(maxRadPoints), angPoints(3,maxAngPoints), &
               angWeights(maxAngPoints), partitionFunction(maxAngPoints), angPointsUnitSph(3,maxAngPoints), &
               angWeightsUnitSph(maxAngPoints), STAT=istat)

      totNumPoints      = 0_KINT

      ! -----------------------------------------------------------
      ! Main loop over the atoms. 
      !
      ! Skip if it's not a symmetry unique atom.
      ! (parallel generation of the points. The tape15 will be combined later)
      ! ------------------------------------------------------------
      
      do iAtom=1, beckeConfig%nAtoms
         if ( beckeConfig%isSymmetryUnique(iAtom) ) then

         ! --------------------------------
         ! Generation of the radial points:
         ! --------------------------------
      
         call CreateM4RadialGrid(radPoints, radWeights, beckeConfig%numRadPoints(iAtom), &
                                 beckeConfig%beckeMapParams(iAtom))
 
         if (present(minimumRadius)) then 
            call fixM4MinimumRadius( beckeConfig%numRadPoints(iAtom), radPoints, radWeights, minimumRadius )
         end if

         ! Determine the symmetry operators of this particular atom:
         
!          call countLocalOperators(beckeConfig%xyzAtoms(:,iAtom), beckeConfig%oper, beckeConfig%transl,     &
!                                   beckeConfig%latticeDimension, beckeConfig%vectors, beckeConfig%kvectors, &
!                                   beckeConfig%nOper, numLocalOper)   
         numLocalOper = 1_KINT   
         allocate(localOper(3,3,numLocalOper))
!          call fillLocalOperators (beckeConfig%xyzAtoms(:,iAtom), beckeConfig%oper, beckeConfig%transl,     &
!                                   beckeConfig%latticeDimension, beckeConfig%vectors, beckeConfig%kvectors, &
!                                   localOper)
         localOper(:,:,1) = beckeConfig%oper (:,:,1)

         ! --------------------------------------------------------------
         ! Generation of the Angular points:
         !
         ! For each radial shell (r), generate the xyz points and weights
         ! --------------------------------------------------------------

         previousAngOrder = -1_KINT

         do iRadShell=1, beckeConfig%numRadPoints(iAtom) 
           
            ! Genetaye the "xyz" points on the unit sphere (only if I don't have them already)
            
            if (beckeConfig%angLOrder(iAtom,iRadShell) /= previousAngOrder) then
               call CreateLebedevPoints (beckeConfig%angLOrder(iAtom,iRadShell), &
                                         numAngPoints, angPointsUnitSph, angWeightsUnitSph, symOperIn=localOper, &
                                         symToleranceIn=beckeConfig%symTolerance,           &
                                         kvectors=beckeConfig%kvectors, allowSpecialGrid=beckeConfig%allowSpecialGrid)
               previousAngOrder = beckeConfig%angLOrder(iAtom,iRadShell)
            end if

            ! From the unit sphere (r=1) to the real one:
            
            angPoints(:,1:numAngPoints) = angPointsUnitSph(:,1:numAngPoints) * radPoints(iRadShell)
            angWeights(1:numAngPoints)  = angWeightsUnitSph(1:numAngPoints)  * radWeights(iRadShell)

            ! From the origin-centered system to the real atom
            
            do ipnt=1, numAngPoints
               angPoints(:,ipnt) = angPoints(:,ipnt) + beckeConfig%xyzAtoms(:,iAtom)
            end do
            
            ! ----------------------------------------------------------------------------------
            ! Calculate the becke-partition weights-correction (size adjusted scheme in ref.[1])
            ! ----------------------------------------------------------------------------------
 
            call CalcPartition( angPoints(:,1:numAngPoints), beckeConfig%xyzAtomsSupercell,             &
                                beckeConfig%qAtomsSupercell, iAtom, partitionFunction(1:numAngPoints),  &
                                beckeConfig%partitionFunction, .true. )

            angWeights(1:numAngPoints) = angWeights(1:numAngPoints) * partitionFunction(1:numAngPoints)

            ! ---------------------
            ! Write points to file:
            ! ---------------------
           
            allocate( pntsAndWeights(4, numAngPoints) )

            actualNumPoints = 0_KINT
            do ipnt=1, numAngPoints
                if (partitionFunction(ipnt) > beckeConfig%partitionFunThresh) then
                  actualNumPoints = actualNumPoints + 1
                  pntsAndWeights(1:3,actualNumPoints) = angPoints(:,ipnt)
                  pntsAndWeights(4,actualNumPoints)   = angWeights(ipnt)
               end if
            end do            
   
            if (actualNumPoints>0_KINT) then
               do ipnt=1, actualNumPoints
                  write(iupto1,"(4e25.14)") pntsAndWeights(:,ipnt)
               end do
               totNumPoints = totNumPoints + actualNumPoints
            end if

            deallocate( pntsAndWeights )
            
         end do

         deallocate(localOper)

         end if
      end do
      
      deallocate(angPoints, angWeights, radPoints, radWeights, partitionFunction, &
                 angPointsUnitSph, angWeightsUnitSph, STAT=istat)

   end subroutine


!    subroutine CountLocalOperators(point, oper, transl, latticeDimension, vectors, kvectors, nOper, numLocalOper)

!       real(KREAL),     intent(in)   :: point(3)
!       real(KREAL),     intent(in)   :: oper(:,:,:)       ! oper(3,3,nOper)
!       real(KREAL),     intent(in)   :: transl(:,:)       ! oper(3,nOper)
!       integer(KINT),   intent(in)   :: latticeDimension
!       real(KREAL),     intent(in)   :: vectors(:,:)      
!       real(KREAL),     intent(in)   :: kvectors(:,:)     
!       integer(KINT),   intent(out)  :: numLocalOper
      
!       integer(KINT)         :: nOper, iOper
!       logical, external     :: symequ

!       nOper = size(oper(1,1,:))

!       numLocalOper = 0_KINT

!       do iOper = 1, nOper
!          if (symequ(point,point,oper(:,:,iOper), transl(1,iOper),&
!              latticeDimension, vectors, kvectors)) then  
!                numLocalOper = numLocalOper + 1
!          end if
!       end do
!    end subroutine


!    subroutine FillLocalOperators(point, oper, transl, latticeDimension, vectors, kvectors, localOper)

!       real(KREAL),     intent(in)   :: point(3)
!       real(KREAL),     intent(in)   :: oper(:,:,:)          ! oper(3,3,nOper)
!       real(KREAL),     intent(in)   :: transl(:,:)          ! transl(3,nOper)
!       integer(KINT),   intent(in)   :: latticeDimension
!       real(KREAL),     intent(in)   :: vectors(:,:)      
!       real(KREAL),     intent(in)   :: kvectors(:,:)     
!       real(KREAL),     intent(out)  :: localOper(:,:,:)     ! localOper(3,3,numLocalOper)

!       integer(KINT)         :: nOper, numLocalOper, iOper, iLocalOper
!       logical, external     :: symequ

!       nOper = size(oper(1,1,:))
!       numLocalOper = size(localOper(1,1,:))

!       iLocalOper = 0_KINT

!       do iOper = 1, nOper
!          if (symequ(point,point,oper(:,:,iOper), transl(1,iOper),&
!              latticeDimension, vectors, kvectors)) then  
!                iLocalOper = iLocalOper + 1
!                localOper(:,:,iLocalOper) = oper(:,:,iOper)
!          end if
!       end do
!    end subroutine

end module
