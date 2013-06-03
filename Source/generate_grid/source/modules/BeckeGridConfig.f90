module BeckeGridConfigModule

   use Vartypes
   use M4RadialGridModule


   implicit none

   private

   public :: BeckeConfigType
   public :: New, Delete
   public :: PrintInfo

   interface New
      module procedure NewBeckeConfig
   end interface

   interface Delete
      module procedure DeleteBeckeConfig
   end interface

   interface PrintInfo
      module procedure PrintBeckeConfigInfo
   end interface



   type BeckeConfigType

      ! geometry info:                                                                
      real(KREAL)                :: symTolerance
      real(KREAL)                :: partitionFunThresh
      integer(KINT)              :: nAtoms
      integer(KINT)              :: nAtomsSupercell
      integer(KINT)              :: nOper
      real(KREAL),   allocatable :: xyzAtoms(:,:)               ! size ( 3,nAtoms )
      real(KREAL),   allocatable :: qAtoms(:)                   ! size ( nAtoms )
      real(KREAL),   allocatable :: xyzAtomsSupercell(:,:)      ! size ( 3,nAtoms*nCells )
      real(KREAL),   allocatable :: qAtomsSupercell(:)          ! size ( nAtoms*nCells )
      logical,       allocatable :: isSymmetryUnique(:)         ! size ( nAtoms )
      real(KREAL),   allocatable :: oper(:,:,:)                 ! size ( 3,3,nOper )
      real(KREAL),   allocatable :: transl(:,:)                 ! size ( 3,nOper )
      real(KREAL)                :: stdRotation(3,3) 
      integer(KINT)              :: latticeDimension 
      real(KREAL)                :: vectors(3,3)
      real(KREAL)                :: kvectors(3,3)
      
      ! variables for M04 mapping:
      integer(KINT), allocatable :: numRadPoints(:)         ! size ( nAtoms )
      real(KREAL),   allocatable :: beckeMapParams(:)       ! size ( nAtoms )
      
      ! variables for lebedev:                                                        
      integer(KINT), allocatable :: angLOrder(:,:)          ! size ( nAtoms,radPoints )

      character(LCHARS)          :: partitionFunction
      character(LCHARS)          :: quality 
      logical                    :: allowSpecialGrid

   end type


   integer(KINT), parameter             :: maxRadPoints = 50000_KINT    ! 1.3^6 * 1000 \approx 5000 (just to be sure *10)
   real(KREAL)  , parameter             :: zero    = 0.0_KREAL
   real(KREAL)  , parameter             :: one     = 1.0_KREAL

contains


   subroutine NewBeckeConfig (self, xyzAtoms, qAtoms, angLOrder, nRadialPoints, pruning)
      
      type(BeckeConfigType),    intent(inout)   :: self
      real(KREAL),              intent(in)      :: xyzAtoms(:,:)    ! xyzAtoms(3,nAtoms)
      real(KREAL),              intent(in)      :: qAtoms(:)        ! qAtoms(nAtoms)
      integer(KINT),            intent(in)      :: angLOrder        
      integer(KINT),            intent(in)      :: nRadialPoints    
      logical,                  intent(in)      :: pruning 

      integer(KINT)     :: nAtoms, istat, iAtom, ipnt
      real(KREAL)       :: oper(3,3,1) 

      nAtoms = size(qAtoms)

      oper = 0.0_KREAL
      oper(1,1,1) = 1.0_KREAL
      oper(2,2,1) = 1.0_KREAL
      oper(3,3,1) = 1.0_KREAL
   
      self%noper    = 1_KINT
      self%nAtoms   = nAtoms
      self%nAtomsSupercell = nAtoms
      self%stdRotation = reshape((/one, zero, zero, zero, one, zero, zero, zero, one/), shape(self%stdRotation))

      allocate(self%numRadPoints(self%nAtoms), self%angLOrder(self%nAtoms,maxRadPoints)                  , &
               self%beckeMapParams(self%nAtoms), self%xyzAtoms(3,self%nAtoms), self%qAtoms(self%nAtoms)  , &
               self%xyzAtomsSupercell(3,self%nAtomsSupercell), self%qAtomsSupercell(self%nAtomsSupercell), & 
               self%isSymmetryUnique(self%nAtoms), self%oper(3,3,self%nOper)                             , &
               self%transl(3,self%nOper), STAT=istat)
      
      self%xyzAtoms           = xyzAtoms
      self%qAtoms             = qAtoms
      self%xyzAtomsSupercell  = xyzAtoms
      self%qAtomsSupercell    = qAtoms
      self%isSymmetryUnique   = .true.
      self%oper               = oper
      self%transl             = 0.0_KREAL
      self%symTolerance       = 1.0e-7_KREAL
      self%partitionFunThresh = 1.0E-8_KREAL
      self%angLOrder          = angLOrder
      self%partitionFunction  = 'YUKAWALIKE'
      self%allowSpecialGrid   = .false.
      self%quality            = 'Normal'
      self%latticeDimension   = 0_KINT
      self%vectors            = 0_KINT
      self%kvectors           = 0_KINT

      do iAtom=1, nAtoms
         
         self%beckeMapParams(iAtom) = M4RadialGridParam%mapParam( int(self%qAtoms(iAtom),KINT  ))

         ! Magical fourmula:
         self%numRadPoints(iAtom)  = nRadialPoints
        
         ! Pruning :
         if (pruning) then
            do ipnt=1, self%numRadPoints(iAtom)
               if (ipnt <= int(self%numRadPoints(iAtom)*0.3)) self%angLOrder(iAtom,ipnt) = 11_KINT
            end do
         end if
      end do

   end subroutine 


         
   subroutine DeleteBeckeConfig(self)
      
      type(BeckeConfigType), intent(inout)    :: self

      integer(KINT)     :: istat

      deallocate( self%numRadPoints, self%angLOrder, self%beckeMapParams,   &
                  self%xyzAtoms, self%qAtoms, self%xyzAtomsSupercell, self%qAtomsSupercell,  &
                  self%isSymmetryUnique, self%oper, self%transl, STAT=istat)

   end subroutine

   subroutine PrintBeckeConfigInfo (self)
      type (BeckeConfigType),       intent(in)  :: self

      write(*,9000) self%quality
      write(*,*) " "
      write(*,9010) minval(self%angLOrder(:,1)), maxval(self%angLOrder)
      write(*,9020) minval(self%numRadPoints), maxval(self%numRadPoints)
      
      9000 format(' Becke grid quality: ',t30,a15)
      9010 format(' Lebedev angular grid order range: 't45,'from ', i3,' to ',i3)
      9020 format(' Nr. of radial points range: 't45,'from ', i3,' to ',i3)

      if (.true.) then
         write(*,*) 'BeckeGrid symTolerance', self%symTolerance 
         write(*,*) 'BeckeGrid partitionFunThresh', self%partitionFunThresh 
         write(*,*) 'BeckeGrid nAtoms', self%nAtoms 
         write(*,*) 'BeckeGrid nAtomsSupercell', self%nAtomsSupercell 
         write(*,*) 'BeckeGrid nOper', self%nOper 
         write(*,*) 'BeckeGrid numRadPoints', self%numRadPoints 
         write(*,*) 'BeckeGrid beckeMapParams', self%beckeMapParams 
        ! write(*,*) 'BeckeGrid angLOrder', self%angLOrder 
         write(*,*) 'BeckeGrid xyzAtoms', self%xyzAtoms 
         write(*,*) 'BeckeGrid qAtoms', self%qAtoms 
         write(*,*) 'BeckeGrid xyzAtomsSupercell', self%xyzAtomsSupercell 
         write(*,*) 'BeckeGrid qAtomsSupercell', self%qAtomsSupercell 
         write(*,*) 'BeckeGrid isSymmetryUnique', self%isSymmetryUnique 
         write(*,*) 'BeckeGrid oper', self%oper 
         write(*,*) 'BeckeGrid transl', self%transl 
         write(*,*) 'BeckeGrid stdRotation', self%stdRotation 
         write(*,*) 'BeckeGrid partitionFunction', self%partitionFunction 
         write(*,*) 'BeckeGrid allowSpecialGrid', self%allowSpecialGrid 
         write(*,*) 'BeckeGrid quality', self%quality 
         write(*,*) 'BeckeGrid latticeDimension', self%latticeDimension 
         write(*,*) 'BeckeGrid vectors', self%vectors 
         write(*,*) 'BeckeGrid kvectors', self%kvectors 
      end if

   end subroutine 

end module
