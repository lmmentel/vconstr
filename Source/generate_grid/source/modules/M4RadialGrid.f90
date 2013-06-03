module M4RadialGridModule

   use vartypes
   use LegendreGridModule
   implicit none
  
   private
   public :: CreateM4RadialGrid, fixM4MinimumRadius
   public :: M4RadialGridParam
   
   ! -------------------------------------------------------------------
   ! Purpose: Create a radial grid [0:infinity] using the radial mapping 
   !          function M4 proposed in: 
   !          Treutler O., Ahlrichs R., J. Chem. Phys. 102, 346 (1995)  
   !          - Legendre mapping of the intervall [-1,1]
   !          - The weights already contain the jacobian r**2
   ! -------------------------------------------------------------------

   real(KREAL) , parameter              :: one     = 1.0_KREAL
   real(KREAL) , parameter              :: two     = 2.0_KREAL

   type M4RadialGridParamType
      real(KREAL)                          :: mapParam(120) = (/&
        ! H     He
        0.80, 0.9,                                              &
        ! Li    Be    B     C     N     O     F     Ne
        1.80, 1.40, 1.30, 1.10, 0.90, 0.90, 0.90, 0.90,         &
        ! Na    Mg    Al    Si    P     S     Cl    Ar
        1.40, 1.30, 1.30, 1.20, 1.10, 1.00, 1.00, 1.00,         &
        ! K     Ca    Sc    Ti    V     Cr    Mn    Fe    Co
        1.50, 1.40, 1.30, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20,   &
        ! Ni    Cu    Zn    Ga    Ge    As    Se    Br    Kr
        1.10, 1.10, 1.10, 1.10, 1.00, 0.90, 0.90, 0.90, 0.90,   &
        ! Rb    Sr    Y     Zr    Nb    Mo    Tc    Ru    Rh
        1.40, 1.40, 1.10, 1.30, 1.00, 1.20, 0.90, 0.90, 0.90,   &
        ! Pd    Ag    Cd    In    Sn    Sb    Te    I     Xe
        1.00, 0.90, 1.00, 1.00, 1.30, 1.20, 1.20, 0.90, 1.00,   &
        ! Cs    Ba    
        1.7,  1.5,                                              &
        ! La    Ce    Pr    Nd    Pm    Sm    Eu    Gd
        1.50, 1.30, 1.30, 1.40, 1.80, 1.40, 1.20, 1.30,         &
        ! Tb    Dy    Ho    Er    Tm    Yb    Lu
        1.30, 1.40, 1.10, 1.10, 1.20, 1.60, 1.40,               &
        ! Hf    Ta    W     Re    Os    Ir    Pt    Au    Hg
        1.30, 1.20, 1.00, 1.00, 0.90, 1.30, 1.20, 1.00, 1.00,   &
        ! Tl    Pb    Bi    Po    At    Rn
        1.20, 1.20, 1.10, 1.20, 1.10, 2.10,                     &
        ! Fr    Ra
        2.20, 1.80,                                             &
        ! Ac    Th    Pa    U     Np    Pu    Am    Cm
        1.70, 1.30, 1.40, 1.20, 1.20, 1.30, 1.40, 1.40,         &
        ! Bk    Cf    Es    Fm    Md    No    Lw
        1.70, 1.90, 1.90, 2.00, 2.00, 1.60, 2.00,               &
        ! other usless elements:
        2.0, 1.6, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,       &
        2.0, 2.0, 2.0, 2.0, 2.0, 0.1, 0.1 /)
   end type

   type(M4RadialGridParamType), save :: M4RadialGridParam


contains
   

   subroutine CreateM4RadialGrid(radPoints, weights, numRadPoints, mapParam)

      real(KREAL),              intent(out) :: radPoints (:)        
      real(KREAL),              intent(out) :: weights (:)          
      integer(KINT),            intent(in)  :: numRadPoints 
      real(KREAL),   optional,  intent(in)  :: mapParam

      real(KREAL)                   :: mapParamLoc
      real(KREAL),     allocatable  :: legPoints(:), legWeights(:)
      integer(KINT)                 :: i_pnt
   
      
      if ( present(mapParam) ) then
         mapParamLoc = mapParam
      else
         mapParamLoc = 1.0_KREAL
      end if

      allocate(legPoints(numRadPoints), legWeights(numRadPoints))
      
      call legpnt (numRadPoints, legPoints, legWeights)

      do i_pnt=1, numRadPoints
         radPoints(i_pnt) = M4Map(mapParamLoc,legPoints(i_pnt))
         weights(i_pnt)   = legWeights(i_pnt) * M4MapJacobian(mapParamLoc,legPoints(i_pnt)) * radPoints(i_pnt)**2
      end do

      deallocate(legPoints, legWeights)
   end subroutine


   subroutine fixM4MinimumRadius( numRadPoints, radPoints, radWeights, minimumRadius )
      ! Transform the radial grid from [0:infinity] to [minimumRadius:infinity]
      integer(KINT),    intent(in)      :: numRadPoints
      real(KREAL),      intent(inout)   :: radPoints(:)
      real(KREAL),      intent(inout)   :: radWeights(:)
      real(KREAL),      intent(in)      :: minimumRadius
      
      integer(KINT)         :: i_pnt

      ! minimumRadius-shift of the points (and weights adjustment)
      
      do i_pnt=1, numRadPoints
         radWeights(i_pnt) = (radWeights(i_pnt) * (radPoints(i_pnt) + minimumRadius)**2 )/radPoints(i_pnt)**2
         radPoints(i_pnt) = radPoints(i_pnt) + minimumRadius
      end do

      ! Add to the first point the weight of the range [0:minimumRadius]
      
      radWeights(1)=radWeights(1) + minimumRadius**3 
   end subroutine 
   

   real(KREAL) function M4Map(rm,x)
      real(KREAL), intent(in) :: rm, x
      M4Map = (rm/log(two)) * (one + x) * log(two/(one-x)) 
   end function
   
   real(KREAL) function M4MapJacobian(rm,x)
      real(KREAL), intent(in) :: rm, x
      M4MapJacobian = (rm/log(two)) * (log(two/(one-x)) + (one+x)/(one-x) )
   end function

end module
