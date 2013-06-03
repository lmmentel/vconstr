module LebedevAngGridModule

   use vartypes
   use LegendreGridModule

   implicit none

   private
   save
   
   public :: CreateLebedevPoints
   public :: MinAcceptableLebedevOrder

   real(KREAL) , parameter   :: pi      =   3.1415926535897931_KREAL
   real(KREAL) , parameter   :: fourpi  =   4.0_KREAL*pi
   real(KREAL) , parameter   :: zero    =   0.0_KREAL
   real(KREAL) , parameter   :: one     =   1.0_KREAL
   real(KREAL) , parameter   :: mone    =  -1.0_KREAL
   
   ! maxNumPointPerOper = 48 (number of octahedral oper) * 36 (max number of symmetry 
   ! unique points (corresponding to LOrder=59 ))
   ! If higher L-orders are implemented, this number should be changed as well.
   integer(KINT),   parameter :: maxNumPointPerOper = 1728_KINT 
   integer(KINT),   parameter :: Orders(15) = (/5_KINT, 7_KINT, 11_KINT, 13_KINT, 17_KINT, 19_KINT, 21_KINT, &
                                               23_KINT, 29_KINT, 31_KINT, 35_KINT, 41_KINT, 47_KINT, 53_KINT, 59_KINT/)

contains 

   subroutine CreateLebedevPoints( LOrder, nPoints, xyz, weights, symOperIn, octaOperIn, symToleranceIn, kvectors, &
                                   allowSpecialGrid )
      ! =======================================================================
      ! Purpose: Create a lebedev spherical grid of order LOrder on 
      !          the unit sphere (r=1)  
      !          ref: V. Lebedev, D. Laikov, Doklady Math. 59 (1999) 477–481 
      !           
      !          - possible LOrders: 5,7,11,13,17,19,21,23,29,31,35,41,47,53,59
      !          - symOperIn and octaOperIn are optional. If symOperIn is not 
      !            present, it will assume NOSYM and the full spherical grid 
      !            will be created.
      ! =======================================================================

      integer(KINT),                    intent(in)  :: LOrder 
      integer(KINT),                    intent(out) :: nPoints
      real(KREAL),                      intent(out) :: xyz(:,:)     ! size(3,numMaxPoints)
      real(KREAL),                      intent(out) :: weights(:)   ! size(numMaxPoints)
      real(KREAL),  optional, target,   intent(in)  :: symOperIn(:,:,:)
      real(KREAL),  optional, target,   intent(in)  :: octaOperIn(:,:,:)  
      real(KREAL),  optional,           intent(in)  :: symToleranceIn
      real(KREAL),  optional,           intent(in)  :: kvectors(:,:)     ! needed for periodic systems
      logical,      optional,           intent(in)  :: allowSpecialGrid  ! special grid with three-fold symmetry
 
      real(KREAL), allocatable  :: xyzSupp1(:,:), weightsSupp1(:)      
      real(KREAL), allocatable  :: xyzSupp2(:,:), weightsSupp2(:)      
      real(KREAL), allocatable  :: translSupp(:,:)
      real(KREAL)               :: matSupp(3,3)
      integer(KINT)             :: i_glob, i_oper, i_pnt, numSymOper, intSupp 
      real(KREAL),  pointer     :: octaOper(:,:,:)  => null()
      real(KREAL),  pointer     :: symOper(:,:,:)   => null()
      real(KREAL)               :: symTolerance
      real(KREAL) :: det, axis(3), axisOld(3), angle, oneoneoneDirection(3), rotAxis(3)
      real(KREAL) :: tempRotateMatrix1(3,3), tempRotateMatrix2(3,3), tempRotateMatrix3(3,3)
      logical :: specialGridNeeded
      real(KREAL) :: axis1(3), axis2(3), axisOrig(3)

      if (present(symOperIn)) then
         symOper => symOperIn
         allocate( translSupp(3,size(symOperIn(1,1,:))) )
      else 
         allocate( symOper(3,3,1) )
         allocate( translSupp(3,1) )
         symOper = reshape((/one ,  zero,  zero,  zero,  one ,  zero,  zero,  zero,  one /), shape(symOper))
      end if
      translSupp = 0.0_KREAL
      numSymOper = size(symOper(1,1,:))
     
      if (present(octaOperIn)) then
         octaOper => octaOperIn
      else 
         allocate(octaOper(3,3,48))
         call fillOctahedralOperator( octaOper )
      end if

      if (present(symToleranceIn)) then
         symTolerance = symToleranceIn
      else 
         symTolerance = 1.0E-7_KREAL
      end if

      
      allocate( xyzSupp1(3,maxNumPointPerOper*numSymOper), weightsSupp1(maxNumPointPerOper*numSymOper), &
                xyzSupp2(3,maxNumPointPerOper*numSymOper), weightsSupp2(maxNumPointPerOper*numSymOper))

      ! ----------------------------------------------
      ! Create the symmetry (octahedral) unique points
      ! ----------------------------------------------
 
      call CreateLebedevUniquePoints(nPoints, xyzSupp1,  weightsSupp1, LOrder, .false., octaOper) 

      ! --------------------------------------------------------
      ! Generate the full octahedral grid from the "unique part"
      ! --------------------------------------------------------
      i_glob=0_KINT
      do i_oper=1, 48 ! number of octahedral operators
         do i_pnt=1, nPoints
            i_glob=i_glob + 1_KINT
            xyzSupp2(:,i_glob)   = matmul(octaOper(:,:,i_oper),xyzSupp1(:,i_pnt))
            weightsSupp2(i_glob) = weightsSupp1(i_pnt)
         end do
      end do
      nPoints = i_glob
      xyzSupp1(:,1:nPoints)     = xyzSupp2(:,1:nPoints)
      weightsSupp1(1:nPoints)   = weightsSupp2(1:nPoints)

      ! ========================================
      ! Special grid for "three-fold symmetries"
      ! ========================================

      if (present (allowSpecialGrid)) then
         if (allowSpecialGrid .and. present(symOperIn)) then
            stop
         end if
      end if

      ! --------------------------------------------------------
      ! Multiply the octahedral grid for all the local operators
      ! --------------------------------------------------------
      i_glob=0_KINT
      do i_oper=1, numSymOper 
         do i_pnt=1, nPoints
            i_glob=i_glob + 1_KINT
            xyzSupp2(:,i_glob)   = matmul(symOper(:,:,i_oper),xyzSupp1(:,i_pnt))
            weightsSupp2(i_glob) = weightsSupp1(i_pnt)/real(numSymOper)
         end do
      end do
      nPoints = i_glob
      xyzSupp1(:,1:nPoints)     = xyzSupp2(:,1:nPoints)
      weightsSupp1(1:nPoints)   = weightsSupp2(1:nPoints)



      ! -------------------------
      ! Merge overlapping points:
      ! -------------------------
      i_glob=0
      do i_pnt=1, nPoints
         intSupp=findOverlappingPoints(xyzSupp1(:,i_pnt), xyzSupp2(:,:), i_glob, symTolerance)
         ! if unique point keep, if non unique merge
        if ( intSupp == -1_KINT) then
            i_glob=i_glob+1
            xyzSupp2(:,i_glob)   = xyzSupp1(:,i_pnt)
            weightsSupp2(i_glob) = weightsSupp1(i_pnt)
         else
            weightsSupp2(intSupp)   = weightsSupp2(intSupp) + weightsSupp1(i_pnt)
         end if
      end do
      nPoints = i_glob
      xyzSupp1(:,1:nPoints)     = xyzSupp2(:,1:nPoints)
      weightsSupp1(1:nPoints)   = weightsSupp2(1:nPoints)

      ! --------------------------------------------------
      ! Reduce the points according to the local symmetry:
      ! --------------------------------------------------
      translSupp    = 0.0_KREAL
      matSupp       = 0.0_KREAL
      intSupp       = 0_KINT
      if (present(kvectors)) then 
         call pntred (intSupp, matSupp, kvectors, numsymOper, symOper, translSupp, nPoints, &
                      xyzSupp1, weightsSupp1)
      else
         call pntred (intSupp, matSupp, matSupp, numsymOper, symOper, translSupp, nPoints, &
                      xyzSupp1, weightsSupp1)
      end if


      if (nPoints > size(xyz(1,:))) then
         write(*,*) "LebedevAngGrid: xyz and weights arrays too small"
         stop
      end if

      xyz(:,1:nPoints) = xyzSupp1(:,1:nPoints)
      weights(1:nPoints) = weightsSupp1(1:nPoints)


      ! Octahedral operators deallocation:
      if ( .not. present(octaOperIn)) then
         deallocate(octaOper)
      end if
      if ( .not. present(symOperIn)) then
         deallocate(symOper)
      end if

      deallocate( xyzSupp1, weightsSupp1, xyzSupp2, weightsSupp2, translSupp )
   end subroutine 



   subroutine CreateLebedevUniquePoints(nPoints, xyzUnique,  wUnique, LOrder, weightsNorm, octaOper)
      ! =============================================================================================
      ! purpose: This elegant subroutine generates the symmetry unique lebedev points and weights. 
      !          If weightsNorm=true the weights are the 'standard' ones (see original lebedev article), 
      !          if it's false the weights are divided by the number of the point's symmetric-copies. 
      !          The reason in that CalcAngPoints will use all the operators even if a point il laying 
      !          on a symmetry plane (the overlapping points will be merged and the weights combined). 
      ! =============================================================================================

      integer(KINT), parameter   :: KR=KREAL        !for "formatting" reasons i need a shorter name for KREAL

      integer(KINT), intent(in)  :: LOrder
      logical,       intent(in)  :: weightsNorm
      real(KREAL),   intent(in)  :: octaOper(:,:,:) ! octaOper(3,3,numOctaOper) (numOctaOper=48)
      integer(KINT), intent(out) :: nPoints
      real(KREAL),   intent(out) :: xyzUnique(:,:)  !3, maxNumPoints
      real(KREAL),   intent(out) :: wUnique(:)
      
      integer(KINT)              :: numOctaOper, i_pnt
      real(KREAL), allocatable   :: xyzw(:,:)   !this guy will contain the "tabulated" lebedev weights and points

      
      numOctaOper = size(octaOper(1,1,:))
      
      select case(LOrder)    
      case(5_KINT)
         nPoints=2_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.8377580409572781E+0_KR /)
         xyzw(:,2)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.9424777960769379E+0_KR /)

      case(7_KINT)
         nPoints=3_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.5983986006837702E+0_KR /)
         xyzw(:,2)=(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.4787188805470161E+0_KR /)
         xyzw(:,3)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.4039190554615448E+0_KR /)

      case(11_KINT)
         nPoints=4_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.1595729601823387E+0_KR /)
         xyzw(:,2)=(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.2836852625463799E+0_KR /)
         xyzw(:,3)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.2650718801466388E+0_KR /)
         xyzw(:,4)=(/ 0.3015113445777636E+0_KR, 0.3015113445777636E+0_KR, 0.9045340337332909E+0_KR, 0.2535056108973113E+0_KR /)
   
      case(13_KINT)
         nPoints=5_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.1450663274384897E+0_KR /)
         xyzw(:,2)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.1500915881570819E+0_KR /)
         xyzw(:,3)=(/ 0.3696028464541502E+0_KR, 0.3696028464541502E+0_KR, 0.8525183117012676E+0_KR, 0.1396193607909271E+0_KR /)
         xyzw(:,4)=(/ 0.6943540066026664E+0_KR, 0.6943540066026664E+0_KR, 0.1890635528853950E+0_KR, 0.1492445168690702E+0_KR /)
         xyzw(:,5)=(/ 0.3742430390903412E+0_KR, 0.9273306571511725E+0_KR, 0.0000000000000000E+0_KR, 0.1484377866929852E+0_KR /)

      case(17_KINT)
         nPoints=6_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.4810746585139659E-1_KR /)
         xyzw(:,2)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.1230717352816702E+0_KR /)
         xyzw(:,3)=(/ 0.1851156353447362E+0_KR, 0.1851156353447362E+0_KR, 0.9651240350865941E+0_KR, 0.1031917340883304E+0_KR /)
         xyzw(:,4)=(/ 0.6904210483822922E+0_KR, 0.6904210483822922E+0_KR, 0.2159572918458484E+0_KR, 0.1249450968725133E+0_KR /)
         xyzw(:,5)=(/ 0.3956894730559419E+0_KR, 0.3956894730559419E+0_KR, 0.8287699812525923E+0_KR, 0.1205802490285279E+0_KR /)
         xyzw(:,6)=(/ 0.4783690288121502E+0_KR, 0.8781589106040661E+0_KR, 0.0000000000000000E+0_KR, 0.1218309173855214E+0_KR /)

      case(19_KINT)
         nPoints=7_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.7535190013117138E-2_KR /)
         xyzw(:,2)=(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.9265184700375431E-1_KR /)
         xyzw(:,3)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.9061000833610514E-1_KR /)
         xyzw(:,4)=(/ 0.6764410400114264E+0_KR, 0.6764410400114264E+0_KR, 0.2912988822095268E+0_KR, 0.8942676055004592E-1_KR /)
         xyzw(:,5)=(/ 0.4174961227965453E+0_KR, 0.4174961227965453E+0_KR, 0.8070898183595826E+0_KR, 0.8487112439121475E-1_KR /)
         xyzw(:,6)=(/ 0.1574676672039082E+0_KR, 0.1574676672039082E+0_KR, 0.9748886436771732E+0_KR, 0.9518264418191037E-1_KR /)
         xyzw(:,7)=(/ 0.1403553811713183E+0_KR, 0.4493328323269557E+0_KR, 0.8822700112603227E+0_KR, 0.8785259467896815E-1_KR /)
        
      case(21_KINT)
         nPoints=8_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.6967855090540039E-1_KR /)
         xyzw(:,2)=(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.7629461771935279E-1_KR /)
         xyzw(:,3)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.8021962308552601E-1_KR /)
         xyzw(:,4)=(/ 0.2551252621114134E+0_KR, 0.2551252621114134E+0_KR, 0.9326425903126906E+0_KR, 0.6513636946550791E-1_KR /)
         xyzw(:,5)=(/ 0.6743601460362766E+0_KR, 0.6743601460362766E+0_KR, 0.3007935951377015E+0_KR, 0.7939343745253054E-1_KR /)
         xyzw(:,6)=(/ 0.4318910696719410E+0_KR, 0.4318910696719410E+0_KR, 0.7917955593934921E+0_KR, 0.7793248373075363E-1_KR /)
         xyzw(:,7)=(/ 0.2613931360335988E+0_KR, 0.9652324219764484E+0_KR, 0.0000000000000000E+0_KR, 0.6882781368562169E-1_KR /)
         xyzw(:,8)=(/ 0.4990453161796037E+0_KR, 0.1446630744325115E+0_KR, 0.8544158046846588E+0_KR, 0.7500092515800830E-1_KR /)

      case(23_KINT)
         nPoints=9_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1)=(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.2239755062103847E-1_KR /)
         xyzw(:,2)=(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.7184075893484736E-1_KR /)
         xyzw(:,3)=(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.7003719860124849E-1_KR /)
         xyzw(:,4)=(/ 0.6712973442695226E+0_KR, 0.6712973442695226E+0_KR, 0.3141969941825863E+0_KR, 0.7048105416807013E-1_KR /)
         xyzw(:,5)=(/ 0.2892465627575439E+0_KR, 0.2892465627575439E+0_KR, 0.9125090968674737E+0_KR, 0.6482032680351046E-1_KR /)
         xyzw(:,6)=(/ 0.4446933178717437E+0_KR, 0.4446933178717437E+0_KR, 0.7774932193147671E+0_KR, 0.6935092759371100E-1_KR /)
         xyzw(:,7)=(/ 0.1299335447650067E+0_KR, 0.1299335447650067E+0_KR, 0.9829723027072532E+0_KR, 0.5160728216651316E-1_KR /)
         xyzw(:,8)=(/ 0.3457702197611283E+0_KR, 0.9383192181375916E+0_KR, 0.0000000000000000E+0_KR, 0.6348336993464156E-1_KR /)
         xyzw(:,9)=(/ 0.1590417105383530E+0_KR, 0.8360360154824589E+0_KR, 0.5251185724436420E+0_KR, 0.6949515747104322E-1_KR /)

      case(29_KINT)
         nPoints=12_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1) =(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.1073910939755579E-1_KR /)
         xyzw(:,2) =(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.4522786682091873E-1_KR /)
         xyzw(:,3) =(/ 0.3515640345570105E+0_KR, 0.3515640345570105E+0_KR, 0.8676436245440834E+0_KR, 0.4335131988095388E-1_KR /)
         xyzw(:,4) =(/ 0.6566329410219612E+0_KR, 0.6566329410219612E+0_KR, 0.3710341783848209E+0_KR, 0.4529953680846059E-1_KR /)
         xyzw(:,5) =(/ 0.4729054132581005E+0_KR, 0.4729054132581005E+0_KR, 0.7434520429875557E+0_KR, 0.4494651051683867E-1_KR /)
         xyzw(:,6) =(/ 0.9618308522614784E-1_KR, 0.9618308522614784E-1_KR, 0.9907056213794081E+0_KR, 0.2955737808697618E-1_KR /)
         xyzw(:,7) =(/ 0.2219645236294178E+0_KR, 0.2219645236294178E+0_KR, 0.9494543172264431E+0_KR, 0.3906825715891940E-1_KR /)
         xyzw(:,8) =(/ 0.7011766416089545E+0_KR, 0.7011766416089545E+0_KR, 0.1292386727105144E+0_KR, 0.4586782837866035E-1_KR /)
         xyzw(:,9) =(/ 0.2644152887060663E+0_KR, 0.9644089148792060E+0_KR, 0.0000000000000000E+0_KR, 0.3747725210708425E-1_KR /)
         xyzw(:,10)=(/ 0.5718955891878961E+0_KR, 0.8203264198277593E+0_KR, 0.0000000000000000E+0_KR, 0.4524925035017432E-1_KR /)
         xyzw(:,11)=(/ 0.2510034751770465E+0_KR, 0.8000727494073951E+0_KR, 0.5448677372580774E+0_KR, 0.4488130226921316E-1_KR /)
         xyzw(:,12)=(/ 0.1233548532583327E+0_KR, 0.4127724083168531E+0_KR, 0.9024425295330004E+0_KR, 0.4262905240772150E-1_KR /)

      case(31_KINT)
         nPoints=13_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1) =(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.3778452231568862E-1_KR /)
         xyzw(:,2) =(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.3833531885829462E-1_KR /)
         xyzw(:,3) =(/ 0.7068965463912316E+0_KR, 0.7068965463912316E+0_KR, 0.2438330166935553E-1_KR, 0.2037140121187405E-1_KR /)
         xyzw(:,4) =(/ 0.4794682625712025E+0_KR, 0.4794682625712025E+0_KR, 0.7349968505877456E+0_KR, 0.3777075881540511E-1_KR /)
         xyzw(:,5) =(/ 0.1927533154878019E+0_KR, 0.1927533154878019E+0_KR, 0.9621290551360144E+0_KR, 0.3758592063240899E-1_KR /)
         xyzw(:,6) =(/ 0.6930357961327123E+0_KR, 0.6930357961327123E+0_KR, 0.1985013112233654E+0_KR, 0.3747506154911825E-1_KR /)
         xyzw(:,7) =(/ 0.3608302115520091E+0_KR, 0.3608302115520091E+0_KR, 0.8600018121275470E+0_KR, 0.3420018485683569E-1_KR /)
         xyzw(:,8) =(/ 0.6498486161496169E+0_KR, 0.6498486161496169E+0_KR, 0.3941998886058389E+0_KR, 0.3812025862193428E-1_KR /)
         xyzw(:,9) =(/ 0.1932945013230339E+0_KR, 0.9811407828432572E+0_KR, 0.0000000000000000E+0_KR, 0.3779900890017291E-1_KR /)
         xyzw(:,10)=(/ 0.3800494919899303E+0_KR, 0.9249661526986790E+0_KR, 0.0000000000000000E+0_KR, 0.3621583529945751E-1_KR /)
         xyzw(:,11)=(/ 0.2899558825499574E+0_KR, 0.7934537856582315E+0_KR, 0.5351230477182762E+0_KR, 0.3717581834486352E-1_KR /)
         xyzw(:,12)=(/ 0.9684121455103957E-1_KR, 0.8280801506686862E+0_KR, 0.5521820743493993E+0_KR, 0.3815175284444799E-1_KR /)
         xyzw(:,13)=(/ 0.1833434647041659E+0_KR, 0.9074658265305127E+0_KR, 0.3780091898744867E+0_KR, 0.3559031656705768E-1_KR /)

      case(35_KINT)
         nPoints=16_KINT
         allocate(xyzw(4,nPoints))
         xyzw(1:3,1) = (/1.0_KREAL/sqrt(3.0_KREAL), 1.0_KREAL/sqrt(3.0_KREAL), 1.0_KREAL/sqrt(3.0_KREAL) /)
         xyzw(1:3,2) = (/1.0_KREAL,                 0.0_KREAL,                 0.0_KREAL                 /)
         xyzw(1:3,3) = (/0.70710678118654752_KREAL, 0.70710678118654752_KREAL, 0.0_KREAL                 /)
         xyzw(1:3,4) = (/0.69093463105113458_KREAL, 0.69093463105113458_KREAL, 0.21264682275657227_KREAL /)
         xyzw(1:3,5) = (/0.64566647095194987_KREAL, 0.64566647095194987_KREAL, 0.40771266423415123_KREAL /)
         xyzw(1:3,6) = (/0.49143426555639500_KREAL, 0.49143426555639500_KREAL, 0.71901649861049314_KREAL /)
         xyzw(1:3,7) = (/0.39272598223217649_KREAL, 0.39272598223217649_KREAL, 0.83158439485090411_KREAL /)
         xyzw(1:3,8) = (/0.28612891787658218_KREAL, 0.28612891787658218_KREAL, 0.91447279057911405_KREAL /)
         xyzw(1:3,9) = (/0.17748365242374568_KREAL, 0.17748365242374568_KREAL, 0.96798714156989403_KREAL /)
         xyzw(1:3,10)= (/0.07568095866244468_KREAL, 0.07568095866244468_KREAL, 0.99425589512552888_KREAL /)
         xyzw(1:3,11)= (/0.97764280892098723_KREAL, 0.21027253307334757_KREAL, 0.0_KREAL                 /)
         xyzw(1:3,12)= (/0.88181328936054412_KREAL, 0.47159868819488597_KREAL, 0.0_KREAL                 /)
         xyzw(1:3,13)= (/0.09921769971362576_KREAL, 0.33443631695748371_KREAL, 0.93718098463607886_KREAL /)
         xyzw(1:3,14)= (/0.20548237125466495_KREAL, 0.45023303874296735_KREAL, 0.86894603165434486_KREAL /)
         xyzw(1:3,15)= (/0.10680182513533723_KREAL, 0.59051570309804130_KREAL, 0.79992785583600399_KREAL /)
         xyzw(1:3,16)= (/0.31042840327515130_KREAL, 0.55501523681448068_KREAL, 0.77174626228042463_KREAL /)
         xyzw(4,1)   = fourpi*0.25123173709441058e-2_KREAL
         xyzw(4,2)   = fourpi*0.52659765761428065e-3_KREAL
         xyzw(4,3)   = fourpi*0.25482199909403521e-2_KREAL
         xyzw(4,4)   = fourpi*0.25304038224001323e-2_KREAL
         xyzw(4,5)   = fourpi*0.25132671684706878e-2_KREAL
         xyzw(4,6)   = fourpi*0.25017251210647733e-2_KREAL
         xyzw(4,7)   = fourpi*0.24453733047996786e-2_KREAL
         xyzw(4,8)   = fourpi*0.23026944325620758e-2_KREAL
         xyzw(4,9)   = fourpi*0.20142782609526094e-2_KREAL
         xyzw(4,10)  = fourpi*0.14624950815475142e-2_KREAL
         xyzw(4,11)  = fourpi*0.19109513147305082e-2_KREAL
         xyzw(4,12)  = fourpi*0.24174423575419847e-2_KREAL
         xyzw(4,13)  = fourpi*0.22366077071364263e-2_KREAL
         xyzw(4,14)  = fourpi*0.24169300107381179e-2_KREAL
         xyzw(4,15)  = fourpi*0.25122368647336706e-2_KREAL
         xyzw(4,16)  = fourpi*0.24966440519292456e-2_KREAL
    
      case(41_KINT)
         nPoints=20_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1) =(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.3889444129321297E-2_KR /)
         xyzw(:,2) =(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.2327768981109099E-1_KR /)
         xyzw(:,3) =(/ 0.7040954938227469E+0_KR, 0.7040954938227469E+0_KR, 0.9219040707689825E-1_KR, 0.2352161488565241E-1_KR /)
         xyzw(:,4) =(/ 0.6807744066455244E+0_KR, 0.6807744066455244E+0_KR, 0.2703560883591648E+0_KR, 0.2335852785125307E-1_KR /)
         xyzw(:,5) =(/ 0.6372546939258752E+0_KR, 0.6372546939258752E+0_KR, 0.4333738687771544E+0_KR, 0.2327328064484758E-1_KR /)
         xyzw(:,6) =(/ 0.5044419707800358E+0_KR, 0.5044419707800358E+0_KR, 0.7007685753735730E+0_KR, 0.2320651712444717E-1_KR /)
         xyzw(:,7) =(/ 0.4215761784010967E+0_KR, 0.4215761784010967E+0_KR, 0.8028368773352738E+0_KR, 0.2285159031614609E-1_KR /)
         xyzw(:,8) =(/ 0.3317920736472123E+0_KR, 0.3317920736472123E+0_KR, 0.8830787279341326E+0_KR, 0.2198567789717927E-1_KR /)
         xyzw(:,9) =(/ 0.2384736701421887E+0_KR, 0.2384736701421887E+0_KR, 0.9414141582204025E+0_KR, 0.2032246835488661E-1_KR /)
         xyzw(:,10)=(/ 0.1459036449157763E+0_KR, 0.1459036449157763E+0_KR, 0.9784805837626939E+0_KR, 0.1740112129664928E-1_KR /)
         xyzw(:,11)=(/ 0.6095034115507196E-1_KR, 0.6095034115507196E-1_KR, 0.9962781297540164E+0_KR, 0.1227022042213690E-1_KR /)
         xyzw(:,12)=(/ 0.6116843442009876E+0_KR, 0.7911019296269020E+0_KR, 0.0000000000000000E+0_KR, 0.2333777588926988E-1_KR /)
         xyzw(:,13)=(/ 0.3964755348199858E+0_KR, 0.9180452877114540E+0_KR, 0.0000000000000000E+0_KR, 0.2142759707326609E-1_KR /)
         xyzw(:,14)=(/ 0.1724782009907724E+0_KR, 0.9850133350280019E+0_KR, 0.0000000000000000E+0_KR, 0.1634032422273241E-1_KR /)
         xyzw(:,15)=(/ 0.5610263808622060E+0_KR, 0.3518280927733519E+0_KR, 0.7493106119041159E+0_KR, 0.2315814309130472E-1_KR /)
         xyzw(:,16)=(/ 0.4742392842551980E+0_KR, 0.2634716655937950E+0_KR, 0.8400474883590504E+0_KR, 0.2265288026067282E-1_KR /)
         xyzw(:,17)=(/ 0.5984126497885380E+0_KR, 0.1816640840360209E+0_KR, 0.7803207424799203E+0_KR, 0.2324565639630277E-1_KR /)
         xyzw(:,18)=(/ 0.3791035407695563E+0_KR, 0.1720795225656878E+0_KR, 0.9092134750923736E+0_KR, 0.2153755923392349E-1_KR /)
         xyzw(:,19)=(/ 0.2778673190586244E+0_KR, 0.8213021581932511E-1_KR, 0.9571020743100725E+0_KR, 0.1954339052477729E-1_KR /)
         xyzw(:,20)=(/ 0.5033564271075117E+0_KR, 0.8999205842074876E-1_KR, 0.8593798558907212E+0_KR, 0.2264760481825463E-1_KR /)

      case(47_KINT)
         nPoints=25_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1) =(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.2755732301591147E-2_KR /)
         xyzw(:,2) =(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.1805075719815613E-1_KR /)
         xyzw(:,3) =(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.1786862935843413E-1_KR /)
         xyzw(:,4) =(/ 0.5087204410502360E-1_KR, 0.5087204410502360E-1_KR, 0.9974086776528230E+0_KR, 0.8542773952204923E-2_KR /)
         xyzw(:,5) =(/ 0.1228198790178831E+0_KR, 0.1228198790178831E+0_KR, 0.9847997535723012E+0_KR, 0.1245727470691386E-1_KR /)
         xyzw(:,6) =(/ 0.2026890814408786E+0_KR, 0.2026890814408786E+0_KR, 0.9580366759833914E+0_KR, 0.1483092903445044E-1_KR /)
         xyzw(:,7) =(/ 0.2847745156464294E+0_KR, 0.2847745156464294E+0_KR, 0.9153179504831548E+0_KR, 0.1629355113817948E-1_KR /)
         xyzw(:,8) =(/ 0.3656719078978026E+0_KR, 0.3656719078978026E+0_KR, 0.8559019286978865E+0_KR, 0.1716404656862801E-1_KR /)
         xyzw(:,9) =(/ 0.4428264886713469E+0_KR, 0.4428264886713469E+0_KR, 0.7796213195276351E+0_KR, 0.1763047477532942E-1_KR /)
         xyzw(:,10)=(/ 0.5140619627249735E+0_KR, 0.5140619627249735E+0_KR, 0.6866444472641543E+0_KR, 0.1782722592255887E-1_KR /)
         xyzw(:,11)=(/ 0.6306401219166803E+0_KR, 0.6306401219166803E+0_KR, 0.4523119203136584E+0_KR, 0.1786154692073831E-1_KR /)
         xyzw(:,12)=(/ 0.6716883332022612E+0_KR, 0.6716883332022612E+0_KR, 0.3125213050016533E+0_KR, 0.1789446746456066E-1_KR /)
         xyzw(:,13)=(/ 0.6979792685336881E+0_KR, 0.6979792685336881E+0_KR, 0.1601558034988290E+0_KR, 0.1798943864849984E-1_KR /)
         xyzw(:,14)=(/ 0.1446865674195309E+0_KR, 0.9894775374955985E+0_KR, 0.0000000000000000E+0_KR, 0.1162942390613896E-1_KR /)
         xyzw(:,15)=(/ 0.3390263475411216E+0_KR, 0.9407768787937587E+0_KR, 0.0000000000000000E+0_KR, 0.1571097913473697E-1_KR /)
         xyzw(:,16)=(/ 0.5335804651263506E+0_KR, 0.8457493051936533E+0_KR, 0.0000000000000000E+0_KR, 0.1752211795927858E-1_KR /)
         xyzw(:,17)=(/ 0.6944024393349413E-1_KR, 0.2355187894242326E+0_KR, 0.9693858634984321E+0_KR, 0.1416341927904775E-1_KR /)
         xyzw(:,18)=(/ 0.2269004109529460E+0_KR, 0.4102182474045730E+0_KR, 0.8833103605221128E+0_KR, 0.1691124051527118E-1_KR /)
         xyzw(:,19)=(/ 0.8025574607775339E-1_KR, 0.6214302417481605E+0_KR, 0.7793481057026609E+0_KR, 0.1790654133178910E-1_KR /)
         xyzw(:,20)=(/ 0.1467999527896572E+0_KR, 0.3245284345717394E+0_KR, 0.9344148270524022E+0_KR, 0.1585276984465826E-1_KR /)
         xyzw(:,21)=(/ 0.1571507769824727E+0_KR, 0.5224482189696630E+0_KR, 0.8380641334583124E+0_KR, 0.1749926303261150E-1_KR /)
         xyzw(:,22)=(/ 0.2365702993157246E+0_KR, 0.6017546634089558E+0_KR, 0.7628406246046698E+0_KR, 0.1782868505766069E-1_KR /)
         xyzw(:,23)=(/ 0.7714815866765733E-1_KR, 0.4346575516141163E+0_KR, 0.8972853361328333E+0_KR, 0.1681841177508118E-1_KR /)
         xyzw(:,24)=(/ 0.3062936666210730E+0_KR, 0.4908826589037616E+0_KR, 0.8156092232039754E+0_KR, 0.1751376156594036E-1_KR /)
         xyzw(:,25)=(/ 0.3822477379524787E+0_KR, 0.5648768149099500E+0_KR, 0.7313007936597657E+0_KR, 0.1779290960066995E-1_KR /)

      case(53_KINT)
         nPoints=30_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1) =(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.1807413785064742E-2_KR /)
         xyzw(:,2) =(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.1414687180198969E-1_KR /)
         xyzw(:,3) =(/ 0.4292963545341347E-1_KR, 0.4292963545341347E-1_KR, 0.9981553450238465E+0_KR, 0.6217877052165790E-2_KR /)
         xyzw(:,4) =(/ 0.1051426854086404E+0_KR, 0.1051426854086404E+0_KR, 0.9888832243546856E+0_KR, 0.9246323068805974E-2_KR /)
         xyzw(:,5) =(/ 0.1750024867623087E+0_KR, 0.1750024867623087E+0_KR, 0.9688902204347074E+0_KR, 0.1117041368444565E-1_KR /)
         xyzw(:,6) =(/ 0.2477653379650257E+0_KR, 0.2477653379650257E+0_KR, 0.9366027304071631E+0_KR, 0.1242606437075843E-1_KR /)
         xyzw(:,7) =(/ 0.3206567123955957E+0_KR, 0.3206567123955957E+0_KR, 0.8912679426476061E+0_KR, 0.1323615416834776E-1_KR /)
         xyzw(:,8) =(/ 0.3916520749849983E+0_KR, 0.3916520749849983E+0_KR, 0.8325967237023519E+0_KR, 0.1373226348846247E-1_KR /)
         xyzw(:,9) =(/ 0.4590825874187624E+0_KR, 0.4590825874187624E+0_KR, 0.7605829053152514E+0_KR, 0.1400383013450966E-1_KR /)
         xyzw(:,10)=(/ 0.5214563888415861E+0_KR, 0.5214563888415861E+0_KR, 0.6754009691084143E+0_KR, 0.1412114215519805E-1_KR /)
         xyzw(:,11)=(/ 0.6253170244654199E+0_KR, 0.6253170244654199E+0_KR, 0.4668589056957432E+0_KR, 0.1414017439086521E-1_KR /)
         xyzw(:,12)=(/ 0.6637926744523170E+0_KR, 0.6637926744523170E+0_KR, 0.3446136542374379E+0_KR, 0.1415165938221183E-1_KR /)
         xyzw(:,13)=(/ 0.6910410398498301E+0_KR, 0.6910410398498301E+0_KR, 0.2119541518501843E+0_KR, 0.1420360447706885E-1_KR /)
         xyzw(:,14)=(/ 0.7052907007457760E+0_KR, 0.7052907007457760E+0_KR, 0.7162440144995555E-1_KR, 0.1426266143312456E-1_KR /)
         xyzw(:,15)=(/ 0.1236686762657990E+0_KR, 0.9923235654314901E+0_KR, 0.0000000000000000E+0_KR, 0.8574497021019509E-2_KR /)
         xyzw(:,16)=(/ 0.2940777114468387E+0_KR, 0.9557815124965484E+0_KR, 0.0000000000000000E+0_KR, 0.1188044552909464E-1_KR /)
         xyzw(:,17)=(/ 0.4697753849207649E+0_KR, 0.8827859807011816E+0_KR, 0.0000000000000000E+0_KR, 0.1350168526987324E-1_KR /)
         xyzw(:,18)=(/ 0.6334563241139567E+0_KR, 0.7737784472573748E+0_KR, 0.0000000000000000E+0_KR, 0.1419120342265561E-1_KR /)
         xyzw(:,19)=(/ 0.5974048614181342E-1_KR, 0.2029128752777523E+0_KR, 0.9773727228453100E+0_KR, 0.1060210174688767E-1_KR /)
         xyzw(:,20)=(/ 0.1375760408473636E+0_KR, 0.4602621942484054E+0_KR, 0.8770584618658027E+0_KR, 0.1351206188837047E-1_KR /)
         xyzw(:,21)=(/ 0.3391016526336286E+0_KR, 0.5030673999662036E+0_KR, 0.7949422999642084E+0_KR, 0.1393079241308106E-1_KR /)
         xyzw(:,22)=(/ 0.1271675191439820E+0_KR, 0.2817606422442134E+0_KR, 0.9510201693743899E+0_KR, 0.1202158743917833E-1_KR /)
         xyzw(:,23)=(/ 0.2693120740413512E+0_KR, 0.4331561291720157E+0_KR, 0.8601434616017620E+0_KR, 0.1358001491783288E-1_KR /)
         xyzw(:,24)=(/ 0.1419786452601918E+0_KR, 0.6256167358580814E+0_KR, 0.7671021862205583E+0_KR, 0.1415975035780934E-1_KR /)
         xyzw(:,25)=(/ 0.6709284600738255E-1_KR, 0.3798395216859157E+0_KR, 0.9226161107308090E+0_KR, 0.1284997745583855E-1_KR /)
         xyzw(:,26)=(/ 0.7057738183256172E-1_KR, 0.5517505421423520E+0_KR, 0.8310175524134743E+0_KR, 0.1393560572068188E-1_KR /)
         xyzw(:,27)=(/ 0.2783888477882155E+0_KR, 0.6029619156159187E+0_KR, 0.7476206108340857E+0_KR, 0.1410940347341234E-1_KR /)
         xyzw(:,28)=(/ 0.1979578938917407E+0_KR, 0.3589606329589096E+0_KR, 0.9121183784091215E+0_KR, 0.1297354423382698E-1_KR /)
         xyzw(:,29)=(/ 0.2087307061103274E+0_KR, 0.5348666438135476E+0_KR, 0.8187485362810218E+0_KR, 0.1391410610029917E-1_KR /)
         xyzw(:,30)=(/ 0.4055122137872836E+0_KR, 0.5674997546074373E+0_KR, 0.7165918454670238E+0_KR, 0.1409670383749578E-1_KR /)

      case(59_KINT)
         nPoints=36_KINT
         allocate(xyzw(4,nPoints))
         xyzw(:,1) =(/ 0.1000000000000000E+1_KR, 0.0000000000000000E+0_KR, 0.0000000000000000E+0_KR, 0.1388821750423976E-2_KR /)
         xyzw(:,2) =(/ 0.0000000000000000E+0_KR, 0.7071067811865475E+0_KR, 0.7071067811865475E+0_KR, 0.1156763661782805E-1_KR /)
         xyzw(:,3) =(/ 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.5773502691896258E+0_KR, 0.1147706707566113E-1_KR /)
         xyzw(:,4) =(/ 0.3712636449657089E-1_KR, 0.3712636449657089E-1_KR, 0.9986206817999193E+0_KR, 0.4637520929383973E-2_KR /)
         xyzw(:,5) =(/ 0.9140060412262223E-1_KR, 0.9140060412262223E-1_KR, 0.9916107397220139E+0_KR, 0.7042182692930802E-2_KR /)
         xyzw(:,6) =(/ 0.1531077852469906E+0_KR, 0.1531077852469906E+0_KR, 0.9762766063946851E+0_KR, 0.8627187438744667E-2_KR /)
         xyzw(:,7) =(/ 0.2180928891660612E+0_KR, 0.2180928891660612E+0_KR, 0.9512470674805785E+0_KR, 0.9701663550202072E-2_KR /)
         xyzw(:,8) =(/ 0.2839874532200175E+0_KR, 0.2839874532200175E+0_KR, 0.9158068862086683E+0_KR, 0.1043203031916077E-1_KR /)
         xyzw(:,9) =(/ 0.3491177600963764E+0_KR, 0.3491177600963764E+0_KR, 0.8696169151819541E+0_KR, 0.1091601979985500E-1_KR /)
         xyzw(:,10)=(/ 0.4121431461444309E+0_KR, 0.4121431461444309E+0_KR, 0.8125737222999156E+0_KR, 0.1121809491106090E-1_KR /)
         xyzw(:,11)=(/ 0.4718993627149127E+0_KR, 0.4718993627149127E+0_KR, 0.7447294696321065E+0_KR, 0.1138616251879345E-1_KR /)
         xyzw(:,12)=(/ 0.5273145452842337E+0_KR, 0.5273145452842337E+0_KR, 0.6662422537361044E+0_KR, 0.1146025009059901E-1_KR /)
         xyzw(:,13)=(/ 0.6209475332444019E+0_KR, 0.6209475332444019E+0_KR, 0.4783809380769523E+0_KR, 0.1147148804964644E-1_KR /)
         xyzw(:,14)=(/ 0.6569722711857291E+0_KR, 0.6569722711857291E+0_KR, 0.3698308664594258E+0_KR, 0.1147399478559670E-1_KR /)
         xyzw(:,15)=(/ 0.6841788309070143E+0_KR, 0.6841788309070143E+0_KR, 0.2525839557007183E+0_KR, 0.1150184041631593E-1_KR /)
         xyzw(:,16)=(/ 0.7012604330123631E+0_KR, 0.7012604330123631E+0_KR, 0.1283261866597231E+0_KR, 0.1154527292189332E-1_KR /)
         xyzw(:,17)=(/ 0.1072382215478166E+0_KR, 0.9942333548213224E+0_KR, 0.0000000000000000E+0_KR, 0.6505581557685621E-2_KR /)
         xyzw(:,18)=(/ 0.2582068959496968E+0_KR, 0.9660896432961190E+0_KR, 0.0000000000000000E+0_KR, 0.9212586853640415E-2_KR /)
         xyzw(:,19)=(/ 0.4172752955306717E+0_KR, 0.9087801316819105E+0_KR, 0.0000000000000000E+0_KR, 0.1063521204175644E-1_KR /)
         xyzw(:,20)=(/ 0.5700366911792503E+0_KR, 0.8216192370614335E+0_KR, 0.0000000000000000E+0_KR, 0.1134884348397456E-1_KR /)
         xyzw(:,21)=(/ 0.9827986018263947E+0_KR, 0.1771774022615325E+0_KR, 0.5210639477011254E-1_KR, 0.8150269576507465E-2_KR /)
         xyzw(:,22)=(/ 0.9624249230326228E+0_KR, 0.2475716463426288E+0_KR, 0.1115640957156485E+0_KR, 0.9343135395662094E-2_KR /)
         xyzw(:,23)=(/ 0.9402007994128811E+0_KR, 0.3354616289066489E+0_KR, 0.5905888853235372E-1_KR, 0.1005124658581385E-1_KR /)
         xyzw(:,24)=(/ 0.9320822040143202E+0_KR, 0.3173615246611977E+0_KR, 0.1746551677578629E+0_KR, 0.1018093606152102E-1_KR /)
         xyzw(:,25)=(/ 0.9043674199393299E+0_KR, 0.4090268427085357E+0_KR, 0.1217235051095989E+0_KR, 0.1066054174603431E-1_KR /)
         xyzw(:,26)=(/ 0.8912407560074747E+0_KR, 0.3854291150669224E+0_KR, 0.2390278479381724E+0_KR, 0.1075216275547464E-1_KR /)
         xyzw(:,27)=(/ 0.8676435628462708E+0_KR, 0.4932221184851285E+0_KR, 0.6266250624154199E-1_KR, 0.1106243828651345E-1_KR /)
         xyzw(:,28)=(/ 0.8581979986041619E+0_KR, 0.4785320675922435E+0_KR, 0.1857505194547337E+0_KR, 0.1107228969613374E-1_KR /)
         xyzw(:,29)=(/ 0.8396753624049856E+0_KR, 0.4507422593157064E+0_KR, 0.3029466973528983E+0_KR, 0.1112159279420600E-1_KR /)
         xyzw(:,30)=(/ 0.8165288564022188E+0_KR, 0.5632123020762100E+0_KR, 0.1267774800684282E+0_KR, 0.1133655307687399E-1_KR /)
         xyzw(:,31)=(/ 0.8015469370783529E+0_KR, 0.5434303569693900E+0_KR, 0.2494112162362238E+0_KR, 0.1132241512838555E-1_KR /)
         xyzw(:,32)=(/ 0.7773563069070351E+0_KR, 0.5123518486419871E+0_KR, 0.3649832260597654E+0_KR, 0.1133825034038340E-1_KR /)
         xyzw(:,33)=(/ 0.7661621213900394E+0_KR, 0.6394279634749102E+0_KR, 0.6424549224220785E-1_KR, 0.1150830253434940E-1_KR /)
         xyzw(:,34)=(/ 0.7553584143533510E+0_KR, 0.6269805509024392E+0_KR, 0.1906018222779231E+0_KR, 0.1147507934820083E-1_KR /)
         xyzw(:,35)=(/ 0.7344305757559503E+0_KR, 0.6031161693096310E+0_KR, 0.3112275947149608E+0_KR, 0.1144521609262729E-1_KR /)
         xyzw(:,36)=(/ 0.7043837184021765E+0_KR, 0.5693702498468441E+0_KR, 0.4238644781522338E+0_KR, 0.1144263581397218E-1_KR /)

      case default
         write(*,*) 'LebedevAngGrid: Octahedral order not implemented'
         stop
      end select
       
      do i_pnt=1, nPoints
         xyzUnique(:,i_pnt) = xyzw(1:3,i_pnt)
         wUnique(i_pnt)     = xyzw(4,i_pnt) / countMulteplicity(xyzUnique(:,i_pnt),weightsNorm, octaOper)
      end do

      deallocate(xyzw)
   end subroutine


   real(KREAL) function CountMulteplicity(point, weightsNorm, oper)
      ! ----------------------------------------------------------
      ! Count the number of symmetry-equivalent octahedral points 
      ! of a given point. If weightsNorm=true then just return 1.0
      ! ----------------------------------------------------------

      real(KREAL), intent(in)   :: point(:)
      logical    , intent(in)   :: weightsNorm
      real(KREAL)               :: oper (:,:,:)     ! oper(3,3,nOper)

      real(KREAL)     :: octa(3,48)
      integer(KINT)   :: i_oper, double, nOper
      
      ! If the normalized points are needed, just return one:
      if (weightsNorm) then
         countMulteplicity=1.0_KREAL
         return
      end if
   
      nOper = size(oper(1,1,:))

      do i_oper=1, nOper
         octa(:,i_oper) = matmul(oper(:,:,i_oper), point(:))
      end do

      double=0_KINT

      do i_oper=1, nOper
         if (sum((point(:)-octa(:,i_oper))**2)< 1.0e-5_KREAL) double = double+1
      end do

      countMulteplicity=real(double,KREAL)
   end function


   integer(KINT) function FindOverlappingPoints(point, set, nPoints, tolerance)
      ! =====================================================================
      ! If point is in set return the match's index, otherwise return -1_KINT
      ! =====================================================================

      integer(KINT),    intent(in)  :: nPoints
      real(KREAL),      intent(in)  :: point(3)
      real(KREAL),      intent(in)  :: set(3,nPoints)     ! set(3,nPoints)
      real(KREAL),      intent(in)  :: tolerance
      
      integer(KINT)                 :: i_pnt
      real(KREAL)                   :: distance

      findOverlappingPoints = -1_KINT

      if (nPoints==0_KINT) return

      do i_pnt=1, nPoints
         distance = sum((point(:)-set(:,i_pnt))**2)
         if (distance < tolerance) then
            findOverlappingPoints = i_pnt
            return
         end if
      end do
   end function


   subroutine fillOctahedralOperator( octaOper )
      real(KREAL),   intent(out)    :: octaOper (:,:,:)     ! size(3,3,48) 
      
      octaOper(:,:,1)  = reshape((/one ,  zero,  zero,  zero,  one ,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,2)  = reshape((/mone,  zero,  zero,  zero,  mone,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,3)  = reshape((/zero,  one ,  zero,  zero,  zero,  one ,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,4)  = reshape((/zero,  one ,  zero,  mone,  zero,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,5)  = reshape((/one ,  zero,  zero,  zero,  zero,  mone,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,6)  = reshape((/zero,  zero,  one ,  zero,  one ,  zero,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,7)  = reshape((/zero,  zero,  one ,  zero,  mone,  zero,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,8)  = reshape((/mone,  zero,  zero,  zero,  zero,  one ,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,9)  = reshape((/mone,  zero,  zero,  zero,  mone,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,10) = reshape((/zero,  zero,  mone,  mone,  zero,  zero,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,11) = reshape((/zero,  one ,  zero,  zero,  zero,  mone,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,12) = reshape((/zero,  zero,  one ,  mone,  zero,  zero,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,13) = reshape((/zero,  mone,  zero,  zero,  zero,  mone,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,14) = reshape((/zero,  mone,  zero,  zero,  zero,  one ,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,15) = reshape((/zero,  mone,  zero,  one ,  zero,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,16) = reshape((/zero,  zero,  mone,  zero,  mone,  zero,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,17) = reshape((/mone,  zero,  zero,  zero,  zero,  mone,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,18) = reshape((/zero,  zero,  mone,  zero,  one ,  zero,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,19) = reshape((/one ,  zero,  zero,  zero,  zero,  one ,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,20) = reshape((/zero,  zero,  mone,  one ,  zero,  zero,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,21) = reshape((/zero,  one ,  zero,  one ,  zero,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,22) = reshape((/one ,  zero,  zero,  zero,  mone,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,23) = reshape((/mone,  zero,  zero,  zero,  one ,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,24) = reshape((/zero,  mone,  zero,  mone,  zero,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,25) = reshape((/zero,  zero,  one ,  one ,  zero,  zero,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,26) = reshape((/zero,  zero,  mone,  mone,  zero,  zero,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,27) = reshape((/zero,  mone,  zero,  zero,  zero,  mone,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,28) = reshape((/zero,  mone,  zero,  one ,  zero,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,29) = reshape((/mone,  zero,  zero,  zero,  zero,  one ,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,30) = reshape((/zero,  zero,  mone,  zero,  mone,  zero,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,31) = reshape((/zero,  zero,  mone,  zero,  one ,  zero,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,32) = reshape((/one ,  zero,  zero,  zero,  zero,  mone,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,33) = reshape((/one ,  zero,  zero,  zero,  one ,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,34) = reshape((/zero,  zero,  one ,  one ,  zero,  zero,  zero,  mone,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,35) = reshape((/zero,  mone,  zero,  zero,  zero,  one ,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,36) = reshape((/zero,  zero,  mone,  one ,  zero,  zero,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,37) = reshape((/zero,  one ,  zero,  zero,  zero,  one ,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,38) = reshape((/zero,  one ,  zero,  zero,  zero,  mone,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,39) = reshape((/zero,  one ,  zero,  mone,  zero,  zero,  zero,  zero,  mone/), shape(octaOper(:,:,1)))
      octaOper(:,:,40) = reshape((/zero,  zero,  one ,  zero,  one ,  zero,  one ,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,41) = reshape((/one ,  zero,  zero,  zero,  zero,  one ,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,42) = reshape((/zero,  zero,  one ,  zero,  mone,  zero,  mone,  zero,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,43) = reshape((/mone,  zero,  zero,  zero,  zero,  mone,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,44) = reshape((/zero,  zero,  one ,  mone,  zero,  zero,  zero,  one ,  zero/), shape(octaOper(:,:,1)))
      octaOper(:,:,45) = reshape((/zero,  mone,  zero,  mone,  zero,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,46) = reshape((/mone,  zero,  zero,  zero,  one ,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,47) = reshape((/one ,  zero,  zero,  zero,  mone,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
      octaOper(:,:,48) = reshape((/zero,  one ,  zero,  one ,  zero,  zero,  zero,  zero,  one /), shape(octaOper(:,:,1)))
   end subroutine 


   integer(KINT) function MinAcceptableLebedevOrder (l)
      integer(KINT), intent(in) :: l
      
      integer(KINT) :: i

      MinAcceptableLebedevOrder = -1_KINT
      do i=1, size(orders)
         if ( orders(i) >= l ) then
            MinAcceptableLebedevOrder = orders(i)
            return
         end if
      end do

   end function


   subroutine createSpecialGrid (refAxis, lOrder, grid, weights, nPoints)
      real(KREAL),   intent(in)  :: refAxis(:)    !size(3)
      integer(KINT), intent(in)  :: lOrder 
      real(KREAL),   intent(out) :: grid(:,:)     !size(3,maxNPoints)
      real(KREAL),   intent(out) :: weights(:)    !size(maxNPoints)
      integer(KINT), intent(out) :: nPoints

      integer(KINT), parameter :: maxNPoints = 10000_KINT 
      real(KREAL) :: thetaGrid(maxNPoints), thetaWeights(maxNPoints)
      real(KREAL) :: phiGrid(maxNPoints), phiWeights(maxNPoints)
      real(KREAL) :: axisOrig(3), rotMatrix(3,3)
      integer(KINT) :: numThetaPoints, numPhiPoints, ipnt, jpnt, iGlobal, totalNPoints, iRot
      real(KREAL) :: thetaRot, phiRot

      axisOrig = (/0.0_KREAL, 0.0_KREAL, 1.0_KREAL/)

      select case(LOrder)    
      case(5_KINT)
         numPhiPoints   = 3_KINT    ! 1
         numThetaPoints = 8_KINT    ! 3
      case(7_KINT)
         numPhiPoints   = 3_KINT    ! 2
         numThetaPoints = 8_KINT    ! 4
      case(11_KINT)
         numPhiPoints   = 3_KINT    ! 2
         numThetaPoints = 8_KINT    ! 5
      case(13_KINT)
         numPhiPoints   = 3_KINT    ! 3
         numThetaPoints = 8_KINT    ! 7
      case(17_KINT)
         numPhiPoints   = 3_KINT
         numThetaPoints = 8_KINT
      case(19_KINT)
         numPhiPoints   = 3_KINT
         numThetaPoints = 9_KINT
      case(21_KINT)
         numPhiPoints   = 4_KINT
         numThetaPoints = 10_KINT
      case(23_KINT)
         numPhiPoints   = 4_KINT
         numThetaPoints = 10_KINT
      case(29_KINT)
         numPhiPoints   = 5_KINT
         numThetaPoints = 13_KINT
      case(31_KINT)
         numPhiPoints   = 5_KINT
         numThetaPoints = 14_KINT
      case(35_KINT)
         numPhiPoints   = 6_KINT
         numThetaPoints = 15_KINT
      case(41_KINT)
         numPhiPoints   = 6_KINT
         numThetaPoints = 17_KINT
      case(47_KINT)
         numPhiPoints   = 7_KINT
         numThetaPoints = 20_KINT
      case(53_KINT)
         numPhiPoints   = 8_KINT
         numThetaPoints = 23_KINT
      case(59_KINT)
         numPhiPoints   = 10_KINT
         numThetaPoints = 26_KINT
      case default
         write(*,*) 'LebedevAngGrid: Octahedral order not implemented'
         stop
      end select

      call legpnt (numThetaPoints, thetaGrid, thetaWeights)
      thetaWeights = thetaWeights * pi / 2.0_KREAL
      thetaGrid = (thetaGrid + 1.0_KREAL) * pi / 2.0_KREAL

!       do ipnt=1, numThetaPoints
!          thetaGrid(ipnt)     = (pi) * (ipnt-0.5_KREAL) / (numThetaPoints)
!          thetaWeights(ipnt)  = (pi) / numThetaPoints 
!       end do

      !call legpnt (numPhiPoints, phiGrid, phiWeights)
      !phiWeights = phiWeights * pi / 6.0_KREAL
      !phiGrid = (phiGrid + 1.0_KREAL) * pi / 3.0_KREAL

      do ipnt=1, numPhiPoints
         phiGrid(ipnt)     = (2*pi/6) * ipnt / (numPhiPoints)
         phiWeights(ipnt)  = (2*pi/6) / numPhiPoints ! I need a three because I generate only the symmetry unique wedge
      end do

      thetaRot = acos( dot_product(axisOrig,refAxis) / &
                       (sqrt(dot_product(axisOrig,axisOrig)) * sqrt(dot_product(refAxis,refAxis))) )
      phiRot = 0.0_KREAL
      if (abs(refAxis(2)) > 1.0E-10_KREAL) phiRot = atan(refAxis(1)/refAxis(2))

      call calcRotationMat(rotMatrix, thetaRot, phiRot)

      iGlobal = 0_KINT
      do ipnt=1, numThetaPoints
         do jpnt=1, numPhiPoints
            do iRot=0,5
            iGlobal = iGlobal+1
               grid(1,iGlobal)  = sin(thetaGrid(ipnt)) * cos(phiGrid(jpnt)+iRot*2.0_KREAL*pi/6.0_KREAL)
               grid(2,iGlobal)  = sin(thetaGrid(ipnt)) * sin(phiGrid(jpnt)+iRot*2.0_KREAL*pi/6.0_KREAL)
               grid(3,iGlobal)  = cos(thetaGrid(ipnt))
               weights(iGlobal) = sin(thetaGrid(ipnt)) * thetaWeights(ipnt) * phiWeights(jpnt)
               grid(:,iGlobal)  = matmul(rotMatrix,grid(:,iGlobal))
            end do
         end do
      end do

      nPoints = iGlobal

   end subroutine 



   subroutine calcRotationMat(rotMatrix, theta, phi)
      real(KREAL), intent(out) :: rotMatrix(:,:) 
      real(KREAL), intent(in) :: theta
      real(KREAL), intent(in) :: phi

      real(KREAL) :: firstRot(3,3), secondRot(3,3)

      firstRot    = 0.0_KREAL
      secondRot   = 0.0_KREAL

      firstRot(1,1) = 1.0_KREAL
      firstRot(2,2) = cos(theta)
      firstRot(2,3) = sin(theta)
      firstRot(3,2) = -sin(theta)
      firstRot(3,3) = cos(theta)

      secondRot(1,1) = cos(phi)
      secondRot(1,2) = sin(phi)
      secondRot(2,1) = -sin(phi)
      secondRot(2,2) = cos(phi) 
      secondRot(3,3) = 1.0_KREAL

      rotMatrix = matmul(secondRot, firstRot)

   end subroutine 


subroutine pntred (ldim, avec, bvec, noper, oper, transl, np, point, weight)

   use Vartypes
   implicit real(KREAL) (a-h, o-z)
   implicit integer(KINT) (i-n)



   real(KREAL), intent(in)      :: avec(3,3), bvec(3,3), oper(3,3,noper), transl(3,noper)
   integer(KINT), intent(inout) :: np
   real(KREAL), intent(inout)   :: point(3,np), weight(np)

!  ======================================================================
!  purpose: reduction of a set of points and weights to the generators:
!           the symmetry unique subset. The weights of the reduced set
!           do NOT represent the equivalent points. The sum of weights
!           of the reduced set should, therefore, be 1/noper times the
!           total (unreduced) sum-of-weights, assuming that the input set
!           is symmetric under the operators
!           This is explicitly tested (at the expense of more CPU time)
!
!  in-out : np      nr. of points
!           point   coordinates of the points
!           weight  weights of the points
!
!  remark * the set of operators is assumed to constitute a group
!         * it is assumed that the first operator is the identity
!  ======================================================================

   real(KREAL), parameter :: tol = 1e-6_KREAL, big = 1e20_KREAL

   real(KREAL), parameter :: zero=0d0, one=1d0, two=2d0, three=3d0, four=4d0, half=0.5d0,          &
          third=one/three


   integer(KINT) :: idequi(np)
   real(KREAL)   :: proj(3), aux(3), dmin, projmin(3)

   if (noper<=1) return

   sumw1  = sum(weight)
   idequi = 0

   ip_: do ip = 1, np
      if (idequi(ip)==-1) cycle ip_

!     ------------------------------------------
!     ip mapped onto itself by the unit operator
!     NB: assumed to be the FIRST operator
!     ------------------------------------------

      idequi(ip) = 1

!     --------------------
!     scan other operators
!     --------------------

      iop_: do iop = 2, noper

         proj = matmul(oper(:,:,iop),point(:,ip)) + transl(:,iop)

!        ------------------------------------------------------------
!        check equivalency with other points. It would be feasible
!        to start the loop at ip, because earlier points have already
!        been scanned. However, to verify that the set of points
!        actually is symmetric, we have to go over all the points and
!        because earlier points cannot be related to this point
!        ------------------------------------------------------------

         dmin    = big
         projmin = zero

         jp_: do jp = 1, np

!------------------------------------------
!***      jp_: do jp = ip, np
!***        if (idequi(jp) == -1) cycle jp_
!------------------------------------------

            aux = point(:,jp) - proj

            if (ldim==0) then
               d = abs(aux(1)) + abs(aux(2)) + abs(aux(3))

            else if (ldim==1) then
               s1 = aux(1)*bvec(1,1)
               s1 = s1 - nint(s1)
               d  = abs(s1) + abs(aux(2)) + abs(aux(3))

            else if (ldim==2) then
               s1 = aux(1)*bvec(1,1) + aux(2)*bvec(2,1)
               s2 = aux(1)*bvec(1,2) + aux(2)*bvec(2,2)
               s1 = s1 - nint(s1)
               s2 = s2 - nint(s2)
               d  = abs(s1) + abs(s2) + abs(aux(3))

            else if (ldim==3) then
               s1 = aux(1)*bvec(1,1) + aux(2)*bvec(2,1) + aux(3)*bvec(3,1)
               s2 = aux(1)*bvec(1,2) + aux(2)*bvec(2,2) + aux(3)*bvec(3,2)
               s3 = aux(1)*bvec(1,3) + aux(2)*bvec(2,3) + aux(3)*bvec(3,3)
               s1 = s1 - nint(s1)
               s2 = s2 - nint(s2)
               s3 = s3 - nint(s3)
               d  = abs(s1) + abs(s2) + abs(s3)
            end if

            if (d<symtol) then
               if (abs(weight(ip)-weight(jp))>tol) then
                  write(*,*) 'different weights. PNTRED'
                  stop
               end if

               if (jp==ip) then
                  idequi(ip) = idequi(ip) + 1
               else
                  idequi(jp) = -1
               end if
               cycle iop_
            end if

            if (d<dmin) then
               dmin    = d
               projmin = proj
            end if

         end do jp_

         write (iuout,9000) ip, point(:,ip), oper(:,:,iop), dmin, projmin
         write (iuout,9010) (i,point(:,i),i = 1,np)
         write (iuout,9020) noper, (j,oper(:,:,j),j = 1,noper)
         write(*,*) 'equivalent point not found. PNTRED'

      end do iop_
   end do ip_

!  -----------------
!  reduced point set
!  -----------------

   mp = 0
   do ip = 1, np
      if (idequi(ip)==0) then
         write(*,*) " fuckkkkkk"
         stop 
      end if
      if (idequi(ip)>0) then
         mp          = mp + 1
         point(:,mp) = point(:,ip)
         weight(mp)  = weight(ip)/idequi(ip)
      end if
   end do
   np = mp

   sumw2 = sum(weight(1:np))*noper
   if (abs(sumw1-sumw2)>tol) then
      write(*,*) 'sum-of-weights changed. PNTRED'
      stop
   end if

 9000 format(' mapping failed in PNTRED for point ',i4,' at',3e12.4/' under operator',3(t20,       &
             3e12.4/),' dmin=',e12.4/' p   =',3e12.4)
 9010 format(' input list'/(1x,i6,4x,3e12.4))
 9020 format(' all operators:',i6/(i4,3(t10,3e12.4/)))
end subroutine pntred



end module
