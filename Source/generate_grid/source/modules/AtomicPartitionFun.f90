module AtomicPartitionFun

   use vartypes

   implicit none
   private
   
   public CalcPartition
   public CalcPartitionYukawaLike, CalcPartitionBecke, CalcPartitionFermi

   ! TODO comment this module
   ! TODO implement the voronoi partition (the real voronoi cells)

contains 

  
   subroutine CalcPartition( points, xyzAtom, qAtoms, indexAtom, partitionFunction, functionType, sizeAdjusted )
      
      real(KREAL),      intent(in)  ::  points(:,:)                 ! points(3,numPoints)
      real(KREAL),      intent(in)  ::  xyzAtom(:,:)                ! xyzAtom(3,nAtoms) 
      real(KREAL),      intent(in)  ::  qAtoms(:)                   ! qAtoms(nAtoms) 
      integer(KINT),    intent(in)  ::  indexAtom
      real(KREAL),      intent(out) ::  partitionFunction(:)        ! partitionFunction(numPoints)
      character*(*),    intent(in)  :: functionType
      logical,          intent(in)  :: sizeAdjusted
   
      character(LCHARS) :: functionTypeLoc
      real(KREAL)       :: qAtomsLoc(size(qAtoms))

      functionTypeLoc = functionType

      if (sizeAdjusted) then
         qAtomsLoc = qAtoms
      else 
         qAtomsLoc = 1.0_KREAL
      end if

      select case (functionTypeLoc)
         case ('FUZZYVORONOIBECKE')
            call CalcPartitionBecke( points, xyzAtom, qAtomsLoc, indexAtom, partitionFunction  )
         case ('FUZZYVORONOIFERMI')
            call CalcPartitionFermi( points, xyzAtom, qAtomsLoc, indexAtom, partitionFunction  )
         case ('YUKAWALIKE')
            call CalcPartitionYukawaLike( points, xyzAtom, qAtomsLoc, indexAtom, partitionFunction  )
         case default 
            partitionFunction = 0.0_KREAL
            write(*,*) "Error: Partition function not recognize: ", functionTypeLoc
            stop
      end select 

   end subroutine 



   subroutine CalcPartitionBecke( points, xyzAtom, qAtoms, indexAtom, partitionFunction  )
      real(KREAL),      intent(in)  ::  points(:,:)                 ! points(3,numPoints)
      real(KREAL),      intent(in)  ::  xyzAtom(:,:)                ! xyzAtom(3,nAtoms) 
      real(KREAL),      intent(in)  ::  qAtoms(:)                   ! qAtoms(nAtoms) 
      integer(KINT),    intent(in)  ::  indexAtom
      real(KREAL),      intent(out) ::  partitionFunction(:)        ! partitionFunction(numPoints)
     
      integer(KINT)             ::  nAtoms, nPoints
      real(KREAL),  allocatable ::  psum(:), pjatm(:), djatm(:), hypr(:) ! wbecke workspace 
      real(KREAL),  allocatable ::  distAtoms(:,:), cmid(:,:)

      nPoints = size(points(1,:))
      nAtoms  = size(xyzAtom(1,:))

      allocate ( distAtoms(nAtoms,nAtoms), cmid(nAtoms,nAtoms), psum(nPoints), pjatm(nPoints), djatm(nPoints),  &
                 hypr(nPoints) )

      call calcDinstanceAtomsLoc(xyzAtom, qAtoms, distAtoms, cmid)

      partitionFunction = 1.0_KREAL

      call wbecke (nPoints, points, partitionFunction, nAtoms, indexAtom, xyzAtom, &
                   distAtoms, cmid, psum, pjatm, djatm, hypr)
      

      deallocate ( distAtoms, cmid, psum, pjatm, djatm, hypr )
   end subroutine 



   subroutine CalcPartitionYukawaLike( points, xyzAtom, qAtoms, indexAtom, partitionFunction  )
      real(KREAL),      intent(in)  ::  points(:,:)                 ! points(3,numPoints)
      real(KREAL),      intent(in)  ::  xyzAtom(:,:)                ! xyzAtom(3,nAtoms) 
      real(KREAL),      intent(in)  ::  qAtoms(:)                   ! qAtoms(nAtoms) 
      integer(KINT),    intent(in)  ::  indexAtom
      real(KREAL),      intent(out) ::  partitionFunction(:)        ! partitionFunction(numPoints)

      real(KREAL), allocatable  :: pfun(:,:), r(:,:)              ! size(nPoints, nAtoms)
      real(KREAL), allocatable  :: norm(:)                        ! size(nPoints)
      integer(KINT)             :: nPoints, nAtoms, i_atom, i_pnt, supp(1), i_q
      real(KREAL) :: radius

      nPoints = size(points(1,:))
      nAtoms  = size(xyzAtom(1,:))

      allocate (pfun(nPoints, nAtoms), norm(nPoints), r(nPoints, nAtoms))

      norm = 0.0_KREAL

      do i_pnt=1, nPoints
         do i_atom=1, nAtoms
            r(i_pnt, i_atom) = sqrt(sum((points(:,i_pnt) - xyzAtom(:,i_atom))**2))
         end do
      end do

      ! Yukawa like partition
      do i_pnt=1, nPoints
         do i_atom=1, nAtoms
            i_q = nint(qAtoms(i_atom))
            radius = partitionFunEffectiveRadius(i_q)
            pfun(i_pnt, i_atom) = radius*exp(-2.0_KREAL*(r(i_pnt, i_atom))) / r(i_pnt, i_atom)**3 
            norm(i_pnt) = norm(i_pnt) + pfun(i_pnt, i_atom)
         end do
      end do

      do i_pnt=1, nPoints
         if (norm(i_pnt) > 1.0E-200_KREAL ) then
            partitionFunction(i_pnt) = pfun(i_pnt, indexAtom) / norm(i_pnt)
         else
            write(*,*) "WARNING: small weight in partition (AtomicPartitionFun)"
            partitionFunction(i_pnt) = 0.0_KREAL
            supp = minloc(r(i_pnt, :))
            if ( supp(1) == indexAtom  ) partitionFunction(i_pnt) = 1.0_KREAL
         end if
      end do

      deallocate (pfun, norm, r)
   end subroutine 


   subroutine CalcPartitionFermi( points, xyzAtom, qAtoms, indexAtom, partitionFunction )
      real(KREAL),      intent(in)  ::  points(:,:)                 ! points(3,numPoints)
      real(KREAL),      intent(in)  ::  xyzAtom(:,:)                ! xyzAtom(3,nAtoms) 
      real(KREAL),      intent(in)  ::  qAtoms(:)                   ! qAtoms(nAtoms) 
      integer(KINT),    intent(in)  ::  indexAtom
      real(KREAL),      intent(out) ::  partitionFunction(:)        ! partitionFunction(numPoints)

      integer(KINT)             :: nPoints, nAtoms, i_atom, j_atom, i_pnt
      real(KREAL), allocatable  :: pfun(:,:), distAtomAtom(:,:), distAtomPoint(:,:)
      real(KREAL), allocatable  :: norm(:), mu(:)                                   ! size(nPoints)
      real(KREAL),  parameter   :: alfa = 3.0_KREAL ! it's just a reasonable number, is not super-optimized

      nPoints = size(points(1,:))
      nAtoms  = size(xyzAtom(1,:))

      allocate ( pfun(nPoints, nAtoms), norm(nPoints), distAtomAtom(nAtoms, nAtoms), distAtomPoint(nAtoms,nPoints), & 
                 mu(nPoints) )
   
      norm = 0.0_KREAL
      pfun = 1.0_KREAL
      
      distAtomAtom = 0.0_KREAL
      do i_atom=1, nAtoms
         do j_atom=i_atom+1, nAtoms
            distAtomAtom(i_atom, j_atom) = sqrt(sum((xyzAtom(:,i_atom) - xyzAtom(:,j_atom))**2))
            distAtomAtom(j_atom, i_atom) = distAtomAtom(i_atom, j_atom)
         end do
      end do
      
      do i_atom=1, nAtoms
         do i_pnt=1, nPoints
            distAtomPoint(i_atom, i_pnt) = sqrt(sum((points(:,i_pnt)-xyzAtom(:,i_atom))**2))
         end do
      end do


      do i_atom=1, nAtoms
         do j_atom=1, nAtoms
            if (i_atom /= j_atom) then

            mu(:) = (distAtomPoint(i_atom,:) - distAtomPoint(j_atom,:))/ distAtomAtom(i_atom, j_atom)

            pfun(:, i_atom) = pfun(:, i_atom) / ( ((1.0_KREAL+mu(:))/(1.0_KREAL-mu(:)))**alfa + 1.0_KREAL )
            end if
         end do
         pfun(:, i_atom) = pfun(:, i_atom) * qAtoms(i_atom)
         norm(:) = norm(:) + pfun(:, i_atom)
      end do

      do i_pnt=1, nPoints
         if (norm(i_pnt) > 1.0E-200_KREAL ) then
            partitionFunction(i_pnt) = pfun(i_pnt, indexAtom) / norm(i_pnt)
         else
            partitionFunction(i_pnt) = 0.0_KREAL
            write(*,*) "WARNING: normalization in AtomicPartitionFun too small"
         end if
      end do

      deallocate( pfun, norm, distAtomAtom, distAtomPoint, mu )
   end subroutine
 



   subroutine CalcDinstanceAtomsLoc(xyzAtoms, qAtoms, distAtoms, cmid)
      ! ---------------------------------------------------------------------------
      ! Calculate the distance between atom-pairs and cmid (for becke partitioning)
      ! This piece of code comes from rpntb.
      ! ---------------------------------------------------------------------------

      real(KREAL),     intent(in)   ::  xyzAtoms(:,:)       ! xyzAtoms(3,nAtoms) 
      real(KREAL),     intent(in)   ::  qAtoms(:)           ! qAtoms(nAtoms)

      real(KREAL), intent(out)      :: distAtoms(:,:)       ! distAtoms(nAtoms, nAtoms)
      real(KREAL), intent(out)      :: cmid(:,:)            ! cmid(nAtoms, nAtoms)

      integer(KINT)                 :: nAtoms
      integer(KINT)                 :: iatm, jatm, nuc1, nuc2
      real(KREAL)                   :: rs1, rs2, xmid, ratio

      nAtoms = size(xyzAtoms(1,:))

      do iatm = 1, nAtoms
         do jatm = 1, nAtoms
            distAtoms(jatm,iatm) = sqrt((xyzAtoms(1,jatm)-xyzAtoms(1,iatm))**2+                                &
            (xyzAtoms(2,jatm)-xyzAtoms(2,iatm))**2+(xyzAtoms(3,jatm)-xyzAtoms(3,iatm))**2)

            nuc1 = nint(qAtoms(iatm))
            nuc2 = nint(qAtoms(jatm))

            rs1=1 
            rs2=1
!             call elinfo (nuc1, 'Sphere', rs1)
!             call elinfo (nuc2, 'Sphere', rs2)

            ratio = rs1/rs2

            xmid            = (ratio-1.0_KREAL)/(ratio+1.0_KREAL)
            cmid(iatm,jatm) = xmid/(xmid**2-1.0_KREAL)

            if (cmid(iatm,jatm)>0.5_KREAL) cmid(iatm,jatm) = 0.5_KREAL
            if (cmid(iatm,jatm)<(-0.5_KREAL)) cmid(iatm,jatm) = - 0.5_KREAL
         end do
      end do

   end subroutine
 


   subroutine wbecke (npnts, pnts, wghts, natm, iatm, xyzatm, datm, cmid, psum, pjatm, djatm, hypr)
   !  ======================================================================
   !  purpose:  correct weight of a number of integration points centered
   !            around atom iatm.
   !
   !  input  :  npnts  : the number of integration points.
   !            pnts   : the integration pooins.
   !            natm   : the total number of atoms.
   !            iatm   : atom.
   !            xyzatm : the coordinates of the atoms.
   !            datm   : datm(i,j) gives the distance between atom i and j.
   !            cmid   : factor for atomic size adjustment.
   !  in-out :  wghts  : the integration weights.
   !
   !  ======================================================================

      real(KREAL),  parameter ::    zero    = 0.0_KREAL
      real(KREAL),  parameter ::    half    = 0.5_KREAL
      real(KREAL),  parameter ::    one     = 1.0_KREAL
      real(KREAL),  parameter ::    onehlf  = 1.5_KREAL
   
      integer(KINT),    intent(in)      :: npnts
      real(KREAL),      intent(in)      :: pnts(3,npnts)
      real(KREAL),      intent(inout)   :: wghts(npnts)
      integer(KINT),    intent(in)      :: natm
      integer(KINT),    intent(in)      :: iatm
      real(KREAL),      intent(in)      :: xyzatm(3,natm)
      real(KREAL),      intent(in)      :: datm(natm,natm)
      real(KREAL),      intent(in)      :: cmid(natm,natm)
      real(KREAL),      intent(inout)   :: psum(npnts)
      real(KREAL),      intent(inout)   :: pjatm(npnts)
      real(KREAL),      intent(inout)   :: djatm(npnts)
      real(KREAL),      intent(inout)   :: hypr(npnts)

      real(KREAL)       :: distAtomsPoints(npnts,natm), rm1
      integer(KINT)     :: jatm, katm

   !  ==============================================
   !  psum is the sum of the nuclear weights.
   !  it is an accumulated sum, so we start at zero.
   !  ==============================================
      
      psum(:) = zero

      do jatm=1, natm
         distAtomsPoints(:,jatm) = sqrt((pnts(1,:)-xyzatm(1,jatm))**2+(pnts(2,:)-xyzatm(2,jatm))**2+ &
                                      (pnts(3,:)-xyzatm(3,jatm))**2)
      end do

      do jatm = 1, natm
   !     -------------------------------------------------------
   !     *** pjatm is the nuclear weight function for atom jatm.
   !     *** it is an accumulated product, so we start at one.
   !     -------------------------------------------------------
         pjatm(:) = one
         do katm = 1, natm
            if (katm/=jatm) then
               rm1 = 1.0_KREAL/datm(jatm,katm)
   !           ------------------------------
   !           *** the hyperbolic coordinate.
   !           ------------------------------
               hypr(:) = (distAtomsPoints(:,jatm)-distAtomsPoints(:,katm))*rm1
               hypr(:) = hypr(:) + cmid(jatm,katm)*(one-hypr(:)**2)
   !           ---------------------
   !           *** three iterations.
   !           ---------------------
               hypr(:) = (onehlf-half*hypr(:)**2)*hypr(:)
               hypr(:) = (onehlf-half*hypr(:)**2)*hypr(:)
               hypr(:) = (onehlf-half*hypr(:)**2)*hypr(:)
               
               pjatm(:) = pjatm(:)*(one-hypr(:))*half
            end if
         end do 
            
         psum(:) = psum(:) + pjatm(:)
         if (jatm==iatm) then
   !        -------------------------------
   !        *** include the nuclear weight.
   !        -------------------------------
            wghts(:) = wghts(:)*pjatm(:)
         end if
      end do 
   !  --------------------------
   !  *** normalize the weights.
   !  --------------------------
      wghts(:) = wghts(:)/psum(:)
      return
   end subroutine 


   real(KREAL) function partitionFunEffectiveRadius(charge)
      integer(KINT),    intent(in) :: charge

      ! here I define the radii of the atoms used by the partition functions. 
      ! Now I'll use a very basic scheme: 0.x for hydrogen and 1 for all other atoms.

      if (charge == 1_KINT) then
         partitionFunEffectiveRadius = 0.3_KREAL
         return
      else
         partitionFunEffectiveRadius = 1.0_KREAL
         return
      end if
   end function


end module

