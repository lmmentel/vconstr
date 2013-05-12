module IntegralsModule
  use varModule
  use dictionaryModule
  use ioModule
  implicit none 

  private

  type energyType
      real(DP) :: nuclearRepulsion
      real(DP) :: electronic
      real(DP) :: total
      real(DP) :: sz
      real(DP) :: szz
      real(DP) :: core
      real(DP) :: scf
      real(DP) :: electronNucleus
      real(DP) :: electronElectron
      real(DP) :: potential
      real(DP) :: kinetic
  end type

  interface print_energies
      module procedure print_energy_type
  end interface

  public :: energyType, print_energies

! ======================
! one electron integrals 
! ======================
  interface readS
      module procedure readOverlap
  end interface
  interface readH
      module procedure readCoreHamiltonian
  end interface
  interface readT
      module procedure readKineticEnergyIntegrals
  end interface
  interface readXdip
      module procedure readXdipoleIntegrals
  end interface
  interface readYdip
      module procedure readYdipoleIntegrals
  end interface
  interface readZdip
      module procedure readZdipoleIntegrals
  end interface

  public :: readS, readH, readT, readXdip, readYdip, readZdip

! ======================
! two electron integrals 
! ======================
  interface ReadTwoIntAO
    module procedure ReadTwoIntAO
  end interface

  interface ReadTwoIntMO
    module procedure ReadTwoIntMO
  end interface

  interface addr
    module procedure address
  end interface 
 
  interface print4i
    module procedure print4indexVec
  end interface 

  public ReadTwoIntAO, ReadTwoIntMO, addr, print4i

! ===============================
! orbitals and occupation numbers 
! ===============================
  interface readMOs
      module procedure readHFMolecularOrbitals
  end interface
  interface readNOs
      module procedure readNaturalOrbitals
  end interface
  interface readOcc
      module procedure readOccupationNumbers
  end interface
  interface readEnergies
      module procedure readEnergiesCommonBlock
  end interface

  public :: readMOs, readNOs, readOcc, readEnergies

! ===============================
! orbitals and denities together 
! ===============================
  interface readOneEintegrals
      module procedure readOneElectronIntegralsInAO
  end interface
  interface getOrbitalsAndDensities
      module procedure reshapeOrbitals
  end interface

  public :: readOneEintegrals, getOrbitalsAndDensities

  contains

! subroutine readCoreHamiltonianMatrixInMO()
! end subroutine readCoreHamiltonianMatrixInMO

 subroutine reshapeOrbitals(Vpaomo, Vaomo, Pmomo, Pnomo, nb, rnel, nele, dictionaryFile)
  character(*), intent(in)  :: dictionaryFile
  integer,      intent(in)  :: nele, nb 
  real(DP),     intent(out) :: rnel 
  real(DP),     intent(out) :: Vpaomo(:), Vaomo(:), Pmomo(:), Pnomo(:)

  real(DP), allocatable :: Mpaomo(:,:), Maomo(:,:)
  integer :: i, j, ij

  allocate(Mpaomo(nb, nb), Maomo(nb, nb))
  Mpaomo = reshape(Vpaomo, (/nb, nb/))
  Maomo  = reshape(Vaomo,  (/nb, nb/))

  call readOrbitalsAndComposeDensities(Mpaomo, Maomo, Pmomo, Pnomo, nb, rnel, nele, dictionaryFile)
  ij = 0
  do j = 1, nb
    do i = 1, nb
      ij = ij + 1 
      Vpaomo(ij) = Mpaomo(i,j)
      Vaomo(ij)  = Maomo(i,j)
    enddo
  enddo
  deallocate(Mpaomo, Maomo)
 end subroutine reshapeOrbitals

 subroutine readOrbitalsAndComposeDensities(Vpaomo, Vaomo, Pmomo, Pnomo, nb, rnel, nele, dictionaryFile)
  use dictionaryModule
  use ioModule
  character(*), intent(in)  :: dictionaryFile
  integer,      intent(in)  :: nele, nb 
  real(DP),     intent(out) :: rnel 
  real(DP),     intent(out) :: Vpaomo(:,:), Vaomo(:,:), Pmomo(:), Pnomo(:)

  type(dictionaryType)  :: dictionary
  real(DP), allocatable :: VaomoInv(:,:), Vaono(:,:), Occ(:), Vmono(:,:), Pnono(:)
  integer  :: i, ier
  real(DP) :: rhfnel, det

!  Vpaomo = reshape(Vpaomo, (/nb, nb/))
!  Vaomo  = reshape(Vaomo,  (/nb, nb/))
  call new(dictionary, dictionaryfile)
  allocate(VaomoInv(nb, nb), Vaono(nb, nb), Occ(nb), Vmono(nb, nb), Pnono(nb*(nb+1)/2))
  
! read HF MO's over AO's

  call readMOs(Vaomo, dictionary)

  Pmomo = 0.0_DP
  forall (i = 1:nele/2) 
    Pmomo(i*(i+1)/2) = 2.0_DP
  end forall
  rhfnel = sum(Pmomo)

! HF MOs are in primitive AOs so just copy the Vaomo to Vmopao
! this will need to be changes if the required Vpamo matrix is in 
! fully uncontracted gaussian primitives
 
  Vpaomo = Vaomo

! get NO's over AO's and Occupation numbers

  call readNOs(Vaono, dictionary)
  call readOcc(Occ, dictionary)

  forall (i=1:size(Occ))
    Pnono(i*(i+1)/2) = Occ(i) 
  end forall
  rnel = sum(Occ)

  write(*,'(/"Sum of MO occupation numbers = ",f6.3)') rhfnel
  write(*,'( "Sum of NO occupation numbers = ",f6.3/)') rnel

! transform NO's from AO to MO

  VaomoInv = Vaomo

  call inverseMatrix(VaomoInv)
  Vmono = matmul(VaomoInv, Vaono)
!  call minvr(VaomoInv, 1.0d-7, det, ier, nb)
!  call matml(Vaono, VaomoInv, Vmono, nb)
  call matPrint(Vmono, "NO's in MO basis")

! transform the natural orbital occupation numbers from NO to MO basis
! store the resulting Pnomo in pnono
!  call matTrans(Vmono, Pnono, Vmono)

  call tmtdag(Pnono, nb, Pnomo, nb, Vmono, 0.5_DP)
  call matPrint(Pnono, 'CI-density in MO basis')
  call matPrint(Pmomo, 'HF MO-density in MO basis')
  deallocate(VaomoInv, Vaono, Occ, Vmono,  Pnono)
  call delete(dictionary)
 end subroutine readOrbitalsAndComposeDensities


 subroutine readOneElectronIntegralsInAO(T, Ven, nuclearRepulsion, dictionaryFile)
  use dictionaryModule
  character(*), intent(in)  :: dictionaryFile
  real(DP),     intent(out) :: T(:), Ven(:) 
  real(DP),     intent(out) :: nuclearRepulsion

  type(dictionaryType)  :: dictionary
  type(energyType)      :: energy
  real(DP), allocatable :: Htemp(:)

  call new(dictionary, dictionaryfile)

  allocate(Htemp(size(T)))
! read core hamiltonian matrix elements
  call readH(Htemp, dictionary)
! read kinetic energy operator matrix elements
  call readT(T, dictionary)
! calculate the electron-nuclear attraction matrix elements
  Ven = Htemp - T 
! read energies from the dictionary file
  call readEnergies(energy, dictionary)
  nuclearRepulsion = energy%nuclearRepulsion
  deallocate(Htemp)
  call delete(dictionary) 
 end subroutine readOneElectronIntegralsInAO

 subroutine readOverlap(S, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: S(:)

  call readFrom(dictionary, 'overlap', S)
 end subroutine readOverlap

 subroutine readCoreHamiltonian(H, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: H(:)

  call readFrom(dictionary, 'core hamiltonian', H)
 end subroutine readCoreHamiltonian

 subroutine readKineticEnergyIntegrals(T, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: T(:)

  call readFrom(dictionary, 'kinetic energy', T)
 end subroutine readKineticEnergyIntegrals

 subroutine readXdipoleIntegrals(X, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: X(:)

  call readFrom(dictionary, 'X dipole', X)
 end subroutine readXdipoleIntegrals

 subroutine readYdipoleIntegrals(Y, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Y(:)

  call readFrom(dictionary, 'Y dipole', Y)
 end subroutine readYdipoleIntegrals

 subroutine readZdipoleIntegrals(Z, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Z(:)

  call readFrom(dictionary, 'Z dipole', Z)
 end subroutine readZdipoleIntegrals

 subroutine readHFMolecularOrbitals(Vaomo, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Vaomo(:,:)

  call readFrom(dictionary, 'HF MOs', Vaomo)
 end subroutine readHFMolecularOrbitals

 subroutine readNaturalOrbitals(Vaono, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Vaono(:,:)

  call readFrom(dictionary, 'natural orbitals', Vaono)
 end subroutine readNaturalOrbitals

 subroutine readOccupationNumbers(Occ, dictionary)
  use dictionaryModule
  type(dictionaryType)  :: dictionary
  real(DP), intent(out) :: Occ(:)

  call readFrom(dictionary, 'occupation numbers', Occ)
 end subroutine readOccupationNumbers

 subroutine readEnergiesCommonBlock(self, dictionary)
  use dictionaryModule
  type(energyType)      :: self
  type(dictionaryType)  :: dictionary
  real(DP), allocatable :: Etemp(:)

! extract the values form the common block and store them in the 
! custom type energyType
!COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
!                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP

  allocate(Etemp(115))
  call readFrom(dictionary, 'energies', Etemp)
  self%nuclearRepulsion = Etemp( 1)
  self%electronic       = Etemp( 2)
  self%total            = Etemp( 3)
  self%sz               = Etemp( 4)
  self%szz              = Etemp( 5)
  self%core             = Etemp( 6)
  self%scf              = Etemp( 7)
  self%electronNucleus  = Etemp(11)
  self%electronElectron = Etemp(12)
  self%potential        = Etemp(13)
  self%kinetic          = Etemp(14)
  deallocate(Etemp)
 end subroutine readEnergiesCommonBlock

 subroutine print_energy_type(self)
  type(energyType) :: self

  write(*,'(/4x,a)') repeat('=',len(trim('Energies from dictionary file'))) 
  write(*,'(4x,a)') trim('Energies from dictionary file') 
  write(*,'(4x,a/)') repeat('=',len(trim('Energies from dictionary file'))) 

  write(*,'(1x,a40," = ",f12.8)') 'Nuclear repulsion energy', self%nuclearRepulsion
  write(*,'(1x,a40," = ",f12.8)') 'Electronic energy', self%electronic      
  write(*,'(1x,a40," = ",f12.8)') 'Total energy', self%total           
  write(*,'(1x,a40," = ",f12.8)') 'Total spin Sz', self%sz              
  write(*,'(1x,a40," = ",f12.8)') 'Square of spin ', self%szz             
  write(*,'(1x,a40," = ",f12.8)') 'Frozen core energy', self%core            
  write(*,'(1x,a40," = ",f12.8)') 'SCF total energy', self%scf             
  write(*,'(1x,a40," = ",f12.8)') 'Electron-nucleus attraction', self%electronNucleus 
  write(*,'(1x,a40," = ",f12.8)') 'Electron-electron repulsion', self%electronElectron
  write(*,'(1x,a40," = ",f12.8)') 'Potential energy', self%potential       
  write(*,'(1x,a40," = ",f12.8)') 'Kinetic energy', self%kinetic         
 end subroutine print_energy_type

!==============================================
! procedure for handling two electron integrals
!==============================================

   function address(i,j,k,l)                                       
                                                                        
   integer(IK) :: i,j,k,l,ij,kl,address                                                
                                                                        
    ij = max(i,j)*(max(i,j)-1)/2 + min(i,j)                                
    kl = max(k,l)*(max(k,l)-1)/2 + min(k,l)                                
                                                                        
    address = max(ij,kl)*(max(ij,kl)-1)/2 + min(ij,kl)                     
                                                                        
  end function

  subroutine ReadTwoIntAO(TwoIntAO)
    use VarModule
! read two-electron integrals over AO's
    real(DP),         intent(inout)   :: TwoIntAO(:)

    real(DP),          allocatable  :: buffer(:)
    integer(IPgamess), allocatable  :: indexBuffer(:)

    real(DP)          :: temp

    integer(IK)       :: readStatus
    integer(IK)       :: m
    integer(IK)       :: i,j,k,l
    integer(IK)       :: bufLength, twoIntIndexBufSize, twoIntBufferSize
    integer(IPgamess) :: label, label1, label2
    integer(IPgamess) :: length, nintmx, labsiz, twoeao, iw 
    logical           :: largeLabels    

    nintmx = 15000                               ! gamess parameter controling read buffer
    labsiz = 1                                   ! gamess parameter controling the index range, for maxao > 256 set to 2  
    twoeao = twoein

    if (labsiz /= 1_IPgamess .and. labsiz /= 2_IPgamess) then
      write(*,*) 'RdTwoIntAO:  CONFUSION IN LABSIZ! '
      stop
    endif

    largeLabels = (labsiz == 2_IPgamess)

    twoIntBufferSize = int(nintmx, kind=IK)

    twoIntIndexBufSize = twoIntBufferSize
    if (largeLabels) then
      if (IPgamess == 4) twoIntIndexBufSize = 2*twoIntIndexBufSize
    else
      if (IPgamess == 8) twoIntIndexBufSize = (twoIntIndexBufSize + 1) / 2
    endif

    rewind(twoeao)

    allocate(buffer(twoIntBufferSize))
    allocate(indexBuffer(twoIntIndexBufSize))

    length = 1_IPgamess
    do while (length > 0_IPgamess)
      Read(twoeao,iostat=readStatus) length,indexBuffer,buffer 
      if (readStatus /= 0) then
        if (readStatus == 2) then
          write(*,*) 'RdTwoIntAO: ENCOUNTERED UNEXPECTED END WHILE READING TWO-ELECTRON FILE'
        else
          write(*,*) 'RdTwoIntAO: ENCOUNTERED I/O ERROR WHILE READING TWO-ELECTRON FILE'
        endif
        stop
      endif

      bufLength = abs(int(length, kind=IK))
      if (bufLength > twoIntBufferSize) stop

      do m=1, bufLength
        if (IPgamess == 4) then                  ! 32-bit Gamess integers
          if (largeLabels) then       
            label1 = indexbuffer(2*m-1) 
            label2 = indexBuffer(2*m)
            i = ishft(label1, -16_IPgamess)
            j = iand( label1,  65535_IPgamess)
            k = ishft(label2, -16_IPgamess)
            l = iand( label2,  65535_IPgamess)
          else
            label = indexBuffer(m)
            i =      ishft(label, -24_IPgamess)
            j = iand(ishft(label, -16_IPgamess), 255_IPgamess)
            k = iand(ishft(label,  -8_IPgamess), 255_IPgamess)
            l = iand(      label,                255_IPgamess)
          endif  
        else                                     ! 64-bit Gamess integers
          if (largeLabels) then
            label = indexBuffer(m)
            i = int(     ishft(label,   -48_IPgamess),                  kind=IK)
            j = int(iand(ishft(label,   -32_IPgamess), 65535_IPgamess), kind=IK)
            k = int(iand(ishft(label,   -16_IPgamess), 65535_IPgamess), kind=IK)
            l = int(iand(      label,                  65535_IPgamess), kind=IK)  
          else
            if (mod(m,2) == 0) then
              label = indexBuffer(m/2)
              i = int(iand(ishft(label, -24_IPgamess ), 255_IPgamess), kind=IK)
              j = int(iand(ishft(label, -16_IPgamess ), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label,  -8_IPgamess ), 255_IPgamess), kind=IK)
              l = int(iand(      label,                 255_IPgamess), kind=IK)
            else
              label = indexBuffer(m/2+1)
              i = int(     ishft(label, -56_IPgamess),                kind=IK)
              j = int(iand(ishft(label, -48_IPgamess), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label, -40_IPgamess), 255_IPgamess), kind=IK)
              l = int(iand(ishft(label, -32_IPgamess), 255_IPgamess), kind=IK)
            endif      
          endif
        endif

        temp = buffer(m)
        TwoIntAO(addr(i,j,k,l)) = temp
 
      enddo
    enddo
   
    deallocate(buffer)
    deallocate(indexBuffer)
  end subroutine

  subroutine ReadTwoIntMO(TwoIntMO)
    use VarModule
    real(DP),         intent(inout)   :: TwoIntMO(:)

    real(DP),          allocatable  :: buffer(:)
    integer(IPgamess), allocatable  :: indexBuffer(:)

    real(DP)          :: temp

    integer(IK)       :: readStatus
    integer(IK)       :: m
    integer(IK)       :: i,j,k,l
    integer(IK)       :: bufLength, twoIntIndexBufSize, twoIntBufferSize
    integer(IPgamess) :: label, label1, label2
    integer(IPgamess) :: length, nintmx, labsiz, twoemo, iw, ipu, is, transInt
    integer(IPgamess) :: nft11,nft12,nft13,nft14,nft15,nft16
    logical           :: largeLabels    

    nintmx = 15000
    labsiz = 1 
    twoemo = twoein

    if (labsiz /= 1_IPgamess .and. labsiz /= 2_IPgamess) then
      write(*,*) 'RdTwoIntAO:  CONFUSION IN LABSIZ! '
      stop
    endif

    largeLabels = (labsiz == 2_IPgamess)

    twoIntBufferSize = int(nintmx, kind=IK)

    twoIntIndexBufSize = twoIntBufferSize
    if (largeLabels) then
      if (IPgamess == 4) twoIntIndexBufSize = 2*twoIntIndexBufSize
    else
      if (IPgamess == 8) twoIntIndexBufSize = (twoIntIndexBufSize + 1) / 2
    endif

    rewind(twoemo)
    read(twoemo)

    allocate(buffer(twoIntBufferSize))
    allocate(indexBuffer(twoIntIndexBufSize))

    length = 1_IPgamess
    do while (length > 0_IPgamess)
      Read(twoemo,iostat=readStatus) length,indexBuffer,buffer 
      if (readStatus /= 0) then
        if (readStatus == 2) then
          write(*,*) 'RdTwoIntMO: ENCOUNTERED UNEXPECTED END WHILE READING TWO-ELECTRON FILE'
        else
          write(*,*) 'RdTwoIntMO: ENCOUNTERED I/O ERROR WHILE READING TWO-ELECTRON FILE'
        endif
        stop
      endif

      bufLength = abs(int(length, kind=IK))
      if (bufLength > twoIntBufferSize) stop

      do m=1, bufLength
        if (IPgamess == 4) then                  ! 32-bit Gamess integers
          if (largeLabels) then       
            label1 = indexbuffer(2*m-1) 
            label2 = indexBuffer(2*m)
            i = ishft(label1, -16_IPgamess)
            j = iand( label1,  65535_IPgamess)
            k = ishft(label2, -16_IPgamess)
            l = iand( label2,  65535_IPgamess)
          else
            label = indexBuffer(m)
            i =      ishft(label, -24_IPgamess)
            j = iand(ishft(label, -16_IPgamess), 255_IPgamess)
            k = iand(ishft(label,  -8_IPgamess), 255_IPgamess)
            l = iand(      label,                255_IPgamess)
          endif  
        else                                     ! 64-bit Gamess integers
          if (largeLabels) then
            label = indexBuffer(m)
            i = int(     ishft(label,   -48_IPgamess),                  kind=IK)
            j = int(iand(ishft(label,   -32_IPgamess), 65535_IPgamess), kind=IK)
            k = int(iand(ishft(label,   -16_IPgamess), 65535_IPgamess), kind=IK)
            l = int(iand(      label,                  65535_IPgamess), kind=IK)  
          else
            if (mod(m,2) == 0) then
              label = indexBuffer(m/2)
              i = int(iand(ishft(label, -24_IPgamess ), 255_IPgamess), kind=IK)
              j = int(iand(ishft(label, -16_IPgamess ), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label,  -8_IPgamess ), 255_IPgamess), kind=IK)
              l = int(iand(      label,                 255_IPgamess), kind=IK)
            else
              label = indexBuffer(m/2+1)
              i = int(     ishft(label, -56_IPgamess),                kind=IK)
              j = int(iand(ishft(label, -48_IPgamess), 255_IPgamess), kind=IK)
              k = int(iand(ishft(label, -40_IPgamess), 255_IPgamess), kind=IK)
              l = int(iand(ishft(label, -32_IPgamess), 255_IPgamess), kind=IK)
            endif      
          endif
        endif

        temp = buffer(m)
        TwoIntMO(addr(i,j,k,l)) = temp
 
      enddo
    enddo
   
    deallocate(buffer)
    deallocate(indexBuffer)
  end subroutine

  subroutine print4indexVec(fi, nb)
 !  use dmftHelpModule
  integer(IK) :: i,j,k,l,ij,kl,ir,iw,ipu,is,nb, transInt
  real(DP), intent(in)  :: fi(:)

  write(*,*)
   ij = 0
   do i = 1,nb
     do j = 1,i
       ij = ij+1
       kl = 0
       do k = 1,nb
         do l = 1,k
           kl = kl + 1
             if(ij >= kl .and. abs(fi(addr(i,j,k,l))) > 1.0e-10) write(*,'(2x,4i5,2x,f20.14)') i,j,k,l,fi(addr(i,j,k,l))
         enddo
       enddo 
     enddo 
   enddo

  end subroutine 


end module IntegralsModule
