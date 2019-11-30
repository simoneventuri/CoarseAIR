! -*-F90-*-
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
!
! Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
! 
! This program is free software; you can redistribute it and/or modify it under the terms of the 
! Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
! 
! This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU Lesser General Public License for more details. 
! 
! You should have received a copy of the GNU Lesser General Public License along with this library; 
! if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
! 
!---------------------------------------------------------------------------------------------------------------
!===============================================================================================================

! ==============================================================================================================
! Forward Propagating Uncertainties on Rates for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! This Program ...
!
!   Mandatory Input Arguments: - Nb of Translational Temperatures
!                              - Nb of Internal Temperatures
!                              - Total Number of Processors, NProc
!                              - Processor Identifier for Parallelization, iProc (1 <= iProc <= NProc)
!                              - Nb of Binned Molecules = Input%NBinnedMolecules
!                                   - Nb of Bins for BinnedMolecule(iMol),          iMol=1:NBinnedMolecules
!
! ==============================================================================================================
Program FwdPropagation

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, Six, Ten
  use Logger_Class              ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class               ,only:  Error
  use Input_Class               ,only:  Input_Type
  use LevelsContainer_Class     ,only:  LevelsContainer_Type
  use BinsContainer_Class       ,only:  BinsContainer_Type
  use Collision_Class           ,only:  Collision_Type
  use System_Class              ,only:  System_Type
  use System_Factory_Class      ,only:  System_Factory_Type
  use Molecule_Class            ,only:  Molecule_Type
  use Molecule_Factory_Class    ,only:  Molecule_Factory_Type
  use KonigPostprocessing_Class ,only:  KonigPostprocessing_Type
  use Timing_Module
  use fgsl                      ,only:   fgsl_rng_type, fgsl_rng, fgsl_rng_env_setup, fgsl_rng_alloc, fgsl_rng_set, fgsl_ran_gaussian                       !!! FGSL

  implicit none
  
  Type(Input_Type)                                           :: Input
  Type(LevelsContainer_Type)      ,dimension(:) ,allocatable :: LevelsContainer_Orig
  Type(LevelsContainer_Type)      ,dimension(:) ,allocatable :: LevelsContainer_Cut
  Type(BinsContainer_Type)        ,dimension(:) ,allocatable :: BinnedMolecule
  Type(Collision_Type)                                       :: Collision
  Type(System_Factory_Type)                                  :: System_Factory
  Class(System_Type)                            ,allocatable :: System_Adjst
  Type(Molecule_Factory_Type)                                :: Molecule_Factory
  Class(Molecule_Type)                          ,allocatable :: Molecule_Adjst
  Type(KonigPostprocessing_Type)                             :: KonigPostprocessing
  
  integer                                                     :: iMol
  integer       ,dimension(:),   allocatable                  :: MolOK
  integer                                                     :: iMol_to_Pair
  logical                                                     :: flag_Binned
  integer       ,dimension(:),   allocatable                  :: Molecule_to_Bin
  integer                                                     :: iBinnedMol
  integer                                                     :: iBins
  character(6)                                                :: iBins_char
  integer                                                     :: iBinsIn
  integer                                                     :: iBinsFn
  character(10)                                               :: SystemTemp
  character(10)                                               :: PES_Model
  character(10)                                               :: Molecules_Name
  real(rkp)                                                   :: Ttra
  integer                                                     :: iTtra
  real(rkp)     ,dimension(:)   ,allocatable                  :: Ttra_Vec
  character(5)  ,dimension(:)   ,allocatable                  :: Ttra_Vec_char
  real(rkp)                                                   :: Tint
  integer                                                     :: iTint
  character(2)                                                :: iTintChar
  character(:)                  ,allocatable                  :: iTintCharTemp
  real(rkp)                                                   :: Ecut
  integer                                                     :: NTtra
  integer                                                     :: NTraj
  integer                                                     :: iArr
  integer       ,dimension(:)   ,allocatable                  :: NTotArr_vec
  integer       ,dimension(:)   ,allocatable                  :: NArr_vec
  integer                                                     :: iP, iPTemp, jP
  character(20) ,dimension(:)   ,allocatable                  :: ReactionIn
  character(20) ,dimension(:)   ,allocatable                  :: ReactionFn
  integer       ,dimension(:)   ,allocatable                  :: TypeVec
  integer       ,dimension(:,:) ,allocatable                  :: AtomsToPairs
  integer       ,dimension(:)   ,allocatable                  :: OtherAtoms
    
  logical                                                     :: FirstTime = .true.
  integer                                                     :: iFwdProp
  character(:)                  ,allocatable                  :: ThermoFileOrig 
  character(:)                  ,allocatable                  :: ThermoFileNew
  character(:)                  ,allocatable                  :: OutputPath
  integer                                                     :: iTimeNodes
  character(:)                  ,allocatable                  :: KonigOKFile 
  logical                                                     :: KonigOKFlag
  integer       ,dimension(:)   ,allocatable                  :: iWrite
  
  integer                                                     :: iComp
  character(2)                                                :: iComp_char

  logical                                                     :: exist_flag
  character(:)                  ,allocatable                  :: FileName
  integer                                                     :: Status
  integer                                                     :: Unit
  integer                                                     :: err
  character(2)                                                :: iArr_char
  character(:)                  ,allocatable                  :: format_char
  character(:)                  ,allocatable                  :: format_charRead
  character(:)                  ,allocatable                  :: format_charReadBins
  real(rkp)                                                   :: StartTime, EndTime
  
  logical                                                     :: i_Debug      = .False.
  logical                                                     :: i_Debug_Deep = .False.
  
  type(fgsl_rng_type)                                         :: t_fgsl                                                                !!! FGSL                                    
  type(fgsl_rng)                                              :: r_fgsl 
  
  
  call CPU_Time( StartTime )
  
  call Logger%Initialize(             "FwdPropagation.log",  &                                          ! Opening the Log File using
              Status          =       'REPLACE',             &                                          ! replacing any previous log file
              Position        =       'REWIND',              &                                          ! rewinding to the top of the file
              Procedure       =       'FwdPropagation',      &                                          ! loading the calling procedure name
              Indentation     =       2           )                                                     ! and setting the initial indentation level


! ==============================================================================================================
!  READING THE PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 1, Input%NTtra_char )
  read(Input%NTtra_char, "(I10)", iostat=Status) Input%NTtra
  if (Status/=0) call Error( "Error reading the argument Input%NTtra" )
  allocate(Ttra_Vec(max(Input%NTtra,3)), Stat=Status )
  if (Status/=0) call Error( "Error allocating Ttra_Vec" )
  if ( Logger%On() ) call Logger%Write( "Allocated Ttra_Vec.")
  allocate(Ttra_Vec_char(max(Input%NTtra,3)), Stat=Status )
  if (Status/=0) call Error( "Error allocating Ttra_Vec_char" )
  if ( Logger%On() ) call Logger%Write( "Allocated Ttra_Vec_char.")
  
  call getarg( 2, Input%NTint_char )
  read(Input%NTint_char, "(I10)", iostat=Status) Input%NTint
  if (Status/=0) call Error( "Error reading the argument Input%NTint" )
  
  call getarg( 3, Input%NProc_char )
  if (trim(adjustl(Input%NProc_char)) .eq. '') then
    Input%NProc_char  = '1'
    Input%NProc       =  1
  else 
    read( Input%NProc_char, '(I3)' ) Input%NProc
  end if
  
  call getarg( 4, Input%iProc_char )
  if (trim(adjustl(Input%iProc_char)) .eq. '') then
    Input%iProc_char  = '1'
    Input%iProc       =  1
  else 
    read( Input%iProc_char, '(I3)' ) Input%iProc
  end if
! ==============================================================================================================
  

  
! ==============================================================================================================
!   INITIALIZING INPUT OBJECT
! ==============================================================================================================
  if ( Logger%On() ) call Logger%Write( "-> Calling Input%Initialize" )
  call Input%Initialize( )
  if ( Logger%On() ) call Logger%Write( "-> Calling Input%FwdProagation" )
  call Input%FwdProagation( )
  if ( Logger%On() ) call Logger%Write( "-> Done reading Input for FwdProagation" )
! ============================================================================================================== 


! ==============================================================================================================                 !!! FGSL
!   INITIALIZING RANDOM GENERATOR
! ==============================================================================================================
  t_fgsl = fgsl_rng_env_setup()
  r_fgsl = fgsl_rng_alloc (t_fgsl)
  call fgsl_rng_set(r_fgsl, int(Input%FwdSeed,8) )
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING PAIR OBJECT
! ==============================================================================================================
  if ( Logger%On() ) call Logger%Write( "Initializing the Collision%Pairs object", NewLine=.True. )
  if ( Logger%On() ) call Logger%Write( "-> Calling Collision%InitializePairs" )
  call Collision%InitializePairs( Input )
  if ( Logger%On() ) call Logger%Write( "-> Done with Collision%InitializePairs" )
! ==============================================================================================================
  

! ==============================================================================================================
!   ASSIGNING MOLECULE to PAIRS
! ==============================================================================================================
  if ( Logger%On() ) call Logger%Write( "Finding Corrispondences between Pairs and Molecules", NewLine=.True. )
  do iMol = 1,Input%NMolecules
    if ( Logger%On() ) call Logger%Write( "Molecule Nb", iMol )
    if ( Logger%On() ) call Logger%Write( "-> Calling Molecule_Factory%Define_Molecule" )
    call Molecule_Factory%Define_Molecule( Input, Molecule_Adjst, iMol )
    if ( Logger%On() ) call Logger%Write( "-> Calling Molecule%AssignMoleculesToPairs" )
    call Molecule_Adjst%AssignMoleculesToPairs( Collision, Input, iMol )
    if ( Logger%On() ) call Logger%Write( "-> Done with Molecule%AssignMoleculesToPairs" )
    deallocate(Molecule_Adjst)
    if ( Logger%On() ) call Logger%Write( "Molecule Deallocated" )
  end do
! ==============================================================================================================


! ==============================================================================================================
!    DEFINING ARRANGEMENTS
! ==============================================================================================================
  if ( Logger%On() ) call Logger%Write( "Finding Equal Pairs in the System", NewLine=.True. )
  if ( Logger%On() ) call Logger%Write( "-> Calling System_Factory%Define_System" )
  call System_Factory%Define_System( Input, System_Adjst )
  if ( Logger%On() ) call Logger%Write( "-> Calling System%AssignPairsArrangements" )
  call System_Adjst%AssignPairsArrangements( Collision, Input )
  if ( Logger%On() ) call Logger%Write( "-> Done with System%AssignPairsArrangements" )
! ==============================================================================================================


! ==============================================================================================================
!   ALLOCATING LEVELS CONTAINERS
! ==============================================================================================================
  if ( Logger%On() ) call Logger%Write( "Allocating the Molecules Original Energy Levels (containers for q.n.s etc) based on the Nb of Molecules" )
  allocate(LevelsContainer_Orig(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating LevelsContainer_Orig" )
  if ( Logger%On() ) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Orig" )

  if ( Logger%On() ) call Logger%Write( "Allocating the Molecules Cut Energy Levels (containers for q.n.s etc) based on the Nb of Molecules" )
  allocate(LevelsContainer_Cut(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating LevelsContainer_Cut" )
  if ( Logger%On() ) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Cut" )
  
  if ( Logger%On() ) call Logger%Write( "Allocating a counter for the Molecules" )
  allocate(MolOK(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating MolOK" )
  if ( Logger%On() ) call Logger%Write( "Allocated ", Input%NMolecules, " MolOK" )
  MolOK = 0
! ==============================================================================================================


! ==============================================================================================================
!   LOOKING FOR BINNED MOLECULES
! ==============================================================================================================
  if (Input%NBinnedMolecules /= 0) then 
    if ( Logger%On() ) call Logger%Write( "Allocating the Binned-Molecules (containers for bins etc.) based on the Nb of Molecules" )
    
    allocate(BinnedMolecule(Input%NBinnedMolecules), Stat=Status )
    if (Status/=0) call Error( "Error allocating BinnedMolecule" )
    if ( Logger%On() ) call Logger%Write( "Allocated ", Input%NBinnedMolecules, " BinnedMolecule" )
    
    allocate(Molecule_to_Bin(Input%NBinnedMolecules), Stat=Status )
    if (Status/=0) call Error( "Error allocating Molecule_to_Bin" )
    if ( Logger%On() ) call Logger%Write( "Allocated ", Input%NBinnedMolecules, " BinnedMolecule" )
    
    if ( Logger%On() ) call Logger%Write( "Finding Corrispondences between Pairs and Molecules", NewLine=.True. )
    do iBinnedMol = 1,Input%NBinnedMolecules
      if ( Logger%On() ) call Logger%Write( "Binned Molecule Nb", iBinnedMol )
      

      !==============================================================================================================
      ! INITIALIZING BINNED MOLECULE
      !==============================================================================================================
      if (i_Debug) call Logger%Write( "Initialize Bins for the Molecule Nb", iMol )
      call BinnedMolecule(iMol)%InitializeBins( Input, iMol, i_Debug=.True.  ) 
    
      iMol = BinnedMolecule(iBinnedMol)%To_Molecule
      ! ============================================================================================================== 
      

      !==============================================================================================================
      ! ASSIGNING BINNING MOLECULES TO PAIRS
      !==============================================================================================================
      if ( Logger%On() ) call Logger%Write( "-> Calling Molecule_Factory%Define_Molecule" )
      call Molecule_Factory%Define_Molecule( Input, Molecule_Adjst, iMol )
      if ( Logger%On() ) call Logger%Write( "-> Calling Molecule%AssignBinnedMoleculesToPairs" )
      call Molecule_Adjst%AssignBinnedMoleculesToPairs( Collision, Input, iMol )
      deallocate(Molecule_Adjst)
      if ( Logger%On() ) call Logger%Write( "Molecule Deallocated" )
      ! ============================================================================================================== 
      
      
    end do
    if ( Logger%On() ) call Logger%Write( "-> Done with Molecule%AssignBinnedMoleculesToPairs" )
    
  end if
! ==============================================================================================================
  

  if ( Logger%On() ) call Logger%Write( "Iterating on all the Molecules in the Chemical System" )
  do iMol = 1,Input%NMolecules
    if ( Logger%On() ) call Logger%Write( "Molecule Nb", iMol )
    
    
    ! ==============================================================================================================
    !   READING MOLECULE ENERGY LEVELS
    ! ==============================================================================================================
    if ( Logger%On() ) call Logger%Write( "Reading the Original Levels List file for Molecule Nb", iMol )
    FileName = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
               trim(adjustl(Input%LevelsFileName(iMol)))
    if ( Logger%On() ) call Logger%Write( "Reading File: ", FileName )
    call LevelsContainer_Orig(iMol)%Initialize( Input, iMol, FileName )

    if ( Logger%On() ) call Logger%Write( "Reading the Cut Levels List file for Molecule Nb", iMol )
    FileName  = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/levels_cut.inp'
    if ( Logger%On() ) call Logger%Write( "Reading File: ", FileName )
    call LevelsContainer_Cut(iMol)%Initialize( Input, iMol, FileName )
    ! ==============================================================================================================
    
    
    ! ==============================================================================================================
    !   CHECKING IF CURRENT MOLECULE IS BINNED
    ! ==============================================================================================================
    flag_Binned           = .FALSE.
    if (Input%NBinnedMolecules /= 0) then 
    
      Molecule_to_Bin(iMol) =    0
      do iBinnedMol = 1,Input%NBinnedMolecules
        if ( Input%BinnedMolecules_Name(iBinnedMol) .eq. Input%Molecules_Name(iMol) ) then
          flag_Binned = .true.
          Molecule_to_Bin(iMol) = iBinnedMol
          if ( Logger%On() ) call Logger%Write( "Current Molecule Nb", iMol, "(", trim(adjustl(Input%Molecules_Name(iMol))), ") is the Binned Molecule Nb", iBinnedMol )
        end if
      end do
      
      if (flag_Binned) then
      
      
        ! ==============================================================================================================
        !   SORTING LEVELS AND COMPUTING BINS PARTITION FUNCTIONS
        ! ==============================================================================================================
        if (trim(Input%BSortMethod(iMol)) .eq. "State-Specific" ) then 
          if (i_Debug) call Logger%Write( "Sorting the Molecule Nb", iMol, " by State Specific" )
          call BinnedMolecule(iMol)%SortByStateSpecific( Input, LevelsContainer_Cut(iMol), iMol, .false., i_Debug=i_Debug_Deep )
        elseif (trim(Input%BSortMethod(iMol)) .eq. "Vib-Specific" ) then 
          if (i_Debug) call Logger%Write( "Sorting the Molecule Nb", iMol, " by Vibrational Specific" )
          call BinnedMolecule(iMol)%SortByVibSpecific( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_Deep )
        elseif (trim(Input%BSortMethod(iMol)) .eq. "RoVib-CG" ) then 
          if (i_Debug) call Logger%Write( "Sorting the Molecule Nb", iMol, " by increasing Ro-Vibrational Energy" )
          call BinnedMolecule(iMol)%SortByRoVibEnergy( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_Deep )
        elseif (trim(Input%BSortMethod(iMol)) .eq. "Hybrid" ) then
          if (i_Debug) call Logger%Write( "Sorting the Molecule Nb", iMol, " by Hybrid Mthd" )
          call BinnedMolecule(iMol)%SortByHybrid( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_Deep )
        else
          if (i_Debug) call Logger%Write( "ERROR: Sorting Method Unknown! Input%BSortMethod(i)", trim(Input%BSortMethod(iMol))) !!! @TODO: Adding New Sorting Methods
          stop ( "ERROR: Sorting Method Unknown!" )    
        end if
        
        if (i_Debug) call Logger%Write( "Reading Bins Energies and Partition Functions" )
        call BinnedMolecule(iMol)%ReadPartFunEnergy( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_Deep )
        ! ==============================================================================================================
      
      
      end if

    end if
    ! ==============================================================================================================


  end do


  !==============================================================================================================
  ! COUNTING THE RATES "CONTAINERS" REQUIRED
  !==============================================================================================================
  !
  ! Example: System=CO2; CO Binned Molecule (10Bins) and O2 Molecule -> Pair1=CO, Pair2=CO, Pair3=O2.
  ! 
  !                           [ Diss | 16 -> CO Bin1 | 17 -> CO Bin1 | 16 -> CO Bin2 | ... | 16 -> CO Bin1 | 17 -> CO Bin1 | 16 -> CO Bin2 | ... | 16 -> O2 | 17 -> O2 ]
  !                       |      |         \         +       |       +        /      +  /        \         +       |       +        /      +  /        |          |
  !                       |      |          \                |               /         /          \                |               /         /         |          |
  !   NTotArr_vec =   [   0   |  1   |        NBins of CO * 2 ( 16 -> and 17 -> )          |        NBins of CO * 2 ( 16 -> and 17 -> )          |          2          ] 
  ! 
  !
  !                           [ Diss | 16 -> CO Bin1 | 17 -> CO Bin1 | 16 -> CO Bin2 | ... | 16 -> CO Bin1 | 17 -> CO Bin1 | 16 -> CO Bin2 | ... | 16 -> O2 | 17 -> O2 ]
  !                       |      |         \         +       |       +        /      +  /        \         +       |       +        /      +  /        |          |
  !                       |      |          \                |               /         /          \                |               /         /         |          |
  !   NTArr_vec   =   [   0   |  1   |                 NBins of CO                         |                   NBins of CO                       |          1          ] 
  !
  allocate(NArr_vec(Collision%NPairs+2), Stat=Status )
  if (Status/=0) call Error( "Error allocating NArr_vec" )
  if ( Logger%On() ) call Logger%Write( "Allocated NArr_vec" )
  NArr_vec = 0
  
  allocate(NTotArr_vec(Collision%NPairs+2), Stat=Status )
  if (Status/=0) call Error( "Error allocating NTotArr_vec" )
  if ( Logger%On() ) call Logger%Write( "Allocated NTotArr_vec" )
  NTotArr_vec = 0

  NTotArr_vec(1) = 0
  NArr_vec(1)    = 0
  NTotArr_vec(2) = 1
  NArr_vec(2)    = 1
  do iP = 3,Collision%NPairs+2
    if (Collision%Pairs(iP-2)%To_BinnedMolecule /= 0) then
      NTotArr_vec(iP) = NTotArr_vec(iP-1) + int(Input%NBins(Collision%Pairs(iP-2)%To_BinnedMolecule)*2.0_rkp)
      NArr_vec(iP)    = NArr_vec(iP-1)    + Input%NBins(Collision%Pairs(iP-2)%To_BinnedMolecule)
    else
      NTotArr_vec(iP) = NTotArr_vec(iP-1) + 2
      NArr_vec(iP)    = NArr_vec(iP-1)    + 1
    end if
  end do
  
  if ( Logger%On() ) call Logger%Write( "NTotArr_vec = ", NTotArr_vec )
  if ( Logger%On() ) call Logger%Write( "NArr_vec    = ", NArr_vec    )
  ! ==============================================================================================================
  

  format_charRead        =  "(1X,A5,5X,A5,5X,A5,5X,f10.2,f10.2,f10.2,I10,*(es20.10))"
  format_charReadBins    =  "(5X,A5,4X,A5,7X,A5,I10,f10.2,f10.2,f10.2,I10,400es20.10)"
  

  do iMol = 1,Input%NMolecules
    if ( Logger%On() ) call Logger%Write( "Molecule Nb", iMol )
    
    
    ! ==============================================================================================================
    !   CREATING FORMAT for ARRHENIUS COEFFICIENTS OUTPUT FILE
    ! ==============================================================================================================
    select case (Input%NAtoms)       ! Setting the index of the two atoms associated to each pair
    case (3)
      allocate(AtomsToPairs(2,3))
      AtomsToPairs(:,1) = [1,2]
      AtomsToPairs(:,2) = [1,3]
      AtomsToPairs(:,3) = [2,3]
    end select
    
    select case (Input%NAtoms)       ! Setting the index of the two atoms associated to each pair
    case (3)
      allocate(OtherAtoms(3))
      OtherAtoms(1) = 3
      OtherAtoms(2) = 2
      OtherAtoms(3) = 1
    end select
    
    do iPTemp = 1,Collision%NPairs
      if (Collision%Pairs(iPTemp)%To_Molecule == iMol) then
        iP=iPTemp
        exit
      end if
    end do
    if (flag_Binned .eqv. .false.) then
      allocate(ReactionIn(1))
      ReactionIn(1) = trim(adjustl(Input%Molecules_Name(iMol))) // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(iP))))
    else
      allocate(ReactionIn(Input%NBins(Collision%Pairs(iP)%To_Molecule)))
      iArr = 1
      do iBins = 1,Input%NBins(Collision%Pairs(iP)%To_Molecule)
        write(iBins_char,'(I6)') iBins
        ReactionIn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(iP)%To_Molecule))) // "(" // trim(adjustl(iBins_char)) // ")" // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(iP))))
        iArr = iArr+1
      end do
    end if
    
    allocate(ReactionFn(NArr_vec(Collision%NPairs+2)), Stat=Status )
    if (Status/=0) call Error( "Error allocating ReactionFn" )
    if ( Logger%On() ) call Logger%Write( "Allocated ReactionFn with dimension = (", NArr_vec(Collision%NPairs+2),  ")" )
    
    allocate(TypeVec(NArr_vec(Collision%NPairs+2)), Stat=Status )
    if (Status/=0) call Error( "Error allocating TypeVec" )
    if ( Logger%On() ) call Logger%Write( "Allocated ReactionFn with dimension = (", NArr_vec(Collision%NPairs+2),  ")" )
    TypeVec = 0
    
    iArr = 1
    ReactionFn(iArr) = trim(adjustl(Input%AtomsName(1))) // "+" // trim(adjustl(Input%AtomsName(2))) // "+" // trim(adjustl(Input%AtomsName(3)))
    iArr = 2
    do jP = 1,Collision%NPairs 
      if (Collision%Pairs(jP)%To_BinnedMolecule == 0) then
        if (Collision%Pairs(jP)%To_Molecule == 0) then
          if (trim(adjustl(Input%AtomsName(AtomsToPairs(1,jP)))) == trim(adjustl(Input%AtomsName(AtomsToPairs(2,jP))))) then
            ReactionFn(iArr) = trim(adjustl(Input%AtomsName(AtomsToPairs(1,jP)))) // "2+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
          else
            ReactionFn(iArr) = trim(adjustl(Input%AtomsName(AtomsToPairs(1,jP)))) // trim(adjustl(Input%AtomsName(AtomsToPairs(2,jP)))) // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
          end if
          TypeVec(iArr) = 2
        else
          ReactionFn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(jP)%To_Molecule))) // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
          if (MolOK(Collision%Pairs(jP)%To_Molecule) == 0) then
            if (Collision%Pairs(jP)%To_Molecule == Collision%Pairs(iP)%To_Molecule ) then
              TypeVec(iArr) = 1
            else 
              TypeVec(iArr) = 2
            end if
          else
            TypeVec(iArr) = -1
          end if
          iArr = iArr+1
          MolOK(Collision%Pairs(jP)%To_Molecule) = 1
        end if
      else
        do iBins = 1,Input%NBins(Collision%Pairs(jP)%To_Molecule)
          write(iBins_char,'(I6)') iBins
          ReactionFn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(jP)%To_Molecule))) // "(" // trim(adjustl(iBins_char)) // ")" // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
          if (MolOK(Collision%Pairs(jP)%To_Molecule) == 0) then
            if (Collision%Pairs(jP)%To_Molecule == Collision%Pairs(iP)%To_Molecule ) then
              TypeVec(iArr) = 1
            else 
              TypeVec(iArr) = 2
            end if
          else
            TypeVec(iArr) = -1
          end if
          iArr = iArr+1
        end do
        MolOK(Collision%Pairs(jP)%To_Molecule) = 1
      end if
    end do
    deallocate(AtomsToPairs)
    ! ==============================================================================================================
  
    
    if (flag_Binned .eqv. .false.) then
      if ( Logger%On() ) call Logger%Write( "Current Molecule Nb", iMol, "(", Input%Molecules_Name(iMol), ") is NOT binned" )
      
      
      ! ==============================================================================================================
      !   ALLOCATING RATES and ARRHENIUS COEFFICIENTS
      ! ==============================================================================================================
      allocate(LevelsContainer_Orig(iMol)%RateConst_Final(NArr_vec(Collision%NPairs+2)), Stat=Status )
      if (Status/=0) call Error( "Error allocating LevelsContainer_Orig(iMol)%RateConst_Final" )
      if ( Logger%On() ) call Logger%Write( "Allocated LevelsContainer_Orig(iMol)%RateConst_Final with dimension = (", NArr_vec(Collision%NPairs+2), ")" )
      LevelsContainer_Orig(iMol)%RateConst_Final = Zero
      
      allocate(LevelsContainer_Orig(iMol)%RateConst_Sigma(NArr_vec(Collision%NPairs+2)), Stat=Status )
      if (Status/=0) call Error( "Error allocating LevelsContainer_Orig(iMol)%RateConst_Sigma" )
      if ( Logger%On() ) call Logger%Write( "Allocated LevelsContainer_Orig(iMol)%RateConst_Sigma with dimension = (", NArr_vec(Collision%NPairs+2), ")" )
      LevelsContainer_Orig(iMol)%RateConst_Sigma = Zero
      
      allocate(LevelsContainer_Orig(iMol)%RateConst_Arr(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
      if (Status/=0) call Error( "Error allocating LevelsContainer_Orig(iMol)%RateConst_Arr" )
      if ( Logger%On() ) call Logger%Write( "Allocated LevelsContainer_Orig(iMol)%RateConst_Arr with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
      LevelsContainer_Orig(iMol)%RateConst_Arr = Zero
      
      allocate(LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
      if (Status/=0) call Error( "Error allocating LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr" )
      if ( Logger%On() ) call Logger%Write( "Allocated LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
      LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr = Zero
      
      allocate(LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
      if (Status/=0) call Error( "Error allocating LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd" )
      if ( Logger%On() ) call Logger%Write( "Allocated LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
      LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd = Zero
      
      allocate(LevelsContainer_Orig(iMol)%CArr(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
      if (Status/=0) call Error( "Error allocating LevelsContainer_Orig%Bin(iBins)%CArr" )
      if ( Logger%On() ) call Logger%Write( "Allocated LevelsContainer_Orig%Bin(iBins)%CArr with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
      LevelsContainer_Orig(iMol)%CArr = Zero
      ! ==============================================================================================================
      
      
      ! ==============================================================================================================
      !   READING RATES
      ! ==============================================================================================================
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/Rate-Constants.dat'
      if ( Logger%On() ) call Logger%Write( "Reading File: ", FileName )
      open( File=FileName, Unit=101, status='OLD', iostat=Status )
      read(101,*)
      
      
      ! ==============================================================================================================
      !   READING RATES St DEVIATIONS
      ! ==============================================================================================================
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/Rate-Constants-Sigma.dat'
      if ( Logger%On() ) call Logger%Write( "Reading File: ", FileName )
      open( File=FileName, Unit=111, status='OLD', iostat=Status )
      read(111,*)

      
      do iTint = 1,Input%NTint
        if ( Logger%On() ) call Logger%Write( "Internal Temperature Nb", iTint )                                                  !!! WARNING: Possible Bug when running Rates not Thermal 

        LevelsContainer_Orig(iMol)%RateConst_Arr   = Zero
  
        Ttra_Vec = Zero
        do iTtra = 1,Input%NTtra
          if ( Logger%On() ) call Logger%Write( "Translational Temperature Nb", iTtra )
        
        
          read(121,format_charRead, iostat=Status) SystemTemp, PES_Model, Molecules_Name, Ttra, Tint, Ecut, NTraj, LevelsContainer_Orig(iMol)%RateConst_Final
          read(122,format_charRead, iostat=Status) SystemTemp, PES_Model, Molecules_Name, Ttra, Tint, Ecut, NTraj, LevelsContainer_Orig(iMol)%RateConst_Sigma
          
          
          ! ==============================================================================================================
          !   CREATING VECTORS OF TEMPERATURES and RATES
          ! ==============================================================================================================
          Input%Tint = Tint
          Input%Ttra = Ttra
          
          if (Input%NTtra == 1) then
            Ttra_Vec(1) = Input%Ttra / Two
            Ttra_Vec(2) = Input%Ttra
            Ttra_Vec(3) = Input%Ttra * Two
            LevelsContainer_Orig(iMol)%RateConst_Arr(:,2) = LevelsContainer_Orig(iMol)%RateConst_Final(:)
            LevelsContainer_Orig(iMol)%RateConst_Arr(:,1) = LevelsContainer_Orig(iMol)%RateConst_Arr(:,2) / Ten
            LevelsContainer_Orig(iMol)%RateConst_Arr(:,3) = LevelsContainer_Orig(iMol)%RateConst_Arr(:,2) * Ten
            LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,2) = LevelsContainer_Orig(iMol)%RateConst_Sigma(:)
            LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,1) = LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,2) / Ten
            LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,3) = LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,2) * Ten
          elseif (Input%NTtra == 2) then
            Ttra_Vec(iTtra) = Input%Ttra
            LevelsContainer_Orig(iMol)%RateConst_Arr(:,iTtra) = LevelsContainer_Orig(iMol)%RateConst_Final(:)
            LevelsContainer_Orig(iMol)%RateConst_Arr(:,3)     = LevelsContainer_Orig(iMol)%RateConst_Arr(:,2) * Ten
            LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,iTtra) = LevelsContainer_Orig(iMol)%RateConst_Sigma(:)
            LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,3)     = LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,2) * Ten
            Ttra_Vec(3) = Ttra_Vec(2) * Two
          else
            Ttra_Vec(iTtra) = Input%Ttra
            LevelsContainer_Orig(iMol)%RateConst_Arr(:,iTtra) = LevelsContainer_Orig(iMol)%RateConst_Final(:)
            LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(:,iTtra) = LevelsContainer_Orig(iMol)%RateConst_Sigma(:)
          end if
          ! ==============================================================================================================
          
          
        end do   ! iTtra
        
      end do   ! iTint

    else  
      if ( Logger%On() ) call Logger%Write( "Current Molecule Nb", iMol, "(", Input%Molecules_Name(iMol), ") is binned" )
      
      
      ! ==============================================================================================================
      !   ALLOCATING RATES and ARRHENIUS COEFFICIENTS
      ! ==============================================================================================================
      do iBins= 1,Input%NBins(iMol)

        allocate(BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(NArr_vec(Collision%NPairs+2)), Stat=Status )
        if (Status/=0) call Error( "Error allocating BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final" )
        if ( Logger%On() ) call Logger%Write( "Allocated BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final with dimension = (", NArr_vec(Collision%NPairs+2), ")" )
        BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final = Zero
        
        allocate(BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma(NArr_vec(Collision%NPairs+2)), Stat=Status )
        if (Status/=0) call Error( "Error allocating BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma" )
        if ( Logger%On() ) call Logger%Write( "Allocated BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma with dimension = (", NArr_vec(Collision%NPairs+2), ")" )
        BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma = Zero
        
        allocate(BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
        if (Status/=0) call Error( "Error allocating BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr" )
        if ( Logger%On() ) call Logger%Write( "Allocated BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
        BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr = Zero
        
        allocate(BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
        if (Status/=0) call Error( "Error allocating BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr" )
        if ( Logger%On() ) call Logger%Write( "Allocated BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
        BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr = Zero
        
        allocate(BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr_Rnd(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
        if (Status/=0) call Error( "Error allocating BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr_Rnd" )
        if ( Logger%On() ) call Logger%Write( "Allocated BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr_Rnd with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
        BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr_Rnd = Zero
        
        allocate(BinnedMolecule(iMol)%Bin(iBins)%CArr(NArr_vec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
        if (Status/=0) call Error( "Error allocating BinnedMolecule(iMol)%Bin(iBins)%CArr" )
        if ( Logger%On() ) call Logger%Write( "Allocated BinnedMolecule(iMol)%Bin(iBins)%CArr with dimension = (", NArr_vec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
        BinnedMolecule(iMol)%Bin(iBins)%CArr = Zero
        
      end do
      ! ==============================================================================================================
      
      
      ! ==============================================================================================================
      !   READING BINS RATES
      ! ==============================================================================================================
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/Bins-Rate-Constants.dat'
      if ( Logger%On() ) call Logger%Write( "Reading File: ", FileName )
      open( File=FileName, Unit=121, status='OLD', iostat=Status )
      read(121,*)

      
      ! ==============================================================================================================
      !   READING BINS RATES St DEVIATIONS
      ! ==============================================================================================================
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/Bins-Rate-Constants-Sigma.dat'
      if ( Logger%On() ) call Logger%Write( "Reading File: ", FileName )
      open( File=FileName, Unit=122, status='OLD', iostat=Status )
      read(122,*)
      

      do iTint = 1,Input%NTint                                                                                                    !!! WARNING: Possible Bug when running Rates not Thermal 

        do iBins=1,Input%NBins(iMol)
          BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr   = Zero
        end do
        
        Ttra_Vec = Zero
        do iTtra = 1,Input%NTtra
        
          iBins = 1
          do while (iBins < Input%NBins(iMol))
          
          
            read(121,format_charReadBins, iostat=Status) SystemTemp, PES_Model, Molecules_Name, iBins, Ttra, Tint, Ecut, NTraj, BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final
            read(122,format_charReadBins, iostat=Status) SystemTemp, PES_Model, Molecules_Name, iBins, Ttra, Tint, Ecut, NTraj, BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma
            

            ! ==============================================================================================================
            !   CREATING VECTORS OF TEMPERATURES and RATES
            ! ==============================================================================================================

            Input%Tint = Tint
            Input%Ttra = Ttra
            if (Input%NTtra == 1) then
              Ttra_Vec(1) = Input%Ttra / Two
              Ttra_Vec(2) = Input%Ttra
              Ttra_Vec(3) = Input%Ttra * Two
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,2) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(:)
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,1) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,2) / Ten
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,3) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,2) * Ten
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,2) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma(:)
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,1) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,2) / Ten
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,3) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,2) * Ten
            elseif (Input%NTtra == 2) then
              Ttra_Vec(iTtra) = Input%Ttra
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,iTtra) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(:)
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,3)     = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,2) * Ten
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,iTtra) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma(:)
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,3)     = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,2) * Ten
              Ttra_Vec(3) = Ttra_Vec(2) * Two
            else
              Ttra_Vec(iTtra) = Input%Ttra
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(:,iTtra) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(:)
              BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma_Arr(:,iTtra) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Sigma(:)
            end if
            ! ==============================================================================================================
            
            
          end do   ! iBins
                        
        end do   ! iTtra
        
      end do ! iTint
      
      close(121)
      close(122)
      
      
      ! ==============================================================================================================
      !   COPYING THERMO FILE
      ! ==============================================================================================================
      do iTtra = 1,Input%NTtra
  
        write(Ttra_Vec_char(iTtra),'(I5)') int(Ttra_Vec(iTtra))
        call system( 'mkdir -p ./T_' // trim(adjustl(Ttra_Vec_char(iTtra))) )
        call system( 'mkdir -p ./T_' // trim(adjustl(Ttra_Vec_char(iTtra))) // '/output' )
      
        FileName = trim(adjustl(Input%OutputDir)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                   trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(Ttra_Vec_char(iTtra))) // '.dat'
                   
        open( File=FileName, Unit=127, status='OLD', iostat=Status )
        if (Status/=0) call Error( "Error opening file: " // FileName )   
        read(127,*)                    
      
        ThermoFileOrig = trim(adjustl(Input%KonigDtbPath)) // '/thermo/' // trim(adjustl(Input%Molecules_Name(iMol))) // '_Format'
        ThermoFileNew  = './database/thermo/' // trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Ttra_Vec_char(iTtra)))
        call system( 'scp ' // trim(adjustl(ThermoFileOrig)) // ' ' // trim(adjustl(ThermoFileNew)) )
        open( File=ThermoFileNew, Unit=142, status='OLD', position="append", action="write", iostat=Status ) 
        write(142,'(A)') 'NB_ENERGY_LEVELS = ' // trim(adjustl(Input%NBins_char(iMol)))
        
        do iBins=1,Input%NBins(iMol)
          read(127,'(2E14.6)')   BinnedMolecule(iMol)%Bin(iBins)%Q,  BinnedMolecule(iMol)%Bin(iBins)%ToteinteV
          write(142,'(2(E12.6, 2X))')  BinnedMolecule(iMol)%Bin(iBins)%Q,  BinnedMolecule(iMol)%Bin(iBins)%ToteinteV
        end do
        
        close(127)
        close(142)
      
      end do
      ! ==============================================================================================================  
      
      
    end if
      
  end do   ! iMol
  

  ! ==============================================================================================================
  !   CREAGING VECTORS for EVALUATING KONIG'S OUTPUT
  ! ==============================================================================================================
  allocate(KonigPostprocessing%Ttra(Input%NTtra), Stat=Status )
  if (Status/=0) call Error( "Error allocating KonigPostprocessing%Ttra" )
  if ( Logger%On() ) call Logger%Write( "Allocated KonigPostprocessing%Ttra with dimension = (", Input%NTtra, ")" )
  
  allocate(iWrite(Input%NTtra), Stat=Status )
  if (Status/=0) call Error( "Error allocating iWrite" )
  if ( Logger%On() ) call Logger%Write( "Allocated iWrite with dimension = (", Input%NTtra, ")" )
  iWrite(:) = 0
  
  do iTtra = 1,Input%NTtra
  
    call KonigPostprocessing%EvenlySpacedVector( Input%FwdTimeMin, Input%FwdTimeMax, Input%NTimeNodes, Input%TimeScale, KonigPostprocessing%Ttra(iTtra)%TimeGrid)
  
    if (trim(adjustl(Input%FWD_MolFractions)) == 'yes') then
      call KonigPostprocessing%Allocate2DOutputToAnalyze( Input%NComp, Input, Input%NTimeNodes, Input%NMolFractionsBins, KonigPostprocessing%Ttra(iTtra)%XAtNodes, KonigPostprocessing%Ttra(iTtra)%XSum, &
                                       KonigPostprocessing%Ttra(iTtra)%XSqSum, KonigPostprocessing%Ttra(iTtra)%XHist )

      call KonigPostprocessing%EvenlySpacedVector( Input%MolFractionsMin, Input%MolFractionsMax, Input%NMolFractionsBins+1, 'lin', KonigPostprocessing%Ttra(iTtra)%XBinsExtremes)
    end if
    
    if (trim(adjustl(Input%FWD_Temperatures)) == 'yes') then
      call KonigPostprocessing%Allocate2DOutputToAnalyze( Input%NBinnedMolecules, Input, Input%NTimeNodes, Input%NTemperaturesBins, KonigPostprocessing%Ttra(iTtra)%TIntAtNodes, KonigPostprocessing%Ttra(iTtra)%TIntSum, &
                                         KonigPostprocessing%Ttra(iTtra)%TIntSqSum, KonigPostprocessing%Ttra(iTtra)%TIntHist )
 
      call KonigPostprocessing%EvenlySpacedVector( Input%TemperaturesMin, Input%TemperaturesMax, Input%NTemperaturesBins+1, 'lin', KonigPostprocessing%Ttra(iTtra)%TIntBinsExtremes)
    end if
    
    if (trim(adjustl(Input%FWD_Populations)) == 'yes') then
      call KonigPostprocessing%Allocate2DOutputToAnalyze( Input%NBins(1), Input, Input%NTimeNodes, Input%NPopulationsBins, KonigPostprocessing%Ttra(iTtra)%PopAtNodes, KonigPostprocessing%Ttra(iTtra)%PopSum, &
                                       KonigPostprocessing%Ttra(iTtra)%PopSqSum, KonigPostprocessing%Ttra(iTtra)%PopHist )

      call KonigPostprocessing%EvenlySpacedVector( Input%PopulationsMin, Input%PopulationsMax, Input%NPopulationsBins+1, 'log', KonigPostprocessing%Ttra(iTtra)%PopBinsExtremes)

    end if
  end do
  ! ==============================================================================================================
  

  ! ==============================================================================================================
  !   LOOP ON FWD PROPAGATION
  ! ==============================================================================================================
  FirstTime   = .true.
  do iFwdProp = 1,Input%NFwdPropProc
  
    write(*,*) iFwdProp
        
    ! ==============================================================================================================
    !   OPENING KINETIC FILE
    ! ==============================================================================================================
    FileName = './database/kinetics/'// trim(adjustl(Input%System))
    if ( Logger%On() ) call Logger%Write( "Writing File: ", FileName )
    open( File=FileName, Unit=102, status='REPLACE', iostat=Status )
    write(102,*) 'Units=cm^3/s'
  
  
    do iMol = 1,Input%NMolecules
      if ( Logger%On() ) call Logger%Write( "Molecule Nb", iMol )
        
      if (iMol == Collision%Pairs(1)%To_BinnedMolecule) then
      
      
        if (flag_Binned .eqv. .false.) then
          if ( Logger%On() ) call Logger%Write( "Current Molecule Nb", iMol, "(", Input%Molecules_Name(iMol), ") is NOT binned" )
                   
          LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd = Zero
          LevelsContainer_Orig(iMol)%CArr              = Zero
          
          do iBinsFn = 1,NArr_vec(Collision%NPairs+2)
          
          
            ! ==============================================================================================================
            !   RANDOMLY GENERATING RATES
            ! ==============================================================================================================
            do iTtra = 1,Input%NTtra
            
              if (LevelsContainer_Orig(iMol)%RateConst_Arr(iBinsFn,iTtra) == Zero) then
              
                LevelsContainer_Orig(iMol)%RateConst_Arr(iBinsFn,iTtra) = Input%MinForRates
              
              end if
            
              if (LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(iBinsFn,iTtra) /= Zero) then
                                                                                                                                  !!! FGSL
                LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd(iBinsFn,iTtra) = LevelsContainer_Orig(iMol)%RateConst_Arr(iBinsFn,iTtra) + fgsl_ran_gaussian( r_fgsl, LevelsContainer_Orig(iMol)%RateConst_Sigma_Arr(iBinsFn,iTtra) )

              else 
            
                LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd(iBinsFn,iTtra) = LevelsContainer_Orig(iMol)%RateConst_Arr(iBinsFn,iTtra)
                
              end if
          
            end do   ! iTtra
            ! ==============================================================================================================
            
            
            ! ==============================================================================================================      !!! GSL
            !   COMPUTING ARRHENIUS COEFFICIENTS
            ! ==============================================================================================================
            call fit(int(max(Input%NTtra,3)), Ttra_Vec, LevelsContainer_Orig(iMol)%RateConst_Arr_Rnd(iBinsFn,:), &
                     LevelsContainer_Orig(iMol)%CArr(iBinsFn,1), LevelsContainer_Orig(iMol)%CArr(iBinsFn,2), LevelsContainer_Orig(iMol)%CArr(iBinsFn,3), err)
            
            if (err /= 0) then
              LevelsContainer_Orig(iMol)%CArr(iBinsFn,:) = Zero
            end if
            ! ==============================================================================================================
            
            
            ! ==============================================================================================================
            !   WRITING ARRHENIUS COEFFICIENTS
            ! ==============================================================================================================
            if (iBinsFn == 1) then
              format_char =  "'" // trim(adjustl(ReactionIn(1))) // "=" // trim(adjustl(ReactionFn(iBinsFn))) // ":', SP,  ES11.4, ',', SP,  ES11.4, ',', SP,  ES11.4, ',2'"
            else
              format_char =  "'" // trim(adjustl(ReactionIn(1))) // "=" // trim(adjustl(ReactionFn(iBinsFn))) // ":', SP,  ES11.4, ',', SP,  ES11.4, ',', SP,  ES11.4, ',6'"
            end if
            format_char = '(' // trim(adjustl(format_char)) // ')'
            
            if ( TypeVec(iBinsFn) == 1 ) then
              write(102,trim(adjustl(format_char))) LevelsContainer_Orig(iMol)%CArr(iBinsFn,:)
            elseif ( TypeVec(iBinsFn) == 1 ) then
              write(102,trim(adjustl(format_char))) LevelsContainer_Orig(iMol)%CArr(iBinsFn,:)
            elseif (TypeVec(iBinsFn) == 0) then
              write(102,trim(adjustl(format_char))) LevelsContainer_Orig(iMol)%CArr(iBinsFn,:)
            elseif (TypeVec(iBinsFn) == 2) then
              write(102,trim(adjustl(format_char))) LevelsContainer_Orig(iMol)%CArr(iBinsFn,:)
            end if
            ! ==============================================================================================================
            
            
          end do   ! iBinsFn
          
        
        else  
          if ( Logger%On() ) call Logger%Write( "Current Molecule Nb", iMol, "(", Input%Molecules_Name(iMol), ") is binned" )
        

          do iBinsIn = 1,Input%NBins(iMol)
          
            BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr_Rnd = Zero
            BinnedMolecule(iMol)%Bin(iBinsIn)%CArr              = Zero
          
            do iBinsFn = 1,NArr_vec(Collision%NPairs+2)
              
              
              ! ==============================================================================================================
              !   RANDOMLY GENERATE BINS RATES
              ! ============================================================================================================== 
              do iTtra = 1,Input%NTtra
                  
                if (BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr(iBinsFn,iTtra) < Input%MinForRates) then 
                
                  BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr_Rnd(iBinsFn,iTtra) = Input%MinForRates           
                
                else

                  do while (BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr_Rnd(iBinsFn,iTtra) <= Input%MinForRates)
                
                    if (BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr(iBinsFn,iTtra) /= Zero) then
                                                                                                                                      !!! FGSL !!!!!!!!!!!!!!!!!!!!!
                      BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr_Rnd(iBinsFn,iTtra) = dexp( dlog(BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr(iBinsFn,iTtra)) + fgsl_ran_gaussian(r_fgsl, dlog(BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr(iBinsFn,iTtra)) * 0.1d0) )

                    else 
                  
                      BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr_Rnd(iBinsFn,iTtra) = BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr(iBinsFn,iTtra)
                      
                    end if
                    
                  end do
                  
                  
                end if
              
              end do   ! iTtra
              ! =============================================================================================================
            
              
              
              ! ==============================================================================================================    !!! GSL
              !   COMPUTING ARRHENIUS COEFFICIENTS
              ! ==============================================================================================================
              call fit(int(max(Input%NTtra,3)), Ttra_Vec, BinnedMolecule(iMol)%Bin(iBinsIn)%RateConst_Arr_Rnd(iBinsFn,:), &
                       BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,1), BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,2), BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,3), err)
                       
              if (err /= 0) then
                BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,:) = Zero
              end if
              ! ==============================================================================================================
            

              ! ==============================================================================================================
              !   WRITING ARRHENIUS COEFFICIENTS
              ! ==============================================================================================================
              if (iBinsFn == 1) then
                format_char =  "'" // trim(adjustl(ReactionIn(iBinsIn))) // "=" // trim(adjustl(ReactionFn(iBinsFn))) // ":', SP,  ES11.4, ',', SP,  ES11.4, ',', SP,  ES11.4, ',2'"
              else
                format_char =  "'" // trim(adjustl(ReactionIn(iBinsIn))) // "=" // trim(adjustl(ReactionFn(iBinsFn))) // ":', SP,  ES11.4, ',', SP,  ES11.4, ',', SP,  ES11.4, ',6'"
              end if
              format_char = '(' // trim(adjustl(format_char)) // ')'
              
              if ( ( BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,1) < One ) .and. ( BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,1) > 1.e-100_rkp ) ) then
              
                if ( (TypeVec(iBinsFn) == 1) .and. (iBinsIn > iBinsFn - NArr_vec(iP+1)) .and. (trim(adjustl(Input%Kinetics_Exo)) == 'yes') ) then
                  write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,:)
                elseif ( (TypeVec(iBinsFn) == 1) .and. (iBinsIn < iBinsFn - NArr_vec(iP+1) ) .and. (trim(adjustl(Input%Kinetics_Endo)) == 'yes') ) then
                  write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,:)
                elseif ( (TypeVec(iBinsFn) == 0) .and. (trim(adjustl(Input%Kinetics_Diss)) == 'yes') ) then
                  write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,:)
                elseif ( (TypeVec(iBinsFn) == 2) .and. (trim(adjustl(Input%Kinetics_Exch)) == 'yes') ) then
                  write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBinsIn)%CArr(iBinsFn,:)
                end if
                
              end if
              ! ==============================================================================================================
                
                
            end do   ! iBinsFn
            
          end do   ! iBinsIn
          
        end if
        
      end if
      
    end do   ! iMol
    
    close(102)
        
        
    do iTtra = 1,Input%NTtra
        
      ! ==============================================================================================================
      !   CALLING KONIG
      ! ==============================================================================================================
      OutputPath = './T_' // trim(adjustl(Ttra_Vec_char(iTtra))) // '/output/'
      
      call system( 'cd ' // trim(adjustl(OutputPath)) // ' && rm -rf *.dat && ' // trim(adjustl(Input%KonigRunCMD)) // ' ' // &
                                                                   trim(adjustl(Input%KonigOrigDir)) // '/input/' // trim(adjustl(Input%KonigInputFileName)) // '_FWD_' // trim(adjustl(Ttra_Vec_char(iTtra))) // ' >/dev/null ' )
      ! ==============================================================================================================

  
      KonigOKFlag = .false.
      KonigOKFile = trim(adjustl(OutputPath ))// 'Tint.dat'
      inquire( file=trim(adjustl(KonigOKFile)), exist=KonigOKFlag)
      if (KonigOKFlag) then
        
        iWrite(iTtra) = iWrite(iTtra) + 1
    
      
        ! ==============================================================================================================
        !   READING KONIG OUTPUT
        ! ==============================================================================================================
        if (FirstTime) then
          call KonigPostprocessing%ReadComponents( OutputPath )
          FirstTime = .false.
        end if

        
        call KonigPostprocessing%ReadBox(  OutputPath, Input, iTtra )
        
        
        if (trim(adjustl(Input%FWD_Temperatures)) == 'yes') then
          call KonigPostprocessing%ReadTemperature( OutputPath, Input, iTtra )
        end if
        
        
        if (trim(adjustl(Input%FWD_Populations)) == 'yes') then
          call KonigPostprocessing%ReadPop(  OutputPath, Input, iTtra )
        end if
        ! ==============================================================================================================


        ! ==============================================================================================================
        !   UPLOADING KONIG OUTPUT DISTRIBUTION
        ! ==============================================================================================================
        if (trim(adjustl(Input%FWD_MolFractions)) == 'yes') then
          do iComp = 1,KonigPostprocessing%NComponents
            call KonigPostprocessing%LinearInterpAtAbscissaGrid( KonigPostprocessing%Ttra(iTtra)%Time, KonigPostprocessing%Ttra(iTtra)%X(:,iComp), KonigPostprocessing%Ttra(iTtra)%TimeGrid,  KonigPostprocessing%Ttra(iTtra)%XAtNodes(:,iComp,iWrite(iTtra)) )
          end do
        end if
        
        if (trim(adjustl(Input%FWD_Temperatures)) == 'yes') then
          do iBinnedMol=1,1
            call KonigPostprocessing%LinearInterpAtAbscissaGrid( KonigPostprocessing%Ttra(iTtra)%Time, KonigPostprocessing%Ttra(iTtra)%Tint(:,iBinnedMol), KonigPostprocessing%Ttra(iTtra)%TimeGrid,  KonigPostprocessing%Ttra(iTtra)%TIntAtNodes(:,iBinnedMol,iWrite(iTtra)) )
            ! KonigPostprocessing%Ttra(iTtra)%TIntSum(:)   =  KonigPostprocessing%Ttra(iTtra)%TIntSum(:)   +  KonigPostprocessing%Ttra(iTtra)%TIntAtNodes(:,iWrite(iTtra))  
            ! KonigPostprocessing%Ttra(iTtra)%TIntSqSum(:) =  KonigPostprocessing%Ttra(iTtra)%TIntSqSum(:) +  KonigPostprocessing%Ttra(iTtra)%TIntAtNodes(:,iWrite(iTtra))**2
          end do
        end if
        
        if (trim(adjustl(Input%FWD_Populations)) == 'yes') then
          do iBins = 1,Input%NBins(1)
            call KonigPostprocessing%LinearInterpAtAbscissaGrid( KonigPostprocessing%Ttra(iTtra)%Time, KonigPostprocessing%Ttra(iTtra)%Pop(:,iBins), KonigPostprocessing%Ttra(iTtra)%TimeGrid,  KonigPostprocessing%Ttra(iTtra)%PopAtNodes(:,iBins,iWrite(iTtra)) )
          end do
        end if
        ! ==============================================================================================================
        
        
        ! ==============================================================================================================
        !   DEALLOCATING KONIG OUTPUT
        ! ==============================================================================================================
        call KonigPostprocessing%DeallocateKonigPost( Input, iTtra )
        ! ==============================================================================================================


      end if

    end do   ! iTtra

  end do   ! iFwdProp
  
  
  do iTtra = 1,Input%NTtra
  
    Input%NFwdProp = iWrite(iTtra)
  
  
    ! ==============================================================================================================
    !   UPLOADING KONIG OUTPUT HISTOGRAM
    ! ==============================================================================================================  
    if (trim(adjustl(Input%FWD_MolFractions)) == 'yes') then
      do iComp = 1,KonigPostprocessing%NComponents
        do iTimeNodes = 1,Input%NTimeNodes
          call KonigPostprocessing%Histogram( KonigPostprocessing%Ttra(iTtra)%XBinsExtremes,  KonigPostprocessing%Ttra(iTtra)%XAtNodes(iTimeNodes,iComp,1:iWrite(iTtra)),  KonigPostprocessing%Ttra(iTtra)%XHist(:,iComp,iTimeNodes))
        end do
      end do
    end if
    
    if (trim(adjustl(Input%FWD_Temperatures)) == 'yes') then
      do iBinnedMol=1,1
        do iTimeNodes = 1,Input%NTimeNodes
          call KonigPostprocessing%Histogram( KonigPostprocessing%Ttra(iTtra)%TIntBinsExtremes,  KonigPostprocessing%Ttra(iTtra)%TIntAtNodes(iTimeNodes,iBinnedMol,1:iWrite(iTtra)),  KonigPostprocessing%Ttra(iTtra)%TIntHist(:,iBinnedMol,iTimeNodes))
        end do
      end do
    end if
    
    if (trim(adjustl(Input%FWD_Populations)) == 'yes') then
      do iBins = 1,Input%NBins(1)
        do iTimeNodes = 1,Input%NTimeNodes
          call KonigPostprocessing%Histogram( KonigPostprocessing%Ttra(iTtra)%PopBinsExtremes,  KonigPostprocessing%Ttra(iTtra)%PopAtNodes(iTimeNodes,iBins,1:iWrite(iTtra)),  KonigPostprocessing%Ttra(iTtra)%PopHist(:,iBins,iTimeNodes))
        end do
      end do
    end if
    ! ==============================================================================================================
    
    
    ! ==============================================================================================================
    !   WRITING KONIG OUTPUT DISTRIBUTION
    ! ==============================================================================================================
    OutputPath = './T_' // trim(adjustl(Ttra_Vec_char(iTtra))) // '/output/'
    
    FileName = 'TimeGrid.dat'
    call KonigPostprocessing%WriteTimeGrid( OutputPath, FileName, Input, KonigPostprocessing%Ttra(iTtra)%TimeGrid )
    
    if (trim(adjustl(Input%FWD_MolFractions)) == 'yes') then
      do iComp = 1,KonigPostprocessing%NComponents
        write(iComp_char,'(I2)') iComp
        ! FileName = 'TIntAtGrid.dat'
        ! call WriteOutputSums( OutputPath, FileName, Input, KonigPostprocessing%Ttra(iTtra)%TimeGrid,  KonigPostprocessing%Ttra(iTtra)%TIntSum,  KonigPostprocessing%Ttra(iTtra)%TIntSqSum, i_Debug )
        FileName = 'MolFracHist_Comp' //  trim(adjustl(iComp_char)) // '.dat'
        call KonigPostprocessing%WriteOutputHist( OutputPath, FileName, Input,  KonigPostprocessing%Ttra(iTtra)%XBinsExtremes,  KonigPostprocessing%Ttra(iTtra)%XHist(:,iComp,:) )
      end do
    end if
    
    if (trim(adjustl(Input%FWD_Temperatures)) == 'yes') then
      do iBinnedMol=1,1
        ! FileName = 'TIntAtGrid.dat'
        ! call WriteOutputSums( OutputPath, FileName, Input, KonigPostprocessing%Ttra(iTtra)%TimeGrid,  KonigPostprocessing%Ttra(iTtra)%TIntSum,  KonigPostprocessing%Ttra(iTtra)%TIntSqSum, i_Debug )
        FileName = trim(adjustl(Input%BinnedMolecules_Name(iBinnedMol))) // 'TIntHist.dat'
        call KonigPostprocessing%WriteOutputHist( OutputPath, FileName, Input,  KonigPostprocessing%Ttra(iTtra)%TIntBinsExtremes,  KonigPostprocessing%Ttra(iTtra)%TIntHist(:,iBinnedMol,:) )
      end do
    end if
    
    if (trim(adjustl(Input%FWD_Populations)) == 'yes') then
      do iBins = 1,Input%NBins(1)
        write(iBins_char,'(I6)') iBins
        ! FileName = 'TIntAtGrid.dat'
        ! call WriteOutputSums( OutputPath, FileName, Input, KonigPostprocessing%Ttra(iTtra)%TimeGrid,  KonigPostprocessing%Ttra(iTtra)%TIntSum,  KonigPostprocessing%Ttra(iTtra)%TIntSqSum, i_Debug )
        FileName = 'PopHist_Bin' // trim(adjustl(iBins_char)) // '.dat'
        call KonigPostprocessing%WriteOutputHist( OutputPath, FileName, Input,  KonigPostprocessing%Ttra(iTtra)%PopBinsExtremes,  KonigPostprocessing%Ttra(iTtra)%PopHist(:,iBins,:) )
      end do
    end if
    ! ==============================================================================================================

    
  end do
        
  
  call CPU_Time( EndTime )
  t_total   =   t_total + EndTime - StartTime
  if ( Logger%On() ) then
    call Logger%Write( "-> Total    : ", t_total, Fr="es15.8" )
  end if


  if ( Logger%On() ) call Logger%Write( "Normal termination" )

End Program FwdPropagation
