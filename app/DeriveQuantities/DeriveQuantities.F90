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
!  Computing Derived Quantities from Rates for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! This Program ...
!
!   Mandatory Input Arguments: - Nb of Translational Temperatures
!                              - Nb of Internal Temperatures
!                              - Nb of Trajectories
!                              - For iMol = 1,NMol
!                                 - BinStart(iMol)
!                                 - FinalStart(iMol)
!
! ==============================================================================================================
Program DeriveQuantities

  use Parameters_Module        ,only:  rkp, Zero, Half, One, Two, Six, Ten
  use Logger_Class             ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class              ,only:  Error
  use Input_Class              ,only:  Input_Type
  use LevelsContainer_Class    ,only:  LevelsContainer_Type
  use BinsContainer_Class      ,only:  BinsContainer_Type
  use Collision_Class          ,only:  Collision_Type
  use System_Class             ,only:  System_Type
  use System_Factory_Class     ,only:  System_Factory_Type
  use Molecule_Class           ,only:  Molecule_Type
  use Molecule_Factory_Class   ,only:  Molecule_Factory_Type
  use Rates_Class 
  use Arrhenius_Class
  use KONIG_Class 
  use HEGEL_Class              
  use Timing_Module

  implicit none
  
  Type(Input_Type)                                      :: Input
  Type(LevelsContainer_Type) ,dimension(:) ,allocatable :: LevelsContainer_Orig
  Type(LevelsContainer_Type),dimension(:,:),allocatable :: LevelsContainerT
  Type(BinsContainer_Type)   ,dimension(:) ,allocatable :: BinnedMolecule
  type(Collision_Type)                                  :: Collision
  type(System_Factory_Type)                             :: System_Factory
  class(System_Type)                       ,allocatable :: System_Adjst
  type(Molecule_Factory_Type)                           :: Molecule_Factory
  class(Molecule_Type)                     ,allocatable :: Molecule_Adjst
  type(Rates_Type)                                      :: Rates
  type(Arrhenius_Type)                                  :: Arrhenius
  
  integer                                               :: iMol, jMol
  integer                                               :: iBins
  character(6)                                          :: iBinsChar
  integer                                               :: NTraj
  integer                                               :: NTraj_temp
  integer   ,dimension(:), allocatable                  :: NTraj_Tot
  real(rkp) ,dimension(:), allocatable                  :: TtraVec
  character(10)                                         :: TintChar
  integer                                               :: iTint
  character(10)                                         :: TtraChar
  integer                                               :: iTtra
  character(2)                                          :: iTtraChar
  character(:)                             ,allocatable :: FileName
  integer                                               :: Status
  real(rkp)                                             :: StartTime, EndTime
  
  logical                                               :: i_Debug_DQ        = .True.
  logical                                               :: i_Debug_DQ_Medium = .True.
  logical                                               :: i_Debug_DQ_Deep   = .True.
  
  
  call CPU_Time( StartTime )
  
  if (i_Debug_DQ) call Logger%Initialize(             "DeriveQuantities.log",  &                                        ! Opening the Log File using
                              Status          =       'REPLACE',               &                                        ! replacing any previous log file
                              Position        =       'REWIND',                &                                        ! rewinding to the top of the file
                              Procedure       =       'DeriveQuantities',      &                                        ! loading the calling procedure name
                              Indentation     =       2           )                                                     ! and setting the initial indentation level

! ==============================================================================================================
!   INITIALIZING INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_DQ) call Logger%Write( "Calling Input%Initialize" )
  call Input%Initialize( i_Debug=i_Debug_DQ )
  if (i_Debug_DQ) call Logger%Write( "Calling Input%DeriveQuantities" )
  call Input%DeriveQuantities( i_Debug=i_Debug_DQ )
  if (i_Debug_DQ) call Logger%Write( "Done reading Input for DeriveQuantities" )
  Input%TaskType = 7
  
  allocate(NTraj_Tot(Input%NTtra), Stat=Status )
  if (Status/=0) call Error( "Error allocating NTraj_Tot" )
  if (i_Debug_DQ) call Logger%Write( "Allocated NTraj_Tot.")
  NTraj_Tot = 0
  
  allocate(TtraVec(max(Input%NTtra,3)), Stat=Status )
  if (Status/=0) call Error( "Error allocating TtraVec" )
  if (i_Debug_DQ) call Logger%Write( "Allocated TtraVec.")
  TtraVec = Zero
! ==============================================================================================================


! ==============================================================================================================
!  READING THE PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 1, Input%NTraj_char )
  read( Input%NTraj_char, '(I10)' ) NTraj_temp
  Input%NTraj = NTraj_temp
  if (i_Debug_DQ) call Logger%Write( "Total Nb of Converged Trajectories:     Input%NTraj = ", Input%NTraj )
  
  call getarg( 2, Input%iPES_char )
  read( Input%iPES_char, '(I6)', iostat=Status ) Input%iPES
  if (Status/=0) call Error( "Error reading the argument Input%iPES" )
  if (i_Debug_DQ) call Logger%Write( "Current PES:    Input%iPES = ", Input%iPES )
  
  do iMol = 1,Input%NInitMolecules
  
    call getarg( 3 + (iMol-1)*3, Input%NBins_char(iMol) )
    read( Input%NBins_char(iMol), '(I6)', iostat=Status) Input%NBins(iMol)
    if (Status/=0) call Error( "Error reading the argument Input%NBins(iMol)" )
    if (i_Debug_DQ) call Logger%Write( "Nb of Levels/Bins for Molecule Nb", iMol, ":     Input%NBins(iMol) = ", Input%NBins(iMol) )
  
    call getarg( 4 + (iMol-1)*3, Input%BinStart_char(iMol) )
    read( Input%BinStart_char(iMol), '(I6)', iostat=Status) Input%BinStart(iMol)
    if (Status/=0) call Error( "Error reading the argument Input%BinStart(iMol)" )
    if (i_Debug_DQ) call Logger%Write( "Starting Level/Bin for Molecule Nb", iMol, ":     Input%BinStart(iMol) = ", Input%BinStart(iMol) )
    
    call getarg( 5 + (iMol-1)*3, Input%BinFinal_char(iMol) )
    read( Input%BinFinal_char(iMol), '(I6)', iostat=Status ) Input%BinFinal(iMol) 
    if (Status/=0) call Error( "Error reading the argument Input%BinFinal(iMol)" )
    if (i_Debug_DQ) call Logger%Write( "Final Level/Bin for Molecule Nb", iMol, ":     Input%BinFinal(iMol) = ", Input%BinFinal(iMol) )
  
  end do
! ==============================================================================================================
  
  
! ==============================================================================================================
!   INITIALIZING PAIR OBJECT
! ==============================================================================================================
  if (i_Debug_DQ) call Logger%Write( "Initializing the Collision%Pairs object. Calling Collision%InitializePairs", NewLine=.True. )
  call Collision%InitializePairs( Input, i_Debug=i_Debug_DQ )
  if (i_Debug_DQ) call Logger%Write( "Done with Collision%InitializePairs" )
! ==============================================================================================================

  
! ==============================================================================================================
!   ASSIGNING MOLECULE to PAIRS
! ==============================================================================================================
  if (i_Debug_DQ) call Logger%Write( "Finding Corrispondences between Pairs and Molecules", NewLine=.True. )
  do iMol = 1,Input%NMolecules
    if (i_Debug_DQ) call Logger%Write( "Molecule Nb", iMol )
    if (i_Debug_DQ) call Logger%Write( "Calling Molecule_Factory%Define_Molecule" )
    call Molecule_Factory%Define_Molecule( Input, Molecule_Adjst, iMol, i_Debug=i_Debug_DQ )
    if (i_Debug_DQ) call Logger%Write( "Calling Molecule%AssignMoleculesToPairs" )
    call Molecule_Adjst%AssignMoleculesToPairs( Collision, Input, iMol, i_Debug=i_Debug_DQ )
    if (i_Debug_DQ) call Logger%Write( "Done with Molecule%AssignMoleculesToPairs" )
    deallocate(Molecule_Adjst)
    if (i_Debug_DQ) call Logger%Write( "Molecule Deallocated" )
  end do
! ==============================================================================================================


! ==============================================================================================================
!    DEFINING ARRANGEMENTS
! ==============================================================================================================
  if (i_Debug_DQ) call Logger%Write( "Finding Equal Pairs in the System", NewLine=.True. )
  if (i_Debug_DQ) call Logger%Write( "Calling System_Factory%Define_System" )
  call System_Factory%Define_System( Input, System_Adjst, i_Debug=i_Debug_DQ )
  if (i_Debug_DQ) call Logger%Write( "Calling System%AssignPairsArrangements" )
  call System_Adjst%AssignPairsArrangements( Collision, Input, i_Debug=i_Debug_DQ )
  if (i_Debug_DQ) call Logger%Write( "Done with System%AssignPairsArrangements" )
! ==============================================================================================================


! ==============================================================================================================
!   ALLOCATING LEVELS CONTAINERS
! ==============================================================================================================
  if (i_Debug_DQ) call Logger%Write( "Allocating the Molecules Original Energy Levels (containers for q.n.s etc) based on the Nb of Molecules" )
  allocate(LevelsContainer_Orig(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating LevelsContainer_Orig" )
  if (i_Debug_DQ) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Orig" )

  if (i_Debug_DQ) call Logger%Write( "Allocating the Molecules Cut Energy Levels (containers for q.n.s etc) based on the Nb of Molecules" )
  allocate(LevelsContainerT(Input%NMolecules,Input%NTtra), Stat=Status )
  if (Status/=0) call Error( "Error allocating LevelsContainerT" )
  if (i_Debug_DQ) call Logger%Write( "Allocated LevelsContainerT with dimension ( ", Input%NMolecules, ",", Input%NTtra, " )" )
! ==============================================================================================================


! ==============================================================================================================
!   LOOKING FOR MOLECULES
! ==============================================================================================================
  if (i_Debug_DQ) call Logger%Write( "Allocating the Binned-Molecules (containers for bins etc.) based on the Nb of Molecules" )
  
  allocate(BinnedMolecule(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating BinnedMolecule" )
  if (i_Debug_DQ) call Logger%Write( "Allocated ", Input%NMolecules, " BinnedMolecule" )
  
  
  if (i_Debug_DQ) call Logger%Write( "Finding Corrispondences between Pairs and Molecules", NewLine=.True. )
  do iMol = 1,Input%NMolecules
    if (i_Debug_DQ) call Logger%Write( "Binned Molecule Nb", iMol )
    
    
    !==============================================================================================================
    ! INITIALIZING BINNED MOLECULE
    !==============================================================================================================
    if (i_Debug_DQ) call Logger%Write( "Initialize Bins for the Molecule Nb", iMol )
    call BinnedMolecule(iMol)%InitializeBins( Input, iMol, i_Debug=i_Debug_DQ  ) 
    ! ============================================================================================================== 
    
    
    !==============================================================================================================
    ! ASSIGNING BINNING MOLECULES TO PAIRS
    !==============================================================================================================
    if (i_Debug_DQ) call Logger%Write( "Calling Molecule_Factory%Define_Molecule" )
    call Molecule_Factory%Define_Molecule( Input, Molecule_Adjst, iMol, i_Debug=i_Debug_DQ )
    if (i_Debug_DQ) call Logger%Write( "Calling Molecule%AssignBinnedMoleculesToPairs" )
    call Molecule_Adjst%AssignBinnedMoleculesToPairs( Collision, Input, iMol, i_Debug=i_Debug_DQ )
    deallocate(Molecule_Adjst)
    if (i_Debug_DQ) call Logger%Write( "Molecule Deallocated" )
    ! ============================================================================================================== 
    
  end do
  if (i_Debug_DQ) call Logger%Write( "Done with Molecule%AssignBinnedMoleculesToPairs" )
! ==============================================================================================================
  

  if (i_Debug_DQ) call Logger%Write( "Iterating on all the Molecules in the Chemical System" )
  do iMol = 1,Input%NMolecules
    if (i_Debug_DQ) call Logger%Write( "Molecule Nb", iMol )
    
    
    do iTtra = 1,Input%NTtra
    
      ! ==============================================================================================================
      !   READING MOLECULE ENERGY LEVELS
      ! ==============================================================================================================
      if (i_Debug_DQ) call Logger%Write( "Reading the Original Levels List file for Molecule Nb", iMol )
      FileName = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // trim(adjustl(Input%LevelsFileName(iMol)))
      if (i_Debug_DQ) call Logger%Write( "Reading File: ", FileName )
      call LevelsContainerT(iMol,iTtra)%Initialize( Input, iMol, FileName, i_Debug=i_Debug_DQ )

      if (i_Debug_DQ) call Logger%Write( "Reading the Cut Levels List file for Molecule Nb", iMol )
      FileName  = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/levels_cut.inp'
      if (i_Debug_DQ) call Logger%Write( "Reading File: ", FileName )
      call LevelsContainerT(iMol,iTtra)%Initialize( Input, iMol, FileName, i_Debug=i_Debug_DQ )
      ! ==============================================================================================================
      
      
      ! ==============================================================================================================
      !   READING THE MAPPING STS CG
      ! ==============================================================================================================
      if (i_Debug_DQ) call Logger%Write( "Reading File: ", FileName )
      if (Input%BinStsFlg) call BinnedMolecule(iMol)%ReadLevelToBin( Input, LevelsContainerT(iMol,iTtra), iMol, i_Debug=i_Debug_DQ )
      ! ==============================================================================================================
      
    end do


  end do


  !==============================================================================================================
  ! COUNTING THE RATES "CONTAINERS" REQUIRED
  !==============================================================================================================
  if ((Input%ReadFormRatesFlg)) call Rates%DifferentiatingFinalStates( Input, Collision, i_Debug=i_Debug_DQ )  
  ! ===============================================================================================================
  
  
  ! ==============================================================================================================
  !   CREATING FORMAT and HEADER FOR OVERALL RATES OUTPUT FILE
  ! ==============================================================================================================
  if (Input%WriteOverallRatesFlg) call Rates%CreatingDescriptorsOverall( Input, Collision, i_Debug=i_Debug_DQ )
  ! ==============================================================================================================

  
  ! ==============================================================================================================
  !   ALLOCATING COARSE-GRAINED RATES 
  ! ==============================================================================================================
  if (Input%BinStsFlg) call Rates%AllocatingCGRates( Input, Collision, LevelsContainerT, i_Debug=i_Debug_DQ )
  ! ==============================================================================================================
  

  if (Input%RunKonigFlg) call KONIGInquiringVariables( Input, i_Debug=i_Debug_DQ )  
  if (Input%RunHEGELFlg) call HEGELInquiringVariables( Input, i_Debug=i_Debug_DQ )  

  
  do iMol = 1,Input%NMolecules
    if (i_Debug_DQ) call Logger%Write( "Molecule Nb", iMol )
    
    if ((iMol == Collision%Pairs(1)%To_Molecule)) then

      
      ! ==============================================================================================================
      !   CREATING HEADER FOR RATES OUTPUT FILE
      ! ==============================================================================================================
      if (Input%ReadFormRatesFlg) call Rates%CreatingDescriptors( Input, Collision, i_Debug=i_Debug_DQ, i_Debug_Deep=i_Debug_DQ_Deep )
      ! ==============================================================================================================


      ! ==============================================================================================================
      !   CREATING FORMAT for ARRHENIUS COEFFICIENTS OUTPUT FILE
      ! ==============================================================================================================
      if ((Input%ReadFormRatesFlg) .and. ((Input%WriteArrFlg) .or. (Input%WriteKonigFlg) .or. (Input%WriteHegelFlg))) call Arrhenius%CreatingFormat( Input, Collision, Rates, iMol, i_Debug=i_Debug_DQ )
      ! ==============================================================================================================


      ! ==============================================================================================================
      !   ALLOCATING OVERALL RATES 
      ! ==============================================================================================================
      if (Input%WriteOverallRatesFlg) call Rates%AllocatingOverall( Input, BinnedMolecule(iMol), i_Debug=i_Debug_DQ )
      ! ==============================================================================================================


      ! ==============================================================================================================
      !   OPENING KINETIC FILE for KONIG / HEGEL
      ! ==============================================================================================================
      if (Input%RunKonigFlg) call KONIGOpeningKinetics( Input, i_Debug=i_Debug_DQ )
      if (Input%RunHEGELFlg) call HEGELOpeningKinetics( Input, i_Debug=i_Debug_DQ )
      ! ==============================================================================================================


      do iBins = Input%BinStart(iMol),Input%BinFinal(iMol)
        write(iBinsChar,'(I0)') iBins
      

        ! ==============================================================================================================
        !   ALLOCATING RATES
        ! ==============================================================================================================
        if ((Input%ReadFormRatesFlg)) call Rates%Allocating( Input, Collision, BinnedMolecule(iMol), iBins, i_Debug=i_Debug_DQ_Deep )
        ! ==============================================================================================================


        do iTint = 1,Input%NTint      
                                                                                                     
          
          ! ==============================================================================================================
          !   WRITING ARRHENIUS COEFFICIENTS
          ! ==============================================================================================================
          if (Input%WriteArrFlg) call ArrOpeningFiles( Input, iMol, iBins, iTint, i_Debug=i_Debug_DQ )
          ! ==============================================================================================================
          
          
          ! ==============================================================================================================
          !   ALLOCATING ARRHENIUS COEFFICIENTS
          ! ==============================================================================================================
          if ((Input%WriteArrFlg) .or. (Input%WriteKonigFlg) .or. (Input%WriteHegelFlg)) call ArrAllocating( Input, Collision, BinnedMolecule(iMol), Rates, iBins, i_Debug=i_Debug_DQ_Deep )
          ! ==============================================================================================================
          

          do iTtra = 1,Input%NTtra
            
            if (Input%NTint==1) then
              Input%Tint = Input%TtraVec(iTtra)
            else
              Input%Tint = Input%TintVec(iTint)
            end if
            Input%Ttra = Input%TtraVec(iTtra)
            
            write(TtraChar,'(I0)') int(Input%Ttra)
            if (i_Debug_DQ) call Logger%Write( "Input%Ttra = ", Input%Ttra, "; TtraChar = ", TtraChar )            
            write(TintChar,'(I0)') int(Input%Tint)
            if (i_Debug_DQ) call Logger%Write( "Input%Tint = ", Input%Tint, "; TintChar = ", TintChar )
            
            
            if (iBins == Input%BinStart(iMol)) then
              
              do jMol = 1,Input%NMolecules
              
                ! ==============================================================================================================
                !   SORTING LEVELS AND COMPUTING BINS PARTITION FUNCTIONS
                ! ==============================================================================================================
                if (trim(Input%BSortMethod(jMol)) .eq. "State-Specific" ) then 
                  if (i_Debug_DQ) call Logger%Write( "Sorting the Molecule Nb", jMol, " by State Specific" )
                  call BinnedMolecule(jMol)%SortByStateSpecific( Input, LevelsContainerT(jMol,iTtra), jMol, .false., i_Debug=i_Debug_DQ_Deep )
                 elseif (trim(Input%BSortMethod(iMol)) .eq. "From-File" ) then 
                  if (i_Debug_DQ) call Logger%Write( "Sorting the Molecule Nb", iMol, " from File's Mapping'" )
                  call BinnedMolecule(jMol)%SortFromFile( Input, LevelsContainerT(jMol,iTtra), jMol, i_Debug=i_Debug_DQ_Deep )
                elseif (trim(Input%BSortMethod(jMol)) .eq. "Vib-Specific" ) then 
                  if (i_Debug_DQ) call Logger%Write( "Sorting the Molecule Nb", jMol, " by Vibrational Specific" )
                  call BinnedMolecule(jMol)%SortByVibSpecific( Input, LevelsContainerT(jMol,iTtra), jMol, i_Debug=i_Debug_DQ_Deep )
                elseif (trim(Input%BSortMethod(jMol)) .eq. "RoVib-CG" ) then 
                  if (i_Debug_DQ) call Logger%Write( "Sorting the Molecule Nb", jMol, " by increasing Ro-Vibrational Energy" )
                  call BinnedMolecule(jMol)%SortByRoVibEnergy( Input, LevelsContainerT(jMol,iTtra), jMol, i_Debug=i_Debug_DQ_Deep )
                elseif (trim(Input%BSortMethod(jMol)) .eq. "Hybrid" ) then
                  if (i_Debug_DQ) call Logger%Write( "Sorting the Molecule Nb", jMol, " by Hybrid Mthd" )
                  call BinnedMolecule(jMol)%SortByHybrid( Input, LevelsContainerT(jMol,iTtra), jMol, i_Debug=i_Debug_DQ_Deep )
                else
                  if (i_Debug_DQ) call Logger%Write( "ERROR: Sorting Method Unknown! Input%BSortMethod(i)", trim(Input%BSortMethod(jMol))) !!! @TODO: Adding New Sorting Methods
                  stop ( "ERROR: Sorting Method Unknown!" )    
                end if
              
                
                if (i_Debug_DQ) call Logger%Write( "Reading Bins Energies and Partition Functions" )
                call BinnedMolecule(jMol)%ReadPartFunEnergy( Input, LevelsContainerT(jMol,iTtra), jMol, i_Debug=i_Debug_DQ_Deep )
                ! ==============================================================================================================
                
              end do
              

              ! ==============================================================================================================
              if ((Input%WriteFormRatesFlg)) then
              
                call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/Rates/')
                call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/Rates/T_' // trim(adjustl(TtraChar)) // '_' // trim(adjustl(TintChar)) )             
                 
              end if
              ! ==============================================================================================================
              
              
              ! ==============================================================================================================
              !   ALLOCATING COARSE-GRAINED RATES 
              ! ==============================================================================================================
              if (Input%BinStsFlg) call Rates%ComputingQ( Input, LevelsContainerT(iMol,iTtra), Input%Tint, iTtra, i_Debug=i_Debug_DQ_Deep )
              ! ==============================================================================================================

            end if
            ! ==============================================================================================================
          
      
            ! ==============================================================================================================
            !   READING LEVELS / BINS FORMATTED RATES
            ! ==============================================================================================================
            if (Input%ReadFormRatesFlg) call Rates%Reading( Input, BinnedMolecule(iMol), iMol, LevelsContainerT(iMol,iTtra), iBins, iBinsChar, TintChar, iTtra, TtraChar, NTraj, i_Debug=i_Debug_DQ_Deep )
            NTraj_Tot(iTtra) = NTraj_Tot(iTtra) + NTraj
            ! ==============================================================================================================


            ! ==============================================================================================================
            !   WRITING BINS RATES
            ! ==============================================================================================================
            if (Input%WriteFormRatesFlg) call Rates%Writing( Input, Collision, BinnedMolecule(iMol), iMol, iBins, iBinsChar, TtraChar, TintChar, NTraj, i_Debug=i_Debug_DQ_Deep )
            ! ==============================================================================================================
            
            
            ! ==============================================================================================================
            !   CREATING VECTORS OF TEMPERATURES and RATES
            ! ==============================================================================================================
            if ((Input%WriteArrFlg) .or. (Input%WriteKonigFlg) .or. (Input%WriteHegelFlg)) call  ArrCreatingVectors( Input, BinnedMolecule(iMol), iBins, TtraVec, iTtra, i_Debug=i_Debug_DQ_Deep )
            ! ==============================================================================================================

            
            !! ==============================================================================================================
            !!   MERGING RATES FOR INTERNAL EXCHANGE WITH THE ONES FOR INELASTIC PROCESSES
            !! ==============================================================================================================
            !if ((.not. Input%ConsiderExcFlg) .and. ((Input%WriteArrFlg) .or. (Input%RunKonigFlg) .or. (Input%WriteAllTsFlg) .or. (Input%WriteUnformRatesFlg))) &
            !                                                                          call ArrMergeIntExch( Input, Collision, BinnedMolecule(iMol), Rates, iBins, iTtra, i_Debug=i_Debug_DQ )
            !! ==============================================================================================================


            ! ==============================================================================================================
            !   COMPUTE THERMAL RATES
            ! ==============================================================================================================
            if (Input%WriteOverallRatesFlg) call Rates%ComputingOverall( Input, Collision, BinnedMolecule(iMol), iBins, iTtra, i_Debug=i_Debug_DQ_Deep )
            ! ==============================================================================================================
          

          end do   ! iTtra

          ! ==============================================================================================================
          !   COMPUTING and WRITING ARRHENIUS COEFFICIENTS
          ! ==============================================================================================================
          if ((Input%WriteArrFlg) .or. (Input%WriteKonigFlg) .or. (Input%WriteHegelFlg)) call Arrhenius%ComputingWriting( Input, Collision, BinnedMolecule(:), iMol, LevelsContainerT(:,1), Rates, iBins, iBinsChar, TtraVec, i_Debug=i_Debug_DQ_Medium ) 
          ! ==============================================================================================================


          ! ==============================================================================================================
          !   DeALLOCATING RATES
          ! ==============================================================================================================
          if ((Input%ReadFormRatesFlg)) call Rates%DeAllocating( Input, BinnedMolecule(iMol), iBins, i_Debug=i_Debug_DQ_Medium )
          ! ==============================================================================================================


          ! ==============================================================================================================
          !   DeALLOCATING ARRHENIUS COEFFICIENTS
          ! ==============================================================================================================
          if ((Input%WriteArrFlg) .or. (Input%WriteKonigFlg) .or. (Input%WriteHegelFlg)) call ArrDeAllocating( Input, BinnedMolecule(iMol), iBins, i_Debug=i_Debug_DQ_Medium )
          ! ==============================================================================================================

          
        end do   ! iTint
                        
      end do   ! iBins 

      
    ! ==============================================================================================================
    !   CREATING OUTPUT FOLDER FOR KONIG
    ! ==============================================================================================================
    if (Input%RunKonigFlg) call KONIGCreatingDirs( Input, iMol, i_Debug=i_Debug_DQ )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   WRITING OVERALL RATE
    ! ==============================================================================================================
    if (Input%WriteOverallRatesFlg) call Rates%WritingOverall( Input, BinnedMolecule, iMol, Input%TintVec, TintChar, Input%TtraVec, NTraj_Tot, i_Debug=i_Debug_DQ )
    ! ==============================================================================================================
    
    
    ! ==============================================================================================================
    !   COMPUTE CG RATES
    ! ==============================================================================================================
    if (Input%BinStsFlg) call Rates%WritingCGRates( Input, Collision, Input%TintVec, TintChar, Input%TtraVec, i_Debug=i_Debug_DQ )
    ! ==============================================================================================================

  
    end if


    ! ==============================================================================================================
    !   WRITE TABLES for HEGEL
    ! ==============================================================================================================
    if (Input%RunHegelFlg) call HEGELComputeTables( Input, LevelsContainerT(iMol,1), BinnedMolecule(iMol), Collision%Species, i_Debug=i_Debug_DQ )
    ! ==============================================================================================================
    
  end do   ! iMol
  
  
  ! ==============================================================================================================
  !   FINISH TO WRITE ARRHENIUS FOR KONIG / HEGEL
  ! ==============================================================================================================
  if (Input%WriteKonigFlg) call KONIGFinishingWritingKinetics( Input, i_Debug=i_Debug_DQ )
  if (Input%WriteHegelFlg) call HEGELFinishingWritingKinetics( Input, i_Debug=i_Debug_DQ )
  ! ==============================================================================================================  
  
  
  if (Input%RunKonigFlg) then
    if (i_Debug_DQ) call Logger%Write( "Input%TtraVecChar = ", Input%TtraVecChar )
    
    ! ==============================================================================================================
    !   COPYING THERMO FILE
    ! ==============================================================================================================
    call KONIGWritingThermo( Input, BinnedMolecule, i_Debug=i_Debug_DQ )
    ! ============================================================================================================== 
    
    ! ==============================================================================================================
    !   COMPUTING AND WRITING INITIAL POPULATIONS
    ! ==============================================================================================================
    call KONIGComputingWritingPop0( Input, BinnedMolecule, i_Debug=i_Debug_DQ )
    ! ==============================================================================================================   
      
    ! ==============================================================================================================
    !   CALLING EXTERNAL CODE
    ! ==============================================================================================================
    call KONIGRunningExtCode( Input, i_Debug=i_Debug_DQ )
    ! ==============================================================================================================   
    
  end if
  
  
  call CPU_Time( EndTime )
  t_total   =   t_total + EndTime - StartTime
  if (i_Debug_DQ) call Logger%Write( "Total    : ", t_total, Fr="es15.8" )

  if (i_Debug_DQ) call Logger%Write( "Normal termination" )

End Program DeriveQuantities
