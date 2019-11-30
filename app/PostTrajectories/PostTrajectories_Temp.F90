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
!         Postprocessing Trajectories for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! This Program ...
!
!   Mandatory Input Arguments: - Translational Temperature [K]
!                              - Internal Temperature [K]
!                              - Nb of Trajectories, for checking convergence
!                              - Current Initial Levels/Bins
!
! ==============================================================================================================
Program PostTrajectories

  use Parameters_Module        ,only:  rkp, Zero, Half, One, Two, Six
  use Logger_Class             ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class              ,only:  Error
  use Input_Class              ,only:  Input_Type
  use LevelsContainer_Class    ,only:  LevelsContainer_Type, Compute_PartitionRatios
  use BinsContainer_Class      ,only:  BinsContainer_Type
  use Degeneracy_Class         ,only:  Degeneracy_Type
  use Degeneracy_Factory_Class ,only:  Degeneracy_Factory_Type
  use System_Class             ,only:  System_Type
  use System_Factory_Class     ,only:  System_Factory_Type
  use Molecule_Class           ,only:  Molecule_Type
  use Molecule_Factory_Class   ,only:  Molecule_Factory_Type
  use Collision_Class          ,only:  Collision_Type

  use Processes_Factory_Class  ,only:  Processes_Factory_Type
  use Processes_Class          ,only:  Processes_Type

  use Timing_Module

  implicit none
  
  Type(Input_Type)                                      :: Input

  Type(LevelsContainer_Type) ,dimension(:) ,allocatable :: LevelsContainer_Orig
  Type(LevelsContainer_Type) ,dimension(:) ,allocatable :: LevelsContainer_Cut

  Type(Degeneracy_Factory_Type)                         :: Degeneracy_Factory
  class(Degeneracy_Type)                   ,allocatable :: Degeneracy_Post
  type(Collision_Type)                                  :: Collision

  type(System_Factory_Type)                             :: System_Factory
  class(System_Type)                     ,allocatable   :: System_Adjst

  Type(BinsContainer_Type)   ,dimension(:) ,allocatable :: Molecule

  type(Molecule_Factory_Type)                           :: Molecule_Factory
  class(Molecule_Type)                   ,allocatable   :: Molecule_Adjst
 
  type(Processes_Factory_Type)                          :: Processes_Factory
  class(Processes_Type)                  ,allocatable   :: SetOfFinalProcesses

  real(rkp)                                             :: Tran
  integer                                               :: TranInt
  character(10)                                         :: Tran_char
  character(10)                                         :: TranInt_char
  character(30)                                         :: Velocity_char
  integer                                               :: iMol
  real(rkp)                                             :: Velocity
  integer                                               :: NTraj_temp
  character(10)                                         :: Tint_char
  character(:)                             ,allocatable :: FileName
  integer                                               :: Status
  integer                                               :: Unit
  real(rkp)                                             :: StartTime, EndTime

  logical                                  ,parameter   :: i_Debug_PT        = .True.
  logical                                  ,parameter   :: i_Debug_PT_Medium = .True.
  logical                                  ,parameter   :: i_Debug_PT_Deep   = .True.

  
  call CPU_Time( StartTime )
  
  if (i_Debug_PT) call Logger%Initialize(             "PostTrajectories.log",   &                                       ! Opening the Log File using
                              Status          =       'REPLACE',                &                                       ! replacing any previous log file
                              Position        =       'REWIND',                 &                                       ! rewinding to the top of the file
                              Procedure       =       'PostTrajectories',       &                                       ! loading the calling procedure name
                              Indentation     =       2           )                                                     ! and setting the initial indentation level


! ==============================================================================================================
!  READING THE PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 1, Tran_char )
  read(Tran_char, "(d20.10)", iostat=Status) Tran
  TranInt = int(Tran)
  write(TranInt_char, "(I10)", iostat=Status) TranInt
  if (Status/=0) call Error( "Error writing the argument TranInt_char" )
  if (i_Debug_PT) call Logger%Write( "TranInt_char = ", TranInt_char, "; Tran = ", Tran  )
  
  call getarg( 2, Input%Tint_char )
  read(Input%Tint_char, "(d20.10)", iostat=Status) Input%Tint
  if (Status/=0) call Error( "Error reading the argument Input%Tint" )
  write(Tint_char,"(I10)") int(Input%Tint)
  if (i_Debug_PT) call Logger%Write( "Input%Tint_char = ", Input%Tint_char, "; Input%Tint = ", Input%Tint, "Tint_char = ", Tint_char )
  
  call getarg( 3, Input%NTraj_char )
  read( Input%NTraj_char, '(I10)' ) NTraj_temp
  if (i_Debug_PT) call Logger%Write( "Input%NTraj_char = ", Input%NTraj_char, "; NTraj_temp = ", NTraj_temp )

  ! ==============================================================================================================
  !   READING VELOCITY CONVERSION FACTOR for CROSS SECTIONS RATES
  ! ==============================================================================================================
  call getarg( 4, Velocity_char )
  read( Velocity_char, "(d20.10)" ) Velocity
  if (i_Debug_PT) call Logger%Write( "Velocity_char = ", Velocity_char, "; Velocity for Conversion Factor = ", Velocity )
  ! ==============================================================================================================

! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Calling Input%Initialize" )
  call Input%Initialize( i_Debug=i_Debug_PT )

  Input%NTraj = NTraj_temp
  if (i_Debug_PT) call Logger%Write( "Total Nb of Converged Trajectories:     Input%NTraj = ", Input%NTraj )
  
  if ( trim(adjustl(Input%TtraModel)) .eq. "Boltzmann" ) then
    Input%Ttra      = Tran
    Input%Ttra_char = Tran_char
    if (i_Debug_PT) call Logger%Write( "Input%Ttra_char = ", Input%Ttra_char, "; Input%Ttra = ", Input%Ttra  )
  elseif ( trim(adjustl(Input%TtraModel)) .eq. "Uniform" ) then
    Input%Erel      = Tran
    Input%Erel_char = Tran_char
    if (i_Debug_PT) call Logger%Write( "Input%Erel_char = ", Input%Erel_char, "; Input%Erel = ", Input%Erel )
  end if
! ==============================================================================================================


! ==============================================================================================================
!  READING THE REMAINING PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 5, Input%PESoI_char )
  read( Input%PESoI_char,     '(I6)' ) Input%PESoI
  if (i_Debug_PT) call Logger%Write( "PES of Interest:     Input%PESoI = ", Input%PESoI)
  
  do iMol = 1,Input%NInitMolecules
  
    call getarg( 5 + iMol, Input%BinOI_char(iMol) )
    read( Input%BinOI_char(iMol),     '(I6)' ) Input%BinOI(iMol)
    if (i_Debug_PT) call Logger%Write( "Level/Bin of Interest for Molecule Nb", iMol, ":     Input%BinOI(iMol) = ", Input%BinOI(iMol))
    
  end do
! ==============================================================================================================
  
  
! ==============================================================================================================
!   INITIALIZING PAIR OBJECT
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Initializing the Collision%Pairs object", NewLine=.True. )
  if (i_Debug_PT) call Logger%Write( "Calling Collision%InitializePairs" )
  call Collision%InitializePairs( Input, i_Debug=i_Debug_PT )
  if (i_Debug_PT) call Logger%Write( "Done with Collision%InitializePairs" )
! ==============================================================================================================
  
  
! ==============================================================================================================
!   ASSIGNING MOLECULE to PAIRS
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Finding Corrispondences between Pairs and Molecules", NewLine=.True. )
  do iMol = 1,Input%NMolecules
    if (i_Debug_PT) call Logger%Write( "Molecule Nb", iMol )
    if (i_Debug_PT) call Logger%Write( "Calling Molecule_Factory%Define_Molecule" )
    call Molecule_Factory%Define_Molecule( Input, Molecule_Adjst, iMol, i_Debug=i_Debug_PT )
    if (i_Debug_PT) call Logger%Write( "Calling Molecule%AssignMoleculesToPairs" )
    call Molecule_Adjst%AssignMoleculesToPairs( Collision, Input, iMol, i_Debug=i_Debug_PT )
    if (i_Debug_PT) call Logger%Write( "Done with Molecule%AssignMoleculesToPairs" )
    deallocate(Molecule_Adjst)
    if (i_Debug_PT) call Logger%Write( "Molecule Deallocated" )
  end do
! ==============================================================================================================


! ==============================================================================================================
!    DEFINING ARRANGEMENTS
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Finding Equal Pairs in the System", NewLine=.True. )
  if (i_Debug_PT) call Logger%Write( "Calling System_Factory%Define_System" )
  call System_Factory%Define_System( Input, System_Adjst, i_Debug=i_Debug_PT )
  if (i_Debug_PT) call Logger%Write( "Calling System%AssignPairsArrangements" )
  call System_Adjst%AssignPairsArrangements( Collision, Input, i_Debug=i_Debug_PT )
  if (i_Debug_PT) call Logger%Write( "Done with System%AssignPairsArrangements" )
! ==============================================================================================================


! ==============================================================================================================
!   ALLOCATING LEVELS CONTAINERS
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Allocating the Molecules Original Energy Levels (containers for q.n.s etc) based on the Nb of Molecules" )
  allocate(LevelsContainer_Orig(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating LevelsContainer_Orig" )
  if (i_Debug_PT) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Orig" )

  if (i_Debug_PT) call Logger%Write( "Allocating the Molecules Cut Energy Levels (containers for q.n.s etc) based on the Nb of Molecules" )
  allocate(LevelsContainer_Cut(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating LevelsContainer_Cut" )
  if (i_Debug_PT) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Cut" )
  ! ==============================================================================================================
  
  
  ! ==============================================================================================================
  !   LOOKING FOR BINNED MOLECULES
  ! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Allocating the Molecules (containers for Levels/Bins) based on the Nb of Molecules" )
  
  allocate(Molecule(Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating Molecule" )
  if (i_Debug_PT) call Logger%Write( "Allocated ", Input%NMolecules, " Molecules" )
      
  if (i_Debug_PT) call Logger%Write( "Finding Corrispondences between Pairs and Molecules", NewLine=.True. )
  do iMol = 1,Input%NMolecules
    if (i_Debug_PT) call Logger%Write( "-> Molecule Nb", iMol )
    
    
    !==============================================================================================================
    ! INITIALIZING BINNED MOLECULE
    !==============================================================================================================
    if (i_Debug_PT) call Logger%Write( "Initialize Bins for the Molecule Nb", iMol )
    call Molecule(iMol)%InitializeBins( Input, iMol, i_Debug=i_Debug_PT  ) 
    ! ============================================================================================================== 
    
    
    !==============================================================================================================
    ! ASSIGNING BINNING MOLECULES TO PAIRS
    !==============================================================================================================
    if (i_Debug_PT) call Logger%Write( "Calling Molecule_Factory%Define_Molecule" )
    call Molecule_Factory%Define_Molecule( Input, Molecule_Adjst, iMol, i_Debug=i_Debug_PT )
    if (i_Debug_PT) call Logger%Write( "Calling Molecule%AssignBinnedMoleculesToPairs" )
    call Molecule_Adjst%AssignBinnedMoleculesToPairs( Collision, Input, iMol, i_Debug=i_Debug_PT )
    deallocate(Molecule_Adjst)
    if (i_Debug_PT) call Logger%Write( "Molecule Deallocated" )
    ! ============================================================================================================== 

  end do
  if (i_Debug_PT) call Logger%Write( "Done with Molecule%AssignBinnedMoleculesToPairs" )

  
  if (i_Debug_PT) call Logger%Write( "Iterating on all the Molecules in the Chemical System" )
  do iMol = 1,Input%NMolecules
    if (i_Debug_PT) call Logger%Write( "-> Molecule Nb", iMol )
    

    ! ==============================================================================================================
    !   READING MOLECULE ENERGY LEVELS
    ! ==============================================================================================================
    if (i_Debug_PT) call Logger%Write( "Reading the Original Levels List file for Molecule Nb", iMol )
    FileName = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // trim(adjustl(Input%LevelsFileName(iMol)))
    if (i_Debug_PT) call Logger%Write( "Reading File: ", FileName )
    call LevelsContainer_Orig(iMol)%Initialize( Input, iMol, FileName, i_Debug=i_Debug_PT )
    
    if (i_Debug_PT) call Logger%Write( "Reading the Cut Levels List file for Molecule Nb", iMol )
    FileName  = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/levels_cut.inp'
    if (i_Debug_PT) call Logger%Write( "Reading File: ", FileName )
    call LevelsContainer_Cut(iMol)%Initialize( Input, iMol, FileName, i_Debug=i_Debug_PT )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   COMPUTING PARTITION FUNCTIONS AND PARTITION FUNCTIONS RATIOS
    ! ==============================================================================================================
    if (i_Debug_PT) call Logger%Write( "Computing Partition Functions and related Ratios for Molecule Nb", iMol )
    call Compute_PartitionRatios( Input, iMol, LevelsContainer_Orig(iMol), LevelsContainer_Cut(iMol), i_Debug=i_Debug_PT )
    ! ==============================================================================================================
  
  
    ! if ( trim(adjustl(Input%TtraModel)) .eq. "Boltzmann" ) then
    !   ! ==============================================================================================================
    !   !   READING VELOCITY CONVERSION FACTOR for CROSS SECTIONS RATES
    !   ! ==============================================================================================================
    !   if (i_Debug_PT) call Logger%Write( "Reading the Velocity for Conversion Factor" )
    !   FileName = trim(adjustl(Input%OutputDir))// '/Velocity_' // trim(adjustl(Input%Ttra_char)) // '.dat'
    !   if (i_Debug_PT) call Logger%Write( "Reading File: ", FileName )
    !   open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
    !   if (Status/=0) call Error( "Postprocessing: Error opening file: " // FileName )
    !   read(Unit,*)
    !   read(Unit,*) gbx
    !   close(Unit)
    !   if (i_Debug_PT) call Logger%Write( "Velocity for Conversion Factor = ", gbx )
    !   ! ==============================================================================================================
    ! end if
   
    
    ! ==============================================================================================================
    !   SORTING LEVELS AND COMPUTING BINS PARTITION FUNCTIONS
    ! ==============================================================================================================
    if (trim(Input%BSortMethod(iMol)) .eq. "State-Specific" ) then 
      if (i_Debug_PT) call Logger%Write( "Sorting the Molecule Nb", iMol, " by State Specific" )
      call Molecule(iMol)%SortByStateSpecific( Input, LevelsContainer_Cut(iMol), iMol, .false., i_Debug=i_Debug_PT_Deep )
    elseif (trim(Input%BSortMethod(iMol)) .eq. "From-File" ) then 
      if (i_Debug_PT) call Logger%Write( "Sorting the Molecule Nb", iMol, " from File's Mapping'" )
      call Molecule(iMol)%SortFromFile( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_PT_Deep )
    elseif (trim(Input%BSortMethod(iMol)) .eq. "Vib-Specific" ) then 
      if (i_Debug_PT) call Logger%Write( "Sorting the Molecule Nb", iMol, " by Vibrational Specific" )
      call Molecule(iMol)%SortByVibSpecific( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_PT_Deep )
    elseif (trim(Input%BSortMethod(iMol)) .eq. "RoVib-CG" ) then 
      if (i_Debug_PT) call Logger%Write( "Sorting the Molecule Nb", iMol, " by increasing Ro-Vibrational Energy" )
      call Molecule(iMol)%SortByRoVibEnergy( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_PT_Deep )
    else
      if (i_Debug_PT) call Logger%Write( "ERROR: Sorting Method Unknown! Input%BSortMethod(i)", trim(Input%BSortMethod(iMol))) !!! @TODO: Adding New Sorting Methods
      stop ( "ERROR: Sorting Method Unknown!" )    
    end if
    
    if (i_Debug_PT) call Logger%Write( "Reading Bins Energies and Partition Functions" )
    call Molecule(iMol)%ReadPartFunEnergy( Input, LevelsContainer_Cut(iMol), iMol, i_Debug=i_Debug_PT_Deep )
    ! ==============================================================================================================


  end do
  

  call Processes_Factory%Define_Processes( Input, Collision, LevelsContainer_Orig, SetOfFinalProcesses, i_Debug=i_Debug_PT )

  call SetOfFinalProcesses%Convert_CrossSect_To_Rates( Input, Collision, LevelsContainer_Orig, Molecule, [Velocity], i_Debug=i_Debug_PT, i_Debug_Deep=i_Debug_PT_Medium )

  !
  call CPU_Time( EndTime )
  t_total   =   t_total + EndTime - StartTime
  if (i_Debug_PT) then
    call Logger%Write( "Total    : ", t_total, Fr="es15.8" )
  end if
  

  if (i_Debug_PT) call Logger%Write( "Normal termination" )

End Program PostTrajectories
