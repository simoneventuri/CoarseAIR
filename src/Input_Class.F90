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

Module Input_Class

  use Parameters_Module      ,only:  rkp, Zero, One, autime_to_sec
  use Logger_Class           ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class            ,only:  Error
  
  implicit none

  private
  public    ::    Input_Type

  Type      ::    Input_Type

    character(150)                            ::    InputDir                            !< Input Directory
    character(150)                            ::    InputFile                           !< Input File for CoarseAIR
    character(150)                            ::    OutputDir                           !< Output Directory
    character(150)                            ::    LevelOutputDir                        !< Bin Output Directory
    
    character(3)                              ::    GenLev
    character(3)                              ::    RunPrep
    character(3)                              ::    RunTraj
    character(3)                              ::    RunPost
    character(3)                              ::    CompQuant

    integer                                   ::    NNode                               !< Total Nb of Nodes
    character(5)                              ::    NNode_char
    integer                                   ::    iNode                               !< Node Identifier
    character(5)                              ::    iNode_char
    integer                                   ::    iProc                               !< Processor Identifier
    character(5)                              ::    iProc_char
    integer                                   ::    NProc            = 1                !< Total Nb of Processors
    character(5)                              ::    NProc_char

    integer                                   ::    TaskType         = 0                !< = 1, Plotting PESs
                                                                                        !< = 2, Computing and Merging Levels
                                                                                        !< = 3, Preprocessing Levels
                                                                                        !< = 4, Running Trajectories
                                                                                        !< = 5, Computing Statistics
                                                                                        !< = 6, Postprocessing Trajectories
                                                                                        !< = 7, Derivanging Quantities
    integer                                   ::    NTint
    integer                                   ::    NTtra
    character(2)                              ::    NTtra_char
    real(rkp)      ,dimension(:) ,allocatable ::    TtraVec
    character(10)  ,dimension(:) ,allocatable ::    TtraVecChar
    integer                                   ::    NErel
    character(2)                              ::    NErel_char
    real(rkp)      ,dimension(:) ,allocatable ::    ErelVec
    character(10)  ,dimension(:) ,allocatable ::    ErelVecChar
    character(10)  ,dimension(:) ,allocatable ::    ErelVecIntChar
    character(2)                              ::    NTint_char
    character(10)  ,dimension(:) ,allocatable ::    TtraVecIntChar
    real(rkp)      ,dimension(:) ,allocatable ::    TintVec

    integer        ,dimension(:) ,allocatable ::    NBins                               !< Vector with the Total Nbs of Bins for each of the Binned Molecules
    integer        ,dimension(:) ,allocatable ::    BinOI                               !< Vector with the Current Nbs of Bins for each of the Binned Molecules
    character(6)   ,dimension(:) ,allocatable ::    NBins_char
    character(6)   ,dimension(:) ,allocatable ::    BinOI_char
    integer        ,dimension(:) ,allocatable ::    BinStart
    character(6)   ,dimension(:) ,allocatable ::    BinStart_char
    integer        ,dimension(:) ,allocatable ::    BinFinal
    character(6)   ,dimension(:) ,allocatable ::    BinFinal_char
    

! GENERATE LEVELS
! =================
    integer        ,dimension(:) ,allocatable ::    ComputeLevels
    integer                                   ::    iPair
    real(rkp)      ,dimension(2)              ::    xExtremes
    integer                                   ::    NGrid
    real                                      ::    DeltaGrid        = -1.e0_rkp
    real(rkp)                                 ::    EStart
    real(rkp)                                 ::    DeltaE
    real(rkp)                                 ::    EEpsilon
    integer                                   ::    ivqnMin           = 0
    integer                                   ::    ijqnMin           = 0
    integer                                   ::    ivqnMax           = 100
    integer                                   ::    ijqnMax           = 300
    real(rkp)      ,dimension(2,3)            ::    BoundayConditions
    character(3)                              ::    Levels_for_UQ
    character(150) ,dimension(:) ,allocatable ::    GeneratedLevelsFile
    logical                                   ::    PrintLevelsFlg = .false.
    logical                                   ::    SortLevelsFlg  = .true.

! CUT INPUTS
! =================
    character(150)                            ::    DtbPath                             !< Database Folder Path
    character(:)               ,allocatable   ::    System                              !< Chemical System of Interest
    integer                                   ::    NMolecules                          !< Nb of Molecules in the System
    character(5) ,dimension(:) ,allocatable   ::    Molecules_Name                      !< Molecules Names
    integer      ,dimension(:) ,allocatable   ::    Molecule_NAtoms                     !< Nbs of Atoms in each of the Molecules
    character(150),dimension(:) ,allocatable  ::    LevelsFileName                      !< Names of the Molecules Energy Levels File
    real(rkp)    ,dimension(:) ,allocatable   ::    Ecut                                !< Minum Energy required for considering Level [Eh]
    real(rkp)    ,dimension(:) ,allocatable   ::    Tcut                                !< Minum Time required for considering Level [ps]

! SORT INPUTS
! =================
    character(50) ,dimension(:)  ,allocatable ::    BSortMethod                         !< Methods for sorting the Levels of each molecule
    integer       ,dimension(:,:),allocatable ::    BSortInfo                           !< Info about such sorting (Ex: Starting/Final Level, Nb of Bins, etc.)
    character(150),dimension(:,:),allocatable ::    BinDataFile
    
! DRIVER INPUTS
! =================
    logical                                   ::    ParamsFlg   =   .false.
    logical                                   ::    ProgressFlg =   .false.
    logical                                   ::    PaQSolFlg   =   .false.
    logical                                   ::    PESEvoFlg   =   .false.
    logical                                   ::    PaQEvoFlg   =   .false.
    logical                                   ::    XXEvoFlg    =   .false.
    logical                                   ::    XXEvoSnglFileFlg = .false.
    logical                                   ::    EnEvoFlg    =   .false.
    real(rkp)                                 ::    TimeMax     =   Zero                !< CPU time limit.
    integer                                   ::    NAtoms      =   0                   !< Nb of atoms
    character(5)  ,dimension(:) ,allocatable  ::    AtomsName                           !< Names of each atoms
    real(rkp)     ,dimension(:) ,allocatable  ::    AtomsMass                           !< Masses of each atoms [a.u.]
    character(150),dimension(:) ,allocatable  ::    FileForAtomMass                     !< Where to Read atom masses
    character(:)               ,allocatable   ::    ProportionalAllocation              !< Flag for using Proportional Allocation for Stratified Sampling when Sampling into the Bins
    integer                                   ::    NTraj       =   0                   !< Nb of trajectories
    character(10)                             ::    NTraj_char  =   ''
    integer                                   ::    NTrajBatch  =   0
    character(10)                             ::    NTrajNTrajBatch_char =   ''
    integer                                   ::    NBatch      =   0 
    logical                                   ::    Restart     =   .False.             !< Restart indicator
    integer                                   ::    NSpecies    =   0                   !< Nb of species
    character(6)  ,dimension(:) ,allocatable  ::    SpeciesName                         !< Name of the species involved in the collision. For now, dim=2
    integer                                   ::    iseed       =   0                   !< Initial random Nb seed
    integer                                   ::    PESiseed    =   0                   !< Random Nb seed for PES parameters sampling      
    character(3)                              ::    Randomize   =   'no'                !< Generate a CPU-clock-dependent seed?
    real(rkp)                                 ::    RndNb       =   0.d0                !< CPU-clock-dependent seed
    real(rkp)                                 ::    dinit       =   Zero                !< Initial distance between target and projectile
    logical                                   ::    BwrdIntegrationFlg = .false.        !< Flag for Backward Integration
    integer                                   ::    BwrdIntegrationFreq = 0             !< Backward Integration Frequency
!   ODE parmeters
    real(rkp)                                 ::    dt          =   Zero                !< initial time step for integrator
    integer                                   ::    ncall       =   0                   !< Nb of calls to integrator between end checking
    real(rkp)                                 ::    eps         =   Zero                !< error control parameter
    integer                                   ::    NSteps      =   0                   !< Nb of extrapolation steps to use
    character(3)                              ::    IncStpSzChar=   'yes'
    logical                                   ::    IncStpSzFlg =   .True.
    real(rkp)                                 ::    Relax       =   Zero                !< stepsize increase damping factor.
    real(rkp)                                 ::    rmax        =   Zero                !< if any atom pair distance exceeds rmax, the trajectory is assumed finished.
    real(rkp)                                 ::    tmax        =   Zero                !< if any trajectory lives longer the tmax, it will be stoped.
!   Relative kinetic energy and impact parameter variables
    character(:)               ,allocatable   ::    TtraModel                           !< Model for Translational Mode
    real(rkp)                                 ::    Ttra        =   Zero                !< Translational Temperature [K] 
    character(10)                             ::    Ttra_char
    real(rkp)                                 ::    Erel        =   Zero                !< Translational Energy [Eh]
    character(10)                             ::    Erel_char
    character(:)               ,allocatable   ::    TintModel                           !< Model for Internal Modes
    real(rkp)                                 ::    Tint        =   Zero                !< Internal temperature [K]
    character(10)                             ::    Tint_char
    character(:)               ,allocatable   ::    ImpParModel                         !< Model for Sampling the Impact Parameter
    character(:)               ,allocatable   ::    ImpParStrataType
    real(rkp)                                 ::    Etot        =   Zero                !< Total energy for trajectories
    real(rkp)                                 ::    TotalAngularMomentum = Zero         !< Imposed total angular momentum. Read from input only if 'bstart=2345432'. Old name: xjtot
    integer                                   ::    nring       =   0                   !< Nb of stratified sampling used for the impact parameters
    integer      ,dimension(:) ,allocatable   ::    ntring                              !< Nb of trajectories in each range of impact parameter. Dim=(nring). Old: DIM=maxr=10
    real(rkp)    ,dimension(:) ,allocatable   ::    bring                               !< Impact parameter value in each range of impact parameter. Dim=(nring). Old: DIM=maxtj
    integer                                   ::    NPESs                               !< Nb of Potential Energy Surfaces (PESs) / Nb of Samples from the Stochastci PES
    integer                                   ::    iPES        =   0                   !< Current PES
    character(6)                              ::    iPES_char
    integer                                   ::    PESoI                               !< Potential Energy Surfaces (PES) of Interest
    character(6)                              ::    PESoI_char
    logical                                   ::    PESRndIniCondFlg = .False.
    logical                                   ::    StochPESFlg      = .False.
    logical                                   ::    SampleParamsStochPES = .False.
    integer      ,dimension(2)                ::    NHiddenLayerNeurons
    character(150),dimension(:),allocatable   ::    PES_Model                           !< Potential Energy Surface's Model
    character(150)                            ::    PES_ParamsFldr
    character(150)                            ::    PES_ParamsFile = 'NONE'
    real(rkp)     ,dimension(:),allocatable   ::    PES_Degeneracy                      !< Potential Energy Surface's Degeneracy
    real(rkp)     ,dimension(:),allocatable   ::    PES_Degeneracy_Vec                  !< Potential Energy Surface's Degeneracy Vector
    character(150),dimension(:),allocatable   ::    Diatomic_Model                      !< Model for Diatomic Potential
    character(150),dimension(:),allocatable   ::    DiatPot_ParamsFile
    real(rkp)                                 ::    eintShift   =   Zero                !< Shift for Internal Energy Levels
    integer                                   ::    NInitMolecules

! Inputs for LevelsContainer_Class
! =================================
    integer                                   ::    jIn         =   0
    integer                                   ::    vIn         =   0
    integer                                   ::    iodd
    character(:)               ,allocatable   ::    ioddChar
    character(:)               ,allocatable   ::    jInMethod
    character(:)               ,allocatable   ::    vInMethod
    integer                                   ::    UnitArbDist =  31

! ANALYSIS INPUTS
! ===============
    integer                                   ::    nquad       =  50                   !<
    real(rkp)                                 ::    rsmal                               !<
    real(rkp)                                 ::    rcent                               !<
    real(rkp)                                 ::    rfact                               !<
    real(rkp)                                 ::    tthres                              !<
    real(rkp)                                 ::    tthresau                            !<

! STATISTICS INPUTS
! =================
    logical                                   ::    StatReadsBinaryFlg    = .False.
    logical                                   ::    StatWritesBinaryFlg   = .False.
    integer                                   ::    NCond
    integer      ,dimension(:) ,allocatable   ::    iskip
    logical      ,dimension(:) ,allocatable   ::    PresEvOdd
    character(:)                ,allocatable  ::    AssignmentMethod                    !< Method for states assigment
    integer                                   ::    NTrajectoriesToAnalyze = -1         !< Max Nb of trajectories to analyze. -1 means analyze all trajectories on file.
    logical                                   ::    IdenticalDiatoms       = .false.

! POSTPROCESSING TRJS INPUTS
! =================
    character(:)               ,allocatable   ::    ConsiderExc
    logical                                   ::    ConsiderExcFlg
    logical                                   ::    MergeExchsFlg         = .false.
    logical                                   ::    MergeExchToInelFlg    = .false.
    character(3)                              ::    ConsiderQB  = 'yes'

! MERGE PESs
! =================
    integer                                   ::    NPESsToMerge
    real(rkp)      ,dimension(:) ,allocatable ::    PESWeight
    character(150) ,dimension(:) ,allocatable ::    ResultsFile

! RUNNING EXTERNAL CODE
! =================
    logical                                   ::    MergeIntExchFlg        = .false.
    character(3)                              ::    MergeIntExchChar
    logical                                   ::    WriteOverallRatesFlg   = .false.
    character(3)                              ::    WriteOverallRatesChar
    logical                                   ::    BinStsFlg              = .false.
    character(3)                              ::    ReadAllTsChar
    logical                                   ::    ReadAllTsFlg           = .false.
    character(3)                              ::    WriteAllTsChar
    logical                                   ::    WriteAllTsFlg          = .false.
    logical                                   ::    ReadFormRatesFlg       = .false.
    character(3)                              ::    ReadFormRatesChar
    logical                                   ::    ReadUnformRatesFlg     = .false.
    character(3)                              ::    ReadUnformRatesChar    
    character(150)                            ::    RatesFolderPath        = 'default'
    logical                                   ::    WriteFormRatesFlg      = .false.
    character(3)                              ::    WriteFormRatesChar
    logical                                   ::    WriteUnformRatesFlg    = .false.
    character(3)                              ::    WriteUnformRatesChar
    logical                                   ::    WriteArrFlg            = .false.
    character(3)                              ::    WriteArrChar
    character(3)                              ::    BinStsChar
    integer      ,dimension(:) ,allocatable   ::    NBinsCG
    character(6) ,dimension(:) ,allocatable   ::    NBinsCG_char
    character(150)                            ::    SystemCGPath
    real(rkp)                                 ::    MaxForRates
    real(rkp)                                 ::    MinForRates            = 1.d-30
    real(rkp)                                 ::    TInit                  = 300.d0              !< Initial Temperature
    character(10)                             ::    TInit_char
    character(150)                            ::    Kinetics_ExoChar       
    logical                                   ::    Kinetics_ExoFlg        = .true.
    character(150)                            ::    Kinetics_EndoChar      
    logical                                   ::    Kinetics_EndoFlg       = .false.
    character(150)                            ::    Kinetics_InelChar      
    logical                                   ::    Kinetics_InelFlg       = .false.
    character(150)                            ::    Kinetics_DissChar      
    logical                                   ::    Kinetics_DissFlg       = .false.
    character(150)                            ::    Kinetics_ExchChar      
    logical                                   ::    Kinetics_ExchFlg       = .false.
    character(150)                            ::    Kinetics_COCDissChar   
    logical                                   ::    Kinetics_COCDissFlg    = .false.
    character(150)                            ::    Kinetics_O2DissChar    
    logical                                   ::    Kinetics_O2DissFlg     = .false.
    character(5) ,dimension(:) ,allocatable   ::    CompName
    real(rkp)    ,dimension(:) ,allocatable   ::    CompMolFrac
    character(6)                              ::    NProcExtCodeChar       = '    16'
    integer                                   ::    NProcForExtCode        = 16
    real(rkp)                                 ::    DissCorrFactor         = 1.d0
    integer                                   ::    RunExtCodeIntFlg

! RUNNING KONIG
! =================
    logical                                   ::    WriteKONIGFlg          = .false.
    character(3)                              ::    WriteKONIGChar
    logical                                   ::    RunKONIGFlg            = .false.
    character(3)                              ::    RunKONIGChar
    character(200)                            ::    KONIGRunCMD
    character(200)                            ::    KONIGDtbPath
    character(200)                            ::    KONIGOrigDir
    character(100)                            ::    KONIGInputFileName = 'input'

! RUNNING HEGEL
! =================
    logical                                   ::    WriteHEGELFlg          = .false.
    character(3)                              ::    WriteHEGELChar
    logical                                   ::    RunHEGELFlg            = .false.
    character(3)                              ::    RunHEGELChar
    character(200)                            ::    HEGELRunCMD
    character(200)                            ::    HEGELDtbPath
    character(200)                            ::    HEGELOrigDir
    character(100)                            ::    HEGELInputFileName = 'input'

! FWD PROPAGATION
! =================
    integer                                   ::    NFwdProp
    integer                                   ::    NFwdPropProc
    integer                                   ::    FwdSeed
    real(rkp)                                 ::    FwdTimeMin
    real(rkp)                                 ::    FwdTimeMax
    integer                                   ::    NTimeNodes
    character(3)                              ::    TimeScale
    character(:)               ,allocatable   ::    FWD_MolFractions
    integer                                   ::    NComp
    real(rkp)                                 ::    MolFractionsMin
    real(rkp)                                 ::    MolFractionsMax
    integer                                   ::    NMolFractionsBins
    character(:)               ,allocatable   ::    FWD_Temperatures
    real(rkp)                                 ::    TemperaturesMin
    real(rkp)                                 ::    TemperaturesMax
    integer                                   ::    NTemperaturesBins
    character(:)               ,allocatable   ::    FWD_Populations
    real(rkp)                                 ::    PopulationsMin
    real(rkp)                                 ::    PopulationsMax
    integer                                   ::    NPopulationsBins
    
! PLOTTING PES
! ==================    
    logical                                   ::    PlotPES_ReadPntsFlg        = .False.
    logical                                   ::    PlotPES_GridFlg            = .False.
    logical                                   ::    PlotPES_DoubleGridFlg      = .False.
    logical                                   ::    PlotPES_TripleGridFlg      = .False.
    logical                                   ::    PlotPES_GridForScatterFlg  = .False.
    logical                                   ::    PlotPES_StatsFlg           = .False.
    logical                                   ::    PlotPES_WriteParamsFlg     = .False.
    logical                                   ::    PlotPES_OnlyTriatFlg       = .False.
    logical                                   ::    PlotPES_VargasPaperFlg     = .False.
    logical                                   ::    PlotPES_Rot3rdFlg          = .False.
    logical                                   ::    PlotPES_IsoTriFlg          = .False.
    integer       ,dimension(2)               ::    NGridPlot
    real(rkp)     ,dimension(2)               ::    MinGridPlot
    real(rkp)     ,dimension(2)               ::    MaxGridPlot
    integer                                   ::    NAnglesPlot
    real(rkp)     ,dimension(:) ,allocatable  ::    AnglesPlot
    character(20) ,dimension(:) ,allocatable  ::    AnglesPlotChar
    character(20)                             ::    POTorFR
    character(20)                             ::    PESOrDiatFlg
    character(20)                             ::    YAxisVar
    character(20)                             ::    UnitFR
    character(20)                             ::    UnitDist
    character(20)                             ::    UnitPot
    character(3)                              ::    RandomPointsChar = ' no'
    logical                                   ::    RandomPointsFlg  = .false.
    character(20)                             ::    PointsDistr
    character(20) ,dimension(5)               ::    DistrParChar
    real(rkp)     ,dimension(5)               ::    DistrPar
    real(rkp)                                 ::    EnergyCutOff     = 1.d10
    logical                                   ::    PESZeroRefFlg    = .False.
    
  contains
    private
    procedure     ,public                     ::    InitializeIO        =>    InitializeIOpaths
    procedure     ,public                     ::    Initialize          =>    InitializeCGQCTInput
    procedure     ,public                     ::    WriteBashInputVariables
    procedure     ,public                     ::    GenerateLevels
    procedure     ,public                     ::    DeriveQuantities
    procedure     ,public                     ::    FwdProagation
    procedure     ,public                     ::    PlotPES
    procedure                                 ::    ReadAtomMass
  End Type

  integer ,parameter                        ::    DefaultLogLevel = LogLevel_INFO
  INTEGER                                   ::    i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE        ::    seed
  real(rkp)                                 ::    RndNb
  logical                                   ::    i_Debug_Global = .True.!.False.

  contains

Subroutine InitializeIOpaths( This, i_Debug )

  class(Input_Type)                           ,intent(inout)  ::    This
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    len
  integer                                                     ::    status
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeIOpaths" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Write( "Inquiring for $COARSEAIR_INPUT_DIR Environment Variable" )
  call get_environment_variable ("COARSEAIR_INPUT_DIR", This%InputDir, len, status, .true.)
  if (status .ge. 2) then
    call Logger%Write( "get_environment_variable failed for $InputDir: status = ", status )
    stop
  elseif (status .eq. 1) then
    call Logger%Write( "$COARSEAIR_INPUT_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_INPUT_DIR: length = ", len, " truncated to 150")
    len = 150
  end if
  if (len .eq. 0) then
    call Logger%Write( "$COARSEAIR_INPUT_DIR exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "Input Directory Path, This%InputDir = ", This%InputDir)
  end if

  if (i_Debug_Loc) call Logger%Write( "Inquiring for $COARSEAIR_INPUT_FILE Environment Variable" )
  call get_environment_variable ("COARSEAIR_INPUT_FILE", This%InputFile, len, status, .true.)
  if (status .ge. 2) then
    call Logger%Write( "get_environment_variable failed for $InputFile: status = ", status )
    stop
  elseif (status .eq. 1) then
    call Logger%Write( "$COARSEAIR_INPUT_FILE does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_INPUT_FILE: length = ", len, " truncated to 150")
    len = 150
  end if
  if (len .eq. 0) then
    call Logger%Write( "$COARSEAIR_INPUT_FILE exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "Input File Path, This%InputFile = ", This%InputFile)
  end if

  if (i_Debug_Loc) call Logger%Write( "Inquiring for $COARSEAIR_OUTPUT_DIR Environment Variable" )
  call get_environment_variable ("COARSEAIR_OUTPUT_DIR", This%OutputDir, len, status, .true.)
  if (status .ge. 2) then
    call Logger%Write( "get_environment_variable failed for $OutputDir: status = ", status )
    stop
  elseif (status .eq. 1) then
    call Logger%Write( "$COARSEAIR_OUTPUT_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_OUTPUT_DIR: length = ", len, " truncated to 150")
    len = 150
  end if
  if (len .eq. 0) then
    call Logger%Write( "$COARSEAIR_OUTPUT_DIR exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "Output Directory Path, This%OutputDir = ", This%OutputDir)
  end if

  if (i_Debug_Loc) call Logger%Write( "Inquiring for $COARSEAIR_BIN_OUTPUT_DIR Environment Variable" )
  call get_environment_variable ("COARSEAIR_BIN_OUTPUT_DIR", This%LevelOutputDir, len, status, .true.)
  if (status .ge. 2) then
    call Logger%Write( "get_environment_variable failed for $LevelOutputDir: status = ", status )
    stop
  elseif (status .eq. 1) then
    call Logger%Write( "$COARSEAIR_BIN_OUTPUT_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_BIN_OUTPUT_DIR: length = ", len, " truncated to 150")
    len = 150
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_BIN_OUTPUT_DIR exists, but has no value")
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "Bin Output Directory Path, This%OutputDir = ", This%LevelOutputDir)
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine InitializeCGQCTInput( This, i_Debug)

  class(Input_Type)                         ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, Unit
  integer                                                   ::    Status
  character(:)  ,allocatable                                ::    FileName
  integer                                                   ::    NLevels
  character(100)                                            ::    Line

  character(150)                                            ::    i_case
  character(150)                                            ::    line_input
  character(150)                                            ::    line_input_temp
  integer                                                   ::    i_eq

  integer                                                   ::    iSp
  character(2)                                              ::    iSp_char
  character(:) ,allocatable                                 ::    Sp_case

  integer                                                   ::    iMol
  character(2)                                              ::    iMol_char
  character(:) ,allocatable                                 ::    Mol_case
  character(:) ,allocatable                                 ::    Molecule_NAtoms_case
  character(:) ,allocatable                                 ::    LevelsFileName_case
  character(:) ,allocatable                                 ::    Ecut_case
  character(:) ,allocatable                                 ::    Tcut_case
  character(:) ,allocatable                                 ::    DiatPot_case
  
  integer                                                   ::    iAt
  character(2)                                              ::    iAt_char
  character(:) ,allocatable                                 ::    At_case
  character(:) ,allocatable                                 ::    AtMassFile_case
  character(:) ,allocatable                                 ::    AtMass_case1, AtMass_case2

  integer                                                   ::    iring
  character(2)                                              ::    iring_char
  character(:) ,allocatable                                 ::    ntring_case
  character(:) ,allocatable                                 ::    bring_case

  integer                                                   ::    iCond
  character(2)                                              ::    iCond_char
  character(150)                                            ::    iskip_case

  integer                                                   ::    iRotQuaNum
  character(2)                                              ::    iRotQuaNum_char
  character(150)                                            ::    PresEvOdd_case

  character(:) ,allocatable                                 ::    BSortMethod_case
  character(:) ,allocatable                                 ::    BinDataFile_case
  character(:) ,allocatable                                 ::    MinLevel_case, MaxLevel_case
  character(:) ,allocatable                                 ::    Minvqn_case, Maxvqn_case

  integer                                                   ::    iiPES, IHL
  character(2)                                              ::    iiPES_char, iHL_char
  character(:) ,allocatable                                 ::    PESModel_case, PESDegeneracy_case, NHL_case           
                       
  character(20)                                             ::    iTint_char, iTtra_char, iErel_char
  integer                                                   ::    iTint, iTtra, iErel
  character(:) ,allocatable                                 ::    Tint_case, Tint_case_bis, Ttra_case, Ttra_case_bis, Erel_case
  character(:) ,allocatable                                 ::    DiatPot_File_case

  character(150)                                            ::    NBins_case
  
  integer      ,dimension(64)                               ::    SeedVec
  
  character(150)                                            ::    PES_Model_Temp 

  logical                                                   ::    ExFlg
  character(150)                                            ::    MolFldr              

  integer                                                   ::    len
  logical                                                   ::    i_Debug_Loc = .true.


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeCGQCTInput" )
  !i_Debug_Loc   =     Logger%On()


  This%NTtra            = 1
  This%NErel            = 1
  This%NTint            = 1
  This%NSpecies         = 0
  This%NMolecules       = 0
  This%NInitMolecules   = 1
  This%NAtoms           = 0
  This%nring            = 0
  This%NCond            = 0
  This%NPESs            = 1
  allocate( This%PES_Model(This%NPESs), stat=Status )
  if (Status/=0) call Error( "Error allocating This%PES_Model" )
  allocate( This%PES_Degeneracy(This%NPESs), stat=Status )
  if (Status/=0) call Error( "Error allocating This%PES_Degeneracy" )
  This%PES_Degeneracy = One
  allocate( This%PES_Degeneracy_Vec(This%NPESs), stat=Status )
  if (Status/=0) call Error( "Error allocating This%PES_Degeneracy_Vec" )
  This%PES_Degeneracy_Vec = One
  
  This%ImpParStrataType = 'Circles'
  
  This%MaxForRates      = Huge(Zero)

!   --------------------------------------------------------------------------------------------------------------
!   --------------------------------------------------------------------------------------------------------------
  call InitializeIOpaths( This, i_Debug_Loc )
  FileName    = trim(This%InputFile)
  if (i_Debug_Loc) call Logger%Write( "Reading the Input Parameter file" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )

  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )

!   --------------------------------------------------------------------------------------------------------------
  do
    read(Unit,'(A150)',iostat=Status) line_input
    if (Status /= 0) then
      exit
    else
      line_input_temp=trim(adjustl(line_input))
      if ( (line_input_temp(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
        i_eq = 1
        do
          if (line_input(i_eq:i_eq) == '=') exit
          i_eq = i_eq + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_eq-2))))
        !write(*,*) (adjustl((TRIM(i_case))))

        select case (adjustl((TRIM(i_case))))

          case("System")
            This%System = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "System:      This%System  = ", This%System )
          
          case("Nb of Species")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NSpecies
            if (i_Debug_Loc) call Logger%Write( "Nb of Species:      This%NSpecies  = ", This%NSpecies )
            allocate( This%SpeciesName(This%NSpecies), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%SpeciesName" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%SpeciesNames to This%NSpecies = ", This%NSpecies )
            
          case("Nb of Molecules")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NMolecules
            if (i_Debug_Loc) call Logger%Write( "Nb of Molecules:      This%NMolecules  = ", This%NMolecules )
            allocate( This%Molecules_Name(This%NMolecules), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%Molecules_Name" )
            allocate( This%LevelsFileName(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%LevelsFileName" )
            allocate( This%Ecut(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%Ecut" )
            allocate( This%Tcut(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%Tcut" )
            allocate( This%Molecule_NAtoms(This%NMolecules), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%Molecule_NAtoms" )
            allocate( This%BSortMethod(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%SortMethod" )
            allocate( This%BinDataFile(This%NMolecules,2), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinDataFile" )
            This%BinDataFile(:,1) = 'LevToBin.dat'
            This%BinDataFile(:,1) = '            '
            allocate( This%NBins(This%NMolecules), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%NBins" )
            allocate( This%NBins_char(This%NMolecules), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%NBins_char" )
            This%NBins = 0
            allocate(This%BinOI_char(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinOI" )
            allocate(This%BinOI(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinOI_Char" )
            allocate(This%BinStart(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinStart" )
            allocate(This%BinStart_char(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinStart_char" )
            allocate(This%BinFinal(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinFinal" )
            allocate(This%BinFinal_char(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BinFinal_char" )        
            allocate(This%BSortInfo(This%NMolecules,3), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%BSortInfo" )
            This%BSortInfo      = 0
            This%BSortInfo(:,1) = 1
            allocate(This%Diatomic_Model(This%NMolecules), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%Diatomic_Model" ) 
            allocate( This%DiatPot_ParamsFile(This%NMolecules), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%DiatPot_ParamsFile" )
            This%DiatPot_ParamsFile = 'NONE'


        case("Nb of Initial Molecules")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NInitMolecules
            if (i_Debug_Loc) call Logger%Write( "Nb of Initial Molecules:      This%NInitMolecules  = ", This%NInitMolecules )
            
          case("Nb of Atoms")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NAtoms
            if (i_Debug_Loc) call Logger%Write( "Nb of Atoms:      This%TimeMax  = ", This%NAtoms )
            allocate( This%AtomsName(This%NAtoms), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%AtomsName" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%AtomsName to This%NAtoms = ", This%NAtoms )
            allocate( This%AtomsMass(This%NAtoms), Stat=Status )
            This%AtomsMass = Zero
            if (Status/=0) call Error( "Error allocating This%AtomsMass" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%AtomsMass to This%NAtoms = ", This%NAtoms )
            allocate( This%FileForAtomMass(This%NAtoms), Stat=Status )
            This%FileForAtomMass = 'NONE'
            if (Status/=0) call Error( "Error allocating This%FileForAtomMass" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%FileForAtomMass to This%NAtoms = ", This%NAtoms )

            
            
          case("Database Path")
            This%DtbPath = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Database Path:      This%DtbPath  = ", This%DtbPath )
            if (trim(This%DtbPath) == 'FromEnvVariable') then
              if (i_Debug_Loc) call Logger%Write( "Inquiring for $COARSEAIR_DTB_DIR Environment Variable" )
              call get_environment_variable ("COARSEAIR_DTB_DIR", This%DtbPath, len, status, .true.)
              This%Levels_for_UQ = 'yes'
              if (status .ge. 2) then
                if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $DtbPath: status = ", status )
                stop
              elseif (status .eq. 1) then
                if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_DTB_DIR does not exist" )
                stop
              elseif (status .eq. -1) then
                if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_DTB_DIR: length = ", len, " truncated to 150")
                len = 150
              end if
              if (len .eq. 0) then
                if (i_Debug_Loc) call Logger%Write( "$COARSEAIR_DTB_DIR exists, but has no value")
                stop
              end if
              if (status .eq. 0) then
                if (i_Debug_Loc) call Logger%Write( "Input Directory Path, This%DtbPath = ", This%DtbPath)
              end if
            end if
            
            
          case("Nb of Translational Energies")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NErel
            if (i_Debug_Loc) call Logger%Write( "Nb of Translational Energies = ",  This%NErel )
            allocate( This%ErelVec(This%NErel), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%ErelVec" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%ErelVec to This%NErel = ", This%NErel )
            This%ErelVec = Zero
            allocate( This%ErelVecChar(This%NErel), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%ErelVecChar" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%ErelVecChar to This%NErel = ", This%NErel )
            allocate( This%ErelVecIntChar(This%NErel), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%ErelVecIntChar" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%ErelVecIntChar to This%NErel = ", This%NErel )
            
            
          case("Nb of Translational Temperatures")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NTtra
            if (i_Debug_Loc) call Logger%Write( "Nb of Translational Temperatures = ",  This%NTtra )
            allocate( This%TtraVec(This%NTtra), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%TtraVec" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%TtraVec to This%NTtra = ", This%NTtra )
            This%TtraVec = Zero
            allocate( This%TtraVecChar(This%NTtra), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%TtraVecChar" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%TtraVecChar to This%NTtra = ", This%NTtra )
            allocate( This%TtraVecIntChar(This%NTtra), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%TtraVecIntChar" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%TtraVecIntChar to This%NTtra = ", This%NTtra )
            
          case("Nb of Internal Temperatures")
            This%NTint_char = trim(adjustl(line_input(i_eq+2:150)))
            READ(This%NTint_char, '(I2)') This%NTint
            if (i_Debug_Loc) call Logger%Write( "Nb of Internal Temperatures = ",  This%NTint )
            allocate( This%TintVec(This%NTint), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%TintVec" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%TintVec to This%NTint = ", This%NTint )
            This%TintVec = Zero
      
      
          case("Generating Energy Levels?")
            This%GenLev = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Generating Energy Levels?:      This%GenLev  = ", This%GenLev )
            
          case("Preprocessing Energy Levels List?")
            This%RunPrep = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Preprocessing Energy Levels List?:      This%System  = ", This%RunPrep )
            
          case("Running Trajectories?")
            This%RunTraj = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Running Trajectories?:      This%System  = ", This%RunTraj )
            
          case("Postprocessing Trajectories?")
            This%RunPost = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Postprocessing Trajectories?:      This%System  = ", This%RunPost )
            
          case("Generating Derivated Quantities?")
            This%CompQuant = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Generating Derivated Quantities?:      This%System  = ", This%CompQuant )
            
            
!          case("Internal Energy Levels Shift")
!            line_input = line_input(i_eq+2:150)
!            READ(line_input, '(d20.10)') This%eintShift
            
            
          case("Using Proportional Allocation?")
            This%ProportionalAllocation = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Using Proportional Allocation?:      This%ProportionalAllocation  = ", trim(This%ProportionalAllocation) )
            
          case("Nb of Total Trajectories")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%NTraj
            if (i_Debug_Loc) call Logger%Write( "Nb of Total Trajectories:      This%NTraj  = ", This%NTraj )
            
          case("Nb of Trajectories per Level / Bin")
              line_input = line_input(i_eq+2:150)
              READ(line_input, '(I20)') This%NTraj
              if (i_Debug_Loc) call Logger%Write( "Nb of Trajectories per Level / Bin:      This%NTraj  = ", This%NTraj )

          case("Nb of Trajectories per Mini-Batch")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%NTrajBatch
            if (i_Debug_Loc) call Logger%Write( "Nb of Trajectories per Mini-Batch:      This%NTrajBatch  = ", This%NTrajBatch )
            
          
          case("Generate a CPU-Clock Dependent Seed?")
            This%Randomize = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Generate a CPU-Clock Dependent Seed?:      This%Randomize  = ",  This%Randomize )

          case("Initial Random Nb Seed")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%iseed
            if (i_Debug_Loc) call Logger%Write( "Initial random Nb seed:      This%iseed  = ", This%iseed )
            if ((This%Randomize == 'yes') .or. (This%Randomize == 'YES')) then
              CALL RANDOM_SEED(size = n)
              ALLOCATE(seed(n))
              CALL SYSTEM_CLOCK(COUNT=clock)
              seed = clock + 37 * (/ (i - 1, i = 1, n) /)
              CALL RANDOM_SEED(PUT = seed)
              CALL RANDOM_NUMBER(RndNb)
              DEALLOCATE(seed)
            else
              RndNb = Zero
            end if
            This%RndNb = RndNb
            if (i_Debug_Loc) call Logger%Write( "RndNb:      RndNb  = ",  This%RndNb )
            
            
          case("Initial Random Nb Seed for PES Parameters Sampling")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%PESiseed
            if (i_Debug_Loc) call Logger%Write( "Initial Random Nb Seed for PES Parameters Sampling:      This%PESiseed  = ", This%PESiseed )
            
            
          case("Initial Distance between Target and Projectile")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%dinit
            if (i_Debug_Loc) call Logger%Write( "Initial distance between target and projectile:      This%dinit  = ", This%dinit )
            

          case("Nb of PESs")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%NPESs
            if (i_Debug_Loc) call Logger%Write( "Nb of PESs:      This%NPESs  = ", This%NPESs )          
            deallocate( This%PES_Model )
            allocate( This%PES_Model(This%NPESs), stat=Status )
            if (Status/=0) call Error( "Error allocating This%PES_Model" )
            deallocate( This%PES_Degeneracy )
            allocate( This%PES_Degeneracy(This%NPESs), stat=Status )
            if (Status/=0) call Error( "Error allocating This%PES_Degeneracy" )
            This%PES_Degeneracy = One
            deallocate( This%PES_Degeneracy_Vec )
            allocate( This%PES_Degeneracy_Vec(This%NPESs), stat=Status )
            if (Status/=0) call Error( "Error allocating This%PES_Degeneracy_Vec" )
            This%PES_Degeneracy_Vec = Zero
            
          case("PES Parameters Folder Name")
            This%PES_ParamsFldr = line_input(i_eq+2:150)
            This%PES_ParamsFldr = trim(adjustl(This%PES_ParamsFldr))

          case("PES Parameters File Name")
            This%PES_ParamsFile = line_input(i_eq+2:150)
            This%PES_ParamsFile = trim(adjustl(This%PES_ParamsFile))
            
          case("Nb of PESs Samples")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%NPESs
            This%StochPESFlg = .True.
            if (i_Debug_Loc) call Logger%Write( "Found a Stochastic PES! This%StochPESFlg set to .TRUE.!" )
            if (i_Debug_Loc) call Logger%Write( "Nb of Samples from the Stochastic PES:       This%NPES = ", This%NPESs )        
            deallocate( This%PES_Degeneracy )
            allocate( This%PES_Degeneracy(This%NPESs), stat=Status )
            if (Status/=0) call Error( "Error allocating This%PES_Degeneracy" )
            This%PES_Degeneracy = One
            PES_Model_Temp = This%PES_Model(1)
            deallocate( This%PES_Model )
            allocate( This%PES_Model(This%NPESs), stat=Status )
            if (Status/=0) call Error( "Error re-allocating This%PES_Model" )
            This%PES_Model = PES_Model_Temp
            
            
          case("Randomize PES Initial Conditions?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%PESRndIniCondFlg = .true.
            end if 
            if (i_Debug_Loc) call Logger%Write( "Randomize PES Initial Conditions?:       This%PESRndIniCondFlg = ", This%PESRndIniCondFlg )
            
            
          case("Sample Parameters for Stochastic PES?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%SampleParamsStochPES = .true.
            end if 
            if (i_Debug_Loc) call Logger%Write( "Sample Parameters for Stochastic PES?:       This%SampleParamsStochPES = ", This%SampleParamsStochPES )
            

          case("Relative Translational Energy Model")
            This%TtraModel = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Relative Translational Energy Model:      This%TtraModel  = ", This%TtraModel )

!          case("Relative Translational Energy [Eh]")
!            line_input = line_input(i_eq+2:150)
!            READ(line_input, '(d20.10)') This%Ttra
!            if (i_Debug_Loc) call Logger%Write( "Initial Translational Temperature [K]:       This%Ttra = ", This%Ttra )

!          case("Initial Translational Temperature [K]")
!            line_input = line_input(i_eq+2:150)
!            READ(line_input, '(d20.10)') This%Ttra
!            if (i_Debug_Loc) call Logger%Write( "Initial Translational Temperature [K]:       This%Ttra = ", This%Ttra )

          case("Fixed Total Energy [Eh]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%Etot
            if (i_Debug_Loc) call Logger%Write( "Fixed Total Energy [Eh]:       This%Etot = ", This%Etot )
            
            
          case("Initial Impact Parameter Model")
            This%ImpParModel = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Initial Impact Parameter Model:      This%ImpParModel  = ", This%ImpParModel )
            
          case("Strata Type")
            This%ImpParStrataType = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Strata Type:      This%ImpParStrataType  = ", This%ImpParStrataType )

          case("Initial Impact Parameter")
            This%nring  =  1
            allocate( This%bring(This%nring) )
            if (Status/=0) call Error( "Error allocating This%bring" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%bring to This%nring = ", This%nring )
            allocate( This%ntring(This%nring) )
            if (Status/=0) call Error( "Error allocating This%ntring" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%ntring to This%nring = ", This%nring )
            This%ntring(1) =  -1
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%bring(1)
            if (i_Debug_Loc) call Logger%Write( "Initial impact parameter:       This%bring(1) = ", This%bring(1) )

          case("Nb of Stratified Sampling Ranges")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%nring
            if (i_Debug_Loc) call Logger%Write( "Nb of Stratified Sampling Ranges:       This%nring = ", This%nring )
            allocate( This%bring(This%nring) )
            if (Status/=0) call Error( "Error allocating This%bring" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%bring to This%nring = ", This%nring )
            allocate( This%ntring(This%nring) )
            if (Status/=0) call Error( "Error allocating This%ntring" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%ntring to This%nring = ", This%nring )

          case("Total Angular Momentum")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TotalAngularMomentum
            if (i_Debug_Loc) call Logger%Write( "Total Angular Momentum:       This%TotalAngularMomentum = ", This%TotalAngularMomentum )
            This%nring = 1
            allocate( This%bring( This%nring) ); This%bring  = Zero
            if (Status/=0) call Error( "Error allocating This%bring" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%bring to This%nring = ", This%nring )
            allocate( This%ntring(This%nring) ); This%ntring = -1
            if (Status/=0) call Error( "Error allocating This%ntring" )
            if (i_Debug_Loc) call Logger%Write( "-> Allocating This%ntring to This%nring = ", This%nring )


          case("Method for Initial Vibrational Quantum Nb")
            This%vInMethod = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Method for Initial Vibrational Quantum Nb:      This%vInMethod  = ", This%vInMethod )

          case("Method for Initial Rotational Quantum Nb")
            This%jInMethod = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Method for Initial Rotational Quantum Nb:      This%jInMethod  = ", This%jInMethod )

          case("All / Only Even / Only Odd values for Rotational Quantum Nb")
            This%ioddChar = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "All / Only Even / Only Odd values for Rotational Quantum Nb:      This%ioddChar  = ", This%ioddChar )
            if (trim(adjustl(This%ioddChar)) == "All") then
              This%iodd = 0
            elseif (trim(adjustl(This%ioddChar)) == "Only Even") then
              This%iodd = 2
            elseif (trim(adjustl(This%ioddChar)) == "Only Odd") then
              This%iodd = 1
            end if
            if (i_Debug_Loc) call Logger%Write( "All / Only Even / Only Odd values for Rotational Quantum Nb:      This%iodd  = ", This%iodd )

          case("Initial Vibrational Quantum Nb")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%vIn
            if (i_Debug_Loc) call Logger%Write( "Initial Vibrational Quantum Nb:      This%vIn  = ", This%vIn )

          case("Initial Rotational Quantum Nb")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%jIn
            if (i_Debug_Loc) call Logger%Write( "Initial Rotational Quantum Nb:      This%jIn  = ", This%jIn )

          case("Unit for the Arbitraty Distribution")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I20)') This%UnitArbDist
            if (i_Debug_Loc) call Logger%Write( "Initial Rotational Quantum Nb:      This%UnitArbDist  = ", This%UnitArbDist )
            
          
          case("Nb of Quadrature Pts for Action Integral")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%nquad
            if (i_Debug_Loc) call Logger%Write( "Nb of Quadrature Pts for Action Integral:       This%nquad = ", This%nquad )

          case("Lower Bound for Position of Centrifugal Barrier Max")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%rsmal
            if (i_Debug_Loc) call Logger%Write( "Lower Bound for Position of Centrifugal Barrier Max:       This%rsmal = ", This%rsmal )

          case("Upper Bound for Position of Centrifugal Barrier Max")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%rcent
            if (i_Debug_Loc) call Logger%Write( "Upper Bound for Position of Centrifugal Barrier Max:      This%rcent  = ", This%rcent )

          case("Factor for q.b. Levels")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%rfact
            if (i_Debug_Loc) call Logger%Write( "Factor for q.b. Levels:      This%rfact  = ", This%rfact )

          case("Time threshold for q.b. Levels [Eh]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%tthres
            if (i_Debug_Loc) call Logger%Write( "Time threshold for q.b. Levels [Eh]:      This%tthres    = ", This%tthres )
            This%tthresau  =   This%tthres * 1.0E-12_rkp / autime_to_sec
            if (i_Debug_Loc) call Logger%Write( "Time threshold for q.b. Levels [au]:      This%tthresau  = ", This%tthresau )

          case("Time threshold for q.b. Levels [Au]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%tthresau
            if (i_Debug_Loc) call Logger%Write( "Time threshold for q.b. Levels [au]:      This%tthresau  = ", This%tthresau )
            
            
          case("Write the Initial Parameters File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%ParamsFlg = .true.
            end if 
            if (i_Debug_Loc) call Logger%Write( "Write the Initial Parameters File?:                 This%ParamsFlg    = ", This%ParamsFlg )
            
          case("Write the Progress File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%ProgressFlg = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the Progress File?:                           This%ProgressFlg  = ", This%ProgressFlg )
            
          case("Write the PaQ Solution File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%PaQSolFlg = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the PaQ Solution File?:                       This%PaQSolFlg    = ", This%PaQSolFlg )
            
          case("Write the PaQ Evolution File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%PaQEvoFlg = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the PaQ Evolution File?:                      This%PaQEvoFlg    = ", This%PaQEvoFlg )
            
          case("Write the PES Evolution File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%PESEvoFlg = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the PES Evolution File?:                      This%PESEvoFlg    = ", This%PESEvoFlg )
            
          case("Write the Spatial Coordinates Evolution File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%XXEvoFlg         = .true.
              This%XXEvoSnglFileFlg = .false.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the Spatial Coordinates Evolution File?:      This%XXEvoFlg     = ", This%XXEvoFlg )

          case("Write the Spatial Coordinates Evolution in a Single File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%XXEvoFlg         = .true.
              This%XXEvoSnglFileFlg = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the Spatial Coordinates Evolution in a Single File?:      This%XXEvoSnglFileFlg     = ", This%XXEvoSnglFileFlg )
            
          case("Write the Energy Evolution File?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%EnEvoFlg = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Write the Energy Evolution File?:                   This%EnEvoFlg     = ", This%EnEvoFlg )

          
          case("CPU Time Limit ")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TimeMax
            if (i_Debug_Loc) call Logger%Write( "CPU Time Max:      This%TimeMax  = ", This%TimeMax )

          case("Restart?")
            if ((trim(line_input(i_eq+2:150)) .eq. "yes") .or. (trim(line_input(i_eq+2:150)) .eq. "YES")) then
              This%Restart=.true.
            elseif ((trim(line_input(i_eq+2:150)) .eq. "no") .or. (trim(line_input(i_eq+2:150)) .eq. "NO")) then
              This%Restart=.false.
            else
              This%Restart=.false.
              if (i_Debug_Loc) call Logger%Write( "Not Received YES / NO Command for This%Restart. This%Restart set to ", This%Restart )
            end if
            if (i_Debug_Loc) call Logger%Write( "Restart Driver?:      This%Restart  = ", This%Restart )
            
            
          case("Backward Integrating?")
            if ((adjustl(trim(line_input(i_eq+2:150))) == 'YES') .or. (adjustl(trim(line_input(i_eq+2:150))) == 'yes')) then
              This%BwrdIntegrationFlg = .TRUE.
            end if
            if (i_Debug_Loc) call Logger%Write( "Backward Integrating?:      This%BwrdIntegrationFlg  = ", This%BwrdIntegrationFlg )

          case("Backward Integration Frequency")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BwrdIntegrationFreq
            if (i_Debug_Loc) call Logger%Write( "Backward Integration Frequency:      This%BwrdIntegrationFreq  = ", This%BwrdIntegrationFreq )
            
            
          case("Initial time step")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%dt
            if (i_Debug_Loc) call Logger%Write( "Initial time step:      This%dt  = ", This%dt )
            
          case("Error control parameter")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%eps
            if (i_Debug_Loc) call Logger%Write( "Error control parameter:      This%eps  = ", This%eps )
            
          case("Stepsize increase damping factor")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%Relax
            if (i_Debug_Loc) call Logger%Write( "Stepsize increase damping factor:      This%Relax  = ", This%Relax )
            
          case("Nb of extrapolation steps")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NSteps
            if (i_Debug_Loc) call Logger%Write( "Nb of extrapolation steps:      This%NSteps  = ", This%NSteps )
            
          case("Try to increase stepsize?")
            This%IncStpSzChar = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Try to increase stepsize?:      This%IncStpSzChar  = ", This%IncStpSzChar )
            if (((trim(adjustl(This%IncStpSzChar))) .eq. 'no') .or. ((trim(adjustl(This%IncStpSzChar))) .eq. 'NO')) This%IncStpSzFlg = .false.
            if (i_Debug_Loc) call Logger%Write( "Try to increase stepsize?:      This%IncStpSzFlg  = ", This%IncStpSzFlg )
            
          case("Nb of integration between end checking")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%ncall
            if (i_Debug_Loc) call Logger%Write( "Nb of integration between end checking:      This%ncall  = ", This%ncall )
            
          case("Max distance between atom pair [bohr]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%rmax
            if (i_Debug_Loc) call Logger%Write( "Max distance between atom pair [bohr]:      This%rmax  = ", This%rmax )

          case("Max trajectory time [a.u.]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%tmax
            if (i_Debug_Loc) call Logger%Write( "Max trajectory time [a.u.]:       This%tmax = ", This%tmax )
            
            
          case("Distinguish Exchanges?")
            call Error("'Distinguish Exchanges?' is OBSOLETE. Please, change INPUT-FILE with 'Distinguish Exchanges between Each Other?' and/or 'Distinguish Exchanges from Inelastic?'")
            ! This%ConsiderExc = line_input(i_eq+2:150)
            ! if (i_Debug_Loc) call Logger%Write( "Distinguish Exchanges?:      This%ConsiderExc  = ", This%ConsiderExc )
            ! if (((trim(adjustl(This%ConsiderExc))) .eq. 'yes') .or. ((trim(adjustl(This%ConsiderExc))) .eq. 'YES')) This%ConsiderExcFlg = .true.
            ! if (i_Debug_Loc) call Logger%Write( "Distinguish Exchanges?:      This%ConsiderExcFlg  = ", This%ConsiderExcFlg )

        case("Distinguish Exchanges between Each Other?")
            line_input = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Distinguish Exchanges between Each Other?:      This%MergeExchsFlg  = ", This%MergeExchsFlg )
            if (((trim(adjustl(line_input))) .eq. 'no') .or. ((trim(adjustl(line_input))) .eq. 'NO')) This%MergeExchsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Distinguish Exchanges between Each Other?:      This%MergeExchsFlg  = ", This%MergeExchsFlg )

        case("Distinguish Exchanges from Inelastic?")
            line_input = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Distinguish Exchanges from Inelastic?:      This%MergeExchToInelFlg  = ", This%MergeExchToInelFlg )
            if (((trim(adjustl(line_input))) .eq. 'no') .or. ((trim(adjustl(line_input))) .eq. 'NO')) This%MergeExchToInelFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Distinguish Exchanges from Inelastic?:      This%MergeExchToInelFlg  = ", This%MergeExchToInelFlg )

!          case("Distinguish Quasi-Bound?")
!            This%ConsiderQb = line_input(i_eq+2:150)
!            if (i_Debug_Loc) call Logger%Write( "Distinguish Quasi-Bound?:      This%ConsiderQb  = ", This%ConsiderQb )
          
          
          case("Nb of Final Conditions on the Trajectories")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NCond
            if (i_Debug_Loc) call Logger%Write( "Nb of Conditions to be Distinguished:      This%NCond  = ", This%NCond )
            allocate( This%iskip(This%NCond), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%iskip" )
            allocate( This%PresEvOdd(This%NCond), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%PresEvOdd" )
            This%PresEvOdd = .false.
            
            
          case("Identical Diatoms?")
            line_input = line_input(i_eq+2:150)
            if ((adjustl(trim(line_input)) == 'yes') .or. (adjustl(trim(line_input)) == 'YES')) then
              This%IdenticalDiatoms = .true.
            end if
            if (i_Debug_Loc) call Logger%Write( "Identical Diatoms?:      This%ConsiderQb  = ", This%ConsiderQb )
            
            
          case("Final States Assigment Mthd")
            This%AssignmentMethod = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Final States Assigment Mthd:      This%AssignmentMethod  = ", This%AssignmentMethod )


          case("Writing Trajectory Files in Binary Format?")
            line_input = line_input(i_eq+2:150)
            if ((adjustl(trim(line_input)) == 'yes') .or. (adjustl(trim(line_input)) == 'YES')) then
              This%StatWritesBinaryFlg = .True.
            end if


!          case("Nb Max of Trajectories per Level / Bin to use for Stastics")
!            line_input = line_input(i_eq+2:150)
!            READ(line_input, '(I20)') This%NTrajectoriesToAnalyze
!            if (i_Debug_Loc) call Logger%Write( "Nb Max of Trajectories to use for Stastics:      This%NTrajectoriesToAnalyze   = ", This%NTrajectoriesToAnalyze  )
      
      
          case("KONIG Inital Temperature [K]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TInit
            if (i_Debug_Loc) call Logger%Write( "KONIG Inital Temperature [K]:      This%TInit  = ", This%TInit )
            write(This%TInit_char,"(I10)") floor(This%TInit)
          
          case("Temperature for Distribution Function Computation [K]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TInit
            if (i_Debug_Loc) call Logger%Write( "KONIG Inital Temperature [K]:      This%TInit  = ", This%TInit )
            write(This%TInit_char,"(I10)") floor(This%TInit) 


         case("Flag for Running External Code")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%RunExtCodeIntFlg
            if (i_Debug_Loc) call Logger%Write( "KONIG Inital Temperature [K]:      This%RunExtCodeIntFlg  = ", This%RunExtCodeIntFlg )
      
        end select
        
        
        do iMol = 1,This%NMolecules
          write(iMol_char, "(I2)") iMol

          Mol_case = "Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(Mol_case)) then
            This%Molecules_Name(iMol) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Molecule Nb", iMol, " Name:      This%Molecules_Name(i)  = ", This%Molecules_Name(iMol) )
          end if
          
          
          LevelsFileName_case = "Levels File Name, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(LevelsFileName_case)) then
            This%LevelsFileName(iMol) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Levels Original File Name, Molecule Nb", iMol, ":      This%LevelsFileName(i)  = ", This%LevelsFileName(iMol) )
          end if


          Ecut_case = "Min Energy for Cutting Levels, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(Ecut_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%Ecut(iMol)
            if (i_Debug_Loc) call Logger%Write( "Min Energy for Cut, Molecule Nb", iMol, ":      This%Ecut(i)  = ", This%Ecut(iMol), " Eh" )
          end if
          
          Tcut_case = "Min Time   for Cutting Levels, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(Tcut_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%Tcut(iMol)
            if (i_Debug_Loc) call Logger%Write( "Min Time for Cut, Molecule Nb", iMol, ":      This%Tcut(i)  = ", This%Tcut(iMol), " ps" )
          end if
          

          Molecule_NAtoms_case = "Nb of Atoms, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(Molecule_NAtoms_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%Molecule_NAtoms(iMol)
            if (i_Debug_Loc) call Logger%Write( "Nb of Atoms, Molecule Nb", iMol, ":      This%Molecule_NAtoms(i)  = ", This%Molecule_NAtoms(iMol) )
          end if
          
          
          BSortMethod_case = "Sorting Method, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(BSortMethod_case)) then
            This%BSortMethod(iMol) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Sorting Method, Molecule Nb", iMol, ":      This%BSortMethod(i)  = ", This%BSortMethod(iMol) )
          end if
          
          MinLevel_case = "Starting Level, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(MinLevel_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BSortInfo(iMol,1)
            if (i_Debug_Loc) call Logger%Write( "Starting Level, Molecule Nb", iMol, ":      This%BSortInfo(iMol,1)  = ", This%BSortInfo(iMol,1) )
          end if
          
          MaxLevel_case = "Final    Level, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(MaxLevel_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BSortInfo(iMol,2)
            if (i_Debug_Loc) call Logger%Write( "Final Level, Molecule Nb", iMol, ":      This%BSortInfo(iMol,2)  = ", This%BSortInfo(iMol,2) )
          end if


          MinLevel_case = "Starting Bin, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(MinLevel_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BSortInfo(iMol,2)
            if (i_Debug_Loc) call Logger%Write( "Starting Level, Molecule Nb", iMol, ":      This%BSortInfo(iMol,2)  = ", This%BSortInfo(iMol,2) )
          end if
          
          MaxLevel_case = "Final    Bin, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(MaxLevel_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BSortInfo(iMol,3)
            if (i_Debug_Loc) call Logger%Write( "Final Level, Molecule Nb", iMol, ":      This%BSortInfo(iMol,3)  = ", This%BSortInfo(iMol,3) )
          end if


          Minvqn_case = "Starting vqn, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(Minvqn_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BSortInfo(iMol,1)
            if (i_Debug_Loc) call Logger%Write( "Starting vqn, Molecule Nb", iMol, ":      This%BSortInfo(iMol,1)  = ", This%BSortInfo(iMol,1) )
          end if
          
          Maxvqn_case = "Final    vqn, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(Maxvqn_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%BSortInfo(iMol,2)
            if (i_Debug_Loc) call Logger%Write( "Final vqn, Molecule Nb", iMol, ":      This%BSortInfo(iMol,2)  = ", This%BSortInfo(iMol,2) )
          end if

          BinDataFile_case = "File with Bins Energy Minima, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(BinDataFile_case)) then
            This%BinDataFile(iMol,1) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( " File with Bins Energy Minima, Molecule Nb", iMol, ":      This%BinDataFile(i,1)  = ", This%BinDataFile(iMol,1) )
          end if
          
          BinDataFile_case = "File with Level-Bin Mapping, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(BinDataFile_case)) then
            This%BinDataFile(iMol,1) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( " File with Bins Energy Minima, Molecule Nb", iMol, ":      This%BinDataFile(i,1)  = ", This%BinDataFile(iMol,1) )
          end if
          
          BinDataFile_case = "File with Bins QNs Bounds, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(BinDataFile_case)) then
            This%BinDataFile(iMol,2) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( " File with QNs Bounds, Molecule Nb", iMol, ":      This%BinDataFile(i,2)  = ", This%BinDataFile(iMol,2) )
          end if

          NBins_case = "Nb of Bins, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(NBins_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NBins(iMol)
            WRITE(This%NBins_char(iMol), '(I6)') This%NBins(iMol)
            if (i_Debug_Loc) call Logger%Write( "Nb of Bins, Molecule Nb", iMol, ":      This%NBins(i)  = ", This%NBins(iMol), " = ", This%NBins_char(iMol) )
            
            This%BSortInfo(iMol,1) = This%NBins(iMol)
            This%BSortInfo(iMol,2) = 0
          end if
          
          DiatPot_case = "Diatomic Potential Model, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(DiatPot_case)) then
            This%Diatomic_Model(iMol) = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Diatomic Potential Model:      This%Diatomic_Model(iMol)  = ", This%Diatomic_Model(iMol) )
          end if

          DiatPot_File_case = "Diatomic Potential Parameters File Name, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(DiatPot_File_case)) then
            This%DiatPot_ParamsFile(iMol) = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Diatomic Potential Parameters File Name:      This%DiatPot_ParamsFile(iMol)  = ", This%DiatPot_ParamsFile(iMol) )
          end if
          
        end do
        
        
        do iAt = 1,This%NAtoms
          write(iAt_char, "(I2)") iAt

          At_case = "Atom" // iAt_char
          if (adjustl(trim(i_case)) == TRIM(At_case)) then
            This%AtomsName(iAt) = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Atom Nb", iAt, " Name:      This%AtomsName(i)  = ", This%AtomsName(iAt) )
          end if

          AtMass_case1 = "Mass [a.u.] of Atom" // iAt_char
          AtMass_case2 = "Mass [a.u.], Atom" // iAt_char
          if ((adjustl(trim(i_case)) == TRIM(AtMass_case1)) .or. (adjustl(trim(i_case)) == TRIM(AtMass_case2)) ) then
            line_input = line_input(i_eq+2:150)
            if (trim(adjustl(line_input)) == 'FromDataFile') then
                This%AtomsMass(iAt) = - 1.d0
            else
                READ(line_input, '(d20.10)') This%AtomsMass(iAt)
            end if
            if (i_Debug_Loc) call Logger%Write( "Mass of Atom Nb", iAt, ":      This%AtomsMass(i)  = ", This%AtomsMass(iAt) )
          end if

          AtMassFile_case = "Mass File, Atom" // iAt_char
          if (adjustl(trim(i_case)) == TRIM(AtMassFile_case)) then
            This%FileForAtomMass(iAt) = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Atom Nb", iAt, " File for Mass:      This%FileForAtomMass(i)  = ", This%FileForAtomMass(iAt) )
          end if

        end do
        
        
        do iSp = 1,This%NSpecies
          write(iSp_char, "(I2)") iSp

          Sp_case = "Species" // iSp_char
          if (adjustl(trim(i_case)) == TRIM(Sp_case)) then
            This%SpeciesName(iSp) = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Species Nb", iSp, " Name:      This%SpeciesName(i)  = ", This%SpeciesName(iSp) )
          end if

        end do
        
        
        do iTint = 1,This%NTint
          write(iTint_char, "(I2)") iTint

          Tint_case     = "Internal Temperature [K]" // iTint_char
          Tint_case_bis = "Internal Temperature"     // iTint_char
          if ( (adjustl(trim(i_case)) == TRIM(Tint_case)) .or. (adjustl(trim(i_case)) == TRIM(Tint_case_bis))) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TintVec(iTint)
            if (i_Debug_Loc) call Logger%Write( "Internal Temperature Nb", iTint, " = ", This%TintVec(iTint) )
          end if

        end do
        
        
        do iTtra = 1,This%NTtra
          write(iTtra_char, "(I2)") iTtra

          Ttra_case     = "Translational Temperature [K]" // iTtra_char
          Ttra_case_bis = "Translational Temperature"     // iTtra_char
          if ( (adjustl(trim(i_case)) == TRIM(Ttra_case)) .or. (adjustl(trim(i_case)) == TRIM(Ttra_case_bis)) ) then
            This%TtraVecChar(iTtra) = adjustl(trim(line_input(i_eq+2:150)))
            READ(This%TtraVecChar(iTtra), '(d20.10)') This%TtraVec(iTtra)
            WRITE(This%TtraVecIntChar(iTtra), '(I10)') int(This%TtraVec(iTtra))
            if (i_Debug_Loc) call Logger%Write( "Translational Temperature Nb", iTtra, " = ", This%TtraVec(iTtra), "; This%TtraVecChar(iTtra) = ", This%TtraVecChar(iTtra), "; This%TtraVecIntChar(iTtra) = ", This%TtraVecIntChar(iTtra) )
          end if

        end do
        
        
        do iErel = 1,This%NErel
          write(iErel_char, "(I2)") iErel

          Erel_case = "Translational Energy [Eh]" // iErel_char
          if (adjustl(trim(i_case)) == TRIM(Erel_case)) then
            This%ErelVecChar(iErel) = adjustl(trim(line_input(i_eq+2:150)))
            READ(This%ErelVecChar(iErel), '(d20.10)') This%ErelVec(iErel)
            WRITE(This%ErelVecIntChar(iErel), '(I10)') int(This%ErelVec(iErel))
            if (i_Debug_Loc) call Logger%Write( "Relative Translational Energy Nb", iErel, " = ", This%ErelVec(iErel), "; This%ErelVecChar(iErel) = ", This%ErelVecChar(iErel), "; This%ErelVecIntChar(iErel) = ", This%ErelVecIntChar(iErel) )
          end if

        end do
        
        
        do iring = 1,This%nring
          write(iring_char, "(I2)") iring

          ntring_case = trim("Nb of Trajectories in Range") // iring_char
          if (adjustl(trim(i_case)) == TRIM(ntring_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%ntring(iring)
            if (i_Debug_Loc) call Logger%Write( "Nb of Trajectories in Range Nb", iring, " :      This%ntring(i)  = ", This%ntring(iring) )
          end if

          bring_case = trim("Impact Parameter for Range") // iring_char
          if (adjustl(trim(i_case)) == TRIM(bring_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%bring(iring)
            if (i_Debug_Loc) call Logger%Write( "Impact Parameter for Range Nb", iring, " :      This%bring(i)  = ", This%bring(iring) )
          end if

        end do
        
        
        do iCond = 1,This%NCond
          write(iCond_char, "(I2)") iCond

          iskip_case = "Iskip for Condition" // iCond_char
          if (adjustl(trim(i_case)) == TRIM(iskip_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)')  This%iskip(iCond)
            if (i_Debug_Loc) call Logger%Write( "Iskip for Condition ", iCond, ":      This%iskip(i)  = ", This%iskip(iCond) )
          end if

          PresEvOdd_case = "Preserving Even / Oddness of Condition" // iCond_char
          if (adjustl(trim(i_case)) == TRIM(PresEvOdd_case)) then
            line_input = line_input(i_eq+2:150)
            if ((adjustl(trim(line_input)) == 'yes') .or. (adjustl(trim(line_input)) == 'YES')) then
              This%PresEvOdd(iCond) = .true.
              if (i_Debug_Loc) call Logger%Write( "Preserving Even / Oddness of Condition", iCond )
            end if
          end if

        end do
        
        
        do iiPES = 1,This%NPESs
          write(iiPES_char, "(I2)") iiPES

          PESModel_case = "Model, PES" // iiPES_char
          if (adjustl(trim(i_case)) == TRIM(PESModel_case)) then
            This%PES_Model(iiPES) = trim(line_input(i_eq+2:150))
            if (i_Debug_Loc) call Logger%Write( "Model, PES ", iiPES, ":      This%PES_Model(iiPES)  = ", This%PES_Model(iiPES) )
          end if
          
          PESDegeneracy_case = "Degeneracy, PES" // iiPES_char
          if (adjustl(trim(i_case)) == TRIM(PESDegeneracy_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)')  This%PES_Degeneracy(iiPES)
            if (i_Debug_Loc) call Logger%Write( "Degeneracy, PES ", iiPES, ":      This%PES_Degeneracy(iiPES)  = ", This%PES_Degeneracy(iiPES) )
            This%PES_Degeneracy_Vec(iiPES) = This%PES_Degeneracy_Vec(max(1,iiPES-1)) + This%PES_Degeneracy(iiPES)
          end if

        end do
  
        
        do iHL = 1,2
          write(iHL_char, "(I2)") iHL
          
          NHL_case = "Nb Neurons, Hidden Layer" // iHL_char
          if (adjustl(trim(i_case)) == TRIM(NHL_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)')  This%NHiddenLayerNeurons(iHL)
            if (i_Debug_Loc) call Logger%Write( "Nb Neurons, Hidden Layer ", iHL, ":      This%NHiddenLayerNeurons(iHL)  = ", This%NHiddenLayerNeurons(iHL) )
          end if
          
        end do
                    
      end if

    end if


  end do

  close(Unit)

  
  This%Tcut=This%Tcut*1.e-12_rkp/autime_to_sec


  if (This%StochPESFlg) then
    deallocate( This%PES_Degeneracy_Vec )
    allocate( This%PES_Degeneracy_Vec(This%NPESs), stat=Status )
    if (Status/=0) call Error( "Error allocating This%PES_Degeneracy_Vec" )
    This%PES_Degeneracy_Vec(1) = One / This%NPESs
    do i = 2,This%NPESs
      This%PES_Degeneracy_Vec(i) = This%PES_Degeneracy_Vec(i-1) + One / This%NPESs
    end do
  else
    This%PES_Degeneracy_Vec = This%PES_Degeneracy_Vec / sum(This%PES_Degeneracy,1)
  end if
  if (i_Debug_Loc) call Logger%Write( "Degeneracy Vector:   This%PES_Degeneracy_Vec = ", This%PES_Degeneracy_Vec )
  This%PES_Degeneracy = This%PES_Degeneracy / sum(This%PES_Degeneracy,1)


  if ( trim(This%TtraModel) .eq. "Uniform" ) then
    if (i_Debug_Loc) call Logger%Write( "Relative Translational Energy Model Selected: ", This%TtraModel, " and the related Fixed Relative Translational Energy [Eh] has been defined: This%Ttra = ", This%Erel )
  elseif ( trim(This%TtraModel) .eq. "Boltzmann" ) then
    if (i_Debug_Loc) call Logger%Write( "Relative Translational Energy Model Selected: ", This%TtraModel, " and the related Initial Translational Temperature [K] has been defined: This%Ttra = ", This%Ttra )
  elseif ( trim(This%TtraModel) .eq. "FixedTotEn" ) then
    if (i_Debug_Loc) call Logger%Write( "Relative Translational Energy Model Selected: ", This%TtraModel, " and the related Fixed Total Energy [Eh] has been defined: This%Ttra = ", This%Etot )
  else
    if (i_Debug_Loc) call Logger%Write( "WARNING: No Relative Translational Energy Model Selected" )
  end if

  do iMol = 1,This%NMolecules
  
    if (This%NBins(iMol) == 0) then
      if (i_Debug_Loc) call Logger%Write( "Found the ", iMol, "-th Molecule with 0 Levels / Bins. Going to read the overall Nb of Levels.")
            
      FileName = trim(adjustl(This%OutputDir)) // '/' // trim(adjustl(This%System)) // '/' // trim(adjustl(This%Molecules_Name(iMol))) // '/NLevels.inp '
      if (i_Debug_Loc) call Logger%Write( "Opening File: ", FileName )
      open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )      
        if (Status==0) then 
          read(Unit,'(I10)',iostat=Status) This%NBins(iMol)
      close(Unit)      
        else
      close(Unit)
          FileName = trim(adjustl(This%OutputDir)) // '/' // trim(adjustl(This%System)) // '/' // trim(adjustl(This%Molecules_Name(iMol))) // '/levels_cut.inp '
          INQUIRE( FILE=trim(adjustl(FileName)), EXIST=ExFlg )
          if (ExFlg) then
            open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
              if ((Status/=0) .and. (i_Debug_Loc)) call Logger%Write( trim(adjustl(This%OutputDir)) // "/" // trim(adjustl(This%System)) // "/" // trim(adjustl(This%Molecules_Name(iMol))) // "/levels_cut.inp NOT FOUND. Nb of Levels set to 0 (For Now)!" )
              NLevels = 0
              do
                read(Unit,'(A100)',iostat=Status) Line
                if (Line(2:2) /= '#') then
                  NLevels = NLevels + 1
                end if
                if (Status /= 0) exit
              end do
            close(Unit)
            NLevels = NLevels - 1
            This%NBins(iMol) = NLevels   
          end if
        end if
      
      WRITE(This%NBins_char(iMol), '(I6)') This%NBins(iMol)
      if (i_Debug_Loc) call Logger%Write( "Nb of Bins, Molecule Nb", iMol, ":      This%NBins(i)  = ", This%NBins(iMol), " = ", This%NBins_char(iMol) )

    end if
    
  end do

  do iAt = 1,This%NAtoms
    if (This%AtomsMass(iAt) == Zero) then
      call Error( "ERROR! Found an Atom Mass Equal to 0.0; the Mass has not been specified!" )
    elseif (This%AtomsMass(iAt) < Zero) then
      if (i_Debug_Loc) call Logger%Write( "Found a Negative Mass for Atom ", iAt, "; I am going to read the Mass from Data-File: ", adjustl(trim(This%FileForAtomMass(iAt))) )
      call This%ReadAtomMass(iAt, i_Debug=i_Debug_Loc)
    end if
  end do 

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine WriteBashInputVariables( This )

  class(Input_Type)                           ,intent(inout)  ::    This
  
  character(:)  ,allocatable                                  ::    FileName
  integer                                                     ::    iTint, iTtra, iErel, iMol
  integer                                                     ::    Unit!, NewUnit
  integer                                                     ::    len
  integer                                                     ::    iTemp
  integer                                                     ::    status

  FileName = trim(adjustl(This%OutputDir)) // '/InputForBash.inp'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  
    write(Unit,'(A20)') trim(adjustl(This%System))
          
    if ( trim(This%TtraModel) .eq. "Uniform" ) then
      write(Unit,'(I20)') 0
      write(Unit,'(I20)') This%NErel
      do iErel = 1,This%NErel
        write(Unit,'(f20.5)') This%ErelVec(iErel)
      end do
    elseif ( trim(This%TtraModel) .eq. "Boltzmann" ) then
      write(Unit,'(I20)') 1
      write(Unit,'(I20)') This%NTtra
      do iTtra = 1,This%NTtra
        write(Unit,'(f20.1)') This%TtraVec(iTtra)
      end do
    end if
    
    write(Unit,'(I20)') This%NTint
    if (This%NTint > 1) then
      do iTint = 1,This%NTint
        write(Unit,'(f20.1)') This%TintVec(iTint)
      end do
    end if
    
    write(Unit,'(I20)') This%NInitMolecules
    do iMol = 1,This%NInitMolecules
      write(Unit,'(A20)') trim(adjustl(This%Molecules_Name(iMol)))
      write(Unit,'(A20)') trim(adjustl(This%BSortMethod(iMol)))
      write(Unit,'(I20)') This%BSortInfo(iMol,1)
      write(Unit,'(I20)') This%BSortInfo(iMol,2)
      write(Unit,'(I20)') This%BSortInfo(iMol,3)
    end do
    
    if ((This%GenLev == 'yes') .or. (This%GenLev == 'YES')) then
      write(Unit,'(I20)') 1
    else
      write(Unit,'(I20)') 0
    end if
    
    if ((This%RunPrep == 'yes') .or. (This%RunPrep == 'YES')) then
      write(Unit,'(I20)') 1
    else
      write(Unit,'(I20)') 0
    end if
    
    if ((This%RunTraj == 'yes') .or. (This%RunTraj == 'YES')) then
      write(Unit,'(I20)') 1
    else
      write(Unit,'(I20)') 0
    end if
    
    if ((This%RunPost == 'yes') .or. (This%RunPost == 'YES')) then
      write(Unit,'(I20)') 1
    else
      write(Unit,'(I20)') 0
    end if
    
    if ((This%CompQuant == 'yes') .or. (This%CompQuant == 'YES')) then
      write(Unit,'(I20)') 1
    else
      write(Unit,'(I20)') 0
    end if
    
    if (This%StochPESFlg) then
      write(Unit,'(I20)') 1
    else
      write(Unit,'(I20)') 0
    end if
    
    write(Unit,'(I20)') This%NPESs
    
    write(Unit,'(I20)') This%PESiseed

    write(Unit,'(I20)') This%RunExtCodeIntFlg

  close(Unit)

End Subroutine


Subroutine GenerateLevels( This, i_Debug)

  class(Input_Type)                         ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, Unit
  integer                                                   ::    Status
  character(:)  ,allocatable                                ::    FileName
  integer                                                   ::    NLevels
  character(100)                                            ::    Line

  character(150)                                            ::    i_case
  character(150)                                            ::    line_input
  character(150)                                            ::    line_input_temp
  integer                                                   ::    i_eq

  integer                                                   ::    iMol
  character(2)                                              ::    iMol_char
  
  character(:) ,allocatable                                 ::    GeneratedLevelsFile_case
  character(:) ,allocatable                                 ::    ComputeLevels_case

  integer                                                   ::    iPsi
  character(2)                                              ::    iPsi_char
  character(150)                                            ::    Psi_case
  character(150)                                            ::    PsiValue_case

  integer                                                   ::    len
  logical                                                   ::    i_Debug_Loc = .true.


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "GenerateLevels" )
  !i_Debug_Loc   =     Logger%On()


  allocate(This%ComputeLevels(This%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%ComputeLevels" )
  This%ComputeLevels=0
  
  allocate(This%GeneratedLevelsFile(This%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%GeneratedLevelsFile" )


!   --------------------------------------------------------------------------------------------------------------
!   --------------------------------------------------------------------------------------------------------------
  FileName = adjustl(trim(This%InputDir)) // '/GenerateLevels.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading the Input file for GenerateLevels" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )

  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )

!   --------------------------------------------------------------------------------------------------------------
  do
    read(Unit,'(A150)',iostat=Status) line_input
    if (Status /= 0) then
      exit
    else
      line_input_temp=trim(adjustl(line_input))
      if ( (line_input_temp(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
        i_eq = 1
        do
          if (line_input(i_eq:i_eq) == '=') exit
          i_eq = i_eq + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_eq-2))))

        select case (adjustl((TRIM(i_case))))


          case("Left Extreme in the Integration Grid")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%xExtremes(1)
            if (i_Debug_Loc) call Logger%Write( "Left Extreme in the Integration Grid:      This%xExtremes(1)  = ",  This%xExtremes(1) )

          case("Right Extreme in the Integration Grid")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%xExtremes(2)
            if (i_Debug_Loc) call Logger%Write( "Right Extreme in the Integration Grid:      This%xExtremes(1)  = ",  This%xExtremes(2) )

          case("Nb of Points in the Integration Grid")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NGrid
            if (i_Debug_Loc) call Logger%Write( "Nb of Points in the Integration Grid:      This%NGrid  = ",  This%NGrid )

          case("Spatial Delta in the Integration Grid")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%DeltaGrid
            if (i_Debug_Loc) call Logger%Write( "Spatial Delta in the Integration Grid:      This%DeltaGrid  = ",  This%DeltaGrid )

          case("Min Value for Energy in the Energy Level Search")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%EStart
            if (i_Debug_Loc) call Logger%Write( "Min Value for Energy in the Energy Level Search:      This%EStart  = ",  This%EStart )

          case("Delta Energy for the Energy Level Search")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%DeltaE
            if (i_Debug_Loc) call Logger%Write( "Delta Energy for the Energy Level Search:      This%EStart  = ",  This%DeltaE )

          case("Energy Tolerance for the Energy Level Search")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%EEpsilon
            if (i_Debug_Loc) call Logger%Write( "Energy Tolerance for the Energy Level Search:      This%EStart  = ",  This%EEpsilon )

          case("Max Vibrational Quantum Nb")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%ivqnMax
            if (i_Debug_Loc) call Logger%Write( "Max Vibrational Quantum Nb:      This%ivqnMax  = ",  This%ivqnMax )

          case("Max Rotational Quantum Nb")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%ijqnMax
            if (i_Debug_Loc) call Logger%Write( "Max Rotational Quantum Nb:      This%ijqnMax  = ",  This%ijqnMax )

          case("Min Vibrational Quantum Nb")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%ivqnMin
            if (i_Debug_Loc) call Logger%Write( "Min Vibrational Quantum Nb:      This%ivqnMin  = ",  This%ivqnMin )

          case("Min Rotational Quantum Nb")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%ijqnMin
            if (i_Debug_Loc) call Logger%Write( "Min Rotational Quantum Nb:      This%ijqnMin  = ",  This%ijqnMin )

          case("Print Levels?")
            line_input = line_input(i_eq+2:150)
            if ((adjustl(trim(line_input)) == 'YES') .or. (adjustl(trim(line_input)) == 'yes')) then
              This%PrintLevelsFlg = .True.
            end if
            if (i_Debug_Loc) call Logger%Write( "New Levels' File Name:      This%PrintLevelsFlg  = ",  This%PrintLevelsFlg )

          case("Sort Levels?")
            line_input = line_input(i_eq+2:150)
            if ((adjustl(trim(line_input)) == 'YES') .or. (adjustl(trim(line_input)) == 'yes')) then
              This%SortLevelsFlg = .True.
            end if
            if (i_Debug_Loc) call Logger%Write( "New Levels' File Name:      This%SortLevelsFlg  = ",  This%SortLevelsFlg )


        end select
        

        do iMol = 1,This%NMolecules
          write(iMol_char, "(I2)") iMol
          
          ComputeLevels_case = "Compute Levels, Molecule " // adjustl(trim(iMol_char))
          if (adjustl(trim(i_case)) == TRIM(ComputeLevels_case)) then
            line_input = line_input(i_eq+2:150)
            if ((adjustl(trim(line_input)) == 'YES') .or. (adjustl(trim(line_input)) == 'yes')) then
              This%ComputeLevels(iMol) = iMol
              if (i_Debug_Loc) call Logger%Write( "Molecule Nb", iMol, ":      Computing Energy Levels." )
            end if
          end if
          
          GeneratedLevelsFile_case = "New Levels File Name, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(GeneratedLevelsFile_case)) then
            This%GeneratedLevelsFile(iMol) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "New Levels File Name, Molecule Nb", iMol, ":      This%GeneratedLevelsFile  = ",  This%GeneratedLevelsFile(iMol) )
          end if

        end do
        
        
        do iPsi = 1,3
          write(iPsi_char, "(I2)") iPsi

          Psi_case = "Location of Boundary Condition, Point Nb" // iPsi_char
          if (adjustl(trim(i_case)) == TRIM(Psi_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%BoundayConditions(2,iPsi)
            if (i_Debug_Loc) call Logger%Write( "Boundary Conditions are set for Point Nb", This%BoundayConditions(2,iPsi), " of the Integration Grid." )
          end if

          PsiValue_case = "Value of Boundary Condition, Point Nb" // iPsi_char
          if (adjustl(trim(i_case)) == TRIM(PsiValue_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%BoundayConditions(1,iPsi)
            if (i_Debug_Loc) call Logger%Write( "Boundary Conditions on Point Nb", This%BoundayConditions(2,iPsi), " of the Integration Grid are This%BoundayConditions(1,iPsi) ", This%BoundayConditions(1,iPsi)  )
          end if

        end do


      end if

    end if


  end do

  close(Unit)

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine DeriveQuantities( This, i_Debug)

  class(Input_Type)                         ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, Unit
  integer                                                   ::    Status
  character(:)  ,allocatable                                ::    FileName
  integer                                                   ::    NLevels
  character(100)                                            ::    Line

  character(150)                                            ::    i_case
  character(150)                                            ::    line_input
  character(150)                                            ::    line_input_temp
  integer                                                   ::    i_eq

  integer                                                   ::    iMol
  character(2)                                              ::    iMol_char

  integer                                                   ::    iComp
  character(3)                                              ::    iComp_char
  character(:) ,allocatable                                 ::    iCompName_case
  character(:) ,allocatable                                 ::    iCompMolFrac_case

  integer                                                   ::    iiPES
  character(2)                                              ::    iiPES_char
  character(150)                                            ::    iiPES_case
  character(150)                                            ::    iiPESWeight_case

  character(150)                                            ::    NBinsCG_case

  integer                                                   ::    len
  logical                                                   ::    i_Debug_Loc = .true.


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "DeriveQuantities" )
  !i_Debug_Loc   =     Logger%On()


  This%NCond           = 0
  This%NPESsToMerge    = 1


  allocate(This%NBinsCG(This%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%NBinsCG" )
  
  allocate(This%NBinsCG_char(This%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%NBinsCG_char" )


!   --------------------------------------------------------------------------------------------------------------
!   --------------------------------------------------------------------------------------------------------------
  FileName = adjustl(trim(This%InputDir)) // '/DeriveQuantities.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading the Input file for DeriveQuantities" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )

  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )

!   --------------------------------------------------------------------------------------------------------------
  do
    read(Unit,'(A150)',iostat=Status) line_input
    if (Status /= 0) then
      exit
    else
      line_input_temp=trim(adjustl(line_input))
      if ( (line_input_temp(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
        i_eq = 1
        do
          if (line_input(i_eq:i_eq) == '=') exit
          i_eq = i_eq + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_eq-2))))

        select case (adjustl((TRIM(i_case))))


          case("Writing Arrhenius for KONIG?")
            This%WriteKONIGChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteKONIGChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteKONIGChar))) .eq. 'YES')) This%WriteKONIGFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Arrhenius for KONIG?:      This%WriteKONIGFlg  = ", This%WriteKONIGFlg )
    
          case("Writing Arrhenius for HEGEL?")
            This%WriteHEGELChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteHEGELChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteHEGELChar))) .eq. 'YES')) This%WriteHEGELFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Arrhenius for HEGEL?:      This%WriteHEGELFlg  = ", This%WriteHEGELFlg )
            
          case("Running KONIG?")
            This%RunKONIGChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%RunKONIGChar))) .eq. 'yes') .or. ((trim(adjustl(This%RunKONIGChar))) .eq. 'YES')) This%RunKONIGFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Running KONIG?:      This%RunKONIGFlg  = ", This%RunKONIGFlg )

          case("Running HEGEL?")
            This%RunHEGELChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%RunHEGELChar))) .eq. 'yes') .or. ((trim(adjustl(This%RunHEGELChar))) .eq. 'YES')) This%RunHEGELFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Running HEGEL?:      This%RunHEGELFlg  = ", This%RunHEGELFlg )

        
          case("Nb of PESs to Merge")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NPESs
            if (i_Debug_Loc) call Logger%Write( "Nb of PESs to Merge:      This%NPESs  = ", This%NPESs )
            allocate( This%ResultsFile(This%NPESs), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%ResultsFile" )
            allocate( This%PESWeight(This%NPESs), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%PESWeight" )
            This%PESWeight = One
            
          case("Min Value for Bin Rate Coefficient")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MinForRates
            if (i_Debug_Loc) call Logger%Write( "Min Value for Bin Rate Coefficient:      This%MinForRates  = ", This%MinForRates )
            
            
          case("Reading Formatted Rates?")
            This%ReadFormRatesChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%ReadFormRatesChar))) .eq. 'yes') .or. ((trim(adjustl(This%ReadFormRatesChar))) .eq. 'YES')) This%ReadFormRatesFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Reading Rates?:      This%ReadFormRatesFlg  = ", This%ReadFormRatesFlg )

          case("Reading Unformatted Rates?")
            This%ReadUnformRatesChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%ReadUnformRatesChar))) .eq. 'yes') .or. ((trim(adjustl(This%ReadUnformRatesChar))) .eq. 'YES')) This%ReadUnformRatesFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Reading Rates?:      This%ReadUnformRatesFlg  = ", This%ReadUnformRatesFlg )
            
          case("Reading Unique Rate-File for Multiple Temperatures?")
            This%ReadAllTsChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%ReadAllTsChar))) .eq. 'yes') .or. ((trim(adjustl(This%ReadAllTsChar))) .eq. 'YES')) This%ReadAllTsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Read a Single File for Rates for Multiple Temperatures?:      This%ReadAllTsFlg  = ", This%ReadAllTsFlg )

          case("Path to Rates Folder")
            This%RatesFolderPath = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Path to Rates Folder:      This%RatesFolderPath  = ", This%RatesFolderPath )
            
            

          case("Writing Formatted Merged Rates?")
            This%WriteFormRatesChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteFormRatesChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteFormRatesChar))) .eq. 'YES')) This%WriteFormRatesFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Merged Rates?:      This%WriteFormRatesFlg  = ", This%WriteFormRatesFlg )

          case("Writing Unformatted Merged Rates?")
            This%WriteUnformRatesChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteUnformRatesChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteUnformRatesChar))) .eq. 'YES')) This%WriteUnformRatesFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Merged Rates?:      This%WriteUnformRatesFlg  = ", This%WriteUnformRatesFlg )
            
          case("Writing Unique Rate-File for Multiple Temperatures?")
            This%WriteAllTsChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteAllTsChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteAllTsChar))) .eq. 'YES')) This%WriteAllTsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Write a Single File for Rates for Multiple Temperatures?:      This%WriteAllTsFlg  = ", This%WriteAllTsFlg )

          case("Writing Arrhenius Coefficients?")
            This%WriteArrChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteArrChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteArrChar))) .eq. 'YES')) This%WriteArrFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Arrhenius Coefficients?:      This%WriteArrFlg  = ", This%WriteArrFlg )

          case("Merge Internal Exchanges in Arrhenius Rates?")
            This%MergeIntExchChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%MergeIntExchChar))) .eq. 'yes') .or. ((trim(adjustl(This%MergeIntExchChar))) .eq. 'YES')) This%MergeIntExchFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Merge Internal Exchanges in Arrhenius Rates?:      This%MergeIntExchFlg  = ", This%MergeIntExchFlg )

          case("Nb of Processors for External Code")
            This%NProcExtCodeChar = trim(adjustl(line_input(i_eq+2:150)))
            READ(This%NProcExtCodeChar, '(I6)') This%NProcForExtCode
            if (i_Debug_Loc) call Logger%Write( "Nb Processors for External Code:      This%NProcForExtCode  = ", This%NProcForExtCode )


          case("Computing Dissociation Thermal Rates?")
            This%WriteOverallRatesChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%WriteOverallRatesChar))) .eq. 'yes') .or. ((trim(adjustl(This%WriteOverallRatesChar))) .eq. 'YES')) This%WriteOverallRatesFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Dissociation Thermal Rates?:      This%WriteOverallRatesFlg  = ", This%WriteOverallRatesFlg )
            
          

          case("Binning the StS Rates?")
            This%BinStsChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%BinStsChar))) .eq. 'yes') .or. ((trim(adjustl(This%BinStsChar))) .eq. 'YES')) This%BinStsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Binning the StS Rates?:      This%BinStsFlgf  = ", This%BinStsFlg )

          case("Path to the Folder of the Coarse-Grained System")
            This%SystemCGPath = trim(adjustl(line_input(i_eq+2:150)))
            if (i_Debug_Loc) call Logger%Write( "Path to the Folder of the Coarse-Grained System:      This%SystemCGPath  = ", This%SystemCGPath )
            
            
          case("Writing Exothermic Processes?")
            This%Kinetics_ExoChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_ExoChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_ExoChar))) .eq. 'YES')) This%Kinetics_ExoFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Exothermic Processes?:      This%Kinetics_ExoFlg  = ", This%Kinetics_ExoFlg )

          case("Writing Endothermic Processes?")
            This%Kinetics_EndoChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_EndoChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_EndoChar))) .eq. 'YES')) This%Kinetics_EndoFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Endothermic Processes?:      This%Kinetics_EndoFlg  = ", This%Kinetics_EndoFlg )
            
            
          case("Writing Inelastic Processes?")
            This%Kinetics_InelChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_InelChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_InelChar))) .eq. 'YES')) This%Kinetics_InelFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Inelastic Processes?:      This%Kinetics_InelFlg  = ", This%Kinetics_InelFlg )

          case("Writing Dissociation Processes?")
            This%Kinetics_DissChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_DissChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_DissChar))) .eq. 'YES')) This%Kinetics_DissFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Dissociation Processes?:      This%Kinetics_DissFlg  = ", This%Kinetics_DissFlg )
            
          case("Dissociation Rates Correction Factor")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%DissCorrFactor
            if (i_Debug_Loc) call Logger%Write( "Dissociation Rates Correction Factor:      This%DissCorrFactor  = ", This%DissCorrFactor )            

          case("Writing Exchange Processes?")
            This%Kinetics_ExchChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_ExchChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_ExchChar))) .eq. 'YES')) This%Kinetics_ExchFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing Exchange Processes?:      This%Kinetics_ExchFlg  = ", This%Kinetics_ExchFlg )


          case("Writing CO+C Dissociation?")
            This%Kinetics_COCDissChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_COCDissChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_COCDissChar))) .eq. 'YES')) This%Kinetics_COCDissFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing CO+C Dissociation?:      This%Kinetics_COCDissFlg  = ", This%Kinetics_COCDissFlg )

          case("Writing O2+M Dissociation?")
            This%Kinetics_O2DissChar = trim(adjustl(line_input(i_eq+2:150)))
            if (((trim(adjustl(This%Kinetics_O2DissChar))) .eq. 'yes') .or. ((trim(adjustl(This%Kinetics_O2DissChar))) .eq. 'YES')) This%Kinetics_O2DissFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Writing O2+M Dissociation?:      This%Kinetics_O2DissFlg  = ", This%Kinetics_O2DissFlg )


          case("KONIG Inital Temperature [K]")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TInit
            if (i_Debug_Loc) call Logger%Write( "KONIG Inital Temperature [K]:      This%TInit  = ", This%TInit )
            write(This%TInit_char,"(I10)") floor(This%TInit)
            
    
          case("Path to KONIG Executable")
            This%KONIGRunCMD = trim(adjustl(line_input(i_eq+2:150)))
            if (i_Debug_Loc) call Logger%Write( "Path to KONIG Executable:      This%KONIGRunCMD  = ", This%KONIGRunCMD )

          case("Path to HEGEL Executable")
            This%HEGELRunCMD = trim(adjustl(line_input(i_eq+2:150)))
            if (i_Debug_Loc) call Logger%Write( "Path to HEGEL Executable:      This%HEGELRunCMD  = ", This%HEGELRunCMD )
            

          case("Nb of Chemical Components")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NComp
            if (i_Debug_Loc) call Logger%Write( "Nb of Chemical Components:      This%NComp  = ", This%NComp )
            allocate( This%CompName(This%NComp), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%CompName" )
            allocate( This%CompMolFrac(This%NComp), Stat=Status )
            if (Status/=0) call Error( "Error allocating This%CompMolFrac" )

        end select
        
        
        do iMol = 1,This%NMolecules
          write(iMol_char, "(I2)") iMol

          NBinsCG_case = "Nb of Bins for Binning StS, Molecule" // iMol_char
          if (adjustl(trim(i_case)) == TRIM(NBinsCG_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NBinsCG(iMol)
            WRITE(This%NBinsCG_char(iMol), '(I6)') This%NBinsCG(iMol)
            if (i_Debug_Loc) call Logger%Write( "Nb of Bins for Binning StS, Molecule Nb", iMol, ":      This%NBinsCG(i)  = ", This%NBinsCG(iMol), " = ", This%NBinsCG_char(iMol) )
          end if

        end do
        
        
        do iiPES = 1,This%NPESs
          write(iiPES_char, "(I2)") iiPES

          iiPES_case = "Results File Path, PES" // iiPES_char
          if (adjustl(trim(i_case)) == TRIM(iiPES_case)) then
            This%ResultsFile(iiPES) = line_input(i_eq+2:150)
          end if
          
          iiPESWeight_case = "Weight, PES" // iiPES_char
          if (adjustl(trim(i_case)) == TRIM(iiPESWeight_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%PESWeight(iiPES)
            if (i_Debug_Loc) call Logger%Write( "Weight, PES Nb", iiPES, ":      This%PESWeight(iiPES)  = ", This%PESWeight(iiPES) )
          end if
          
        end do
        
        
        do iComp = 1,This%NComp
          write(iComp_char, "(I3)") iComp

          iCompName_case = "Comp " // adjustl(trim(iComp_char)) // 'Name'
          if (adjustl(trim(i_case)) == TRIM(iCompName_case)) then
            This%CompName(iComp) = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Component ", iComp, " Name:      This%CompName(iComp)  = ", This%CompName(iComp) )
          end if

          iCompMolFrac_case = "Comp " // adjustl(trim(iComp_char)) // 'Initial Mole Fraction'
          if (adjustl(trim(i_case)) == TRIM(iCompMolFrac_case)) then
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)')  This%CompMolFrac(iComp)
            if (i_Debug_Loc) call Logger%Write( "Initial Mole Fraction for Component ", iComp, ":      This%CompMolFrac(iComp)  = ", This%CompMolFrac(iComp) )
          end if

        end do


      end if

    end if


  end do

  close(Unit)
  
  if (.not. allocated(This%PESWeight)) then
    allocate( This%PESWeight(1), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%PESWeight" )
    This%PESWeight = 1
    if (i_Debug_Loc) call Logger%Write( "This%PESWeight Allocated with Lenght 1 and Set to This%PESWeight=1" )
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine FwdProagation( This, i_Debug)

  class(Input_Type)                         ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, Unit
  integer                                                   ::    Status
  character(:)  ,allocatable                                ::    FileName
  integer                                                   ::    NLevels
  character(100)                                            ::    Line

  character(150)                                            ::    i_case
  character(150)                                            ::    line_input
  character(150)                                            ::    line_input_temp
  integer                                                   ::    i_eq

  integer                                                   ::    len
  logical                                                   ::    i_Debug_Loc = .true.


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "FwdProagation" )
  !i_Debug_Loc   =     Logger%On()

  This%NTtra           = 1
  This%NTint           = 1
  This%NSpecies        = 0
  This%NMolecules      = 0
  This%NAtoms          = 0
  This%nring           = 0

!   --------------------------------------------------------------------------------------------------------------
!   --------------------------------------------------------------------------------------------------------------
  FileName = adjustl(trim(This%InputDir)) // '/FwdProagation.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading the Input file for FwdProagation" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )

  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )

!   --------------------------------------------------------------------------------------------------------------
  do
    read(Unit,'(A150)',iostat=Status) line_input
    if (Status /= 0) then
      exit
    else
      line_input_temp=trim(adjustl(line_input))
      if ( (line_input_temp(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
        i_eq = 1
        do
          if (line_input(i_eq:i_eq) == '=') exit
          i_eq = i_eq + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_eq-2))))

        select case (adjustl((TRIM(i_case))))

          case("Nb of Iterations for Forward Progation")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NFwdProp
            if (i_Debug_Loc) call Logger%Write( "Nb of Iterations for Forward Progation:      This%NFwdProp  = ", This%NFwdProp )
            This%NFwdPropProc = floor(real(This%NFwdProp) / real(This%NProc))
            if (i_Debug_Loc) call Logger%Write( "Nb of Fwd Propagation Iterations per Processor:     NFwdProp = ", This%NFwdPropProc )
            
          case("Seed for Forward Progation")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%FwdSeed
            if (i_Debug_Loc) call Logger%Write( "Seed for Forward Progation:      This%FwdSeed  = ", This%FwdSeed )
            This%FwdSeed = int(This%FwdSeed * This%iProc * This%iNode)
            if (i_Debug_Loc) call Logger%Write( "The random Nb seed has been multiplied by the PROCESSOR Nb (", This%iProc ,") and by  the NODE Nb (", This%iNode ,"). Initial random Nb seed now is: This%FwdSeed  = ", This%FwdSeed )
            
          case("KONIG Input File")
            This%KONIGInputFileName = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Run KONIG Command:      This%KONIGInputFileName  = ", This%KONIGInputFileName )

          case("HEGEL Input File")
            This%HEGELInputFileName = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Run HEGEL Command:      This%HEGELInputFileName  = ", This%HEGELInputFileName )
            
          case("Time Min")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%FwdTimeMin
            if (i_Debug_Loc) call Logger%Write( "Time Min:      This%FwdTimeMin  = ", This%FwdTimeMin )

          case("Time Max")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%FwdTimeMax
            if (i_Debug_Loc) call Logger%Write( "Time Max:      This%FwdTimeMax  = ", This%FwdTimeMax )

          case("Nb of Time Nodes")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NTimeNodes
            if (i_Debug_Loc) call Logger%Write( "Nb of Time Steps:      This%NTimeNodes  = ", This%NTimeNodes )

          case("Time Scale")
            This%TimeScale = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Time Scale:      This%TimeScale  = ", This%TimeScale )

          case("Forward Propagating Mole Fractions?")
            This%FWD_MolFractions = line_input(i_eq+2:150)
            
          case("Mole Fractions Min")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MolFractionsMin
            if (i_Debug_Loc) call Logger%Write( "Mole Fractions Min:      This%MolFractionsMin  = ", This%MolFractionsMin )

          case("Mole Fractions Max")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MolFractionsMax
            if (i_Debug_Loc) call Logger%Write( "Mole Fractions Max:      This%MolFractionsMax  = ", This%MolFractionsMax )

          case("Nb of Mole Fractions Bins")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NMolFractionsBins
            if (i_Debug_Loc) call Logger%Write( "Nb of Mole Fractions Bins:      This%NMolFractionsBins  = ", This%NMolFractionsBins )
            
          case("Forward Propagating Temperatures?")
            This%FWD_Temperatures = line_input(i_eq+2:150)

          case("Temperatures Min")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TemperaturesMin
            if (i_Debug_Loc) call Logger%Write( "Temperatures Min:      This%TemperaturesMin  = ", This%TemperaturesMin )

          case("Temperatures Max")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%TemperaturesMax
            if (i_Debug_Loc) call Logger%Write( "Temperatures Max:      This%TemperaturesMax  = ", This%TemperaturesMax )

          case("Nb of Temperatures Bins")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NTemperaturesBins
            if (i_Debug_Loc) call Logger%Write( "Nb of Temperatures Bins:      This%NTemperaturesBins  = ", This%NTemperaturesBins )

          case("Forward Propagating Bins Populations?")
            This%FWD_Populations = line_input(i_eq+2:150)

          case("Bins Populations Min")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%PopulationsMin
            if (i_Debug_Loc) call Logger%Write( "Temperatures Min:      This%PopulationsMin  = ", This%PopulationsMin )

          case("Bins Populations Max")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%PopulationsMax
            if (i_Debug_Loc) call Logger%Write( "Temperatures Max:      This%PopulationsMax  = ", This%PopulationsMax )

          case("Nb of Bins Populations Bins")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NPopulationsBins
            if (i_Debug_Loc) call Logger%Write( "Nb of Populations Bins:      This%NPopulationsBins  = ", This%NPopulationsBins )
        
        end select

      end if

    end if


  end do

  close(Unit)


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine PlotPES( This, i_Debug)

  class(Input_Type)                         ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, Unit
  integer                                                   ::    Status
  character(:)  ,allocatable                                ::    FileName
  integer                                                   ::    NLevels
  character(100)                                            ::    Line

  character(150)                                            ::    i_case
  character(150)                                            ::    line_input
  character(150)                                            ::    line_input_temp
  integer                                                   ::    i_eq
  
  integer                                                   ::    iPar
  character(2)                                              ::    iParChar
  character(20)                                             ::    iPar_case
  
  integer                                                   ::    iA
  character(150)                                            ::    iAngles_case
  character(3)                                              ::    iA_char
  integer                                                   ::    tempInt
  
  integer                                                   ::    len
  logical                                                   ::    i_Debug_Loc = .true.


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES" )
  !i_Debug_Loc   =     Logger%On()


  This%NAnglesPlot = 1

!   --------------------------------------------------------------------------------------------------------------
!   --------------------------------------------------------------------------------------------------------------
  FileName = adjustl(trim(This%InputDir)) // '/PlotPES.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading the Input file for FwdProagation" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )

  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )

!   --------------------------------------------------------------------------------------------------------------
  do
    read(Unit,'(A150)',iostat=Status) line_input
    if (Status /= 0) then
      exit
    else
      line_input_temp=trim(adjustl(line_input))
      if ( (line_input_temp(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
        i_eq = 1
        do
          if (line_input(i_eq:i_eq) == '=') exit
          i_eq = i_eq + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_eq-2))))

        select case (adjustl((TRIM(i_case))))
        
          case("Reading Grid Points?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_ReadPntsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Reading Grid Points?:      This%PlotPES_ReadPntsFlg  = ", This%PlotPES_ReadPntsFlg )
            
          case("Computing Grid Points?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_GridFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Grid Points?:      This%PlotPES_GridFlg  = ", This%PlotPES_GridFlg )
            
          case("Computing Double-Grid Points?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_DoubleGridFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Double-Grid Points?:      This%PlotPES_DoubleGridFlg  = ", This%PlotPES_DoubleGridFlg )
          
          case("Computing Triple-Grid Points?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_TripleGridFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Triple-Grid Points?:      This%PlotPES_TripleGridFlg  = ", This%PlotPES_TripleGridFlg )
          
          case("Computing Points for Scatter Comparison?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_GridForScatterFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Points for Scatter Comparison?:      This%PlotPES_GridForScatterFlg  = ", This%PlotPES_GridForScatterFlg )
            
          case("Computing Stochastic PES Statistics?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_StatsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Stochastic PES Statistics?:      This%PlotPES_StatsFlg  = ", This%PlotPES_StatsFlg )
          
          case("Writing Sampled Parameters of Stochastic PES?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_WriteParamsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Stochastic PES Statistics?:      This%PlotPES_WriteParamsFlg  = ", This%PlotPES_WriteParamsFlg )
            
          case("Plotting Only Triatomic Component of PES?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_OnlyTriatFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Stochastic PES Statistics?:      This%PlotPES_OnlyTriatFlg  = ", This%PlotPES_OnlyTriatFlg )
          
          case("Plotting Cuts with Varga et al.'s style?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_VargasPaperFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Plotting Cuts with Varga et al.'s style?:      This%PlotPES_VargasPaperFlg  = ", This%PlotPES_VargasPaperFlg )

          case("Computing Rotating 3rd Atom?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_Rot3rdFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing Rotating 3rd Atom?:      This%PlotPES_Rot3rdFlg  = ", This%PlotPES_Rot3rdFlg )
            

         case("Computing isosceles triangle?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PlotPES_IsoTriFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Computing isosceles triangle?:      This%PlotPES_IsoTriFlg  = ", This%PlotPES_IsoTriFlg )

        
        
          case("Potential or Force")
            This%POTorFR = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Potential or Force:      This%POTorFR  = ", This%POTorFR)

          case("PES or Diatomic")
            This%PESOrDiatFlg = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "PES or Diatomic:      This%PESOrDiatFlg  = ", This%PESOrDiatFlg )

          case("Y Axis Variable")
            This%YAxisVar = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Y Axis Variable:      This%YAxisVar  = ", This%YAxisVar )

          case("Unit of measure for Force")
            This%UnitFR = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Unit of measure for Force:      This%UnitFR  = ", This%UnitFR )

          case("Unit of measure for Distance")
            This%UnitDist = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Unit of measure for Distance:      This%UnitDist  = ", This%UnitDist )
            
          case("Unit of measure for Potential")
            This%UnitPot = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Unit of measure for Potential:      This%UnitPot  = ", This%UnitPot )
          
          
          case("Random Points?")
            This%RandomPointsChar = line_input(i_eq+2:150)
            if (i_Debug_Loc) call Logger%Write( "Random Points?:      This%RandomPointsChar = ", This%RandomPointsChar )
            if (trim(adjustl(This%RandomPointsChar)) == 'yes') This%RandomPointsFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Random Points?:      This%RandomPointsFlg  = ", This%RandomPointsFlg )
            
          case("Points Distribution")
            This%PointsDistr = trim(adjustl(line_input(i_eq+2:150)))
            if (i_Debug_Loc) call Logger%Write( "Points Distribution:      This%PointsDistr  = ", This%PointsDistr )
            
            
            
          case("Nb of Grid Points")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NGridPlot(1)
            if (i_Debug_Loc) call Logger%Write( "Nb of Grid Points:      This%NGridPlot(1)  = ", This%NGridPlot(1) )

          case("Nb of Grid Points, Axis 1")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NGridPlot(1)
            if (i_Debug_Loc) call Logger%Write( "Nb of Grid Points, Pair 1:      This%NGridPlot(1)  = ", This%NGridPlot(1) )

          case("Nb of Grid Points, Axis 2")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NGridPlot(2)
            if (i_Debug_Loc) call Logger%Write( "Nb of Grid Points, Pair 2:      This%NGridPlot(2)  = ", This%NGridPlot(2) )


          case("Grid Min, Axis 1")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MinGridPlot(1)
            if (i_Debug_Loc) call Logger%Write( "Grid Min, Pair 1:      This%MinGridPlot(1)  = ", This%MinGridPlot(1) )

          case("Grid Min, Axis 2")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MinGridPlot(2)
            if (i_Debug_Loc) call Logger%Write( "Grid Min, Pair 2:      This%MinGridPlot(2)  = ", This%MinGridPlot(2) )


          case("Grid Max, Axis 1")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MaxGridPlot(1)
            if (i_Debug_Loc) call Logger%Write( "Grid Max, Pair 1:      This%MaxGridPlot(1)  = ", This%MaxGridPlot(1) )

          case("Grid Max, Axis 2")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%MaxGridPlot(2)
            if (i_Debug_Loc) call Logger%Write( "Grid Max, Pair 2:      This%MaxGridPlot(2)  = ", This%MaxGridPlot(2) )


          case("Nb of Angles")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(I10)') This%NAnglesPlot
            if (i_Debug_Loc) call Logger%Write( "Nb of Angles:      This%NAnglesPlot  = ", This%NAnglesPlot )
            allocate( This%AnglesPlot(This%NAnglesPlot), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%Angles" )
            allocate( This%AnglesPlotChar(This%NAnglesPlot), Stat=Status)
            if (Status/=0) call Error( "Error allocating This%AnglesChar" )

          
            
          case("Energy Cut-Off (From Dissociation Energy)")
            line_input = line_input(i_eq+2:150)
            READ(line_input, '(d20.10)') This%EnergyCutOff
            if (i_Debug_Loc) call Logger%Write( "Energy Cut-Off (From Dissociation Energy):      This%EnergyCutOff  = ", This%EnergyCutOff )
            
          case("Setting the Energy Reference to 0.0?")
            line_input = line_input(i_eq+2:150)
            if ((trim(adjustl(line_input)) == 'yes') .or. (trim(adjustl(line_input)) == 'YES')) This%PESZeroRefFlg = .true.
            if (i_Debug_Loc) call Logger%Write( "Setting the Energy Reference to 0.0?:      This%PESZeroRefFlg  = ", This%PESZeroRefFlg )
        

        end select
        
        
        do iPar = 1,5
          write(iParChar, "(I2)") iPar
          
          iPar_case = "Parameter " // adjustl(trim(iParChar))
          if (adjustl(trim(i_case)) == TRIM(iPar_case)) then
            This%DistrParChar(iPar) = adjustl(trim(line_input(i_eq+2:150)))
            READ(This%DistrParChar(iPar), '(d20.10)')  This%DistrPar(iPar)
            if (i_Debug_Loc) call Logger%Write( "Parameter ", iPar, ":      This%DistrParChar(iPar)  = ", This%DistrParChar(iPar), " = This%DistrPar(iPar)  = ", This%DistrPar(iPar) )
          end if
          
        end do
        
        
        do iA = 1,This%NAnglesPlot
          write(iA_char, "(I3)") iA
          
          iAngles_case = "Angle " // adjustl(trim(iA_char))
          if (adjustl(trim(i_case)) == TRIM(iAngles_case)) then
            This%AnglesPlotChar(iA) = adjustl(trim(line_input(i_eq+2:150)))
            READ(This%AnglesPlotChar(iA), '(d20.10)')  This%AnglesPlot(iA)
            tempInt = int(This%AnglesPlot(iA))
            WRITE(This%AnglesPlotChar(iA), '(I10)') tempInt
            if (i_Debug_Loc) call Logger%Write( "Angle ", iA, ":      This%AnglesPlotChar(iA)  = ", This%AnglesPlotChar(iA), " = This%AnglesPlot(iA)  = ", This%AnglesPlot(iA) )
          end if
          
        end do
                  

      end if

    end if


  end do

  close(Unit)


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine ReadAtomMass(This, iAtom, i_Debug)

  class(Input_Type)                         ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iAtom
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    i, Unit
  integer                                                   ::    Status
  character(:)  ,allocatable                                ::    FileName
  character(20)                                             ::    AtomMassChar
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadAtomMass" )
  !i_Debug_Loc   =     Logger%On()

  !   --------------------------------------------------------------------------------------------------------------
  !   --------------------------------------------------------------------------------------------------------------
  if ( adjustl(trim(This%FileForAtomMass(iAtom))) == 'NONE') then
    FileName = adjustl(trim(This%DtbPath)) // '/Masses/' // adjustl(trim(This%AtomsName(iAtom))) // '.dat'
  elseif ( adjustl(trim(This%FileForAtomMass(iAtom))) == 'Local') then
    FileName = trim(adjustl(This%OutputDir)) // '/' // trim(adjustl(This%System)) // '/Masses/' // adjustl(trim(This%AtomsName(iAtom))) // '.dat'
  else
    FileName = adjustl(trim(This%FileForAtomMass(iAtom))) // '/' // adjustl(trim(This%AtomsName(iAtom))) // '.dat'
  end if
  
  if (i_Debug_Loc) call Logger%Write( "Reading the Mass for " // adjustl(trim(This%AtomsName(iAtom))) // " Atom" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )

  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )

    read(Unit,*)
    read(Unit,'(A20)') AtomMassChar
    READ(AtomMassChar, '(d20.10)')  This%AtomsMass(iAtom)
    if (i_Debug_Loc) call Logger%Write( "Found Mass: Mass = " // adjustl(trim(AtomMassChar)) // " [a.u.]" )

  close(Unit)

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


End Module