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

! @BUG: Problem with the impact parameter b when impose ang. molemtum: The value should be read from the inputi
! @BUG: Variable Ein is not defiend
! @BUG: In the original code, for an impose total angular momentum, the impact parameter was not set ('CoordVeloc_ImpTotAngMom')
! @TODO: The traj object should be initialized fron the Collision object to pass the variables (param) which are in common between Collision?TRaj
! This was initially inicong.f
Module Collision_Class

#include "../qct.inc"

  use Parameters_Module           ,only:  rkp, Zero, Half, One, Two, Pi, Kelvin_To_Hartree, Rugc
  use Logger_Class                ,only:  Logger, LogLevel_INFO
  use Error_Class                 ,only:  Error

  use Transformation_Class        ,only:  Transformation_Type
  use Atom_Class                  ,only:  Atom_Type
  use AtomsPair_Class             ,only:  AtomsPair_Type
  use Species_Class               ,only:  Species_Type
  use PESsContainer_Class         ,only:  PESsContainer_Type
  use PES_Class                   ,only:  PES_Type, DiatPotContainer_Type, PESEvoFile, PESEvoFlg
  use ImpactParameter_Class       ,only:  ImpactParameter_Type
  use StateInitDiatomAtom_Module  ,only:  Params
  use MoleculesContainer_Class    ,only:  MoleculesContainer_Type
  use Output_Module               ,only:  OutputFile_Type
  use Timing_Module

  implicit none

  private
  public    ::    Collision_Type


  Type      ::    Collision_Type
    logical                                               ::    Initialized   = .False.
    integer                                               ::    NTrajTot      =   0
    integer                                               ::    NAtoms        =   0          ! Number of atoms involved in the collision
    integer                                               ::    NPairs        =   0          ! Number of pairs of atoms (Old name was 'maxpr'): NAtoms*(NAtoms-1)/2
    integer                                               ::    NPESs         =   0          ! Number of PESs
    integer                                               ::    NSpecies      =   0          ! Number of species involved in the collision: Always equal to 2
    integer                                               ::    NAtoMaxSpe    =   0          ! Maximum number of atoms per species: Always equal to 2
    integer                                               ::    NEqtVar       =   0          ! Number of equations per variable (Coordinates/Momenta): 3*(NAtoms-1)
    integer                                               ::    NEqtTot       =   0          ! Total number of equations (for both Coordinates/Momenta): 2 * NEqtVar
    integer                                               ::    iSeed
    integer                    ,dimension(:) ,allocatable ::    iSeedVec
    integer                    ,dimension(:) ,allocatable ::    jSeedVec
    logical                                               ::    ImposedTotalAngularMomentum  ! Indicator of an imposed total angular momentum. Old name: ijtot
    integer                                               ::    icoord        =   0          ! Indicator of the coordinate system: jacobi coordinates => True, cartesiancoordinates => False.
                                                                                             !   Old name: icoord=0 => cartesian coordinates, icoord=1 => jacobi coordinates
    character(:)                            ,allocatable  ::    TtraModel                    ! Model for Translational Energy Description
    real(rkp)                                             ::    Ttra                         ! Translational temperature [K] or initial relative translational energy [hartree]
    real(rkp)                                             ::    Etot                         ! Total energy for trajectories
    real(rkp)                                             ::    Erel                         ! Fixed relative translational kinetic energy [hartree]
    real(rkp)                                             ::    EMu                          ! Mean Energy
    real(rkp)                                             ::    ESD                          ! Energy Standard Deviation
    real(rkp)                                             ::    TotalAngularMomentum         ! Imposed total angular momentum. Read from input only if 'bstart=2345432'. Old name: xjtot
    real(rkp)                                             ::    Dinit                        ! Initial separation of the the target and projectile
    real(rkp)                                             ::    RedMass                      ! Reduced mass of the target/projectile
    type(Transformation_Type)                             ::    Transformation

    integer   ,dimension(:)   ,allocatable                ::    NAtomsPerSpecies
    integer   ,dimension(:,:) ,allocatable                ::    AtomsToSpecies                  ! Index mapping from atoms to species. Dim=(NAtomsMaxPerSpecies=2,NPairs)

    type(AtomsPair_Type)      ,dimension(:) ,allocatable  ::    Pairs                        ! List of Pairs object:
    integer                   ,dimension(:)  ,allocatable ::    OppositePair
    integer                   ,dimension(:,:),allocatable ::    Atoms_To_Pair                ! Index mapping from a set of 2 atoms to the associated pairs. Dim=(NAtoms,NAtoms) @TODO
    integer   ,dimension(:,:) ,allocatable                ::    Pair_To_Atoms                    ! Index mapping from atoms to pairs.   Dim=(NAtomsPerPair=2,NPairs)

    type(Atom_Type)           ,dimension(:) ,allocatable  ::    Atoms                        ! List of Atoms object: Atoms of the target species are listed first, 
                                                                                             !   then the atoms of the projectile species. Dim=(NAtoms)
    type(Species_Type)        ,dimension(:) ,allocatable  ::    Species
    class(PESsContainer_Type) ,dimension(:) ,allocatable  ::    PESsContainer
    type(OutputFile_Type)                                 ::    OutputFile
    type(ImpactParameter_Type),dimension(:) ,allocatable  ::    ImpactPara
    
    integer                   ,dimension(:) ,allocatable  ::    DiffPairs
    !integer                   ,dimension(0:51)            ::    Arrangements
    integer                   ,dimension(0:101)           ::    Arrangements
!   ********************************************
    real(rkp)                 ,dimension(:) ,allocatable  ::    mMiMn                         !< Opposite of mass ratio: -Mi/Mn, with Mn mass of the last atom. Dim=(NAtoms-1). Old name xmss
    logical                                               ::    NormalKineticEnergy = .False. ! In the old code, this indicator used to be true if 'tpq(1,1) = 1.02030201d0'
!   ********************************************
    real(rkp)                                             ::    EqVelocity
    class(MoleculesContainer_Type) ,dimension(:) ,allocatable  ::    MoleculesContainer

    integer                                               ::    NTrajOverall
    integer                                               ::    iProc

    procedure(SetterInitialState)       ,pointer  ,nopass ::    SetInitialState
    procedure(ComputePairDistance)      ,pointer  ,nopass ::    Compute_Rp_CartCoord

  contains
    private
    procedure ,public   ::    Initialize => InitializeCollision
    procedure ,public   ::    InitializeAtomsPairsSpecies
    procedure ,public   ::    InitializeMolecules
    procedure ,public   ::    InitializeTrajectories
    procedure ,public   ::    SetInitialConditions
    procedure ,public   ::    ApplyTransformation
    procedure ,public   ::    ApplyAntiTransformation

    procedure ,public   ::    Compute_PES =>  Compute_PES_1d

    generic   ,public   ::    Hamiltonian =>  Hamiltonian_0d, Hamiltonian_1d
    procedure           ::    Hamiltonian_0d
    procedure           ::    Hamiltonian_1d

    generic   ,public   ::    Compute_Rp  =>  Compute_Rp_0d, Compute_Rp_1d
    procedure           ::    Compute_Rp_0d
    procedure           ::    Compute_Rp_1d

    procedure           ::    InitializeCollisionOutput
    procedure           ::    InitializeImpactParameter
    procedure           ::    SetIniCondForward
    
    procedure           ::    ComputeCoordinatesVelocities                      !< Computes the cartesian coordinates of the target and projectile wrt their centers of mass.
    
    procedure           ::    SetPESIndx
    procedure           ::    ShiftCoordinates
    procedure           ::    SetPaQ
    procedure   ,nopass ::    SetErrorFactor
    procedure           ::    inia2
    procedure           ::    SetInitialPaQ
  End Type

  Abstract Interface
    Subroutine SetterInitialState( Species, Qijk, dQijk, iTraj, iPES, NTrajOverall, iProc, i_Debug, i_Debug_Deep )
      use Species_Class         ,only:  Species_Type
      use Parameters_Module     ,only:  rkp
      type(Species_Type)  ,dimension(:)     ,intent(in)   ::    Species         !< Species objects which always have 2 elements [Target,Projectile]. Dim=(NSpecies)=(2)
      real(rkp)           ,dimension(:,:,:) ,intent(out)  ::    Qijk            !< Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). 
                                                                                !     Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
      real(rkp)           ,dimension(:,:,:) ,intent(out)  ::    dQijk           !< Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile).                   !     Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
      integer                               ,intent(in)   ::    iTraj
      integer                               ,intent(in)   ::    iPES
      integer                               ,intent(in)   ::    NTrajOverall
      integer                               ,intent(in)   ::    iProc
      logical                     ,optional ,intent(in)   ::    i_Debug
      logical                     ,optional ,intent(in)   ::    i_Debug_Deep
    End Subroutine

    PURITY Subroutine ComputePairDistance( mMiMn, x, Rxyz, Rp )
      use Parameters_Module     ,only:  rkp
      real(rkp) ,dimension(:)    CONTIGUOUS ,intent(in)   ::    mMiMn           !< Opposite of mass ratio: -Mi/Mn, with Mn mass of the last atom. Dim=(NAtoms-1). Old name xmss
      real(rkp) ,dimension(:,:)  CONTIGUOUS ,intent(in)   ::    x               !< Position. Dim=(NEqtVar,NTraj)
      real(rkp) ,dimension(:,:)  CONTIGUOUS ,intent(out)  ::    Rxyz            !< Cartesian coordinates of last atom. Dim=(NSpace,NTraj)
      real(rkp) ,dimension(:,:)  CONTIGUOUS ,intent(out)  ::    Rp              !< Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
    End Subroutine

  End Interface

  integer   ,parameter    ::    NSpace         = 3
  logical   ,parameter    ::    Formatted      = .True.
  integer                 ::    iSpeTar        = 1
  integer                 ::    iSpePro        = 2
  real(rkp)               ::    MostProbEr 
  
  logical   ,parameter    ::    i_Debug_Global = .False.
  logical   ,parameter    ::    i_Debug_Deep   = .False.
  
  contains

!________________________________________________________________________________________________________________________________!
Subroutine InitializeCollision( This, Input, i_Debug, i_Debug_Deep )
! This procedure initializes the collision parameters. The input parameters are loaed into the Collision object. These parameters are then kept constant for all trajectories. This procedure is independent on the number of trajectories.
  
  use Input_Class                 ,only:  Input_Type
  use Parameters_Module           ,only:  rkp, Zero, One, Two, Pi, Kelvin_To_Hartree, Rugc
  use RandomVector_Module         ,only:  RandSetVec
  use RandomVectorBis_Module      ,only:  RandSetVecBis
  use PES_Factory_Class           ,only:  PES_Factory_Type
  use Global_Module               ,only:  UnitTraj
  use StateInitDiatomAtom_Module  ,only:  SetInitialState_DiatomAtom,   InitializeLevels_DiatomAtom
  use StateInitDiatomDiatom_Module,only:  SetInitialState_DiatomDiatom, InitializeLevels_DiatomDiatom


  class(Collision_Type)                     ,intent(out)    ::    This
  type(Input_Type)                          ,intent(inout)  ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  integer                                                   ::    iA                              ! Index of atoms
  integer                                                   ::    iS                              ! Index of species
  integer                                                   ::    iP                              ! Index of pairs of atoms
  integer                                                   ::    i, j,n
  integer                                                   ::    iPES, iPESSample
  integer                                                   ::    iMol
  integer                                                   ::    Unit
  integer                                                   ::    Status
  integer                                                   ::    nsum
  integer                                                   ::    Opposite
  real(rkp)                                                 ::    Vp                              ! Most probable speed: the speed most likely to be possessed by species
  real(rkp)                                                 ::    Va                              ! Average or mean speed: expected value of the speed distribution
  real(rkp)                                                 ::    Vrms                            ! The root mean square speed is the second-order moment of speed:
  real(rkp)                                                 ::    xmu                             ! Reduced mass for relative translation
  type(PES_Factory_Type)                                    ::    PES_Factory
  character(:)              ,allocatable                    ::    FileName
  integer   ,dimension(64)                                  ::    seed
  integer                                                   ::    iSeed
  integer                                                   ::    NTrajTemp
  integer                                                   ::    NTrajPES
  real(rkp) ,dimension(2)                                   ::    array
  integer                                                   ::    TaskType_Loc
  logical                                                   ::    i_Debug_Loc
  
  !!! Input%TaskType = 1 in Plotting PESs
  !!! Input%TaskType = 2 in Computing Levels
  !!! Input%TaskType = 3 in Preprocessing Levels
  !!! Input%TaskType = 4 in Running Trajectories
  !!! Input%TaskType = 5 in Computing Statistics
  !!! Input%TaskType = 6 in Postprocessing Trajectories

  TaskType_Loc = 4; if ( Input%TaskType > 0 ) TaskType_Loc = Input%TaskType

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeCollision" )
  !i_Debug_Loc   =     Logger%On()


  ! ==============================================================================================================
  ! The type of coordinate system is defined by the variable This%icoord which can take the following values:
  !  * icoord = 0:  Cartesian coordinates are used
  !  * icoord = 1:  Cartesian jacobi coordinates are used in which the target coordinates will be assumed to
  !                 be the vector from the center-of-mass of the target and the center-of-mass of the projectile,
  !                 and the coordinates supplied by the 'SetInitialState' procedure will be overwriten by the
  !                 vector between the "pseudo"-atoms.
  !  * icoord = 2:  Cartesian jacobi coordinates are used in which "target" and "projectile" are exchanged compared to icoord = 1
  ! Note that when jacobi coordinates are used, the number of atoms read in will be one more than the
  ! number of jacobi coordinates and the last mass will not be used.
  ! Note also that if This%icoord=0, NAtomsTarg+NAtomsProj must equal the total number of atoms NAtoms,
  ! while if This%icoord > 0, NAtomsTarg+NAtomsProj must equal the number of coordinate vectors, which is NAtoms-1.
  ! Also note that if This%icoord=1, NAtomsTarg=1 also and if This%icoord=2, NAtomsProj=1 also.
  ! ==============================================================================================================


  !  1. INITIALIZING RANDOM GENERATOR: Initializing a Vector of Random Gen.s of dimension equal to the number of different PESs to run
  !                                    Each PES, on each Processor, on each Node, will have a different seed. 
  ! 
  !  2. IMPORTING LOCAL VARIABLES FROM INPUT OBJECT
  ! 
  !  3. INITIALIZING ATOM, PAIR, AND SPECIES OBJECTS
  !
  !  4. INITIALIZING PES OBJECTS
  !
  !  5. COMPUTING TRANFORMATION MATRIX
  !
  !  6. COMPUTING NORMALIZED and REDUCED MASSES
  !
  !  7. INITIALING MOLECULES (WITH LEVELS AND BINS)
  !
  !  8. SETTING PARAMETERS SPECIFYING THE RELATIVE KINETIC ENERGY AND INITIAL IMPACT PARAMETER
  !
  !  9. INITIALIZING MOLECULE STATES
  !
  ! 10. SETTING THE AVERAGE VELOCITY
  !
  ! 11. FIND MOST PROBABLE TRANSLATIONAL ENERGY AT THE CURRENT ENERGY
  ! 
  ! 12. INITIALIZING OUTPUT FILES


  if (TaskType_Loc == 4) then   !!! <----- ONLY if RUNNING TRAJECTORIES
    ! ==============================================================================================================
    !  1. INITIALIZING RANDOM GENERATOR
    ! ==============================================================================================================
      allocate( This%iSeedVec(Input%NPESs+1), stat=Status )
      if (Status/=0) call Error( "Error allocating This%iSeedVec" )
      if (i_Debug_Loc) call Logger%Write( "Allocated This%iSeedVec with dimension Input%NPES = ", Input%NPESs, " + 1" )
      This%iSeedVec = 1
      
      if (Input%PESRndIniCondFlg) then
        This%iSeedVec = (/ (i, i = 1,Input%NPESs) /)
      end if
      
      This%iSeedVec = int( Input%iSeed * (This%iSeedVec + (Input%iProc-1)*(Input%NPESs+1) + (Input%iNode-1)*200.d0*(Input%NPESs+1)) + &
                           (Input%RndNb)*200.d0*100.d0*(Input%NPESs+2) )
      if (i_Debug_Loc) call Logger%Write( "Initialized This%iSeedVec. First value = ", This%iSeedVec(1) )
      
      call RandSetVec( This%iSeedVec, .False. )

      allocate( This%jSeedVec(Input%NPESs+1), stat=Status )
      if (Status/=0) call Error( "Error allocating This%jSeedVec" )
      if (i_Debug_Loc) call Logger%Write( "Allocated This%jSeedVec with dimension Input%NPES = ", Input%NPESs, " + 1" )
      This%jSeedVec = 1
      
      if (Input%PESRndIniCondFlg) then
        This%jSeedVec = (/ (i, i = 1,Input%NPESs) /)
      end if
      
      This%jSeedVec = int( Input%iSeed * (This%jSeedVec + (Input%iProc-1)*(Input%NPESs+1) + (Input%iNode-1)*200.d0*(Input%NPESs+1)) + &
                           (Input%RndNb)*200.d0*100.d0*(Input%NPESs+2) )
      if (i_Debug_Loc) call Logger%Write( "Initialized This%jSeedVec. First value = ", This%jSeedVec(1) ) 
      call RandSetVecBis( This%jSeedVec, .False. )
    ! ==============================================================================================================  
  end if


  ! ==============================================================================================================
  !  2. SETTING DIMENSIONS: This%NTrajTot, This%NAtoms, This%NEqtTot, This%NPairs
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting dimensions" )
  This%Initialized  =   .True.
  This%NTrajTot     =   Input%NTrajBatch                                                          ! Getting the total number of trajectories
  This%NPESs        =   Input%NPESs                                                               ! Nb of PESs
  This%Dinit        =   Input%Dinit                                                               ! Setting the initial separation of the the target and projectile  
  if (i_Debug_Loc) then
    call Logger%Write( "-> Total number of trajectories:                        This%NTrajTot = ", This%NTrajTot )
    call Logger%Write( "-> Total number of PESs:                                This%NPESs    = ", This%NPESs )
    call Logger%Write( "-> Initial separation of the the target and projectile: This%Dinit    = ", This%Dinit )
  end if
  ! ==============================================================================================================


  ! if (TaskType_Loc == 4) then
  !   !==============================================================================================================
  !   !    consider back integrating
  !   !==============================================================================================================
  !     if ( Input%imode < 0 ) then
  !       if (i_Debug_Loc) call Logger%Write( "Backward integration" )
  !       imode1      =   abs(Input%imode)
  !       if (i_Debug_Loc) call Logger%Write( "-> Every",imode1,"-th trajectory" )
  !       if (i_Debug_Loc) call Logger%Write( "-> Opening unit 8" )
  !       read(5,*) LongString
  !       FileName  =   LongString
  !       open( Unit=UnitTraj, file=FileName, Form='unformatted', Access='direct', Status='old', Recl=rkp*(This%NEqtTot+2) )
  !       allocate( Rvec(This%NEqtTot+2) )
  !       read(UnitTraj,rec=1) Rvec
  !       temp1         =   Rvec(1)
  !       temp2         =   Rvec(2)
  !       This%NTrajTot      =   int(temp2+0.001d0) / imode1
  !       This%NTrajTot      =   max(This%NTrajTot,1)
  !     else
  !       if (i_Debug_Loc) call Logger%Write( "Opening file unit: UnitTraj = ", UnitTraj )
  !       allocate( FileName , source = 'fort.8' )
  !       open( Unit=UnitTraj, Form='unformatted', Access='direct', Status='unknown', Recl=rkp*(This%NEqtTot+2), File=FileName, Iostat=Status )
  !       if (Status/=0) call Error( "Error opening file " // FileName )
  !     end if
  !   !==============================================================================================================
  ! end if
  

  ! ==============================================================================================================
  !   3. INITIALING PAIRS
  ! ==============================================================================================================  
  if (i_Debug_Loc) call Logger%Write( "Entering This%InitializeAtomsPairsSpecies" )
  call This%InitializeAtomsPairsSpecies( Input, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) call Logger%Write( "Done with This%InitializeAtomsPairsSpecies" )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   4. SETTING THE PES (POTENTIAL ENERGY SURFACE) OBJECT (Allocating the abstract PES to the system PES)
  ! ==============================================================================================================  
  if (Input%StochPESFlg) then
    if (i_Debug_Loc) call Logger%Write( "Selected a Stochastic PES" )
    
    allocate(This%PESsContainer(Input%NPESs), stat=Status)
    if (Status/=0) call Error( "Error allocating This%PESsContainer" )
    if (i_Debug_Loc) call Logger%Write( "Allocated PESsContainer" )
      
    if (i_Debug_Loc) call Logger%Write( "Setting the potential energy surface object. Calling PES_Factory%Construct_PES: System = ", Input%System )
    
    if (Input%NTrajBatch < Input%NPESs) call Error( "Nb of Trajectories < Nb of PES Samples! Please, change Input File!" )
    do iPESSample = 1,Input%NPESs
      !associate( PES => This%PESsContainer(iPES)%PES )
        call PES_Factory%Construct_PES( Input, This%Atoms, iPESSample, This%PESsContainer(iPESSample)%PES, i_Debug=i_Debug_Loc )
        seed  = Input%PESiseed + iPESSample
        iSeed = Input%PESiseed + iPESSample
        call random_seed(put = seed)
        if (Input%SampleParamsStochPES) then
          call This%PESsContainer(iPESSample)%PES%SampleParamPost()
        else
          call This%PESsContainer(iPESSample)%PES%ReadParamPost(iSeed)
        end if 
        
        if (TaskType_Loc == 1) then
          if (Input%PlotPES_WriteParamsFlg) then
            call This%PESsContainer(iPESSample)%PES%WriteParamSample(iSeed, Input%LevelOutputDir)
          end if
        end if
      !end associate
    end do
    
  else
    if (i_Debug_Loc) call Logger%Write( "Selected a Deterministic PES" )
      
    allocate(This%PESsContainer(Input%NPESs),stat=Status)
    if (Status/=0) call Error( "Error allocating This%PESsContainer" )
    if (i_Debug_Loc) call Logger%Write( "Allocated PESsContainer" )
      
    if (i_Debug_Loc) call Logger%Write( "Setting the potential energy surface object. Calling PES_Factory%Construct_PES: System = ", Input%System )
    do iPES = 1,Input%NPESs
      !associate( PES => This%PESsContainer(iPES)%PES )
        call PES_Factory%Construct_PES( Input, This%Atoms, iPES, This%PESsContainer(iPES)%PES, i_Debug=i_Debug_Loc )
      !end associate
    end do
    if (i_Debug_Loc) call Logger%Write( "-> Done constructing the PES object" )
    
  end if    
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the PES object" )
  ! ==============================================================================================================


  if (TaskType_Loc == 4) then   !!! <----- ONLY if RUNNING TRAJECTORIES
    ! ==============================================================================================================
    !   5. CONSTRUCTING TRANFORMATION MATRIX FROM/TO CARTESIAN TO/FROM GENERALIZED COORDINATES
    ! ==============================================================================================================
      if (i_Debug_Loc) call Logger%Write( "Constructing tranformation matrix" )
      if (i_Debug_Loc) call Logger%Write( "-> Calling Transformation%Initialize" )
      call This%Transformation%Initialize( This%Atoms(:)%Mass, i_Debug=i_Debug_Deep )
      if (i_Debug_Loc) call Logger%Write( "-> Done with Transformation%Initialize" )
    ! ==============================================================================================================
  end if


  ! ==============================================================================================================
  !   6. SETTING REDUCED AND OPPOSITE MASSES
  ! ==============================================================================================================
    This%NEqtVar      =   3 * (This%NAtoms-1)                                                       ! Setting the number of equations per variable (Coordinates/momenta): 3*(NAtoms-1)
    This%NEqtTot      =   2 * This%NEqtVar                                                          ! Setting the total number of equations (for both Coordinates/Momenta): 2*NEqtVar

    allocate( This%mMiMn(This%NAtoms-1), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%mMiMn" )
    This%mMiMn(1:This%NAtoms-1) = - Input%AtomsMass(1:This%NAtoms-1) / Input%AtomsMass(This%NAtoms)
    do iPES = 1,Input%NPESs
      associate( PES => This%PESsContainer(iPES)%PES )
        allocate( PES%mMiMn( size(This%mMiMn) ) )
        PES%mMiMn = This%mMiMn
        if (i_Debug_Loc) call Logger%Write( "This%PES(iPES)%mMiMn = ", PES%mMiMn )
      end associate
    end do


    This%RedMass    =     This%Species(iSpeTar)%Mass * This%Species(iSpePro)%Mass / ( This%Species(iSpeTar)%Mass + This%Species(iSpePro)%Mass )
    if ( This%icoord == 0 ) then
      xmu     =   This%RedMass
      if (i_Debug_Loc) call Logger%Write( "cartesian atomic coordinates will be used, however make cm stationary at origin to eliminate 3 coordinates and 3 momenta." )
    else if ( This%icoord == 1 ) then
      xmu     =   Input%AtomsMass(1)
      if (i_Debug_Loc) call Logger%Write( "cartesian jacobi coordinates will be used, with the first 3 coordinates the vector between centers of mass." )
      if ( This%Species(iSpeTar)%NAtoms /= 1 ) then
        call Logger%Write( "This%Species(iSpeTar)%NAtoms = ",This%Species(iSpeTar)%NAtoms, " which is incompatible with This%icoord=1." )
        stop
      end if
    else if ( This%icoord == 2 ) then
      xmu     =     Input%AtomsMass(This%NAtoms-1)
      if (i_Debug_Loc) call Logger%Write( "cartesian jacobi coordinates will be used, with the last 3 coordinates the vector between centers of mass." )
      if ( This%Species(iSpePro)%NAtoms /= 1 ) then
        call Logger%Write( "This%Species(iSpePro)%NAtoms = ",This%Species(iSpePro)%NAtoms," which is incompatible with This%icoord=2." )
        stop
      end if
    else
      if (i_Debug_Loc) call Logger%Write( "unknown This%icoord option in InitializeCollision: ", This%icoord )
      stop
    end if
    if (i_Debug_Loc) call Logger%Write( "Reduced mass for relative translation: xmu = ",xmu, Fr="es15.8" )

    nsum = This%Species(iSpeTar)%NAtoms + This%Species(iSpePro)%NAtoms
    if ( This%icoord == 0 ) then
      if ( This%NAtoms /= nsum ) then
        call Logger%Write( "the number of atoms do not match")
        stop
      end if
    else
      This%NormalKineticEnergy   =   .True.                               !  if This%NormalKineticEnergy=True, then the kinetic energy is "normal", ie Ekin=0.5 sum p**2/m.
      if ( This%NAtoms-1 /= nsum ) then
        call Logger%Write( "the number of atoms do not match")
        stop
      end if
      do i = 1,nsum
        This%mMiMn(i)  =   One / Input%AtomsMass(i)
      end do
      This%Species(iSpeTar)%Mass   =   xmu * Two
      This%Species(iSpePro)%Mass   =   xmu * Two
    end if
    if (i_Debug_Loc) call Logger%Write( "This%NormalKineticEnergy = ", This%NormalKineticEnergy )
  ! ==============================================================================================================


  if ( (TaskType_Loc == 1) .or. (TaskType_Loc == 2) .or. (TaskType_Loc == 3) .or. (TaskType_Loc > 5) ) then   !!! <----- If NOT RUNNING TRAJECTORIES or COMPUTING STATISTICS
    ! ==============================================================================================================
    !   7. INITIALING MOLECULES (WITH LEVELS AND BINS)
    ! ==============================================================================================================
    call This%InitializeMolecules( Input, i_Debug=i_Debug_Loc )
    ! ==============================================================================================================
  end if


  if (TaskType_Loc == 4) then   !!! <----- ONLY if RUNNING TRAJECTORIES

    ! ==============================================================================================================
    !   8. SETTING PARAMETERS SPECIFYING THE RELATIVE KINETIC ENERGY AND INITIAL IMPACT PARAMETER
    ! ==============================================================================================================
    allocate(This%ImpactPara(Input%NPESs), stat=Status)
    if (Status/=0) call Error( "Error allocating This%ImpactPara" )
    if (i_Debug_Loc) call Logger%Write( "Allocated ImpactPara" )

    if (i_Debug_Loc) call Logger%Write( "Calling InitializeImpactParameter" )
    NTrajTemp = 0
    do iPES = 1,Input%NPESs
      !if (iPES < Input%NPESs) then
        NTrajPES  = int ( Input%PES_Degeneracy(iPES) * Input%NTrajBatch )
        NTrajTemp = NTrajTemp + NTrajPES
      !else
      !  NTrajPES = Input%NTrajBatch - NTrajTemp
      !  if (i_Debug_Loc) call Logger%Write( "-> For the Last PES, allocated ", NTrajPES,  " trajectories" )
      !end if
      call This%InitializeImpactParameter( Input, iPES, NTrajPES, i_Debug=i_Debug_Loc )
    end do
    
    Input%NTrajOverall = NTrajTemp
    This%NTrajOverall  = Input%NTrajOverall
    This%iProc         = Input%iProc
    if (i_Debug_Loc) call Logger%Write( "-> This%NTrajOverall = ", This%NTrajOverall )
    if (i_Debug_Loc) call Logger%Write( "-> This%iProc        = ", This%iProc )

    if (i_Debug_Loc) call Logger%Write( "-> Done with InitializeImpactParameter" )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   9. INITIALIZING MOLECULE STATES
    ! ==============================================================================================================
    select case (This%NAtoms)
      case (3)
        if (i_Debug_Loc) call Logger%Write( "Calling InitializeLevels_DiatomAtom" )
        if ( This%ImpactPara(1)%NRings /= 0) then  ! Case when we are Initializing from Analysis
          call InitializeLevels_DiatomAtom( Input, This%Species, i_Debug=i_Debug_Loc )
        end if
        This%SetInitialState        =>    SetInitialState_DiatomAtom
        This%Compute_Rp_CartCoord   =>    Compute_Rp_CartCoord_3Atoms
      case (4)
        ! call Error( "The 4 atoms case is not fully implemented")
        if (i_Debug_Loc) call Logger%Write( "Calling InitializeLevels_DiatomDiatom" )
        if ( This%ImpactPara(1)%NRings /= 0) then  ! Case when we are Initializing from Analysis
          call InitializeLevels_DiatomDiatom( Input, This%Species, i_Debug=i_Debug_Loc )
        end if
        This%SetInitialState  =>    SetInitialState_DiatomDiatom
    end select
    if (i_Debug_Loc) call Logger%Write( "-> This%icoord = ", This%icoord )
    ! ==============================================================================================================

  end if


  if (TaskType_Loc > 3) then   !!! <----- ONLY if RUNNING TRAJECTORIES, COMPUTING STATISTICS or POSTPROCESSING CROSS SECTIONS
    ! ==============================================================================================================
    !   10. SETTING THE AVERAGE VELOCITY
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Setting parameters specifying the relative kinetic energy and initial impact parameter" )
    This%TtraModel                    =       Input%TtraModel
    This%Ttra                         =       Input%Ttra
    This%Etot                         =       Input%Etot
    This%Erel                         =       Input%Erel
    
    if ( This%TtraModel .eq. "Boltzmann" ) then

      if (i_Debug_Loc) call Logger%Write( "-> Temperature of the Boltzmann distribution for the initial relative translation energy: T [K] = ", This%Ttra, Fr = "es15.8" )
      Vp         =   sqrt( Two * This%Ttra * Kelvin_To_Hartree / xmu )
      Va         =   Two / sqrt(Pi) * Vp
      Vrms       =   sqrt(1.5_rkp)  * Vp
      This%EqVelocity =   Va * 4.134137e16_rkp * 0.5291771e-8_rkp**3  
      if (i_Debug_Loc) then
        call Logger%Write( "-> Most probable velocity:  Vp   [a.u.]   = ", Vp,              Fr="es15.8" )
        call Logger%Write( "-> Mean velocity:           Va   [a.u.]   = ", Va,              Fr="es15.8" )
        call Logger%Write( "-> Mean velocity:           Va   [cm^3/s] = ", This%EqVelocity, Fr="es15.8" )
        call Logger%Write( "-> Root mean square speed:  Vrms [a.u.]   = ", Vrms,            Fr="es15.8" )
      end if
      if (i_Debug_Loc) call Logger%Write( "Write the Velocity Scaling Factor" )  
      FileName = trim(adjustl(Input%OutputDir)) // '/Velocity_' // trim(adjustl(Input%Ttra_char)) // '.dat'
      open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,*) "#     Mean velocity [cm^3/s]   "
      if (Status/=0) call Error( "Error opening file: " // FileName )                   
        write(Unit,'(d20.10)') This%EqVelocity
      close(Unit) 

    elseif ( This%TtraModel .eq. "Gaussian" ) then
      if (i_Debug_Loc) call Logger%Write( "-> Initial relative translational kinetic energy Normally Distributed " )

      This%EMu = Input%EMu
      This%ESD = Input%ESD
      if (i_Debug_Loc) call Logger%Write( "-> Mean [Eh] = ", This%EMu, "; SD [Eh] = ", This%ESD )

    elseif ( This%TtraModel .eq. "FixedTotEn" ) then
      if (i_Debug_Loc) call Logger%Write( "-> Initial relative translational kinetic energy from total energy: Etot [Eh] = ", This%Etot, Fr="es15.8" )
  
    elseif (  This%TtraModel .eq. "Uniform" ) then
      if (i_Debug_Loc) call Logger%Write( "-> Fixed relative translational kinetic energy: Erel [Eh] = ", This%Erel, Fr="es15.8" )
  
    else
      call Logger%Write( "WARNING: No Relative Translational Energy Model Selected" )
    end if

    ! ==============================================================================================================
  end if


  if (TaskType_Loc == 4) then   !!! <----- ONLY if RUNNING TRAJECTORIES
    ! ==============================================================================================================
    !   11. FIND MOST PROBABLE TRANSLATIONAL ENERGY AT THE CURRENT ENERGY
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Finding the most probable Energy from Boltzmann Distribution at T=",This%Ttra,"K. Calling FindEnergyMaxwellianMax."  )
    !call FindEnergyMaxwellianMax( This%Ttra, MostProbEr )
    if (i_Debug_Loc) call Logger%Write( "Most Probable Energy from Boltzmann Distribution at T=",This%Ttra,"K is MostProbEr = ",MostProbEr )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   12. INITIALIZING OUTPUT FILES
    ! ==============================================================================================================
    if ( This%ImpactPara(1)%NRings /= 0 ) then   ! Only if Driver
      if (i_Debug_Loc) call Logger%Write( "Initializing the output files" )
      if (i_Debug_Loc) call Logger%Write( "-> Calling This%InitializeCollisionOutput" )
      call This%InitializeCollisionOutput( i_Debug_Loc )
      if (i_Debug_Loc) call Logger%Write( "-> Done with This%InitializeCollisionOutput" )
    end if
    ! ============================================================================================================== 

  end if


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine InitializeAtomsPairsSpecies( This, Input, i_Debug )

  use Input_Class                 ,only:  Input_Type
  use Parameters_Module           ,only:  Atom_To_Species_3At, Atom_To_Species_4At, NAtomsPerSpecies_3At, NAtomsPerSpecies_4At, AtomsToSpecies_3At, AtomsToSpecies_4At, &
                                          OppositePair_3At, OppositePair_4At, Atoms_To_Pair_3At, Atoms_To_Pair_4At, Pair_To_Atoms_3At, Pair_To_Atoms_4At

  class(Collision_Type)                     ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    iA                ! Index of atoms
  integer                                                   ::    iS                ! Index of species
  integer                                                   ::    iP                ! Index of pairs of atoms
  integer                                                   ::    i, j
  integer                                                   ::    Opposite

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeAtomsPairsSpecies" )
  !i_Debug_Loc   =     Logger%On()

  ! ==============================================================================================================
  !   3.1. SETTING THE ATOM OBJECTS (Defining Atom%Idx, Atom%Mass, Atom%Name and with This%Atoms(:)%To_Species)
  ! ==============================================================================================================
  ! The order of the individual elements within the vector of Atom object is important since both 'Pair' and
  ! 'Species' objects defined some index mapping variables between their constitutive atoms and the associated
  ! atoms in the variable 'This%Atoms'. Therefore, the following convention is used to order atoms:
  !  * for AB+C  collisions, the atoms are ordered [A,B,C]
  !  * for AB+CD collisions, the atoms are ordered [A,B,C,D].
  ! where AB is the target species and C/CD is the projectile species. The ordering thus first list the atoms
  ! from the target species, and then the ones from the projectile species.
  ! The index mapping between the Atoms to the associated Species is established when initializing the Species.
  ! ==============================================================================================================
    This%NAtoms       =   Input%NAtoms                                                              ! Getting the number of atoms involved in the collision    
    if (i_Debug_Loc) then
      call Logger%Write( "-> Number of atoms: This%NAtoms   = ", This%NAtoms   )
    end if

    if (i_Debug_Loc) call Logger%Write( "Setting Atoms objects" )
    if (i_Debug_Loc) call Logger%Write( "-> Number of Atoms: This%NAtoms = ", This%NAtoms )
    allocate( This%Atoms(This%NAtoms) )
    do iA = 1,This%NAtoms
      associate( Atom => This%Atoms(iA) )
        call Atom%Initialize( Input%AtomsMass(iA), Input%AtomsName(iA), i_Debug=i_Debug_Loc )
        Atom%Idx = iA
        if (i_Debug_Loc) call Logger%Write( "-> Atom%Idx = ", Atom%Idx, "Atom%Name = ", Atom%Name, "Atom%Mass = ", Atom%Mass, Fr="es15.8" )
      end associate
    end do
    select case (This%NAtoms)
      case (3); This%Atoms(:)%To_Species = Atom_To_Species_3At
      case (4); This%Atoms(:)%To_Species = Atom_To_Species_4At
    end select
    if (i_Debug_Loc) call Logger%Write( "-> Done setting Atom objects" )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   3.2. SETTING THE PAIRS OBJECTS (Defining Pair_To_Atoms and constructing Diatomic Potential for each of the pairs)
  ! ==============================================================================================================
  ! For AB+C collisions:
  ! The atoms are A, B and C. This is the order to the associated elements in the variables This%Atoms.
  ! For the atom-atom pairs the order is the A-B, A-C, B-C, that is:
  ! * 1) A-B = 1-2 : corresponds to the intramolecular bond of the target molecular species AB.
  ! * 2) A-C = 1-3 : is between the 1st atom A of the target species and the atom of the projectile atomic species C.
  ! * 3) B-C = 2-3 : is between the 2nd atom B of the target species and the atom of the projectile atomic species C.
  ! For AB+CD collisions:   (ordered taken for NASA N4 potential)
  ! * 1) A-B = 1-2 : corresponds to the intramolecular bond of the target molecular species AB.
  ! * 2) C-D = 3-4 : corresponds to the intramolecular bond of the projectile molecular species CD.
  ! * 3) A-C = 1-3 :
  ! * 4) A-D = 1-4 :
  ! * 5) B-C = 2-3 :
  ! * 6) B-D = 2-4 :
  ! Note that this ordering of the atom-atom pairs should be consitent with the ordering of the vector containing
  ! the bond length. This variables is computed by the PES object (in the 'potcl' procedure).
  ! When new PES object are implemented, one need to be sure that the pair order is consistent.
  ! ==============================================================================================================
    This%NPairs       =   This%NAtoms * ( This%NAtoms - 1 ) / 2                                     ! Setting the number of pairs of atoms (Identical pairs are counted several times)
    if (i_Debug_Loc) then
      call Logger%Write( "-> Number of atoms: This%NPairs   = ", This%NPairs   )
    end if

    if (i_Debug_Loc) call Logger%Write( "Setting Pairs objects" )
    if (i_Debug_Loc) call Logger%Write( "-> Number of Pairs: This%NPairs = ", This%NPairs )
    allocate( This%Pairs(This%NPairs) )
    allocate( This%OppositePair(This%NPairs) )
    allocate( This%Pair_To_Atoms(2, This%NPairs) )
    allocate( This%Atoms_To_Pair(This%NAtoms,This%NAtoms) )

    select case(This%Natoms)
      case(3)
        This%OppositePair  = OppositePair_3At

        This%Pair_To_Atoms = Pair_To_Atoms_3At

        do iP = 1,This%NPairs
           associate( Pair => This%Pairs(iP) )
             if (i_Debug_Loc) call Logger%Write( "-> Calling Pair%Initialize: iP = ", iP, "This%Pair_To_Atoms = ", This%Pair_To_Atoms(:,iP), Fi="i2" )
             call Pair%Initialize( Input, iP, This%Pair_To_Atoms(:,iP), This%Atoms, This%OppositePair(iP), i_Debug=i_Debug_Loc )
           end associate
        end do

        This%Atoms_To_Pair = Atoms_To_Pair_3At

      case(4) 
        This%OppositePair  = OppositePair_4At

        This%Pair_To_Atoms = Pair_To_Atoms_4At

        do iP = 1,This%NPairs
          associate( Pair => This%Pairs(iP) )
            if (i_Debug_Loc) call Logger%Write( "-> Calling Pair%Initialize: iP = ", iP, "This%Pair_To_Atoms = ", This%Pair_To_Atoms(:,iP), Fi="i2" )
            call Pair%Initialize( Input, iP, This%Pair_To_Atoms(:,iP), This%Atoms, This%OppositePair(iP), i_Debug=i_Debug_Loc )
          end associate
        end do 

        This%Atoms_To_Pair = Atoms_To_Pair_4At
    end select
    if (i_Debug_Loc) call Logger%Write( "-> Done setting Pair objects" )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   3.3 SETTING THE SPECIES OBJECTS (Defining NAtomsPerSpecies, AtomsToSpecies and This%NAtoMaxSpe)
  ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Setting Species objects" )
    This%NSpecies   =   size(Input%SpeciesName)
    if (i_Debug_Loc) call Logger%Write( "-> Number of Species: This%NSpecies = ", This%NSpecies )
    allocate( This%Species(This%NSpecies) )

    allocate( This%NAtomsPerSpecies(This%NSpecies) )
    allocate( This%AtomsToSpecies(2,This%NSpecies) )
    select case (This%NAtoms)
      case (3)
        This%NAtomsPerSpecies = NAtomsPerSpecies_3At
        This%AtomsToSpecies   = AtomsToSpecies_3At
      case (4)
        This%NAtomsPerSpecies = NAtomsPerSpecies_4At
        This%AtomsToSpecies   = AtomsToSpecies_4At
    end select
    This%NAtoMaxSpe     =   0
    do iS = 1,This%NSpecies
      associate( Species => This%Species(iS), AtoToSpe => This%AtomsToSpecies(1:This%NAtomsPerSpecies(iS),iS) )
        if (i_Debug_Loc) call Logger%Write( "-> Calling Species%Initialize: iS = ", iS, "Name = ", trim(Input%SpeciesName(iS)), "This%AtomsToSpecies = ", AtoToSpe(:), F2="i3", F4="a6", F6="i2" )
        call Species%Initialize( Input, iS, AtoToSpe, trim(Input%SpeciesName(iS)), This%Atoms, i_Debug=i_Debug_Loc )
        This%NAtoMaxSpe   =   max( This%NAtoMaxSpe, Species%NAtoms )
      end associate
    end do
    if (i_Debug_Loc) call Logger%Write( "-> Maximum number of atoms per species: This%NAtoMaxSpe = ", This%NAtoMaxSpe )
    if (i_Debug_Loc) call Logger%Write( "-> Done setting Species objects" )
  ! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!--------------------------------------------------------------------------------------------------------------------------------!
Subroutine InitializeMolecules( This, Input, i_Debug )

  use Input_Class                     ,only:  Input_Type
  use Molecule_Factory_Class          ,only:  Molecule_Factory_Type

  class(Collision_Type)                     ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(inout)  ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iMol
  integer                                                   ::    Status
  type(Molecule_Factory_Type)                               ::    Molecule_Factory
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeMolecules")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  

  allocate(This%MoleculesContainer(Input%NMolecules), stat=Status)
  if (i_Debug_Loc) call Logger%Write( "Allocate This%MoleculesContainer with Dimension (", Input%NMolecules, ")"  )

  do iMol = 1,Input%NMolecules
    if (i_Debug_Loc) call Logger%Write( "Molecule Nb", iMol )

    if (i_Debug_Loc) call Logger%Write( "Calling Molecule_Factory%Define_Molecule" )
    call Molecule_Factory%Define_Molecule( Input, This%NPairs, This%Pairs, This%Atoms, iMol, This%MoleculesContainer(iMol)%Molecule, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) call Logger%Write( "Done with Molecule_Factory%Define_Molecule" )

  end do

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine InitializeImpactParameter( This, Input, iPES, NTrajPES, i_Debug )
! This procedure initializes the parameters related to the impact parameter. It should only be called for the Driver.

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  rkp, Zero, One, Two, Pi
  use Global_Module         ,only:  UnitTraj

  class(Collision_Type)                                 ,intent(inout)  ::    This
  type(Input_Type)                                      ,intent(in)     ::    Input
  integer                                               ,intent(in)     ::    iPES
  integer                                               ,intent(in)     ::    NTrajPES
  logical                                     ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeImpactParameter" )
  !i_Debug_Loc   =     Logger%On()
  

! ==============================================================================================================
!     TREATING THE CASE WHEN THERE IS NO PARAMETER TO INITIALIZE
! ==============================================================================================================
  if ( Input%nring == 0) then  ! Case when we are Initializing from Analysis
    if (i_Debug_Loc) then
      call Logger%Write( "Input%nring = ", Input%nring  )
      call Logger%Write( "-> Noting to do => Exiting" )
      call Logger%Exiting
    end if
    return
  end if
! ==============================================================================================================


! ==============================================================================================================
!     SETTING PARAMETERS SPECIFYING THE RELATIVE KINETIC ENERGY
! ==============================================================================================================
! The variables 'This%TtraModel' is the model to use for describing the kinetic energy.
! It can take the following values:
!  * Boltzmann:  Then, the initial relative translational energy will be determined from a Boltzmann distribution at the temperature abs(This%Ttra) in Kelvin.
!  * FixedTotEn: Then, the initial relative translational energy will be determined from fixed total energy This%Etot, in hartree.
!  * Uniform:    Then, the initial relative translational energy corresponds to the value stored in This%Erel, in hartree.
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting parameters specifying the relative kinetic energy and initial impact parameter" )
  This%TtraModel                    =       Input%TtraModel
  This%Ttra                         =       Input%Ttra
  This%Etot                         =       Input%Etot
  This%Erel                         =       Input%Erel
  if (i_Debug_Loc) then
    
    if ( This%TtraModel .eq. "Boltzmann" ) then
      call Logger%Write( "-> Relative translational kinetic energy from Maxwell-Boltzmann distribution at translational temperature: T [K] = ",This%Ttra, Fr="es15.8" )
    
    elseif ( This%TtraModel .eq. "Gaussian" ) then
      if (i_Debug_Loc) call Logger%Write( "-> Initial relative translational kinetic energy Normally Distributed " )

      This%EMu = Input%EMu
      This%ESD = Input%ESD
      if (i_Debug_Loc) call Logger%Write( "-> Mean [Eh] = ", This%EMu, "; SD [Eh] = ", This%ESD )

    elseif (  This%TtraModel .eq. "FixedTotEn" ) then
      call Logger%Write( "-> Relative translational kinetic energy from total energy: Etot [Eh] = ", This%Etot, Fr="es15.8" )
    
    elseif (  This%TtraModel .eq. "Uniform" ) then
      call Logger%Write( "-> Fixed relative translational kinetic energy: Erel [Eh] = ", This%Erel, Fr="es15.8" )
      !This%Ttra = This%Erel
    
    else
      call Logger%Write( "WARNING: No Model Selected for the Relative Translational Energy!" )
    end if
  end if
! ==============================================================================================================


! ==============================================================================================================
!     SETTING PARAMETERS SPECIFYING THE RELATIVE KINETIC ENERGY AND INITIAL IMPACT PARAMETER
! ==============================================================================================================
! If "ImposedTotalAngularMomentum" 
!  * = true,  the Impact parameter will be chosen so that the trajectory has a particular value of the total angular momentum, "This%TotalAngularMomentum". 
!  * = false, the Impact parameter will be
!      - The same for all the trajectories                                                                    (if bstart >= 0)
!      - abs(bstart)*sqrt(eta), where eta is a random number on [0,1[                                         (if bstart <  0)
!      - Obtained from stratified sampling will be used. This requires additional input - nring,ntring,bring  (if bstart = 1234321)
! ============================================================================================================== 
  This%TotalAngularMomentum         =       Input%TotalAngularMomentum
  This%ImposedTotalAngularMomentum  =       (Input%nring == 1) .and. (Input%bring(1) == Zero)
  if (i_Debug_Loc) then
    call Logger%Write( "-> Indicator of imposed total angular momentum: This%ImposedTotalAngularMomentum = ", This%ImposedTotalAngularMomentum )
    if (This%ImposedTotalAngularMomentum) call Logger%Write( "-> Total angular momentum: This%TotalAngularMomentum = ", This%TotalAngularMomentum, Fr="es15.8" )
  end if
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE NB OF TRAJECTORIES THAT ARE STARTING FROM EACH OF THE IMPACT PARAMETERS RINGS
! ==============================================================================================================
  if ( .Not. This%ImposedTotalAngularMomentum ) then
    if (i_Debug_Loc) call Logger%Write( "Calling This%ImpactPara%Initialize" )
    call This%ImpactPara(iPES)%Initialize( Input, NTrajPES, i_Debug=i_Debug_Loc )
  end if
! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine InitializeCollisionOutput( This, i_Debug )
! This procedure initializes the collision parameters. The input parameters are loaed into the Collision object. 
! These parameters are then kept constant for all trajectories. This procedure is independent on the number of trajectories.

  class(Collision_Type)                     ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  character(:)  ,allocatable                                ::    FileName
  character(:)  ,allocatable                                ::    Form
  real(rkp) ,dimension(:)         ,allocatable              ::    Rvec
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeCollisionOutput" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine InitializeTrajectories( This, Input, Trajectories, NTrajLoc, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Trajectories_Class    ,only:  Trajectories_Type
  use Parameters_Module     ,only:  One

  class(Collision_Type)                     ,intent(in)     ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Trajectories_Type)                   ,intent(out)    ::    Trajectories
  integer                                   ,intent(in)     ::    NTrajLoc
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Trajectory" )
  !i_Debug_Loc   =     Logger%On()


  if (i_Debug_Loc) call Logger%Write( "Calling Trajectories%Initialize" )
  call Trajectories%Initialize( NTrajLoc, This%NEqtTot, This%NPairs, This%NPESs )                                   ! Must be called first

! Setting the variables from the input object
  Trajectories%tmax                    =     Input%tmax                                                             ! Setting the maximum trajectory time
  Trajectories%RpiMax                  =     One / Input%rmax                                                       ! Setting the inverse of the maximum distance between atom pair [1/bohr]
  Trajectories%BwrdIntegrationFlg      =     Input%BwrdIntegrationFlg                                               ! Setting the backward integration indicator
  
  if (i_Debug_Loc) then
    call Logger%Write( "-> Maximum trajectory time:                         Trajectories%tmax   =             ", Trajectories%tmax    , Fr="es15.8" )
    call Logger%Write( "-> Inverse of max. distance between pairs [1/bohr]: Trajectories%RpiMax =             ", Trajectories%RpiMax  , Fr="es15.8" )
    call Logger%Write( "-> Backward integration indicator:                  Trajectories%BwrdIntegrationFlg = ", Trajectories%BwrdIntegrationFlg )
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetInitialConditions( This, iTraj, Input, Traj, i_Debug, i_Debug_Deep )

  use Input_Class           ,only:  Input_Type
  use Trajectories_Class    ,only:  Trajectories_Type

  class(Collision_Type)                     ,intent(in)     ::    This      !
  integer                                   ,intent(in)     ::    iTraj     ! Index of the trajectory to be initialized
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  logical                                                   ::    i_Debug_Loc
  logical                                                   ::    i_Debug_Loc_Deep
  real(rkp)                                                 ::    StartTime, EndTime

  call CPU_Time(StartTime)

  i_Debug_Loc      = i_Debug_Global; if ( present(i_Debug) )      i_Debug_Loc      = i_Debug
  i_Debug_Loc_Deep = .False.;        if ( present(i_Debug_Deep) ) i_Debug_Loc_Deep = i_Debug_Deep
  if (i_Debug_Loc) call Logger%Entering( "SetInitialConditions" )
  !i_Debug_Loc   =     Logger%On()

  if (.not. Traj%BwrdIntegrationFlg) then
    ! ==============================================================================================================
    !     CASE OF A FORWARD INTEGRATION OF TRAJECTORIES
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Calling SetinitCondForward" )
    call This%SetIniCondForward( iTraj, Input, Traj, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Loc_Deep )
  else if ( Traj%BwrdIntegrationFlg ) then
    ! ==============================================================================================================
    !     CASE OF A BACKWARD  INTEGRATION OF TRAJECTORIES
    ! ==============================================================================================================
    ! call This%SetIniCondBackward( Input, Trajectory, PaQ, erra, dt, irec, t, smax, smin, irej, iste, i_Debug )

  end if

  call CPU_Time( EndTime )
  t_inicon  =   t_inicon + EndTime - StartTime
  i_inicon  =   i_inicon + 1
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetIniCondForward( This, iTraj, Input, Traj, i_Debug, i_Debug_Deep )
! This procedure set the initial values of the momenta and coordinates.
! This version is for arbitrary numbers of atoms, but it calls the two following procedures which depend on number of atoms:
!  * SetInitialState
!  * preini: InitializeInitialStates_3Atoms InitializeInitialStates_4Atoms
! This version performs most of the set up in scaler mode so that it makes no difference how many trajectories are run at a time.
! The momenta are stored before the coordinates.
! The variable 'imode' controls the choice of initial condiations.
!  * 0 : read in initial state. everything else is sampled using the crude monte carlo method
!  * -n: back integrate every n'th trajectory previously run.

  use Input_Class                 ,only:  Input_Type
  use Trajectories_Class          ,only:  Trajectories_Type
  use Parameters_Module           ,only:  Zero
  use StateInitDiatomAtom_Module  ,only:  SetInitialState_DiatomAtom
  use StateInitDiatomDiatom_Module,only:  SetInitialState_DiatomDiatom
  
  class(Collision_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTraj         ! Index of the trajectory to be initialized
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  real(rkp) ,dimension( NSpace, This%NAtoms )               ::    Q             ! Target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension( NSpace, This%NAtoms )               ::    dQdt          ! Time derivatives of the target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(NSpace,This%NAtoMaxSpe,This%NSpecies)::    Qijk          ! Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). 
                                                                                !    Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp) ,dimension(NSpace,This%NAtoMaxSpe,This%NSpecies)::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). 
                                                                                !    Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  logical                                                   ::    i_Debug_Loc
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetIniCondForward" )
  !i_Debug_Loc   =     Logger%On()


  Traj%NTrajInit      =   Traj%NTrajInit + 1                                    ! Incrementing the number of initialized trajectories so far. Required to set the trajectory index
  Traj%Idx(iTraj)     =   Traj%NTrajInit                                        ! Setting the global index of current trajectory (Required for the impact parameter selection)
  Traj%t(iTraj)       =   Zero                                                  ! Initializing to zero the trajectory time
  Traj%dt(iTraj)      =   Input%dt                                              ! Setting the initial timestep to the input value
  Traj%smin(iTraj)    =   huge(Zero)
  Traj%smax(iTraj)    =   -huge(Zero)
  Traj%irej(iTraj)    =   0
  Traj%iste(iTraj)    =   0
  
  Params              = Zero 
  
! ==============================================================================================================
!     SETTING PES INDEXES
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the PES Index for the current Trajectory" )
  call This%SetPESIndx( Input, iTraj, Traj, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) call Logger%Write( "Done with This%SetPESIndx" )
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING ATOMS (SPECIES COMPONENTS) INITIAL POSITIONS AND VELOCITIES (Qijk and dQijk)
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Computing initial position and velocity. Calling This%SetInitialState" )
  call This%SetInitialState( This%Species, Qijk, dQijk, iTraj, Traj%iPES(iTraj), This%NTrajOverall, This%iProc, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Deep )
  if (i_Debug_Loc) call Logger%Write( "Done with This%SetInitialState" )
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING SPECIES INITIAL POSITIONS AND VELOCITIES (Q and dQdt)
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "This%ImposedTotalAngularMomentum = ", This%ImposedTotalAngularMomentum )
  if ( .Not. This%ImposedTotalAngularMomentum ) then
    call CoordVeloc_FreeTotAngMom( This, Input, iTraj, Traj, Q, dQdt, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) call Logger%Write( "Done with CoordVeloc_FreeTotAngMom" )
  else
    call CoordVeloc_ImpTotAngMom(  This, Input, Qijk, dQijk, Traj%iPES(iTraj), Q, dQdt, Traj%b(iTraj), i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) call Logger%Write( "Done with CoordVeloc_ImpTotAngMom" )
  end if
! ==============================================================================================================


! ==============================================================================================================
!     SHIFTING COORDINATES AND VELOCITIES OF INDIVIDUAL ATOMS
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Shifting coordinates and velocities of individual atoms. Calling This%ShiftCoordinates" )
  call This%ShiftCoordinates( Q, dQdt, Qijk, dQijk, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) call Logger%Write( "Done with This%ShiftCoordinates" )
! ==============================================================================================================


! ==============================================================================================================
!     SETTING THE COORDINATES/VELOCITIES IN FINAL ARRAY AND PREPARE TO TRANSFORM VELOCITIES TO MOMENTA
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the Coordinates/Velocities in Final Array. Calling This%SetPaQ" )
  call This%SetPaQ( Qijk, dQijk, Traj%PaQ(:,iTraj) )
  if (i_Debug_Loc) then
    call Logger%Write( "Done with This%This%SetPaQ" )
    call Logger%Write( "-> Traj%PaQ(1:6, iTraj) = ", Traj%PaQ(1:6, iTraj), Fr="*(es15.8,3x)")
    call Logger%Write( "-> Traj%PaQ(7:12,iTraj) = ", Traj%PaQ(7:12,iTraj), Fr="*(es15.8,3x)")
  end if
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING THE INITIAL VALUE OF HAMILTONIAN
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Computing the initial value of Hamiltonian. Call This%Hamiltonian" )  
  call This%Hamiltonian( iTraj, Traj )
  if (i_Debug_Loc) then
    call Logger%Write( "Done with This%Hamiltonian" )
    call Logger%Write( "-> Traj%H(iTraj) = ", Traj%H(iTraj), Fr="*(es15.8,3x)")
  end if
! ==============================================================================================================
 

! ==============================================================================================================
!     SETTING THE RELATIVE ERROR FACTORS FOR INTEGRATOR
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the relative error factors for integrator. Calling This%SetErrorFactor" )
  Traj%erra(:,iTraj) = This%SetErrorFactor( Input%eps, Traj%PaQ(:,iTraj) )
  if (i_Debug_Loc) then
    call Logger%Write( "Done with This%SetErrorFactor" )
    call Logger%Write( "-> Traj%erra(1:6, iTraj) = ", Traj%erra(1:6, iTraj), Fr="*(es15.8,3x)")
    call Logger%Write( "-> Traj%erra(6:12,iTraj) = ", Traj%erra(7:12,iTraj), Fr="*(es15.8,3x)")
  end if
! ==============================================================================================================


! ==============================================================================================================
!     write entry to file specifying initial conditions
! ==============================================================================================================
  if (i_Debug_Loc) then
    call Logger%Write( "Write entry to file specifying initial conditions. Calling SetInitialPaQ" )
    call Logger%Write( "Before SetInitialPaQ" )
    call Logger%Write( "-> Traj%PaQ(1:6, iTraj) = ", Traj%PaQ(1:6, iTraj), Fr="*(es15.8,3x)")
    call Logger%Write( "-> Traj%PaQ(7:12,iTraj) = ", Traj%PaQ(7:12,iTraj), Fr="*(es15.8,3x)")
  end if
  call SetInitialPaQ( This, iTraj, Traj, dQijk )
  if (i_Debug_Loc) then
    call Logger%Write( "Done with SetInitialPaQ" )
    call Logger%Write( "Traj%H0(iTraj)        = ", Traj%H0(iTraj)        , Fr="es15.8" )
    call Logger%Write( "Traj%PaQ0(1:6, iTraj) = ", Traj%PaQ0(1:6, iTraj) , Fr="es15.8" )
    call Logger%Write( "Traj%PaQ0(7:12,iTraj) = ", Traj%PaQ0(7:12,iTraj) , Fr="es15.8" )
  end if
! ==============================================================================================================
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine CoordVeloc_FreeTotAngMom( This, Input, iTraj, Traj, Q, dQdt, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero
  use Trajectories_Class    ,only:  Trajectories_Type
  use RandomVector_Module   ,only:  RanddwsVec

  Type(Collision_Type)                      ,intent(in)     ::    This                   !
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iTraj                  ! Index of the trajectory to be initialized
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    Q                      ! Target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    dQdt                   ! Time derivatives of the target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    Er                     ! Relative translational energy
  real(rkp)                                                 ::    RandNum                ! Random number used to compute b and Erel
  real(rkp)                                                 ::    bmin
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "CoordVeloc_FreeTotAngMom" )
  !i_Debug_Loc   =     Logger%On()
  
  
! ==============================================================================================================
!     SETTING THE RELATIVE TRANSLATIONAL ENERGY
! ==============================================================================================================
  if ( This%TtraModel .eq. "Boltzmann" ) then
    RandNum = RanddwsVec(Traj%iPES(iTraj))
    if (i_Debug_Loc) call Logger%Write( "Translational Energy from Boltzmann Distribution at T=",This%Ttra,"K. Calling RelativeTranslationEnergy for Generating Random Energy."  )
    call RelativeTranslationEnergy_FromBoltzmann( This%Ttra, Er, RandNum )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )
  
  elseif ( This%TtraModel .eq. "Gaussian" ) then
      RandNum = RanddwsVec(Traj%iPES(iTraj))
      if (i_Debug_Loc) call Logger%Write( "-> Initial relative translational kinetic energy Normally Distributed " )
      call RelativeTranslationEnergy_FromGaussian( This%EMu, This%ESD, Er, RandNum )
      if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )

  elseif ( This%TtraModel .eq. "FixedTotEn" ) then
    Er = This%Etot                                                                                                                ! Original version: Er = This%Etot - Ein where Ein is undefined !!!
    if (i_Debug_Loc) call Logger%Write( "Translational Energy from Fixed Total Energy" )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )
  
  elseif (  This%TtraModel .eq. "Uniform" ) then
    Er = This%Erel
    if (i_Debug_Loc) call Logger%Write( "Translational Energy from Fixed Relative Translational Energy" )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )
  
  else
    if (i_Debug_Loc) call Logger%Write( "WARNING: No Relative Translational Energy Model Selected" )
  end if
! ==============================================================================================================


! ==============================================================================================================
!     SETTING THE IMPACT PARAMETER
! ==============================================================================================================  
  if (i_Debug_Loc) call Logger%Write( "Computing the impact parameter. Calling This%ImpactPara%GetValue" )           
  call This%ImpactPara(Traj%iPES(iTraj))%GetValue( Traj%iTrajPES(Traj%iPES(iTraj)), bmin, Traj%b(iTraj) )
  Traj%ib(iTraj) = Traj%iTrajPES(Traj%iPES(iTraj))

  if (i_Debug_Loc) call Logger%Write( "-> Impact parameter: Traj%b(iTraj) = ", Traj%b(iTraj), Fr="es15.8" )
  if ( Traj%b(iTraj) < Zero ) then
    RandNum = RanddwsVec(Traj%iPES(iTraj))

    if (trim(adjustl(Input%ImpParStrataType)) == 'Circles') then                                                                  ! Original David
      if (i_Debug_Loc) call Logger%Write( "Impact parameter is negative. bNew = - sqrt(RandNum) * bOld " )                          
      Traj%b(iTraj) =  - (RandNum)**(1.0/2.0) * Traj%b(iTraj)  
    
    elseif (trim(adjustl(Input%ImpParStrataType)) == 'Rings') then                                                                !!! SIMONE
      if (i_Debug_Loc) call Logger%Write( "Impact parameter is negative. bNew = sqrt( bmin**2 + RandNum * (Traj%b(iTraj)**2 - bmin**2) ) " )
       Traj%b(iTraj) =  sqrt( bmin**2 + RandNum * (Traj%b(iTraj)**2 - bmin**2) )
    
    elseif (trim(adjustl(Input%ImpParStrataType)) == 'Segments') then                                                             !!! SIMONE
      if (i_Debug_Loc) call Logger%Write( "Impact parameter is negative. bmin + RandNum * (- Traj%b(iTraj) - bmin) " )
       Traj%b(iTraj) = bmin + RandNum * (- Traj%b(iTraj) - bmin)
    
    end if
    if (i_Debug_Loc) call Logger%Write( "-> Impact parameter: Traj%b(iTraj) = ", Traj%b(iTraj), Fr="es15.8" )
  end if
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING THE COORDINATES AND VELOCITIES OF THE TARGET/PROJECTILE
! ==============================================================================================================
  call This%ComputeCoordinatesVelocities( Traj%b(iTraj), Er, Q, dQdt, i_Debug=i_Debug_Loc )
! ==============================================================================================================


  Params(6)     = Er
  Params(10)    = Traj%b(iTraj)


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine CoordVeloc_ImpTotAngMom( This, Input, Qijk, dQijk, iPES, Q, dQdt, b, i_Debug )
! This procedure computes the the coordinates and velocities of the target/projectile for the case where
! the total angular momentum is imposed to the given value This%TotalAngularMomentum.
! The angular momenta is computed by randomly choosing the direction of the total angular momentum vector.
! The steps are:
! * compute angular momentum of target
! * compute angular momentum of projectile
! * Computing the sum of the target/projectile angular momentum components
! * Computing the components of the orbital angular momentum

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero, One, Two, Half, TwoPi
  use RandomVector_Module   ,only:  RanddwsVec

  Type(Collision_Type)                      ,intent(in)     ::    This          !
  type(Input_Type)                          ,intent(in)     ::    Input
  real(rkp) ,dimension(:,:,:)               ,intent(in)     ::    Qijk          ! Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). ...
                                                                                !    Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp) ,dimension(:,:,:)               ,intent(in)     ::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species ...
  integer                                   ,intent(in)     ::    iPES
                                                                                !    (dim-3: 1=target, 2:projectile). Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    Q             ! Target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    dQdt          ! Time derivatives of the target coordinates for each spatial direction (dim-1) and for each atom (dim-2). ...
                                                                                !    Dim=(NSpace=3,NAtoms)
  real(rkp)                                 ,intent(out)    ::    b             ! Impact parameter
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  real(rkp)                                                 ::    Er              ! Relative translational energy
  real(rkp)                                                 ::    Term            ! Random number term
  real(rkp)                                                 ::    RandNumA        ! Random number used to compute b
  real(rkp)                                                 ::    RandNumB        ! Random number used to compute b
  real(rkp)                                                 ::    xjphi, xjcthet  ! Random direction of the total angular momentum vector
  real(rkp) ,dimension( NSpace )                            ::    TarAngMom       ! Components of the target angular momentum. Dim=(NSpace=3)
  real(rkp) ,dimension( NSpace )                            ::    ProAngMom       ! Components of the projectile angular momentum. Dim=(NSpace=3)
  real(rkp) ,dimension( NSpace )                            ::    SumAngMom       ! Components of the sum of the internal angular momenta of the target/projectile. Dim=(NSpace=3)
  real(rkp) ,dimension( NSpace )                            ::    TotAngMom       ! Components of the total angular momentum vector. Dim=(NSpace=3)
  real(rkp) ,dimension( NSpace )                            ::    OrbAngMom       ! Components of the orbital angular momentum vector. Dim=(NSpace=3)
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "CoordVeloc_ImpTotAngMom" )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!     COMPUTING THE COMPONENTS OF THE ORBITAL ANGULAR MOMENTUM
! ==============================================================================================================
  RandNumA = RanddwsVec(iPES)
  RandNumB = RanddwsVec(iPES)
  xjphi       =   RandNumA  * TwoPi
  xjcthet     =   ( RandNumB - Half ) * Two
  Term        =   sqrt(One-xjcthet**2)
  TotAngMom   =   [ Term*cos(xjphi) , Term*sin(xjphi) , xjcthet ] * This%TotalAngularMomentum                                     ! Computing the components of the total angular momentum
  
  TarAngMom   =   This%Species(iSpeTar)%AngularMomentum( Qijk(:,:,1), dQijk(:,:,1) )                                              ! Computing the angular momentum of the target
  ProAngMom   =   This%Species(iSpePro)%AngularMomentum( Qijk(:,:,2), dQijk(:,:,2) )                                              ! Computing the angular momentum of the projectile
  SumAngMom   =   ProAngMom + TarAngMom                                                                                           ! Computing the sum of the target/projectile angular momentum components

  OrbAngMom   =   TotAngMom - SumAngMom      
                                                                                                                                  ! Computing the components of the orbital angular momentum vector
  if (i_Debug_Loc) then
    call Logger%Write( "-> Random numbers:                             [xjphi,xjcthet]   = ", xjphi, xjcthet           , Fr="es15.8" )
    call Logger%Write( "-> Components of target angular momentum:      TarAngMom         = ", TarAngMom                , Fr="es15.8" )
    call Logger%Write( "-> Magnitude  of target angular momentum:      norm2(TarAngMom)  = ", norm2(TarAngMom)         , Fr="es15.8" )
    call Logger%Write( "-> Components of projectile angular momentum:  ProAngMom         = ", ProAngMom                , Fr="es15.8" )
    call Logger%Write( "-> Magnitude  of projectile angular momentum:  norm2(ProAngMom)  = ", norm2(ProAngMom)         , Fr="es15.8" )
    call Logger%Write( "-> Sum of internal angular momenta:            SumAngMom         = ", SumAngMom                , Fr="es15.8" )
    call Logger%Write( "-> Magnitude  of total angular momentum:       This%TotAngMom    = ", This%TotalAngularMomentum, Fr="es15.8" )
    call Logger%Write( "-> Components of total angular momentum:       TotAngMom         = ", TotAngMom                , Fr="es15.8" )
    call Logger%Write( "-> Components of orbital angular momentum:     OrbAngMom         = ", OrbAngMom                , Fr="es15.8" )
  end if
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING THE RELATIVE TRANSLATIONAL ENERGY
! ==============================================================================================================
  if ( This%TtraModel .eq. "Boltzmann" ) then
    RandNumA = RanddwsVec(iPES)
    if (i_Debug_Loc) call Logger%Write( "Translational Energy from Boltzmann Distribution at T=",This%Ttra,"K. Calling RelativeTranslationEnergy for Generating Random Energy."  )
    call RelativeTranslationEnergy_FromBoltzmann( This%Ttra, Er, RandNumA )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )
  
  elseif ( This%TtraModel .eq. "Gaussian" ) then
    RandNumA = RanddwsVec(iPES)
    if (i_Debug_Loc) call Logger%Write( "-> Initial relative translational kinetic energy Normally Distributed " )
    call RelativeTranslationEnergy_FromGaussian( This%EMu, This%ESD, Er, RandNumA )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )

  elseif ( This%TtraModel .eq. "FixedTotEn" ) then
    Er = This%Etot                                                                                                                ! Original version: Er = This%Etot - Ein where Ein is undefined !!!
    if (i_Debug_Loc) call Logger%Write( "Translational Energy from Fixed Total Energy" )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )
  
  elseif (  This%TtraModel .eq. "Uniform" ) then
    Er = This%Erel
    if (i_Debug_Loc) call Logger%Write( "Translational Energy from Fixed Relative Translational Energy" )
    if (i_Debug_Loc) call Logger%Write( "-> Translational Energy: Er [Eh] = ", Er, Fr="es15.8" )
  
  else
    if (i_Debug_Loc) call Logger%Write( "WARNING: No Relative Translational Energy Model Selected" )
  end if
  if ( b < Zero ) then
    RandNumA = RanddwsVec(iPES)
    b         =   -sqrt(RandNumA) * b
  end if
! ==============================================================================================================

  Params(6)     = Er
  Params(10)    = b
  Params(11:13) = TarAngMom
  Params(14:16) = ProAngMom
  Params(17:19) = TotAngMom

! ==============================================================================================================
!     COMPUTING THE COORDINATES AND VELOCITIES OF THE TARGET/PROJECTILE
! ==============================================================================================================
  call This%inia2( Er, OrbAngMom, Q, dQdt )
! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine FindEnergyMaxwellianMax( Ttra, Er )
! This procedure computes the relative translational energy from Boltzmann distribution.

  use Parameters_Module     ,only:  Zero, One, Kelvin_To_Hartree, Pi

  real(rkp)                         ,intent(in)     ::    Ttra                          !< Translational temperature [K]
  real(rkp)                         ,intent(out)    ::    Er                            !< Relative translational energy [?]

  real(rkp) ,parameter                              ::    Tolerence =   1.0E-8_rkp
  integer                                           ::    iter
  real(rkp)                                         ::    betai, err
  real(rkp)                                         ::    FuncDerRx, FuncDerLx, FuncDerNew
  real(rkp)                                         ::    ErRx, ErLx, ErNew
  real(rkp)                                         ::    sqrKT, sqrtB
  
  betai = Ttra * Kelvin_To_Hartree
  sqrKT = sqrt(Pi) * betai**2
  sqrtB = Two  * (One/betai)**(1.5d0)
  
  iter = 0
  ErRx = 1.d-3
  ErLx = 1.5d-1
  FuncDerRx  = exp(-ErRx/betai) * (One - Two * ErRx / betai) / sqrt(ErRx/(betai)) / sqrKT
  FuncDerLx  = exp(-ErLx/betai) * (One - Two * ErLx / betai) / sqrt(ErLx/(betai)) / sqrKT
  if (FuncDerRx*FuncDerLx >= Zero) call Error( "Error in Finding the Maximum of the Maxwellian Energy Distribution (FindEnergyMaxwellianMax)." )   
  do 
    iter = iter + 1
    
    ErNew      = (ErRx + ErLx) / Two
    FuncDerNew = exp(-ErNew/betai) * (One - Two * ErNew / betai) / sqrt(ErNew/(betai)) / sqrKT
    
    if (FuncDerNew > Tolerence) then
      ErRx      = ErNew
      FuncDerRx = FuncDerNew
    elseif (FuncDerNew < - Tolerence) then
      ErLx      = ErNew
      FuncDerLx = FuncDerNew
    else
      Er        = ErNew
      exit
    end if

  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine RelativeTranslationEnergy_FromBoltzmann( Ttra, Er, RandNum )
! This procedure computes the relative translational energy from Boltzmann distribution.

  use Parameters_Module     ,only:  Zero, One, Kelvin_To_Hartree, Pi

  real(rkp)                         ,intent(in)     ::    Ttra                          !< Translational temperature [K]
  real(rkp)                         ,intent(out)    ::    Er                            !< Relative translational energy [?]
  real(rkp)                         ,intent(inout)  ::    RandNum                       !< Random number

  real(rkp) ,parameter                              ::    Tolerence =   1.0E-5_rkp
  integer                                           ::    k, iter
  real(rkp)                                         ::    betai, err
  real(rkp)                                         ::    erri
  real(rkp)                                         ::    Func, FuncInt
  real(rkp)                                         ::    sqrtPi, sqrtB
  real(rkp)                                         ::    sqrtEb, expEb
  
  betai    =   Ttra * Kelvin_To_Hartree
  Er       = One
  iter     =   0
  
!!! SIMONE: Newton method !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  iter = 0
!  do 
!    iter = iter + 1

!    Func     = exp( -Er ) * Er                
!    FuncInt  = One - exp( -Er ) * (One + Er)   

!    err      = - (FuncInt - RandNum) / Func 
!    Er       = Er + err 

!    if ( abs(err) <= Tolerence ) exit
!  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! DAVID'S VERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RandNum  =   RandNum - One
  do
    iter   =   iter + 1
    
    do k = 1,4
      Er   =   ( One + Er*(One+Er) + RandNum*exp(Er) ) / Er                              
    end do

    err    =   Zero
    erri   =   ( ( One + Er*(One+Er) + RandNum*exp(Er) ) / Er ) - Er
    err    =   max(err,erri)
    
    Er     =   Er + erri
    
    if ( err <= Tolerence ) exit
  end do

  Er       =   betai * Er
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine RelativeTranslationEnergy_FromGaussian( EMu, ESD, Er, RandNum )
! This procedure computes the relative translational energy from Boltzmann distribution.

  use Parameters_Module     ,only:  Zero, One, Kelvin_To_Hartree, Pi

  real(rkp)                         ,intent(in)     ::    EMu                           !< Mean Energy [Eh]
  real(rkp)                         ,intent(in)     ::    ESD                           !< Energy Standard Deviation [Eh]
  real(rkp)                         ,intent(out)    ::    Er                            !< Relative translational energy [?]
  real(rkp)                         ,intent(inout)  ::    RandNum                       !< Random number

  real(rkp) ,parameter                              ::    Tolerance =   1.0E-13_rkp
  real(rkp)                                         ::    X, Func, FuncInt
  real(rkp)                                         ::    err
  integer                                           ::    iter
  
  Er   = EMu
  iter = 0
  do
    iter = iter + 1
    X        = (Er - EMu) / (ESD * sqrt(Two))
    Func     = exp(-(X ** 2)) / (ESD * sqrt(Two * Pi))
    FuncInt  = Half * (One + erf(X))

    err      = -(FuncInt - RandNum) / Func
    Er       = Er + err
    if ( abs(err) <= Tolerance ) exit
  end do  

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ComputeCoordinatesVelocities( This, b, Er, Q, dQdt, i_Debug )
! This procedure computes the initial coordinates and velocities for the collision of two atoms.

  use Parameters_Module     ,only:  Zero, One, Two

  class(Collision_Type)                     ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    b           !< Initial impact parameter. Dim=(NTraj)
  real(rkp)                                 ,intent(in)     ::    Er          !< Initial relative translational energy. Dim=(NTraj)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    Q           !< Target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    dQdt        !< Time derivatives of the target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer   ,parameter                                      ::    x = 1
  integer   ,parameter                                      ::    y = 2
  integer   ,parameter                                      ::    z = 3

  integer                                                   ::    i
  real(rkp)                                                 ::    Ratio          ! Ratio of the initial impact parameter over the initial separation distance
  real(rkp)                                                 ::    vEr, vErx, vErz
  real(rkp)                                                 ::    fa
  real(rkp)                                                 ::    fb
  real(rkp)                                                 ::    Qx1
  real(rkp)                                                 ::    Qx2
  real(rkp)                                                 ::    theta
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ComputeCoordinatesVelocities" )
  !i_Debug_Loc   =     Logger%On()


  fa            = - This%Species(iSpePro)%Mass / ( This%Species(iSpeTar)%Mass + This%Species(iSpePro)%Mass )
  fb            =   This%Species(iSpeTar)%Mass / ( This%Species(iSpeTar)%Mass + This%Species(iSpePro)%Mass )
  Qx1           =   fa * This%Dinit
  Qx2           =   fb * This%Dinit


  if ( b > This%Dinit ) then
    call Logger%Write( "Current value of impact parameter exceeds initial separation" )
    call Logger%Write( " -> Trajectory index:   i = ", i )
    call Logger%Write( " -> Impact parameter:   b = ", b, Fr="es15.8" )
    call Logger%Write( " -> Initial separation: D = ", This%Dinit, Fr="es15.8" )
    error stop
  end if

  vEr         =   sqrt( Two * Er / This%RedMass )           ! Computing sqrt(2*Er/mu)
  Ratio       =   b / This%Dinit                            ! Computing b/d
  if (i_Debug_Loc) then
    call Logger%Write( " -> System Recudced Mass:              This%RedMass:  = ", This%RedMass, Fr="es15.8" )
    call Logger%Write( " -> Projectile-Target Relative Speed : vEr            = ", vEr,          Fr="es15.8" )
  end if

  vErx        = - vEr * sqrt( One - Ratio**2 )
  vErz        = - vEr * Ratio

  Q(x,1)    =   Qx1                         ! Computing Q(x,1)
  Q(y,1)    =   Zero                        ! Computing ya'
  Q(z,1)    =   Zero                        ! Computing za'
  Q(x,2)    =   Qx2                         ! Computing xb'
  Q(y,2)    =   Zero                        ! Computing yb'
  Q(z,2)    =   Zero                        ! Computing zb'


  dQdt(x,1) =   fa * vErx                   ! Computing xa' dot                     ! fa and fb are used for destributing the relative velocity to the atoms velocities. ...
  dQdt(y,1) =   Zero                        ! Computing ya' dot                     ! ... This is because the frame is fixed to the center of mass.
  dQdt(z,1) =   fa * vErz                   ! Computing za' dot
  dQdt(x,2) =   fb * vErx                   ! Computing xb' dot
  dQdt(y,2) =   Zero                        ! Computing yb' dot
  dQdt(z,2) =   fb * vErz                   ! Computing zb' dot
  
  if (i_Debug_Loc) then
    call Logger%Write( "-> Q(A)    = ", Q(:,1), Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Q(B)    = ", Q(:,2), Fr="*(es15.8,3x)" )
    call Logger%Write( "-> dQdt(A) = ", dQdt(:,1), Fr="*(es15.8,3x)" )
    call Logger%Write( "-> dQdt(B) = ", dQdt(:,2), Fr="*(es15.8,3x)" )
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ShiftCoordinates( This, Q, dQ, Qijk, dQijk,  i_Debug )
! This procedure shifts the coordinates and velocities of individual atoms in the target and projectile in order for the center-of-mass to ride on the pseudo atoms.

  class(Collision_Type)                     ,intent(in)     ::    This
  real(rkp) ,dimension(:,:)                 ,intent(in)     ::    Q             ! coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(:,:)                 ,intent(in)     ::    dQ            ! Time derivatives of the coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(:,:,:)               ,intent(inout)  ::    Qijk          ! Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile).
                                                                                ! Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp) ,dimension(:,:,:)               ,intent(inout)  ::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile).
                                                                                ! Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    ia            ! Index of atoms
  integer                                                   ::    ic            ! Index of spatial direction

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ShiftCoordinates" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) then
    call Logger%Write( "-> Coordinate indicator: This%icoord = ", This%icoord )
    call Logger%Write( "Raw target coord/Veloc.:" )
    call Logger%Write( "-> Qijk( :,1,1) = ", Qijk( :,1,1), Fr="es15.8" )
    call Logger%Write( "-> Qijk( :,2,1) = ", Qijk( :,2,1), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,1,1) = ", dQijk(:,1,1), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,2,1) = ", dQijk(:,2,1), Fr="es15.8" )
    call Logger%Write( "Raw projectile coord/Veloc.:" )
    call Logger%Write( "-> Qijk( :,1,2) = ", Qijk( :,1,2), Fr="es15.8" )
    call Logger%Write( "-> Qijk( :,2,2) = ", Qijk( :,2,2), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,1,2) = ", dQijk(:,1,2), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,2,2) = ", dQijk(:,2,2), Fr="es15.8" )
    call Logger%Write( "Q:" )
    call Logger%Write( "-> Q(:,1) = ", Q(:,1), Fr="es15.8" )
    call Logger%Write( "-> Q(:,2) = ", Q(:,2), Fr="es15.8" )
    call Logger%Write( "dQ:" )
    call Logger%Write( "-> dQ(:,1) = ", dQ(:,1), Fr="es15.8" )
    call Logger%Write( "-> dQ(:,2) = ", dQ(:,2), Fr="es15.8" )
  end if

  if ( This%icoord == 0 ) then

    do ia = 1,This%Species(iSpeTar)%NAtoms
      do ic = 1,NSpace
         Qijk(ic,ia,1)    =    Qijk(ic,ia,1) +  Q(ic,1)
        dQijk(ic,ia,1)    =   dQijk(ic,ia,1) + dQ(ic,1)
      end do
    end do
    do ia = 1,This%Species(iSpePro)%NAtoms
      do ic = 1,NSpace
         Qijk(ic,ia,2)    =    Qijk(ic,ia,2) +  Q(ic,2)
        dQijk(ic,ia,2)    =   dQijk(ic,ia,2) + dQ(ic,2)
      end do
    end do

  else if ( This%icoord == 1 ) then

    do ic = 1,NSpace
       Qijk(ic,1,1)     =    Q(ic,2) -  Q(ic,1)
      dQijk(ic,1,1)     =   dQ(ic,2) - dQ(ic,1)
    end do

  else

    do ic = 1,NSpace
      Qijk( ic,1,2)     =    Q(ic,2) -  Q(ic,1)
      dQijk(ic,1,2)     =   dQ(ic,2) - dQ(ic,1)
    end do

  end if

  if (i_Debug_Loc) then
    call Logger%Write( "-> Shifted target coord/Veloc.:" )
    call Logger%Write( "-> Qijk( :,1,1) = ", Qijk( :,1,1), Fr="es15.8" )
    call Logger%Write( "-> Qijk( :,2,1) = ", Qijk( :,2,1), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,1,1) = ", dQijk(:,1,1), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,2,1) = ", dQijk(:,2,1), Fr="es15.8" )
    call Logger%Write( "-> Shifted projectile coord/Veloc.:" )
    call Logger%Write( "-> Qijk( :,1,2) = ", Qijk( :,1,2), Fr="es15.8" )
    call Logger%Write( "-> Qijk( :,2,2) = ", Qijk( :,2,2), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,1,2) = ", dQijk(:,1,2), Fr="es15.8" )
    call Logger%Write( "-> dQijk(:,2,2) = ", dQijk(:,2,2), Fr="es15.8" )
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetPaQ( This, Qijk, dQijk, PaQ )
! This procedure stores the coordinates and velocities in the final array and prepare to transform velocities into momenta.

  use Parameters_Module     ,only:  Zero

  class(Collision_Type)                     ,intent(in)     ::    This
  real(rkp) ,dimension(:,:,:)               ,intent(in)     ::    Qijk          ! Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). 
                                                                                !    Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp) ,dimension(:,:,:)               ,intent(in)     ::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). 
                                                                                !    Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp) ,dimension(:)                   ,intent(out)    ::    PaQ           ! Dim=(NEqtTot)

  integer                                                   ::    i             ! Index of elements of the 'PaQ=[P,Q]' variable
  integer                                                   ::    iP            ! Index of elements associated to momenta P within the 'PaQ=[P,Q]' variable
  integer                                                   ::    iQ            ! Index of elements associated to coordinates Q within the 'PaQ=[P,Q]' variable
  integer                                                   ::    iA            ! Index of atoms
  integer                                                   ::    iS            ! Index of spatial directions

  iP            =   1                   ! Setting the starting index of momenta P within the 'PaQ=[P,Q]' variable
  iQ            =   1 + This%NEqtVar    ! Setting the starting index of coordinates Q within the 'PaQ=[P,Q]' variable

  if ( This%icoord == 0 ) then

    do iA = 1,This%Species(iSpeTar)%NAtoms
      do iS = 1,NSpace
        PaQ(iP)   =   Zero
        PaQ(iQ)   =   Qijk(iS,iA,1)     ! This line causes an error
        iP        =   iP + 1
        iQ        =   iQ + 1
      end do
    end do

    do iA = 1,This%Species(iSpePro)%NAtoms-1
      do iS = 1,NSpace
        PaQ(iP)   =   Zero
        PaQ(iQ)   =   Qijk(iS,iA,2)
        iP        =   iP + 1
        iQ        =   iQ + 1
      end do
    end do

    iP              =   1

    do iA = 1,This%Species(iSpeTar)%NAtoms
      do iS = 1,NSpace
        do i=1,size(PaQ)/2
          PaQ(i)  =   PaQ(i) + This%Transformation%Tqp(iP,i) * dQijk(iS,iA,1)
        end do
      iP            =   iP + 1
      end do
    end do

    do iA = 1,This%Species(iSpePro)%NAtoms-1
      do iS = 1,NSpace
        do i=1,size(PaQ)/2
          PaQ(i)  =   PaQ(i) + This%Transformation%Tqp(iP,i) * dQijk(iS,iA,2)
        end do
      iP            =   iP + 1
      end do
    end do

  else

    do iA = 1,This%Species(iSpeTar)%NAtoms
      do iS = 1,NSpace
        PaQ(iP)   =   dQijk(iS,iA,1) * This%Species(iSpeTar)%Atoms(iA)%Mass
        PaQ(iQ)   =   Qijk( iS,iA,1)
        iP        =   iP + 1
        iQ        =   iQ + 1
      end do
    end do

    do iA = 1,This%Species(iSpePro)%NAtoms
      do iS = 1,NSpace
        PaQ(iP)   =   dQijk(iS,iA,2) * This%Species(iSpePro)%Atoms(iA)%Mass
        PaQ(iQ)   =   Qijk( iS,iA,2)
        iP        =   iP + 1
        iQ        =   iQ + 1
      end do
    end do

  end if

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine inia2( This, Erel, xoam, Q, dQdt )
! This procedure determines initial coordinates and velocities for the collision of Two atoms.

  use Parameters_Module     ,only:  Zero, One, Two, Four, Half

  class(Collision_Type)                     ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    Erel          !< Relative translational energy
  real(rkp) ,dimension(:)                   ,intent(in)     ::    xoam          !< Orbital angular momentum components. Dim=(NSpace=3)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    Q             !< coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms=2)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    dQdt          !< Time derivatives of the coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms=2)

  real(rkp)                                                 ::    Discriminant
  real(rkp)                                                 ::    Mt_over_SumMtMp
  real(rkp)                                                 ::    Mp_over_SumMtMp
  real(rkp)                                                 ::    MtInv
  real(rkp)                                                 ::    px1, px2, py1, py2, et1, et2
  real(rkp)                                                 ::    xjx1, xjy1, xjz1, xjx2, xjy2, xjz2
  real(rkp)                                                 ::    ddd, dv1, dv2, cm

  Mt_over_SumMtMp   =   This%Species(iSpeTar)%Mass / ( This%Species(iSpeTar)%Mass + This%Species(iSpePro)%Mass )
  Mp_over_SumMtMp   =   This%Species(iSpePro)%Mass / ( This%Species(iSpeTar)%Mass + This%Species(iSpePro)%Mass )
  MtInv             =   One / This%Species(iSpeTar)%Mass

!   set b stationary at origin; set z of a equal to Zero

!   pz (stored in dQdt)
    dQdt(3,1)    =   sqrt( xoam(1)**2 + xoam(2)**2 ) / This%Dinit
    dQdt(3,1)    = - dQdt(3,1)

    dQdt(3,1)    = - dQdt(3,1)
!
!   coordinates of a
    Q(1,1)       =  - xoam(2) / dQdt(3,1)
    Q(2,1)       =    xoam(1) / dQdt(3,1)
    Q(3,1)       =    Zero
!
!   coefficients relating px to py
    Q(1,2)       =   - xoam(2) / xoam(1)
    Q(2,2)       =   - xoam(3) * dQdt(3,1) / xoam(1)
!
!   sum of px**2 and py**2
    Q(3,2)       =   Two * This%Species(iSpeTar)%Mass * Erel - dQdt(3,1)**2
!
!   quadratic equation coefficients
    dQdt(1,2)    =   Q(1,2)**2 + One
    dQdt(2,2)    =   Two * Q(1,2) * Q(2,2)
    dQdt(3,2)    =   Q(2,2)**2 - Q(3,2)
!
!   solve for py
    Discriminant    =   dQdt(2,2)**2 - Four * dQdt(1,2) * dQdt(3,2)
    if ( Discriminant < Zero ) then
      call Logger%Write( "Discriminant was negative" )
      call Logger%Write( "-> This%Dinit = ", This%Dinit, Fr="es15.8" )
      call Logger%Write( "-> Erel = ", Erel, Fr="es15.8" )
      call Logger%Write( "-> xoam = ", xoam, Fr="es15.8" )
      write(*,*) ' inia2 Error 1'
      error stop
    end if

    Discriminant    =   sqrt(Discriminant)
    py1       =   ( - dQdt(2,2) + Discriminant ) / ( Two * dQdt(1,2) )
    py2       =   ( - dQdt(2,2) - Discriminant ) / ( Two * dQdt(1,2) )

!   solve for px
    px1       =   Q(1,2) * py1 + Q(2,2)
    px2       =   Q(1,2) * py2 + Q(2,2)

!   check roots
    et1       =   ( px1**2 + py1**2 + dQdt(3,1)**2 ) * MtInv * Half
    et2       =   ( px2**2 + py2**2 + dQdt(3,1)**2 ) * MtInv * Half

    xjx1      =    Q(2,1) * dQdt(3,1) - Q(3,1) * py1
    xjy1      =    Q(3,1) * px1 - Q(1,1) * dQdt(3,1)
    xjz1      =    Q(1,1) * py1 - Q(2,1) * px1
    xjx2      =    Q(2,1) * dQdt(3,1) - Q(3,1) * py2
    xjy2      =    Q(3,1) * px2 - Q(1,1) * dQdt(3,1)
    xjz2      =    Q(1,1) * py2 - Q(2,1) * px2

    ddd       =   sqrt(Q(1,1)**2 + Q(2,1)**2 + Q(3,1)**2)
    dv1       =   Q(1,1) * px1 + Q(2,1) * py1 + Q(3,1) * dQdt(3,1)
    dv2       =   Q(1,1) * px2 + Q(2,1) * py2 + Q(3,1) * dQdt(3,1)
    if ( dv1 * dv2 > Zero ) then
      call Logger%Write( "don''t know how to choose root" )
      write(*,*) ' inia2 Error 2'
      error stop
    else
      if ( dv1 < Zero ) then; dQdt(1:2,1) = [px1,py1]
      else;                   dQdt(1:2,1) = [px2,py2]
      end if
    end if

    Q(:,2)     =   Zero
    dQdt(:,2)  =   Zero

!   change from momenta to velocities
    dQdt(1,1)    =   dQdt(1,1) * MtInv
    dQdt(2,1)    =   dQdt(2,1) * MtInv
    dQdt(3,1)    =   dQdt(3,1) * MtInv

!   now shift cm to origin
    cm            =   Q(1,1) * Mt_over_SumMtMp + Q(1,2) * Mp_over_SumMtMp
    Q(1,1)     =   Q(1,1) - cm
    Q(1,2)     =   Q(1,2) - cm

    cm            =   Q(2,1) * Mt_over_SumMtMp + Q(2,2) * Mp_over_SumMtMp
    Q(2,1)     =   Q(2,1) - cm
    Q(2,2)     =   Q(2,2) - cm

    cm            =   Q(3,1) * Mt_over_SumMtMp + Q(3,2) * Mp_over_SumMtMp
    Q(3,1)     =   Q(3,1) - cm
    Q(3,2)     =   Q(3,2) - cm

    cm            =   dQdt(1,1) * Mt_over_SumMtMp + dQdt(1,2) * Mp_over_SumMtMp
    dQdt(1,1)  =   dQdt(1,1) - cm
    dQdt(1,2)  =   dQdt(1,2) - cm

    cm            =   dQdt(2,1) * Mt_over_SumMtMp + dQdt(2,2) * Mp_over_SumMtMp
    dQdt(2,1)  =   dQdt(2,1) - cm
    dQdt(2,2)  =   dQdt(2,2) - cm

    cm            =   dQdt(3,1) * Mt_over_SumMtMp + dQdt(3,2) * Mp_over_SumMtMp
    dQdt(3,1)  =   dQdt(3,1) - cm
    dQdt(3,2)  =   dQdt(3,2) - cm

!   end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
 Subroutine Hamiltonian_0d( This, iTraj, Traj )
! This procedure computes the Hamiltonian for the current values of the coordinates and momenta for a given trajectory.

  use Trajectories_Class    ,only:  Trajectories_Type
  use Parameters_Module     ,only:  Zero, Half, One

  class(Collision_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTraj     ! Index of the trajectory where the computation is to be done
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj

  integer                                                   ::    i, j, jmin
  real(rkp)                                                 ::    V         ! Potential. Dim=(NTraj)
  real(rkp) ,dimension( This%NEqtTot )                      ::    Eki       ! Contribution of a given spatial direction to the kinetic energyDim=(NEqtTot)
  real(rkp) ,dimension( This%NEqtVar )                      ::    Q         ! Generalized coordinates Qi [bohr]. Dim=(NEqtVar)
  real(rkp) ,dimension( This%NPairs )                       ::    Rp        ! Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension( NSpace )                            ::    Rxyz      ! Cartesian coordinates of last atom. Not used. Dim=(NSpace)
  real(rkp) ,dimension( This%NEqtVar+NSpace )               ::    QAll


! ==============================================================================================================
!   COMPUTING THE KINETIC PART OF THE HAMILTONIAN (TIMES TWO)
! ==============================================================================================================
  Traj%H(iTraj)     =   Zero
  if ( .Not. This%NormalKineticEnergy ) then
    Eki = Zero
    do i = 1,This%NEqtVar
      jmin = mod(i,3)
      if ( jmin == 0 ) jmin = 3
      do j = jmin,This%NEqtVar,3
        Eki(i) = Eki(i) + This%Transformation%Tpq(j,i) * Traj%PaQ(j,iTraj)
      end do
    end do
    do i = 1,This%NEqtVar
      Traj%H(iTraj) = Traj%H(iTraj) + Traj%PaQ(i,iTraj) * Eki(i)
    end do
  else
    do i = 1,This%NEqtVar
      j = i / 3
      if ( 3*j /= i ) j = j + 1
      Traj%H(iTraj) =  Traj%H(iTraj) + Traj%PaQ(i,iTraj)**2 * This%mMiMn(j)
    end do
  end if
! ==============================================================================================================


! ==============================================================================================================
!   COMPUTING THE POTENTIAL PART OF THE HAMILTONIAN
! ==============================================================================================================
  Q                 =   Traj%PaQ(This%NEqtVar+1:,iTraj)                                                         ! Extracting the generalized coordinates Qi. Required for the argument to 'Compute_PES' to be contiguous
  call This%Compute_Rp( Q, Rxyz, Rp )                                                                           ! Computing the cartesian coordinates of last atom (not used) and the atom-atom distances
  QAll = [Q, Rxyz]
  V                 =   This%PESsContainer(Traj%iPES(iTraj))%PES%Potential( Rp, QAll )                          ! Computing the intramolecular potential
  Traj%Rpi(:,iTraj) =   One / Rp                                                                                ! Computing the inverse of the atom-atom distances
  Traj%H(iTraj)     =   Traj%H(iTraj) * Half + V                                                                ! Adding the potential part to the Hamiltonian
! ==============================================================================================================
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine Hamiltonian_1d( This, Traj )
! This procedure computes the Hamiltonian for the current values of the coordinates and momenta for all the trajectories 1 to Traj%NTraj.
! The following components are used:
! * as inputs:  NTraj, PaQ
! * as outputs: H, Rpi

  use Trajectories_Class    ,only:  Trajectories_Type
  use Parameters_Module     ,only:  Zero, Half, One

  class(Collision_Type)                     ,intent(in)     ::    This
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj

  integer                                                   ::    i, j, jmin, k, iTraj
  integer                                                   ::    N
  real(rkp) ,dimension( Traj%NTraj )                        ::    V         ! Potential. Dim=(NTraj)
  real(rkp) ,dimension( This%NEqtTot )                      ::    Eki       ! Contribution of a given spatial direction to the kinetic energy. Dim=(NTraj,NEqtTot)
  real(rkp) ,dimension( This%NEqtVar, Traj%NTraj )          ::    Q         ! Generalized coordinates. A separate variable is required so that the argument to 'Potential_From_R' are contiguous. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension( This%NPairs,  Traj%NTraj )          ::    Rp        ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension( NSpace,       Traj%NTraj )          ::    Rxyz      ! Cartesian coordinates of last atom. Not used. Dim=(NSpace,NTraj)
  real(rkp) ,dimension( This%NEqtVar+NSpace )               ::    QAll      ! Cartesian coordinates of all the atoms. Dim=(NEqtVar+NSpace)

  N     =   Traj%NTraj

! ==============================================================================================================
!   COMPUTING THE KINETIC PART OF THE HAMILTONIAN (TIMES TWO)
! ==============================================================================================================
  if ( .Not. This%NormalKineticEnergy ) then
    do k = 1,N
      Traj%H(k) = Zero
      Eki       = Zero
      do i = 1,This%NEqtVar
        jmin = mod(i,3)
        if ( jmin == 0 ) jmin = 3
        do j = jmin,This%NEqtVar,3
          Eki(i)    =   Eki(i) + This%Transformation%Tpq(j,i) * Traj%PaQ(j,k)
        end do
      end do
      do i = 1,This%NEqtVar
        Traj%H(k)   =   Traj%H(k) + Traj%PaQ(i,k) * Eki(i)
      end do
    end do
  else
    Traj%H(1:N)       =   Zero
    do i = 1,This%NEqtVar
      j      =   i / 3
      if ( 3*j /= i ) j = j + 1
      Traj%H(1:N)   =   Traj%H(1:N) + Traj%PaQ(i,1:N)**2 *This%mMiMn(j)
    end do
  end if
! ==============================================================================================================


! ==============================================================================================================
!   COMPUTING THE POTENTIAL PART OF THE HAMILTONIAN
! ==============================================================================================================
  Q                 =   Traj%PaQ(This%NEqtVar+1:,1:N)                                                           ! Extracting the generalized coordinates Qi. Required for the argument to 'Compute_PES' to be contiguous
  call This%Compute_Rp( Q, Rxyz, Rp )                                                                           ! Computing the cartesian coordinates of last atom (not used) and the atom-atom distances
  do iTraj = 1,size(Rp,2)
    QAll     = [Q(:,iTraj), Rxyz(:,iTraj)]
    V(iTraj) = This%PESsContainer(Traj%iPES(iTraj))%PES%Potential( Rp(:,iTraj), QAll(:) )                       ! Computing the intramolecular potential
  end do
  Traj%Rpi(:,1:N)   =   One / Rp                                                                                ! Computing the inverse of the atom-atom distances
  Traj%H(1:N)       =   Traj%H(1:N) * Half + V                                                                  ! Adding the potential part to the Hamiltonian
! ==============================================================================================================

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________!
!PURITY Subroutine Compute_PES_0d( This, iPES, Q, mdVdQ, V, Rpi )
!! This procedure computes:
!!  - the distances of atom-atom pairs
!!  - the PES and
!!  - the opposite of the PES derivatives wrt to the atom-atom distances.

!  class(Collision_Type)                   ,intent(in)     ::    This
!  integer                                 ,intent(in)     ::    iPES
!  real(rkp) ,dimension(:)  CONTIGUOUS     ,intent(in)     ::    Q           !< Generalized coordinates. Dim=(NEqtVar)
!  real(rkp) ,dimension(:)  CONTIGUOUS     ,intent(out)    ::    mdVdQ       !< Negative of the PES derivatives wrt to the generalized coordinates Qi [hartree/bohr]. Dim=(NEqtVar)
!  real(rkp)                               ,intent(out)    ::    V           !< Potential
!  real(rkp) ,dimension(:)  CONTIGUOUS     ,intent(out)    ::    Rpi         !< Inverse of the distances of atom-atom pairs [1/bohr]. Dim=(NPairs)
!  
!  integer                                     ,parameter  ::    N = 1       ! Number of trajectories
!  integer   ,dimension(N)                                 ::    iPES_       ! PES Index. Dim=(NEqtVar,NTraj=1)
!  real(rkp) ,dimension(size(Q),N)                         ::    x_          ! Position. Dim=(NEqtVar,NTraj=1)
!  real(rkp) ,dimension(size(Q),N)                         ::    mdVdQ_      ! Potential derivatives wrt positions. Dim=(NEqtVar,NTraj=1)
!  real(rkp) ,dimension(N)                                 ::    V_          ! Potential. Dim=(NTraj=1)
!  real(rkp) ,dimension(size(Rpi),N)                       ::    Rpi_        ! Inverse of the distances of atom-atom pairs [1/bohr]. Dim=(NPairs,NTraj=1)

!  iPES_(N)    =   iPES
!  x_(:,N)     =   Q
!  call This%Compute_PES( iPES_, x_, mdVdQ_, V_, Rpi_ )
!  mdVdQ       =   mdVdQ_(:,N)
!  V           =   V_(N)
!  Rpi         =   Rpi_(:,N)
!  
!End Subroutine
!!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine Compute_PES_1d( This, iPES, Q, mdVdQ, V, Rpi )
! This procedure computes:
!  - the distances of atom-atom pairs
!  - the PES and
!  - the negative of the PES derivatives wrt to the atom-atom distances.

  use Parameters_Module       ,only:  One

  class(Collision_Type)                   ,intent(in)     ::    This
  integer     ,dimension(:)               ,intent(in)     ::    iPES
  real(rkp)   ,dimension(:,:)  CONTIGUOUS ,intent(in)     ::    Q           !< Generalized coordinates. Dim=(NEqtVar,NTraj)
  real(rkp)   ,dimension(:,:)  CONTIGUOUS ,intent(out)    ::    mdVdQ       !< Negative of the PES derivatives wrt to the generalized coordinates Qi [hartree/bohr]. Dim=(NEqtVar,NTraj)
  real(rkp)   ,dimension(:)    CONTIGUOUS ,intent(out)    ::    V           !< Potential [hartree]. Dim=(NTraj)
  real(rkp)   ,dimension(:,:)  CONTIGUOUS ,intent(out)    ::    Rpi         !< Inverse of the distances of atom-atom pairs [1/bohr]. Dim=(NPairs,NTraj)

  integer                                                 ::    i, iTraj    ! Index of trajectories
  real(rkp)   ,dimension(NSpace     ,size(V,1))           ::    Rxyz        ! Cartesian coordinates of last atom. Dim=(NSpace,NTraj)
  real(rkp)   ,dimension(size(Rpi,1),size(V,1))           ::    Rp          ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp)   ,dimension( This%NEqtVar+NSpace )           ::    QAll        ! Cartesian coordinates of all the atoms. Dim=(NEqtVar+NSpace)
  real(rkp)   ,dimension(size(Rpi,1),size(V,1))           ::    dVdRi       ! Derivative of the PES wrt atom-atom distances [hartree/bohr]. Dim=(NPairs,NTraj)
  logical     ,dimension(size(V,1))                       ::    CartCoordFlg

! ==============================================================================================================
!   COMPUTING THE CARTESIAN COORDINATES OF LAST ATOM, THE ATOM-ATOM DISTANCES AND ITS INVERSE
! ==============================================================================================================
  call This%Compute_Rp( Q, Rxyz, Rp )
  Rpi = One/Rp
! ==============================================================================================================

! ==============================================================================================================
!   Computing the potential and derivatives wrt atom-atom distances
! ==============================================================================================================
  
  do iTraj = 1,size(V,1)

    QAll = [Q(:,iTraj),  Rxyz(:,iTraj)]
   
    call This%PESsContainer(iPES(iTraj))%PES%Compute( Rp(:,iTraj), QAll(:), V(iTraj), dVdRi(:,iTraj), mdVdQ(:,iTraj) )

    ! Need to change sign here. PES container subroutine returns the derivatives with respect to Cartesian coordinates
    ! (for Hamilton's equations ones needs; dp_i/dt= -dV/dQ_i) 
    mdVdQ(:,iTraj) = -mdVdQ(:,iTraj)

  end do
  dVdRi = dVdRi * Rpi

  do iTraj = 1,size(V,1)
    if ( .not. This%PESsContainer(iPES(iTraj))%PES%CartCoordFlg ) then
      if ( This%NAtoms == 3 ) then
        mdVdQ(1,iTraj)  =   -dVdRi(1,iTraj) * (Q(1,iTraj)-Q(4,iTraj)) - dVdRi(2,iTraj) * (Q(1,iTraj)-Rxyz(1,iTraj)) * (one-This%mMiMn(1)) + dVdRi(3,iTraj) * (Q(4,iTraj)-Rxyz(1,iTraj))*This%mMiMn(1)
        mdVdQ(2,iTraj)  =   -dVdRi(1,iTraj) * (Q(2,iTraj)-Q(5,iTraj)) - dVdRi(2,iTraj) * (Q(2,iTraj)-Rxyz(2,iTraj)) * (one-This%mMiMn(1)) + dVdRi(3,iTraj) * (Q(5,iTraj)-Rxyz(2,iTraj))*This%mMiMn(1)
        mdVdQ(3,iTraj)  =   -dVdRi(1,iTraj) * (Q(3,iTraj)-Q(6,iTraj)) - dVdRi(2,iTraj) * (Q(3,iTraj)-Rxyz(3,iTraj)) * (one-This%mMiMn(1)) + dVdRi(3,iTraj) * (Q(6,iTraj)-Rxyz(3,iTraj))*This%mMiMn(1)
        mdVdQ(4,iTraj)  =    dVdRi(1,iTraj) * (Q(1,iTraj)-Q(4,iTraj)) + dVdRi(2,iTraj) * (Q(1,iTraj)-Rxyz(1,iTraj)) * This%mMiMn(2)       - dVdRi(3,iTraj) * (Q(4,iTraj)-Rxyz(1,iTraj))*(one-This%mMiMn(2))
        mdVdQ(5,iTraj)  =    dVdRi(1,iTraj) * (Q(2,iTraj)-Q(5,iTraj)) + dVdRi(2,iTraj) * (Q(2,iTraj)-Rxyz(2,iTraj)) * This%mMiMn(2)       - dVdRi(3,iTraj) * (Q(5,iTraj)-Rxyz(2,iTraj))*(one-This%mMiMn(2))
        mdVdQ(6,iTraj)  =    dVdRi(1,iTraj) * (Q(3,iTraj)-Q(6,iTraj)) + dVdRi(2,iTraj) * (Q(3,iTraj)-Rxyz(3,iTraj)) * This%mMiMn(2)       - dVdRi(3,iTraj) * (Q(6,iTraj)-Rxyz(3,iTraj))*(one-This%mMiMn(2))
      else
        write(*,*) 'Subroutine Compute_PES_1d from Collision_Class.F90. dV/dR -> dV/dQ still to be Implemented for 4 Atoms Systems!'
      end if
    end if
  end do

! ==============================================================================================================

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine Compute_Rp_1d( This, Q, Rxyz, Rp )

  class(Collision_Type)                   ,intent(in)     ::    This
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(in)     ::    Q           !< Generalized coordinates. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(out)    ::    Rxyz        !< Cartesian coordinates of last atom. Dim=(NSpace,NTraj)
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(out)    ::    Rp          !< Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
 
  if ( This%NAtoms == 3 ) then
    call Compute_Rp_CartCoord_3Atoms( This%mMiMn, Q, Rxyz, Rp )
  else if (This%NAtoms == 4) then 
    call Compute_Rp_CartCoord_4Atoms( This%mMiMn, Q, Rxyz, Rp )
  end if
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine Compute_Rp_0d( This, Q, Rxyz, Rp )
  class(Collision_Type)                   ,intent(in)     ::    This
  real(rkp) ,dimension(:)    CONTIGUOUS   ,intent(in)     ::    Q           !< Generalized coordinates. Dim=(NEqtVar)
  real(rkp) ,dimension(:)    CONTIGUOUS   ,intent(out)    ::    Rxyz        !< Cartesian coordinates of last atom. Dim=(NSpace)
  real(rkp) ,dimension(:)    CONTIGUOUS   ,intent(out)    ::    Rp          !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  
  integer                                     ,parameter  ::    N=1         ! Number of trajectories
  real(rkp) ,dimension(size(Q),N)                         ::    Q_          ! Generalized coordinates. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension(NSpace,N)                          ::    Rxyz_       ! Cartesian coordinates of last atom. Dim=(NSpace,NTraj)
  real(rkp) ,dimension(size(Rp),N)                        ::    Rp_         ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  
  Q_(:,N)     =   Q(:)
  if ( This%NAtoms == 3 ) then
    call Compute_Rp_CartCoord_3Atoms( This%mMiMn, Q_, Rxyz_, Rp_ )
  else if (This%NAtoms == 4) then 
    call Compute_Rp_CartCoord_4Atoms( This%mMiMn, Q_, Rxyz_, Rp_ )
  end if
  !call This%Compute_Rp( Q_, Rxyz_, Rp_ )
  Rxyz(:)     =   Rxyz_(:,N)
  Rp(:)       =   Rp_(:,N)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine Compute_Rp_CartCoord_3Atoms( mMiMn, Q, Rxyz, Rp )
! This procedure computes:
!  - the cartesian coordinates of the atoms: xyz
!  - the distances of atom-atom pairs: Rpi

  real(rkp) ,dimension(:)    CONTIGUOUS   ,intent(in)     ::    mMiMn       !< Opposite of mass ratio: -Mi/Mn, with Mn mass of the last atom. Dim=(NAtoms-1). Old name xmss
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(in)     ::    Q           !< Generalized coordinates. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(out)    ::    Rxyz        !< Cartesian coordinates of last atom. Dim=(NSpace,NTraj)
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(out)    ::    Rp          !< Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  
  integer                                                 ::    i           ! Index of trajectories

! ==============================================================================================================
!   Computing the cartesian coordinates of last atom
! ==============================================================================================================
  Rxyz(1:3,:) = mMiMn(1)*Q(1:3,:) + mMiMn(2)*Q(4:6,:)
! ==============================================================================================================

! ==============================================================================================================
!   Computing the atom-atom distances
! ==============================================================================================================
  do i = 1,size(Q,2)
    ! Pair 1 (atoms 1 and 2)
    Rp(1,i) = sqrt( ( Q(1,i) - Q(   4,i) )**2 + ( Q(2,i) - Q(   5,i) )**2 + ( Q(3,i) - Q(   6,i) )**2 )
    ! Pair 2 (atoms 2 and 3)
    Rp(2,i) = sqrt( ( Q(1,i) - Rxyz(1,i) )**2 + ( Q(2,i) - Rxyz(2,i) )**2 + ( Q(3,i) - Rxyz(3,i) )**2 )
    ! Pair 3 (atoms 1 and 3)
    Rp(3,i) = sqrt( ( Q(4,i) - Rxyz(1,i) )**2 + ( Q(5,i) - Rxyz(2,i) )**2 + ( Q(6,i) - Rxyz(3,i) )**2 )
  end do
! ==============================================================================================================

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

!________________________________________________________________________________________________________________________________!
PURITY Subroutine Compute_Rp_CartCoord_4Atoms( mMiMn, Q, Rxyz, Rp )
! This procedure computes:
!  - the cartesian coordinates of the atoms: xyz
!  - the distances of atom-atom pairs: Rpi

  real(rkp) ,dimension(:)    CONTIGUOUS   ,intent(in)     ::    mMiMn       !< Opposite of mass ratio: -Mi/Mn, with Mn mass of the last atom. Dim=(NAtoms-1). Old name xmss
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(in)     ::    Q           !< Generalized coordinates. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(out)    ::    Rxyz        !< Cartesian coordinates of last atom. Dim=(NSpace,NTraj)
  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(out)    ::    Rp          !< Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  
  integer                                                 ::    i           ! Index of trajectories

! ==============================================================================================================
!   Computing the cartesian coordinates of last atom
! ==============================================================================================================
  Rxyz(1:3,:) = mMiMn(1)*Q(1:3,:) + mMiMn(2)*Q(4:6,:) + mMiMn(3)*Q(7:9,:)
! ==============================================================================================================

! ==============================================================================================================
!   Computing the atom-atom distances (see pair information given above)
! ==============================================================================================================
  do i = 1,size(Q,2)
    ! Pair 1 (atoms 1 and 2)
    Rp(1,i) = sqrt((Q(1,i) - Q(4,i))**2 + (Q(2,i) - Q(5,i))**2 + (Q(3,i) - Q(6,i))**2)
    ! Pair 2 (atoms 1 and 3)
    Rp(2,i) = sqrt((Q(1,i) - Q(7,i))**2 + (Q(2,i) - Q(8,i))**2 + (Q(3,i) - Q(9,i))**2)
    ! Pair 3 (atoms 1 and 4)
    Rp(3,i) = sqrt((Q(1,i) - Rxyz(1,i))**2 + (Q(2,i) - Rxyz(2,i))**2 + (Q(3,i) - Rxyz(3,i))**2)
    ! Pair 4 (atoms 2 and 3)
    Rp(4,i) = sqrt((Q(4,i) - Q(7,i))**2 + (Q(5,i) - Q(8,i))**2 + (Q(6,i) - Q(9,i))**2)
    ! Pair 5 (atoms 2 and 4)
    Rp(5,i) = sqrt((Q(4,i) - Rxyz(1,i))**2 + (Q(5,i) - Rxyz(2,i) )**2 + (Q(6,i) - Rxyz(3,i))**2)
    ! Pair 6 (atoms 3 and 4)
    Rp(6,i) = sqrt((Q(7,i) - Rxyz(1,i))**2 + (Q(8,i) - Rxyz(2,i) )**2 + (Q(9,i) - Rxyz(3,i))**2)
  end do
! ==============================================================================================================

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

!________________________________________________________________________________________________________________________________!
Pure Function SetErrorFactor( Eps, PaQ ) result(ErrFac)
! This procedure sets the relative error factors for integrator
! @ERROR: Error in initial code for the relative error factors !!!!

  use Parameters_Module     ,only:  Zero

  real(rkp)                                 ,intent(in)     ::    Eps
  real(rkp) ,dimension(:)                   ,intent(in)     ::    PaQ
  real(rkp) ,dimension( size(PaQ) )                         ::    ErrFac

  integer                                                   ::    i
  integer                                                   ::    NEqtVar       ! Number of equations per variable (Coordinates/Momenta): 3*(NAtoms-1)
  real(rkp)                                                 ::    ErrP          ! Relative error factor for the coordinates
  real(rkp)                                                 ::    ErrQ          ! Relative error factor for the momenta

  NEqtVar   =   size(PaQ) / 2                                                   ! Getting the number of equations per variable (Coordinates/momenta): 3*(NAtoms-1)
  ErrP      =   Zero
  ErrQ      =   Zero

  do i = 1,NEqtVar
    ErrP    =   ErrP + PaQ(i)**2
    ErrQ    =   ErrQ + PaQ(i)**2
  end do
  
  ErrP      =   Eps * sqrt( ErrP / NEqtVar )
  ErrQ      =   Eps * sqrt( ErrQ / NEqtVar )

  do i = 1,NEqtVar
    ErrFac(i)         =   ErrP
    ErrFac(NEqtVar+i) =   ErrQ
  end do

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetInitialPaQ( This, iTraj, Traj, dQijk )
! Write entry to file specifying initial conditions

  use Trajectories_Class    ,only:  Trajectories_Type

  class(Collision_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTraj         ! Index of the trajectory to be initialized
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  real(rkp) ,dimension(:,:,:)               ,intent(in)     ::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). 

  integer                                                   ::    iA            ! Index of atoms
  integer                                                   ::    iS            ! Index of spatial directions
  integer                                                   ::    iE, jE        ! Index of equations

  Traj%H0(iTraj)      =   Traj%H(iTraj)

  iE     =   1
  if ( This%icoord == 0 ) then
    
    do iA = 1,This%Species(iSpeTar)%NAtoms
      do iS = 1,NSpace
        Traj%PaQ0(iE,iTraj) =   dQijk(iS,iA,1)
        iE                  =   iE + 1
      end do
    end do
    
    do iA = 1,This%Species(iSpePro)%NAtoms-1
      do iS = 1,NSpace
        Traj%PaQ0(iE,iTraj) =   dQijk(iS,iA,2)
        iE                  =   iE + 1
      end do
    end do
  
  else
    do jE = 1,This%NEqtVar
      Traj%PaQ0(iE,iTraj)   =   Traj%PaQ(jE,iTraj)
      iE                    =   iE + 1
    end do
  end if
  do jE = This%NEqtVar+1,This%NEqtTot
    Traj%PaQ0(iE,iTraj)     =   Traj%PaQ(jE,iTraj)
    iE                      =   iE + 1
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function ApplyTransformation( This, iTraj, Traj ) result(PaQ)
! If Collision%NormalKineticEnergy=False, then the P will be modified. The Q value always remain unchanged

  use Trajectories_Class    ,only:  Trajectories_Type
  use Parameters_Module     ,only:  Zero

  class(Collision_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTraj
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  real(rkp)   ,dimension( This%NEqtTot )                    ::    PaQ     ! Dim=(NEqtTot)

  integer                                                   ::    i, j

  PaQ     =   Traj%PaQ(:,iTraj)

  if ( .Not. This%NormalKineticEnergy ) then
    PaQ(1:This%NEqtVar)   =   Zero
    
    do i = 1,This%NEqtVar
      do j = 1,This%NEqtVar
        PaQ(i) = PaQ(i) + This%Transformation%Tpq(i,j) * Traj%PaQ(j,iTraj)
      end do
    end do
  
  end if

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function ApplyAntiTransformation( This, iTraj, Traj ) result(PaQ)
! If Collision%NormalKineticEnergy=False, then the P will be modified. The Q value always remain unchanged

  use Trajectories_Class    ,only:  Trajectories_Type
  use Parameters_Module     ,only:  Zero

  class(Collision_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTraj
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  real(rkp)   ,dimension( This%NEqtTot )                    ::    PaQ     ! Dim=(NEqtTot)

  integer                                                   ::    i, j

  PaQ     =   Traj%PaQ(:,iTraj)

  if ( .Not. This%NormalKineticEnergy ) then
    PaQ(1:This%NEqtVar)   =   Zero
    
    do i = 1,This%NEqtVar
      do j = 1,This%NEqtVar
        PaQ(i) = PaQ(i) + This%Transformation%Tqp(i,j) * Traj%PaQ(j,iTraj)
      end do
    end do
  
  end if

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetPESIndx( This, Input, iTraj, Traj, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Trajectories_Class    ,only:  Trajectories_Type
  use RandomVector_Module   ,only:  RanddwsVec
  
  class(Collision_Type)                     ,intent(in)     ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iTraj                  ! Index of the trajectory to be initialized
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    RandNum     
  real(rkp)                                                 ::    NSamplesPerPES          
  integer                                                   ::    iPES
  integer                                                   ::    NTrajPES
  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    NProcTot
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetPESIndx" )
  !i_Debug_Loc   =     Logger%On()

  !RandNum = RanddwsVec(Input%NPESs+1)
  
  iPES     = 1
  NTrajPES = 0
  do
    NTrajPES = NTrajPES + int ( Input%PES_Degeneracy(iPES) * Input%NTrajBatch )
    if (NTrajPES >= iTraj) exit
    iPES = iPES+1
    if (iPES==Input%NPESs) exit
  end do

  Traj%iPES(iTraj) = iPES
  Params(2)        = real(Traj%iPES(iTraj),rkp)
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module