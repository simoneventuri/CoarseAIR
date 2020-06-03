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

Module Integrator_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error, CheckVariable
  use ODE_Solver            ,only:  ODE_Solver_Type
  use TrajectoryPoint_Class ,only:  TrajectoryPoint_Type, ComputeCoordAndVeloc
  use File_Class            ,only:  File_Type

  implicit none

  private
  public    ::    Integrator_Type

  Type    ::    Integrator_Type
    type(ODE_Solver_Type)                   ::    ODE

!   Integration statistics variables
    logical                                    ::    BwrdIntegrationFlg
    integer                                    ::    BwrdIntegrationFreq
    integer                                    ::    N         =   0           ! At the end of the simulation N should corresponds to NTrajTot
    integer                                    ::    ireja     =   0           !< Total number of repeated steps
    integer                                    ::    istea     =   0           !< Total number of good steps
    real(rkp)                                  ::    smaxa     =   Zero        !< Average maximum time step
    real(rkp)                                  ::    smina     =   Zero        !< Average minimum time step
    real(rkp)                                  ::    Srata     =   Zero        !< Average ratio of max/min time step
    real(rkp)                                  ::    freja     =   Zero        !< Average number of repeated steps
    real(rkp)                                  ::    fstea     =   Zero        !< Average number of good steps
!   DoBackwardStuff
    real(rkp)                                  ::    htdifa    =   Zero
    real(rkp)                                  ::    htdifr    =   Zero
    real(rkp)                                  ::    qtdifa    =   Zero
    real(rkp)                                  ::    qtdifr    =   Zero
    integer                                    ::    itdif     =   0
    integer                                    ::    nogood    =   0
    integer                                    ::    iwnogood  =   0
    integer                                    ::    UnitTrajUnconv            ! File unit number where the unconverged trajectories are written
    integer                                    ::    NTrajTimeMax = 0          ! Number of trajectories whose final time exceeds the maximum time limit
!   Analysis
    type(TrajectoryPoint_Type)                 ::    TrajTemp
    type(TrajectoryPoint_Type)                 ::    TrajIni
    type(TrajectoryPoint_Type)                 ::    TrajFin
    integer                                    ::    NumNotCvg =   0           ! Number of trajectories whose energy is not conserved

    type(File_Type)                            ::    UnconvTrajsFile
    type(File_Type)                            ::    AnalysisOutputFile
    type(File_Type)                            ::    ProgressFile
    type(File_Type)                            ::    PaQSolFile
    type(File_Type) ,dimension(:) ,allocatable ::    PaQEvoFile
    type(File_Type) ,dimension(:) ,allocatable ::    XXEvoFile
    integer                                    ::    NTrajOverall
    integer                                    ::    iProc

    logical                                    ::    CMFrameFlg = .True.

  contains
    private
    procedure ,public   ::    Initialize            =>  InitializeIntegrator
    procedure ,public   ::    Integrate
    procedure           ::    AnalysisTrajPoints
    procedure           ::    UpdateStatistics      =>  UpdateIntegrationStatistics
    procedure           ::    NormalizeStatistics   =>  NormalizeIntegrationStatistics
    procedure           ::    OutputStatistics      =>  OutputIntegrationStatistics
    procedure           ::    WriteUnconvergedTrajectories
  
    procedure           ::    InitializePrinting
    procedure           ::    InitializePrintingEvo
    procedure           ::    PrintingParams
    procedure           ::    PrintingEvo
    
  End Type

  logical   ,parameter    ::    i_Debug_Global = .true.!.false.

  integer   ,parameter    ::    NSpace = 3
  
!   @TODO: To be stored in Trajectories_Type

  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeIntegrator( This, Input, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero
  use ODE_System            ,only:  HamiltonODESystem

  class(Integrator_Type)                    ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  character(:)  ,allocatable                                ::    FileName
  character(:)  ,allocatable                                ::    Format

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeIntegrator" )
  !i_Debug_Loc   =     Logger%On()

! ==============================================================================================================
!      INITIALIZING INTEGRATION STATISTICS VARIABLES
! ==============================================================================================================
  This%BwrdIntegrationFlg  =  Input%BwrdIntegrationFlg    ! Setting the backward integration indicator
  This%N          =   0
  This%ireja      =   0           !< Total number of repeated steps
  This%istea      =   0           !< Total number of good steps
  This%Smaxa      =   Zero        !< Average maximum time step
  This%Smina      =   Zero        !< Average minimum time step
  This%Srata      =   Zero        !< Average ratio of max/min time step
  This%freja      =   Zero        !< Average number of repeated steps
  This%fstea      =   Zero        !< Average number of good steps
  This%htdifa     =   Zero
  This%htdifr     =   Zero
  This%qtdifa     =   Zero
  This%qtdifr     =   Zero
  This%itdif      =   0
  This%nogood     =   0
  This%iwnogood   =   0
  This%NTrajTimeMax = 0
  This%NumNotCvg  =   0
! ==============================================================================================================


! ==============================================================================================================
!      INITIALIZING THE ODE SOLVER OBJECT
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Calling ODE%Initialize" )
  call This%ODE%Initialize( eps         =   Input%eps,          &
                            NSteps      =   Input%NSteps,       &
                            IncStpSzFlg =   Input%IncStpSzFlg,  &
                            Relax       =   Input%Relax,        &
                            NanCheck    =   .False.,            &
                            EvaluateRHS =   HamiltonODESystem,  &
                            i_Debug     =   i_Debug_Loc         )
  if (i_Debug_Loc) call Logger%Write( "-> Done initializing ODE" )
! ==============================================================================================================


! ==============================================================================================================
!      INITIALIZING THE ANALYSIS OUTPUT FILE
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Initializing the analysis output file: AnalysisOutputFile" )
  associate( File => This%AnalysisOutputFile )
    if (i_Debug_Loc) call Logger%Write( "-> Calling File%Open" )
    FileName = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/trajectories.out'
    if (i_Debug_Loc) call Logger%Write( " File%Name = ", FileName )
    call File%Open( FileName )
    if (i_Debug_Loc) call Logger%Write( "-> File%Status = ", File%Status )
    if (File%Status/=0) call Error( "Error opening file " // 'fort.7' )
    if (i_Debug_Loc) call Logger%Write( "-> Writing header" )
    if (File%Status/=0) call Error( "Error writing header in file " // 'fort.7' )
    if ( allocated(Format) ) deallocate(Format)
    allocate( Format , source = "( a8,3x,a4,3x,3(a11,3x),*(a14,3x) )" )
    if (Input%NInitMolecules.eq.1) then 
       write(File%Unit,'(A)',iostat=File%Status)  &
       '#    iTraj, iPES,       bmax,        b_i,           j_i,           v_i,         arr_i,           j_f,           v_f,         arr_f'
    elseif (Input%NInitMolecules.eq.2) then
       write(File%Unit,'(A)',iostat=File%Status)  &
       '#    iTraj, iPES,       bmax,        b_i,          j1_i,          v1_i,          j2_i,          v2_i,         arr_i,          j1_f,          v1_f,          j2_f,          v2_f,         arr_f'
    else
       call Error( "Error from 'InitializeIntegrator':: case Input%NInitMolecules > 2 NOT available!" )  
    endif
    if (File%Status/=0) call Error( "Error writing header in file " // 'fort.7' )
    if (i_Debug_Loc) call Logger%Write( "-> File%Status = ", File%Status )
  end associate
! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Integrate( This, Input, Collision, i_Debug, i_Debug_ODE )
! This procedures controls the integration of trajectories

  use, intrinsic :: IEEE_ARITHMETIC
  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Parameters_Module           ,only:  Zero, One
  use Trajectories_Class          ,only:  Trajectories_Type
  use Global_Module               ,only:  UnitTraj
  use ODE_System                  ,only:  Collision_ODE => Collision
  use String_Module               ,only:  Convert_Ratio, Convert_To_Real, Convert_To_String
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)            ,target   ,intent(inout)  ::    Collision
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_ODE

  logical                                                   ::    Restart
  integer                                                   ::    iCall
  integer                                                   ::    i
  integer                                                   ::    Status
  character(160)                                            ::    LongString
  real(rkp)                                                 ::    StartTime1, EndTime1
  real(rkp)                                                 ::    StartTime2, EndTime2
  real(rkp)                                                 ::    TrajPerSecond
  real(rkp)                                                 ::    PercentageDone
  real(rkp)                                                 ::    t_Integration
  real(rkp)                                                 ::    temp1, temp2
  character(:)  ,allocatable                                ::    FileName, String
  logical   ,dimension(:)         ,allocatable              ::    Converged               ! Dim=(NTraj)
  integer   ,dimension(:)         ,allocatable              ::    IndexConverged
  real(rkp) ,dimension(:)         ,allocatable              ::    buff                    ! Dim=(NEqtTot+2)
  type(Trajectories_Type)                                   ::    Traj
  type(Trajectories_Type)                                   ::    TrajTemp
  type(Trajectories_Type)                                   ::    TrajTemp0
  type(Trajectories_Type)                                   ::    TrajTempF
  character(*)                                  ,parameter  ::    ProcName = 'Integrate'
  integer ,parameter                                        ::    IterMax = 100000
  integer                                                   ::    Iter
  integer                                                   ::    iTraj                   ! Index of trajectory
  integer                                                   ::    NTrajTemp
  integer                                                   ::    NTrajTemp_Converged
  integer                                                   ::    NBatch
  integer                                                   ::    iBatch
  integer                                                   ::    iTrajStart
  integer                                                   ::    iTrajEnd
  integer                                                   ::    NTrajBatch
  integer                                                   ::    NTrajBatch_Converged      
  integer                                                   ::    NTrajBatch_Remaining
  integer                                                   ::    NTrajTot
  integer                                                   ::    NTrajTot_Converged      
  integer                                                   ::    NTrajTot_Remaining
  integer                                                   ::    IdxOverall

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Integrate" )
  ! i_Debug_Loc   =     Logger%On()


  NTrajTot   = Input%NTraj                                                                                                   
  if (i_Debug_Loc) call Logger%Write( "Total Nb of Trajectories:                NTrajTot   = ", NTrajTot )
  NTrajBatch = Input%NTrajBatch                                                                                   
  if (i_Debug_Loc) call Logger%Write( "Total Nb of Trajectories per Mini-Batch: NTrajBatch = ", NTrajBatch )
  NBatch     = Input%NBatch
  if (i_Debug_Loc) call Logger%Write( "Nb of Mini-Batches:                      NBatch     = ", NBatch )

  ! ==============================================================================================================
  !     INITIALIZING THE TRAJECTORIES OBJECT
  ! ==============================================================================================================
  ! The trajectories object is initialized using the local number of mini-batch trajectories, NTrajBatch.
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Initializing the trajectories object. Calling Collision%InitializeTrajectories" )
  call Collision%InitializeTrajectories( Input, Traj, NTrajBatch, i_Debug=i_Debug )
  
  if ((Input%PaQEvoFlg) .or. (Input%XXEvoFlg)) then
    if (i_Debug_Loc) call Logger%Write( "Initializing the TEMPORARY trajectories object. Calling Collision%InitializeTrajectories" )
    call Collision%InitializeTrajectories( Input, TrajTemp, NTrajBatch, i_Debug=i_Debug )
  end if
  
  if (i_Debug_Loc) call Logger%Write( "Done with Collision%InitializeTrajectories" )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   INITIALIZING THE TRAJECTORY-POINT OBJECTS ASSOCIATED TO THE INITIAL, FINAL and TEMPORARY CONDITIONS
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Initializing the Trajectory-s objects associated to the initial and final conditions" )
  if (i_Debug_Loc) call Logger%Write( "Calling This%TrajIni%Initialize" )
  call This%TrajIni%Initialize( Input, Collision, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) call Logger%Write( "Calling This%TrajFin%Initialize" )
  call This%TrajFin%Initialize( Input, Collision, i_Debug=i_Debug_Loc )
  
  if ((Input%PaQEvoFlg) .or. (Input%XXEvoFlg)) then
    if (i_Debug_Loc) call Logger%Write( "Calling This%TrajTemp%Initialize" )
    call This%TrajTemp%Initialize( Input, Collision, i_Debug=i_Debug_Loc )
  end if
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   INITIALIZATION OF LOCAL VARIABLES
  ! ==============================================================================================================
  allocate( Converged(NTrajBatch)       ); Converged      = .True.                           ! Initialization to true so that 'SetInitialConditions' is called at the 1st iteration
  allocate( IndexConverged(NTrajBatch)  ); IndexConverged = 0
  allocate( buff(Collision%NEqtTot+2) )
  
  NTrajTot_Converged  =   0     
  Collision_ODE       =>  Collision                                                          ! For the RHS porcedure oin the ODE solver
  Restart             =   .False.                                                            ! Initializing the restart indicator
  if ( Input%TimeMax < Zero ) then
    Restart   =   .True.
  end if
  
  if (i_Debug_Loc) call Logger%Write( "Initializing the Output Files. Calling This%InitializePrinting" )
  call This%InitializePrinting( Input, Collision, Traj, NTrajBatch, i_Debug=i_Debug_Loc )
  ! ==============================================================================================================
  
  
  ! ==============================================================================================================
  !   IN CASE OF BACKINTEGRATION
  ! ==============================================================================================================
  if ( Traj%BwrdIntegrationFlg ) then
    if (i_Debug_Loc) call Logger%Write( "Backward integration" )
    This%BwrdIntegrationFreq = Input%BwrdIntegrationFreq
    if (i_Debug_Loc) call Logger%Write( "-> Frequency: This%BwrdIntegrationFreq = ", This%BwrdIntegrationFreq )
    read(5,*) LongString
    allocate( FileName , source = LongString )
    open( Unit=UnitTraj, file=FileName, Form='unformatted', Access='direct', Status='old', Recl=8*(Collision%NEqtTot+2) )
    read(UnitTraj,rec=1) buff
    temp1         =   buff(1)
    temp2         =   buff(2)
    NTrajTot      =   int(temp2+0.001d0) / This%BwrdIntegrationFreq
    NTrajTot      =   max(NTrajTot,1)
  end if
  ! ==============================================================================================================


  This%iProc        = Input%iProc
  if (i_Debug_Loc) call Logger%Write( "Processor Idx:                           This%iProc        = ", This%iProc )  
  This%NTrajOverall = Input%NTrajOverall
  if (i_Debug_Loc) call Logger%Write( "Nb of Overall Trahectories in this Node: This%NTrajOverall = ", This%NTrajOverall )


  ! ==============================================================================================================                
  !  LOOP OVER MINI-BATCHES OF TRAJECTORIES
  ! *** ========================================================================================================== ! Initializing the number of converged trajectories
  iTrajStart          =   0
  call CPU_Time( StartTime1 )  
  if (i_Debug_Loc) call Logger%Write( "Loop over Batches of Trajectories", NewLine=.True. )
  MiniBatches_Loop: do iBatch = 1,NBatch
    iTrajEnd   = iTrajStart + NTrajBatch
    if (i_Debug_Loc) call Logger%Write( "Trajectories Mini-Batch Nb:               iBatch     = ", iBatch,     NewLine=.True. )
    if (i_Debug_Loc) call Logger%Write( "Mini Batch goes from Trajectory Nb ...:   iTrajStart = ", iTrajStart, NewLine=.True. )
    if (i_Debug_Loc) call Logger%Write( "... to Trajectory Nb:                     iTrajEnd   = ", iTrajEnd,   NewLine=.True. )
  
    
    ! ==============================================================================================================
    !  SETTING THE TRAJECTORIES INITIAL CONDITIONS
    ! ==============================================================================================================
    if (.Not. Restart) then
      if (i_Debug_Loc) call Logger%Write( "Number of Trajectories to be Initialized: ", Convert_Ratio(count(Converged(1:NTrajBatch)),NTrajBatch) )

      do iTraj = 1,NTrajBatch

        if (i_Debug_Loc) call Logger%Write( "Setting Initial Condition for iTraj = ", iTraj )
        call Collision%SetInitialConditions( iTraj, Input, Traj, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_ODE )

        if ((Input%PaQEvoFlg) .or. (Input%XXEvoFlg) .or. (Input%ParamsFlg)) then
          if (i_Debug_Loc) call Logger%Write( "Printing the Initial Parameters. Calling This%PrintingParams" )
          call This%PrintingParams( i_Debug=i_Debug_Loc )
        
          if (i_Debug_Loc) call Logger%Write( "Initializing the Evolution Files. Calling This%InitializePrintingEvo" )
          call This%InitializePrintingEvo( Input, Collision, Traj, TrajTemp, iTraj, i_Debug=i_Debug_Loc )
        end if
      end do
      
      if (i_Debug_Loc) call Logger%Write( "Done with Collision%SetInitialConditions", NewLine=.True. )
    end if
    ! ==============================================================================================================
  
  
    IndexConverged       = 0
    Converged            = .True.
    NTrajTemp            = NTrajBatch
    NTrajBatch_Converged = 0
    ! ==============================================================================================================                
    !  LOOP OVER TRAJECTORIES
    ! *** ========================================================================================================== 
    Trajectories_Loop: do
      flush(Logger%Unit)
      call CPU_Time( StartTime2 )
      Traj%NTraj = NTrajTemp
      TrajTemp0     = Traj
      TrajTemp0%PaQ = Traj%PaQ0


      ! ==============================================================================================================
      !  INTEGRATION LOOP
      ! *** ===========================================================================================================
      if (i_Debug_Loc) call Logger%Write( "Starting integration", NewLine=.True. )

     
      Iter =   0
      Intagration_Loop: do 
        Iter = Iter + 1
   
        if ( .Not. Restart ) then
          Outer: do iCall = 1,Input%ncall

            !if ( Input%BwrdIntegrationFlg ) then
            !  do iTraj = 1,Traj%NTraj
            !    if ( Traj%t(iTraj) + Traj%dt(iTraj) < Zero ) Traj%dt(iTraj) = - Traj%t(iTraj)
            !    if ( Traj%dt(iTraj) == Zero ) exit Outer
            !  end do
            !end if

            !do iTraj = 1,Traj%NTraj
            !  if (i_Debug_Loc) call CheckVariable( Traj%dt(iTraj), ProcName=ProcName, VarName="Traj%dt(iTraj)" )
            !end do

            call This%ODE%Integrate( Traj, i_Debug=i_Debug_ODE )                                                                    ! Integrating the trajectories
            
          end do Outer
        end if
            
        call Collision%Hamiltonian( Traj )

        Converged(1:Traj%NTraj) = Traj%AssessConvergence()
        if ((Input%PaQEvoFlg) .or. (Input%XXEvoFlg)) call This%PrintingEvo( Input, Collision, Traj, TrajTemp, Converged, i_Debug=i_Debug_Loc )


        if ( any(Converged(1:Traj%NTraj)) ) exit Intagration_Loop                                   ! If at least one trajectory is converged, then exiting the integration loop
      end do Intagration_Loop                                                                                                                    
      if ( Iter == IterMax ) call Error( "Error: In Integrate, the maximum number of iteration has been reached" )
      ! *** ===========================================================================================================
      
      
      ! ==============================================================================================================
      do iTraj = 1,Traj%NTraj
        if ( .Not. Converged(iTraj) ) cycle                                                                                         ! If current trajectory is not converged, then cycle
        
        Traj%PaQ(:,iTraj) = Collision%ApplyTransformation( iTraj, Traj )
        
        !if ( Input%BwrdIntegrationFlg ) call Traj%DoBackwardStuff( iTraj, PaQ(:), itemp )
        IdxOverall = (This%NTrajOverall) * (This%iProc-1) + Traj%Idx(iTraj)


        if (Input%PaQSolFlg) write(This%PaQSolFile%Unit,This%PaQSolFile%Format) IdxOverall, Traj%t(iTraj), Traj%H0(iTraj), TrajTemp0%PaQ(:,iTraj), Traj%H(iTraj), Traj%PaQ(:,iTraj)
        
        
        call This%UpdateStatistics( iTraj, Traj )
        
        call This%AnalysisTrajPoints( iTraj, iTrajStart, Collision, Traj, i_Debug=i_Debug_Loc )
        
      end do
      ! ==============================================================================================================


      ! ==============================================================================================================
      !  CHECK TO SEE IF WE SHOULD PICK UP MORE TRAJECTORIES
      ! ==============================================================================================================
      ! When counting the number of converged local trajectories NTrajBatch_Converged, one must specified the current
      ! number of local trajectories in the iDone variable, iDone(1:NTrajBatch), in order not to count trajectories
      ! which have already be done
      ! ==============================================================================================================
      NTrajTemp_Converged  =   count( Converged(1:Traj%NTraj) )                                               ! Getting the local number of converged trajectories  
      NTrajBatch_Converged =   NTrajTemp_Converged + NTrajBatch_Converged                                     ! Updating the total number of converged trajectories
      NTrajBatch_Remaining =   NTrajBatch - NTrajBatch_Converged                                              ! Setting the number of trajectories still to be computed
      NTrajTot_Converged   =   NTrajTot_Converged + NTrajBatch_Converged
      NTrajTot_Remaining   =   NTrajTot - NTrajTot_Converged
      NTrajTemp            =   NTrajBatch_Remaining

      
      if (i_Debug_Loc) then
        call Logger%Write( "-> NTrajTemp_Converged  = ", NTrajTemp_Converged )
        call Logger%Write( "-> NTrajBatch_Converged = ", NTrajBatch_Converged   )
        call Logger%Write( "-> NTrajBatch_Remaining = ", NTrajBatch_Remaining   )
        call Logger%Write( "-> NTrajTemp            = ", NTrajTemp           )
        call Logger%Write( "-> Number of Local trajectories:  Converged = ", Convert_Ratio(NTrajBatch_Converged,NTrajBatch), "Remaining = ", Convert_Ratio(NTrajBatch_Remaining,NTrajBatch) )
        call Logger%Write( "-> Number of Total trajectories:  Converged = ", Convert_Ratio(NTrajTot_Converged,NTrajTot),     "Remaining = ", Convert_Ratio(NTrajTot_Remaining,NTrajTot) )
      end if
      
      
      if ( NTrajBatch_Remaining <= 0 ) exit Trajectories_Loop                                                 ! If all the trajectory have been proceesed, then exit the main loop
      if ( NTrajBatch_Remaining < NTrajBatch ) then
        call Traj%Reorder( Converged, NTrajBatch_Remaining )                                                  ! If there are still some trajectories to add, then reordering
      end if
      ! ==============================================================================================================
      
      
    end do Trajectories_Loop
    ! *** ==========================================================================================================
    !      OUTPUT PROGRESS INFORMATION
    ! ==============================================================================================================
    call CPU_Time( EndTime2 )
    t_Integration   =   EndTime2 - StartTime2
    TrajPerSecond   =   NTrajBatch_Converged / t_Integration
    PercentageDone  =   100.0*NTrajTot_Converged / NTrajTot
    if (i_Debug_Loc) call Logger%Write(  "Progress: ", Convert_Ratio(NTrajTot_Converged,NTrajTot)//" trajectories", "=> ", PercentageDone, Fr="f6.2,'%'" )
    if (This%ProgressFile%Active) write(This%ProgressFile%Unit,This%ProgressFile%Format) NTrajTot_Converged, PercentageDone, TrajPerSecond, NTrajTot_Remaining/TrajPerSecond, t_Integration
    ! ==============================================================================================================
    
    
    ! ==============================================================================================================
    !      OUTPUT PROGRESS INFORMATION
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Resetting the Trajectory object; Calling Traj%ReSet")
    call Traj%ReSet(NTrajBatch)
    ! ==============================================================================================================
    
    
    iTrajStart = iTrajStart + NTrajBatch
    if (i_Debug_Loc) call Logger%Write( "Mini-Batch Nb ", iBatch ," is Done!", NewLine=.True. )
  end do MiniBatches_Loop
  if (i_Debug_Loc) call Logger%Write( "All Trjs have been processed!", NewLine=.True. )
  ! *** ==========================================================================================================


  call CPU_Time( EndTime1 )
  t_Integration   =   EndTime1 - StartTime1
  TrajPerSecond   =   NTrajTot_Converged / t_Integration
  if (i_Debug_Loc) call Logger%Write( "Total Integration Time: ", t_Integration, Fr="es15.8" )
  if (i_Debug_Loc) call Logger%Write( "Trjs per Second:  ",       TrajPerSecond, Fr="es15.8" )


  ! ==============================================================================================================
  !  NORMALIZING THE STATISTICS
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Calling This%NormalizeStatistics" )
  call This%NormalizeStatistics( )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !  WRITING IN LOGGER the STATISTICS ABOUT TRAJECTORIES 
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Calling This%OutputStatistics" )
  call This%OutputStatistics( i_Debug=i_Debug_Loc )
  ! ==============================================================================================================


  nullify( Collision_ODE )

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine NormalizeIntegrationStatistics( This )
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  
  This%Smina      =   This%Smina         / This%N
  This%Smaxa      =   This%Smaxa         / This%N
  This%Srata      =   This%Srata         / This%N
  This%freja      =   dfloat(This%ireja) / This%N
  This%fstea      =   dfloat(This%istea) / This%N
  
  if ( This%BwrdIntegrationFlg ) then
    This%htdifa   =   This%htdifa / This%itdif
    This%htdifr   =   This%htdifr / This%itdif
    This%qtdifa   =   This%qtdifa / This%itdif
    This%qtdifr   =   This%qtdifr / This%itdif
  end if
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine OutputIntegrationStatistics( This, i_Debug )

  class(Integrator_Type)                    ,intent(in)     ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "OutputIntegrationStatistics" )
  ! i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Write( "Averaged integration statistics:" )
  if (i_Debug_Loc) call Logger%Write( "-> Average minimum time step:           This%Smina  = ", This%Smina  , Fr="es15.8" )
  if (i_Debug_Loc) call Logger%Write( "-> Average maximum time step:           This%Smaxa  = ", This%Smaxa  , Fr="es15.8" )
  if (i_Debug_Loc) call Logger%Write( "-> Average ratio of max/min time step:  This%Srata  = ", This%Srata  , Fr="es15.8" )
  if (i_Debug_Loc) call Logger%Write( "-> Average number of repeated steps:    This%freja  = ", This%freja  , Fr="es15.8" )
  if (i_Debug_Loc) call Logger%Write( "-> Average number of good steps:        This%fstea  = ", This%fstea  , Fr="es15.8" )
  
  if ( This%BwrdIntegrationFlg ) then
    if (i_Debug_Loc) call Logger%Write( "Averaged back integrated trajectories:" )
    if (i_Debug_Loc) call Logger%Write( "-> Average abs energy difference:       This%htdifa = ", This%htdifa , Fr="es15.8" )
    if (i_Debug_Loc) call Logger%Write( "-> Average rel energy difference:       This%htdifr = ", This%htdifr , Fr="es15.8" )
    if (i_Debug_Loc) call Logger%Write( "-> Average rms p&q difference:          This%qtdifa = ", This%qtdifa , Fr="es15.8" )
    if (i_Debug_Loc) call Logger%Write( "-> Average rel p&q difference:          This%qtdifr = ", This%qtdifr , Fr="es15.8" )
  end if
  
  if ( This%NTrajTimeMax > 0 ) then
    if (i_Debug_Loc) call Logger%Write( "<WARNING> There are some trajectories whose final time exceeds the maximum time limit" )
    if (i_Debug_Loc) call Logger%Write( "<WARNING> -> Number of trajectories with t>tmax: This%NTrajTimeMax = ", This%NTrajTimeMax )
    if (i_Debug_Loc) call Logger%Write( "<WARNING> -> File storing the non-converged trajectories: Unconverged_trajectories.out" )
  end if
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine UpdateIntegrationStatistics( This, iTraj, Traj )

  use Trajectories_Class    ,only:  Trajectories_Type
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iTraj     ! Index of the trajectory being processed
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  
  This%N        =   This%N     + 1
  This%Smaxa    =   This%Smaxa + Traj%Smax(iTraj)
  This%Smina    =   This%Smina + Traj%Smin(iTraj)
  This%Srata    =   This%Srata + Traj%Smax(iTraj) / Traj%Smin(iTraj)
  This%ireja    =   This%ireja + Traj%irej(iTraj)
  This%istea    =   This%istea + Traj%iste(iTraj)
  
  call This%WriteUnconvergedTrajectories( iTraj, Traj )
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteUnconvergedTrajectories( This, iTraj, Traj )
  
  use Trajectories_Class    ,only:  Trajectories_Type
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iTraj     ! Index of the trajectory being processed
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  
  integer                                                   ::    Status
  character(*)  ,parameter                                  ::    FileName = 'Unconverged_trajectories.out' ! Name of the file where the unconverged trajectories are written
  
  if ( Traj%t(iTraj) <= Traj%tmax ) return
  
  This%NTrajTimeMax   =   This%NTrajTimeMax + 1
  if ( This%NTrajTimeMax == 1 ) then
    
    This%UnconvTrajsFile%Name =  FileName
    open( NewUnit=This%UnconvTrajsFile%Unit, File=This%UnconvTrajsFile%Name, Action="WRITE", iostat=This%UnconvTrajsFile%Status )
    if (This%UnconvTrajsFile%Status /= 0) call Error( "Error opening file: " // This%UnconvTrajsFile%Name )
    This%UnconvTrajsFile%Format = "(2x,2(i15,3x),(3x,es15.8:))"
    write(This%UnconvTrajsFile%Unit,"(A)") "#          %NTraj               Idx                    t"
    write(This%UnconvTrajsFile%Unit,This%UnconvTrajsFile%Format) This%NTrajTimeMax, Traj%Idx(iTraj), Traj%t(iTraj)
  
  else 
  
    write(This%UnconvTrajsFile%Unit,This%UnconvTrajsFile%Format) This%NTrajTimeMax, Traj%Idx(iTraj), Traj%t(iTraj)
  
  end if
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine AnalysisTrajPoints( This, iTraj, iTrajStart, Collision, Traj, i_Debug )

  use Collision_Class       ,only:  Collision_Type
  use Trajectories_Class    ,only:  Trajectories_Type
  use Parameters_Module     ,only:  Half

  class(Integrator_Type)                    ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iTraj
  integer                                   ,intent(in)     ::    iTrajStart
  type(Collision_Type)                      ,intent(in)     ::    Collision
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  logical                                       ,parameter  ::    WriteData = .True.
  real(rkp)                                     ,parameter  ::    HamAbsTol = 1.0E-5_rkp                                          !< Absolute tolerance for comparion the inital and final Hamiltonial values
  integer                                                   ::    iM          ! Index of molecuels
  integer                                                   ::    ib
!  real(rkp)                                                 ::    StartTime, EndTime
  real(rkp)                                                 ::    dot3
  real(rkp)                                                 ::    dEkin
  character(3)                                              ::    econs
  real(rkp), dimension(2)                                   ::    iTypeIni
  real(rkp), dimension(2)                                   ::    iTypeFin
  real(rkp)                                                 ::    bmin
  real(rkp)                                                 ::    bmax       
  integer                                                   ::    IdxOverall
  character(*)                                  ,parameter  ::    Format = "( i8,3x,i4,3x,3(es11.5,3x),2(4(i3,3x)),2(es14.8,3x),i9 )"

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AnalysisTrajPoints" )

!  call CPU_Time( StartTime )
  

  ! ==============================================================================================================
  !     SETTING THE TRAJECTORYPOINT OBJECTS ASSOCIATED TO THE INITIAL AND FINAL CONDITIONS
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the TrajectoryPoint objects associated to the initial and final conditions" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling This%TrajIni%SetData and This%TrajFin%SetData" )
  call This%TrajIni%SetData( Zero,          Traj%H0(iTraj), Traj%PaQ0(:,iTraj) )
  call This%TrajFin%SetData( Traj%t(iTraj), Traj%H (iTraj), Traj%PaQ (:,iTraj) )
  if (i_Debug_Loc) then
    call Logger%Write( "-> Initial Condition: t         = ", This%TrajIni%t,         Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Initial Condition: H         = ", This%TrajIni%H,         Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Initial Condition: PaQ(1:6)  = ", This%TrajIni%PaQ(1:6),  Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Initial Condition: PaQ(7:12) = ", This%TrajIni%PaQ(7:12), Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Final Condition:   t         = ", This%TrajFin%t,         Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Final Condition:   H         = ", This%TrajFin%H,         Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Final Condition:   PaQ(1:6)  = ", This%TrajFin%PaQ(1:6),  Fr="*(es15.8,3x)" )
    call Logger%Write( "-> Final Condition:   PaQ(7:12) = ", This%TrajFin%PaQ(7:12), Fr="*(es15.8,3x)" )
  end if
  ! ==============================================================================================================


  econs     =   ''
  if ( abs(This%TrajIni%H-This%TrajFin%H) > HamAbsTol ) then
    if (i_Debug_Loc) then
      call Logger%Write( "<WARNING> Energy not conserved!" )
      call Logger%Write( "-> Hf - Hi = ", This%TrajIni%H - This%TrajFin%H, Fr="es15.8" )
    end if
    econs          = 'not'
    This%NumNotCvg = This%NumNotCvg + 1
  end if


  ! ==============================================================================================================
  !     ANALYZING THE INITIAL TRAJECTORY POINT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Processing the INITIAL Trajectory Point. Calling This%TrajIni%Analyze" )
  call This%TrajIni%Analyze( Collision, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) then
    associate( Tini => This%TrajIni )
      call Logger%Write( "Properties of the INITIAL Trajectory Point: TrajIni" )
      call Logger%Write( "-> t          = ", Tini%t          , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> H          = ", Tini%H          , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> PaQ(1:6)   = ", Tini%PaQ(1:6)   , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> PaQ(7:12)  = ", Tini%PaQ(7:12)  , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> Vrel       = ", Tini%Vrel       , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> xkin       = ", Tini%xkin       , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> b          = ", Tini%b          , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> NMolecules = ", Tini%NMolecules )
      do iM = 1,Tini%NMolecules
        associate( Molecule => Tini%Molecules(iM) )
          call Logger%Write( " --- iM = ", iM, "viba = ", Molecule%viba, "AngMom = ", Molecule%AngMom, "Eint = ", Molecule%Eint, "itype = ", Molecule%itype, " ---> ", F2="i3", F4="f6.2", F6="f6.2", F8="es15.8" )
        end associate
      end do
    end associate
  end if
  ! ==============================================================================================================


  ! ==============================================================================================================
  !     ANALYZING THE FINAL TRAJECTORY POINT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( " " )
  if (i_Debug_Loc) call Logger%Write( "Processing the FINAL Trajectory Point. Calling This%TrajFin%Analyze" )
  call This%TrajFin%Analyze( Collision, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) then
    associate( Tfin => This%TrajFin )
      call Logger%Write( "Properties of the FINAL trajectory point: TrajFin" )
      call Logger%Write( "-> t          = ", Tfin%t          , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> H          = ", Tfin%H          , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> PaQ(1:6)   = ", Tfin%PaQ(1:6)   , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> PaQ(7:12)  = ", Tfin%PaQ(7:12)  , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> Vrel       = ", Tfin%Vrel       , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> xkin       = ", Tfin%xkin       , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> b          = ", Tfin%b          , Fr="*(es15.8,3x)" )
      call Logger%Write( "-> NMolecules = ", Tfin%NMolecules )
      do iM = 1,Tfin%NMolecules
        associate( Molecule => Tfin%Molecules(iM) )
          call Logger%Write( " ---> iM = ", iM, "viba = ", Molecule%viba, "AngMom = ", Molecule%AngMom, "Eint = ", Molecule%Eint, "itype = ", Molecule%itype, F2="i3", F4="f6.2", F6="f6.2", F8="es15.8" )
        end associate
      end do
    end associate
  end if
  ! ==============================================================================================================

  dot3    =   dot_product( This%TrajIni%Vrel , This%TrajFin%Vrel ) / sqrt( sum(This%TrajIni%Vrel**2) * sum(This%TrajFin%Vrel**2) )
  dEkin   =   This%TrajIni%xkin - This%TrajFin%xkin   ! change in internal energy: recall ein0 + This%TrajIni%xkin = etot0 = etotf = einf + This%TrajFin%xkin,   so einf - ein0 = This%TrajIni%xkin - This%TrajFin%xkin

  if ((Traj%t(iTraj) <= Traj%tmax) .and. (WriteData)) then

    if (size(This%TrajIni%Molecules).gt.2) call Error( "Error: in 'AnalysisTrajPoints', maximum number of molecules MUST be 2!")

    ! ==============================================================================================================
    !     SETTING THE ARRANGEMENT QUANTUM NUMBERS
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( " " )
    if (i_Debug_Loc) call Logger%Write( "Setting the arrangement quantum numbers" )
    do iM = 1,size(This%TrajIni%Molecules) 
       iTypeIni(iM) = dfloat(This%TrajIni%Molecules(iM)%itype)
       iTypeFin(iM) = dfloat(This%TrajFin%Molecules(iM)%itype)
    enddo
    if (i_Debug_Loc) then
      call Logger%Write( "-> iTypeIni = ", iTypeIni )
      call Logger%Write( "-> iTypeFin = ", iTypeFin )
    end if
    do iM = 1,size(This%TrajIni%Molecules)
       if ( iTypeIni(iM) > 100.0_rkp ) call Error( "Error: Bad arrangement quantum number of initial point")
       if ( iTypeFin(iM) > 100.0_rkp ) call Error( "Error: Bad arrangement quantum number of final point")
    enddo
    ! ==============================================================================================================


    ! ==============================================================================================================
    !     GETTING THE MAXIMUM IMPACT PARAMETER USED FOR CURRENT TRAJECTORY
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Getting the maximum impact parameter used for current trajectory. Calling Collision%ImpactPara%GetValue" )
    !call Collision%ImpactPara(Traj%iPES(iTraj))%GetValue( Traj%Idx(iTraj), bmax ); bmax = abs(bmax)                              ! Original by David and Bruno
    ib = Traj%ib(iTraj) - 1

    call Collision%ImpactPara(Traj%iPES(iTraj))%GetValue( ib, bmin, bmax ); bmax = abs(bmax)                      
    if (i_Debug_Loc) call Logger%Write( "-> Traj%Idx(iTraj) = ", Traj%Idx(iTraj), "; ib = ", ib ,"bmax = ", bmax, Fr="es15.8" )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !     WRITING TRAJECTORY'S FINAL STATE
    ! ==============================================================================================================
    ! -> for 3 ATOMS one writes 16*iPair + p, with  
    !   - p = 0 bound
    !   - p = 1 quasi-bound
    !   - p = 2 pre-dissociated
    ! -> for 4 ATOMS one writes 16*iPair + p1 + 4*p2, with p1 and p2 refer to the first and second molecule 
    !    and may have the sames values as 'p' in the 3 ATOMS case 
    ! (convention taken from David Schwenke's QCT code)
    IdxOverall = (This%NTrajOverall) * (This%iProc-1) + Traj%Idx(iTraj)
    if ( Collision%NAtoms == 3 ) then

      write(This%AnalysisOutputFile%Unit, 1)                                                                        &
                                              IdxOverall,                                                       ',',&
                                              Traj%iPES(iTraj),                                                 ',',&
                                              bmax,                                                             ',',& 
                                              This%TrajIni%b,                                                   ',',& 
                                              This%TrajIni%Molecules(1)%AngMom,                                 ',',&
                                              This%TrajIni%Molecules(1)%viba,                                   ',',& 
                                              This%TrajIni%Molecules(1)%iPair*16.0_rkp + iTypeIni(1) + Half,    ',',&
                                              This%TrajFin%Molecules(1)%AngMom,                                 ',',&
                                              This%TrajFin%Molecules(1)%viba,                                   ',',& 
                                              This%TrajFin%Molecules(1)%iPair*16.0_rkp + iTypeFin(1) + Half
      flush(This%AnalysisOutputFile%Unit)

    else

      write(This%AnalysisOutputFile%Unit, 1)                                                                                            & 
                                              IdxOverall,                                                                           ',',&
                                              Traj%iPES(iTraj),                                                                     ',',& 
                                              bmax,                                                                                 ',',& 
                                              This%TrajIni%b,                                                                       ',',& 
                                              This%TrajIni%Molecules(1)%AngMom,                                                     ',',&
                                              This%TrajIni%Molecules(1)%viba,                                                       ',',&
                                              This%TrajIni%Molecules(2)%AngMom,                                                     ',',&
                                              This%TrajIni%Molecules(2)%viba,                                                       ',',&
                                              This%TrajIni%Molecules(1)%iPair*16.0_rkp + iTypeIni(1) + 4.0_rkp*iTypeIni(2) + Half,  ',',&
                                              This%TrajFin%Molecules(1)%AngMom,                                                     ',',&
                                              This%TrajFin%Molecules(1)%viba,                                                       ',',&
                                              This%TrajFin%Molecules(2)%AngMom,                                                     ',',&
                                              This%TrajFin%Molecules(2)%viba,                                                       ',',&
                                              This%TrajFin%Molecules(1)%iPair*16.0_rkp + iTypeFin(1) + 4.0_rkp*iTypeFin(2) + Half
      flush(This%AnalysisOutputFile%Unit)                       
    end if
    ! ==============================================================================================================

  end if


  !1 format(i8, 3x, i4, 3x, 2(es12.5,3x), *(es15.8,3x) )
  !2 format(i8, 3x, i4, 3x, 2(es12.5,3x), *(es15.8,3x) )
  1 format(I10, A, I5, 2(A, es11.5), *(A, es14.8) )
  


!  if (i_Debug_Loc) then
!    call CPU_Time( EndTime )
!    call Logger%Write( "Time in procedure ", EndTime - StartTime )
!  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine InitializePrinting( This, Input, Collision, Traj, NTrajBatch, i_Debug )
! This procedures controls the integration of trajectories

  use, intrinsic :: IEEE_ARITHMETIC
  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Parameters_Module           ,only:  Zero, One
  use PES_Class                   ,only:  PESEvoFile, PESEvoFlg
  use StateInitDiatomAtom_Module  ,only:  ParamsFile
  use Trajectories_Class          ,only:  Trajectories_Type
  use String_Module               ,only:  Convert_Ratio, Convert_To_Real, Convert_To_String
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  integer                                   ,intent(in)     ::    NTrajBatch
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, iTraj
  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializePrinting" )
  ! i_Debug_Loc   =     Logger%On()

  
  ! ==============================================================================================================
  !  OPENING AND WRITING HEADER FOR THE INITIAL PARAMETERS
  ! ==============================================================================================================
  ParamsFile%Active = Input%ParamsFlg
  if (ParamsFile%Active) then
    if (i_Debug_Loc) call Logger%Write( "Opening the Initial Parameters file (Params.out)? ", ParamsFile%Active )
    ParamsFile%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Params.out'
    if (i_Debug_Loc) call Logger%Write( "Initial Parameters file: ", ParamsFile%Name )
    open( NewUnit=ParamsFile%Unit, File=ParamsFile%Name, Action="WRITE", iostat=ParamsFile%Status )
    if (ParamsFile%Status /= 0) call Error( "Error opening file: " // ParamsFile%Name )
    ParamsFile%Format = "(i7,',',i4,*(',',es15.8:))"
    write(ParamsFile%Unit,'(A)') "# iTraj, iPES, Angle(1), Angle(2), Angle(3), ETran, Ekin, rBond, rdotBond, b, TarAngMom(1), TarAngMom(2), TarAngMom(3), ProAngMom(1), ProAngMom(2), ProAngMom(3), TotAngMom(1), TotAngMom(2), TotAngMom(3)"
  end if
  ! ==============================================================================================================
  
  
  ! ==============================================================================================================
  !     OPENING AND WRITING HEADER FOR THE PROGRESS FILE
  ! ==============================================================================================================
  This%ProgressFile%Active = Input%ProgressFlg
  if (This%ProgressFile%Active) then
    if (i_Debug_Loc) call Logger%Write( "Opening the progress file (Progress.out)? ", This%ProgressFile%Active )
    This%ProgressFile%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Progress.out'
    open( NewUnit=This%ProgressFile%Unit, File=This%ProgressFile%Name, Action='WRITE', Form='FORMATTED', iostat=This%ProgressFile%Status )
    if (This%ProgressFile%Status /= 0) call Error( "Error opening file: " // This%ProgressFile%Name )
    This%ProgressFile%Format = "(2x,i15,3x,f15.2,3(3x,es15.8))"
    write(This%ProgressFile%Unit,"('#',1x,5(a15,3x))") [ character(15) :: "Ntraj conv.", "Pencentage done", "Traj./second", "Time left [s]", "Time spend [s]"]
  end if
  ! ==============================================================================================================
  
  
  ! ==============================================================================================================
  !  OPENING AND WRITING HEADER FOR THE PAQs INITIAL AND FINAL CONDITIONS
  ! ==============================================================================================================
  This%PaQSolFile%Active = Input%PaQSolFlg
  if (This%PaQSolFile%Active) then
    if (i_Debug_Loc) call Logger%Write( "Opening the solution file (PaQSol.out)? ", This%PaQSolFile%Active )
    This%PaQSolFile%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/PaQSol.out'
    open( NewUnit=This%PaQSolFile%Unit, File=This%PaQSolFile%Name, Action="WRITE", iostat=This%PaQSolFile%Status )
    if (This%PaQSolFile%Status /= 0) call Error( "Error opening file: " // This%PaQSolFile%Name )
    This%PaQSolFile%Format = "(2x,i7,*(',',es15.8:))"
    write(This%PaQSolFile%Unit,"('# Traj. index, t_fin', *(',',A) )") [ character(12) :: "H_ini", ("PaQ_ini("//Convert_To_String(i)//")",i=1,Collision%NEqtTot), "H_fin", ("PaQ_fin("//Convert_To_String(i)//")",i=1,Collision%NEqtTot) ]
  end if
  ! ==============================================================================================================
  
  
  ! ==============================================================================================================
  !  OPENING AND WRITING HEADER FOR THE PES EVOLUTION
  ! ==============================================================================================================
  if (Input%PESEvoFlg) then
    PESEvoFlg = Input%PESEvoFlg
    allocate( PESEvoFile(NTrajBatch) )
    PESEvoFile(:)%Active = Input%PESEvoFlg
    do iTraj = 1,Traj%NTraj      
      if (i_Debug_Loc) call Logger%Write( "Opening the evolution file (PESEvo-*.out)? ", PESEvoFile(iTraj)%Active )
      if (PESEvoFile(iTraj)%Active) then
        PESEvoFile(iTraj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/PESEvo-'// Convert_To_String(iTraj) // '.out'
        if (i_Debug_Loc) call Logger%Write( "Evolution File for iTraj=", iTraj,": ", PESEvoFile(iTraj)%Name )
        open( NewUnit=PESEvoFile(iTraj)%Unit, File=PESEvoFile(iTraj)%Name, Action="WRITE", iostat=PESEvoFile(iTraj)%Status )
        if (PESEvoFile(iTraj)%Status /= 0) call Error( "Error opening file: " // PESEvoFile(iTraj)%Name )
        PESEvoFile(iTraj)%Format = "(2x,*(3x,es15.8:))"
        !write(PESEvoFile(iTraj)%Unit,"('#',3x,*(3x,a15:))") [ character(15) :: "time", ("R("//Convert_To_String(i)//")",i=1,Collision%NAtoms), "V", ("dV("//Convert_To_String(i)//")",i=1,Collision%NAtoms) ]
        write(PESEvoFile(iTraj)%Unit,"('#',3x,*(3x,a15:))") [ character(15) :: ("R("//Convert_To_String(i)//")",i=1,Collision%NAtoms), "V", ("dV("//Convert_To_String(i)//")",i=1,Collision%NAtoms) ]
      end if
    end do
  end if
  ! ==============================================================================================================

  
  if (Input%PaQEvoFlg) then
    allocate( This%PaQEvoFile(NTrajBatch) )
    This%PaQEvoFile(:)%Active = Input%PaQEvoFlg
  end if
  
  if (Input%XXEvoFlg) then
    allocate( This%XXEvoFile(NTrajBatch) )
    if (i_Debug_Loc) call Logger%Write( "Allocated ", NTrajBatch, " This%XXEvoFile" )
    This%XXEvoFile(:)%Active = Input%XXEvoFlg
    
    if (Input%XXEvoSnglFileFlg) then
      do iTraj = 1,NTrajBatch

        if (i_Debug_Loc) call Logger%Write( "Opening the evolution file (This%XXEvoFile-*.out)? ", This%XXEvoFile(iTraj)%Active )
        if (This%XXEvoFile(iTraj)%Active) then
          This%XXEvoFile(iTraj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/XXEvo-'// Convert_To_String(iTraj) // '.out'
          if (i_Debug_Loc) call Logger%Write( "Evolution File for iTraj=", iTraj, ": ", This%XXEvoFile(iTraj)%Name )
          
          open( NewUnit=This%XXEvoFile(iTraj)%Unit, File=This%XXEvoFile(iTraj)%Name, Action="WRITE", iostat=This%XXEvoFile(iTraj)%Status )
          write(*,*) This%XXEvoFile(iTraj)%Unit, This%XXEvoFile(iTraj)%Status 
          if (This%XXEvoFile(iTraj)%Status /= 0) call Error( "Error opening file: " // This%XXEvoFile(iTraj)%Name )
          if (Collision%NAtoms.eq.3) then 
             This%XXEvoFile(iTraj)%Format = "(3x, es20.8, 18(',', es20.8), ',', I10 )"
             write(This%XXEvoFile(iTraj)%Unit,"('#', 2x, a20, 18(',', a20:), ',', a10)") 'Time [s]', 'x(1)', 'y(1)', 'z(1)', 'x(2)', 'y(2)', 'z(2)', 'x(3)', 'y(3)', 'z(3)', 'vx(1)', 'vy(1)', 'vz(1)', 'vx(2)', 'vy(2)', 'vz(2)', 'vx(3)', 'vy(3)', 'vz(3)', 'iP'
          else if (Collision%Natoms.eq.4) then
             This%XXEvoFile(iTraj)%Format = "(3x, es20.8, 24(',', es20.8), ',', I10, I10 )"
             write(This%XXEvoFile(iTraj)%Unit,"('#', 2x, a20, 24(',', a20:), ',', a10)") 'Time [s]', 'x(1)', 'y(1)', 'z(1)', 'x(2)', 'y(2)', 'z(2)', 'x(3)', 'y(3)', 'z(3)', 'x(4)', 'y(4)', 'z(4)', 'vx(1)', 'vy(1)', 'vz(1)', 'vx(2)', 'vy(2)', 'vz(2)', 'vx(3)', 'vy(3)', 'vz(3)', 'vx(4)', 'vy(4)', 'vz(4)', 'iP', 'jP'
          endif
        endif 
      end do
    end if

  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine InitializePrintingEvo( This, Input, Collision, Traj, TrajTemp, iTraj, i_Debug )
! This procedures controls the integration of trajectories

  use, intrinsic :: IEEE_ARITHMETIC
  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Parameters_Module           ,only:  Zero, One
  use String_Module               ,only:  Convert_Ratio, Convert_To_Real, Convert_To_String
  use Trajectories_Class          ,only:  Trajectories_Type
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  type(Trajectories_Type)                   ,intent(inout)  ::    TrajTemp
  integer                                   ,intent(in)     ::    iTraj
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, iP, iTrajj
  real(rkp) ,dimension(NSpace,Collision%NAtoms)             ::    xx          ! Dim=(NSpace,NAtoms)
  real(rkp) ,dimension(NSpace,Collision%NAtoms)             ::    xxdot       ! Dim=(NSpace,NAtoms)
  logical                                                   ::    i_Debug_Loc
  integer   ,dimension(3,4)                                 ::    iPToAtom        


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializePrintingEvo" )
  ! i_Debug_Loc   =     Logger%On()

    
  iPToAtom(1,:) = [1,2,3,4]
  iPToAtom(2,:) = [1,3,2,4]
  iPToAtom(3,:) = [2,3,1,4]


  TrajTemp              = Traj
  TrajTemp%PaQ(:,iTraj) = Traj%PaQ0(:,iTraj)
  call This%TrajTemp%SetData( Traj%t(iTraj), Traj%H(iTraj), TrajTemp%PaQ(:,iTraj) )


  ! ==============================================================================================================
  !  OPENING AND WRITING HEADER FOR THE PAQ EVOLUTION
  ! ==============================================================================================================
  if (Input%PaQEvoFlg) then
    if (This%PaQEvoFile(iTraj)%Active) then
      if (i_Debug_Loc) call Logger%Write( "Opening the evolution file (PaQEvo-*.out)? ", This%PaQEvoFile(iTraj)%Active )
      This%PaQEvoFile(iTraj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/PaQEvo-'// Convert_To_String(iTraj) // '.out'
      if (i_Debug_Loc) call Logger%Write( "Evolution File for iTraj=", iTraj,": ", This%PaQEvoFile(iTraj)%Name )
      open( NewUnit=This%PaQEvoFile(iTraj)%Unit, File=This%PaQEvoFile(iTraj)%Name, Action="WRITE", iostat=This%PaQEvoFile(iTraj)%Status )
      if (This%PaQEvoFile(iTraj)%Status /= 0) call Error( "Error opening file: " // This%PaQEvoFile(iTraj)%Name )
      This%PaQEvoFile(iTraj)%Format = "(2x,*(3x,es15.8:))"
      write(This%PaQEvoFile(iTraj)%Unit,"('#',3x,*(3x,a15:))") [ character(15) :: "time", "H", ("PaQ("//Convert_To_String(i)//")",i=1,Collision%NEqtTot) ]
      
      if (.not. This%CMFrameFlg) then
        call This%TrajTemp%SetData( Traj%t(iTraj), Traj%H(iTraj), TrajTemp%PaQ(:,iTraj) )
        TrajTemp%PaQ(:,iTraj) = Collision%ApplyAntiTransformation( iTraj, TrajTemp )
      end if

      call This%TrajTemp%SetData( Traj%t(iTraj), Traj%H(iTraj), TrajTemp%PaQ(:,iTraj) )
      write(This%PaQEvoFile(iTraj)%Unit,This%PaQEvoFile(iTraj)%Format) Zero, Traj%H0(iTraj), TrajTemp%PaQ(:,iTraj)
    end if
  end if
  ! ==============================================================================================================
        
        
  ! ==============================================================================================================
  !  OPENING AND WRITING HEADER FOR THE COORDINATES EVOLUTION
  ! ==============================================================================================================
  if (Input%XXEvoFlg) then
    if (This%XXEvoFile(iTraj)%Active) then
      
      call ComputeCoordAndVeloc( This%TrajTemp, Collision, xx, xxdot )   

      if (Input%XXEvoSnglFileFlg) then
        if (i_Debug_Loc) call Logger%Write( "Opening the evolution file (XXEvo-*.out)? ", This%XXEvoFile(iTraj)%Active )

        if (size(xx,2).eq.3) then

           write(This%XXEvoFile(Traj%Idx(iTraj))%Unit, This%XXEvoFile(Traj%Idx(iTraj))%Format) Traj%t(iTraj), &
                                               xx(1,iPToAtom(1,1)), xx(2,iPToAtom(1,1)), xx(3,iPToAtom(1,1)), &
                                               xx(1,iPToAtom(1,2)), xx(2,iPToAtom(1,2)), xx(3,iPToAtom(1,2)), &
                                               xx(1,iPToAtom(1,3)), xx(2,iPToAtom(1,3)), xx(3,iPToAtom(1,3)), &
                                               xxdot(1,iPToAtom(1,1)), xxdot(2,iPToAtom(1,1)), xxdot(3,iPToAtom(1,1)), &
                                               xxdot(1,iPToAtom(1,2)), xxdot(2,iPToAtom(1,2)), xxdot(3,iPToAtom(1,2)), &
                                               xxdot(1,iPToAtom(1,3)), xxdot(2,iPToAtom(1,3)), xxdot(3,iPToAtom(1,3)), 1

        else 

           write(This%XXEvoFile(Traj%Idx(iTraj))%Unit, This%XXEvoFile(Traj%Idx(iTraj))%Format) Traj%t(iTraj), &
                                               xx(1,iPToAtom(1,1)), xx(2,iPToAtom(1,1)), xx(3,iPToAtom(1,1)), &
                                               xx(1,iPToAtom(1,2)), xx(2,iPToAtom(1,2)), xx(3,iPToAtom(1,2)), &
                                               xx(1,iPToAtom(1,3)), xx(2,iPToAtom(1,3)), xx(3,iPToAtom(1,3)), &
                                               xx(1,iPToAtom(1,4)), xx(2,iPToAtom(1,4)), xx(3,iPToAtom(1,4)), &
                                               xxdot(1,iPToAtom(1,1)), xxdot(2,iPToAtom(1,1)), xxdot(3,iPToAtom(1,1)), &
                                               xxdot(1,iPToAtom(1,2)), xxdot(2,iPToAtom(1,2)), xxdot(3,iPToAtom(1,2)), &
                                               xxdot(1,iPToAtom(1,3)), xxdot(2,iPToAtom(1,3)), xxdot(3,iPToAtom(1,3)), & 
                                               xxdot(1,iPToAtom(1,4)), xxdot(2,iPToAtom(1,4)), xxdot(3,iPToAtom(1,4)), 1, 2

        endif

      else
        if (i_Debug_Loc) call Logger%Write( "Opening the evolution file (XXEvo-*.out)? ", This%XXEvoFile(iTraj)%Active )
      
        iTrajj = Traj%Idx(iTraj)

        if (size(xx,2).eq.3) then
        
          iP = 1
          call system('mkdir -p '   // trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Traj-'// Convert_To_String(iTraj) )
          This%XXEvoFile(iTrajj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Traj-'// Convert_To_String(iTraj) // &
                                       '/XXEvo.vtk.' // Convert_To_String(int(Traj%t(iTraj)))
          open( NewUnit=This%XXEvoFile(iTrajj)%Unit, File=This%XXEvoFile(iTrajj)%Name, Action="WRITE", iostat=This%XXEvoFile(iTrajj)%Status )
            if (This%XXEvoFile(iTrajj)%Status /= 0) call Error( "Error opening file: " // This%XXEvoFile(iTrajj)%Name )
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "# vtk DataFile Version 3.0"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "vtk output"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "ASCII"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "DATASET POLYDATA"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "FIELD FieldData 1"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "TIME 1 1 double"
            write(This%XXEvoFile(iTrajj)%Unit,"(1es15.8)")     Traj%t(iTraj)
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINTS 3 double"
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,1)), xx(2,iPToAtom(iP,1)), xx(3,iPToAtom(iP,1))
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,2)), xx(2,iPToAtom(iP,2)), xx(3,iPToAtom(iP,2))
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,3)), xx(2,iPToAtom(iP,3)), xx(3,iPToAtom(iP,3))
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "LINES 2 12"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 0 1"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 1 0"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 2 2"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 2 2"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINT_DATA 3"
          close(This%XXEvoFile(iTrajj)%Unit)

        else

          iP = 1
          call system('mkdir -p '   // trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Traj-'// Convert_To_String(iTraj) )
          This%XXEvoFile(iTrajj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Traj-'// Convert_To_String(iTraj) // &
                                       '/XXEvo.vtk.' // Convert_To_String(int(Traj%t(iTraj)))
          open( NewUnit=This%XXEvoFile(iTrajj)%Unit, File=This%XXEvoFile(iTrajj)%Name, Action="WRITE", iostat=This%XXEvoFile(iTrajj)%Status )
            if (This%XXEvoFile(iTrajj)%Status /= 0) call Error( "Error opening file: " // This%XXEvoFile(iTrajj)%Name )
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "# vtk DataFile Version 3.0"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "vtk output"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "ASCII"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "DATASET POLYDATA"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "FIELD FieldData 1"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "TIME 1 1 double"
            write(This%XXEvoFile(iTrajj)%Unit,"(1es15.8)")     Traj%t(iTraj)
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINTS 4 double"
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,1)), xx(2,iPToAtom(iP,1)), xx(3,iPToAtom(iP,1))
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,2)), xx(2,iPToAtom(iP,2)), xx(3,iPToAtom(iP,2))
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,3)), xx(2,iPToAtom(iP,3)), xx(3,iPToAtom(iP,3))
            write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,4)), xx(2,iPToAtom(iP,4)), xx(3,iPToAtom(iP,4))
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "LINES 2 12"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 0 1"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 1 0"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 2 3"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 3 2"
            write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINT_DATA 4"
          close(This%XXEvoFile(iTrajj)%Unit)

        end if
        
      end if
      
!      This%XXEvoFile(iTraj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/XXEvo-'// Convert_To_String(iTraj) // '.out'
!      if (i_Debug_Loc) call Logger%Write( "Evolution File for iTraj=", iTraj,": ", This%XXEvoFile(iTraj)%Name )
!      open( NewUnit=This%XXEvoFile(iTraj)%Unit, File=This%XXEvoFile(iTraj)%Name, Action="WRITE", iostat=This%XXEvoFile(iTraj)%Status )
!      if (This%XXEvoFile(iTraj)%Status /= 0) call Error( "Error opening file: " // This%XXEvoFile(iTraj)%Name )
!      This%XXEvoFile(iTraj)%Format = "(2x,*(3x,es15.8:))"
!      write(This%XXEvoFile(iTraj)%Unit,"('#',3x,*(3x,a15:))") [ character(15) :: "time", "x(1,1)",  "x(2,1)",  "x(3,1)",  "x(1,2)",  "x(2,2)",  "x(3,2)",  "x(1,3)",  "x(2,3)",  "x(3,3)", &
!                                                                                               "dx(1,1)", "dx(2,1)", "dx(3,1)", "dx(1,2)", "dx(2,2)", "dx(3,2)", "dx(1,3)", "dx(2,3)", "dx(3,3)"  ]
!      TrajTemp%PaQ(:,Traj%Idx(iTraj)) = Collision%ApplyTransformation( Traj%Idx(iTraj), Traj )
!      call This%TrajIni%SetData( Traj%t(Traj%Idx(iTraj)), Traj%H0(Traj%Idx(iTraj)), TrajTemp%PaQ(:,Traj%Idx(iTraj)) )
!      call ComputeCoordAndVeloc( This%TrajIni, Collision, xx, xxdot )
!      write(This%XXEvoFile(Traj%Idx(iTraj))%Unit,This%XXEvoFile(Traj%Idx(iTraj))%Format)  Traj%t(iTraj), xx(1:3,1),    xx(1:3,2),    xx(1:3,3), xxdot(1:3,1), xxdot(1:3,2), xxdot(1:3,3)  
    end if
  end if
  ! ==============================================================================================================
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PrintingParams( This, i_Debug )
! This procedures controls the integration of trajectories

  use, intrinsic :: IEEE_ARITHMETIC
  use Input_Class                 ,only:  Input_Type
  use Parameters_Module           ,only:  Zero, One
  use StateInitDiatomAtom_Module  ,only:  ParamsFile, Params
  
  class(Integrator_Type)                    ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i, iTraj
  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PrintingParams" )
  ! i_Debug_Loc   =     Logger%On()

  ! ==============================================================================================================
  !  WRITING THE INITIAL PARAMETERS
  ! ==============================================================================================================
  if (ParamsFile%Active) then
    write(ParamsFile%Unit,ParamsFile%Format) int(Params(1)), int(Params(2)), Params(3:19)
  end if
  ! ==============================================================================================================
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PrintingEvo( This, Input, Collision, Traj, TrajTemp, Converged, i_Debug )
! This procedures controls the integration of trajectories

  use, intrinsic :: IEEE_ARITHMETIC
  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Parameters_Module           ,only:  Zero, One
  use Trajectories_Class          ,only:  Trajectories_Type
  use String_Module               ,only:  Convert_Ratio, Convert_To_Real, Convert_To_String
  use Sorting_Module              ,only:  hpsort

  class(Integrator_Type)                    ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  type(Trajectories_Type)                   ,intent(in)     ::    Traj
  type(Trajectories_Type)                   ,intent(inout)  ::    TrajTemp
  logical   ,dimension(:)                   ,intent(in)     ::    Converged
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iTraj
  real(rkp) ,dimension(NSpace,Collision%NAtoms)             ::    xx          ! Dim=(NSpace,NAtoms)
  real(rkp) ,dimension(NSpace,Collision%NAtoms)             ::    xxdot       ! Dim=(NSpace,NAtoms)
  logical                                                   ::    i_Debug_Loc
  real(rkp) ,dimension(3)                                   ::    Rbs, Rbs_                   
  integer   ,dimension(3)                                   ::    iRbs 
  integer                                                   ::    i, j, iP, jP, iTrajj, NEqtVar 
  integer   ,dimension(3,4)                                 ::    iPToAtom        
  real(rkp) ,dimension(12)                                  ::    PaQ 

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PrintingEvo" )
  ! i_Debug_Loc   =     Logger%On()
  
  iPToAtom(1,:) = [1,2,3,4]
  iPToAtom(2,:) = [1,3,2,4]
  iPToAtom(3,:) = [2,3,1,4]
  

  do iTraj = 1,Traj%NTraj
    if (.not. Converged(iTraj)) then
      
      if (This%CMFrameFlg) then
        TrajTemp%PaQ(:,iTraj) = Collision%ApplyTransformation( iTraj, Traj )
      else
        TrajTemp%PaQ(:,iTraj) = Traj%PaQ(:,iTraj)
      end if
      call This%TrajTemp%SetData( Traj%t(iTraj), Traj%H(iTraj), TrajTemp%PaQ(:,iTraj) )

      if (Input%PaQEvoFlg) then
        if (This%PaQEvoFile(Traj%Idx(iTraj))%Active) then
          write(This%PaQEvoFile(Traj%Idx(iTraj))%Unit,This%PaQEvoFile(Traj%Idx(iTraj))%Format)  Traj%t(iTraj), Traj%H(iTraj), TrajTemp%PaQ(:,iTraj)
        end if
      end if

      
      if (Input%XXEvoFlg) then
        if (This%XXEvoFile(Traj%Idx(iTraj))%Active) then
        
          call ComputeCoordAndVeloc( This%TrajTemp, Collision, xx, xxdot )          
          
          iTrajj = Traj%Idx(iTraj)
          
          iP = 0
          do i = 1,2
            do j = i+1,3
              iP      = iP + 1
              Rbs(iP) = ( xx(1,i) - xx(1,j) )**2 + ( xx(2,i) - xx(2,j) )**2 + ( xx(3,i) - xx(3,j) )**2
            end do
          end do
          Rbs_ = Rbs
          call hpsort( Rbs, iRbs )
          iP = iRbs(1)

          if (Input%XXEvoSnglFileFlg) then
            if (i_Debug_Loc) call Logger%Write( "Opening the evolution file (XXEvo-*.out)? ", This%XXEvoFile(Traj%Idx(iTraj))%Active )

            if (size(xx,2).eq.3) then

               jP = 1
               write(This%XXEvoFile(Traj%Idx(iTraj))%Unit, This%XXEvoFile(Traj%Idx(iTraj))%Format) Traj%t(iTraj), &
                                                  xx(1,iPToAtom(jP,1)),    xx(2,iPToAtom(jP,1)),    xx(3,iPToAtom(jP,1)), &
                                                  xx(1,iPToAtom(jP,2)),    xx(2,iPToAtom(jP,2)),    xx(3,iPToAtom(jP,2)), &
                                                  xx(1,iPToAtom(jP,3)),    xx(2,iPToAtom(jP,3)),    xx(3,iPToAtom(jP,3)), &
                                               xxdot(1,iPToAtom(jP,1)), xxdot(2,iPToAtom(jP,1)), xxdot(3,iPToAtom(jP,1)), &
                                               xxdot(1,iPToAtom(jP,2)), xxdot(2,iPToAtom(jP,2)), xxdot(3,iPToAtom(jP,2)), &
                                               xxdot(1,iPToAtom(jP,3)), xxdot(2,iPToAtom(jP,3)), xxdot(3,iPToAtom(jP,3)), iP

            else 

                jP = iRbs(2)

                write(This%XXEvoFile(Traj%Idx(iTraj))%Unit, This%XXEvoFile(Traj%Idx(iTraj))%Format) Traj%t(iTraj), &
                                               xx(1,iPToAtom(1,1)), xx(2,iPToAtom(1,1)), xx(3,iPToAtom(1,1)), &
                                               xx(1,iPToAtom(1,2)), xx(2,iPToAtom(1,2)), xx(3,iPToAtom(1,2)), &
                                               xx(1,iPToAtom(1,3)), xx(2,iPToAtom(1,3)), xx(3,iPToAtom(1,3)), &
                                               xx(1,iPToAtom(1,4)), xx(2,iPToAtom(1,4)), xx(3,iPToAtom(1,4)), &
                                               xxdot(1,iPToAtom(1,1)), xxdot(2,iPToAtom(1,1)), xxdot(3,iPToAtom(1,1)), &
                                               xxdot(1,iPToAtom(1,2)), xxdot(2,iPToAtom(1,2)), xxdot(3,iPToAtom(1,2)), &
                                               xxdot(1,iPToAtom(1,3)), xxdot(2,iPToAtom(1,3)), xxdot(3,iPToAtom(1,3)), & 
                                               xxdot(1,iPToAtom(1,4)), xxdot(2,iPToAtom(1,4)), xxdot(3,iPToAtom(1,4)), iP, jP

            endif


            !write(This%XXEvoFile(Traj%Idx(iTraj))%Unit, This%XXEvoFile(Traj%Idx(iTraj))%Format) Traj%t(iTraj), &
            !                                                           xx(1,iPToAtom(1,1)), xx(2,iPToAtom(1,1)), xx(3,iPToAtom(1,1)), &
            !                                                           xx(1,iPToAtom(1,2)), xx(2,iPToAtom(1,2)), xx(3,iPToAtom(1,2)), &
            !                                                           xx(1,iPToAtom(1,3)), xx(2,iPToAtom(1,3)), xx(3,iPToAtom(1,3)), &
            !                                                           xxdot(1,iPToAtom(1,1)), xxdot(2,iPToAtom(1,1)), xxdot(3,iPToAtom(1,1)), &
            !                                                           xxdot(1,iPToAtom(1,2)), xxdot(2,iPToAtom(1,2)), xxdot(3,iPToAtom(1,2)), &
            !                                                           xxdot(1,iPToAtom(1,3)), xxdot(2,iPToAtom(1,3)), xxdot(3,iPToAtom(1,3)), iP

          else

            This%XXEvoFile(iTrajj)%Name = trim(adjustl(Input%LevelOutputDir)) // '/Node_' // trim(adjustl(Input%iNode_char)) // '/Proc_' // trim(adjustl(Input%iProc_char)) // '/Traj-'// Convert_To_String(Traj%Idx(iTraj)) // &
                                          '/XXEvo.vtk.' // Convert_To_String(int(Traj%t(iTraj)))

            if (size(xx,2).eq.3) then

              open( NewUnit=This%XXEvoFile(iTrajj)%Unit, File=This%XXEvoFile(iTrajj)%Name, Action="WRITE", iostat=This%XXEvoFile(iTrajj)%Status )
                if (This%XXEvoFile(iTrajj)%Status /= 0) call Error( "Error opening file: " // This%XXEvoFile(iTrajj)%Name )
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "# vtk DataFile Version 3.0"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "vtk output"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "ASCII"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "DATASET POLYDATA"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "FIELD FieldData 1"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "TIME 1 1 double"
                write(This%XXEvoFile(iTrajj)%Unit,"(1es20.8)")     Traj%t(iTraj)
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINTS 3 double"
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,1)), xx(2,iPToAtom(iP,1)), xx(3,iPToAtom(iP,1))
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,2)), xx(2,iPToAtom(iP,2)), xx(3,iPToAtom(iP,2))
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(iP,3)), xx(2,iPToAtom(iP,3)), xx(3,iPToAtom(iP,3))
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "LINES 2 12"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 0 1"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 1 0"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 2 2"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 2 2"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINT_DATA 3"
              close(This%XXEvoFile(iTrajj)%Unit)      

            else

              open( NewUnit=This%XXEvoFile(iTrajj)%Unit, File=This%XXEvoFile(iTrajj)%Name, Action="WRITE", iostat=This%XXEvoFile(iTrajj)%Status )
                if (This%XXEvoFile(iTrajj)%Status /= 0) call Error( "Error opening file: " // This%XXEvoFile(iTrajj)%Name )
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "# vtk DataFile Version 3.0"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "vtk output"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "ASCII"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "DATASET POLYDATA"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "FIELD FieldData 1"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "TIME 1 1 double"
                write(This%XXEvoFile(iTrajj)%Unit,"(1es20.8)")     Traj%t(iTraj)
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINTS 4 double"
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(1,1)), xx(2,iPToAtom(1,1)), xx(3,iPToAtom(1,1))
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(1,2)), xx(2,iPToAtom(1,2)), xx(3,iPToAtom(1,2))
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(1,3)), xx(2,iPToAtom(1,3)), xx(3,iPToAtom(1,3))
                write(This%XXEvoFile(iTrajj)%Unit,"(3es20.8)") xx(1,iPToAtom(1,4)), xx(2,iPToAtom(1,4)), xx(3,iPToAtom(1,4))
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "LINES 2 12"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 0 1"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 1 0"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 2 3"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "2 3 2"
                write(This%XXEvoFile(iTrajj)%Unit,'(A)')           "POINT_DATA 4"
              close(This%XXEvoFile(iTrajj)%Unit)      

            end if

          end if
          !write(This%XXEvoFile(Traj%Idx(iTraj))%Unit,This%XXEvoFile(Traj%Idx(iTraj))%Format)  Traj%t(iTraj), xx(1:3,1), xx(1:3,2), xx(1:3,3), xxdot(1:3,1), xxdot(1:3,2), xxdot(1:3,3)
        end if
      end if
    
    end if
  end do
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
