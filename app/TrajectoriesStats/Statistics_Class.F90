! -*-F90-*-
!===============================================================================================================
! 
! Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
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

Module Statistics_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use File_Class            ,only:  File_Type

  implicit none

  private
  public    ::    Statistics_Type

  Type  ::    Statistics_Type                                       
    character(:)                    ,allocatable  ::    AssignmentMethod                                                              !< Method for States Assigment.
    integer(rkp)                                  ::    NTrajectoriesToAnalyze                                                        !< Maximum Nb of Trajectories to Analyze. -1 means Analyze all Trajectories on File.
    logical                                       ::    IdenticalDiatoms                                                              !< ???
    integer(rkp)                                  ::    NCond                                                                         !< Number of Overall Conditions (j/v/arrangment) on the Trajectories.
    integer(rkp)    ,dimension(:)   ,allocatable  ::    iskip                                         !Dim[NCond]                     !< Vector of Flags (0/1) for skipping the Conditions (j/v/arrangment) on the Trajectories.
    logical         ,dimension(:)   ,allocatable  ::    PresEvOdd                                     !Dim[NCond]                     !< True if Condition corresponds to Rotational Q.N. that has to keep oddness / evenness. 
    integer(rkp)                                  ::    QuantumNumberMax = 300                                                        !< Maximum Value Allowed for QNs. Above, an Error is Raised. 
    integer(rkp)                                  ::    NTraj                                                                         !< Nb of Trajectories.
    integer(rkp)                                  ::    NRings                                                                        !< Nb of Impact Parameters' Rings.
    real(rkp)       ,dimension(:)   ,allocatable  ::    bSampled                                      !Dim[NTraj]                     !< Vector of Impact Parameters Sampled at the Begininning of Trajectories.
    real(rkp)       ,dimension(:)   ,allocatable  ::    bMax                                          !Dim[NRings] -> Dim[<= NRings]  !< Vector of Impact Parameters' Maxima for Sampling at the Begininning of Trajectories.
    real(rkp)       ,dimension(:,:) ,allocatable  ::    Qini                                          !Dim[NCond,NTraj]               !< Set of Overall Initial Conditions.
    real(rkp)       ,dimension(:,:) ,allocatable  ::    Qfin                                          !Dim[NCond,NTraj]               !< Set of Overall Fin Conditions.
    integer(rkp)    ,dimension(:)   ,allocatable  ::    IniStateCode                                  !Dim[NTraj]                     !< Vector of Identification Nbs for the Trajectories' Initial States.
    integer(rkp)    ,dimension(:)   ,allocatable  ::    SortedIndx_IniStateCode                       !Dim[NTraj]                     !< Sorted Indx for the Vector of Identification Nbs.
    real(rkp)       ,dimension(:)   ,allocatable  ::    RingArea                                      !Dim[<= NRings]                 !< Vector of Impact Parameters' Rings Areas.
    integer(rkp)                                  ::    ifact                                                                         !< Temporary Integer used for computing the Trajectory's IniStateCode.
    integer                                       ::    PESoI                                                                         !< Potential Energy Surfaces (PES) of Interest
    character(6)                                  ::    PESoI_char
    type(File_Type)                               ::    ResidOutFile
    type(File_Type)                               ::    ProbaOutFile
    type(File_Type)                               ::    TrajeOutFile
    type(File_Type)                               ::    StatOutFile
    type(File_Type)                               ::    bSensitivityFile
    logical                                       ::    StatReadsBinaryFlg  = .False.
    logical                                       ::    StatWritesBinaryFlg = .False.
  contains
    private
    procedure ,public   ::    Initialize    =>    InitializeStatistics
    procedure ,public   ::    Process       =>    ProcessStatistics
    procedure           ::    AddFinState
    procedure           ::    ReadInputs
    procedure           ::    ReadInputsUnformatted
    procedure           ::    SetRings
    procedure           ::    IdentifyInitialStates
    procedure           ::    PrepareOutputFiles
    procedure           ::    WriteFinStateProbabilities
  End Type
  
  logical   ,parameter  ::    i_Debug_Global = .True.

  contains

! ***********************************************************************************************************************************************!
Subroutine InitializeStatistics( This, Input, i_Debug, i_Debug_Deep )

  use Input_Class     ,only:  Input_Type

  class(Statistics_Type)                                ,intent(out)    ::    This
  type(Input_Type)                                      ,intent(in)     ::    Input
  logical                                     ,optional ,intent(in)     ::    i_Debug
  logical                                     ,optional ,intent(in)     ::    i_Debug_Deep

  logical                                                               ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeStatistics" )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!   EXTRACTING DATA FROM THE INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Extracting data from the input object" )
  This%NCond                  =   Input%NCond
  This%AssignmentMethod       =   Input%AssignmentMethod
  allocate( This%PresEvOdd(size(Input%PresEvOdd)), source = Input%PresEvOdd )
  This%NTrajectoriesToAnalyze =   Input%NTrajectoriesToAnalyze
  This%IdenticalDiatoms       =   Input%IdenticalDiatoms
  allocate( This%iskip(size(Input%iskip)), source = int(Input%iskip,8) )
  if (i_Debug_Loc) then
    call Logger%Write( "-> Indicator of Identical Diatoms:                 This%IdenticalDiatoms = ", This%IdenticalDiatoms )
    call Logger%Write( "-> Number of Overall Conditions on Trajs:          This%NCond            = ", This%NCond )
    call Logger%Write( "-> Vector of Skipping Flgs for Overall Conditions: This%iskip            = ", This%iskip )
    call Logger%Write( "-> Maximum value allowed for Q.N.s:                This%QuantumNumberMax = ", This%QuantumNumberMax )
    if (adjustl(trim(This%AssignmentMethod)) == 'Histogram') then 
      call Logger%Write( "Method for states assigment = Histogram method" )
    elseif (adjustl(trim(This%AssignmentMethod)) == 'QSS') then 
      call Logger%Write( "Method for states assigment = Quadratic smooth sampling method" )
    else
      call Error( "Error: Unknown assignment method for Fin states" )
    end if
  end if

  This%PESoI      = Input%PESoI 
  This%PESoI_char = Input%PESoI_char
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE IMPACT PARAMETERS RINGS
! ==============================================================================================================
  This%StatReadsBinaryFlg  = Input%StatReadsBinaryFlg
  This%StatWritesBinaryFlg = Input%StatWritesBinaryFlg

  if (i_Debug_Loc) call Logger%Write( "Reading the statistings input data" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling This%ReadInputs" )
  if (This%StatReadsBinaryFlg) then
    call This%ReadInputsUnformatted(i_Debug_Loc)
  else
    call This%ReadInputs(i_Debug_Loc)
  end if
  if (i_Debug_Loc) call Logger%Write( "-> Done reading the statistics input data" )
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE IMPACT PARAMETERS RINGS
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the impact parameters rings" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling This%SetRings" )
  call This%SetRings(i_Debug_Loc)
  if (i_Debug_Loc) call Logger%Write( "-> Done setting the impact parameters rings" )
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE IMPACT PARAMETERS RINGS
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Identifying and sorting initial states" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling This%IdentifyInitialStates" )
  call This%IdentifyInitialStates( i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Deep )
  if (i_Debug_Loc) call Logger%Write( "-> Done identifying and sorting initial states" )
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE IMPACT PARAMETERS RINGS
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Preparing the output files" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling This%PrepareOutputFiles" )
  call This%PrepareOutputFiles(i_Debug_Loc)
  if (i_Debug_Loc) call Logger%Write( "-> Done preparing the output files" )
! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
Subroutine ReadInputs( This, i_Debug )

  use, intrinsic :: iso_fortran_env ,only:  IOStat_End

  class(Statistics_Type)                                ,intent(inout)  ::    This
  logical                                       ,optional ,intent(in)   ::    i_Debug

  logical                                                               ::    i_Debug_Loc
  integer(rkp)                                                          ::    iTraj, Idx, iPES
  logical                                                               ::    Limit                                               ! TRUE if There is a limit on the Nb of Trajs to Analyze
  Integer                                                               ::    iCond
  Integer                                                               ::    UnitWrite, StatusWrite
  type(File_Type)                                                       ::    DataFile
  integer                                                               ::    iTrajExcluded
  character(10)                                                         ::    iTrajExcluded_Char
  integer(rkp)                                                          ::    iTemp1, iTemp2
  real(rkp)                                                             ::    Temp1, Temp2
  real(rkp)             ,dimension(:) ,allocatable                      ::    TempVec1, TempVec2 

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadInputs" )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!     OPENING AND WRITING HEADER FOR THE PROGRESS FILE
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Opening the data file" )
  if (This%PESoI == 0) then
    !DataFile%Name     =   'trajectories.out'                                                                                        !# FOR COMPATIBILITY WITH CG-QCT CODE
    DataFile%Name     =   'trajectories.csv'                                                                                         !#
  else
    !DataFile%Name     =   trim(adjustl( 'trajectories.out.' // trim(adjustl(This%PESoI_char)) ))                                     !# FOR COMPATIBILITY WITH CG-QCT CODE
    DataFile%Name     =   trim(adjustl( 'trajectories.csv.' // trim(adjustl(This%PESoI_char)) ))                                     !#
  end if
  open( NewUnit=DataFile%Unit, File=DataFile%Name, Action='READ', Form='FORMATTED', iostat=DataFile%Status )
  if (DataFile%Status/=0) call Error( "Error opening file: " // DataFile%Name )
  DataFile%Format   =   "( i9,3x, 2(es15.8,3x), 2(4(i3,3x)))"
  read(DataFile%Unit,*)
! ==============================================================================================================


! ==============================================================================================================
!     SETTING THE NUMBER OF TRAJECTORIES TO BE ANALYZED
! ==============================================================================================================
  allocate( TempVec1(This%NCond) )
  allocate( TempVec2(This%NCond) )

  if (i_Debug_Loc) call Logger%Write( "Setting the number of trajectories to be analyzed" )
  Limit         =   This%NTrajectoriesToAnalyze > 0 
  iTraj         =   0
  iTrajExcluded =   0
  do
    read(DataFile%Unit,*,iostat=DataFile%Status) iTemp1, iTemp2, Temp1, Temp2, TempVec1, TempVec2 
    if ( DataFile%Status == IOStat_End ) exit                                                                                     
    !if (DataFile%Status/=0) call Error( "Error reading the data file for statistics: " // DataFile%Name  )                        
    if (DataFile%Status/=0) then
      iTrajExcluded = iTrajExcluded + 1
    else
      iTraj         = iTraj + 1                                                                                  
    end if
    if ( Limit .and. iTraj > This%NTrajectoriesToAnalyze ) exit                                                 
  end do
  This%NTraj = iTraj
  if (i_Debug_Loc) call Logger%Write( "-> Nb of trajectories: This%NTraj = ", This%NTraj )
  if (iTrajExcluded>0) then
    write(iTrajExcluded_Char, '(I10)') iTrajExcluded
    if (i_Debug_Loc) call Logger%Write( "-> Nb of trajectories Excluded: iTrajExcluded = ", iTrajExcluded )
    write(*,*) "          [Statistics_Class.F90]: Nb of Excluded Trajectories " // adjustl(trim(iTrajExcluded_Char))
  end if
! ==============================================================================================================


! ==============================================================================================================
!     ALLOCATE THE DATA
! ==============================================================================================================
  allocate( This%bMax(This%NTraj) )
  allocate( This%bSampled(This%NTraj) )

  allocate( This%Qini(This%NCond,This%NTraj) )
  allocate( This%Qfin(This%NCond,This%NTraj) )
! ==============================================================================================================


! ==============================================================================================================
!   IDENTIFYING AND SORTING INITIAL STATES AND IMPACT PARAMETER RINGS
! ==============================================================================================================
  if (This%StatWritesBinaryFlg) then
    open( NewUnit=UnitWrite, File='./trajectories.bin', Action='WRITE', access="Stream", form="Unformatted", iostat=StatusWrite )
    if (StatusWrite/=0) call Error( "Error writing the binary data file for statistics: " // './trajectories.bin'  ) 
      write(UnitWrite) int(This%NTraj, rkp)
  end if

  if (i_Debug_Loc) call Logger%Write( "Reading the trajectory data: bMax, bSampled, Qini, Qfin" )
  rewind(DataFile%Unit)                                                                                           
  read(DataFile%Unit,*)                                                                                     
  iTraj=0
  do                                                                                     
    read(DataFile%Unit,*,iostat=DataFile%Status) iTemp1, iTemp2, Temp1, Temp2, TempVec1, TempVec2 
    !if (DataFile%Status/=0) call Error( "Error reading the data file for statistics: " // DataFile%Name  )  
    if (DataFile%Status==0) then
      iTraj = iTraj+1 

      Idx                   = iTemp1
      iPES                  = iTemp2
      This%bMax(iTraj)      = Temp1
      This%bSampled(iTraj)  = Temp2
      This%Qini(:,iTraj)    = TempVec1
      This%Qfin(:,iTraj)    = TempVec2

      if (This%StatWritesBinaryFlg) then
        write(UnitWrite) int(Idx,  rkp)
        write(UnitWrite) int(iPES, rkp)
        write(UnitWrite) This%bMax(iTraj)
        write(UnitWrite) This%bSampled(iTraj)
        do iCond=1,This%NCond
          write(UnitWrite) This%Qini(iCond,iTraj)
        end do
        do iCond=1,This%NCond
          write(UnitWrite) This%Qfin(iCond,iTraj)
        end do
      end if

    end if

    if (iTraj==This%NTraj) exit
  end do                                                                                                        
  if (i_Debug_Loc) then
    call Logger%Write( "-> Done reading the trajectory data" )
    call Logger%Write( "-> Last line: iTraj = ", "This%bMax(iTraj) = ", This%bMax(This%NTraj), "This%bSampled(iTraj) = ", This%bSampled(This%NTraj), Fi="i9", Fr="es15.8")
  end if
  call DataFile%Close()

  if (This%StatWritesBinaryFlg) close(UnitWrite)
! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
Subroutine ReadInputsUnformatted( This, i_Debug )

  use, intrinsic :: iso_fortran_env ,only:  IOStat_End

  class(Statistics_Type)                                ,intent(inout)  ::    This
  logical                                       ,optional ,intent(in)   ::    i_Debug

  logical                                                               ::    i_Debug_Loc
  integer(rkp)                                                          ::    iTraj, Idx, iPES
  logical                                                               ::    Limit                                               ! TRUE if There is a limit on the Nb of Trajs to Analyze
  Integer                                                               ::    iCond
  Integer                                                               ::    POSTemp, TotBytes
  Integer                                                               ::    UnitWrite, StatusWrite
  type(File_Type)                                                       ::    DataFile

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadInputs" )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!     OPENING AND WRITING HEADER FOR THE PROGRESS FILE
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the number of trajectories to be analyzed" )
  Limit   =   This%NTrajectoriesToAnalyze > 0 
  iTraj   =   0

  if (i_Debug_Loc) call Logger%Write( "Opening the data file" )
  DataFile%Name     =   'trajectories.bin'
  open( NewUnit=DataFile%Unit, File=DataFile%Name, Action='READ', access="Stream", form="Unformatted", iostat=DataFile%Status )
  if (DataFile%Status/=0) call Error( "Error opening file: " // DataFile%Name )
    
    read(DataFile%Unit, POS=1) iTraj
    if (i_Debug_Loc) call Logger%Write( "-> The Unformatted File Contains ", iTraj, " Trajectories" )
    
    if ( Limit .and. iTraj > This%NTrajectoriesToAnalyze ) then
      This%NTraj = This%NTrajectoriesToAnalyze 
    else
      This%NTraj = int(iTraj, 4)
    end if    
    if (i_Debug_Loc) call Logger%Write( "-> Number of trajectories: This%NTraj = ", This%NTraj )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !     ALLOCATE THE DATA
    ! ==============================================================================================================
    allocate( This%bMax(This%NTraj) )
    allocate( This%bSampled(This%NTraj) )
    allocate( This%Qini(This%NCond,This%NTraj) )
    allocate( This%Qfin(This%NCond,This%NTraj) )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   IDENTIFYING AND SORTING INITIAL STATES AND IMPACT PARAMETER RINGS
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Reading the trajectory data: bMax, bSampled, Qini, Qfin" )
    
    TotBytes = int(2 * rkp) + int(2 * rkp) + int(This%NCond * rkp) + int(This%NCond * rkp)
    do iTraj = 1,This%NTraj

      POSTemp = rkp + (iTraj-1)*TotBytes + 1
      read(DataFile%Unit, POS=PosTemp, iostat=DataFile%Status ) Idx
      
      POSTemp = rkp + (iTraj-1)*TotBytes + 1 + rkp
      read(DataFile%Unit, POS=PosTemp, iostat=DataFile%Status ) iPES
      
      POSTemp = rkp + (iTraj-1)*TotBytes + 1 + int(2*rkp)
      read(DataFile%Unit, POS=PosTemp, iostat=DataFile%Status ) This%bMax(iTraj)
      
      POSTemp = rkp + (iTraj-1)*TotBytes + 1 + int(3*rkp)
      read(DataFile%Unit, POS=PosTemp, iostat=DataFile%Status ) This%bSampled(iTraj)  
      
      do iCond = 1,This%NCond
        PosTemp = rkp + (iTraj-1)*TotBytes + 1 + int(3*rkp) + iCond*rkp
        read(DataFile%Unit, POS=PosTemp, iostat=DataFile%Status ) This%Qini(iCond,iTraj)  
      end do
      
      do iCond = 1,This%NCond
        PosTemp = rkp + (iTraj-1)*TotBytes + 1 + int(3*rkp) + int(This%NCond*rkp) + iCond*rkp
        read(DataFile%Unit, POS=PosTemp, iostat=DataFile%Status ) This%Qfin(iCond,iTraj)  
      end do                                                                              
      
      if (DataFile%Status/=0) call Error( "Error reading the data file for statistics: " // DataFile%Name  )      
      if (i_Debug_Loc) then
        call Logger%Write( "-> Done reading the trajectory data" )
        call Logger%Write( "-> Last line: iTraj = ", "This%bMax(iTraj) = ", This%bMax(This%NTraj), "This%bSampled(iTraj) = ", This%bSampled(This%NTraj), Fi="i9", Fr="es15.8")
      end if
    end do
    
  call DataFile%Close()
  ! ==============================================================================================================


  ! ==============================================================================================================
  open( NewUnit=UnitWrite, File='./NConvTraj.dat', Action='WRITE', form="Formatted", iostat=StatusWrite )
  if (StatusWrite/=0) call Error( "Error writing the Nb of Converged Trajectories: " // './NConvTraj.dat' )
    write(UnitWrite, '(I10)') This%NTraj
  close(UnitWrite)     
  ! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
Subroutine SetRings( This, i_Debug )

  use Parameters_Module     ,only:  Zero, One, Two, Pi
  use Sorting_Module        ,only:  hpsort

  class(Statistics_Type)                                ,intent(inout)  ::    This
  logical                                     ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ,parameter  ::    ImpactParameterTolerence = 1.0E-03_rkp              ! Min difference between Impact Par's Maxs Allowed.
  real(rkp)                                                             ::    DeltabMax                                           ! Min difference between Impact Par's Maxs.
  integer(rkp)                                                          ::    iTraj, iRing                             
  real(rkp)                                                             ::    Temp
  real(rkp) ,dimension(This%NTraj)                                      ::    bMaxTemp
  logical                                                               ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetRings" )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!   SETTING THE NUMBER OF RINGS
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "-> Current Number of rings: This%NRings = ", This%NRings )
  bMaxTemp = This%bMax
  iRing = 0                                                                                                 
  do iTraj = 1,This%NTraj                                                                                                         ! Loop on all trajectories
    if ( iRing > 0 ) then                                                                                                         ! If some rings have already been added, then check that the current potential ring is ...
      DeltabMax = minval( abs( bMaxTemp(iTraj) - bMaxTemp(1:iRing) ) )                                                            ! Computing the minimum difference between current bMax and the one from all the previous rings
      if ( DeltabMax < ImpactParameterTolerence ) cycle                                                                           ! Do not accept the ring if the bMax of current trajectory is too close to the ones already processed
    end if
    iRing           = iRing + 1                                                                                                   ! Incrementing the number of rings
    bMaxTemp(iRing) = bMaxTemp(iTraj)                                                                                                 
  end do
  This%NRings = iRing
  if ( allocated(This%bMax) ) deallocate(This%bMax)                                                                               ! Renewing bMax with only the accepted values
  allocate( This%bMax(iRing), source = bMaxTemp(1:iRing) )  
  if (i_Debug_Loc) call Logger%Write( "-> Updated Number of rings: This%NRings = ", This%NRings )
! ==============================================================================================================


! ==============================================================================================================
!   SORTING THE RINGS
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Sorting the rings in increasing order" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling hpsort" )
  call hpsort( This%bMax )                                                                                                        ! bMax = Sorted Vector of Impact Params
  if (i_Debug_Loc) call Logger%Write( "-> Sorted This%bMax = ", This%bMax, Fr="es15.8" )
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE RINGS AREA
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the rings area" )
  if ( allocated(This%RingArea) ) deallocate(This%RingArea)
  allocate( This%RingArea(This%NRings) )
  Temp = Zero
  do iRing = 1,This%NRings
    This%RingArea(iRing) = Temp
    Temp                 = Pi * This%bMax(iRing)**2
    This%RingArea(iRing) = Temp - This%RingArea(iRing)                                                                            ! RingArea = Vector of Impact Param Ring Area
  end do
  if (i_Debug_Loc) call Logger%Write( "-> This%RingArea = ", This%RingArea, Fr="es15.8" )

  if ( (This%NRings == 1) .and. (This%RingArea(1) == Zero) ) then
    if (i_Debug_Loc) call Logger%Write( "Outputting probabilities, not Cross-sections" )
    This%RingArea(1)     =   One
  end if
  ! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting()


End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
Subroutine IdentifyInitialStates( This, i_Debug, i_Debug_Deep )

  use Parameters_Module     ,only:  Zero, One, Two, Half
  use Sorting_Module        ,only:  hpsort

  class(Statistics_Type)                                ,intent(inout)  ::    This
  logical                                     ,optional ,intent(in)     ::    i_Debug
  logical                                     ,optional ,intent(in)     ::    i_Debug_Deep

  integer(rkp)                                                          ::    IniCond                                             ! Trajectory Inital Condition
  real(rkp)                                                             ::    FinCond                                             ! Trajectory Fin Condition
  integer(rkp)                                                          ::    iTraj, iCond, iPairIni, iPairFin
  integer(rkp)                                                          ::    isum
  logical                                                               ::    i_Debug_Loc
  logical                                                               ::    i_Debug_Deep_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "IdentifyInitialStates" )
  !i_Debug_Loc   =     Logger%On()
  i_Debug_Deep_Loc = .false.; if ( present(i_Debug_Deep) )i_Debug_Deep_Loc = i_Debug_Deep
  

  if (i_Debug_Loc) call Logger%Write( "Identifying and sorting initial states" )

  allocate( This%IniStateCode(This%NTraj) )
  allocate( This%SortedIndx_IniStateCode(This%NTraj) )

  do iTraj = 1,This%NTraj

    This%ifact = 1
    isum       = 0
   
    iPairIni   = int((This%Qini(This%NCond,iTraj) -0.49_rkp)/16)   ! Finding Initial Arrangement
    iPairFin   = int((This%Qfin(This%NCond,iTraj) -0.49_rkp)/16)   ! Finding Final Arrangement
  
    do iCond = 1,This%NCond
    
      This%Qini(iCond,iTraj) = This%Qini(iCond,iTraj) - Half
      IniCond                = nint(This%Qini(iCond,iTraj))
      
      if ( IniCond >= This%QuantumNumberMax ) then
        call Logger%Write( "Error: Quantum number is above maximum allowed value - abort IniCond = ", IniCond )
        call Logger%Write( "Error: -> iTraj = ", iTraj, "IniCond = ", IniCond, "This%QuantumNumberMax = ", This%QuantumNumberMax )
        call Error( "Error: Quantum number is above maximum allowed value" )
      end if
      
      
      FinCond                = This%Qfin(iCond,iTraj) 
      if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", 1 - This%Qfin(:,iTraj)   = ", This%Qfin(:,iTraj) )                                                          
      if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", FinCond              = ", FinCond )                                                          
      
      This%Qfin(iCond,iTraj) = This%Qfin(iCond,iTraj) - Half                                                                      ! QUANTUM RULES ???
      if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", 2- This%Qfin(:,iTraj)   = ", This%Qfin(:,iTraj) )                                                          
      if ( This%PresEvOdd(iCond) .and. (iPairIni == iPairFin) ) then                                                              ! IF ((iCond == Pair/Odd RotQN) and (In/Elastic Collision))
        if ( mod(IniCond,2) == 1 ) then                                                                                           !    => IF (ODD RotQN(Ini))
          This%Qfin(iCond,iTraj) = This%Qfin(iCond,iTraj) - One                                                                   !          => RotQN(Fin) = RotQN(Fin) - 1
          if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", 3 - This%Qfin(:,iTraj)   = ", This%Qfin(:,iTraj) )                                                          
          if ( (nint(This%Qfin(iCond,iTraj)*Half) < 0) .and. (FinCond > Zero) ) then                                              !          => IF RotQN(Fin) == -1   
            This%Qfin(iCond,iTraj) = This%Qfin(iCond,iTraj) + Two                                                                 !                    => RotQN(Fin) = 1
            if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", 4 - This%Qfin(:,iTraj)   = ", This%Qfin(:,iTraj) )                                                          
          end if
        end if                                                                                                                    !
        This%Qfin(iCond,iTraj) = This%Qfin(iCond,iTraj) * Half                                                                    !    => RotQN(Fin) = RotQN(Fin) /2
        if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", 5 - This%Qfin(:,iTraj)   = ", This%Qfin(:,iTraj) )                                                          
      end if
     
      IniCond = IniCond + 2                                                                                                       ! Creating a Traj Identification Nb (IniStateCode) based on the IniCond
      if ( This%iskip(iCond) /= 0 ) then                                                                                         
        isum = isum + IniCond * This%ifact
      else
        isum = isum + 2 * This%ifact
      end if
      This%ifact = This%ifact * This%QuantumNumberMax
     
    end do
    
    This%IniStateCode(iTraj) = isum                                                                                                    
    if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", This%ifact = ", This%ifact, ", This%IniStateCode(iTraj) = ", This%IniStateCode(iTraj) )
    if ( i_Debug_Deep ) call Logger%Write( "-> iTraj = ", iTraj, ", 6 - This%Qfin(:,iTraj)   = ", This%Qfin(:,iTraj) )
    
  end do
  
  call hpsort( This%IniStateCode, This%SortedIndx_IniStateCode )      

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
Subroutine PrepareOutputFiles( This, i_Debug )

  use String_Module         ,only:  Convert_To_String

  class(Statistics_Type)                                ,intent(inout)  ::    This
  logical                                     ,optional ,intent(in)     ::    i_Debug

  logical                                                               ::    i_Debug_Loc
  integer(rkp)                                                          ::    N, i
  character(:)  ,allocatable                                            ::    Ni, Nr
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadInputs" )
  !i_Debug_Loc   =     Logger%On()


! ! ==============================================================================================================
! !   PREPARING THE OUTPUT FILE TO WRITE RESIDUALS ('statistics-residuals.out')
! ! ==============================================================================================================
!   if (i_Debug_Loc) call Logger%Write( "Preparing the output file to write residuals" )
!   This%ResidOutFile%Name = 'statistics-residuals.out'
!   if (i_Debug_Loc) call Logger%Write( "-> Opening file: This%ResidOutFile%Name = ", This%ResidOutFile%Name )
!   open( NewUnit=This%ResidOutFile%Unit, File=This%ResidOutFile%Name, Action='WRITE', Status='REPLACE', Form='FORMATTED', iostat=This%ResidOutFile%Status )
!   if (This%ResidOutFile%Status/=0) call Error( "Error opening file: " // This%ResidOutFile%Name )
!   N  = This%NCond
!   Ni = Convert_To_String(int(N+2,4))
!   Nr = Convert_To_String(int(N,4))
!   This%ResidOutFile%Format = "(1x,"//Ni//"i6,1p"//Nr//"e15.7)"
!   if (i_Debug_Loc) call Logger%Write( "-> Setting the write format to: This%ResidOutFile%Format = ", This%ResidOutFile%Format )
!   if (i_Debug_Loc) call Logger%Write( "-> Writing the header")
!   write(This%ResidOutFile%Unit,"('#',"//Ni//"a6,1p"//Nr//"a15)") ("i("//Convert_To_String(int(i,4))//")",i=1,N), "NRMS", "NCont",("RMS("//Convert_To_String(int(i,4))//")",i=1,N)
!   flush(This%ResidOutFile%Unit)
! ! ==============================================================================================================


! ! ==============================================================================================================
! !   PREPARING THE OUTPUT FILE TO WRITE PROBABILITIES DATA ('statistics-probabilities.out')
! ! ==============================================================================================================
!   if (i_Debug_Loc) call Logger%Write( "Preparing the output file to write probabilities data" )
!   This%ProbaOutFile%Name = 'statistics-probabilities.out'
!   if (i_Debug_Loc) call Logger%Write( "-> Opening file: This%ProbaOutFile%Name = ", This%ProbaOutFile%Name )
!   open( NewUnit=This%ProbaOutFile%Unit, File=This%ProbaOutFile%Name, Action='WRITE', Status='REPLACE', Form='FORMATTED', iostat=This%ProbaOutFile%Status )
!   if (This%ProbaOutFile%Status/=0) call Error( "Error opening file: " // This%ProbaOutFile%Name )
!   Ni  =   Convert_To_String(int(N+1,4))
!   This%ProbaOutFile%Format = "(1x,"//Ni//"(i5,3x),*(e15.7,3x))"
!   if (i_Debug_Loc) call Logger%Write( "-> Setting the write format to: This%ProbaOutFile%Format = ", This%ProbaOutFile%Format )
!   if (i_Debug_Loc) call Logger%Write( "-> Writing the header")
!   write(This%ProbaOutFile%Unit,"('#',"//Ni//"(a5,3x),*(a15,3x))") 'idx', ("i("//Convert_To_String(int(i,4))//")",i=1,N), 'Cross', 'CrossSD'
!   flush(This%ProbaOutFile%Unit)
! ! ==============================================================================================================


! ==============================================================================================================
!   PREPARING THE OUTPUT FILE TO WRITE CROSS SECTIONS ('statistics.out')
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Preparing the output file to write cross sections" )
  This%StatOutFile%Name     =   'statistics.csv'
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: This%StatOutFile%Name = ", This%StatOutFile%Name )
  open( NewUnit=This%StatOutFile%Unit, File=This%StatOutFile%Name, Action='WRITE', Status='REPLACE', Form='FORMATTED', iostat=This%StatOutFile%Status )
  if (This%StatOutFile%Status/=0) call Error( "Error opening file: " // This%StatOutFile%Name )
  N   = This%NCond
  Ni  =   Convert_To_String(int(N*2,4))
  This%StatOutFile%Format = "(1x,"//Ni//"(i5,3x),*(e15.7,3x))"
  if (i_Debug_Loc) call Logger%Write( "-> Setting the write format to: This%StatOutFile%Format = ", This%StatOutFile%Format )
  if (i_Debug_Loc) call Logger%Write( "-> Writing the header")
  write(This%StatOutFile%Unit,"('#',"//Ni//"(a5,3x),*(a15,3x))") ("f("//Convert_To_String(int(i,4))//")",i=1,N),("i("//Convert_To_String(int(i,4))//")",i=1,N), 'Cross', 'CrossSD'
  flush(This%StatOutFile%Unit)
! ==============================================================================================================


! ! ==============================================================================================================
! !   PREPARING THE OUTPUT FILE TO WRITE STRATA CONTRIBUTES TO CROSS SECTIONS ('statistics-bSensitivity.out')
! ! ==============================================================================================================
!   if (i_Debug_Loc) call Logger%Write( "Preparing the output file to write strata contributes to cross sections" )
!   This%bSensitivityFile%Name     =   'statistics-bSensitivity.out'
!   if (i_Debug_Loc) call Logger%Write( "-> Opening file: This%bSensitivityFile%Name = ", This%bSensitivityFile%Name )
!   open( NewUnit=This%bSensitivityFile%Unit, File=This%bSensitivityFile%Name, Action='WRITE', Status='REPLACE', Form='FORMATTED', iostat=This%bSensitivityFile%Status )
!   if (This%bSensitivityFile%Status/=0) call Error( "Error opening file: " // This%bSensitivityFile%Name )
!   Ni  =   Convert_To_String(int(N*2,4))
!   This%bSensitivityFile%Format = "(1x,"//Ni//"(i5,3x),*(e15.7,3x))"
!   if (i_Debug_Loc) call Logger%Write( "-> Setting the write format to: This%bSensitivityFile%Format = ", This%bSensitivityFile%Format )
!   if (i_Debug_Loc) call Logger%Write( "-> Writing the header")
!   write(This%bSensitivityFile%Unit,"('#',(a10,3x))") 'NRings'
!   flush(This%bSensitivityFile%Unit)
!   write(This%bSensitivityFile%Unit,"('$',(i10,3x))") This%NRings
!   flush(This%bSensitivityFile%Unit)
!   write(This%bSensitivityFile%Unit,"('#',"//Ni//"(a5,3x),*(a15,3x))") ("f("//Convert_To_String(int(i,4))//")",i=1,N),("i("//Convert_To_String(int(i,4))//")",i=1,N), &
!                                                                       ("Cross("//Convert_To_String(int(i,4))//")",i=1,This%NRings),("CrossVar("//Convert_To_String(int(i,4))//")",i=1,This%NRings)
!   flush(This%bSensitivityFile%Unit)
! ! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine


Subroutine ProcessStatistics( This, i_Debug, i_Debug_Deep )

  use Parameters_Module     ,only:  Zero

  implicit none

  class(Statistics_Type)                                ,intent(inout)  ::    This
  logical                                     ,optional ,intent(in)     ::    i_Debug
  logical                                     ,optional ,intent(in)     ::    i_Debug_Deep

  integer(rkp)                                                          ::    IniStateCodeTemp                                    !< Current Identification Nb for the Trajectories' Initial States.
  integer(rkp),dimension(:)   ,allocatable                              ::    IniState            !Dim[NCond]                     !< Current Initial State (Set of Initial Conditions).
  integer(rkp)                                                          ::    NIniStateRep                                        !< Nb of times that the Initial State Repeats (for QSS).
  integer(rkp)                                                          ::    NRMS                                                !< Nb of Fin States per Initial State for which all the FinCond >= -0.99d0.
  integer(rkp)                                                          ::    NCont                                               !< Nb of Fin States per Initial State for which at least one FinCond < -0.99d0.
  real(rkp) ,dimension(:)   ,allocatable                                ::    RMS                 !Dim[NCond]                     !< Root Mean Square for (FinCond(iCond) - IniCond(iCond)).
  integer(rkp)                                                          ::    iFinStates                                          !< Nb of Different Fin States per Initial State.
  integer(rkp)   ,dimension(:)   ,allocatable                           ::    FinStateCode        !Dim[NTraj]                     !< Vector of Final States Identification Nbs.
  real(rkp) ,dimension(:,:) ,allocatable                                ::    FinWeight           !Dim[<= NRings,NTraj]           !< Fin Weight for the Rings Areas, in order to Compute Cross Sections.
  real(rkp)                                                             ::    bMaxElastic                                         !< Max Impact Parameter that generates Elastic Collisions.
  integer(rkp)   ,dimension(:)   ,allocatable                           ::    ToFinState          !Dim[NTraj]                     !< Mapping the Trajectory to its Final State.
  integer(rkp)   ,dimension(:)   ,allocatable                           ::    TrajsPerb        !Dim[<= NRings]                 !< Nb of Trajs with same IniConds per Impact Parameter Ring.
  integer(rkp)                                                          ::    MainIter                                            ! Unused ???
  integer(rkp)                                                          ::    IniCondTemp, isum                                         
  integer(rkp)                                                          ::    iCond, iIter, jIter, kIter
  integer(rkp)                                                          ::    ifact
  logical                                                               ::    i_Debug_Loc
  logical                                                               ::    i_Debug_Deep_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ProcessStatistics" )
  !i_Debug_Loc   =     Logger%On()
  i_Debug_Deep_Loc = .false.; if ( present(i_Debug_Deep) )i_Debug_Deep_Loc = i_Debug_Deep
  

  ! ==============================================================================================================
  !   ALLOCATING VARIABLES
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Allocation" )
  allocate( IniState(This%NCond) )
  allocate( TrajsPerb(This%NRings) )
  allocate( FinWeight(This%NRings,This%NTraj) )
  allocate( RMS(This%NCond) )
  allocate( FinStateCode(This%NTraj) )
  allocate( ToFinState(This%NTraj) )
  ! ==============================================================================================================

  ifact            =   This%ifact
  IniStateCodeTemp =   This%IniStateCode(1)
  iIter            =   1
  jIter            =   1
  bMaxElastic      =   Zero

  MainIter  =   0

  ! ==============================================================================================================
  !   LOOPING ON THE TRAJECTORIES INITIAL CONDITIONS
  ! ==============================================================================================================
  do    
    if (i_Debug_Loc) call Logger%Write( "MainIter = ", MainIter, NewLine=.True. )
    if ( MainIter > 100 ) exit

    isum       = IniStateCodeTemp
    This%ifact = ifact
    !if (i_Debug_Loc) call Logger%Write( "isum = ", isum, "This%ifact = ", This%ifact )
    
    if (i_Debug_Loc) call Logger%Write( "Loop on Trajectory's Conditions in reversed order" )
    do iCond = This%NCond,1,-1                                                                                                    ! Reconstructing the IniCond of the Trajectory
      !if (i_Debug_Loc) call Logger%Write( "iCond = ", iCond )
      
      if ( This%QuantumNumberMax == 0 ) call Error( "Error: This%QuantumNumberMax is zero" )
      
      This%ifact   =   This%ifact / This%QuantumNumberMax
      if (i_Debug_Loc) call Logger%Write( "  This%ifact = ", This%ifact )
      if ( This%ifact == 0 ) then
        call Logger%Write( "This%ifact is zero" )
        call Logger%Write( "ifact = ", ifact, "iCond = ", iCond, "This%QuantumNumberMax = ", This%QuantumNumberMax )
        call Error( "Error: ifact is zero" )
      end if
      
      IniCondTemp     = isum / This%ifact
      IniState(iCond) = IniCondTemp - 2
      !if (i_Debug_Loc) call Logger%Write( "iCond = ", iCond , "IniCondTemp = ", IniCondTemp, "IniState(iCond) = ", IniState(iCond) )
      
      isum = isum - IniCondTemp * This%ifact
      
    end do
    
    if (i_Debug_Loc) call Logger%Write( "Initial Condition Nb ", iIter, "; IniState = ", IniState )


    kIter        = jIter                                                                                                          ! Other Trajectories with the same IniCond?
    NIniStateRep = 0
    do
      if ( This%IniStateCode(jIter) /= IniStateCodeTemp ) cycle
      jIter        = jIter        + 1
      NIniStateRep = NIniStateRep + 1
      exit
    end do

    jIter = kIter
    if (adjustl(trim(This%AssignmentMethod)) == 'QSS') NIniStateRep = NIniStateRep * 2**This%NCond


    ! Initializing the counters for this initial state
    iFinStates   =   0
    NRMS         =   0
    NCont        =   0
    TrajsPerb    =   0
    FinWeight    =   Zero
    RMS          =   Zero

   
    ! ==============================================================================================================
    !   LOOPING ON THE FINAL STATES
    ! ==============================================================================================================
    do   
!@TODO: GENERALIZE: The following would not work if we want to compute cross sections from random initial states!
      !if ( This%IniStateCode(jIter) /= IniStateCodeTemp ) exit                         ! Repeating untill either the IniCond of the trajectory is the same OR we finished traj
                                                                                                                                  
      if (i_Debug_Deep_Loc) call Logger%Write( "Calling This%AddFinState: jIter = ",jIter )
      call This%AddFinState( jIter, TrajsPerb, iFinStates, FinWeight, FinStateCode, IniState, RMS, NRMS, NCont, bMaxElastic, ToFinState(This%SortedIndx_IniStateCode(jIter)), i_Debug=i_Debug_Deep_Loc )
     
      jIter = jIter + 1
      if ( jIter > This%NTraj ) exit
    end do
    ! ==============================================================================================================
    
    
    ! ==============================================================================================================
    !   WRITING STATISTICS FOR CURRENT INITIAL STATES TYPE
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Calling WriteFinStateProbabilities" )
    call WriteFinStateProbabilities( This, This%RingArea, TrajsPerb, iFinStates, FinWeight, FinStateCode, This%QuantumNumberMax, IniState, This%PresEvOdd, RMS, NRMS, NCont, ToFinState, This%bSampled, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Deep_Loc )
    ! ==============================================================================================================
    
  
    if ( jIter > This%NTraj ) exit
    IniStateCodeTemp = This%IniStateCode(jIter)
    iIter            = iIter + 1
  end do
  ! ==============================================================================================================
  

  if (i_Debug_Loc) call Logger%Exiting()
  

End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
! This procedures adds the current Fin state to the statistics we have accumulated so far
Subroutine AddFinState( This, jIter, TrajsPerb, iFinStates, FinWeight, FinStateCode, IniState,RMS, NRMS, NCont, bMaxElastic, ToFinState, i_Debug )

  use Parameters_Module     ,only:  Zero, One, Two

  class(Statistics_Type)                                ,intent(inout)  ::    This
  integer(rkp)                                          ,intent(in)     ::    jIter                                               ! Current Iteration
  integer(rkp)    ,dimension(:)                         ,intent(inout)  ::    TrajsPerb                                                  ! ? 
  integer(rkp)                                          ,intent(inout)  ::    iFinStates                 
  real(rkp)       ,dimension(:,:)                       ,intent(inout)  ::    FinWeight    
  integer(rkp)    ,dimension(:)                         ,intent(inout)  ::    FinStateCode  
  integer(rkp)    ,dimension(:)                         ,intent(in)     ::    IniState
  real(rkp)       ,dimension(:)                         ,intent(inout)  ::    RMS                                                 ! RMS = Sum( (IniCond - FinCond)^2 )   
  integer(rkp)                                          ,intent(inout)  ::    NRMS                                                ! Nb of Rings with none of FinCond <= -1.d0
  integer(rkp)                                          ,intent(inout)  ::    NCont                                               ! Nb of Rings with one or more FinCond <= -1.d0
  real(rkp)                                             ,intent(inout)  ::    bMaxElastic
  integer(rkp)                                          ,intent(out)    ::    ToFinState
  logical                                     ,optional ,intent(in)     ::    i_Debug

  integer(rkp)                                                          ::    MainIter
  logical                                                               ::    Found
  integer(rkp)                                                          ::    SortedIndx
  integer(rkp)                                                          ::    NCond
  integer(rkp)                                                          ::    NTraj
  integer(rkp)   ,dimension(2,size(This%Qfin,1))                        ::    FinCond
  integer(rkp)                                                          ::    NRings
  integer(rkp)                                                          ::    FinCondTemp, ArrDiv, ArrDiff, Temp
  integer(rkp)                                                          ::    IniArrDiv, FinArrDiv
  real(rkp)                                                             ::    InElasticFactor, QRulesFac, Diff, Weight
  real(rkp)                                                             ::    bSampled
  real(rkp)      ,dimension(2,size(This%Qfin,1))                        ::    WeightTemp
  integer(rkp)   ,dimension(size(This%Qfin,1))                          ::    ido   
  integer(rkp)                                                          ::    iRings, iRingsPlus, ii, iStates, iCond, jCond,  &
                                                                              iCont, isum, isumTemp, isumIdDiat, iquse,       & 
                                                                              ipart, ipart1, ipart2, ipartMin, ipartMax,      &
                                                                              ifact, iElastic, iu                                             
  logical                                                               ::    i_Debug_Loc   

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AddFinState" )
  !i_Debug_Loc   =     Logger%On()


  SortedIndx = This%SortedIndx_IniStateCode(jIter)
  NCond       = size(This%Qfin(:,SortedIndx),1)
  NRings      = size(This%bMax)
  bSampled    = This%bSampled(SortedIndx)

! ==============================================================================================================
!     DETERMINE WHICH IMPACT PARAMETER RING THIS TRAJECTORY CORRESPONDS TO
! ==============================================================================================================
  if (i_Debug_Loc) then
    call Logger%Write( "Trajectory with bSampled = ", bSampled, Fr="es15.8" )
    call Logger%Write( "This%Qfin(:,SortedIndx) = ", This%Qfin(:,SortedIndx), Fr="es15.8" )
  end if
  iRingsPlus = 1
  if ( NRings > 1 ) then
    do iRings = 1,NRings-1
      if ( (bSampled > This%bMax(iRings+1)) .or. (bSampled <= This%bMax(iRings)) ) cycle
      iRingsPlus = iRings + 1
      if (i_Debug_Loc) call Logger%Write( "iRing for this Trajectory = ", iRingsPlus, Fr="es15.8" )
      exit
    end do
  end if
  TrajsPerb(iRingsPlus) = TrajsPerb(iRingsPlus) + 1
  if (i_Debug_Loc) call Logger%Write( "iRing Nb ", iRingsPlus, " now contains ", TrajsPerb(iRingsPlus), " Trajectories." )
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING THE RMS DEVIATION BETWEEN INITIAL AND Fin QNS
! ==============================================================================================================
  iCont = 0
  do iCond = 1,NCond
    if ( This%Qfin(iCond,SortedIndx) < -0.99d0) iCont = iCont + 1
  end do
  if ( iCont == 0 ) then
    NRMS = NRMS + 1
    IniArrDiv = int( IniState(NCond) / 16 )
    FinArrDiv = int( (This%Qfin(NCond,SortedIndx) + 0.01d0) / 16.d0 )
    InElasticFactor = One
    if ( IniArrDiv == FinArrDiv ) InElasticFactor = Two
    do iCond = 1,NCond
      QRulesFac = One
      if ( This%PresEvOdd(iCond) ) QRulesFac = InElasticFactor
      Diff = This%Qfin(iCond,SortedIndx) * QRulesFac - dfloat( IniState(iCond) )
      RMS(iCond) = RMS(iCond) + Diff**2
    end do
  else
    NCont = NCont + 1
    if (i_Debug_Loc) call Logger%Write( "jIter-th Trajectory has at least one FinCond < -0.99d0. NCont = ", Ncont )
  end if
! ==============================================================================================================


! ==============================================================================================================
!     HISTOGRAM AND QSS METHODS
! ==============================================================================================================
  iElastic = 0
  do iCond = 1,NCond
  
    FinCond(1,iCond) = nint(This%Qfin(iCond,SortedIndx))
    if (i_Debug_Loc) call Logger%Write( "iCond = ", iCond, "; FinCond(1,iCond) = ", FinCond(1,iCond) )
    
    if ( dfloat(FinCond(1,iCond)) <= This%Qfin(iCond,SortedIndx) ) then
      FinCond(2,iCond) = FinCond(1,iCond) + 1
    else
      FinCond(2,iCond) = FinCond(1,iCond) - 1
    end if
    if (i_Debug_Loc) call Logger%Write( "FinCond(2,iCond) = ", FinCond(2,iCond) )
    
    if ((adjustl(trim(This%AssignmentMethod)) == 'Histogram') .or. (This%Qfin(iCond,SortedIndx) < Zero)) then
      WeightTemp(1,iCond)    =   One
      WeightTemp(2,iCond)    =   Zero
      ido(iCond)             =   0
      if ( FinCond(1,iCond) /= IniState(iCond) ) iElastic = iElastic + 1
    else
      WeightTemp(1,iCond)    =   One - ( This%Qfin(iCond,SortedIndx) - dfloat(FinCond(1,iCond)) )**2
      WeightTemp(2,iCond)    =   One - WeightTemp(1,iCond)
      ido(iCond)             =   1
    end if
    
  end do
  if (i_Debug_Loc) call Logger%Write( "FinCond(1,:)    = ", FinCond(1,:) )
  if (i_Debug_Loc) call Logger%Write( "FinCond(2,:)    = ", FinCond(2,:) )
  if (i_Debug_Loc) call Logger%Write( "WeightTemp(1,:) = ", WeightTemp(1,:) )
  if (i_Debug_Loc) call Logger%Write( "WeightTemp(2,:) = ", WeightTemp(2,:) )
  if (i_Debug_Loc) call Logger%Write( "ido(:)          = ", ido )

  if ( iElastic == NCond ) bMaxElastic = max(bMaxElastic,bSampled)                                                                ! Possible bug: bMaxElastic is unused

  MainIter  =  0
  Main_Loop: do

    isum      =   0
    ifact     =   1
    Weight    =   One
    MainIter  =   MainIter + 1
    if ( MainIter > 100 ) stop

    do iCond = 1,NCond

      iu    = max(ido(iCond),1)
      iquse = FinCond(iu,iCond)
      if (i_Debug_Loc) call Logger%Write( "iCond = ", iCond, "; iu = ", iu, "; iquse = ", iquse )

      if ( iCond == NCond ) then

        if (i_Debug_Loc) call Logger%Write( "Arrangement: FinCond(iu,iCond) = ", FinCond(iu,iCond) )
        
        ArrDiv    =  FinCond(iu,iCond)/16
        ArrDiff   =  FinCond(iu,iCond) - ArrDiv*16

        ipart1    =  0
        ipart2    =  0
        if ( NCond == 5 ) then
          ipart1   =  ArrDiff/4
          ipart2   =  ArrDiff - ipart1*4
          ipartMin =  min(ipart1,ipart2)
          ipartMax =  max(ipart1,ipart2)
          ipart    =  ipartMin + ipartMax*4
        end if
        if (i_Debug_Loc) call Logger%Write( "ArrDiv  = ", ArrDiv )
        if (i_Debug_Loc) call Logger%Write( "ArrDiff = ", ArrDiff )
        if (i_Debug_Loc) call Logger%Write( "ipart1  = ", ipart1 )
        if (i_Debug_Loc) call Logger%Write( "ipart2  = ", ipart2 )

        if (NCond == 3) then
        
          !if (ArrDiff == 2) then                                                                                           ! # UnComment for Matching CG-QCT Rates
          !  if (i_Debug_Loc) call Logger%Write( "Only1 1 Molecule in the Collision; Found an Unbound Fin State" )
          !  iquse = 2
          !else
          !  iquse = iquse * This%iskip(iCond)
          !end if                                                                                                           ! # UnComment for Matching CG-QCT Rates
          iquse = iquse * This%iskip(iCond)                                                                                 ! # Comment for Matching CG-QCT Rates

        elseif (NCond == 5) then
        
          ! if (ArrDiff == 10) then                                                                                         ! # UnComment for Matching CG-QCT Rates
          !   if (i_Debug_Loc) call Logger%Write( "2 Molecules in the Collision; Both the Fin States are Unbound" )
          !   iquse = 10
          ! else if (ipart1 == 2) then
          !   if (i_Debug_Loc) call Logger%Write( "2 Molecules in the Collision; The 1st Fin States is Unbound" )
          !   iquse = ipart2
          ! else if (ipart2 == 2) then
          !   if (i_Debug_Loc) call Logger%Write( "2 Molecules in the Collision; The 2nd Fin States is Unbound" )
          !   iquse = ipart1
          ! else
          !   iquse = iquse * This%iskip(iCond)
          ! end if                                                                                                          ! # UnComment for Matching CG-QCT Rates
          iquse = iquse * This%iskip(iCond)                                                                                 ! # Comment for Matching CG-QCT Rates

        end if

      else
        iquse = iquse * This%iskip(iCond)
      end if

      if ( iquse + 2 > This%QuantumNumberMax ) iquse = 0
      if ( (This%iskip(iCond) /= 0) .or. (iCond == NCond) ) isum = isum + (iquse+2) * ifact
      if (i_Debug_Loc) call Logger%Write( "iCond = ", iCond, "; iquse = ", iquse, "; isum = ", isum )
      Weight = Weight * WeightTemp(iu,iCond)
      ifact  = ifact  * This%QuantumNumberMax

    end do


    if ( This%IdenticalDiatoms ) then
    
      isumIdDiat  =   0
      ifact       =   1
      Weight      =   One
      
      do ii = 1,2
        FinCondTemp    =   FinCond(ii,1)                                                                                          ! Switching molecules 
        FinCond(ii,1)  =   FinCond(ii,3)
        FinCond(ii,3)  =   FinCondTemp
        FinCondTemp    =   FinCond(ii,2)
        FinCond(ii,2)  =   FinCond(ii,4)
        FinCond(ii,4)  =   FinCondTemp
        ArrDiv         =   FinCond(ii,5) / 16
        ArrDiff        =   FinCond(ii,5) - 16 * ArrDiv
        Temp           =   ArrDiff / 4
        FinCond(ii,5)  =   Temp + 16 * ArrDiv
      end do

      do iCond = 1,NCond
      
        iu          =   max(ido(iCond),1)
        iquse       =   FinCond(iu,iCond)
        
        if ( iCond == NCond ) then
          if (i_Debug_Loc) call Logger%Write( "Fin Arrangement: FinCond(iu,iCond) = ", FinCond(iu,iCond) )
          ArrDiv    =   FinCond(iu,iCond) / 16
          ArrDiff   =   FinCond(iu,iCond) - 16 * ArrDiv
          ipart1    =   0
          ipart2    =   0
          if ( NCond == 5 ) then
            ipart1  =   ArrDiff / 4
            ipart2  =   ArrDiff - 4 * ipart1
          end if
            
          if (NCond == 3) then
        
            ! if (ArrDiff == 2) then                                                                                        ! # UnComment for Matching CG-QCT Rates
            !   if (i_Debug_Loc) call Logger%Write( "Only1 1 Molecule in the Collision; Found an Unbound Fin State" )       ! #
            !   iquse = 2                                                                                                   ! #
            ! else                                                                                                          ! #
            !   iquse = iquse * This%iskip(iCond)                                                                           ! #
            ! end if                                                                                                        ! # UnComment for Matching CG-QCT Rates
            iquse = iquse * This%iskip(iCond)                                                                               ! # Comment for Matching CG-QCT Rates
            
          elseif (NCond == 5) then
          
            ! if (ArrDiff == 10) then                                                                                       ! # UnComment for Matching CG-QCT Rates
            !   if (i_Debug_Loc) call Logger%Write( "1 Molecule in the Collision; Found an Unbound Fin State" )
            !   iquse = 10
            ! else if (ipart1 == 2) then
            !   if (i_Debug_Loc) call Logger%Write( "2 Molecules in the Collision; Found an Unbound Fin State" )
            !   iquse = ipart2
            ! else if (ipart2 == 2) then
            !   if (i_Debug_Loc) call Logger%Write( "2 Molecules in the Collision; Found an Unbound Fin State" )
            !   iquse = ipart1
            ! else
            !   iquse = iquse * This%iskip(iCond)
            ! end if                                                                                                       ! # UnComment for Matching CG-QCT Rates
            iquse = iquse * This%iskip(iCond)                                                                              ! # Comment for Matching CG-QCT Rates
            
          end if
          
        else
          iquse = iquse * This%iskip(iCond)
        end if

        if ( iquse+2 > This%QuantumNumberMax ) iquse = 0
        if ( (This%iskip(iCond) /= 0) .or. (iCond == NCond) ) isumIdDiat = isumIdDiat + (iquse+2) * ifact
        if (i_Debug_Loc) call Logger%Write( "iCond = ", iCond, "iquse = ", iquse, "isum = ", isum )
        Weight = Weight * WeightTemp(iu,iCond)
        ifact  = ifact  * This%QuantumNumberMax
        
      end do

      if ( isumIdDiat > isum ) then
        isumTemp   =   isum
        isum       =   isumIdDiat
        isumIdDiat =   isumTemp
      end if

    end if

    Found = .False.
    if ( This%IdenticalDiatoms ) then
      do iStates = 1,iFinStates
        if ( (FinStateCode(iStates) /= isum) .and. (FinStateCode(iStates) /= isumIdDiat) ) cycle
        ToFinState = iStates
        Found      = .True.
      end do
    else
      do iStates = 1,iFinStates
        if ( FinStateCode(iStates) /= isum ) cycle
        ToFinState = iStates
        Found      = .True.
        exit
      end do
    end if

    if ( .Not. Found ) then
      iFinStates               = iFinStates + 1
      FinStateCode(iFinStates) = isum
      ToFinState               = iFinStates
    end if
    if (i_Debug_Loc) call Logger%Write( "Identification Nb for Final State:   iFinState = ", ToFinState )

    FinWeight(iRingsPlus,ToFinState) = FinWeight(iRingsPlus,ToFinState) + Weight
    if (i_Debug_Loc) call Logger%Write( "FinWeight(", iRingsPlus, ",", ToFinState, ") = ", FinWeight(iRingsPlus,ToFinState) )
    
    jCond = 1
    do
      !if (i_Debug_Loc) call Logger%Write( "jCond = ", jCond, "; NCond = ", NCond )
      if ( jCond > NCond ) exit Main_Loop
      
      !if (i_Debug_Loc) call Logger%Write( "ido(jCond) = ", ido(jCond) )
      
      if ( ido(jCond) == 1 ) then
        ido(jCond) = 2
        cycle Main_Loop
      end if
      
      if ( ido(jCond) == 0 ) then
        jCond = jCond + 1
        cycle
      end if
      
      if ( ido(jCond) == 2 ) then
        ido(jCond) = 1
        jCond      = jCond + 1
        cycle
      end if
      
      exit
      
    end do

  end do Main_Loop

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
! ***********************************************************************************************************************************************!


! ***********************************************************************************************************************************************!
! This procedure writes out the Fin state probabilities
Subroutine WriteFinStateProbabilities( This, RingArea, TrajsPerb, iFinStates, FinWeight, FinStateCode, QuantumNumberMax, IniState, PresEvOdd, RMS, NRMS, NCont, ToFinState, bSampled, i_Debug, i_Debug_Deep )

  use Parameters_Module     ,only:  Zero, One


  class(Statistics_Type)                    ,intent(inout)  ::    This
  real(rkp) ,dimension(:)                   ,intent(in)     ::    RingArea                                                        ! Vector of Impact Param Ring Area
  integer(rkp), dimension(:)                ,intent(in)     ::    TrajsPerb        
  integer(rkp)                              ,intent(in)     ::    iFinStates                                                      ! Number of Fin States
  real(rkp) ,dimension(:,:)                 ,intent(inout)  ::    FinWeight     
  integer(rkp), dimension(:)                ,intent(in)     ::    FinStateCode       
  integer(rkp)                              ,intent(in)     ::    QuantumNumberMax
  integer(rkp), dimension(:)                ,intent(in)     ::    IniState      
  logical   ,dimension(:)                   ,intent(in)     ::    PresEvOdd  
  real(rkp) ,dimension(:)                   ,intent(inout)  ::    RMS   
  integer(rkp)                              ,intent(in)     ::    NRMS
  integer(rkp)                              ,intent(in)     ::    NCont                                                           ! Number of Traj with at least 1 FinCond < -0.99d0
  integer(rkp),dimension(:)                 ,intent(inout)  ::    ToFinState    
  real(rkp) ,dimension(:)                   ,intent(in)     ::    bSampled        
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  integer(rkp)                                              ::    NRings
  integer(rkp)                                              ::    NCond
  integer(rkp)                                              ::    NTraj
  integer(rkp)                                              ::    IniArr, FinArr
  real(rkp)                                                 ::    Cross, CrossSD, SDPerc, CrossRing, CrossSDRing
  real(rkp)                                                 ::    bSampledMax
  integer(rkp), dimension(size(IniState))                   ::    FinCondVec  
  integer(rkp)                                              ::    jFinStates, iCond, iRings, iTraj, idws, jdws, i2, i0
  integer(rkp)                                              ::    ifact, ifacts, isum, FinCond
  character(:)  ,allocatable                                ::    Fmt
  logical                                                   ::    i_Debug_Loc 
  logical                                                   ::    i_Debug_Deep_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteFinStateProbabilities" )
  !i_Debug_Loc   =     Logger%On()
  i_Debug_Deep_Loc = .false.; if ( present(i_Debug_Deep) )i_Debug_Deep_Loc = i_Debug_Deep


  ! ==============================================================================================================
  !     WRITING THE RESIDUALS FILE
  ! ==============================================================================================================
  if (i_Debug_Deep_Loc) call Logger%Write( "Number of Fin States:   iFinStates = ", iFinStates )
  NRings = size(RingArea)
  NCond  = size(IniState)
  NTraj  = size(ToFinState)

  ! Normalizing FinWeight
  do jFinStates = 1,iFinStates
    if (i_Debug_Deep_Loc) call Logger%Write( "FinWeight(:,jFinStates) Before = ", FinWeight(:,jFinStates) )
    do iRings = 1,NRings
      FinWeight(iRings,jFinStates) = FinWeight(iRings,jFinStates) / max(1,TrajsPerb(iRings) )                                     ! ???  max(1,int(sum(FinWeight(jFinStates,:))))
    end do
    if (i_Debug_Deep_Loc) call Logger%Write( "FinWeight(:,jFinStates) After  = ", FinWeight(:,jFinStates) )
  end do

  ! Computing ifacts and RMS
  ifacts = 1
  do iCond = 1,NCond
    ifacts     = ifacts * QuantumNumberMax
    RMS(iCond) = sqrt( RMS(iCond) / max(1,NRMS) )
  end do
  
  !write(This%ResidOutFile%Unit,This%ResidOutFile%Format) IniState(1:NCond), NRMS, NCont, RMS(1:NCond)
  ! ==============================================================================================================
  

  do jFinStates = 1,iFinStates
    if (i_Debug_Deep_Loc) call Logger%Write( "For Fin state iFinStates = ", jFinStates )
  
  
    ! ==============================================================================================================
    !     WRITING THE CURRENT Fin STATES TYPE IN THE MAIN FILE
    ! ==============================================================================================================
    ! Computing Cross and CrossSD
    Cross   = Zero
    CrossSD = Zero
    do iRings = 1,NRings
      CrossRing   = RingArea(iRings) * FinWeight(iRings,jFinStates)
      Cross       = Cross   + CrossRing
      CrossSDRing = FinWeight(iRings,jFinStates) * (One-FinWeight(iRings,jFinStates)) * RingArea(iRings) * RingArea(iRings) / max(1,TrajsPerb(iRings))
      CrossSD     = CrossSD + CrossSDRing
      if (i_Debug_Deep_Loc) call Logger%Write( "jFinStates = ", jFinStates, ", iRings = ", iRings, ", CrossRing = ", CrossRing, ", CrossVarRing = ", CrossSDRing )
    end do
    CrossSD = sqrt(CrossSD)
    
    ! Reconstructing Fin Arrangements, Fin vqn and Fin jqn
    ifact  = ifacts
    isum   = FinStateCode(jFinStates)
    IniArr = IniState(NCond) / 16
    if (i_Debug_Deep_Loc) call Logger%Write( "isum = ", isum, "; FinStateCode(jFinStates) = ", FinStateCode(jFinStates) )
    do iCond = NCond,1,-1
    
      if ( QuantumNumberMax == 0 ) call Error( "[WriteFinStateProbabilities]: QuantumNumberMax is zero again" )
      ifact = ifact / QuantumNumberMax
      if ( ifact == 0 ) call Error( "[WriteFinStateProbabilities]: ifact is zero again" )
      FinCond = isum / ifact
      if ( iCond == NCond ) FinArr = (FinCond-2) / 16
      isum  = isum - FinCond * ifact
      
      FinCond = FinCond - 2
      if (i_Debug) call Logger%Write( "FinCond = ", FinCond )
      
      ! IF jqnIn is ODD, jqnFn stays ODD
      if ( (PresEvOdd(iCond) ) .and. (IniArr == FinArr)  ) then        
        FinCond = FinCond * 2
        if ( mod(IniState(iCond),2) /= 0 ) FinCond = FinCond + 1   
        if (i_Debug) call Logger%Write( "PresEvOdd(iCond) -> FinCond = ", FinCond )                                                                                 
      end if
    
      FinCondVec(iCond) = FinCond
      
    end do
    
    !write(This%bSensitivityFile%Unit,This%bSensitivityFile%Format) FinCondVec, IniState, [(RingArea(iRings) * FinWeight(iRings,jFinStates), iRings=1,NRings)], [(FinWeight(iRings,jFinStates) * (One-FinWeight(iRings,jFinStates)) * RingArea(iRings) * RingArea(iRings) / max(1,TrajsPerb(iRings)), iRings=1,NRings)]
   
    write(This%StatOutFile%Unit,This%StatOutFile%Format) FinCondVec, IniState, Cross, CrossSD
    ! ==============================================================================================================
    

    ! ==============================================================================================================
    !     WRITING THE CURRENT Fin STATES TYPE IN THE PROBABILITIES FILE
    ! ==============================================================================================================
    ! Finding the b max associated with the current Fin states type
    bSampledMax = Zero
    do iTraj = 1,NTraj
      if ( ToFinState(iTraj) /= jFinStates ) cycle
      ToFinState(iTraj) = 0
      bSampledMax = max(bSampledMax,bSampled(iTraj))
    end do
    
    ! Computing SDPerc
    SDPerc = 1.d2 * CrossSD / Cross
    if (i_Debug_Deep_Loc) call Logger%Write( "   Cross = ", Cross, "; CrossSD = ", CrossSD, "; SDPerc = ", SDPerc, "; bSampledMax = ", bSampledMax, Fr="es15.8" )
    ! ==============================================================================================================


    ! ! ==============================================================================================================
    ! !     WRITING ... ???
    ! ! ==============================================================================================================
    ! i2 = 2                                                                                                             
    ! if ( FinCondVec(1) < 0 ) i2 = 1                                                                                                     
    ! write(This%ProbaOutFile%Unit,This%ProbaOutFile%Format) i2, FinCondVec, Cross, CrossSD
    ! ! ==============================================================================================================
    

  end do

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
! ***********************************************************************************************************************************************!


End Module