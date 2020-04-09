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
Module LevelsContainer_Class

  use Parameters_Module       ,only:  rkp, Half, Zero, One, Two, Three, Four, Five, Ten, Hartree_To_eV, ue, Ukb
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error
  use Level_Class             ,only:  Level_Type
  use Degeneracy_Class        ,only:  Degeneracy_Type

  implicit none

  private
  public    ::    LevelsContainer_Type

  Type      ::    LevelsContainer_Type
    character(:) ,allocatable                         ::    PathToLevelsFile
    integer                                           ::    NStates     =   0   ! Old name: nst
    integer                                           ::    inj         =   0
    integer                                           ::    inv         =   0
    real(rkp)                                         ::    Tint        = Zero
    character(:) ,allocatable                         ::    jInMethod
    integer                                           ::    jIn         = 0
    integer                                           ::    iodd                ! Indicator of states to be included
    character(:) ,allocatable                         ::    vInMethod
    integer                                           ::    vIn
    integer                                           ::    StateIn     =  0   ! Indes of input state: Variable set in SetLevelsData
    integer                                           ::    UnitArbDist
    integer                                           ::    maxjqn      =   0   ! Maximum Rotational Q.N.
    integer                                           ::    maxvqn      =   0   ! Maximum Vibrational Q.N.
    real(rkp)                                         ::    mineint     = 1.0e100_rkp ! Energy of the Lowest Energy Level [Eh]
    real(rkp)                                         ::    mineinteV   = 1.0e100_rkp ! Energy of the Lowest Energy Level [eV]
    real(rkp)                                         ::    CutRatio    = Zero  ! Ratio Between Cut Levels Tot Partition Function and Original Levels Tot Partition Function
    real(rkp)                                         ::    QBoundRatio = Zero  ! Ratio Between Quasi-Bound Levels Tot Partition Function and Original Levels Tot Partition Function
    real(rkp)                                         ::    BoundRatio  = Zero  ! Ratio Between Bound Levels Tot Partition Function and Original Levels Tot Partition Function
    real(rkp)         ,dimension(:,:) ,allocatable    ::    CrossSec            ! Cross Sections
    real(rkp)         ,dimension(:,:) ,allocatable    ::    RateConst           ! Rate Constants
    real(rkp)         ,dimension(:,:) ,allocatable    ::    RateConst_Norm      ! Normalized Rate Constants
    real(rkp)         ,dimension(:)   ,allocatable    ::    RateConst_Final     ! Final Rate Constants
    real(rkp)         ,dimension(:)   ,allocatable    ::    RateConst_Sigma     ! Final Rate Constants St Deviations
    real(rkp)         ,dimension(:,:) ,allocatable    ::    RateConst_Arr
    real(rkp)         ,dimension(:,:) ,allocatable    ::    RateConst_Sigma_Arr
    real(rkp)         ,dimension(:,:) ,allocatable    ::    RateConst_Arr_Rnd
    real(rkp)         ,dimension(:,:) ,allocatable    ::    CArr                ! Arrhenius Coefficients
    integer                                           ::    NBound              ! Number of Bound States of the Levels Container
    integer                                           ::    NQBound             ! Number of Quasi-Bound States of the Levels Container
    real(rkp)                                         ::    Q                   ! Partition Function
    type(Level_Type)  ,dimension(:)   ,allocatable    ::    States              ! Array of interna states
    class(Degeneracy_Type)            ,allocatable    ::    Degeneracy
  contains
    procedure                                         ::    Initialize  => InitializeLevelsContainer
    procedure                                         ::    WriteList
    procedure   ,private                              ::    ReadLevelsData
    procedure   ,private                              ::    SetLevelsData
    procedure   ,private                              ::    SetUniformDistribution
    procedure   ,private                              ::    SetBoltmannDistribution
    procedure   ,private                              ::    SetArbitraryDistribution
  End Type
  
  public                                              ::    Compute_PartitionRatios
  logical                                             ::    i_Debug_Global = .False.

  contains

!________________________________________________________________________________________________________________________________!
Subroutine InitializeLevelsContainer( This, Input, DiatPot, iMol, FileName, NStates, ReCheckFlg, i_Debug )  !This%NStates,maxst,eint,egam,rmin,vmin,vmax,tau,jqn,vqn,ri,ro,rmax,iunit,iodd,inv,inj)
! This procedure reads the data associated to the vibrational-rotational states.
! This procedure was originally called 'statin' and was stored in the 'preini.f' file.
!     input variables:
!     ----------------
!     maxst - the first dimension of eint,tau,jqn, and vqn
!     iunit - the unit number of the file we are to read
!     output variables:
!     ----------------
!     eint(i) - internal energy of i'th quantum state
!     egam(i) - Half width of i'th quantum state
!     rmin(i) - the position of the potential minimum (included centrifugal potential) for i'th quantum state
!     vmin(i) - the value of the potential minimun (inc. cent. pot.)
!     vmax(i) - the value of the local potential maximum (inc. cent. pot.)
!     tau(i) - the vibrational period of the i'th quantum state
!     jqn(i) - the rotational q.n. of the i'th quantum state
!     vqn(i) - the vibrational q.n. of the i'th quantum state
!     ri(i) - inner turning point
!     ro(i) - outter turning point
!     rmax(i) - location of maximum in centrifugal barrier
!     iodd - if 0, include all states. if 1, include only odd states. if 2, include only even states. (rotational quantum number that is)

  use Input_Class              ,only:  Input_Type
  use DiatomicPotential_Class  ,only:  DiatomicPotential_Type

  class(LevelsContainer_Type)               ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  class(DiatomicPotential_Type)             ,intent(in)     ::    DiatPot           ! Intra-molecular diatomic potential object
  integer                                   ,intent(in)     ::    iMol
  character(*)                    ,optional ,intent(in)     ::    FileName         !> Name of the file containing the state data
  integer                         ,optional ,intent(in)     ::    NStates
  logical                         ,optional ,intent(in)     ::    ReCheckFlg
  logical                         ,optional ,intent(in)     ::    i_Debug
 
  integer                                                   ::    iState, jqn
  real(rkp)                                                 ::    VMaxTemp, VMinTemp, rMaxTemp, rMinTemp
  logical                                                   ::    TempFlg
  logical                                                   ::    ReCheckFlg_Loc
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeLevelsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  ReCheckFlg_Loc = .False.; if ( present(ReCheckFlg) ) ReCheckFlg_Loc = ReCheckFlg

! ! ==============================================================================================================
! !   READING FROM THE INPUT FILE THE DIATOMIC QUANTUM STATE INFORMATION OF MOLECULE
! ! ==============================================================================================================
  This%jIn          =   Input%jIn
  This%vIn          =   Input%vIn
  This%Tint         =   Input%Tint
  This%iodd         =   Input%iodd
  if (allocated(Input%jInMethod)) then
     This%jInMethod = Input%jInMethod
  else
     This%jInMethod = 'Null'
  endif
  if (allocated(Input%vInMethod)) then 
     This%vInMethod = Input%vInMethod
  else
     This%vInMethod = 'Null'
  endif
  This%UnitArbDist  =   Input%UnitArbDist
! ! ==============================================================================================================

! ==============================================================================================================
!   SETTING THE INDICATOR OF SELECTED STATES
! ==============================================================================================================
  if (i_Debug_Loc) then
    call Logger%Write( "-> This%iodd = ", This%iodd )
    select case (This%iodd)
      case(0); call Logger%Write( "-> Including all states" )
      case(1); call Logger%Write( "-> Including only states with odd js" )
      case(2); call Logger%Write( "-> Including only states with even js" )
    end select
  end if
! ==============================================================================================================


  if ( present(FileName) ) then
    allocate( This%PathToLevelsFile,  source = adjustl(trim( FileName )) ) 
    if (i_Debug_Loc) call Logger%Write( "Path to File with Pre-Processed Levels = ", This%PathToLevelsFile    )


    ! ==============================================================================================================
    !   READING THE STATE DATA
    ! ==============================================================================================================
      if (i_Debug_Loc) call Logger%Write( "Calling This%ReadLevelsData" )
      call This%ReadLevelsData( Input, iMol, DiatPot%xmui2, i_Debug=i_Debug_Loc )
      if (i_Debug_Loc) call Logger%Write( "-> Done reading level data" )
      if (i_Debug_Loc) call Logger%Write( "-> This%NStates = ", This%NStates )
    ! ==============================================================================================================

    if (ReCheckFlg_Loc) then
      ! ==============================================================================================================
      !  CHECKING LEVEL PROPERTIES
      ! ==============================================================================================================
        if (i_Debug_Loc) call Logger%Write( "Calling DiatPot%CheckMaxAndMin" )
        do jqn = 0,This%maxjqn
          TempFlg = .True.
          do iState = 1,This%NStates
            if (jqn == This%States(iState)%jqn) then
              if (TempFlg) then            
                call DiatPot%CheckMaxAndMin( This%States(iState), iState, i_Debug=i_Debug_Loc)
                TempFlg  = .False.
                VMaxTemp = This%States(iState)%VMaxNew
                VMinTemp = This%States(iState)%VMinNew
                rMaxTemp = This%States(iState)%rMaxNew
                rMinTemp = This%States(iState)%rMinNew 
              else
                This%States(iState)%VMaxNew = VMaxTemp
                This%States(iState)%VMinNew = VMinTemp
                This%States(iState)%rMaxNew = rMaxTemp
                This%States(iState)%rMinNew = rMinTemp
                if  (This%States(iState)%VMaxNew          < This%States(iState)%EInt) then
                  write(*,'(A,I5,A,I1,A,I2,A,I3,A)') '    [LevelsContainer_Class.F90 - InitializeLevelsContainer]: WARNING!   Level Nb ', iState, ' of Molecule Nb ', iMol, '(', This%States(iState)%vqn, ',', This%States(iState)%jqn, ') has Internal Energy larger than the Centrifugal Barrier!'
                elseif  (This%States(iState)%VMaxNew - 1.e-10 < This%States(iState)%EInt) then
                  write(*,'(A,I5,A,I1,A,I2,A,I3,A)') '    [LevelsContainer_Class.F90 - InitializeLevelsContainer]: WARNING!!! Level Nb ', iState, ' of Molecule Nb ', iMol, '(', This%States(iState)%vqn, ',', This%States(iState)%jqn, ') has Internal Energy significantly close to the Centrifugal Barrier!'
                end if
              end if
            end if
          end do
        end do
        if (i_Debug_Loc) call Logger%Write( "-> Done checking diatomic potential maxima and minima" )


        if (i_Debug_Loc) call Logger%Write( "Calling DiatPot%CheckTurningPoints" )
        do iState = 1,This%NStates
          call DiatPot%CheckTurningPoints( This%States(iState), iState, i_Debug=i_Debug_Loc)
        end do
        if (i_Debug_Loc) call Logger%Write( "-> Done checking level turning points" )



      ! ==============================================================================================================
    end if


    ! ==============================================================================================================
    !   SETTING THE STATE DATA DEPENDING ON WHETHER IT IS A REAL OR PSEUDO-STATE
    ! ==============================================================================================================
      if (i_Debug_Loc) call Logger%Write( "Calling This%SetLevelsData" )
      call This%SetLevelsData( i_Debug=i_Debug_Loc )
      if (i_Debug_Loc) call Logger%Write( "-> Done" )
      if (i_Debug_Loc) call Logger%Write( "-> This%NStates = ", This%NStates )
    ! ==============================================================================================================

  else

    This%NStates = NStates
    allocate(This%States(This%NStates))
    if (i_Debug_Loc) call Logger%Write( "Allocated This%States with Dimension (", This%NStates ,") " )

  end if


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadLevelsData( This, Input, iMol, xmui2, i_Debug )
  
  use Input_Class              ,only:  Input_Type
  use Degeneracy_Factory_Class ,only:  Degeneracy_Factory_Type

  class(LevelsContainer_Type)               ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iMol
  real(rkp)                                 ,intent(in)     ::    xmui2
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    Status
  integer                                                   ::    Unit
  integer                                                   ::    Parity
  integer                                                   ::    i       ! index of variables
  integer                                                   ::    iState  ! Index of states
  integer                                                   ::    NStates ! Tot nb of states
  integer                                                   ::    v       ! Vibration quantum number
  integer                                                   ::    j       ! Rotation quantum number
  integer   ,dimension(2)                                   ::    iArray
  real(rkp) ,dimension(9)                                   ::    rArray
  character(1)                                              ::    cdum
  
  Type(Degeneracy_Factory_Type)                             ::    Degeneracy_Factory

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadLevelsData")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!   OPENING THE FILE CONTAINING THE STATE DATA INFOMATION
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Opening state data file" )
  if (i_Debug_Loc) call Logger%Write( "-> This%PathToLevelsFile, = ", This%PathToLevelsFile )
  open( NewUnit=Unit, File=This%PathToLevelsFile, iostat=Status)
  if (Status/=0) call Error( "Error opening file '" // This%PathToLevelsFile // "' for state data" )
! ==============================================================================================================

                        Parity  =   2   ! Setting to even
  if ( This%iodd == 0 ) Parity  =   1   ! Setting to odd


! ==============================================================================================================
!   GETTING THE NUMBER OF STATES TO BE CONSIDERED AND ALLOCATION
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Getting the number of states to be considered" )
  This%NStates    =   0
  do
    read(Unit,*,iostat=Status) cdum
    if (Status<0) exit                                                                                           
    if (Status>0) call Error( "Error reading file '" // This%PathToLevelsFile // "' for state data" )
    if (cdum=='#') cycle
    backspace(Unit)
    read(Unit,*,iostat=Status) iArray
    if (Status<0) exit                                                                                           
    if (Status>0) call Error( "Error reading file '" // This%PathToLevelsFile // "' for state data" )
    v    =   iArray(1)                                                                                                            ! Getting the vibrational quantum number of current state
    j    =   iArray(2)                                                                                                            ! Getting the rotational quantum number of current state
    if ( ( trim(adjustl(This%jInMethod)) .ne. "Fixed" ) .and. ( trim(adjustl(This%vInMethod)) .eq. "Fixed" ) .and. ( v /= This%vIn ) ) cycle! j is sampled but current v /= This%vFixed; going to the next state
    if ( mod(j+This%iodd,Parity) == 0 ) This%NStates = This%NStates + 1                                                           ! Incrementing number of states, accounting for parity if required
  end do
  if (i_Debug_Loc) call Logger%Write( "-> This%NStates = ", This%NStates )
  if ( This%NStates == 0 ) call Error( "Error: No internal levels found in file '" // This%PathToLevelsFile // "'" )
  if (i_Debug_Loc) call Logger%Write( "-> Allocating This%States" )
  if ( allocated(This%States) ) deallocate(This%States)
  allocate( This%States(This%NStates) )
! ==============================================================================================================


! ==============================================================================================================
!   ALLOCATING DEGENERACIES
! ==============================================================================================================
  if ( i_Debug_Loc ) call Logger%Write( "Computing Degeneracies for Molecule Nb", iMol )  
  if ( i_Debug_Loc ) call Logger%Write( "-> Degeneracy_Factory%Define_Degeneracy" )
  call Degeneracy_Factory%Define_Degeneracy( Input, This%Degeneracy, iMol, i_Debug=i_Debug_Loc ) 
! ==============================================================================================================


! ==============================================================================================================
!   GETTING THE STATES DATA
! ==============================================================================================================
  rewind(Unit)
  iState     =   0
  NStates    =   0 
  do
    read(Unit,*,iostat=Status) cdum
    if (Status<0) exit                                                                                          
    if (Status>0) call Error( "Error reading file '" // This%PathToLevelsFile // "' for state data" )
    if (cdum=='#') cycle
    backspace(Unit)
    read(Unit,*,iostat=Status) iArray, rArray
    if (Status<0) exit                                                                                           
    if (Status>0) call Error( "Error reading file '" // This%PathToLevelsFile // "' for state data" )
    v    =   iArray(1)                                                                                                            ! Getting the vibrational quantum number of current state
    j    =   iArray(2)                                                                                                            ! Getting the rotational quantum number of current state
    if ( ( trim(adjustl(This%jInMethod)) .ne. "Fixed" ) .and. ( trim(adjustl(This%vInMethod)) .eq. "Fixed" ) .and. ( v /= This%vIn ) ) cycle! j is sampled but current v /= This%vFixed; going to the next state
    if ( mod(j+This%iodd,Parity) == 0 ) iState = iState + 1                                                                       ! Incrementing number of states, accounting for parity if required
      associate( State => This%States(iState) )
      i                   =     1
      State%vqn           =     iArray(i) ; i = i + 1
      State%jqn           =     iArray(i)
      i                   =     1
      State%eint          =     rArray(i); i = i + 1
      State%einteV        =     State%eint * Hartree_To_eV 
      State%egam          =     rArray(i) ; i = i + 1
      State%rmin          =     rArray(i) ; i = i + 1
      State%rmax          =     rArray(i) ; i = i + 1
      State%vmin          =     rArray(i) ; i = i + 1
      State%vmax          =     rArray(i) ; i = i + 1
      State%tau           =     rArray(i) ; i = i + 1
      State%ri            =     rArray(i) ; i = i + 1
      State%ro            =     rArray(i)
      State%rlim          =     Zero
      State%Vc_R2         =     xmui2 * ( State%jqn + Half )**2
      if (State%vqn > This%maxvqn) then
        This%maxvqn = State%vqn
      end if
      if (State%jqn > This%maxjqn) then
        This%maxjqn = State%jqn
      end if
      if (State%eint < This%mineint) then
        This%mineint    = State%eint
        This%mineinteV  = State%einteV
      end if
      State%g = This%Degeneracy%Compute_Degeneracy_State( State%jqn )
    end associate
    NStates = NStates + 1
  end do
  do iState = 1,NStates 
    This%States(iState)%einteV_scaled = This%States(iState)%einteV - This%mineinteV
  end do
  close(Unit)
! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteList( This, FileName, SortLevelsFlg, i_Debug) 
! This procedure ...
  
  use Error_Class                   ,only:  Error
  use Sorting_Module                

  class(LevelsContainer_Type)                            ,intent(in)     ::    This
  character(*)                                           ,intent(in)     ::    FileName
  logical                                  ,optional     ,intent(in)     ::    SortLevelsFlg
  logical                                  ,optional     ,intent(in)     ::    i_Debug

  integer                                                                ::    Status
  integer                                                                ::    Unit 
  integer                                                                ::    iState
  integer          ,dimension(:)      ,allocatable                       ::    IdxVec
  integer                                                                ::    NStates
  real(rkp)        ,dimension(:)      ,allocatable                       ::    StatesEint
  logical                                                                ::    SortLevelsFlg_Loc 
  logical                                                                ::    i_Debug_Loc 


  SortLevelsFlg_Loc = .False.; if ( present(SortLevelsFlg) ) SortLevelsFlg_Loc = SortLevelsFlg

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteList")
  !i_Debug_Loc   =     Logger%On()
  

  NStates = This%NStates


  allocate(IdxVec(NStates))
  if (SortLevelsFlg_Loc) then
    allocate(StatesEint(NStates))
    if (i_Debug_Loc) call Logger%Write( "-> Sorting the Energy Level before writing them" )
    do iState = 1,NStates
      StatesEint(iState) = This%States(iState)%eint
    end do
    if (SortLevelsFlg_Loc) then
      call hpsort( StatesEint, IdxVec )
    end if     
  else         
    do iState = 1,NStates
      IdxVec(iState) = iState
    end do 
  end if


  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", FileName )
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName ) 
  
    !write(Unit,'(a)') ('######################################################################################################################################################')
    !write(Unit,'(a)') ('# jqn   : the rotational q.n. of the i''th quantum state')
    !write(Unit,'(a)') ('# vqn   : the vibrational q.n. of the i''th quantum state')
    !write(Unit,'(a)') ('# eint  : internal energy of i''th quantum state [Eh]')
    !write(Unit,'(a)') ('# egam  : Half width of i''th quantum state')
    !write(Unit,'(a)') ('# rmin  : the position of the potential minimum (included centrifugal potential) for i''th quantum state')
    !write(Unit,'(a)') ('# vmin  : the value of the potential minimun (inc. cent. pot.)')
    !write(Unit,'(a)') ('# vmax  : the value of the local potential maximum (inc. cent. pot.)')
    !write(Unit,'(a)') ('# tau   : the vibrational period of the i''th quantum state')
    !write(Unit,'(a)') ('# ri    : inner turning point')
    !write(Unit,'(a)') ('# ro    : outter turning point')
    !write(Unit,'(a)') ('# rmax  : location of maximum in centrifugal barrier')
    !write(Unit,'(a)') ('######################################################################################################################################################')
    !write(Unit,'(a)') ('#   vqn  jqn      eint           egam           rmin           rmax           vmin           vmax            tau             ri             ro')
    !write(Unit,'(a)') ('######################################################################################################################################################')
    
    write(Unit,'(A)') ('#================================================================================================================================================')
    write(Unit,'(A)') ('#')                                                                                              
    write(Unit,'(A)') ('#================================================================================================================================================')
    write(Unit,'(A)') ('# vqn, jqn,         E[Eh],      EGam[au],      rMin[a0],      rMax[a0],      VMin[Eh],      VMax[Eh],       Tau[au],       rIn[a0],      rOut[a0]')

    do iState = 1,NStates
    
      write(Unit, 1)  This%States(IdxVec(iState))%vqn,  &
                      ',',                              &
                      This%States(IdxVec(iState))%jqn,  &
                      ',',                              &
                      This%States(IdxVec(iState))%eint, &
                      ',',                              &
                      This%States(IdxVec(iState))%egam, &
                      ',',                              &
                      This%States(IdxVec(iState))%rmin, &
                      ',',                              &
                      This%States(IdxVec(iState))%rmax, &
                      ',',                              &
                      This%States(IdxVec(iState))%Vmin, &
                      ',',                              &
                      This%States(IdxVec(iState))%Vmax, &
                      ',',                              &
                      This%States(IdxVec(iState))%tau,  &
                      ',',                              &
                      This%States(IdxVec(iState))%ri,   &
                      ',',                              &
                      This%States(IdxVec(iState))%ro
                      
    end do
    
  close(Unit)


  1 format(1X, I4, A, I4, 9(A, es14.7) )


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetLevelsData( This, i_Debug )

  class(LevelsContainer_Type)               ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  real(rkp) ,parameter                                      ::    Factor  =   3.157752E5_rkp  ! ?
  logical                                                   ::    Found
  integer                                                   ::    i, iState
  real(rkp)                                                 ::    Beta

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetLevelsData")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!   CASE WHEN INPUT VALUES OF V AND J ARE PROVIDED
! ==============================================================================================================
  if ( ( trim(adjustl(This%vInMethod)) .eq. "Fixed" ) .and. ( trim(adjustl(This%jInMethod)) .eq. "Fixed" ) ) then                 ! Initial v and j were provided
    if (i_Debug_Loc) call Logger%Write( "Initial v AND j were provided" )
    Found       =   .False.
    do i = 1,This%NStates
      associate( State => This%States(i) )
        if ( This%vIn /= State%vqn .or. This%jIn /= State%jqn ) cycle                   ! If current state does not match the input state, going to the next state
        This%StateIn =   i
        Found        =   .True.
        exit
      end associate
    end do
    if (.Not.Found) then
      call Logger%Write( "This%vIn = ", This%vIn, "This%jIn = ", This%jIn )
      call Error( "Error in initial state data: Data not found" )
    end if
! ==============================================================================================================


! ==============================================================================================================
!   CASE WHEN INPUT VALUES OF V AND/OR J ARE NOT PROVIDED
! ==============================================================================================================
  else
    if (i_Debug_Loc) call Logger%Write( "Either j and v both not fixed OR v fixed and j not" )
    
    if ( trim(adjustl(This%vInMethod)) .eq. "Uniform" ) then
      if (i_Debug_Loc) call Logger%Write( "Initial (v,j) or (j) selected from a Uniform Distribution")
      
      if (i_Debug_Loc) call Logger%Write( "-> Calling This%SetUniformDistribution" )
      call This%SetUniformDistribution( iState )
    
    elseif ( ( trim(adjustl(This%vInMethod)) .eq. "Boltzmann" ) .or. ( trim(adjustl(This%vInMethod)) .eq. "MostProbable" ) ) then
      if (i_Debug_Loc) call Logger%Write( "Initial (v,j) or (j) selected from a Boltzmann Distribution @ T [K] = ", This%Tint )
      
      Beta =   Factor / abs(This%Tint)
      
      if (i_Debug_Loc) call Logger%Write( "-> Calling This%SetBoltmannDistribution" )
      call This%SetBoltmannDistribution( Beta, iState )
      
      if ( trim(adjustl(This%vInMethod)) .eq. "MostProbable" ) then
        if (i_Debug_Loc) call Logger%Write( "Initial (v,j) or (j) fixed to the most probable state of the Boltzmann Distribution @ T [K] = ", This%Tint )
        
        This%vIn     =   This%States(iState)%vqn
        This%jIn     =   This%States(iState)%jqn
        This%StateIn =   iState
        if (i_Debug_Loc) call Logger%Write( "Most probable intial state: v = ", This%vIn, "j = ", This%jIn, " corresponding to the ", This%StateIn, "-th State." )
      
      end if
      
    elseif ( trim(adjustl(This%vInMethod)) .eq. "ArbitratyDistribution" ) then
      if (i_Debug_Loc) call Logger%Write( "Initial (v,j) or (j) sampled using the probabilities from unit", This%UnitArbDist )
  
      if (i_Debug_Loc) call Logger%Write( "-> Calling This%SetArbitraryDistribution" )
      call This%SetArbitraryDistribution()
    
    end if
    
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetUniformDistribution( This, iState, i_Debug )

  class(LevelsContainer_Type)               ,intent(inout)  ::    This
  integer                                   ,intent(out)    ::    iState                                          ! Index of the state having the maximum probability
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    i                                               ! Index of states
  real(rkp)                                                 ::    Qtot                                            ! Total partition function
  real(rkp)                                                 ::    Smax                                            ! Maximum contribution to the total partition function
  real(rkp)                                                 ::    xnorm
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetUniformDistribution")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


  Qtot          =   Zero
  iState        =   0
  Smax          =   Zero

  do i = 1,This%NStates
    associate( State => This%States(i) )
      State%rlim  =   real(2*State%jqn+1,kind=rkp)
      Qtot        =   Qtot + State%rlim
      if ( State%rlim > Smax ) then
        iState    =   i
        Smax      =   State%rlim
      end if
    end associate
  end do

  if (i_Debug_Global) then
    call Logger%Write( " Partition function: Qtot = ", Qtot) 
    call Logger%Write( " State of max. probability: iState = ", iState) 
    call Logger%Write( " Contributionof the max. prob. state: Smax = ", Smax)
  end if

  xnorm     =   One / Qtot
  do i = 1,This%NStates
    This%States(i)%rlim =   This%States(i)%rlim * xnorm
  end do
  do i = 2,This%NStates
    This%States(i)%rlim =   This%States(i-1)%rlim + This%States(i)%rlim                                           
  end do

  if (i_Debug_Loc) call Logger%Exiting( Writing=.False. )

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetBoltmannDistribution( This, Beta, iState, i_Debug)
! rlim(i) = [ (2*J(i)+1)*exp(-epsilon(i)/(KT)) + (2*J(i-1)+1)*exp(-epsilon(i-1)/(KT)) ] / sum_k{ (2*J(k)+1)*exp(-epsilon(k)/(KT)) }

  class(LevelsContainer_Type)               ,intent(inout)  ::    This
  real(rkp)                                 ,intent(in)     ::    Beta                                            ! Beta is 1/(Kb*t)
  integer                                   ,intent(out)    ::    iState                                          ! Index of the state having the maximum probability
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    i                                               ! Index of states
  real(rkp)                                                 ::    Qtot                                            ! Total partition function
  real(rkp)                                                 ::    Smax                                            ! Maximum contribution to the total partition function
  real(rkp)                                                 ::    xnorm
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetBoltmannDistribution")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc)  call Logger%Write(" Case of a heteronuclear diatomic species")


  Qtot          =   Zero
  iState        =   0
  Smax          =   Zero

  do i = 1,This%NStates
    associate( State => This%States(i) )
      State%rlim  =   real(2*State%jqn+1,kind=rkp) * exp( - Beta * State%eint )
      Qtot        =   Qtot + State%rlim
      if ( State%rlim > Smax ) then
        iState    =   i
        Smax      =   State%rlim
      end if
    end associate
  end do

  if (i_Debug_Loc)  then
    call Logger%Write(" Partition function: Qtot = ", Qtot)
    call Logger%Write(" State of max. probability: iState = ", iState)
    call Logger%Write(" Contributionof the max. prob. state: Smax = ", Smax)
  end if

  xnorm     =   One / Qtot
  do i = 1,This%NStates
    This%States(i)%rlim =   This%States(i)%rlim * xnorm
  end do
  do i = 2,This%NStates
    This%States(i)%rlim =   This%States(i-1)%rlim + This%States(i)%rlim                                           
  end do

  if (i_Debug_Loc) call Logger%Exiting( Writing=.False. )

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetArbitraryDistribution( This, i_Debug )
! This procedure reads the probabilities for sampling initial states

  class(LevelsContainer_Type)               ,intent(inout)  ::    This
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Unit
  integer                                                   ::    Status
  integer                                                   ::    i                                               ! Index of states
  integer                                                   ::    v                                               ! Vibrational quantum number read in the probability file
  integer                                                   ::    j                                               ! Rotational quantum number read in the probability file
  real(rkp)                                                 ::    Probability                                     ! Probability of a given state read in the probability file
  real(rkp)                                                 ::    Qtot                                            ! Total partition function
  real(rkp)                                                 ::    xnorm
  character(:)  ,allocatable                                ::    FileName
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetArbitraryDistribution")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  This%States(:)%rlim   =   Zero                                                                                  ! Initializing


! ==============================================================================================================
!   OPENING THE FILE CONTAINING THE STATE PROBABILITIES
! ==============================================================================================================
  allocate( FileName , source = 'fort.31')
  if (i_Debug_Loc) call Logger%Write(" Opening the file containing the state probabilities")
  if (i_Debug_Loc) call Logger%Write("  -> FileName = ", FileName)
  open( NewUnit=Unit, File=FileName, iostat=Status)
  if (Status/=0) call Error( "Error opening file '" // FileName // "' for state probabilities" )
! ==============================================================================================================


  Read_Loop: do
              read( Unit, *, iostat=Status ) v, j, Probability
              if (Status<0) exit                                                                                            ! End-of-file
              if (Status>0) call Error( "Error reading file '" // "fort.31" // "' for state data" )
              do i = 1,This%NStates
                if ( (This%States(i)%vqn/=v) .or. (This%States(i)%jqn/=j) ) cycle
                This%States(i)%rlim   =   Probability
                cycle Read_Loop
              end do
            end do Read_Loop

  Qtot      =   Zero
  do i = 1,This%NStates
    if ( This%States(i)%rlim > Zero ) then
      Qtot  =   Qtot + This%States(i)%rlim
    else
      if (i_Debug_Loc) call Logger%Write(" State with vj = ", This%States(i)%vqn, " not found in file ", FileName)
    end if
  end do

  xnorm     =   One / Qtot
  do i = 1,This%NStates
    This%States(i)%rlim   =   This%States(i)%rlim * xnorm
  end do
  do i = 2,This%NStates
    This%States(i)%rlim   =   This%States(i-1)%rlim + This%States(i)%rlim
  end do

  if (i_Debug_Loc) call Logger%Exiting( Writing=.False. )
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


Subroutine Compute_PartitionRatios( Input, iMol, LevelsContainer_Orig, LevelsContainer_Cut, i_Debug )

  use Input_Class             ,only:  Input_Type

  Type(Input_Type)                        ,intent(in)     ::    Input
  integer                                 ,intent(in)     ::    iMol
  Type(LevelsContainer_Type)              ,intent(inout)  ::    LevelsContainer_Orig
  Type(LevelsContainer_Type)              ,intent(inout)  ::    LevelsContainer_Cut
  logical                       ,optional ,intent(in)     ::    i_Debug
  
  Type(LevelsContainer_Type)                              ::    BLevels
  Type(LevelsContainer_Type)                              ::    QBLevels
  
  integer                                                 ::    Status
  integer                                                 ::    iLevels  
  integer                                                 ::    NQBound
  integer                                                 ::    NBound
  logical                                                 ::    i_Debug_Loc
  

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Compute_PartitionRatios")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


  ! Refering the Levels Energies to the Molecule Eenergy Minimum
  if (i_Debug_Loc) call Logger%Write( "Refering the Energies of the Original Molecule Levels and of the Cut Molecule ones to the Original Molecule Lowest Energy Level: ", LevelsContainer_Orig%Q )  
  LevelsContainer_Orig%States%einteV_scaled = LevelsContainer_Orig%States%einteV - LevelsContainer_Orig%mineinteV
  LevelsContainer_Cut%States%einteV_scaled  = LevelsContainer_Cut%States%einteV  - LevelsContainer_Orig%mineinteV


  ! Choosing the Model for Computing "Original" Molecule and "Cut" Molecule Partition Functions, based on the QCT States Sampling Method
  if (trim(adjustl(Input%vInMethod)) .eq. "Uniform") then
    if (i_Debug_Loc) call Logger%Write( "Uniform Model for Internal Modes Description")   
    
    LevelsContainer_Orig%Q = Zero
    do iLevels = 1,LevelsContainer_Orig%NStates
      LevelsContainer_Orig%Q = LevelsContainer_Orig%Q + LevelsContainer_Orig%States(iLevels)%g
    end do
    if (i_Debug_Loc) call Logger%Write( "Partition Function of the Original States Molecule:         LevelsContainer_Orig%Q = ", LevelsContainer_Orig%Q )   
    
    LevelsContainer_Cut%Q = Zero
    do iLevels = 1,LevelsContainer_Cut%NStates
      LevelsContainer_Cut%Q = LevelsContainer_Cut%Q + LevelsContainer_Cut%States(iLevels)%g
    end do
    if (i_Debug_Loc) call Logger%Write( "Partition Function of the Cut States Molecule:         LevelsContainer_Cut%Q = ", LevelsContainer_Cut%Q )    

  elseif (trim(adjustl(Input%vInMethod)) .eq. "Boltzmann") then
    if (i_Debug_Loc) call Logger%Write( "Boltzmann Model for Internal Modes Description")   
    
    LevelsContainer_Orig%Q = Zero
    do iLevels = 1,LevelsContainer_Orig%NStates
      LevelsContainer_Orig%Q = LevelsContainer_Orig%Q + LevelsContainer_Orig%States(iLevels)%g * dexp( - (LevelsContainer_Orig%States(iLevels)%einteV_scaled * Ue) / (UKb * Input%Tint) )
    end do
    if (i_Debug_Loc) call Logger%Write( "Partition Function of the Original States Molecule:         LevelsContainer_Orig%Q = ", LevelsContainer_Orig%Q )   
    
    LevelsContainer_Cut%Q = Zero
    do iLevels = 1,LevelsContainer_Cut%NStates
      LevelsContainer_Cut%Q = LevelsContainer_Cut%Q + LevelsContainer_Cut%States(iLevels)%g * dexp( - (LevelsContainer_Cut%States(iLevels)%eintev_scaled * Ue) / (UKb * Input%Tint) )
    end do
    if (i_Debug_Loc) call Logger%Write( "Partition Function of the Cut States Molecule:         LevelsContainer_Cut%Q = ", LevelsContainer_Cut%Q )     
  
  else
    call Error( "Error: Not possible to compute Partition Function: Internal Modes Model not Specified!" )
  end if
  
  
  ! Separating Quasi-Bound Levels from Bound ones and computing their degeneracies and energies
  if (i_Debug_Loc) call Logger%Write( "Counting Nb of Bound and Quasi-Bound States" )   
  LevelsContainer_Orig%NQBound = 0
  LevelsContainer_Orig%NBound  = 0

  do iLevels = 1,LevelsContainer_Orig%NStates
    if (LevelsContainer_Orig%States(iLevels)%eint > Zero) then
      LevelsContainer_Orig%NQBound = LevelsContainer_Orig%NQBound + 1
    else
      LevelsContainer_Orig%NBound = LevelsContainer_Orig%NBound + 1
    end if
  end do
  if (i_Debug_Loc) call Logger%Write( "Nb of Bound States:               LevelsContainer_Orig%NBound  = ", LevelsContainer_Orig%NBound  )   
  if (i_Debug_Loc) call Logger%Write( "Nb of Quasi-Bound States:         LevelsContainer_Orig%NQBound = ", LevelsContainer_Orig%NQBound )   
  
  allocate(  BLevels%States(LevelsContainer_Orig%NBound) , Stat=Status )
  if (Status/=0) call Error( "Compute_PartitionRatios: Error allocating BLevels%States" )
  if (i_Debug_Loc) call Logger%Write( "BLevels%States, Bound Levels Container, allocated " ) 
  BLevels%NStates = LevelsContainer_Orig%NBound
  
  allocate( QBLevels%States(LevelsContainer_Orig%NQBound) , Stat=Status )
  if (Status/=0) call Error( "Compute_PartitionRatios: Error allocating QBLevels%States" )
  if (i_Debug_Loc) call Logger%Write( "QBLevels%States, Quasi-Bound Levels Container, allocated " ) 
  QBLevels%NStates = LevelsContainer_Orig%NQBound
  
  NQBound = 0
  NBound  = 0
  do iLevels = 1,LevelsContainer_Orig%NStates
    if (LevelsContainer_Orig%States(iLevels)%eint > Zero) then
      NQBound = NQBound + 1
      QBLevels%States(NQBound)%to_level      = iLevels
      QBLevels%States(NQBound)%g             = LevelsContainer_Orig%States(iLevels)%g
      QBLevels%States(NQBound)%einteV_scaled = LevelsContainer_Orig%States(iLevels)%einteV_scaled
    else
      NBound = NBound + 1
      BLevels%States(NBound)%to_level        = iLevels
      BLevels%States(NBound)%g               = LevelsContainer_Orig%States(iLevels)%g
      BLevels%States(NBound)%einteV_scaled   = LevelsContainer_Orig%States(iLevels)%einteV_scaled
    end if
  end do
  if (i_Debug_Loc) call Logger%Write( "Nb of Bound States:         LevelsContainer_Orig%NBound = ", NBound )   
  if (i_Debug_Loc) call Logger%Write( "Nb of Quasi-Bound States:         LevelsContainer_Orig%NQBound = ", NQBound ) 
  
  
  ! Choosing the Model for Computing Bound and Quasi-Bound Partition Functions, based on the QCT States Sampling Method
  if (trim(adjustl(Input%vInMethod)) .eq. "Uniform") then
    if (i_Debug_Loc) call Logger%Write( "Uniform Model for Internal Modes Description")   
    
    BLevels%Q = Zero
    do iLevels = 1,BLevels%NStates
      BLevels%Q = BLevels%Q + BLevels%States(iLevels)%g 
    end do
    if (i_Debug_Loc) call Logger%Write( "Bound Levels Degeneracy,            BLevels%Q = ", BLevels%Q ) 
    
    QBLevels%Q = Zero
    do iLevels = 1,QBLevels%NStates
      QBLevels%Q = QBLevels%Q + QBLevels%States(iLevels)%g
    end do
    if (i_Debug_Loc) call Logger%Write( "Quasi-Bound Levels Degeneracy,            QBLevels%Q = ", QBLevels%Q ) 
  
  elseif (trim(adjustl(Input%vInMethod)) .eq. "Boltzmann") then
    if (i_Debug_Loc) call Logger%Write( "Boltzmann Model for Internal Modes Description")   
    
    BLevels%Q = Zero
    do iLevels = 1,BLevels%NStates
      BLevels%Q = BLevels%Q + BLevels%States(iLevels)%g * dexp( - (BLevels%States(iLevels)%einteV_scaled * Ue) / (UKb * Input%Tint) )
    end do
    if (i_Debug_Loc) call Logger%Write( "Bound Levels Degeneracy,            BLevels%Q = ", BLevels%Q ) 
    
    QBLevels%Q = Zero
    do iLevels = 1,QBLevels%NStates
      QBLevels%Q = QBLevels%Q + QBLevels%States(iLevels)%g * dexp( - (QBLevels%States(iLevels)%einteV_scaled * Ue) / (UKb * Input%Tint) )
    end do
    if (i_Debug_Loc) call Logger%Write( "Quasi-Bound Levels Degeneracy,            QBLevels%Q = ", QBLevels%Q ) 
    
  else 
    call Error( "Error: Not possible to compute Partition Function: Internal Modes Model not Specified!" )
  end if


  ! Computing the Partition Functions Ratios
  LevelsContainer_Orig%CutRatio     = LevelsContainer_Cut%Q / LevelsContainer_Orig%Q
  if (Input%Ecut(iMol) < Zero) then
    LevelsContainer_Orig%BoundRatio   =  BLevels%Q / LevelsContainer_Orig%Q + LevelsContainer_Orig%CutRatio - One
    LevelsContainer_Orig%QBoundRatio  =  QBLevels%Q / LevelsContainer_Orig%Q
  else
    LevelsContainer_Orig%BoundRatio   = Zero
    LevelsContainer_Orig%QBoundRatio  = ( LevelsContainer_Cut%Q - LevelsContainer_Orig%Q )
  end if
  if (i_Debug_Loc) call Logger%Write( "Partition Function Ratio: Levels of Cut Molecule / Levels of Original Molecule:                      LevelsContainer_Orig%CutRatio    =  ", LevelsContainer_Orig%CutRatio    ) 
  if (i_Debug_Loc) call Logger%Write( "Partition Function Ratio: Q-Bound Levels of Cut Molecule / All Levels of Original Molecule:          LevelsContainer_Orig%QBoundRatio =  ", LevelsContainer_Orig%QBoundRatio ) 
  if (i_Debug_Loc) call Logger%Write( "Partition Function Ratio: Bound Levels of Cut Molecule / All Levels of Original Molecule:            LevelsContainer_Orig%BoundRatio  =  ", LevelsContainer_Orig%BoundRatio  ) 
  
  if (i_Debug_Loc) call Logger%Exiting( Writing=.False. )
  
End Subroutine


End Module