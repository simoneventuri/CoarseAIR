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
Module Processes_Class

#include "../qct.inc"

  use Parameters_Module       ,only:  rkp, Zero, One, Six, Ue, UKb
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error
  use Process_Class           ,only:  Process_Type

  implicit none

  private
  public  ::    Processes_Type
  

  Type    ,abstract                                         ::    Processes_Type
    logical                                                 ::    Initialized         !< Indicator whether the object is initialized

    integer                ,dimension(:)      ,allocatable  ::    IdxVec, IdxVecSorted

    logical                                                 ::    MergeExchToInelFlg = .False. ! If True, merging Equal Exchanges together and to Inelastic Processes
                                                                                               !
    logical                                                 ::    MergeExchsFlg      = .False. ! If True, merging Equal Exchanges together
                                                                                               ! 
    
    integer                                                 ::    NProcType               !
                                                                                          ! 
    integer                ,dimension(:)      ,allocatable  ::    ProcTypeVec             !
                                                                                          !                                                               
    integer                ,dimension(3)                    ::    ExcTypeVec              ! Type of Exchange (ex: N2+N  -> ExcTypeVec=[1,1,1]; CO+O  -> ExcTypeVec=[1,1,2]; 
                                                                                          !                       N2+NO -> ExcTypeVec=[1,1,1]; N2+O2 -> ExcTypeVec=[1,2,2]).
                                                                                          !
    integer                ,dimension(3,2)                  ::    NProc_iPOpp       = 0   ! For 3AtomSyst, NProc_iPOpp = [ (NbProc)_Pair1, (NbProc)_Pair2, (NbProc)_Pair3;              0,              0,               0]    
                                                                                          ! For 4AtomSyst, NProc_iPOpp = [ (NbProc)_Pair1, (NbProc)_Pair2, (NbProc)_Pair3; (NbProc)_Pair6, (NbProc)_Pair5, (NbProc)_Pair4 ].
                                                                                          !   where (NbProc)_Pair = NLevels_Pair + 1 (1 is for dissociation)
                                                                                          !
    integer                ,dimension(0:3)                  ::    NProc_iP          = 0   ! For 3AtomSyst, NProc_iP = [0, NProc_iPOpp(1,1), NProc_iPOpp(2,1), NProc_iPOpp(3,1)];
                                                                                          ! For 4AtomSyst, NProc_iP = [0, NProc_iPOpp(1,1)*NProc_iPOpp(1,2)+1, NProc_iPOpp(2,1)*NProc_iPOpp(2,2)+1, NProc_iPOpp(3,1)*NProc_iPOpp(3,2)+1].
                                                                                          !
    integer                                                 ::    NProc_Tot         = 0   ! NProc_Tot = sum(NProc_iP)
                                                                                          !  
    integer                                                 ::    NProc_Tot_Unique  = 0   ! NProc_Tot = sum(NProc_iP)
                                                                                          !
    integer                ,dimension(:)      ,allocatable  ::    Proc_To_LineVec         ! Mapping Line of Variables (coming from CrossSec File / Rates File) to Correspondent Process
                                                                                          !
    Type(Process_Type)     ,dimension(:)      ,allocatable  ::    ProcessesVec            ! Vector of Processes
    Type(Process_Type)     ,dimension(:)      ,allocatable  ::    ProcessesVecTemp        ! Vector of Processes
    Type(Process_Type)     ,dimension(:)      ,allocatable  ::    ProcessesVecCleaned     ! Vector of Processes
                                                                                          !
    integer                                                 ::    NProc             = 0   !                                                                                           
    integer                                                 ::    NProc_Cleaned     = 0   ! 
    real(rkp)              ,dimension(:,:)    ,allocatable  ::    OvProcRates

    integer                                                 ::    Status            = -1
    character(:)                              ,allocatable  ::    OutputDir
    character(:)                              ,allocatable  ::    LevelOutputDir
    logical                                                 ::    FirstIssueFlg
    integer                                                 ::    UnitIssues
    integer                                                 ::    NTraj
    character(17)                                           ::    System
    logical                                                 ::    StochPESFlg
    integer                                                 ::    NPESs
    character(17)                                           ::    PES_Name
    integer                                                 ::    PESoI
    character(:)                              ,allocatable  ::    PESoI_char   
    character(17)                                           ::    IniMolecules
    integer                ,dimension(:)      ,allocatable  ::    InBins
    character(17)                                           ::    InBinsChar
    character(17)                                           ::    InBinsCharName
    integer                                                 ::    InProc
    character(:)                              ,allocatable  ::    InProcChar
    integer                                                 ::    NTTra
    real(rkp)              ,dimension(:)      ,allocatable  ::    TTra
    character(20)          ,dimension(:)      ,allocatable  ::    TTraChar
    integer                                                 ::    NTInt
    real(rkp)              ,dimension(:)      ,allocatable  ::    TInt
    character(20)          ,dimension(:)      ,allocatable  ::    TIntChar
    real(rkp)              ,dimension(:)      ,allocatable  ::    QRatio
    character(37)                                           ::    QRatioChar
    !real(rkp)              ,dimension(:,:,:)  ,allocatable  ::    QRatio
    !character(37)          ,dimension(:,:)    ,allocatable  ::    QRatioChar

  contains
    private
    procedure              ,public                          ::    Initialize                   =>    Initialize_Processes
    procedure              ,public                          ::    Mask4Excahge                 =>    Mask4Excahge_Processes
    procedure              ,public                          ::    ConstructVecOfProcs          =>    ConstructVecOfProcs_Processes
    procedure              ,public                          ::    FindingFinalLevel            =>    FindingFinalLevel_Processes
    procedure              ,public                          ::    Convert_CrossSect_To_Rates   =>    Convert_CrossSect_To_Rates_Processes
    procedure              ,public                          ::    WritingRates
    procedure              ,public                          ::    WritingRates_Binary
    procedure              ,public                          ::    WritingIssue
    procedure              ,public                          ::    ReadingRates
    procedure              ,public                          ::    ReadingRates_Binary
  End Type

  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains


Subroutine Initialize_Processes( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Processes_Type)                             ,intent(out)   :: This
  Type(Input_Type)                                  ,intent(in)    :: Input
  Type(Collision_Type)                              ,intent(inout) :: Collision
  logical                                 ,optional ,intent(in)    :: i_Debug
 
  logical                                                          :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeProcesses_Processes" )
  !i_Debug_Loc   =     Logger%On()
   
  This%Initialized  =   .True.

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


Subroutine Mask4Excahge_Processes( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Processes_Type)                             ,intent(inout)  :: This
  Type(Input_Type)                                  ,intent(in)     :: Input
  Type(Collision_Type)                              ,intent(inout)  :: Collision
  logical                                 ,optional ,intent(in)     :: i_Debug
  
  integer                                                           :: iP, jP
  integer                                                           :: Status
  integer                                                           :: iMol1, iMol, jMol
  logical                                                           :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Mask4Excahge_Processes" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


Subroutine ConstructVecOfProcs_Processes( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Processes_Type)                             ,intent(inout) :: This
  Type(Input_Type)                                  ,intent(in)    :: Input
  Type(Collision_Type)                              ,intent(in)    :: Collision
  logical                                 ,optional ,intent(in)    :: i_Debug
 
  logical                                                          :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ConstructVecOfProcs_Processes" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


Subroutine FindingFinalLevel_Processes( This, Input, Collision, vqn, jqn, Arr, Name, ProcType, ExcType, Pairs, iLevelFin, iLevelFinChar, Idx, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type
  use Process_Class               ,only: Process_Type

  class(Processes_Type)                             ,intent(in)  :: This
  Type(Input_Type)                                  ,intent(in)  :: Input
  Type(Collision_Type)                              ,intent(in)  :: Collision
  integer       ,dimension(:)                       ,intent(in)  :: vqn
  integer       ,dimension(:)                       ,intent(in)  :: jqn
  integer                                           ,intent(in)  :: Arr
  character(:)  ,allocatable                        ,intent(out) :: Name 
  integer                                           ,intent(out) :: ProcType
  integer                                           ,intent(out) :: ExcType
  integer       ,dimension(:)                       ,intent(out) :: Pairs
  integer       ,dimension(:)                       ,intent(out) :: iLevelFin
  character(6)  ,dimension(:)                       ,intent(out) :: iLevelFinChar  
  integer                                           ,intent(out) :: Idx
  logical                                 ,optional ,intent(in)  :: i_Debug
 
  logical                                                        :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "FindingFinalLevel_Processes" )
  !i_Debug_Loc   =     Logger%On()

  Name          = ''
  ProcType      = 0
  ExcType       = 0
  Pairs         = 0
  iLevelFin     = 0
  iLevelFinChar = '      '
  Idx           = 0
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


! ==============================================================================================================
!   READING INITIAL AND FINAL CONDITIONS AND CROSS SECTIONS FROM statistics.out
! ==============================================================================================================
Subroutine Convert_CrossSect_To_Rates_Processes( This, Input, Collision, Velocity, i_Debug, i_Debug_Deep)

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero
  
  class(Processes_Type)                                 ,intent(inout) ::    This
  type(Input_Type)                                      ,intent(in)    ::    Input
  type(Collision_Type)                                  ,intent(inout) ::    Collision
  real(rkp)                  ,dimension(:)              ,intent(in)    ::    Velocity
  logical                                     ,optional ,intent(in)    ::    i_Debug
  logical                                     ,optional ,intent(in)    ::    i_Debug_Deep


  logical                                                   ::    i_Debug_Loc
  logical                                                   ::    i_Debug_Loc_Deep


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Convert_CrossSect_To_Rates_Processes" )
  i_Debug_Loc_Deep = i_Debug_Global; if ( present(i_Debug_Deep) )i_Debug_Loc_Deep = i_Debug_Deep
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   WRITING BINS RATES
! ==============================================================================================================
!______________________________________________________________________________________________________________!
Subroutine WritingRates( This, iTTra, iTInt, i_Debug )

  use Parameters_Module     ,only:  Zero

  class(Processes_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTTra
  integer                                   ,intent(in)     ::    iTInt
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iP, jP
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)                                 ,allocatable ::    FileName
  character(:)                                 ,allocatable ::    TsString
  character(17)                                             ::    PES_Name
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "WritingRates" )
  !i_Debug_Loc   =     Logger%On()  

  allocate(TsString, source = trim(adjustl( 'T_' // trim(adjustl(This%TTraChar(iTTra))) // '_' // trim(adjustl(This%TIntChar(iTInt))) // '/' )) )

  call system('mkdir -p ' // trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/')
  call system('mkdir -p ' // trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString)

  if (This%PESoI == 0) then
    !FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // trim(adjustl(This%InProcChar)) // '.csv' ))
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/' // adjustl(trim(This%InBinsCharName)) // '.csv' ))


    if (This%NPESs == 1) then
      PES_Name = adjustl(trim(This%PES_Name))
    else
      PES_Name = '           Merged'
    end if
  else 
    !FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // trim(adjustl(This%InProcChar)) // '.csv.' // This%PESoI_char ))
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/' // adjustl(trim(This%InBinsCharName)) // '.csv.' // This%PESoI_char ))
    if (This%StochPESFlg) then
      PES_Name = adjustl(trim('SPES_' // This%PESoI_char))
    else
      PES_Name = adjustl(trim(This%PES_Name))
    end if
  end if
  if ( i_Debug_Loc ) call Logger%Write( "Writing File: ", FileName )
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error writing the ASCI data file for Rates: " // FileName  ) 
  
    write(Unit,'(A)')                                          '#           System |               PES |     In. Molecules |      No Processes |   No Trajectories |'
    write(Unit,'(A1,A17, A3,A17, A3,A17, A3,I17, A3,I17, A2)') '#', adjustr(This%System), ' | ', adjustr(PES_Name), ' | ', adjustr(This%IniMolecules), ' | ', This%NProc_Tot, ' | ', This%NTraj, ' |'

    write(Unit,'(A)') '#   In. Lev.s/Bins |                   Part. Func.s Ratios | Trans Temperature |   Int Temperature |'
    write(Unit,'(A1,A17, A3,A37, A3,es17.10, A3,es17.10, A2)') '#', adjustr(This%InBinsChar), ' | ', This%QRatioChar, ' | ', This%TTra(iTTra), ' | ', This%TInt(iTInt), ' |'
    
    write(Unit,'(A)') '#  Process,            Rate,         Rate SD'
    do iP = 1,This%NProc_Cleaned
      jP = This%IdxVecSorted(iP)
      write(Unit,'(I10,A1,es16.10,A1,es16.10)') This%ProcessesVecCleaned(jP)%Idx, ',', This%ProcessesVecCleaned(jP)%Temperature(iTTra)%Rate, ',', sqrt( This%ProcessesVecCleaned(jP)%Temperature(iTTra)%RateSD )
    end do
  
  close(Unit)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   WRITING LEVEL/BIN RATES IN ASCI FORMAT
! ==============================================================================================================
!______________________________________________________________________________________________________________!
Subroutine WritingRates_Binary( This, iTTra, iTInt, i_Debug )

  use Parameters_Module     ,only:  Zero

  class(Processes_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTTra
  integer                                   ,intent(in)     ::    iTInt
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iP, jP
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)                                 ,allocatable ::    FileName
  character(:)                                 ,allocatable ::    TsString
  character(17)                                             ::    PES_Name
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "WritingRates_Binary" )
  !i_Debug_Loc   =     Logger%On()  

  allocate(TsString, source = trim(adjustl( 'T_' // trim(adjustl(This%TTraChar(iTTra))) // '_' // trim(adjustl(This%TIntChar(iTInt))) // '/' )) )

  call system('mkdir -p ' // trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/')
  call system('mkdir -p ' // trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString)

  if (This%PESoI == 0) then
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // This%InProcChar // '.bin' ))
  else 
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // This%InProcChar // '.bin.' // This%PESoI_char ))
  end if
  if ( i_Debug_Loc ) call Logger%Write( "Writing File: ", FileName )
  open( NewUnit=Unit, File=FileName, Action='WRITE', access="Stream", form="Unformatted", iostat=Status )
  if (Status/=0) call Error( "Error writing the binary data file for Rates: " // FileName  ) 

    write(Unit) int(This%NTraj, rkp)
    write(Unit) int(This%NProc_Cleaned, rkp)
    do iP = 1,This%NProc_Cleaned
      jP = This%IdxVecSorted(iP)

      write(Unit) int(This%ProcessesVecCleaned(jP)%Idx, rkp)
      write(Unit) This%ProcessesVecCleaned(jP)%Temperature(iTTra)%Rate
      write(Unit) sqrt( This%ProcessesVecCleaned(jP)%Temperature(iTTra)%RateSD )
    end do
  
  close(Unit)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



! ==============================================================================================================
!   WRITING LEVEL/BIN RATES IN BINARY FORMAT
! ==============================================================================================================
!______________________________________________________________________________________________________________!
Subroutine WritingIssue( This, vqnIn, jqnIn, ArrIn, IssueIn, vqnFin, jqnFin, ArrFin, IssueFin, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero
  
  class(Processes_Type)                      ,intent(inout) ::     This
  integer   ,dimension(:)                    ,intent(in)    ::    vqnIn
  integer   ,dimension(:)                    ,intent(in)    ::    jqnIn
  integer                                    ,intent(in)    ::    ArrIn
  integer                                    ,intent(in)    ::  IssueIn
  integer   ,dimension(:)                    ,intent(in)    ::   vqnFin
  integer   ,dimension(:)                    ,intent(in)    ::   jqnFin
  integer                                    ,intent(in)    ::   ArrFin
  integer                                    ,intent(in)    :: IssueFin
  logical                          ,optional ,intent(in)    ::  i_Debug

  integer                                                ::    iP, jP
  integer                                                ::    Status
  integer                                                ::    Unit
  character(:)                              ,allocatable ::    FileName
  character(:)                              ,allocatable ::    TsString
  character(17)                                          ::    PES_Name
  logical                                                ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "WritingIssue" )
  !i_Debug_Loc   =     Logger%On()  

  if (This%FirstIssueFlg) then
    FileName = adjustl(trim( trim(adjustl(This%LevelOutputDir)) // '/StatIssues.csv' ))
    if ( i_Debug_Loc ) call Logger%Write( "Writing File: ", FileName )
    open( File=FileName, NewUnit=This%UnitIssues, status='unknown', access='append', iostat=Status )
    if (Status/=0) call Error( "Error writing the binary data file for Rates: " // FileName  ) 
    if (size(vqnIn,1) == 1) then
      write(This%UnitIssues, '(A)') '#    vIn(1),    jIn(1),     ArrIn,   IssueIn,   vFin(1),   jFin(1),    ArrFin,  IssueFin'
    else
      write(This%UnitIssues, '(A)') '#    vIn(1),    jIn(1),    vIn(2),    jIn(2),     ArrIn,   IssueIn,   vFin(1),   jFin(1),   vFin(2),   jFin(2),    ArrFin,  IssueFin'
    end if
    This%FirstIssueFlg = .False.
  end if
    
  if (size(vqnIn,1) == 1) then
    write(This%UnitIssues, 1) vqnIn(1), ',', jqnIn(1), ',', ArrIn, ',', IssueIn, ',', vqnFin(1), ',', jqnFin(1), ',', ArrFin, ',', IssueFin
  else 
    write(This%UnitIssues, 1) vqnIn(1), ',', jqnIn(1), ',', vqnIn(2), ',', jqnIn(2), ',', ArrIn, ',', IssueIn, ',', vqnFin(1), ',', jqnFin(1), ',', vqnFin(2), ',', jqnFin(2), ',', ArrFin, ',', IssueFin
  end if

  1 Format(X, I10, *(A, I10))

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   READING LEVEL/BIN RATES IN ASCI FORMAT
! ==============================================================================================================
!______________________________________________________________________________________________________________!
Subroutine ReadingRates( This, iTTra, iTInt, i_Debug )

  use Parameters_Module     ,only:  Zero

  class(Processes_Type)                     ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iTTra
  integer                                   ,intent(in)     ::    iTInt
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iTemp, jP
  real(rkp)                                                 ::    Temp1, Temp2
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)                                 ,allocatable ::    FileName
  character(:)                                 ,allocatable ::    TsString
  character(17)                                             ::    PES_Name
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "ReadingRates" )
  !i_Debug_Loc   =     Logger%On()  

  allocate(TsString, source = trim(adjustl( 'T_' // trim(adjustl(This%TTraChar(iTTra))) // '_' // trim(adjustl(This%TIntChar(iTInt))) // '/' )) )

  if (This%PESoI == 0) then
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // trim(adjustl(This%InProcChar)) // '.csv' ))
    if (This%NPESs == 1) then
      PES_Name = adjustl(trim(This%PES_Name))
    else
      PES_Name = '           Merged'
    end if
  else 
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // trim(adjustl(This%InProcChar)) // '.csv.' // This%PESoI_char ))
    PES_Name = adjustl(trim('SPES_' // This%PESoI_char))
  end if
  if ( i_Debug_Loc ) call Logger%Write( "Reading File: ", FileName )
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error reading the ASCI data file for Rates: " // FileName  ) 
    read(Unit,*,iostat=Status)
    read(Unit,*,iostat=Status)
    read(Unit,*,iostat=Status)
    read(Unit,*,iostat=Status)    
    read(Unit,*,iostat=Status)
    jP=1
    do 
      read(Unit, *, iostat=Status) iTemp, Temp1, Temp2
      if (Status == 0) then
        This%ProcessesVec(jP)%Idx                       = iTemp
        This%ProcessesVec(jP)%Temperature(iTTra)%Rate   = Temp1
        This%ProcessesVec(jP)%Temperature(iTTra)%RateSD = Temp2**2
      else
        exit
      end if
      jP=jP+1
    end do
  close(Unit)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   READING LEVEL/BIN RATES IN BINARY FORMAT
! ==============================================================================================================
!______________________________________________________________________________________________________________!
Subroutine ReadingRates_Binary( This, iTTra, iTInt, i_Debug )

  use Parameters_Module     ,only:  Zero

  class(Processes_Type)                     ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iTTra
  integer                                   ,intent(in)     ::    iTInt
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iP, jP
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)                                 ,allocatable ::    FileName
  character(:)                                 ,allocatable ::    TsString
  character(17)                                             ::    PES_Name
  integer                                                   ::    POSTemp
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "ReadingRates_Binary" )
  !i_Debug_Loc   =     Logger%On()  

  allocate(TsString, source = trim(adjustl( 'T_' // trim(adjustl(This%TTraChar(iTTra))) // '_' // trim(adjustl(This%TIntChar(iTInt))) // '/' )) )

  if (This%PESoI == 0) then
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // This%InProcChar // '.bin' ))
  else 
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/Proc' // This%InProcChar // '.bin.' // This%PESoI_char ))
  end if
  if ( i_Debug_Loc ) call Logger%Write( "Reading File: ", FileName )
  open( NewUnit=Unit, File=FileName, Action='READ', access="Stream", form="Unformatted", iostat=Status )
  if (Status/=0) call Error( "Error reading the binary data file for Rates: " // FileName  ) 

    POSTemp = rkp
    read(Unit, POS=POSTemp) This%NTraj
    POSTemp = int(2*rkp)
    read(Unit, POS=POSTemp) This%NProc
    do iP = 1,This%NProc_Cleaned
      jP = This%IdxVecSorted(iP)

      POSTemp = POSTemp + rkp
      read(Unit, POS=POSTemp) This%ProcessesVec(jP)%Idx
      POSTemp = POSTemp + rkp
      read(Unit, POS=POSTemp) This%ProcessesVec(jP)%Temperature(iTTra)%Rate
      POSTemp = POSTemp + rkp
      read(Unit, POS=POSTemp) This%ProcessesVec(jP)%Temperature(iTTra)%RateSD
                              This%ProcessesVec(jP)%Temperature(iTTra)%RateSD = This%ProcessesVec(jP)%Temperature(iTTra)%RateSD**2
    end do
  
  close(Unit)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



! ! ==============================================================================================================
! !   READING LEVEL/BIN RATES IN BINARY FORMAT
! ! ==============================================================================================================
! !______________________________________________________________________________________________________________!
! Subroutine MergingPESs( ProcessesRead, ProcessesMerged, iTTra, iTInt, Idx, i_Debug )

!   use Parameters_Module     ,only:  Zero

!   class(Processes_Type)                     ,intent(in)     ::    ProcessesRead
!   class(Processes_Type)                     ,intent(inout)  ::    ProcessesMerged
!   integer                                   ,intent(in)     ::    iTTra
!   integer                                   ,intent(in)     ::    iTInt
!   integer                                   ,intent(in)     ::    Idx
!   logical                         ,optional ,intent(in)     ::    i_Debug

!   integer                                                   ::    iProc
!   logical                                                   ::    i_Debug_Loc

!   i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
!   if ( i_Debug_Loc ) call Logger%Entering( "MergingPESs" )
!   !i_Debug_Loc   =     Logger%On()  

!   do iProc=1,ProcessesRead%NProc
!     ProcessesMerged(Idx)%Temperature(iTTra)%Rate   = ProcessesMerged(Idx)%Temperature(iTTra)%Rate   + ProcessesRead(iProc)%Temperature(iTTra)%Rate
!     ProcessesMerged(Idx)%Temperature(iTTra)%RateSD = ProcessesMerged(Idx)%Temperature(iTTra)%RateSD + ProcessesRead(iProc)%Temperature(iTTra)%RateSD
!   end do

!   if ( i_Debug_Loc ) call Logger%Exiting

! End Subroutine
! !--------------------------------------------------------------------------------------------------------------------------------!


End Module