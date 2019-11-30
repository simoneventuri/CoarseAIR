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

    logical                                                 ::    MergeExchToInelFlg = .False. ! If True, merging Equal Exchanges together and to Inelastic Processes
                                                                                               !
    logical                                                 ::    MergeExchsFlg      = .False. ! If True, merging Equal Exchanges together
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
    integer                ,dimension(:)      ,allocatable  ::    ExchMask                ! Mask For Exchange

    integer                                                 ::    NProc_Cleaned     = 0   ! 

    integer                                                 ::    Status            = -1
    character(:)                              ,allocatable  ::    OutputDir
    integer                                                 ::    NTraj
    character(17)                                           ::    System
    integer                                                 ::    NPESs
    character(17)                                           ::    PES_Name
    integer                                                 ::    PESoI
    character(:)                              ,allocatable  ::    PESoI_char   
    character(17)                                           ::    IniMolecules
    integer                ,dimension(:)      ,allocatable  ::    InBins
    character(17)                                           ::    InBinsChar
    integer                                                 ::    InProc
    character(:)                              ,allocatable  ::    InProcChar
    integer                                                 ::    NTTra
    real(rkp)              ,dimension(:)      ,allocatable  ::    TTra
    character(20)          ,dimension(:)      ,allocatable  ::    TTraChar
    integer                                                 ::    NTInt
    real(rkp)              ,dimension(:)      ,allocatable  ::    TInt
    character(20)          ,dimension(:)      ,allocatable  ::    TIntChar
    real(rkp)              ,dimension(:,:,:)  ,allocatable  ::    QRatio
    character(37)          ,dimension(:,:)    ,allocatable  ::    QRatioChar

  contains
    private
    procedure              ,public                          ::    Initialize                   =>    Initialize_Processes
    procedure              ,public                          ::    Mask4Excahge                 =>    Mask4Excahge_Processes
    procedure              ,public                          ::    ConstructVecOfProcs          =>    ConstructVecOfProcs_Processes
    procedure              ,public                          ::    FindingFinalLevel            =>    FindingFinalLevel_Processes
    procedure              ,public                          ::    Convert_CrossSect_To_Rates   =>    Convert_CrossSect_To_Rates_Processes
    procedure              ,public                          ::    InProc_WritingRates
  End Type

  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains


Subroutine Initialize_Processes( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Processes_Type)                             ,intent(out)   :: This
  Type(Input_Type)                                  ,intent(in)    :: Input
  Type(Collision_Type)                              ,intent(in)    :: Collision
  logical                                 ,optional ,intent(in)    :: i_Debug
 
  logical                                                          :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeProcesses_Processes" )
  !i_Debug_Loc   =     Logger%On()
   
  This%Initialized  =   .True.

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


Subroutine Mask4Excahge_Processes( This, Input, i_Debug )

  use Input_Class                 ,only: Input_Type

  class(Processes_Type)                             ,intent(inout)  :: This
  Type(Input_Type)                                  ,intent(in)     :: Input
  logical                                 ,optional ,intent(in)     :: i_Debug
  
  integer                                                           :: iP
  integer                                                           :: Status
  integer                                                           :: iLevel, jLevel, kLevel
  logical                                                           :: i_Debug_Loc

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
  character(*)                                      ,intent(out) :: Name 
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
  type(Collision_Type)                                  ,intent(in)    ::    Collision
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
Subroutine InProc_WritingRates( This, iTTra, iTInt, i_Debug )

  use Parameters_Module     ,only:  Zero

  class(Processes_Type)                     ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    iTTra
  integer                                   ,intent(in)     ::    iTInt
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iP
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)                                 ,allocatable ::    FileName
  character(:)                                 ,allocatable ::    TsString
  character(17)                                             ::    PES_Name
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "InProc_WritingRates" )
  !i_Debug_Loc   =     Logger%On()  

  allocate(TsString, source = trim(adjustl( 'T_' // trim(adjustl(This%TTraChar(iTTra))) // '_' // trim(adjustl(This%TIntChar(iTInt))) // '/' )) )

  call system('mkdir -p ' // trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/')
  call system('mkdir -p ' // trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString)

  if (This%PESoI == 0) then
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/InProc' // This%InProcChar // '.dat' ))
    if (This%NPESs == 1) then
      PES_Name = adjustl(trim(This%PES_Name))
    else
      PES_Name = '           Merged'
    end if
  else 
    FileName = adjustl(trim( trim(adjustl(This%OutputDir)) // '/'// trim(adjustl(This%System)) // '/Rates/' // TsString // '/InProc' // This%InProcChar // '.dat.' // This%PESoI_char ))
    PES_Name = adjustl(trim('SPES_' // This%PESoI_char))
  end if
  if ( i_Debug_Loc ) call Logger%Write( "Writing File: ", FileName )
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  
    write(Unit,'(A)')                                          '#           System |               PES | Initial Molecules |      No Processes |   No Trajectories |'
    write(Unit,'(A1,A17, A3,A17, A3,A17, A3,I17, A3,I17, A2)') '$', adjustr(This%System), ' | ', adjustr(PES_Name), ' | ', adjustr(This%IniMolecules), ' | ', This%NProc_Tot, ' | ', This%NTraj, ' |'

    write(Unit,'(A)') '$     Initial Bins |                       PartFunc Ratios | Trans Temperature |   Int Temperature |'
    write(Unit,'(A1,A17, A3,A37, A3,es17.10, A3,es17.10, A2)') '$', adjustr(This%InBinsChar), ' | ', This%QRatioChar(iTTra, iTInt), ' | ', This%TTra(iTTra), ' | ', This%TInt(iTInt), ' |'
    
    write(Unit,'(A)') '#  Process,            Rate,         Rate SD'
    do iP = 1,This%NProc_Cleaned
      write(Unit,'(I10,A1,es16.10,A1,es16.10)') This%ProcessesVecCleaned(iP)%Idx, ',', This%ProcessesVecCleaned(iP)%Temperature(iTTra)%Rate, ',', sqrt( This%ProcessesVecCleaned(iP)%Temperature(iTTra)%RateSD )
    end do
  
  close(Unit)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module