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

Module Process_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use Temperature_Class     ,only:  Temperature_Type

  implicit none

  private
  public    ::    Process_Type

  Type      ::    Process_Type
    character(:)                          ,allocatable ::    Name                      ! Name of Process
    integer                                            ::    Idx        =   -1         ! Idx  of Process
    integer                                            ::    ProcType   =   -1
    integer                                            ::    ExcType    =   -1
    integer                 ,dimension(:) ,allocatable ::    To_Pair     
    integer                 ,dimension(:) ,allocatable ::    To_Species
    integer                 ,dimension(:) ,allocatable ::    To_Level
    character(6)            ,dimension(:) ,allocatable ::    To_LevelChar
    class(Temperature_Type) ,dimension(:) ,allocatable ::    Temperature
    real(rkp) ,dimension(3)                            ::    ArrCoeffs
  contains
    private
    procedure ,public   ::    Initialize       => Initialize_Process
    procedure ,public   ::    Shelving_1stTime => Shelving_1stTime_Process
    procedure ,public   ::    Shelving         => Shelving_Process
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Process( This, NInitMolecules, NTTran, Idx, Name, ProcType, ExcType, To_Pair, To_Level, To_LevelChar, i_Debug )

  class(Process_Type)                       ,intent(out)    ::    This
  integer                                   ,intent(in)     ::    NInitMolecules
  integer                                   ,intent(in)     ::    NTTran
  integer                                   ,intent(in)     ::    Idx
  character(*)                              ,intent(in)     ::    Name
  integer                                   ,intent(in)     ::    ProcType
  integer                                   ,intent(in)     ::    ExcType
  integer ,dimension(:)                     ,intent(in)     ::    To_Pair
  integer ,dimension(:)                     ,intent(in)     ::    To_Level
  character(6) ,dimension(:)                ,intent(in)     ::    To_LevelChar
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  integer                                                   ::    iMol
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Process")
  !i_Debug_Loc   =     Logger%On()

  allocate( This%Name , source = trim(Name) )

  This%Idx      = Idx

  This%ProcType = ProcType
  This%ExcType  = ExcType

  allocate( This%To_Pair(NInitMolecules), Stat=Status  )
  if (Status/=0) call Error( "Error allocating To_Pair in Initialize_Process" )
  
  allocate( This%To_Level(NInitMolecules), Stat=Status  )
  if (Status/=0) call Error( "Error allocating To_Level in Initialize_Process" )
  
  allocate( This%To_LevelChar(NInitMolecules), Stat=Status  )
  if (Status/=0) call Error( "Error allocating To_LevelChar in Initialize_Process" )
  if (i_Debug_Loc) call Logger%Write( "Allocated To_Pair, To_Level, To_LevelChar with Dimension = (",NInitMolecules,")" )

  do iMol = 1,NInitMolecules
    This%To_Pair(iMol)      = To_Pair(iMol)
    This%To_Level(iMol)     = To_Level(iMol)
    This%To_LevelChar(iMol) = To_LevelChar(iMol)
  end do


  allocate( This%Temperature(NTTran), Stat=Status  )
  if (Status/=0) call Error( "Error allocating Temperature in Initialize_Process" )
  if (i_Debug_Loc) call Logger%Write( "Allocated Temperature with Dimension = (",NTTran,")" )


  if (i_Debug_Loc) then
    call Logger%Write( "This%Name         = ", This%Name )
    call Logger%Write( "This%Idx          = ", This%Idx )
    call Logger%Write( "This%To_Pair      = ", This%To_Pair )
    call Logger%Write( "This%To_Level     = ", This%To_Level )
    call Logger%Write( "This%To_LevelChar = ", This%To_LevelChar )
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

!________________________________________________________________________________________________________________________________!
Subroutine Shelving_1stTime_Process( This, NInitMolecules, NTTran, Idx, Name, ProcType, ExcType, To_Pair, To_Level, To_LevelChar, CorrFactor, CrossSect, Velocity, Rate, ArrCoeffs, i_Debug )

  class(Process_Type)                       ,intent(out)    ::    This
  integer                                   ,intent(in)     ::    NInitMolecules
  integer                                   ,intent(in)     ::    NTTran
  integer                                   ,intent(in)     ::    Idx
  character(*)                              ,intent(in)     ::    Name
  integer                                   ,intent(in)     ::    ProcType
  integer                                   ,intent(in)     ::    ExcType
  integer      ,dimension(:)                ,intent(in)     ::    To_Pair
  integer      ,dimension(:)                ,intent(in)     ::    To_Level
  character(6) ,dimension(:)                ,intent(in)     ::    To_LevelChar
  real(rkp)                       ,optional ,intent(in)     ::    CorrFactor
  real(rkp)    ,dimension(2)      ,optional ,intent(in)     ::    CrossSect
  real(rkp)                       ,optional ,intent(in)     ::    Velocity
  real(rkp)    ,dimension(2)      ,optional ,intent(in)     ::    Rate
  real(rkp)    ,dimension(3)      ,optional ,intent(in)     ::    ArrCoeffs
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    CorrFactorTemp = 1.0
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Shelving_1stTime_Process")
  !i_Debug_Loc   =     Logger%On()

  if (present(CorrFactor)) then
    CorrFactorTemp = CorrFactor
  end if

  call This%Initialize( NInitMolecules, NTTran, Idx, Name, ProcType, ExcType, To_Pair, To_Level, To_LevelChar, i_Debug=i_Debug_Loc )

  if (present(CrossSect)) then
    This%Temperature(1)%CrossSect   = CorrFactorTemp    * CrossSect(1)
    This%Temperature(1)%CrossSectSD = CorrFactorTemp**2 * CrossSect(2)**2

    This%Temperature(1)%Rate        = CorrFactorTemp    * CrossSect(1)    * Velocity
    This%Temperature(1)%RateSD      = CorrFactorTemp**2 * CrossSect(2)**2 * Velocity**2
  end if

  if (present(Rate)) then
    This%Temperature(1)%Rate        = CorrFactorTemp    * Rate(1)
    This%Temperature(1)%RateSD      = CorrFactorTemp**2 * Rate(2)**2
  end if

  if (present(ArrCoeffs)) then
    This%ArrCoeffs = ArrCoeffs
  end if
  

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Shelving_Process( This, iTTran, CorrFactor, CrossSect, Velocity, Rate, i_Debug )

  class(Process_Type)                       ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    iTTran
  real(rkp)                       ,optional ,intent(in)     ::    CorrFactor
  real(rkp)    ,dimension(2)      ,optional ,intent(in)     ::    CrossSect
  real(rkp)                       ,optional ,intent(in)     ::    Velocity
  real(rkp)    ,dimension(2)      ,optional ,intent(in)     ::    Rate
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    CorrFactorTemp = 1.0
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Shelving_Process")
  !i_Debug_Loc   =     Logger%On()

  if (present(CorrFactor)) then
    CorrFactorTemp = CorrFactor
  end if

  if (present(CrossSect)) then
    This%Temperature(iTTran)%CrossSect   = This%Temperature(iTTran)%CrossSect   + CorrFactorTemp    * CrossSect(1)
    This%Temperature(iTTran)%CrossSectSD = This%Temperature(iTTran)%CrossSectSD + CorrFactorTemp**2 * CrossSect(2)**2

    This%Temperature(iTTran)%Rate        = This%Temperature(iTTran)%Rate        + CorrFactorTemp    * CrossSect(1)      * Velocity
    This%Temperature(iTTran)%RateSD      = This%Temperature(iTTran)%RateSD      + CorrFactorTemp**2 * CrossSect(2)**2   * Velocity**2
  end if

  if (present(Rate)) then
    This%Temperature(iTTran)%Rate        = This%Temperature(iTTran)%Rate        + CorrFactorTemp    * Rate(1)
    This%Temperature(iTTran)%RateSD      = This%Temperature(iTTran)%Rate        + CorrFactorTemp**2 * Rate(2)**2
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module