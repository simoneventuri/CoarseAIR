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
Module Nb4Atoms_Processes_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Six
  use Processes_Class       ,only:  Processes_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Nb4Atoms_Processes_Type


  Type    ,extends(Processes_Type)  ::    Nb4Atoms_Processes_Type
  contains
    procedure         ::  Initialize                 =>    Initialize_Nb4Atoms
    procedure         ::  Mask4Excahge               =>    Mask4Excahge_Nb4Atoms
    procedure         ::  ConstructVecOfProcs        =>    ConstructVecOfProcs_Nb4Atoms
    procedure         ::  FindingFinalLevel          =>    FindingFinalLevel_Nb4Atoms
    procedure         ::  Convert_CrossSect_To_Rates =>    Convert_CrossSect_To_Rates_Nb4Atoms
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.
  
  contains
  


! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for 3Atoms System
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Initialize_Nb4Atoms( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Nb4Atoms_Processes_Type)                    ,intent(out)    :: This
  Type(Input_Type)                                  ,intent(in)     :: Input
  Type(Collision_Type)                              ,intent(in)     :: Collision
  logical                                 ,optional ,intent(in)     :: i_Debug
  
  character(10)          ,dimension(3)                              :: CollNameA, CollNameB
  integer                                                           :: iP, iOpp
  integer                                                           :: iMol, iMolOpp
  integer                                                           :: iLevel, jLevel, kLevel
  integer                                                           :: Status
  integer                                                           :: NProc
  integer                                                           :: iP1, iP6
  integer                                                           :: NLevelsOpp
  character(10)                                                     :: InProcChar
  logical                                                           :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Nb4Atoms_Processes" )
  !i_Debug_Loc   =     Logger%On()

  This%Initialized  =   .True.


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                                     Allocating Main Quantities                                    !!!
  ! !---------------------------------------------------------------------------------------------------! !  
  allocate(  This%OutputDir, source = adjustl(trim(Input%OutputDir)) ) 

  This%NTraj        = Input%NTraj
  This%System       = adjustl(trim(Input%System))

  This%NPESs        = Input%NPESs
  This%PES_Name     = adjustl(trim(Input%PES_Model(1)))
  This%PESoI        = Input%PESoI
  allocate( This%PESoI_char, source = adjustl(trim(Input%PESoI_char)) )

  iP1 = Collision%Pairs(1)%To_Molecule
  iP6 = Collision%Pairs(1)%To_Molecule
  This%IniMolecules = adjustl(trim(Collision%MoleculesContainer(iP1)%Molecule%Name)) // '-' // adjustl(trim(Collision%MoleculesContainer(iP6)%Molecule%Name))
  allocate(This%InBins(2)); This%InBins = Input%BinOI
  This%InBinsChar   = adjustl(trim( adjustl(trim(Input%BinOI_char(1))) // ',' // adjustl(trim(Input%BinOI_char(2))) ))
  iOpp              = Collision%Pairs(1)%Opposite
  iMol              = Collision%Pairs(iOpp)%To_Molecule
  NLevelsOpp        = Collision%MoleculesContainer(iMolOpp)%Molecule%LevelsContainer%NStates
  This%InProc       = ( This%InBins(1) - 1) * NLevelsOpp + This%InBins(2)
  write(InProcChar, "(I10)") This%InProc
  allocate( This%InProcChar, source = adjustl(trim(InProcChar)) )

  This%NTTra = Input%NTTra
  allocate(This%TTra(This%NTTra));      This%TTra     = int(Input%TtraVec)
  allocate(This%TTraChar(This%NTTra));  This%TTraChar = Input%TtraVecIntChar
  This%NTInt  = Input%NTInt
  allocate(This%TInt(This%NTInt));      This%TInt      = int(Input%TtraVec)
  allocate(This%TIntChar(This%NTInt));  This%TIntChar  = Input%TtraVecIntChar

  allocate(This%QRatio(2,   This%NTTra, This%NTInt)); This%QRatio     = Zero
  allocate(This%QRatioChar( This%NTTra, This%NTInt)); This%QRatioChar = '                                    '
  ! !---------------------------------------------------------------------------------------------------! !


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                      Finding Pairs Corresponding to Equal Exchanges                               !!!
  ! !---------------------------------------------------------------------------------------------------! !
  CollNameA(1) = trim(adjustl(Collision%Pairs(1)%Name)) // trim(adjustl(Collision%Pairs(6)%Name))
  CollNameB(1) = trim(adjustl(Collision%Pairs(6)%Name)) // trim(adjustl(Collision%Pairs(1)%Name))
  CollNameA(2) = trim(adjustl(Collision%Pairs(2)%Name)) // trim(adjustl(Collision%Pairs(2)%Name))
  CollNameB(2) = trim(adjustl(Collision%Pairs(5)%Name)) // trim(adjustl(Collision%Pairs(5)%Name))
  CollNameA(3) = trim(adjustl(Collision%Pairs(3)%Name)) // trim(adjustl(Collision%Pairs(3)%Name))
  CollNameB(3) = trim(adjustl(Collision%Pairs(4)%Name)) // trim(adjustl(Collision%Pairs(4)%Name))

  This%ExcTypeVec = [1,2,3]
  if   ( trim(adjustl(CollNameA(2))) == trim(adjustl(CollNameA(1))) ) then
    This%ExcTypeVec(2) =  1
  elseif ( trim(adjustl(CollNameB(2))) == trim(adjustl(CollNameA(1))) ) then
    This%ExcTypeVec(2) = -1
  end if

  if   ( trim(adjustl(CollNameA(3))) == trim(adjustl(CollNameA(1))) ) then 
    This%ExcTypeVec(3) =  1
  elseif   ( trim(adjustl(CollNameB(3))) == trim(adjustl(CollNameA(1))) ) then
    This%ExcTypeVec(3) = -1
  elseif   ( trim(adjustl(CollNameA(3))) == trim(adjustl(CollNameA(2))) ) then
    This%ExcTypeVec(3) =  This%ExcTypeVec(2)
  elseif   ( trim(adjustl(CollNameB(3))) == trim(adjustl(CollNameA(2))) ) then
    This%ExcTypeVec(3) = -This%ExcTypeVec(2)
  end if
  if (i_Debug_Loc) call Logger%Write( "Constructed This%ExcTypeVec! This%ExcTypeVec = ", This%ExcTypeVec )   
  ! !---------------------------------------------------------------------------------------------------! !


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                         Finding Nb of Processes (Per Pair and Total)                              !!!
  ! !---------------------------------------------------------------------------------------------------! !
  This%NProc_Tot = 0
  do iP = 1,3
    iOpp                   = Collision%Pairs(iP)%Opposite
    iMol                   = Collision%Pairs(iP)%To_Molecule
    This%NProc_iPOpp(iP,1) =    Collision%MoleculesContainer(iMol)%Molecule%LevelsContainer%NStates  + 1
    This%NProc_iPOpp(iP,2) = Collision%MoleculesContainer(iMolOpp)%Molecule%LevelsContainer%NStates  + 1
    This%NProc_iP(iP)      = This%NProc_iPOpp(iP,1) * This%NProc_iPOpp(iP,2)                         + 1
    This%NProc_Tot         = This%NProc_Tot + This%NProc_iP(iP) 
  end do
  allocate( This%Proc_To_LineVec(0:This%NProc_Tot), Stat=Status  )
  if (Status/=0) call Error( "Error allocating Proc_To_LineVec in FindEqExchanges_Nb3_Processes" )
  if (i_Debug_Loc) call Logger%Write( "Allocated Proc_To_LineVec with Dimension = (",This%NProc_Tot,"+1)" )
  This%Proc_To_LineVec = 0
  ! !---------------------------------------------------------------------------------------------------! !


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                                    Creating a Mask for Exchange                                   !!!
  ! !---------------------------------------------------------------------------------------------------! !
  if ( (Input%MergeExchToInelFlg) .or. (Input%MergeExchsFlg) ) then
    call This%Mask4Excahge( Input, i_Debug=i_Debug_Loc )
  end if
  ! !---------------------------------------------------------------------------------------------------! !


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
! !===================================================================================================! !


! !---------------------------------------------------------------------------------------------------! !
!!!                                                                                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Subroutine ConstructVecOfProcs_Nb4Atoms( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Nb4Atoms_Processes_Type)                    ,intent(inout) :: This
  Type(Input_Type)                                  ,intent(in)    :: Input
  Type(Collision_Type)                              ,intent(in)    :: Collision
  logical                                 ,optional ,intent(in)    :: i_Debug
 
  integer                                                          :: iProc, NProc
  integer                                                          :: iP, iOpp
  integer                                                          :: iMol, iMolOpp
  character(6)                                                     :: MolNameA, MolNameB
  integer                                                          :: iA, jA, kA, lA
  integer                                                          :: iLevel, jLevel
  character(6)                                                     :: iLevelChar, jLevelChar
  character(10)          ,dimension(3)                             :: CollNameA, CollNameB
  integer                                                          :: ProcType, ExcType
  character(:)                         ,allocatable                :: Name
  integer                                                          :: Status
  logical                                                          :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ConstructVecOfProcs_Nb4Atoms" )
  !i_Debug_Loc   =     Logger%On()

  This%NProc_Tot = 0
  do iP = 1,3
    iOpp = Collision%Pairs(iP)%Opposite

    iMol                   = Collision%Pairs(iP)%To_Molecule
    iMolOpp                = Collision%Pairs(iOpp)%To_Molecule

    This%NProc_iPOpp(iP,1) =    Collision%MoleculesContainer(iMol)%Molecule%LevelsContainer%NStates + 1
    This%NProc_iPOpp(iP,2) = Collision%MoleculesContainer(iMolOpp)%Molecule%LevelsContainer%NStates + 1

    This%NProc_iP(iP)      = This%NProc_iPOpp(iP,1) * This%NProc_iPOpp(iP,2) + 1
    This%NProc_Tot         = This%NProc_Tot + This%NProc_iP(iP) 
  end do
  allocate(This%ProcessesVec(0:This%NProc_Tot), Stat=Status)
  if (Status/=0) call Error( "Error allocating This%ProcessesVec in InitializeProcesses_Nb4Atoms_Processes" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%ProcessesVec with Dimension = (",This%NProc_Tot,"+1)" )


  iProc = 0
  do iP = 1,3
    iOpp = Collision%Pairs(iP)%Opposite

    iMol     = Collision%Pairs(iP)%To_Molecule
    MolNameA = Collision%MoleculesContainer(iMol)%Molecule%Name

    iMolOpp  = Collision%Pairs(iOpp)%To_Molecule
    MolNameB = Collision%MoleculesContainer(iMolOpp)%Molecule%Name

    iA = Collision%Pair_To_Atoms(iP,1)
    jA = Collision%Pair_To_Atoms(iP,2)
    kA = Collision%Pair_To_Atoms(iOpp,1)
    lA = Collision%Pair_To_Atoms(iOpp,2)

    do iLevel = 0,This%NProc_iPOpp(iP,1)
    
      do jLevel = 0,This%NProc_iPOpp(iP,2)

        if (iLevel == 0) then
          if (jLevel == 0) then
            ProcType = 0
            ExcType  = 0
          else
            ProcType = 0
            ExcType  = 1
          end if      
        elseif (jLevel == 0) then
          ProcType = 0
          ExcType  = 2
        else
          if (iP == 1) then
            ProcType = 1
            ExcType  = 0
          else
            ProcType = 2
            ExcType  = abs(This%ExcTypeVec(iP))
          end if
        end if
        
        call CreateName_Nb4Atoms( MolNameA, MolNameB, iLevel, jLevel, iLevelChar, jLevelChar, Name )
        call This%ProcessesVec(iProc)%Initialize( 2, Input%NTtra, iProc, Name, ProcType, ExcType, [iP, iOpp], [iLevel, jLevel], [iLevelChar, jLevelChar], i_Debug )

      iProc = iProc + 1
      end do

    end do
    
  end do

  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
! !===================================================================================================! !


! !---------------------------------------------------------------------------------------------------! !
!!!                                                                                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Subroutine FindingFinalLevel_Nb4Atoms( This, Input, Collision, vqn, jqn, Arr, Name, ProcType, ExcType, Pairs, iLevelFin, iLevelFinChar, Idx, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type
  use Process_Class               ,only: Process_Type

  class(Nb4Atoms_Processes_Type)                    ,intent(in)  :: This
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
 
  integer                                                        :: iType, jType, Temp
  integer                                                        :: iP, iOpp
  integer                                                        :: iA, jA, kA, lA
  character(6)                                                   :: iLevelChar, jLevelChar
  integer                                                        :: iLevel, jLevel
  integer                                                        :: iMol, iMolOpp
  integer                                                        :: Status
  integer                                                        :: NProcPre, NProcCurrA, NProcCurrB
  character(6)                                                   :: MolNameA, MolNameB
  logical                                                        :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "FindingFinalLevel_Nb4Atoms" )
  !i_Debug_Loc   =     Logger%On()


  iP       = int(Arr / 16.0_rkp)
  iOpp  = Collision%Pairs(iP)%Opposite
  NProcPre = sum( This%NProc_iP(0:iP-1) )
  Temp     = mod(Arr , 16)
  iType    = mod(Temp, 4)
  jType    = int(Temp, 4)

  iMol     = Collision%Pairs(iP)%To_Molecule
  MolNameA = Collision%MoleculesContainer(iMol)%Molecule%Name

  iMolOpp  = Collision%Pairs(iOpp)%To_Molecule
  MolNameB = Collision%MoleculesContainer(iMolOpp)%Molecule%Name

  iA    = Collision%Pair_To_Atoms(iP,1)
  jA    = Collision%Pair_To_Atoms(iP,2)
  kA    = Collision%Pair_To_Atoms(iOpp,1)
  lA    = Collision%Pair_To_Atoms(iOpp,2)

  if (iType > 1) then
    
    iLevel     = 0
    NProcCurrA = 1
    ProcType   = 0
    if (jType > 1) then
      !!! Double Dissociation
      jLevel     = 0
      NProcCurrB = 1
      ExcType    = 0
    
    else
      !!! 1st Pair Dissociated, 2nd Bound
      jLevel     = Collision%MoleculesContainer(iMolOpp)%Molecule%BinsContainer%qns_to_Bin(vqn(2),jqn(2))
      NProcCurrB = jLevel + 1
      ExcType    = 1
    end if
  
  elseif (jType > 1) then
    !!! 1st Pair Bound, 2nd Dissociated
    iLevel     = Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%qns_to_Bin(vqn(1),jqn(1))
    NProcCurrA = iLevel + 1
    jLevel     = 0
    NProcCurrB = 1
    ProcType   = 0
    ExcType    = 2
  
  else
    !!! Either Inelastic or Exchange
    iLevel     = Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%qns_to_Bin(vqn(1),jqn(1))
    NProcCurrA = iLevel + 1
    jLevel     = Collision%MoleculesContainer(iMolOpp)%Molecule%BinsContainer%qns_to_Bin(vqn(2),jqn(2))
    if (iP == 1) then
      ProcType = 1
      ExcType  = 0
    else
      ProcType = 2
      ExcType  = abs(This%ExcTypeVec(iP))
    end if

  end if

  call CreateName_Nb4Atoms( MolNameA, MolNameB, iLevel, jLevel, iLevelChar, jLevelChar, Name )  
  iLevelFin     = [iLevel,     jLevel]
  iLevelFinChar = [iLevelChar, jLevelChar]
  Idx           = NProcPre + ( (NProcCurrA-1)*This%NProc_iPOpp(iP,2) + NProcCurrB + 1)
  Pairs         = [iP, iOpp]

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! !---------------------------------------------------------------------------------------------------! !
!!!                                    Creating a Mask for Exchange                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Subroutine Mask4Excahge_Nb4Atoms( This, Input, i_Debug )

  use Input_Class                 ,only: Input_Type

  class(Nb4Atoms_Processes_Type)                    ,intent(inout)  :: This
  Type(Input_Type)                                  ,intent(in)     :: Input
  logical                                 ,optional ,intent(in)     :: i_Debug
  
  integer                                                           :: iP
  integer                                                           :: Status
  integer                                                           :: iLevel, jLevel, kLevel
  logical                                                           :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Mask4Excahge_Nb4Atoms" )
  !i_Debug_Loc   =     Logger%On() 

  if (Input%MergeExchToInelFlg) then
    This%MergeExchToInelFlg = .True.

    This%NProc_Tot_Unique = This%NProc_iP(1)
    if (This%ExcTypeVec(2) /= 1) then
      This%NProc_Tot_Unique = This%NProc_Tot_Unique + This%NProc_iP(2)
    end if
    if ( (This%ExcTypeVec(3) /= 1) .and. (This%ExcTypeVec(3) /= 2) ) then
      This%NProc_Tot_Unique = This%NProc_Tot_Unique + This%NProc_iP(3)
    end if
    allocate( This%ExchMask(0:This%NProc_Tot_Unique), Stat=Status  )
    if (Status/=0) call Error( "Error allocating ExchMask in Initialize_Nb4Atoms" )
    if (i_Debug_Loc) call Logger%Write( "Allocated ExchMask with Dimension = (",This%NProc_Tot_Unique,"+1)" )


    iP     = 1
    do iLevel = 1,This%NProc_iPOpp(iP,1)
      do jLevel = 1,This%NProc_iPOpp(iP,2)
        kLevel                                       = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
        This%ExchMask( kLevel )                      = kLevel
      end do
    end do

    iP     = 2
    if (This%ExcTypeVec(iP) == 1) then
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                     = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(1) ) = kLevel
        end do
      end do
    elseif (This%ExcTypeVec(iP) == -1) then
      do jLevel = 1,This%NProc_iPOpp(iP,2)
        do iLevel = 1,This%NProc_iPOpp(iP,1)
          kLevel                                     = (jLevel - 1)*This%NProc_iPOpp(iP,1) + iLevel
          This%ExchMask( kLevel + This%NProc_iP(1) ) = kLevel
        end do
      end do
    else
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(1) )                    = kLevel + This%NProc_iP(1)
        end do
      end do
    end if

    iP     = 3
    if (This%ExcTypeVec(iP) == 1) then
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel
        end do
      end do
    elseif (This%ExcTypeVec(iP) == -1) then
      do jLevel = 1,This%NProc_iPOpp(iP,2)
        do iLevel = 1,This%NProc_iPOpp(iP,1)
          kLevel                                                        = (jLevel - 1)*This%NProc_iPOpp(iP,1) + iLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel
        end do
      end do
    elseif (This%ExcTypeVec(iP) == 2) then
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel + This%NProc_iP(1)
        end do
      end do
    elseif (This%ExcTypeVec(iP) == -2) then
      do jLevel = 1,This%NProc_iPOpp(iP,2)
        do iLevel = 1,This%NProc_iPOpp(iP,1)
          kLevel                                                        = (jLevel - 1)*This%NProc_iPOpp(iP,1) + iLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel + This%NProc_iP(1)
        end do
      end do
    else
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel + This%NProc_iP(2) + This%NProc_iP(1)
        end do
      end do
    end if


  elseif (Input%MergeExchsFlg) then
    This%MergeExchsFlg = .True.

    This%NProc_Tot_Unique = This%NProc_iP(1)
    This%NProc_Tot_Unique = This%NProc_Tot_Unique   + This%NProc_iP(2)
    if (This%ExcTypeVec(3) /= 2) then
      This%NProc_Tot_Unique = This%NProc_Tot_Unique + This%NProc_iP(3)
    end if
    allocate( This%ExchMask(0:This%NProc_Tot_Unique), Stat=Status  )
    if (Status/=0) call Error( "Error allocating ExchMask in Initialize_Nb4Atoms" )
    if (i_Debug_Loc) call Logger%Write( "Allocated ExchMask with Dimension = (",This%NProc_Tot_Unique,"+1)" )

    iP = 1
    do iLevel = 1,This%NProc_iPOpp(iP,1)
      do jLevel = 1,This%NProc_iPOpp(iP,2)
        kLevel                                                          = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
        This%ExchMask( kLevel )                                         = kLevel
      end do
    end do

    iP = 2
    do iLevel = 1,This%NProc_iPOpp(iP,1)
      do jLevel = 1,This%NProc_iPOpp(iP,1)
        kLevel                                                          = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
        This%ExchMask( kLevel + This%NProc_iP(1) )                      = kLevel + This%NProc_iP(1)
      end do
    end do

    iP = 3
    if (This%ExcTypeVec(iP) == 2) then
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel + This%NProc_iP(1)
        end do
      end do
    elseif (This%ExcTypeVec(iP) == -2) then
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (jLevel - 1)*This%NProc_iPOpp(iP,1) + iLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel + This%NProc_iP(1)
        end do
      end do
    else
      do iLevel = 1,This%NProc_iPOpp(iP,1)
        do jLevel = 1,This%NProc_iPOpp(iP,2)
          kLevel                                                        = (iLevel - 1)*This%NProc_iPOpp(iP,2) + jLevel
          This%ExchMask( kLevel + This%NProc_iP(2) + This%NProc_iP(1) ) = kLevel + This%NProc_iP(1) + This%NProc_iP(1)
        end do
      end do
    end if

  end if
  ! !---------------------------------------------------------------------------------------------------! !

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   READING INITIAL AND FINAL CONDITIONS AND CROSS SECTIONS FROM statistics.out
! ==============================================================================================================
Subroutine Convert_CrossSect_To_Rates_Nb4Atoms( This, Input, Collision, Velocity, i_Debug, i_Debug_Deep)

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero
  
  class(Nb4Atoms_Processes_Type)                        ,intent(inout) ::    This
  type(Input_Type)                                      ,intent(in)    ::    Input
  type(Collision_Type)                                  ,intent(in)    ::    Collision
  real(rkp)                  ,dimension(:)              ,intent(in)    ::    Velocity
  logical                                     ,optional ,intent(in)    ::    i_Debug
  logical                                     ,optional ,intent(in)    ::    i_Debug_Deep

  character(:)                    ,allocatable              ::    FileName
  integer                                                   ::    Unit
  integer                                                   ::    Status
  integer                    ,dimension(2)                  ::    vqnIn, jqnIn, vqnFin, jqnFin
  integer                                                   ::    ArrIn, ArrFin
  real(rkp)    ,dimension(2)                                ::    CrossSectTemp
  integer                                                   ::    iProc, Proc_To_Line
  integer                                                   ::    iLine, NLine
  integer                                                   ::    iTtra, NTtra
  integer                                                   ::    iTint, NTint
  character(:)                    ,allocatable              ::    Name 
  integer                                                   ::    ProcType
  integer                                                   ::    ExcType
  integer       ,dimension(2)                               ::    iP
  integer       ,dimension(2)                               ::    iLevelFin
  character(6)  ,dimension(2)                               ::    iLevelFinChar  
  integer                                                   ::    Idx
  logical                                                   ::    i_Debug_Loc
  logical                                                   ::    i_Debug_Loc_Deep


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Convert_CrossSect_To_Rates_Nb3Atoms" )
  i_Debug_Loc_Deep = i_Debug_Global; if ( present(i_Debug_Deep) )i_Debug_Loc_Deep = i_Debug_Deep
  !i_Debug_Loc   =     Logger%On()
  

  if (i_Debug_Loc) call Logger%Write( "Initializing the Processes. Calling Processes%Initialize" )
  call This%Initialize( Input, Collision, i_Debug=i_Debug_Loc )

  !do iTtra = 1,Input%NTtra
  !  if (i_Debug_Loc) call Logger%Write( "  iTtra = ", iTtra )
    NTtra = 1 !Input%NTtra
    iTtra = 1

    !do iTint = 1,Input%NTint
    !  if (i_Debug_Loc) call Logger%Write( "  iTint = ", iTint )
      NTint = 1
      iTInt = 1

      FileName = trim(adjustl(Input%LevelOutputDir))// '/statistics.out'
      if (i_Debug_Loc) call Logger%Write( "  Reading File: ", FileName )
      Open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
        if (Status/=0) then
          if (i_Debug_Loc) call Logger%Write( "The statistics.out File is not present for this Initial Level / Bin." )
        else
          read(Unit,*,iostat=Status)
          NLine = 0
          do
            read(Unit,*,iostat=Status)
            if (Status .ne. 0) exit
            NLine = NLine + 1
          end do
          if (i_Debug_Loc) call Logger%Write( "  Found ", NLine, " Processes in the Statistics File." )

          !!! Creating an Array of Temporary Processes
          allocate( This%ProcessesVecTemp(NLine), Stat=Status  )
          if (Status/=0) call Error( "Error allocating ProcessesVecTemp in Convert_CrossSect_In_Rates_Nb3Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Allocated ProcessesVecTemp with Dimension = (",NLine,")" )

      Rewind(Unit)

          read(Unit,*,iostat=Status)
          This%NProc_Cleaned = 0
          do iLine=1,NLine
            if (i_Debug_Loc_Deep) call Logger%Write( "    iLine = ", iLine )

            !!! Readin Final State !!!!
            read(Unit,*,iostat=Status) jqnFin(1), vqnFin(1), jqnFin(2), vqnFin(2), ArrFin, jqnIn(1), vqnIn(1), jqnIn(2), vqnIn(2), ArrIn, CrossSectTemp(1), CrossSectTemp(2)
            if (i_Debug_Loc_Deep) call Logger%Write( "    vqnIn  = ", vqnIn(1:2)  )
            if (i_Debug_Loc_Deep) call Logger%Write( "    jqnIn  = ", jqnIn(1:2)  )
            if (i_Debug_Loc_Deep) call Logger%Write( "    ArrIn  = ", ArrIn       )
            if (i_Debug_Loc_Deep) call Logger%Write( "    vqnFin = ", vqnFin(1:2) )
            if (i_Debug_Loc_Deep) call Logger%Write( "    jqnFin = ", jqnFin(1:2) )
            if (i_Debug_Loc_Deep) call Logger%Write( "    ArrFin = ", ArrFin      )
            if (i_Debug_Loc_Deep) call Logger%Write( "    CrossSect and CrossSectSD = ", CrossSectTemp(1:2) )   


            !!! Postprocessing Final State for Reconstructing Process !!!!
            call This%FindingFinalLevel( Input, Collision, vqnFin, jqnFin, ArrFin, Name, ProcType, ExcType, iP, iLevelFin, iLevelFinChar, Idx, i_Debug=i_Debug_Loc )

            !!! Filtering Exchanges !!!!
            if ( (This%MergeExchToInelFlg) .or. (This%MergeExchsFlg) ) then
              Idx = This%ExchMask(Idx)
            end if

            !!! New Process ??? !!!!
            Proc_To_Line = This%Proc_To_LineVec(Idx)
            if (Proc_To_Line < 1) then
              !!! New Process! Allocating it !!!!
              call This%ProcessesVecTemp(iLine)%Shelving_1stTime( 2, NTtra, Idx, Name, ProcType, ExcType, iP, iLevelFin, iLevelFinChar, CorrFactor=1.0, CrossSect=CrossSectTemp, Velocity=Velocity(iTtra), i_Debug=i_Debug_Loc )
              This%Proc_To_LineVec(Idx) = iLine
              This%NProc_Cleaned        = This%NProc_Cleaned + 1
            else
              !!! Old Process! Adding Cross Section !!!!
              call This%ProcessesVecTemp(Proc_To_Line)%Shelving( iTtra, CorrFactor=1.0, CrossSect=CrossSectTemp, Velocity=Velocity(iTtra), i_Debug=i_Debug_Loc )
            end if

          end do


          !!! Creating an Array of Unique Processes and Deallocating Temporary One
          allocate( This%ProcessesVecCleaned(This%NProc_Cleaned), Stat=Status  )
          if (Status/=0) call Error( "Error allocating ProcessesVecCleaned in Convert_CrossSect_In_Rates_Nb3Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Allocated ProcessesVecCleaned with Dimension = (",This%NProc_Cleaned,")" )
          
          do iProc=1,This%NProc_Cleaned
            This%ProcessesVecCleaned(iProc) = This%ProcessesVecTemp(iProc)
          end do
          
          deallocate(This%ProcessesVecTemp, Stat=Status)
          if (Status/=0) call Error( "Error deallocating This%ProcessesVecTemp in Convert_CrossSect_In_Rates_Nb3Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Deallocated This%ProcessesVecTemp" )

        end if
      Close(Unit)
      if (i_Debug_Loc) call Logger%Write( "    Done Reading File: ", FileName )

      !!! Writing Rates !!!!
      call This%InProc_WritingRates( iTTra, iTInt, i_Debug=i_Debug_Loc )

    !end do

  !end do


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! !---------------------------------------------------------------------------------------------------! !
!!!                                                                                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Pure Subroutine CreateName_Nb4Atoms( MolNameA, MolNameB, iLevel, jLevel, iLevelChar, jLevelChar, Name )

  character(*)               ,intent(in)  :: MolNameA
  character(*)               ,intent(in)  :: MolNameB
  integer                    ,intent(in)  :: iLevel
  integer                    ,intent(in)  :: jLevel
  character(6)               ,intent(out) :: iLevelChar
  character(6)               ,intent(out) :: jLevelChar
  character(20)              ,intent(out) :: Name

  if (iLevel<10)then
    write(iLevelChar,'(I1)') iLevel
  elseif (iLevel<100) then
    write(iLevelChar,'(I2)') iLevel
  elseif (iLevel<1000) then
    write(iLevelChar,'(I3)') iLevel
  elseif (iLevel<10000) then
    write(iLevelChar,'(I4)') iLevel
  elseif (iLevel<100000) then
    write(iLevelChar,'(I5)') iLevel
  end if

  if (jLevel<10)then
    write(jLevelChar,'(I1)') jLevel
  elseif (jLevel<100) then
    write(jLevelChar,'(I2)') jLevel
  elseif (jLevel<1000) then
    write(jLevelChar,'(I3)') jLevel
  elseif (jLevel<10000) then
    write(jLevelChar,'(I4)') jLevel
  elseif (jLevel<100000) then
    write(jLevelChar,'(I5)') jLevel
  end if

  Name = trim(adjustl(MolNameA)) // "(" // trim(adjustl(iLevelChar)) // ")+" // trim(adjustl(MolNameB)) // "(" // trim(adjustl(jLevelChar)) // ")"


end Subroutine
! !===================================================================================================! !


End Module