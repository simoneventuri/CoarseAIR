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
  Type(Collision_Type)                              ,intent(inout)  :: Collision
  logical                                 ,optional ,intent(in)     :: i_Debug
  
  integer                                                           :: iP, iOpp
  integer                                                           :: iMol, iMolOpp
  integer                                                           :: iLevel, jLevel, kLevel, iTemp
  integer                                                           :: Status
  integer                                                           :: NProc
  integer                                                           :: NLevelsOpp
  character(10)                                                     :: InProcChar
  character(20)                                                     :: Temp1Char, Temp2Char
  logical                                                           :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Nb4Atoms_Processes" )
  !i_Debug_Loc   =     Logger%On()

  This%Initialized  =   .True.


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                                     Allocating Main Quantities                                    !!!
  ! !---------------------------------------------------------------------------------------------------! !  
  allocate(  This%OutputDir,      source = adjustl(trim(Input%OutputDir))      ) 
  allocate(  This%LevelOutputDir, source = adjustl(trim(Input%LevelOutputDir)) ) 
  This%FirstIssueFlg = .True.

  This%NTraj        = Input%NTraj
  This%System       = adjustl(trim(Input%System))

  This%NPESs        = Input%NPESs
  This%PESoI        = Input%PESoI
  This%StochPESFlg  = Input%StochPESFlg
  iTemp             = 1
  if ( (.not. This%StochPESFlg) .and. (This%PESoI>0) ) iTemp = This%PESoI
  This%PES_Name     = adjustl(trim(Input%PES_Model(iTemp)))
  allocate( This%PESoI_char, source = adjustl(trim(Input%PESoI_char)) )

  This%MergeExchsFlg      = Input%MergeExchsFlg
  This%MergeExchToInelFlg = Input%MergeExchToInelFlg

  iMol    = Collision%Pairs(1)%To_Molecule
  if (i_Debug_Loc) call Logger%Write( "Molecule associated to Pair 1 is the Molecule Nb ",iMol )
  iOpp    = Collision%Pairs(1)%Opposite
  iMolOpp = Collision%Pairs(iOpp)%To_Molecule
  if (i_Debug_Loc) call Logger%Write( "Molecule associated to Pair 6 is the Molecule Nb ",iMolOpp )


  This%IniMolecules = adjustl(trim(Collision%MoleculesContainer(iMol)%Molecule%Name)) // '-' // adjustl(trim(Collision%MoleculesContainer(iMolOpp)%Molecule%Name))
  if (i_Debug_Loc) call Logger%Write( "This%IniMolecules = ", This%IniMolecules )
  allocate(This%InBins(2)); This%InBins = Input%BinOI
  if (i_Debug_Loc) call Logger%Write( "This%InBins = ", This%InBins )
  This%InBinsChar     = adjustl(trim( adjustl(trim(Input%BinOI_char(1))) // ',' // adjustl(trim(Input%BinOI_char(2))) ))
  This%InBinsCharName = adjustl(trim( 'i' // adjustl(trim(Input%BinOI_char(1))) // '_j' // adjustl(trim(Input%BinOI_char(2))) ))

  NLevelsOpp        = Collision%MoleculesContainer(iMolOpp)%Molecule%BinsContainer%NBins
  This%InProc       = ( This%InBins(1) - 1) * NLevelsOpp + This%InBins(2)
  if (i_Debug_Loc) call Logger%Write( "This%InProc = ", This%InProc )
  write(InProcChar, "(I10)") This%InProc
  allocate( This%InProcChar, source = adjustl(trim(InProcChar)) )
  if (i_Debug_Loc) call Logger%Write( "This%InProcChar = ", This%InProcChar )

  This%NTTra = Input%NTTra
  allocate(This%TTra(This%NTTra));      This%TTra     = int(Input%TtraVec)
  allocate(This%TTraChar(This%NTTra));  This%TTraChar = Input%TtraVecIntChar
  This%NTInt  = Input%NTInt
  allocate(This%TInt(This%NTInt));      This%TInt      = int(Input%TtraVec)
  allocate(This%TIntChar(This%NTInt));  This%TIntChar  = Input%TtraVecIntChar

  !allocate(This%QRatio(2,   This%NTTra, This%NTInt)); This%QRatio     = Zero
  !allocate(This%QRatioChar( This%NTTra, This%NTInt)); This%QRatioChar = '                                    '
  allocate(This%QRatio(2)); This%QRatio     = Zero
  This%QRatio(1) = Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%Bin(This%InBins(1))%QRatio
  This%QRatio(2) = Collision%MoleculesContainer(iMolOpp)%Molecule%BinsContainer%Bin(This%InBins(2))%QRatio
  write(Temp1Char, '(es20.10)') This%QRatio(1)
  write(Temp2Char, '(es20.10)') This%QRatio(2)
  This%QRatioChar = adjustl(trim( adjustr(trim(Temp1Char)) // ',' // adjustl(trim(Temp2Char)) ))
  ! !---------------------------------------------------------------------------------------------------! !


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                         Finding Nb of Processes (Per Pair and Total)                              !!!
  ! !---------------------------------------------------------------------------------------------------! !
  This%NProc_Tot = 0
  do iP = 1,3
    iOpp                   = Collision%Pairs(iP)%Opposite
    iMol                   = Collision%Pairs(iP)%To_Molecule
    This%NProc_iPOpp(iP,1) =    Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%NBins  + 1
    This%NProc_iPOpp(iP,2) = Collision%MoleculesContainer(iMolOpp)%Molecule%BinsContainer%NBins  + 1
    This%NProc_iP(iP)      = This%NProc_iPOpp(iP,1) * This%NProc_iPOpp(iP,2)                  
    This%NProc_Tot         = This%NProc_Tot + This%NProc_iP(iP)                                  
  end do
  This%NProc_Tot = This%NProc_Tot + 1
  allocate( This%Proc_To_LineVec(0:This%NProc_Tot-1), Stat=Status  )
  if (Status/=0) call Error( "Error allocating Proc_To_LineVec in FindEqExchanges_Nb3_Processes" )
  if (i_Debug_Loc) call Logger%Write( "Allocated Proc_To_LineVec with Dimension = (",This%NProc_Tot,")" )
  This%Proc_To_LineVec = 0
  ! !---------------------------------------------------------------------------------------------------! !


  ! !---------------------------------------------------------------------------------------------------! !
  !!!                                    Creating a Mask for Exchange                                   !!!
  ! !---------------------------------------------------------------------------------------------------! !
  call This%Mask4Excahge( Input, Collision, i_Debug=i_Debug_Loc )
  ! !---------------------------------------------------------------------------------------------------! !


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
! !===================================================================================================! !


! !---------------------------------------------------------------------------------------------------! !
!!!                                    Creating a Mask for Exchange                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Subroutine Mask4Excahge_Nb4Atoms( This, Input, Collision, i_Debug )

  use Input_Class                 ,only: Input_Type
  use Collision_Class             ,only: Collision_Type

  class(Nb4Atoms_Processes_Type)                    ,intent(inout)  :: This
  Type(Input_Type)                                  ,intent(in)     :: Input
  Type(Collision_Type)                              ,intent(inout)  :: Collision
  logical                                 ,optional ,intent(in)     :: i_Debug
  
  character(10)          ,dimension(3)                              :: CollNameA, CollNameB
  integer                                                           :: iP, jP, kP, gP, iPOpp, jPOpp, kPOpp
  integer                                                           :: Status, NTotTemp
  integer                                                           :: iMol1, jMol1, iMol, jMol, kMol, lMol
  integer                                                           :: ProcType
  logical                                                           :: i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Mask4Excahge_Nb4Atoms" )
  !i_Debug_Loc   =     Logger%On() 

  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Creating:
  !!!   - This%ExcTypeVec
  !!!
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Creating:
  !!!   - Collision%Pairs(...)%NProc
  !!!   - Collision%Pairs(...)%To_Pair_Exch
  !!!   - This%NProcType
  !!!
  This%NProcType = 1
  iMol1 = Collision%Pairs(1)%To_Molecule
  Collision%Pairs(1)%NProc = Collision%MoleculesContainer(iMol1)%Molecule%BinsContainer%NBins + 1
  jMol1 = Collision%Pairs(6)%To_Molecule
  Collision%Pairs(6)%NProc = Collision%MoleculesContainer(jMol1)%Molecule%BinsContainer%NBins + 1

  iP=2
  iPOpp = Collision%Pairs(iP)%Opposite
  iMol  = Collision%Pairs(iP)%To_Molecule
  Collision%Pairs(iP)%NProc    = Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%NBins + 1
  jMol  = Collision%Pairs(iPOpp)%To_Molecule
  Collision%Pairs(iPOpp)%NProc = Collision%MoleculesContainer(jMol)%Molecule%BinsContainer%NBins + 1
  if ( (iMol == iMol1) .and. (jMol == jMol1) .and. (This%MergeExchToInelFlg) ) then
    Collision%Pairs(iP)%To_Pair_Exch =  1
  elseif ( (iMol == jMol1) .and. (jMol == iMol1) .and. (This%MergeExchToInelFlg) ) then
    Collision%Pairs(iP)%To_Pair_Exch = -1
  else
    Collision%Pairs(iP)%To_Pair_Exch = iP
    This%NProcType = This%NProcType + 1
  end if

  iP=3
  iPOpp = Collision%Pairs(iP)%Opposite
  iMol  = Collision%Pairs(iP)%To_Molecule
  Collision%Pairs(iP)%NProc    = Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%NBins + 1
  jMol  = Collision%Pairs(iPOpp)%To_Molecule
  Collision%Pairs(iPOpp)%NProc = Collision%MoleculesContainer(jMol)%Molecule%BinsContainer%NBins + 1
  if ( (iMol == iMol1) .and. (jMol == jMol1) .and. (This%MergeExchToInelFlg) ) then
    Collision%Pairs(iP)%To_Pair_Exch =  1
  elseif ( (iMol == jMol1) .and. (jMol == iMol1) .and. (This%MergeExchToInelFlg) ) then
    Collision%Pairs(iP)%To_Pair_Exch = -1
  else
    kP=2
    kMol  = Collision%Pairs(kP)%To_Molecule
    kPOpp = Collision%Pairs(kP)%Opposite
    lMol  = Collision%Pairs(kPOpp)%To_Molecule
    if ( (iMol == kMol) .and. (jMol == lMol) .and. (This%MergeExchsFlg) ) then
      Collision%Pairs(iP)%To_Pair_Exch =   Collision%Pairs(kP)%To_Pair_Exch
    elseif ( (iMol == lMol) .and. (jMol == kMol) .and. (This%MergeExchsFlg) ) then
      Collision%Pairs(iP)%To_Pair_Exch = - Collision%Pairs(kP)%To_Pair_Exch
    else
      Collision%Pairs(iP)%To_Pair_Exch = iP
      This%NProcType = This%NProcType + 1
    end if
  end if

  if (i_Debug_Loc) call Logger%Write( "Constructed Collision%Pairs(...)%To_Pair_Exch!" )   
  if (i_Debug_Loc) call Logger%Write( "                                               Collision%Pairs(1)%To_Pair_Exch = ", Collision%Pairs(1)%To_Pair_Exch )   
  if (i_Debug_Loc) call Logger%Write( "                                               Collision%Pairs(2)%To_Pair_Exch = ", Collision%Pairs(2)%To_Pair_Exch )   
  if (i_Debug_Loc) call Logger%Write( "                                               Collision%Pairs(3)%To_Pair_Exch = ", Collision%Pairs(3)%To_Pair_Exch )   
  if (i_Debug_Loc) call Logger%Write( "This%NProcType = ", This%NProcType )   



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Creating:
  !!!   - This%ProcTypeVec: 
  !!! 
  allocate(This%ProcTypeVec(This%NProcType))
  allocate(This%OvProcRates(This%NProcType,2))

  ProcType            = 1
  This%ProcTypeVec(1) = ProcType
  iMol1 = Collision%Pairs(1)%To_Molecule
  jMol1 = Collision%Pairs(6)%To_Molecule

  iP=2
  iPOpp = Collision%Pairs(iP)%Opposite
  iMol  = Collision%Pairs(iP)%To_Molecule
  jMol  = Collision%Pairs(iPOpp)%To_Molecule
  if (.not. ( (iMol == iMol1) .and. (jMol == jMol1) .and. (This%MergeExchToInelFlg) ) ) then
    kMol  = Collision%Pairs(kP)%To_Molecule
    kPOpp = Collision%Pairs(kP)%Opposite
    lMol  = Collision%Pairs(kPOpp)%To_Molecule
    ProcType                   = ProcType + 1
    This%ProcTypeVec(ProcType) = ProcType 
  end if

  iP=3
  iPOpp = Collision%Pairs(iP)%Opposite
  iMol  = Collision%Pairs(iP)%To_Molecule
  jMol  = Collision%Pairs(iPOpp)%To_Molecule
  if ( (.not. ( (iMol == iMol1) .and. (jMol == jMol1) .and. (This%MergeExchToInelFlg) ) ) .and. &
       (.not. ( (jMol == iMol1) .and. (iMol == jMol1) .and. (This%MergeExchToInelFlg) ) ) ) then
    kP=2
    kMol  = Collision%Pairs(kP)%To_Molecule
    kPOpp = Collision%Pairs(kP)%Opposite
    lMol  = Collision%Pairs(kPOpp)%To_Molecule
    if ( (.not. ( (iMol == kMol) .and. (jMol == lMol) .and. (This%MergeExchsFlg) ) ) .and. &
         (.not. ( (jMol == kMol) .and. (iMol == lMol) .and. (This%MergeExchsFlg) ) ) ) then
      ProcType                   = ProcType + 1
      This%ProcTypeVec(ProcType) = ProcType 
    end if
  end if
  if (i_Debug_Loc) call Logger%Write( "Constructed This%ProcTypeVec!                  This%ProcTypeVec = ", This%ProcTypeVec )   



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Creating:
  !!!   - Collision%Pairs(...)%NPrevProc
  !!!
  iP=1
  Collision%Pairs(iP)%NPrevProc     = 1 
  NTotTemp                          = Collision%Pairs(1)%NProc * Collision%Pairs(6)%NProc + 1
  do iP=2,3
    iPOpp = Collision%Pairs(iP)%Opposite
    if (Collision%Pairs(iP)%To_Pair_Exch == Collision%Pairs(iP-1)%To_Pair_Exch) then
      jP                            = Collision%Pairs(iP)%To_Pair_Exch
      Collision%Pairs(iP)%NPrevProc = Collision%Pairs(jP)%NPrevProc
    else
      Collision%Pairs(iP)%NPrevProc = NTotTemp
      NTotTemp                      = NTotTemp + Collision%Pairs(iP)%NProc * Collision%Pairs(iPOpp)%NProc
    end if
  end do
  if (i_Debug_Loc) call Logger%Write( "Constructed Collision%Pairs(...)%To_Pair_Exch and Collision%Pairs(...)%NProc!" )   
  if (i_Debug_Loc) call Logger%Write( "                                               Collision%Pairs(1)%NPrevProc = ", Collision%Pairs(1)%NPrevProc, "; Collision%Pairs(1)%NProc = ", Collision%Pairs(1)%NProc )   
  if (i_Debug_Loc) call Logger%Write( "                                               Collision%Pairs(2)%NPrevProc = ", Collision%Pairs(2)%NPrevProc, "; Collision%Pairs(2)%NProc = ", Collision%Pairs(2)%NProc )   
  if (i_Debug_Loc) call Logger%Write( "                                               Collision%Pairs(3)%NPrevProc = ", Collision%Pairs(3)%NPrevProc, "; Collision%Pairs(3)%NProc = ", Collision%Pairs(3)%NProc )   


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


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

    This%NProc_iPOpp(iP,1) =    Collision%MoleculesContainer(iMol)%Molecule%BinsContainer%NBins + 1
    This%NProc_iPOpp(iP,2) = Collision%MoleculesContainer(iMolOpp)%Molecule%BinsContainer%NBins + 1

    This%NProc_iP(iP)      = This%NProc_iPOpp(iP,1) * This%NProc_iPOpp(iP,2)
    This%NProc_Tot         = This%NProc_Tot + This%NProc_iP(iP) 
  end do
  This%NProc_Tot = This%NProc_Tot + 1
  allocate(This%ProcessesVec(0:This%NProc_Tot-1), Stat=Status)
  if (Status/=0) call Error( "Error allocating This%ProcessesVec in InitializeProcesses_Nb4Atoms_Processes" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%ProcessesVec with Dimension = (",This%NProc_Tot,")" )


  iProc = 0
  do iP = 1,3
    iOpp = Collision%Pairs(iP)%Opposite

    iMol     = Collision%Pairs(iP)%To_Molecule
    MolNameA = Collision%MoleculesContainer(iMol)%Molecule%Name

    iMolOpp  = Collision%Pairs(iOpp)%To_Molecule
    MolNameB = Collision%MoleculesContainer(iMolOpp)%Molecule%Name

    iA = Collision%Pair_To_Atoms(1,iP)
    jA = Collision%Pair_To_Atoms(2,iP)
    kA = Collision%Pair_To_Atoms(1,iOpp)
    lA = Collision%Pair_To_Atoms(2,iOpp)

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
  character(:)                  ,allocatable        ,intent(out) :: Name 
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

  if (i_Debug_Loc) call Logger%Write( "vqn = ", vqn )
  if (i_Debug_Loc) call Logger%Write( "jqn = ", jqn )
  if (i_Debug_Loc) call Logger%Write( "Arr = ", Arr )

  iP       = int(Arr / 16.0_rkp)
  if (i_Debug_Loc) call Logger%Write( "iP = ", iP )
  if (iP > 0) then
    iOpp     = Collision%Pairs(iP)%Opposite
    NProcPre = Collision%Pairs(iP)%NPrevProc
    if (i_Debug_Loc) call Logger%Write( "Collision%Pairs(iP)%NPrevProc = ", Collision%Pairs(iP)%NPrevProc )
    Temp     = mod(Arr , 16)
    iType    = mod(Temp, 4)
    jType    = int(Temp / 4)
    if (i_Debug_Loc) call Logger%Write( "iType = ", iType )
    if (i_Debug_Loc) call Logger%Write( "jType = ", jType )
    if (i_Debug_Loc) call Logger%Write( "NProcPre = ", NProcPre )

    iMol     = Collision%Pairs(iP)%To_Molecule
    MolNameA = Collision%MoleculesContainer(iMol)%Molecule%Name
    if (i_Debug_Loc) call Logger%Write( "iMol = ", iMol )

    iMolOpp  = Collision%Pairs(iOpp)%To_Molecule
    MolNameB = Collision%MoleculesContainer(iMolOpp)%Molecule%Name
    if (i_Debug_Loc) call Logger%Write( "iMolOpp = ", iMolOpp )

    iA    = Collision%Pair_To_Atoms(1,iP)
    jA    = Collision%Pair_To_Atoms(2,iP)
    kA    = Collision%Pair_To_Atoms(1,iOpp)
    lA    = Collision%Pair_To_Atoms(2,iOpp)
    if (i_Debug_Loc) call Logger%Write( "iA = ", iA )
    if (i_Debug_Loc) call Logger%Write( "jA = ", jA )
    if (i_Debug_Loc) call Logger%Write( "kA = ", kA )
    if (i_Debug_Loc) call Logger%Write( "lA = ", lA )

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
      NProcCurrB = jLevel + 1
      if (iP == 1) then
        ProcType = 1
        ExcType  = 0
      else
        ProcType = 2
        ExcType  = abs(This%ExcTypeVec(iP))
      end if

    end if

    if (i_Debug_Loc) call Logger%Write( "ProcType   = ", ProcType )
    if (i_Debug_Loc) call Logger%Write( "iLevel     = ", iLevel )
    if (i_Debug_Loc) call Logger%Write( "jLevel     = ", jLevel )
    if (i_Debug_Loc) call Logger%Write( "ExcType    = ", ExcType )
    if (i_Debug_Loc) call Logger%Write( "NProcCurrA = ", NProcCurrA )
    if (i_Debug_Loc) call Logger%Write( "NProcCurrB = ", NProcCurrB )
    if (i_Debug_Loc) call Logger%Write( "NProcPre   = ", NProcPre )
    if (i_Debug_Loc) call Logger%Write( "This%NProc_iPOpp(iP,2) = ", This%NProc_iPOpp(iP,2) )

    call CreateName_Nb4Atoms( MolNameA, MolNameB, iLevel, jLevel, iLevelChar, jLevelChar, Name )  
    
    iLevelFin     = [iLevel,     jLevel]
    iLevelFinChar = [iLevelChar, jLevelChar]
    Idx           = NProcPre + ( (NProcCurrA-1)*This%NProc_iPOpp(iP,2) + NProcCurrB )
    Pairs         = [iP, iOpp]

    if (i_Debug_Loc) call Logger%Write( "iLevelFin  = ", iLevelFin )
    if (i_Debug_Loc) call Logger%Write( "Idx        = ", Idx )
    if (i_Debug_Loc) call Logger%Write( "Pairs      = ", Pairs )

  else
    allocate(Name, source = adjustl(trim( Collision%Atoms(1)%Name // '+' // Collision%Atoms(2)%Name // '+' // Collision%Atoms(3)%Name // '+' // Collision%Atoms(4)%Name )) )
    iLevelFin     = [  0,   0]
    iLevelFinChar = ['0', '0']
    Idx           = 0
    Pairs         = [0, 0]
  end if
 
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   READING INITIAL AND FINAL CONDITIONS AND CROSS SECTIONS FROM statistics.csv
! ==============================================================================================================
Subroutine Convert_CrossSect_To_Rates_Nb4Atoms( This, Input, Collision, Velocity, i_Debug, i_Debug_Deep)

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero
  use Sorting_Module
  
  class(Nb4Atoms_Processes_Type)                        ,intent(inout) ::    This
  type(Input_Type)                                      ,intent(in)    ::    Input
  type(Collision_Type)                                  ,intent(inout) ::    Collision
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
  integer                                                   ::    IssueIn, IssueFin
  logical                                                   ::    i_Debug_Loc
  logical                                                   ::    i_Debug_Loc_Deep


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Convert_CrossSect_To_Rates_Nb4Atoms" )
  i_Debug_Loc_Deep = i_Debug_Global; if ( present(i_Debug_Deep) )i_Debug_Loc_Deep = i_Debug_Deep
  !i_Debug_Loc   =     Logger%On()
  

  !do iTtra = 1,Input%NTtra
  !  if (i_Debug_Loc) call Logger%Write( "  iTtra = ", iTtra )
    NTtra = 1 !Input%NTtra
    iTtra = 1

    !do iTint = 1,Input%NTint
    !  if (i_Debug_Loc) call Logger%Write( "  iTint = ", iTint )
      NTint = 1
      iTInt = 1

      FileName = trim(adjustl(Input%LevelOutputDir))// '/statistics.csv'
      if (i_Debug_Loc) call Logger%Write( "  Reading File: ", FileName )
      Open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
        if (Status/=0) then
          if (i_Debug_Loc) call Logger%Write( "The statistics.csv File is not present for this Initial Level / Bin." )
          FileName = adjustl(trim( trim(adjustl(This%LevelOutputDir)) // '/StatIssues.csv' ))
          if ( i_Debug_Loc ) call Logger%Write( "Writing File: ", FileName )
          open( File=FileName, NewUnit=This%UnitIssues, status='unknown', access='append', iostat=Status )
            if (Status/=0) call Error( "Error writing the binary data file for Rates: " // FileName  ) 
            write(This%UnitIssues, '(A)') '#    vIn(:),    jIn(:),     ArrIn,   IssueIn,   vFin(:),   jFin(:),    ArrFin,  IssueFin'          
            write(This%UnitIssues, '(X, I9, *(A, I9))') -1, ',', 0, ',', 0, ',', 0, ',', 0, ',', 0, ',', 0, ',', 0, ',', 0, ',', 0
          close(This%UnitIssues)
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


            !!! Checking if Initial and Final States are Accettable
            call Mask4InProc_Nb4Atoms(  Collision, Input%BinOI,   vqnIn,  jqnIn,  ArrIn,  IssueIn  )!, i_Debug=i_Debug_Loc )
            call Mask4FinProc_Nb4Atoms( Collision, CrossSectTemp, vqnFin, jqnFin, ArrFin, IssueFin )!, i_Debug=i_Debug_Loc )
            if ( (IssueIn < 1) .and. (IssueFin < 1) ) then

              !!! Postprocessing Final State for Reconstructing Process !!!!
              call This%FindingFinalLevel( Input, Collision, vqnFin, jqnFin, ArrFin, Name, ProcType, ExcType, iP, iLevelFin, iLevelFinChar, Idx, i_Debug=i_Debug_Loc )

              !!! New Process ??? !!!!
              Proc_To_Line = This%Proc_To_LineVec(Idx-1)
              if (Proc_To_Line < 1) then
                !!! New Process! Allocating it !!!!
                This%NProc_Cleaned          = This%NProc_Cleaned + 1
                call This%ProcessesVecTemp(This%NProc_Cleaned)%Shelving_1stTime( 2, NTtra, Idx, Name, ProcType, ExcType, iP, iLevelFin, iLevelFinChar, CorrFactor=1.0, CrossSect=CrossSectTemp, Velocity=Velocity(iTtra), i_Debug=i_Debug_Loc )
                This%Proc_To_LineVec(Idx-1) = This%NProc_Cleaned
                
              else
                !!! Old Process! Adding Cross Section !!!!
                call This%ProcessesVecTemp(Proc_To_Line)%Shelving( iTtra, CorrFactor=1.0, CrossSect=CrossSectTemp, Velocity=Velocity(iTtra), i_Debug=i_Debug_Loc )
              end if

            else
              if (i_Debug_Loc_Deep) call Logger%Write( "    Found an Issue with the Trajectory; writing it in a Separate File." )  
              call This%WritingIssue( vqnIn,  jqnIn,  ArrIn,  IssueIn, vqnFin, jqnFin, ArrFin, IssueFin, i_Debug=i_Debug_Loc )

            end if

          end do


          !!! Creating an Array of Unique Processes and Deallocating Temporary One
          allocate( This%ProcessesVecCleaned(This%NProc_Cleaned), Stat=Status  )
          if (Status/=0) call Error( "Error allocating ProcessesVecCleaned in Convert_CrossSect_In_Rates_Nb4Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Allocated ProcessesVecCleaned with Dimension = (",This%NProc_Cleaned,")" )

          allocate( This%IdxVec(This%NProc_Cleaned), Stat=Status  )
          if (Status/=0) call Error( "Error allocating IdxVec in Convert_CrossSect_In_Rates_Nb4Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Allocated IdxVec with Dimension = (",This%NProc_Cleaned,")" )
          allocate( This%IdxVecSorted(This%NProc_Cleaned), Stat=Status  )
          if (Status/=0) call Error( "Error allocating IdxVecSorted in Convert_CrossSect_In_Rates_Nb4Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Allocated IdxVecSorted with Dimension = (",This%NProc_Cleaned,")" )
          

          do iProc=1,This%NProc_Cleaned
            This%IdxVec(iProc)              = This%ProcessesVecTemp(iProc)%Idx
            This%ProcessesVecCleaned(iProc) = This%ProcessesVecTemp(iProc)
          end do

          call hpsort(This%IdxVec, This%IdxVecSorted)

          deallocate(This%ProcessesVecTemp, Stat=Status)
          if (Status/=0) call Error( "Error deallocating This%ProcessesVecTemp in Convert_CrossSect_In_Rates_Nb4Atoms" )
          if (i_Debug_Loc) call Logger%Write( "    Deallocated This%ProcessesVecTemp" )

        end if
      Close(Unit)
      if (i_Debug_Loc) call Logger%Write( "    Done Reading File: ", FileName )
      

      !!! Writing Rates !!!!
      if (Input%PostWritesBinaryFlg) then
        call This%WritingRates_Binary( iTTra, iTInt, i_Debug=i_Debug_Loc )
      else
        call This%WritingRates( iTTra, iTInt, i_Debug=i_Debug_Loc )
      end if


      if ( allocated(This%IdxVec) ) then
        deallocate(This%IdxVec, Stat=Status)
        if (Status/=0) call Error( "Error deallocating This%IdxVec in Convert_CrossSect_In_Rates_Nb4Atoms" )
        if (i_Debug_Loc) call Logger%Write( "    Deallocated This%IdxVec" )
        deallocate(This%IdxVecSorted, Stat=Status)
        if (Status/=0) call Error( "Error deallocating This%IdxVecSorted in Convert_CrossSect_In_Rates_Nb4Atoms" )
        if (i_Debug_Loc) call Logger%Write( "    Deallocated This%IdxVecSorted" )
      end if

    !end do

  !end do


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! !---------------------------------------------------------------------------------------------------! !
!!!                           Creating a Mask for Statistic Results                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Pure Subroutine Mask4InProc_Nb4Atoms( Collision, BinOI, vqnIn, jqnIn, ArrIn, Issue)!, i_Debug )

  use Collision_Class             ,only: Collision_Type

  Type(Collision_Type)                              ,intent(in)     :: Collision
  integer     ,dimension(2)                         ,intent(in)     :: BinOI
  integer     ,dimension(2)                         ,intent(inout)  :: vqnIn
  integer     ,dimension(2)                         ,intent(inout)  :: jqnIn
  integer                                           ,intent(inout)  :: ArrIn
  integer                                           ,intent(out)    :: Issue
  ! logical                                 ,optional ,intent(in)     :: i_Debug
  
  integer                                                           :: Temp
  integer                                                           :: iPIn, iPFin, jPIn, jPFin
  integer                                                           :: iTypeIn, iTypeFin, jTypeIn, jTypeFin
  integer                                                           :: iMolIn, iMolFin, jMolIn, jMolFin
  integer                                                           :: iBinIn, iBinFin, iBinTemp, jBinIn, jBinFin, jBinTemp
  logical                                                           :: i_Debug_Loc

  ! i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  ! if (i_Debug_Loc) call Logger%Entering( "Mask4InProc_Nb4Atoms" )
  ! !i_Debug_Loc   =     Logger%On() 

  Issue = -1

  iPIn     = int(ArrIn / 16.0_rkp)
  jPIn     = Collision%Pairs(iPIn)%Opposite
  Temp     = mod(ArrIn, 16)
  iTypeIn  = mod(Temp, 4)
  jTypeIn  = int(Temp / 4)
  iMolIn   = Collision%Pairs(iPIn)%To_Molecule
  jMolIn   = Collision%Pairs(jPIn)%To_Molecule
  
  iBinIn   = BinOI(iMolIn)
  jBinIn   = BinOI(jMolIn)
  iBinTemp = Collision%MoleculesContainer(iMolIn)%Molecule%BinsContainer%qns_to_Bin(vqnIn(1),jqnIn(1))
  jBinTemp = Collision%MoleculesContainer(jMolIn)%Molecule%BinsContainer%qns_to_Bin(vqnIn(2),jqnIn(2))


  if ( ( iTypeIn > 1 ) .or. ( jTypeIn > 1 ) ) then
    Issue = 1 ! Current Trajectory Starts from Dissociation
  elseif ( ( Collision%Pairs(iPIn)%To_Molecule /= iMolIn ) .or. ( Collision%Pairs(jPIn)%To_Molecule /= jMolIn ) ) then
    Issue = 2 ! Current Trajectory does not have the molecule of interest as Initial Condition
  elseif ( ( iBinTemp == -1 ) .or. ( jBinTemp == -1 ) ) then
    Issue = 3 ! Current Trajectory has an Initial Condition in which the Q.N.s are not listed in the iMol-th Molecule Energy Levels
  elseif ( ( iBinTemp /= iBinIn ) .or. ( jBinTemp /= jBinIn ) ) then
    Issue = 4 ! Current Trajectory Starts from a Level/Bin that is not the one of Interest
  else
    ! Current Trajectory has an Initial Condition that is accettable!
    vqnIn(1) = Collision%MoleculesContainer(iMolIn)%Molecule%BinsContainer%Bin(iBinIn)%vqnFirst
    jqnIn(1) = Collision%MoleculesContainer(iMolIn)%Molecule%BinsContainer%Bin(iBinIn)%jqnFirst
    vqnIn(2) = Collision%MoleculesContainer(jMolIn)%Molecule%BinsContainer%Bin(jBinIn)%vqnFirst
    jqnIn(2) = Collision%MoleculesContainer(jMolIn)%Molecule%BinsContainer%Bin(jBinIn)%jqnFirst
  end if

  
  ! if (i_Debug_Loc) call Logger%Exiting

End Subroutine
! !===================================================================================================! !


! !---------------------------------------------------------------------------------------------------! !
!!!                           Creating a Mask for Statistic Results                                   !!!
! !---------------------------------------------------------------------------------------------------! !
Subroutine Mask4FinProc_Nb4Atoms( Collision, CrossSect, vqnFin, jqnFin, ArrFin, Issue)!, i_Debug )

  use Collision_Class             ,only: Collision_Type

  Type(Collision_Type)                              ,intent(in)     :: Collision
  real(rkp)   ,dimension(2)                         ,intent(in)     :: CrossSect
  integer     ,dimension(2)                         ,intent(inout)  :: vqnFin
  integer     ,dimension(2)                         ,intent(inout)  :: jqnFin
  integer                                           ,intent(inout)  :: ArrFin
  integer                                           ,intent(out)    :: Issue
  ! logical                                 ,optional ,intent(in)     :: i_Debug
  
  integer                                                           :: Temp
  integer                                                           :: iPIn, iPFin, jPIn, jPFin
  integer                                                           :: iTypeIn, iTypeFin, jTypeIn, jTypeFin
  integer                                                           :: iMolIn, iMolFin, jMolIn, jMolFin
  integer                                                           :: iBinIn, iBinFin, iBinTemp, jBinIn, jBinFin, jBinTemp
  logical                                                           :: i_Debug_Loc

  ! i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  ! if (i_Debug_Loc) call Logger%Entering( "Mask4FinProc_Nb4Atoms" )
  ! !i_Debug_Loc   =     Logger%On() 

  Issue = -1

  iPFin    = int(ArrFin / 16)
  if (CrossSect(1) < 1.e-100) then
    Issue = 12
  else
    if (iPFin > 0) then
      jPFin    = Collision%Pairs(iPFin)%Opposite
      Temp     = mod(ArrFin , 16)
      iTypeFin = mod(Temp, 4)
      jTypeFin = int(Temp / 4)
      iMolFin  = Collision%Pairs(iPFin)%To_Molecule
      jMolFin  = Collision%Pairs(jPFin)%To_Molecule


      !iBinTemp  = Collision%MoleculesContainer(iMolFin)%Molecule%BinsContainer%qns_to_Bin(vqnFin(1),0)
      !jBinTemp  = Collision%MoleculesContainer(jMolFin)%Molecule%BinsContainer%qns_to_Bin(vqnFin(2),0)
      iBinTemp  = Collision%MoleculesContainer(iMolFin)%Molecule%BinsContainer%qns_to_Bin(vqnFin(1),jqnFin(1))
      jBinTemp  = Collision%MoleculesContainer(jMolFin)%Molecule%BinsContainer%qns_to_Bin(vqnFin(2),jqnFin(2))


      if (iTypeFin < 2) then
        if ( iBinTemp .eq. -1 ) then
          if ( Collision%MoleculesContainer(iMolFin)%Molecule%BinsContainer%qns_to_Bin(vqnFin(1),jqnFin(1)-1) .eq. -1 ) then
            Issue = 0 ! Current Trajectory has a Final Condition that is very close to the Centrifugal Barrier. Its lifetime will be very short and for this reason the Final State is considered DISSOCIATED
            vqnFin(1) = 0
            jqnFin(1) = 0
            iTypeFin  = 2
          else
            Issue = 11 ! Current Trajectory has a Final Condition that should not exist
          endif
        else
          ! Current Trajectory has a Final Condition that is accettable!
          iBinFin   = iBinTemp
          vqnFin(1) = Collision%MoleculesContainer(iMolFin)%Molecule%BinsContainer%Bin(iBinFin)%vqnFirst
          jqnFin(1) = Collision%MoleculesContainer(iMolFin)%Molecule%BinsContainer%Bin(iBinFin)%jqnFirst
        end if
      end if

      if (jTypeFin < 2) then
        if ( jBinTemp .eq. -1 ) then
          if ( Collision%MoleculesContainer(jMolFin)%Molecule%BinsContainer%qns_to_Bin(vqnFin(2),jqnFin(2)-1) .eq. -1 ) then
            Issue = 0 ! Current Trajectory has a Final Condition that is very close to the Centrifugal Barrier. Its lifetime will be very short and for this reason the Final State is considered DISSOCIATED
            vqnFin(2) = 0
            jqnFin(2) = 0
            jTypeFin  = 2
          else
            Issue = 12 ! Current Trajectory has a Final Condition that should not exist
          endif
        else
          ! Current Trajectory has a Final Condition that is accettable!
          jBinFin   = jBinTemp
          vqnFin(2) = Collision%MoleculesContainer(jMolFin)%Molecule%BinsContainer%Bin(jBinFin)%vqnFirst
          jqnFin(2) = Collision%MoleculesContainer(jMolFin)%Molecule%BinsContainer%Bin(jBinFin)%jqnFirst
        end if
      end if

      ArrFin = int( iPFin*16 + 4*jTypeFin + iTypeFin )
    end if
  end if
  
  ! if (i_Debug_Loc) call Logger%Exiting

End Subroutine
! !===================================================================================================! !


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
  character(:)  ,allocatable ,intent(out) :: Name

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

  allocate( Name, source=trim(adjustl( trim(adjustl(MolNameA)) // "(" // trim(adjustl(iLevelChar)) // ")+" // trim(adjustl(MolNameB)) // "(" // trim(adjustl(jLevelChar)) // ")" )) )


end Subroutine
! !===================================================================================================! !


End Module