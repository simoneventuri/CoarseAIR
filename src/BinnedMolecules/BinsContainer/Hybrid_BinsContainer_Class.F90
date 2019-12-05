! -*-F90-*
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Berkan Bolkan and Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
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
Module Hybrid_BinsContainer_Class

#include "../../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, One, Two
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Hybrid_BinsContainer_Type


  Type    ,extends(BinsContainer_Type)          ::    Hybrid_BinsContainer_Type
    
  contains
    procedure          ::  Initialize            =>    Initialize_Hybrid_BinsContainer    
    procedure          ::  Compute_Level_To_Bin  =>    Compute_Level_To_Bin_Hybrid_BinsContainer  
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.                                                    
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                               DEFERRED PROCEDURES for Hybrid BinsContainer
! **************************************************************************************************************
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Hybrid_BinsContainer( This, Input, LevelsContainer, iMol, WriteFlg, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(Hybrid_BinsContainer_Type)            ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    WriteFlg
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    NBins = 0
  integer                                                     ::    TempN1, TempN2, TempN3
  real(rkp)                                                   ::    Temp
  integer                                                     ::    iLevels
  integer                                                     ::    iBins, jBins
  integer                                                     ::    iBins_temp
  integer                                                     ::    ijqn
  integer                                                     ::    iwrite
  character(:)                                   ,allocatable ::    FileName
  character(:)                                   ,allocatable ::    FileName2
  character(5)                                                ::    iBins_char
  character(6)                                                ::    NBins_char
  integer                                                     ::    Status
  integer                                                     ::    Unit
  logical                                                     ::    VibSpecFlg

  logical                                                     ::    i_Debug_Loc
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Hybrid_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()


  ! Reading from "/E_low.dat" the Bins Energy Minima
  FileName = adjustl(trim(This%BinDataFile(1)))
  if (i_Debug_Loc)  call Logger%Write( "Reading the File for the Bins Energy Minima" )
  if (i_Debug_Loc)  call Logger%Write( "-> Opening file: ", FileName )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )
  read(Unit,*,iostat=Status)
  do while (Status==0)  
    read(Unit,*,iostat=Status) iBins, Temp
    if (Status==0) then
      This%Bin(iBins)%Hybrid    = 1
      This%Bin(iBins)%MineinteV = Temp
      NBins                     = NBins + 1
    end if 
  end do
  close(Unit)
  
  
  ! Reading from "/E_low.dat" the Bins Energy Minima
  FileName = adjustl(trim(This%BinDataFile(2)))
  if (i_Debug_Loc)  call Logger%Write( "Reading the File for the Bins Energy Minima" )
  if (i_Debug_Loc)  call Logger%Write( "-> Opening file: ", FileName )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )
  read(Unit,*,iostat=Status)
  do while (Status==0)  
    read(Unit,*,iostat=Status) iBins, TempN1, TempN2, TempN3
    if (Status==0) then
      if (This%Bin(iBins)%Hybrid == 1) call Error( "Error in the Hybrid Mthd File: Bin already used with specified Energy!" )
      This%Bin(iBins)%Hybrid   = 2
      This%Bin(iBins)%vqnFirst = TempN1
      This%Bin(iBins)%jqnFirst = TempN2
      This%Bin(iBins)%jqnLast  = TempN3
      NBins                    = NBins + 1
    end if
  end do
  close(Unit)
  
  if (NBins /= Input%NBins(iMol)) call Error( "Error: Nb of Bins in Input and Nb of Bins in Hybrid Mthd Files Differ!" )
  
  ! Copying "/E_low.dat" to the Output Folder
  call system('scp ' // adjustl(trim(This%BinDataFile(1))) // ' ' // This%PathToBinMolFldr )
  call system('scp ' // adjustl(trim(This%BinDataFile(2))) // ' ' // This%PathToBinMolFldr )


  call This%Compute_Level_To_Bin( LevelsContainer, i_Debug=i_Debug_Loc )


  ! Finding the first (v,J)s for each of the bins (vectors vqnFirst and jqnFirst)
  ! Assigning the correspondence (v,J)_level -> bin through the matrix qns_to_Bin
  ! Counting the number of levels per bin (vector NLevels)
  ! Writing all the bin vqns and jqns in the related Bin_*.dat file
  This%Bin%NLevels = 0
  if (i_Debug_Loc)  call Logger%Write( "Finding the first (v,J)s for each of the bins (vectors vqnFirst and jqnFirst) ... ")
  do iBins = 1,Input%NBins(iMol)
    
    iwrite = 1
    do iLevels = 1,LevelsContainer%NStates 
         
      if (LevelsContainer%States(iLevels)%To_Bin == iBins) then
         
        if (iwrite == 1) then
          This%Bin(iBins)%MinLevEintev = LevelsContainer%States(iLevels)%einteV
          This%Bin(iBins)%vqnFirst     = LevelsContainer%States(iLevels)%vqn
          This%Bin(iBins)%jqnFirst     = LevelsContainer%States(iLevels)%jqn
          iwrite = 2
        end if
         
       This%qns_to_Bin( LevelsContainer%States(iLevels)%vqn, LevelsContainer%States(iLevels)%jqn ) = iBins
       This%Bin(iBins)%NLevels = This%Bin(iBins)%NLevels + 1
       
      end if
       
    end do
   
  end do
  if (i_Debug_Loc)  then
    call Logger%Write( "        Bin    Min En [eV]    Nb Levels  Levels %")
    do iBins = 1,Input%NBins(iMol)
      call Logger%Write( "         ", iBins, This%Bin(iBins)%MinLevEinteV,  This%Bin(iBins)%NLevels, real(This%Bin(iBins)%NLevels) * 1.e2_rkp / real(LevelsContainer%NStates) )
    end do
  end if 


  if (WriteFlg) then
  
    call This%WriteNLevels( i_Debug=i_Debug_Loc )

    
    call This%WriteBinsData( Input, LevelsContainer, iMol )

  end if

  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_Level_To_Bin_Hybrid_BinsContainer( This, LevelsContainer, i_Debug )

  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(Hybrid_BinsContainer_Type)            ,intent(inout)  ::    This
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    iLevels, iBins, jBins
  logical                                                     ::    VibSpecFlg
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Compute_Level_To_Bin_Hybrid_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()
  

  do iLevels = 1,LevelsContainer%NStates
    
    VibSpecFlg = .false.
    
    do iBins = 1,This%NBins
      if (This%Bin(iBins)%Hybrid == 2) then
        if (LevelsContainer%States(iLevels)%vqn == This%Bin(iBins)%vqnFirst) then
          if ( (LevelsContainer%States(iLevels)%jqn >= This%Bin(iBins)%jqnFirst) .and. (LevelsContainer%States(iLevels)%jqn <= This%Bin(iBins)%jqnLast) ) then
            LevelsContainer%States(iLevels)%To_Bin = iBins
            VibSpecFlg = .true.
            exit
          end if
        end if
      end if
    end do
    
    if (.not. VibSpecFlg) then
      do iBins = 1,This%NBins
        jBins = iBins
        if ((LevelsContainer%States(iLevels)%einteV > This%Bin(iBins)%MineinteV) .and. (LevelsContainer%States(iLevels)%einteV < This%Bin(iBins+1)%MineinteV)) exit
        jBins = This%NBins
      end do 
      LevelsContainer%States(iLevels)%To_Bin = jBins  
    end if
    
  end do


  if (i_Debug_Loc) call Logger%Exiting

end Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module