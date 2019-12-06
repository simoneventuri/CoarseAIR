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
Module Energy_BinsContainer_Class

#include "../../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, One, Two
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Energy_BinsContainer_Type


  Type    ,extends(BinsContainer_Type)          ::    Energy_BinsContainer_Type
    
  contains
    procedure          ::  Initialize            =>    Initialize_Energy_BinsContainer    
    procedure          ::  Compute_Level_To_Bin  =>    Compute_Level_To_Bin_Energy_BinsContainer  
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.                                                    
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                               DEFERRED PROCEDURES for Energy BinsContainer
! **************************************************************************************************************
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Energy_BinsContainer( This, Input, LevelsContainer, iMol, WriteFlg, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(Energy_BinsContainer_Type)            ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    WriteFlg
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  integer                                                     ::    iLevels
  integer                                                     ::    iBins
  integer                                                     ::    iBins_temp
  integer                                                     ::    ijqn
  integer                                                     ::    iwrite
  character(:)                                   ,allocatable ::    FileName
  character(:)                                   ,allocatable ::    FileName2
  character(5)                                                ::    iBins_char
  integer                                                     ::    Status
  integer                                                     ::    Unit
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Energy_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()

  
  ! Reading from "/E_low.dat" the Bins Energy Minima
  FileName = adjustl(trim(This%BinDataFile(1)))
  if (i_Debug_Loc)  call Logger%Write( "Reading the File for the Bins Energy Minima" )
  if (i_Debug_Loc)  call Logger%Write( "-> Opening file: ", FileName )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )
  read(Unit,*)
  do iBins = 1,Input%NBins(iMol)
    read(Unit,*) iBins_temp, This%Bin(iBins)%MineinteV
  end do
  close(Unit)
  
  ! Copying "/E_low.dat" to the Output Folder
  call system('scp ' // adjustl(trim(This%BinDataFile(1))) // ' ' // This%PathToBinMolFldr )
      

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
Subroutine Compute_Level_To_Bin_Energy_BinsContainer( This, LevelsContainer, i_Debug )

  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(Energy_BinsContainer_Type)            ,intent(inout)  ::    This
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    iLevels, iBins, jBins
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Compute_Level_To_Bin_Energy_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()
  
  if (This%NBins == 1) then
    LevelsContainer%States(:)%To_Bin = 1
  else
    do iLevels = 1,LevelsContainer%NStates
      iBins = 1
      do while ( LevelsContainer%States(iLevels)%einteV > This%Bin(iBins+1)%MineinteV )
        iBins = iBins+1
        if ( iBins == This%NBins ) exit 
      end do
      LevelsContainer%States(iLevels)%To_Bin = iBins
    end do
  end if
  

  if (i_Debug_Loc) call Logger%Exiting

end Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module