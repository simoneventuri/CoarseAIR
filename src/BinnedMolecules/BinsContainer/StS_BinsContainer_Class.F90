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
Module StS_BinsContainer_Class

#include "../../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, One, Two
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    StS_BinsContainer_Type


  Type    ,extends(BinsContainer_Type)          ::    StS_BinsContainer_Type
    
  contains
    procedure          ::  Initialize            =>    Initialize_StS_BinsContainer    
    procedure          ::  Compute_Level_To_Bin  =>    Compute_Level_To_Bin_StS_BinsContainer
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.                                                    
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                               DEFERRED PROCEDURES for StS BinsContainer
! **************************************************************************************************************
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Subroutine Initialize_StS_BinsContainer( This, Input, LevelsContainer, iMol, WriteFlg, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(StS_BinsContainer_Type)               ,intent(inout)  ::    This
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
  logical                                                     ::    i_Debug_Deepest
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_StS_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()


  call This%Compute_Level_To_Bin( LevelsContainer, i_Debug=i_Debug_Loc )
  
  ! Scaling Levels
  LevelsContainer%States%einteV_scaled = LevelsContainer%States%einteV - LevelsContainer%MineinteV
  

  ! Creating qnFirst and qns_to_Bin
  This%Bin%NLevels = 1
  do iBins = 1,This%NBins

    This%Bin(iBins)%vqnFirst     = LevelsContainer%States(iBins)%vqn
    This%Bin(iBins)%jqnFirst     = LevelsContainer%States(iBins)%jqn
    
    This%qns_to_Bin( LevelsContainer%States(iBins)%vqn, LevelsContainer%States(iBins)%jqn ) = iBins    
  end do


  if (WriteFlg) then

    call This%WriteNLevels( i_Debug=i_Debug_Loc )


    call This%WriteBinsData( Input, LevelsContainer, iMol )
  
  end if

  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_Level_To_Bin_StS_BinsContainer( This, LevelsContainer, i_Debug )

  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(StS_BinsContainer_Type)               ,intent(inout)  ::    This
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    iLevels, iBins, jBins
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Compute_Level_To_Bin_StS_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()


  do iLevels = 1,LevelsContainer%NStates
    LevelsContainer%States(iLevels)%To_Bin = iLevels
  end do
  

  if (i_Debug_Loc) call Logger%Exiting

end Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



End Module