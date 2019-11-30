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

Module System_Class

#include "../qct.inc"

  use Parameters_Module       ,only:  rkp, Zero, One, Six, Ue, UKb
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error

  implicit none

  private
  public  ::    System_Type
  
  Type    ,abstract                                         ::    System_Type
    logical                                                 ::    Initialized         !< Indicator whether the object is initialized
    character(:)                              ,allocatable  ::    Name                !< Name of current System
  contains
    private
    procedure              ,public                          ::    Initialize                   =>    Initialize_System
    procedure              ,public                          ::    AssignPairsArrangements      =>    AssignPairsArrangements_System
    procedure              ,public                          ::    AdjustArrangements           =>    AdjustArrangements_System

  End Type

  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains



Subroutine Initialize_System( This, Input, i_Debug )

  use Input_Class                 ,only:  Input_Type
  
  class(System_Type)                        ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_System" )
  !i_Debug_Loc   =     Logger%On()
  
  This%Name         =   '<Unknown>'
  This%Initialized  =   .True.
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine



Subroutine AssignPairsArrangements_System( This, Collision, Input, i_Debug )

  use Collision_Class             ,only: Collision_Type
  use Input_Class                 ,only: Input_Type

  class(System_Type)                        ,intent(in)     ::    This
  type(Collision_Type)                      ,intent(inout)  ::    Collision
  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AssignPairsArrangements_System" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting               

End Subroutine



Subroutine AdjustArrangements_System( This, Collision, Input, i_Debug )

  use Collision_Class             ,only: Collision_Type
  use Input_Class                 ,only:  Input_Type

  class(System_Type)                        ,intent(in)     ::    This
  type(Collision_Type)                      ,intent(inout)  ::    Collision
  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc 
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AdjustArrangements_System" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting              
  
End Subroutine


End Module