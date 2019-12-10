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
Module Processes_Factory_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Processes_Factory_Type

  Type      ::    Processes_Factory_Type
  contains
    private
    procedure ,nopass ,public ::  Define_Processes
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains

Subroutine Define_Processes( Input, Collision, Processes, i_Debug )

  use Input_Class                ,only:    Input_Type
  use Collision_Class            ,only:    Collision_Type
  use Processes_Class            ,only:    Processes_Type
  use Nb3Atoms_Processes_Class   ,only:    Nb3Atoms_Processes_Type                                                             
  use Nb4Atoms_Processes_Class   ,only:    Nb4Atoms_Processes_Type                                                               

  type(Input_Type)                                  ,intent(in)     ::    Input
  Type(Collision_Type)                              ,intent(inout)  ::    Collision
  class(Processes_Type)                ,allocatable ,intent(out)    ::    Processes
  logical                                 ,optional ,intent(in)     ::    i_Debug

  logical                                                           ::    i_Debug_Loc
  integer                                                           ::    Status

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Define_Processes")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  select case ( Input%NAtoms )    
    case(3)
      if (i_Debug_Loc) call Logger%Write( "Defining a Nb3Atoms_Processes_Type object" )
      allocate( Nb3Atoms_Processes_Type :: Processes )
    case(4)
      if (i_Debug_Loc) call Logger%Write( "Defining a Nb4Atoms_Processes_Type object" )
      allocate( Nb4Atoms_Processes_Type :: Processes )
    case default
      call Error( "Nb of Atoms not supported (YET): Check Input%NAtoms. " )
  end select

  if (i_Debug_Loc) call Logger%Write( "Calling Processes%Initialize" )
  call Processes%Initialize( Input, Collision, i_Debug=i_Debug_Loc)

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine


End Module