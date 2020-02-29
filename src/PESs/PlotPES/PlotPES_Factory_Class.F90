! -*-F90-*-
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

Module PlotPES_Factory_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    PlotPES_Factory_Type

  Type      ::    PlotPES_Factory_Type
  contains
    private
    procedure ,nopass ,public ::  Define_PlotPEs
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains

Subroutine Define_PlotPEs( Input, PlotPES, i_Debug )

  use Input_Class             ,only:    Input_Type
  use PlotPES_Class           ,only:    PlotPES_Type
  use Nb3_PlotPES_Class       ,only:    Nb3_PlotPES_Type                                                             
  use Nb4_PlotPES_Class       ,only:    Nb4_PlotPES_Type                                                               

  type(Input_Type)                                  ,intent(in)     ::    Input
  class(PlotPES_Type)                  ,allocatable ,intent(out)    ::    PlotPES
  logical                                 ,optional ,intent(in)     ::    i_Debug

  logical                                                           ::    i_Debug_Loc
  integer                                                           ::    Status

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Define_PlotPES")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  select case ( Input%NAtoms )    
    case(3)
      if (i_Debug_Loc) call Logger%Write( "Defining a Nb3_PlotPES_Type object" )
      allocate( Nb3_PlotPES_Type :: PlotPES )
    case(4)
      if (i_Debug_Loc) call Logger%Write( "Defining a Nb4_PlotPES_Type object" )
      allocate( Nb4_PlotPES_Type :: PlotPES )
    case default
      call Error( "System Nb of Atoms not supported (YET): Check Input%NAtoms. " )
  end select


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine

End Module