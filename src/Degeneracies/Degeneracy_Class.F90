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

Module Degeneracy_Class

#include "../qct.inc"

  use Parameters_Module       ,only:  rkp, Zero, One, Six, Ue, UKb
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error

  implicit none

  private
  public  ::    Degeneracy_Type
  public  ::    Compute_PartFunc_State
  
  Type    ,abstract                                         ::    Degeneracy_Type
    logical                                                 ::    Initialized         !< Indicator whether the object is initialized
    character(:)                              ,allocatable  ::    Name                !< Name of current Degeneracy
  contains
    private
    procedure                                ,public        ::    Initialize             =>    Initialize_Degeneracy
    procedure                                ,public        ::    Output                 =>    Output_Degeneracy    
    procedure(Compute)     ,deferred ,nopass ,public        ::    Compute_Degeneracy_State
  End Type
  
  Abstract Interface
  
    Pure Elemental Function Compute( jqn ) result( g )
      use Parameters_Module      ,only:  rkp, Zero, One, Two, Three, Six
      integer                         ,intent(in)     ::    jqn
      real(rkp)                                       ::    g
    End Function

  End interface

  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains

Subroutine Initialize_Degeneracy( This, Input, i_Debug )

  use Input_Class                 ,only:  Input_Type
  
  class(Degeneracy_Type)                    ,intent(out)    ::    This
  Type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Degeneracy")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  This%Name         =   '<Unknown>'
  This%Initialized  =   .True.
  
  if (i_Debug_Loc) write(Logger%Unit,"(4x,'[Initialize_Degeneracy]: Nothing to do here')")
  if (i_Debug_Loc) call Logger%Exiting()
  
End Subroutine


Subroutine Output_Degeneracy( This, Unit )
  class(Degeneracy_Type)                    ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    Unit
  integer                                                   ::    idum
  logical                                                   ::    ldum
  idum = Unit
  ldum = This%Initialized
End Subroutine


Pure Function Compute_PartFunc_State( g, EeV, T ) result( Q )

  real(rkp)                                 ,intent(in)     ::    g
  real(rkp)                                 ,intent(in)     ::    EeV
  real(rkp)                                 ,intent(in)     ::    T
  real(rkp)                                                 ::    Q

  Q = g * exp( - EeV * Ue / (UKb * T) ) 

End Function


End Module