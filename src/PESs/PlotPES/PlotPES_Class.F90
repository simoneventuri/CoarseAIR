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

Module PlotPES_Class
  
#include "../../qct.inc"

  use Parameters_Module           ,only:  rkp, Zero, Half, One, Two, Pi, Kelvin_To_Hartree, Kcm_To_Hartree, Hartree_To_eV, B_To_Ang, KcmAng_To_HartB
  use Logger_Class                ,only:  Logger, LogLevel_INFO
  use Error_Class                 ,only:  Error

  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Transformation_Class        ,only:  R_to_X, X_to_R, dX_To_dR, dR_To_dX

  implicit none

  private
  public  ::    PlotPES_Type
  
  Type    ,abstract     :: PlotPES_Type
    logical                 ::    Initialized         !< Indicator whether the object is initialized
    real(rkp)               ::    RConverter     = One
    real(rkp)               ::    VConverter     = One
    real(rkp)               ::    dVConverter    = One
  contains
    private
    procedure              ,public                          ::    Initialize             =>    PlotPES_Initialize
    procedure              ,public                          ::    IsoTri                 =>    PlotPES_IsoTri
    procedure              ,public                          ::    Rot3rd                 =>    PlotPES_Rot3rd
    procedure              ,public                          ::    ComputeCuts            =>    PlotPES_ComputeCuts
    procedure              ,public                          ::    EvaluatePoints         =>    PlotPES_EvaluatePoints
    procedure              ,public                          ::    PlotsVargasPaper       =>    PlotPES_PlotsVargasPaper
    procedure              ,public                          ::    StochPESStats          =>    PlotPES_StochPESStats
    procedure              ,public                          ::    GridForScatter         =>    PlotPES_GridForScatter
    procedure              ,public                          ::    ReadPoints             =>    PlotPES_ReadPoints
    procedure              ,public                          ::    TripleGrid             =>    PlotPES_TripleGrid
    procedure              ,public                          ::    DoubleGrid             =>    PlotPES_DoubleGrid
    procedure              ,public                          ::    GridForStochPES        =>    PlotPES_GridForStochPES
    procedure              ,public                          ::    Grid                   =>    PlotPES_Grid
  End Type

  integer   ,parameter    ::    NSpace         = 3
  logical   ,parameter    ::    Formatted      = .True.
  integer                 ::    iSpeTar        = 1
  integer                 ::    iSpePro        = 2
  real(rkp)               ::    MostProbEr                            
  real(rkp)               ::    rVMin_Min      = 2.0d0
  real(rkp)               ::    rVMin_Max      = 2.5d0
  
  logical   ,parameter    ::    i_Debug_Global = .False.
  
  contains


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_Initialize( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_Initialize" )
  !i_Debug_Loc   =     Logger%On()


  if (trim(adjustl(Input%UnitDist)) .eq. 'Angstrom') then           
    This%RConverter  = One              / B_To_Ang
    This%dVConverter = This%dVConverter / B_To_Ang
  end if
  
  if (trim(adjustl(Input%UnitPot)) .eq. 'KcalMol') then                                                                           
    This%VConverter  = One              / Kcm_To_Hartree
    This%dVConverter = This%dVConverter / Kcm_To_Hartree
  elseif (trim(adjustl(Input%UnitPot)) .eq. 'ElectronVolt') then                                                                
    This%VConverter  = One              * Hartree_To_eV
    This%dVConverter = This%dVConverter * Hartree_To_eV
  end if
  if (i_Debug_Loc) call Logger%Write( "This%RConverter  = ", This%RConverter )
  if (i_Debug_Loc) call Logger%Write( "This%VConverter  = ", This%VConverter )
  if (i_Debug_Loc) call Logger%Write( "This%dVConverter = ", This%dVConverter )

  
  call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES' )
  if (i_Debug_Loc) call Logger%Write( "Created PlotPES Output Folder" )

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_Grid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_Grid" )
  !i_Debug_Loc   =     Logger%On()

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_GridForStochPES( This, Input, Collision, NPairs, NAtoms,  i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_GridForStochPES" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_DoubleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_DoubleGrid" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_TripleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_TripleGrid" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_ReadPoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_ReadPoints" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_GridForScatter( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_GridForScatter" )
  !i_Debug_Loc   =     Logger%On()
        
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_StochPESStats( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_StochPESStats" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_PlotsVargasPaper( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_PlotsVargasPaper" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_EvaluatePoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_EvaluatePoints" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_ComputeCuts( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_ComputeCuts" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_Rot3rd( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_Rot3rd" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_IsoTri( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_IsoTri" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module