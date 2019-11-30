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

Module Molecule_Class

#include "../qct.inc"

  use Parameters_Module       ,only:  rkp, Zero, One, Six, Ue, UKb
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error

  use Atom_Class              ,only:  Atom_Type
  use DiatomicPotential_Class ,only:  DiatomicPotential_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type
  use BinsContainer_Class     ,only:  BinsContainer_Type

  implicit none

  private
  public  ::    Molecule_Type
  
  Type    ,abstract                                          ::    Molecule_Type
    logical                                                  ::    Initialized               !< Indicator whether the object is initialized
    character(:)                                ,allocatable ::    Name                      !< Name of current Degeneracy
    character(:)                                ,allocatable ::    PathToMolFldr
    character(:)                                ,allocatable ::    PathToMolDtbFldr
    integer                       ,dimension(:) ,allocatable ::    To_Atoms
    class(Atom_Type)              ,dimension(:) ,allocatable ::    Atom                      ! List of Atoms object: Atoms of the target species are listed first, 
    real(rkp)                                                ::    xmu
    real(rkp)                                                ::    xmui
    real(rkp)                                                ::    xmui2
    integer                       ,dimension(:) ,allocatable ::    To_Pairs
    class(DiatomicPotential_Type)               ,allocatable ::    DiaPot                    ! Intra-molecular diatomic potenitla object
    class(LevelsContainer_Type)                 ,allocatable ::    LevelsContainer
    class(BinsContainer_Type)                   ,allocatable ::    BinsContainer
  contains
    private
    procedure ,public        ::  Initialize                     =>    Initialize_Molecule
  End Type

  logical   ,parameter  ::    i_Debug_Global = .False.

  contains



Subroutine Initialize_Molecule( This, Input, NPairs, Pairs, Atoms, iMol, i_Debug )

  use Input_Class                    ,only:  Input_Type
  use Atom_Class                     ,only:  Atom_Type
  use AtomsPair_Class                ,only:  AtomsPair_Type
  use DiatomicPotential_Factory_Class ,only:  DiatomicPotential_Factory_Type

  class(Molecule_Type)                      ,intent(out)    ::    This
  Type(Input_Type)                          ,intent(inout)  ::    Input
  integer                                   ,intent(in)     ::    NPairs
  type(AtomsPair_Type) ,dimension(:)        ,intent(inout)  ::    Pairs                        ! List of Pairs object
  type(Atom_Type)      ,dimension(:)        ,intent(in)     ::    Atoms                        ! List of Atoms object: Atoms of the target species are listed first, 
  integer                                   ,intent(in)     ::    iMol
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Molecule" )
  !i_Debug_Loc   =     Logger%On()
  
  This%Name         =   '<Unknown>'
  This%Initialized  =   .True.
  
  if (i_Debug_Loc) write(Logger%Unit,"(4x,'[Initialize_Molecule]: Nothing to do here')")
  
  if (i_Debug_Loc) call Logger%Exiting()
  
End Subroutine


End Module