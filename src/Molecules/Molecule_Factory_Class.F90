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

Module Molecule_Factory_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Molecule_Factory_Type

  Type      ::    Molecule_Factory_Type
  contains
    private
    procedure ,nopass ,public ::  Define_Molecule
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains

Subroutine Define_Molecule( Input, NPairs, Pairs, Atoms, iMol, Molecule, i_Debug )

  use Input_Class            ,only:    Input_Type
  use Atom_Class             ,only:    Atom_Type
  use AtomsPair_Class        ,only:    AtomsPair_Type
  use Molecule_Class         ,only:    Molecule_Type
  use Nb2_Molecule_Class     ,only:    Nb2_Molecule_Type                                                             
!  use Nb3_Molecule_Class     ,only:    Nb3_Molecule_Type                                                               

  type(Input_Type)                                  ,intent(inout)  ::    Input
  integer                                           ,intent(in)     ::    NPairs
  type(AtomsPair_Type) ,dimension(:)                ,intent(inout)  ::    Pairs                        ! List of Pairs object
  type(Atom_Type)      ,dimension(:)                ,intent(in)     ::    Atoms                        ! List of Atoms object: Atoms of the target species are listed first, 
  integer                                           ,intent(in)     ::    iMol
  class(Molecule_Type)                 ,allocatable ,intent(out)    ::    Molecule
  logical                                 ,optional ,intent(in)     ::    i_Debug

  logical                                                           ::    i_Debug_Loc
  integer                                                           ::    Status

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Define_Molecule")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  select case ( Input%Molecule_NAtoms(iMol) )    
    case(2)
      if (i_Debug_Loc) call Logger%Write( "Defining a Nb2_Molecule_Type object" )
      allocate( Nb2_Molecule_Type :: Molecule )
!    case(3)
!      if (i_Debug_Loc) call Logger%Write( "Defining a Nb3_Molecule_Type object" )
!      allocate( Nb3_Molecule_Type :: Molecule )
    case default
      call Error( "Nb of Atoms per Molecule not supported (YET): Check Input%2Atoms_Molecule_Type.")
  end select

  if (i_Debug_Loc) call Logger%Write( "Calling Molecule%Initialize" )
  call Molecule%Initialize( Input, NPairs, Pairs, Atoms, iMol, i_Debug=i_Debug_Loc )

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine


End Module