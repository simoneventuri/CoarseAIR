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

Module Atom_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger

  implicit none

  private
  public    ::    Atom_Type

  Type      ::    Atom_Type
    character(:)  ,allocatable  ::    Name                      ! Name of atom
    integer                     ::    Idx         =   0         ! Global index of the atoms in the list of all atoms. This list corresponds to the variable Collision%Atoms(1:NAtoms), where NAtoms is the total number of atoms, 
                                                                !    including atoms from both the target and projectile species.
    real(rkp)                   ::    Mass        =   Zero      ! Mass of atom [a.u.]
    integer                     ::    To_Species  =   0         ! Index mapping from current Atom to the element in the list of Species object which corresponds to the species containing current atom
  contains
    private
    procedure ,public   ::    Initialize => Initialize_Atom
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Atom( This, Mass, Name, i_Debug )

  class(Atom_Type)                          ,intent(out)    ::    This
  real(rkp)                                 ,intent(in)     ::    Mass      ! Mass of the Atom [a.u.]
  character(*)                              ,intent(in)     ::    Name      ! Name of Atom
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Atom")
  !i_Debug_Loc   =     Logger%On()

  This%Mass   =   Mass
  allocate( This%Name , source = trim(Name) )

  if (i_Debug_Loc) then
    call Logger%Write( "This%Mass  = ", This%Mass, Fr="es15.8" )
    call Logger%Write( "This%Name  = ", This%Name )
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
