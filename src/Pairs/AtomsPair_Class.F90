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

Module AtomsPair_Class

  use Parameters_Module         ,only:  rkp, Zero, One, Two, Pi, Kelvin_To_Hartree, Rugc
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public    ::    AtomsPair_Type

  Type      ::    AtomsPair_Type
    character(:)                        ,allocatable  ::    Name                      ! Name of the atom-atom pair: A-B where A and B are the atoms names
    integer                                           ::    Idx         =   0         ! Index of the atom-atom pair
    integer                                           ::    Opposite    =   0         ! Index of the 'opposite' pair. For 3 atoms system, it is always 0. For 4 atoms system, it is the index of the pair which does not contain any of
                                                                                      !     the atoms of current pair
    integer             ,dimension(2)                 ::    To_Atoms    =   0         ! Index mapping from current pair to the two associated atoms in the list of all atoms.
    integer                                           ::    To_Species  =   0         ! Index mapping from current pair to the species initially containing the two atoms. 
                                                                                      !     Equal 0 if the pair is connecting two atoms from different species.
    integer                                           ::    To_Molecule               ! Index mapping from current pair to the Molecule Nb ...
    integer                                           ::    To_BinnedMolecule         ! Index mapping from current pair to the Binned Molecule Nb ...
    integer                                           ::    NLevels     =   0
    class(DiatomicPotential_Type)       ,allocatable  ::    DiaPot                    ! Intra-molecular diatomic potenitla object
  contains
    private
    procedure ,public                   ::    Initialize => Initialize_AtomsPair
  End Type

  logical   ,parameter                  ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_AtomsPair( This, Input, iPair, iAtoms, Atoms, Opposite, i_Debug )

  use Input_Class                     ,only:  Input_Type
  use Atom_Class                      ,only:  Atom_Type
  use DiatomicPotential_Factory_Class  ,only:  DiatomicPotential_Factory_Type

  class(AtomsPair_Type)                     ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iPair               ! Index of current pair
  integer         ,dimension(2)             ,intent(in)     ::    iAtoms              ! Index of the two atoms of the current pair
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms               ! Array of Atoms objects used to construct the current pair
  integer                         ,optional ,intent(in)     ::    Opposite
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer         ,dimension(2)                             ::    iA = [1,2]
  integer                                                   ::    jA, kA
  logical                                                   ::    i_Debug_Loc
  character(:)    ,allocatable                              ::    Atom1Name, Atom2Name
  type(DiatomicPotential_Factory_Type)                       ::    DiaPotFactory

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_AtomsPair")
  !i_Debug_Loc   =     Logger%On()

  allocate( Atom1Name, source = adjustl(trim( Atoms(iAtoms(1))%Name)) )
  allocate( Atom2Name, source = adjustl(trim( Atoms(iAtoms(2))%Name)) )
  if (i_Debug_Loc) call Logger%Write( "Atom1Name = ", Atom1Name)
  if (i_Debug_Loc) call Logger%Write( "Atom2Name = ", Atom2Name)

  if (i_Debug_Loc) call Logger%Write( "Ordering Alphabetically the Names of the Atoms Composing the Pair")
  if (Atom1Name <= Atom2Name) then
    jA = 1; kA = 2
    if (i_Debug_Loc) call Logger%Write( "Alphabetically, 1st Atom Name <= 2nd Atom Name" )
  else
    jA = 2; kA = 1
    if (i_Debug_Loc) call Logger%Write( "Alphabetically, 1st Atom Name > 2nd Atom Name; Inverting the 2 Atoms" )
  end if
  
  associate( Atom1 => Atoms(iAtoms(jA)), Atom2 => Atoms(iAtoms(kA)) )
  
    This%Idx          =   iPair
    This%To_Atoms     =   [iAtoms(jA), iAtoms(kA)]
    if (i_Debug_Loc) call Logger%Write( "This%Idx      = ", This%Idx )
    if (i_Debug_Loc) call Logger%Write( "This%To_Atoms = ", This%To_Atoms )

    allocate( This%Name , source = Atom1%Name // '-' // Atom2%Name )
    if (i_Debug_Loc) call Logger%Write( "This%Name = ", This%Name )

    This%To_Species   =   0
    if ( Atoms(jA)%To_Species == Atoms(kA)%To_Species ) This%To_Species = Atoms(jA)%To_Species
    if ( present(Opposite) ) This%Opposite = Opposite

! ==============================================================================================================
!   CONSTRUCTING THE DIATOMIC POTENTIAL ASSOCIATED TO THE TWO ATOMS WITHIN THE CURRENT PAIR
! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Construction the diatomic potential object" )
    call DiaPotFactory%Construct( [Atom1,Atom2], iA, Input, This%DiaPot, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) call Logger%Write( "-> Done" )
! ==============================================================================================================

  end associate

  if (i_Debug_Loc) then
    call Logger%Write( "Index of pair:      This%Idx         = ", This%Idx         )
    call Logger%Write( "Name of pair:       This%Name        = ", This%Name        )
    call Logger%Write( "Mapping to atoms:   This%To_Atoms    = ", This%To_Atoms    )
    call Logger%Write( "Mapping to species: This%To_Species  = ", This%To_Species  )
    call Logger%Write( "Diatomic pot. name: This%DiaPot%Name = ", This%DiaPot%Name )
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module