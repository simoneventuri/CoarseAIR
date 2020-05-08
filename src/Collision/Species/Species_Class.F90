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

Module Species_Class

  use Parameters_Module         ,only:  rkp, Zero
  use Logger_Class              ,only:  Logger
  use Atom_Class                ,only:  Atom_Type
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type
  use LevelsContainer_Class     ,only:  LevelsContainer_Type

  implicit none

  private
  public    ::    Species_Type

  Type      ::    Species_Type
    integer                                         ::    Idx         =   0           ! Index of species in the list of all Species
    character(:)  ,allocatable                      ::    Name                        ! Name of species
    character(:)  ,allocatable                      ::    BSortMethod                 ! 
    real(rkp)                                       ::    Mass        =   Zero        ! Mass of the species [a.u.]
    real(rkp)                                       ::    RedMass     =   Zero        ! Reduced mass of the species [a.u.]
    integer                                         ::    NAtoms      =   0           ! Number of atomss contained in the species
    integer                                         ::    To_Molecule =   0
    integer         ,dimension(:)   ,allocatable    ::    To_Atoms                    ! Index mapping from current Species to the Atoms its contains in the list of all atoms. 
                                                                                      ! This list corresponds to the variable Collision%Atoms(1:NAtoms), where NAtoms is the total number of atoms, 
                                                                                      ! including atoms from both the target and projectile species.
    type(Atom_Type) ,dimension(:)   ,allocatable    ::    Atoms                       ! List of Atoms object, each one corresponds to a given atom of the current species. Dim=(NAtoms)
    class(DiatomicPotential_Type)   ,allocatable    ::    DiatPot                     ! Intra-molecular diatomic potenitla object
    type(LevelsContainer_Type)                      ::    ListStates
  contains
    private
    procedure ,public                   ::    Initialize => InitializeSpecies
    procedure ,public ,non_overridable  ::    AngularMomentum
    procedure ,public ,non_overridable  ::    GetAtomsMass
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeSpecies( This, Input, iSpecies, iAtoms, Name, Atoms, i_Debug )
!!! The 'Atoms' object has the 'inout' attribute since the component corresponding to the index mapping from the atom to the associated species is being set.
! @TODO: Do not defined the 'DiatPot' object as a component. Instead, add a pointer to the Pairs object. This will be usefull for species with more than 2 atoms

  use Input_Class                        ,only:  Input_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type

  class(Species_Type)                       ,intent(out)    ::    This
  class(Input_Type)                         ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iSpecies      ! Index of current species
  integer         ,dimension(:)             ,intent(in)     ::    iAtoms        ! Index of the atoms contained in current species (1 or 2 for now)
  character(*)                              ,intent(in)     ::    Name          ! Name of current species
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms         ! List of Atoms object corresponding to the atoms contained in current species
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    iA, jA            ! Index of atoms
  integer                                                   ::    iMol
  type(DiatomicPotential_Factory_Type)                      ::    DiatPotFactory
  real(rkp)                                                 ::    SumMass
  real(rkp)                                                 ::    ProMass

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeSpecies")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  

  if (i_Debug_Loc) call Logger%Write( "Initializing species: ", trim(Name) )
  This%Idx            =   iSpecies
  allocate( This%Name , source =   trim(Name) )                                                                 ! Setting the species name
  This%NAtoms         =   size( iAtoms )                                                                        ! Setting the number of atoms
  allocate( This%Atoms(   This%NAtoms) )                                                                        ! Setting the array of Atom objects
  allocate( This%To_Atoms(This%NAtoms) )
  do iA = 1,This%NAtoms
    jA                =   iAtoms(iA)                                                                            ! Getting the global index of current atom
    This%Atoms(iA)    =   Atoms(jA)                                                                             ! Setting the array of Atom objects
    This%To_Atoms(iA) =   Atoms(jA)%Idx                                                                         ! Setting the index mapping from current species to the total list of atoms
  end do

  if (This%NAtoms == 2) then
    do iMol = 1,Input%NMolecules
      if ( trim(adjustl(Input%Molecules_Name(iMol))) == trim(adjustl(This%Name)) ) This%To_Molecule = iMol
    end do
    allocate( This%BSortMethod, source = trim(adjustl(Input%BSortMethod(This%To_Molecule))) )
  end if

  SumMass             =   sum(     This%Atoms(:)%Mass )                                                         ! Computing the sum of the atoms mass
  ProMass             =   product( This%Atoms(:)%Mass )                                                         ! Computing the product of the atoms mass
  This%Mass           =   SumMass                                                                               ! Setting the species mass [a.u.]
  This%RedMass        =   ProMass / SumMass                                                                     ! Setting the reduced species mass [a.u.]

! ==============================================================================================================
!   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiatPotFactory%Construct" )
  call DiatPotFactory%Construct( Atoms, This%To_Atoms, Input, This%DiatPot, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
! ==============================================================================================================


  if (i_Debug_Loc) then
    call Logger%Write( "This%Idx          = ", This%Idx           )
    call Logger%Write( "This%Name         = ", This%Name          )
    call Logger%Write( "This%DiatPot%Name = ", This%DiatPot%Name  )
    call Logger%Write( "This%NAtoms       = ", This%NAtoms        )
    call Logger%Write( "This%Mass         = ", This%Mass    , Fr="es15.8" )
    call Logger%Write( "This%RedMass      = ", This%RedMass , Fr="es15.8" )
    call Logger%Write( "This%To_Molecule  = ", This%To_Molecule   )
    call Logger%Write( "This%To_Atoms     = ", This%To_Atoms      )
    call Logger%Write( "This%BSortMethod  = ", This%BSortMethod   )
    do iA = 1,This%NAtoms
    associate( Ato => This%Atoms(iA) )
      call Logger%Write( "-> Ato%Idx = ", Ato%Idx, "Ato%Name = ", Ato%Name, "Ato%Mass = ", Ato%Mass, F2="i1", F4="a5", Fr=",s15.8" )
    end associate
    end do
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function AngularMomentum( This, Q, dQdt ) result(AngMon)
! COMPUTING THE ANGULAR MOMENTUM OF THE TARGET

  class(Species_Type)                       ,intent(in)     ::    This
  real(rkp) ,dimension(:,:)                 ,intent(in)     ::    Q             ! Species coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  real(rkp) ,dimension(:,:)                 ,intent(in)     ::    dQdt          ! Time derivatives of the target coordinates for each spatial direction (dim-1) and for each atom (dim-2). Dim=(NSpace=3,NAtoms)
  
  real(rkp) ,dimension(3)                                   ::    AngMon
  integer                                                   ::    i               ! Index of atoms in the target/projectile species
  real(rkp)                                                 ::    Mass            ! Mass of a given atom of current species
  
  AngMon        =   Zero
  do i = 1,This%NAtoms
    Mass        =   This%Atoms(i)%Mass
    AngMon(1)   =   AngMon(1) + (Q(2,i) * dQdt(3,i) - Q(3,i) * dQdt(2,i)) * Mass
    AngMon(2)   =   AngMon(2) + (Q(3,i) * dQdt(1,i) - Q(1,i) * dQdt(3,i)) * Mass
    AngMon(3)   =   AngMon(3) + (Q(1,i) * dQdt(2,i) - Q(2,i) * dQdt(1,i)) * Mass
  end do
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function GetAtomsMass( This ) result(AtomsMass)

  class(Species_Type)                       ,intent(in)     ::    This
  
  real(rkp) ,dimension(This%NAtoms)                         ::    AtomsMass
  
  AtomsMass(:)  =   This%Atoms(:)%Mass
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
