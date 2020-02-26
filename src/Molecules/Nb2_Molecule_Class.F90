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

Module Nb2_Molecule_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Six
  use Molecule_Class        ,only:  Molecule_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Nb2_Molecule_Type


  Type    ,extends(Molecule_Type)                 ::    Nb2_Molecule_Type
  contains
    procedure         ::  Initialize                     =>    Initialize_Nb2_Molecule
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.
  
  contains

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for 2Atoms Molecule
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Initialize_Nb2_Molecule( This, Input, NPairs, Pairs, Atoms, iMol, i_Debug )

  use Input_Class                    ,only:  Input_Type
  use Atom_Class                     ,only:  Atom_Type
  use AtomsPair_Class                ,only:  AtomsPair_Type
  use DiatomicPotential_Factory_Class ,only:  DiatomicPotential_Factory_Type
  use BinsContainer_Factory_Class     ,only:  BinsContainer_Factory_Type

  class(Nb2_Molecule_Type)                  ,intent(out)    ::    This
  Type(Input_Type)                          ,intent(inout)  ::    Input
  integer                                   ,intent(in)     ::    NPairs
  type(AtomsPair_Type) ,dimension(:)        ,intent(inout)  ::    Pairs                        ! List of Pairs object
  type(Atom_Type)      ,dimension(:)        ,intent(in)     ::    Atoms                        ! List of Atoms object: Atoms of the target species are listed first, 
  integer                                   ,intent(in)     ::    iMol
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    i, iP, i1, i2
  integer                                                   ::    Status
  integer                                                   ::    iA1, iA2
  integer                                                   ::    pos_m, pos_e
  character(:)                    ,allocatable              ::    AtA_Name, AtB_Name
  character(:)                    ,allocatable              ::    MolA_Name, MolB_Name
  integer       ,dimension(6)                               ::    To_Pairs_Temp = 0
  type(DiatomicPotential_Factory_Type)                      ::    DiaPot_Factory
  type(BinsContainer_Factory_Type)                          ::    BinsContainer_Factory
  character(150)                                            ::    FileName
  logical                                                   ::    WriteFlg
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Nb2_Molecule")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  
  This%Initialized  =   .True.
  allocate( This%Name, source = adjustl(trim( Input%Molecules_Name(iMol) )) )  
  if (i_Debug_Loc) call Logger%Write( "This%Name = ", This%Name) 


  ! ==============================================================================================================
  !   7.1. MAPPING MOLECULES TO PAIRS AND VICEVERSA 
  ! ==============================================================================================================
  To_Pairs_Temp = 0
  if (i_Debug_Loc) call Logger%Write( "To_Pairs_Temp = ", To_Pairs_Temp )
  do iP = 1,NPairs

     ! Form molecule name for current pair. Here one must be account for the fact the NO is the same as ON. 
     ! Atom names
     pos_m = index(Pairs(iP)%Name, '-')
     pos_e = len_trim(Pairs(iP)%Name)

     AtA_Name = Pairs(iP)%Name(1:pos_m - 1)
     AtB_Name = Pairs(iP)%Name(pos_m + 1:pos_e)
 
     ! Homonuclear molecule (e.g., N2, O2)
     if (AtA_Name .eq. AtB_Name) then 

        MolA_Name = AtA_Name//'2'
        
        if (This%Name .eq. MolA_Name) then 
           Pairs(iP)%To_Molecule = iMol
           if (i_Debug_Loc) call Logger%Write( "Pair Nb", iP, " Corresponds to Molecule: ", This%Name) 
           To_Pairs_Temp(iP) = 1
         endif

     ! Heteronucler molecule (e.g., NO, CO)
     else

        MolA_Name = AtA_Name//AtB_Name
        MolB_Name = AtB_Name//AtA_Name

        if ( (This%Name .eq. MolA_Name) .or. (This%Name .eq. MolB_Name) ) then 
          Pairs(iP)%To_Molecule = iMol
          if (i_Debug_Loc) call Logger%Write( "Pair Nb", iP, " Corresponds to Molecule: ", This%Name) 
          To_Pairs_Temp(iP) = 1
        endif

     endif

  end do 
  if (i_Debug_Loc) call Logger%Write( "To_Pairs_Temp = ", To_Pairs_Temp )
  if ( sum(To_Pairs_Temp) == 0 ) then
    call Error( "No Pair Contains the Molecule" )
  else
    if (i_Debug_Loc) call Logger%Write( "Found ", sum(To_Pairs_Temp), " Pairs Containing the Molecule ", This%Name )
    allocate( This%To_Pairs( sum(To_Pairs_Temp) ), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%To_Pairs in AssignMoleculesToPairs_Nb2_Molecule" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%To_Pairs with Dimension = (",sum(To_Pairs_Temp),")" )
    i=1
    do iP=1,6
      if ( To_Pairs_Temp(iP) == 1 ) then
        This%To_Pairs(i) = iP
        i = i+1
        if (i_Debug_Loc) call Logger%Write( "Pair Nb ", iP, " Contains Molecule ", This%Name )
      end if
    end do
  end if
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   7.2. ALLOCATING ATOMS AND COMPUTING REDUCED MASS
  ! ==============================================================================================================
  allocate( This%To_Atoms(2), stat=Status )
  if (Status/=0) call Error( "Error allocating This%To_Atoms" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%To_Atoms with dimension 2" )

  allocate( This%Atom(2), stat=Status )
  if (Status/=0) call Error( "Error allocating This%Atom" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Atom with dimension 2" )

  This%To_Atoms(1) = Pairs(This%To_Pairs(1))%To_Atoms(1)
  This%To_Atoms(2) = Pairs(This%To_Pairs(1))%To_Atoms(2)
  if (i_Debug_Loc) call Logger%Write( "First  Atom of the Molecule ", This%Name, " is the Atom Nb ", This%To_Atoms(1) )
  if (i_Debug_Loc) call Logger%Write( "Second Atom of the Molecule ", This%Name, " is the Atom Nb ", This%To_Atoms(2) )

  i1 = Pairs(This%To_Pairs(1))%To_Atoms(1)
  i2 = Pairs(This%To_Pairs(1))%To_Atoms(2)
  This%Atom = Atoms
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Atom(1): This%Atom(1)%Name = ", This%Atom(1)%Name )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Atom(2): This%Atom(2)%Name = ", This%Atom(2)%Name )

  This%xmu = (This%Atom(1)%Mass * This%Atom(2)%Mass) / (This%Atom(1)%Mass + This%Atom(2)%Mass) 
  if (i_Debug_Loc) call Logger%Write( "Computed Reduced Mass; This%xmu = ", This%xmu )
  This%xmui    = One  / This%xmu
  This%xmui2   = Half * This%xmui
  if (i_Debug_Loc) call Logger%Write( "Computed Inverse of Reduced Mass; This%xmui = ", This%xmui, "; This%xmui2 = ", This%xmui2 )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   7.3. CONSTRUCTING THE DIATOMIC POTENTIAL ASSOCIATED TO THE TWO ATOMS WITHIN THE CURRENT MOLECULE
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Construction the diatomic potential object" )
  call DiaPot_Factory%Construct( This%Atom, [1,2], Input, This%DiaPot, i_Debug=i_Debug_Loc ) 
  if (i_Debug_Loc) call Logger%Write( "Diatomic pot. name: This%DiaPot%Name = ", This%DiaPot%Name )
  if (i_Debug_Loc) call Logger%Write( "-> Done" )
  ! ==============================================================================================================


  ! ==============================================================================================================
  !   7.4. CREATING LOCAL FOLDER FOR THE MOLECULE
  ! ==============================================================================================================
  allocate( This%PathToMolDtbFldr, source = adjustl(trim( trim(adjustl(Input%DtbPath))   // '/Molecules/' // This%Name                   // '/' // trim(adjustl(This%DiaPot%Name)) // '/' )) )
  allocate( This%PathToMolFldr,    source = adjustl(trim( trim(adjustl(Input%OutputDir)) // '/'           // trim(adjustl(Input%System)) // '/' // trim(adjustl(This%Name))        // '/' )) )
  if (i_Debug_Loc) call Logger%Write( "Path to Folder for Original      Molecule, This%PathToMolDtbFldr = ", This%PathToMolDtbFldr )
  if (i_Debug_Loc) call Logger%Write( "Path to Folder for Pre-Processed Molecule, This%PathToMolFldr    = ", This%PathToMolFldr    )

  call system('mkdir -p ' // This%PathToMolFldr )
  if (i_Debug_Loc) call Logger%Write( "Created Folder for Pre-Processed Molecule through system mkdir" )
  ! ==============================================================================================================


  if (Input%TaskType > 2) then
    
    
    ! ==============================================================================================================
    !   7.5. CONSTRUCTING THE LEVELS CONTAINER
    ! ==============================================================================================================
    if (Input%TaskType == 3) then
      call system('scp ' // This%PathToMolDtbFldr // '/' // adjustl(trim( Input%LevelsFileName(iMol) )) // ' ' // This%PathToMolFldr )
      if (i_Debug_Loc) call Logger%Write( "Levels List copied to: ", This%PathToMolFldr )

      !!! NOTE: IF WE WANT TO INCLUDE THE OPTION OF CUTTING THE LEVEL LIST, THEN WE SHOULD MODIFY THIS !!!
      !!! WE SHOULD :
      !!!    - Read the Original List, only if E_Level > EMin, defined in Input File
      !!!    - Compute Partition Function Ratio
      !!!    - Write Partitio Ratio in A file
      !!!    - Write List in This%PathToMolFldr as '/levels_cut.inp'
      !!!
      call system('scp ' // This%PathToMolDtbFldr // '/' // adjustl(trim( Input%LevelsFileName(iMol) )) // ' ' // This%PathToMolFldr // '/levels_cut.inp' )
      if (i_Debug_Loc) call Logger%Write( "Levels List copied as: ", This%PathToMolFldr // '/levels_cut.inp' )
      ! ==============================================================================================================
    end if
    FileName  = This%PathToMolFldr // '/levels_cut.inp'

    if (i_Debug_Loc) call Logger%Write( "Calling This%Levels_Container%InitializeLevelsContainer" )
    allocate(This%LevelsContainer); call This%LevelsContainer%Initialize( Input, iMol, FileName=FileName, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) call Logger%Write( "-> This%Levels_Container%InitializeLevelsContainer" )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   7.6 CONSTRUCTING THE BINS CONTAINER
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Calling This%BinsContainer%Initialize" )

    if (Input%TaskType == 3) then
      WriteFlg = .True.
    else
      WriteFlg = .False.
    end if
    call BinsContainer_Factory%Construct( Input, This%LevelsContainer, iMol, This%PathToMolDtbFldr, This%PathToMolFldr, This%BinsContainer, WriteFlg=WriteFlg, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) call Logger%Write( "-> This%BinsContainer%Initialize" )
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   7.7 COMPUTING / READING THE PARTITION FUNCTIONS
    ! ==============================================================================================================
    if (Input%TaskType == 3) then
      if (i_Debug_Loc) call Logger%Write( "Calling This%BinsContainer%ComputePartFunEnergy" )
      call This%BinsContainer%ComputePartFunEnergy( Input, This%LevelsContainer, iMol, i_Debug=i_Debug_Loc )
      if (i_Debug_Loc) call Logger%Write( "-> This%BinsContainer%ComputePartFunEnergy" )
    else 
      if (i_Debug_Loc) call Logger%Write( "Calling This%BinsContainer%ReadPartFunEnergy" )
      call This%BinsContainer%ReadPartFunEnergy( Input, This%LevelsContainer, iMol, i_Debug=i_Debug_Loc )
      if (i_Debug_Loc) call Logger%Write( "-> This%BinsContainer%ReadPartFunEnergy" )
    end if
    ! ==============================================================================================================

  end if

  
  if (i_Debug_Loc) call Logger%Exiting()  

End Subroutine


End Module