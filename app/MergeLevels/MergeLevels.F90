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

! ==============================================================================================================
!           Merging Levels Lists for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! This program merges the Levels Lists Produced by qct-computelevels.x.
!
!   Mandatory Input Arguments: - Translational Temperature [K]
!                              - Internal Temperature [K]
!   Optional Input Arguments:  - Nb of Binned Molecules = Input%NBinnedMolecules
!                                   - Nb of Bins for BinnedMolecule(iMol),          iMol=1:NBinnedMolecules
!                                   - Current Bin for the BinnedMolecule(iMol),     iMol=1:NBinnedMolecules
!
! ==============================================================================================================
Program MergeLevels

  use Parameters_Module        ,only:  rkp, Zero
  use Logger_Class             ,only:  Logger
  use Input_Class              ,only:  Input_Type
  use Collision_Class          ,only:  Collision_Type
  use LevelsGenerator_Module   ,only:  ReadEnergyLevels
  use Error_Class              ,only:  Error

  implicit none

  type(Input_Type)                                  ::    Input
  type(Collision_Type)                              ::    Collision

  character(200)                                    ::    LevelOutputDir
  character(200)                                    ::    FileName
  integer                                           ::    iMol
  integer                                           ::    Status
  real(rkp)                                         ::    StartTime, EndTime

  logical   ,parameter                              ::    i_Debug_ML        = .True.
  logical   ,parameter                              ::    i_Debug_ML_Deep   = .True.
  
  call CPU_Time( StartTime )
  
  if (i_Debug_ML) call Logger%Initialize( "MergeLevels.log",       &                                               ! Opening the Log File using
              Status          =       'REPLACE',                   &                                               ! replacing any previous log file
              Position        =       'REWIND',                    &                                               ! rewinding to the top of the file
              Procedure       =       'MergeLevels',               &                                               ! loading the calling procedure name
              Indentation     =       2           )                                                                ! and setting the initial indentation level


! ==============================================================================================================
!  READING THE PROGRAM ARGUMENT
! ==============================================================================================================            
  call getarg( 1, LevelOutputDir )
  if (i_Debug_ML) call Logger%Write( "LevelOutputDir = ", LevelOutputDir )

  call getarg( 2, Input%NProc_char )
  if (trim(adjustl(Input%NProc_char)) .eq. '') then
    Input%NProc_char  = '1'
    Input%NProc       =  1
  else 
    read( Input%NProc_char, '(I3)' ) Input%NProc
  end if
  if (i_Debug_ML) call Logger%Write( "Input%NProc_char = ", Input%NProc_char, "Input%NProc = ", Input%NProc )
  
  call getarg( 3, Input%NNode_char )
  if (trim(adjustl(Input%NNode_char)) .eq. '') then
    Input%NNode_char  = '1'
    Input%NNode       =  1
  else 
    read( Input%NNode_char, '(I3)' ) Input%NNode
  end if
  if (i_Debug_ML) call Logger%Write( "Input%NNode_char = ", Input%NNode_char, "Input%NNode = ", Input%NNode )


! ==============================================================================================================
!   INITIALIZING INPUT
! ==============================================================================================================
  if (i_Debug_ML) call Logger%Write( "Reading the input. Calling Input%Initialize", NewLine=.True. )
  call Input%Initialize( i_Debug=i_Debug_ML )
  if (i_Debug_ML) call Logger%Write( "Done initializing Input" )
  if (i_Debug_ML) call Logger%Write( "Calling Input%GenerateLevels" )
  call Input%GenerateLevels( i_Debug=i_Debug_ML )
  Input%TaskType = 2
  if (i_Debug_ML) call Logger%Write( "Done reading Input for GenerateLevels" )
  
  Input%LevelOutputDir = trim(adjustl(LevelOutputDir))
  if (i_Debug_ML) call Logger%Write( "This%LevelOutputDir = ", Input%LevelOutputDir )
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING COLLISION OBJECT
! ==============================================================================================================
  if (i_Debug_ML) call Logger%Write( "Initializing the Collision object", NewLine=.True. )
  if (i_Debug_ML) call Logger%Write( "Calling Collision%Initialize" )
  call Collision%Initialize( Input, i_Debug=i_Debug_ML, i_Debug_Deep=i_Debug_ML_Deep )
  if (i_Debug_ML) call Logger%Write( "Done with Collision%Initialize" )
! ==============================================================================================================


  do iMol = 1,Input%NMolecules
    if ( size(Collision%MoleculesContainer(iMol)%Molecule%To_Atoms) == 1 ) cycle
    if (i_Debug_ML) call Logger%Write( "Checking whether to run iMol ", iMol, ", ", Collision%MoleculesContainer(iMol)%Molecule%Name )

    if (iMol == Input%ComputeLevels(iMol)) then 
      if (i_Debug_ML) call Logger%Write( "I am going to run iMol ", iMol, ", ", Collision%MoleculesContainer(iMol)%Molecule%Name, "; Calling ComputeEnergyLevels" )

      call ReadEnergyLevels(Input, Collision, iMol, i_Debug=i_Debug_ML ) 
      
      ! call system('scp ' // &
      !             Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr  // trim(adjustl(Input%GeneratedLevelsFile(iMol))) // &
      !             ' ' //
      !             Collision%MoleculesContainer(iMol)%Molecule%PathToMolDtbFldr // trim(adjustl(Input%GeneratedLevelsFile(iMol))) )

    end if

  end do
  
  if (i_Debug_ML) call Logger%Write( "Normal termination" )

End Program
