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

! ==================================================================================================================
!           Preprocessing Energy Levels for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==================================================================================================================
! The program sorts, cuts and (if required) groups the energy levels of the desired molecules.
!
!   Mandatory Input Arguments: - Translational Temperature [K]
!                              - Internal Temperature [K]
!
! ==================================================================================================================
Program PreprocLevels

  use Parameters_Module     ,only:  rkp, Zero, One, Hartree_To_eV
  use Logger_Class          ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class           ,only:  Error
  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Timing_Module

  implicit none

  type(Input_Type)                                     :: Input
  type(Collision_Type)                                 :: Collision

  integer                                              :: iMol
  integer                                              :: Status
  real(rkp)                                            :: StartTime, EndTime
  
  logical   ,parameter                                 ::     i_Debug_PL        = .True.
  logical   ,parameter                                 ::     i_Debug_PL_Deep   = .True.

  call CPU_Time( StartTime )

  if (i_Debug_PL) call Logger%Initialize( "PreprocLevels.log", &                                ! Opening the Log File using
               Status          =       'REPLACE',           &                                   ! replacing any previous log file
               Position        =       'REWIND',            &                                   ! rewinding to the top of the file
               Procedure       =       'PreprocLevels',     &                                   ! loading the calling procedure name
               Indentation     =       2           )                                            ! and setting the initial indentation level
  
! ==============================================================================================================
!  READING THE PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 1, Input%Ttra_char )
  read(Input%Ttra_char, "(d20.10)", iostat=Status) Input%Ttra
  if (Status/=0) call Error( "Error reading the argument Input%Ttra" )
  if (i_Debug_PL) call Logger%Write( "Input%Ttra_char = ", Input%Ttra_char, "; Input%Ttra = ", Input%Ttra  )

  call getarg( 2, Input%Tint_char )
  read(Input%Tint_char, "(d20.10)", iostat=Status) Input%Tint
  if (Status/=0) call Error( "Error reading the argument Input%Tint" )
  if (i_Debug_PL) call Logger%Write( "Input%Tint_char = ", Input%Tint_char, "; Input%Tint = ", Input%Tint)
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING THE INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_PL) call Logger%Write( "Reading the input", NewLine=.True. )
  if (i_Debug_PL) call Logger%Write( "-> Calling Input%Initialize" )
  call Input%Initialize( i_Debug=i_Debug_PL )
  Input%TaskType = 3
  if (i_Debug_PL) call Logger%Write( "-> Done Input%Initialize" )
! ==============================================================================================================


! ==============================================================================================================
!   Creating Folder for the System in the Output Folder
! ==============================================================================================================  
  if (i_Debug_PL) call Logger%Write( "Creating Folder ", adjustl(trim(Input%System)) )
  call system('mkdir -p ' // adjustl(trim(Input%OutputDir)) // '/' // adjustl(trim(Input%System)) )
! ============================================================================================================== 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 
!!! TO BE RESTORED
!!!
! ! ==============================================================================================================
! !  ALLOCATING THE ORIGINAL-LEVELS and THE CUT-LEVELS CONTAINERS AND THE BINNED MOLECULES
! ! ==============================================================================================================
!   ! Allocating the Molecules Energy Levels (containers for q.n.s etc) based on the Nb of Molecules
!   allocate(LevelsContainer_Orig(Input%NMolecules), Stat=Status )
!   if (Status/=0) call Error( "Error allocating LevelsContainer_Orig" )
!   if (i_Debug_PL) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Orig" )

!   ! Allocating the Molecules Energy Levels (containers for q.n.s etc) based on the Nb of Molecules
!   allocate(LevelsContainer_Cut(Input%NMolecules), Stat=Status )
!   if (Status/=0) call Error( "Error allocating LevelsContainer_Cut" )
!   if (i_Debug_PL) call Logger%Write( "Allocated ", Input%NMolecules, " LevelsContainer_Cut" )

!   ! Allocating the Binned-Molecules (containers for bins etc.) based on the Nb of Binned Molecules
!   allocate(BinnedMolecule(Input%NMolecules), Stat=Status )
!   if (Status/=0) call Error( "Error allocating BinnedMolecule" )
!   if (i_Debug_PL) call Logger%Write( "Allocated ", Input%NMolecules, " BinnedMolecule" )
! ! ==============================================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! ==============================================================================================================
!   INITIALIZING COLLISION OBJECT
! ==============================================================================================================
  if (i_Debug_PL) call Logger%Write( "Initializing the Collision object", NewLine=.True. )
  if (i_Debug_PL) call Logger%Write( "Calling Collision%Initialize" )
  call Collision%Initialize( Input, i_Debug=i_Debug_PL, i_Debug_Deep=i_Debug_PL_Deep )
  if (i_Debug_PL) call Logger%Write( "Done with Collision%Initialize" )
! ==============================================================================================================  



  call CPU_Time( EndTime )
  t_total   =   t_total + EndTime - StartTime
  if (i_Debug_PL) then
    call Logger%Write( "-> Total    : ", t_total, Fr="es15.8" )
  end if
! ==============================================================================================================


  if (i_Debug_PL) call Logger%Write( "Normal termination" )


End Program PreprocLevels
