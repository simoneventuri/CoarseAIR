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
!         Postprocessing Trajectories for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! This Program ...
!
!   Mandatory Input Arguments: - Translational Temperature [K]
!                              - Internal Temperature [K]
!                              - Nb of Trajectories, for checking convergence
!                              - Current Initial Levels/Bins
!
! ==============================================================================================================
Program PostTrajectories

  use Parameters_Module        ,only:  rkp, Zero, Half, One, Two, Six
  use Logger_Class             ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class              ,only:  Error
  use Input_Class              ,only:  Input_Type
  use Collision_Class          ,only:  Collision_Type
  use Processes_Factory_Class  ,only:  Processes_Factory_Type
  use Processes_Class          ,only:  Processes_Type
  use Timing_Module

  implicit none
  
  Type(Input_Type)                                      :: Input
  type(Collision_Type)                                  :: Collision
  type(Processes_Factory_Type)                          :: Processes_Factory
  class(Processes_Type)                  ,allocatable   :: FinProcesses

  real(rkp)                                             :: Tran
  integer                                               :: TranInt
  character(10)                                         :: Tran_char
  character(10)                                         :: TranInt_char
  character(10)                                         :: Velocity_char
  real(rkp)                                             :: Velocity
  integer                                               :: iMol
  integer                                               :: NTraj_temp
  real(rkp)                                             :: TempVal
  character(10)                                         :: TempVal_char
  integer                                               :: Status
  real(rkp)                                             :: StartTime, EndTime

  logical                                  ,parameter   :: i_Debug_PT        = .False.
  logical                                  ,parameter   :: i_Debug_PT_Medium = .False.
  logical                                  ,parameter   :: i_Debug_PT_Deep   = .False.

  
  call CPU_Time( StartTime )
  
  if (i_Debug_PT) call Logger%Initialize(             "PostTrajectories.log",   &                                       ! Opening the Log File using
                              Status          =       'REPLACE',                &                                       ! replacing any previous log file
                              Position        =       'REWIND',                 &                                       ! rewinding to the top of the file
                              Procedure       =       'PostTrajectories',       &                                       ! loading the calling procedure name
                              Indentation     =       2           )                                                     ! and setting the initial indentation level

! ==============================================================================================================
!  READING THE PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 1, Tran_char )
  read(Tran_char, "(d20.10)", iostat=Status) Tran
  TranInt = int(Tran)
  write(TranInt_char, "(I10)", iostat=Status) TranInt
  if (Status/=0) call Error( "Error writing the argument TranInt_char" )
  if (i_Debug_PT) call Logger%Write( "TranInt_char = ", TranInt_char, "; Tran = ", Tran  )

  call getarg( 2, TempVal_char )
  read(TempVal_char, "(d20.10)", iostat=Status) Input%TintTemp
  if (Status/=0) call Error( "Error reading the argument Input%TintTemp" )
  if (i_Debug_PT) call Logger%Write( "TempVal_char = ", TempVal_char, "; Input%TintTemp = ", Input%TintTemp )

  call getarg( 3, Input%NTraj_char )
  read( Input%NTraj_char, '(I10)' ) NTraj_temp
  if (i_Debug_PT) call Logger%Write( "Input%NTraj_char = ", Input%NTraj_char, "; NTraj_temp = ", NTraj_temp )

  ! ==============================================================================================================
  !   READING VELOCITY CONVERSION FACTOR for CROSS SECTIONS RATES
  ! ==============================================================================================================
  call getarg( 4, Velocity_char )
  read( Velocity_char, "(d20.10)" ) Velocity
  if (i_Debug_PT) call Logger%Write( "Velocity_char = ", Velocity_char, "; Velocity for Conversion Factor = ", Velocity )
  ! ==============================================================================================================

! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Calling Input%Initialize" )
  call Input%Initialize( i_Debug=i_Debug_PT )
  Input%TaskType = 6

  Input%NTraj = NTraj_temp
  if (i_Debug_PT) call Logger%Write( "Total Nb of Converged Trajectories:     Input%NTraj = ", Input%NTraj )
  
  if ( trim(adjustl(Input%TtraModel)) .eq. "Boltzmann" ) then
    Input%Ttra      = Tran
    Input%Ttra_char = Tran_char
    if (i_Debug_PT) call Logger%Write( "Input%Ttra_char = ", Input%Ttra_char, "; Input%Ttra = ", Input%Ttra  )
  elseif ( trim(adjustl(Input%TtraModel)) .eq. "Uniform" ) then
    Input%Erel      = Tran
    Input%Erel_char = Tran_char
    if (i_Debug_PT) call Logger%Write( "Input%Erel_char = ", Input%Erel_char, "; Input%Erel = ", Input%Erel )
  end if
! ==============================================================================================================


! ==============================================================================================================
!  READING THE REMAINING PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 5, Input%PESoI_char )
  read( Input%PESoI_char,     '(I6)' ) Input%PESoI
  if (i_Debug_PT) call Logger%Write( "PES of Interest:     Input%PESoI = ", Input%PESoI)
  
  do iMol = 1,Input%NInitMolecules
  
    call getarg( 5 + iMol, Input%BinOI_char(iMol) )
    read( Input%BinOI_char(iMol),     '(I6)' ) Input%BinOI(iMol)
    if (i_Debug_PT) call Logger%Write( "Level/Bin of Interest for Molecule Nb", iMol, ":     Input%BinOI(iMol) = ", Input%BinOI(iMol))
    
  end do
! ==============================================================================================================
  

! ==============================================================================================================
!   INITIALIZING COLLISION OBJECT
! ==============================================================================================================
  if (i_Debug_PT) call Logger%Write( "Initializing the Collision object", NewLine=.True. )
  if (i_Debug_PT) call Logger%Write( "Calling Collision%Initialize" )
  call Collision%Initialize( Input, i_Debug=i_Debug_PT, i_Debug_Deep=i_Debug_PT_Deep )
  if (i_Debug_PT) call Logger%Write( "Done with Collision%Initialize" )
! ==============================================================================================================
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! TO-BE-RESTORED
!!!
!!! NOT DIFFERENTATING BETWEEN level_cut and level_original
!!! ! ==============================================================================================================
!!! !   COMPUTING PARTITION FUNCTIONS AND PARTITION FUNCTIONS RATIOS
!!! ! ==============================================================================================================
!!! FileName = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // trim(adjustl(Input%LevelsFileName(iMol)))
!!! FileName  = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/levels_cut.inp'
!!! if (i_Debug_PT) call Logger%Write( "Computing Partition Functions and related Ratios for Molecule Nb", iMol )
!!! call Compute_PartitionRatios( Input, iMol, LevelsContainer_Orig(iMol), LevelsContainer_Cut(iMol), i_Debug=i_Debug_PT ) (in LevelsConatainer_Class.F90)
!!! ! ==============================================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

! ==============================================================================================================
!   CONVERTING CROSS SECTIONS INTO RATES
! ==============================================================================================================
  call Processes_Factory%Define_Processes( Input, Collision, FinProcesses, i_Debug=i_Debug_PT )
! ==============================================================================================================


! ==============================================================================================================
!   CONVERTING CROSS SECTIONS INTO RATES
! ==============================================================================================================
  call FinProcesses%Convert_CrossSect_To_Rates( Input, Collision, [Collision%EqVelocity], i_Debug=i_Debug_PT, i_Debug_Deep=i_Debug_PT_Medium )
! ==============================================================================================================


  !
  call CPU_Time( EndTime )
  t_total   =   t_total + EndTime - StartTime
  if (i_Debug_PT) then
    call Logger%Write( "Total    : ", t_total, Fr="es15.8" )
  end if
  

  if (i_Debug_PT) call Logger%Write( "Normal termination" )

End Program PostTrajectories
