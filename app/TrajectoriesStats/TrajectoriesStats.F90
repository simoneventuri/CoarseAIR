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
!   Computing Statistics from Trajectories for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! This Program computes cross sections from results of classical trajectories.
!
!   Mandatory Input Arguments: - Translational Temperature [K]
!                              - Internal Temperature [K]
!
! ==============================================================================================================
Program TrajectoriesStats

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Input_Class           ,only:  Input_Type
  use Statistics_Class      ,only:  Statistics_Type
  use Error_Class           ,only:  Error

  implicit none
  
  type(Input_Type)                        ::    Input
  type(Statistics_Type)                   ::    Statistics

  real(rkp)                               ::    Tran
  character(10)                           ::    Tran_char
  integer                                 ::    Status
  integer                                 ::    StatReadsBinary
  character(1)                            ::    StatReadsBinaryChar
  logical                                 ::    i_Debug_TS      = .True.
  logical                                 ::    i_Debug_TS_Deep = .True.
!  real(rkp)                               ::    StartTime, EndTime


!  call CPU_Time( StartTime )
  if (i_Debug_TS) call Logger%Initialize( "TrajectoriesStats.log"     , &                                        ! Opening the Log File using
                         Status          =       'REPLACE'            , &                                        ! replacing any previous log file
                         Position        =       'REWIND'             , &                                        ! rewinding to the top of the file
                         Procedure       =       'TrajectoriesStats'  , &                                        ! loading the calling procedure name
                         Indentation     =       2               )                                               ! and setting the initial indentation level


! ==============================================================================================================
!  READING THE PROGRAM ARGUMENTS
! ==============================================================================================================
  call getarg( 1, Tran_char )
  read(Tran_char, "(d20.10)", iostat=Status) Tran
  if (Status/=0) call Error( "Error reading the argument Tran_char" )
  if (i_Debug_TS) call Logger%Write( "Tran_char = ", Tran_char, "; Tran = ", Tran  )
  
  call getarg( 2, Input%Tint_char )
  read(Input%Tint_char, "(d20.10)", iostat=Status) Input%Tint
  if (Status/=0) call Error( "Error reading the argument Input%Tint" )
  if (i_Debug_TS) call Logger%Write( "Input%Tint_char = ", Input%Tint_char, "; Input%Tint = ", Input%Tint )

  call getarg( 3, StatReadsBinaryChar )
  read(StatReadsBinaryChar, "(I1)", iostat=Status) StatReadsBinary
  if (Status/=0) call Error( "Error reading the argument StatReadsBinary" )
  if (StatReadsBinary == 1) then
    Input%StatReadsBinaryFlg = .True.
    if (i_Debug_TS) call Logger%Write( "Input%StatReadsBinaryFlg = ", Input%StatReadsBinaryFlg, "; Reading Trajectories From Binary Files" )
  else
    if (i_Debug_TS) call Logger%Write( "Input%StatReadsBinaryFlg = ", Input%StatReadsBinaryFlg, "; Reading Trajectories From ASCI Files" )
  end if

  call getarg( 4, Input%PESoI_char )
  read( Input%PESoI_char,     '(I6)' ) Input%PESoI
  if (i_Debug_TS) call Logger%Write( "PES of Interest:     Input%PESoI = ", Input%PESoI)
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING THE INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_TS) call Logger%Write( "Reading the input", NewLine=.True. )
  if (i_Debug_TS) call Logger%Write( "-> Calling Input%Statistics" )
  call Input%Initialize( i_Debug=i_Debug_TS )
  Input%TaskType = 5
  if (i_Debug_TS) call Logger%Write( "-> Done initializing Input" )
  
  if ( trim(adjustl(Input%TtraModel)) .eq. "Boltzmann" ) then
    Input%Ttra      = Tran
    Input%Ttra_char = Tran_char
    if (i_Debug_TS) call Logger%Write( "Input%Ttra_char = ", Input%Ttra_char, "; Input%Ttra = ", Input%Ttra  )
  elseif ( trim(adjustl(Input%TtraModel)) .eq. "Uniform" ) then
    Input%Erel      = Tran
    Input%Erel_char = Tran_char
    if (i_Debug_TS) call Logger%Write( "Input%Erel_char = ", Input%Erel_char, "; Input%Erel = ", Input%Erel  )
  end if
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING THE STATISTICS OBJECT
! ==============================================================================================================
  if (i_Debug_TS) call Logger%Write( "Initializing the Statistics object", NewLine=.True. )
  if (i_Debug_TS) call Logger%Write( "-> Calling Statistics%Initialize" )
  call Statistics%Initialize( Input, i_Debug=i_Debug_TS, i_Debug_Deep=i_Debug_TS_Deep )
  if (i_Debug_TS) call Logger%Write( "-> Done initializing Statistics" )
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING THE STATISTICS OBJECT
! ==============================================================================================================
  if (i_Debug_TS) call Logger%Write( "Processing the data", NewLine=.True. )
  if (i_Debug_TS) call Logger%Write( "-> Calling Statistics%Process" )
  call Statistics%Process( i_Debug=i_Debug_TS, i_Debug_Deep=i_Debug_TS_Deep )
  if (i_Debug_TS) call Logger%Write( "-> Done processing the data" )
! ==============================================================================================================


!  call CPU_Time( EndTime )
!  call Logger%Write( "-> Total time: ", EndTime - StartTime, Fr="es15.8" )
  
  
  if (i_Debug_TS) call Logger%Write( "Normal termination" )


End Program TrajectoriesStats