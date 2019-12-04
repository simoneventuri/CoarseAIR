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
!                  Initializing CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! The program reads QCT.inp and writes some of the input in a new file, in order to pass such variables to the 
! bash files that are used for executing the code.
!
! ==============================================================================================================
Program Initialize

  use Input_Class           ,only:  Input_Type
  use Logger_Class          ,only:  Logger
  
  implicit none

  type(Input_Type)                      ::    Input

  logical                               ::    i_Debug_In = .False.
  
  if (i_Debug_In)  call Logger%Initialize( "Initialize.log",                      &             ! Opening the Log File using
                                            Status          =       'REPLACE',    &             ! replacing any previous log file
                                            Position        =       'REWIND',     &             ! rewinding to the top of the file
                                            Procedure       =       'Initialize', &             ! loading the calling procedure name
                                            Indentation     =       2           )               ! and setting the initial indentation level

! ==============================================================================================================
!   INITIALIZING INPUT
! ==============================================================================================================
  if (i_Debug_In) call Logger%Write( "Reading the input. Calling Input%Initialize", NewLine=.True. )
  call Input%Initialize( i_Debug=i_Debug_In )
  if (i_Debug_In) call Logger%Write( "Done initializing Input" )
! ==============================================================================================================


! ==============================================================================================================
!   WRITING INPUT FOR BASH FILES
! ==============================================================================================================
  if (i_Debug_In) call Logger%Write( "Reading the input. Calling Input%Initialize", NewLine=.True. )
  call Input%WriteBashInputVariables( )
  if (i_Debug_In) call Logger%Write( "Done initializing Input" )
! ==============================================================================================================


  if (i_Debug_In) call Logger%Write( "Normal termination" )
  

End Program