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
!         Plotting the PES for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! 
!
! ==============================================================================================================
Program PlotPES_Program

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use Timing_Module
  
  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use PlotPES_Factory_Class ,only:  PlotPES_Factory_Type
  use PlotPES_Class         ,only:  PlotPES_Type
  
  implicit none

  type(Input_Type)                      ::    Input
  type(Collision_Type)                  ::    Collision

  type(PlotPES_Factory_Type)            ::    PlotPES_Factory
  class(PlotPES_Type)    ,allocatable   ::    PlotPES

  integer                               ::    Status, Unit
  character(150)                        ::    LevelOutputDir
  
  logical                ,parameter     ::    i_Debug_PP        = .True.
  logical                ,parameter     ::    i_Debug_PP_Deep   = .True.
  
  Input%iProc =  1
  Input%iNode =  1

  ! trim(adjustl(Input%OutputDir))
  if (i_Debug_PP) call Logger%Initialize( "PlotPES.log",           &                                      ! Opening the Log File using
                              Status          =       'REPLACE',   &                                      ! replacing any previous log file
                              Position        =       'REWIND',    &                                      ! rewinding to the top of the file
                              Procedure       =       'PlotPES',   &                                      ! loading the calling procedure name
                              Indentation     =       2   )                                               ! and setting the initial indentation level



!  call CPU_Time( StartTime )
  
! ==============================================================================================================
!   INITIALIZING INPUT
! ==============================================================================================================
  call getarg( 1, LevelOutputDir )
  if (i_Debug_PP) call Logger%Write( "LevelOutputDir = ", LevelOutputDir )
  
  if (i_Debug_PP) call Logger%Write( "Reading the input. Calling Input%Initialize", NewLine=.True. )
  call Input%Initialize( i_Debug=i_Debug_PP )
  if (i_Debug_PP) call Logger%Write( "Calling Input%PlotPES" )
  call Input%PlotPES( i_Debug=i_Debug_PP )
  Input%TaskType = 1
  if (i_Debug_PP) call Logger%Write( "Done initializing Input" )
  
  Input%LevelOutputDir = LevelOutputDir
  Input%TaskType = 1
! ==============================================================================================================

  call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/')


! ==============================================================================================================
!   INITIALIZING COLLISION
! =============================================================================================================
  if (i_Debug_PP) call Logger%Write( "Setting the Collision object; Calling Collision%Initialize", NewLine=.True. )
  call Collision%Initialize( Input, i_Debug=i_Debug_PP, i_Debug_Deep=i_Debug_PP_Deep )
  if (i_Debug_PP) call Logger%Write( "Done Initializing Collision" )
! ==============================================================================================================



! ==============================================================================================================
!   INITIALIZING PlotPES OBJECT
! =============================================================================================================
  if (i_Debug_PP) call Logger%Write( "Setting the Collision object; Calling Collision%Initialize", NewLine=.True. )
  call PlotPES_Factory%Define_PlotPEs( Input, PlotPES, i_Debug=i_Debug_PP )
  call PlotPES%Initialize( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
  if (i_Debug_PP) call Logger%Write( "Done Initializing Collision" )
! ==============================================================================================================



! ==============================================================================================================
!   COMPUTING AND PLOTTING PES
! =============================================================================================================
  if (trim(adjustl(Input%PESOrDiatFlg)) .eq. 'Diatomic') then

    if (i_Debug_PP) call Logger%Write( "Calling DiatPot", NewLine=.True. )
    call PlotPES%DiatPot( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
    if (i_Debug_PP) call Logger%Write( "Done with DiatPot", NewLine=.True. ) 
  
  elseif (trim(adjustl(Input%PESOrDiatFlg)) .eq. 'PES') then

    if (Input%PlotPES_ReadPntsFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling ReadPoints", NewLine=.True. )
      call PlotPES%ReadPoints( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with ReadPoints", NewLine=.True. )
      
  !    if (i_Debug_PP) call Logger%Write( "Calling EvaluatePoints", NewLine=.True. )
  !    call EvaluatePoints( Input, Collision, i_Debug=i_Debug_PP )
  !    if (i_Debug_PP) call Logger%Write( "Done with EvaluatePoints", NewLine=.True. )
      
  !    if (i_Debug_PP) call Logger%Write( "Calling ComputeCuts", NewLine=.True. )
  !    call ComputeCuts( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
  !    if (i_Debug_PP) call Logger%Write( "Done with ComputeCuts", NewLine=.True. )
    end if

    if (Input%PlotPES_GridFlg) then
      if (Input%StochPESFlg) then
        if (i_Debug_PP) call Logger%Write( "Calling GridForStochPES", NewLine=.True. )
        call PlotPES%GridForStochPES( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
        if (i_Debug_PP) call Logger%Write( "Done with GridForStochPES", NewLine=.True. )
      else
        if (i_Debug_PP) call Logger%Write( "Calling Grid", NewLine=.True. )
        call PlotPES%Grid( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
        if (i_Debug_PP) call Logger%Write( "Done with Grid", NewLine=.True. )
      end if
    end if
    
    if (Input%PlotPES_DoubleGridFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling DoubleGrid", NewLine=.True. )
      call PlotPES%DoubleGrid( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with DoubleGrid", NewLine=.True. )
    end if
    
    if (Input%PlotPES_TripleGridFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling TripleGrid", NewLine=.True. )
      call PlotPES%TripleGrid( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with TripleGrid", NewLine=.True. )
    end if
    
    if (Input%PlotPES_GridForScatterFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling GridForScatter", NewLine=.True. )
      call PlotPES%GridForScatter( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with GridForScatter", NewLine=.True. )
    end if
    
    if (Input%PlotPES_StatsFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling StochPESStats", NewLine=.True. )
      call PlotPES%StochPESStats( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with StochPESStats", NewLine=.True. )
    end if
    
    if (Input%PlotPES_VargasPaperFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling PlotsVargasPaper", NewLine=.True. )
      call PlotPES%PlotsVargasPaper( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with PlotsVargasPaper", NewLine=.True. )
    end if

    if (Input%PlotPES_Rot3rdFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling PlotPES_Rot3rd", NewLine=.True. )
      call PlotPES%Rot3rd( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with PlotPES_Rot3rd", NewLine=.True. )
    end if

    if (Input%PlotPES_IsoTriFlg) then
      if (i_Debug_PP) call Logger%Write( "Calling PlotPES_IsoTri", NewLine=.True. )
      call PlotPES%IsoTri( Input, Collision, Collision%NPairs, Collision%NAtoms, i_Debug=i_Debug_PP )
      if (i_Debug_PP) call Logger%Write( "Done with PlotPES_IsoTri", NewLine=.True. )
    end if

  end if ! PES/DIAT

  if (i_Debug_PP) call Logger%Write( "Done Plotting PES" )
! ==============================================================================================================]
  

  if (i_Debug_PP) call Logger%Write( "Normal termination" )


End Program
