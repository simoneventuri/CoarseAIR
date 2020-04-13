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
!         Running Trajectories for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! Driver program for classical trajectory program for a general number of atoms.
! The code uses a bulirsch-stoer ODE integrator.
!
! Mandatory Input Arguments: - Path to the LevelOutputDir, where the trajectories of the current (Level/Bin & Node & Processor) will be written 
!                            - Translational Temperature [K]
!                            - Internal Temperature [K]
!                            - Total Number of Processors
!                            - Processor Identifier for Parallelization, iProc (1 <= iProc <= NProc/NNodes)
!                            - Node      Identifier for Parallelization, iNode (1 <= iNode <= NNodes)
!                            - Current Inital Level/Bin
!
! ==============================================================================================================
Program RunTrajectories

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Integrator_Class      ,only:  Integrator_Type
  ! use LevelsContainer_Class ,only:  WriteTempList
  use Error_Class           ,only:  Error
  use Timing_Module

  implicit none

  type(Input_Type)                      ::    Input
  type(Collision_Type)                  ::    Collision
  type(Integrator_Type)                 ::    Integrator

  character(200)                        ::    LevelOutputDir
  integer                               ::    iMol  
  integer                               ::    iBin
  integer                               ::    iPES
  integer                               ::    TotDeg, NTrajDeg, NTrajPES
  integer                               ::    Status, Unit
  integer                               ::    temp1, temp2
  integer ,dimension(:) ,allocatable    ::    Levels_per_Bin
  character(:)  ,allocatable            ::    FileName
  real(rkp)                             ::    Tran
  character(10)                         ::    Tran_char
!  real(rkp)                             ::    StartTime, EndTime
  
  logical                ,parameter     ::    i_Debug_RT        = .True.
  logical                ,parameter     ::    i_Debug_RT_Deep   = .True. 
  logical                ,parameter     ::    i_Debug_RT_ODE    = .True. 
  
!  call CPU_Time( StartTime )
  
  if (i_Debug_RT) call Logger%Initialize( "RunTrajectories.log",           &                                      ! Opening the Log File using
                              Status          =       'REPLACE',           &                                      ! replacing any previous log file
                              Position        =       'REWIND',            &                                      ! rewinding to the top of the file
                              Procedure       =       'RunTrajectories',   &                                      ! loading the calling procedure name
                              Indentation     =       2           )                                               ! and setting the initial indentation level

! ==============================================================================================================
!  READING THE PROGRAM ARGUMENT
! ==============================================================================================================            
  call getarg( 1, LevelOutputDir )
  if (i_Debug_RT) call Logger%Write( "LevelOutputDir = ", LevelOutputDir )
 
  call getarg( 2, Tran_char )
  read(Tran_char, "(d20.10)", iostat=Status) Tran
  if (Status/=0) call Error( "Error reading the argument Tran_char" )
  if (i_Debug_RT) call Logger%Write( "Tran_char = ", Tran_char, "; Tran = ", Tran  )
  
  call getarg( 3, Input%Tint_char )
  read(Input%Tint_char, "(d20.10)", iostat=Status) Input%Tint
  if (Status/=0) call Error( "Error reading the argument Input%Tint" )
  if (i_Debug_RT) call Logger%Write( "Input%Tint_char = ", Input%Tint_char, "; Input%Tint = ", Input%Tint)
  
  call getarg( 4, Input%NProc_char )
  if (trim(adjustl(Input%NProc_char)) .eq. '') then
    Input%NProc_char  = '1'
    Input%NProc       =  1
  else 
    read( Input%NProc_char, '(I3)' ) Input%NProc
  end if
  if (i_Debug_RT) call Logger%Write( "Input%NProc_char = ", Input%NProc_char, "; Input%NProc = ", Input%NProc)
  
  call getarg( 5, Input%iProc_char )
  if (trim(adjustl(Input%iProc_char)) .eq. '') then
    Input%iProc_char  = '1'
    Input%iProc       =  1
  else 
    read( Input%iProc_char, '(I3)' ) Input%iProc
  end if
  if (i_Debug_RT) call Logger%Write( "Input%iProc_char = ", Input%iProc_char, "; Input%iProc = ", Input%iProc)
  
  call getarg( 6, Input%NNode_char )
  if (trim(adjustl(Input%NNode_char)) .eq. '') then
    Input%NNode_char  = '1'
    Input%NNode       =  1
  else 
    read( Input%NNode_char, '(I3)' ) Input%NNode
  end if
  if (i_Debug_RT) call Logger%Write( "Input%NNode_char = ", Input%NNode_char, "; Input%NNode = ", Input%NNode)
  
  call getarg( 7, Input%iNode_char )
  if (trim(adjustl(Input%iNode_char)) .eq. '') then
    Input%iNode_char  = '1'
    Input%iNode       =  1
  else 
    read( Input%iNode_char, '(I3)' ) Input%iNode
  end if
  if (i_Debug_RT) call Logger%Write( "Input%iNode_char = ", Input%iNode_char, "; Input%iNode = ", Input%iNode)
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING INPUT
! ==============================================================================================================
  if (i_Debug_RT) call Logger%Write( "Reading the input. Calling Input%Initialize", NewLine=.True. )
  call Input%Initialize( i_Debug=i_Debug_RT )
  Input%TaskType = 4
  if (i_Debug_RT) call Logger%Write( "Done initializing Input" )
  
  Input%LevelOutputDir =trim(adjustl(LevelOutputDir))
  if (i_Debug_RT) call Logger%Write( "This%LevelOutputDir = ", Input%LevelOutputDir )
  
  if ( trim(adjustl(Input%TtraModel)) .eq. "Boltzmann" ) then
    Input%Ttra      = Tran
    Input%Ttra_char = Tran_char
    if (i_Debug_RT) call Logger%Write( "Input%Ttra_char = ", Input%Ttra_char, "; Input%Ttra = ", Input%Ttra  )
  elseif ( trim(adjustl(Input%TtraModel)) .eq. "Uniform" ) then
    Input%Erel      = Tran
    Input%Erel_char = Tran_char
    if (i_Debug_RT) call Logger%Write( "Input%Erel_char = ", Input%Erel_char, "; Input%Erel = ", Input%Erel  )
  end if
! ==============================================================================================================


! ==============================================================================================================
!  READING THE REMAINING PROGRAM ARGUMENTS
! ==============================================================================================================
  do iMol = 1,Input%NInitMolecules
  
    call getarg( 7 + iMol, Input%BinOI_char(iMol) )
    read( Input%BinOI_char(iMol),     '(I6)' ) Input%BinOI(iMol)
    if (i_Debug_RT) call Logger%Write( "Level/Bin of Interest for Molecule Nb", iMol, ":     Input%BinOI(iMol) = ", Input%BinOI(iMol))
  
  end do
! ==============================================================================================================


! ==============================================================================================================
!  DEFINING NB of TRAJECTORIES
! ==============================================================================================================
  if  ( ( (Input%ProportionalAllocation == "yes") .or. (Input%ProportionalAllocation == "YES") ) ) then
    if (i_Debug_RT) call Logger%Write( "Proportional Allocation has been selected for sampling the energy levels in the bins.")
    
    iMol = 1
  
    allocate(Levels_per_Bin(Input%NBins(iMol)), Stat=Status )
    if (Status/=0) call Error( "Error allocating Levels_per_Bin" )
    Levels_per_Bin = 0
  
    FileName = trim(adjustl(Input%OutputDir)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
               trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qnsFirst.dat'
    open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
    
      if (Status==0) then
        read(Unit,*,iostat=Status)
        do iBin = 1,Input%NBins(iMol)
          read(Unit,*,iostat=Status) temp1, temp2, Levels_per_Bin(iBin)
          if (Status /= 0) exit
        end do
        
    close(Unit)
        if (i_Debug_RT) call Logger%Write( "Total Nb of Energy Levels,                    sum(Levels_per_Bin)              = ", sum(Levels_per_Bin) )
        if (i_Debug_RT) call Logger%Write( "Nb of Energy Levels in the Bin of Interest,   Levels_per_Bin(Input%BinOI(iMol)) = ", Levels_per_Bin(Input%BinOI(iMol)) )
        
        Input%NTraj =  ceiling( real(Input%NTraj) * real(Levels_per_Bin(Input%BinOI(iMol))) / real(sum(Levels_per_Bin)) )
        if (i_Debug_RT) call Logger%Write( "Nb of Trajectories to run for Input Bin,  Input%NTraj = ", Input%NTraj )
      else
        call Error( "Error reading " // trim(adjustl(Input%OutputDir)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
               trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qnsFirst.dat' )
      end if
    deallocate(Levels_per_Bin)

  end if
  
  if (i_Debug_RT) call Logger%Write( "Nb of Processors,  Input%NProc = ", Input%NProc )
  Input%NTraj =  max(ceiling(real(Input%NTraj) / real(Input%NProc,rkp)), 1)

  


  if (i_Debug_RT) call Logger%Write( "Nb of Trajectories to run for Input Processor,  Input%NTraj = ", Input%NTraj )
  
  if ((Input%NTrajBatch > 0) .and. (Input%NTrajBatch < Input%NTraj)) then
    if (i_Debug_RT) call Logger%Write( "Selected to run Trajectories in Mini-Batch Mode; Nb of Trajectories per Mini-Batch:   Input%NTrajBatch = ", Input%NTrajBatch )
    Input%NBatch = ceiling( real(Input%NTraj) / real(Input%NTrajBatch) )
    if (i_Debug_RT) call Logger%Write( "Nb Mini-Batches:   Input%NBatch = ", Input%NBatch )
    Input%NTraj = int(Input%NTrajBatch * Input%NBatch)
    if (i_Debug_RT) call Logger%Write( "Final Nb of Trajectories to run for Input Processor,  Input%NTraj = ", Input%NTraj )
  else
    Input%NBatch     = 1
    Input%NTrajBatch = Input%NTraj
    if (i_Debug_RT) call Logger%Write( "Nb of Trajectories per Mini-Batch set to:   Input%NTrajBatch = ", Input%NTrajBatch )
  end if


  TotDeg   = sum(Input%PES_DegeneracyINT)
  NTrajDeg = ceiling ( real(Input%NTrajBatch, rkp) / real(TotDeg, rkp) )
  NTrajPES = 0
  do iPES=1,Input%NPESs
    NTrajPES = NTrajPES + ceiling ( Input%PES_DegeneracyINT(iPES) * real(NTrajDeg, rkp) )
  end do
  if (NTrajPES /= Input%NTrajBatch) then
    if (i_Debug_RT) call Logger%Write( "Adding ", (NTrajPES - Input%NTrajBatch), " Trajectories to Input%NTrajBatch in order to respect the Statistics on the PES Degeneracies" )
    Input%NTrajBatch = NTrajPES
  end if
! ==============================================================================================================  
  

! ! ==============================================================================================================
! !   COPING ENERGY LEVELS LIST
! ! ==============================================================================================================
!   if (i_Debug_RT) call Logger%Write( "-> Scanning Levels File" )
!   call WriteTempList( Input, i_Debug=i_Debug_RT ) 
!   if (i_Debug_RT) call Logger%Write( "-> Done Scanning Levels File" )
! ! ==============================================================================================================
  
  
! ==============================================================================================================
!   INITIALIZING COLLISION
! =============================================================================================================
  if (i_Debug_RT) call Logger%Write( "Setting the Collision object", NewLine=.True. )
  if (i_Debug_RT) call Logger%Write( "-> Calling Collision%Initialize" )
  call Collision%Initialize( Input, i_Debug=i_Debug_RT, i_Debug_Deep=i_Debug_RT_Deep )
  if (i_Debug_RT) call Logger%Write( "-> Done Initializing Collision" )
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING INTEGRATOR
! =============================================================================================================
  if (i_Debug_RT) call Logger%Write( "Setting the Integrator object", NewLine=.True. )
  if (i_Debug_RT) call Logger%Write( "-> Calling Integrator%Initialize" )
  call Integrator%Initialize( Input, i_Debug=i_Debug_RT )
  if (i_Debug_RT) call Logger%Write( "-> Done Initializing Integrator" )
! ==============================================================================================================


! ==============================================================================================================
!   INTEGRATE TRAJECTORIES
! ==============================================================================================================
  if (i_Debug_RT) call Logger%Write( "Integrating trajectories", NewLine=.true. )
  if (i_Debug_RT) call Logger%Write( "-> Calling Integrator%Integrate" )
  call Integrator%Integrate( Input, Collision, i_Debug=i_Debug_RT, i_Debug_ODE=i_Debug_RT_ODE )
  if (i_Debug_RT) call Logger%Write( "Done Integrating trajectories" )
! ==============================================================================================================


!! ==============================================================================================================
!!   WRITING TIMING INFO
!! ==============================================================================================================
!  ti1 = Zero   ; ic1 = 0
!  ti2 = Zero   ; ic2 = 0
!  ti3 = t_vphase   ; ic3 = i_vphase
!  ti4 = Zero   ; ic4 = 0
!  ti5 = t_inicon   ; ic5 = i_inicon
!  ti6 = t_potcl    ; ic6 = i_potcl
!  ti7 = t_pot      ; ic7 = i_pot
!  ti8 = Zero   ; ic8 = 0
!  ti9 = Zero   ; ic9 = 0

!  call CPU_Time( EndTime )
!  t_total   =   t_total + EndTime - StartTime
!  if (i_Debug_RT) then
!    call Logger%Write( "               Time [s]          Count" )
!    call Logger%Write( "-> ham      : ", ti1, ic1, Fr="es15.8" )
!    call Logger%Write( "-> bsstepfv : ", ti2, ic2, Fr="es15.8" )
!    call Logger%Write( "-> vphase   : ", ti3, ic3, Fr="es15.8" )
!    call Logger%Write( "-> der      : ", ti4, ic4, Fr="es15.8" )
!    call Logger%Write( "-> inicon   : ", ti5, ic5, Fr="es15.8" )
!    call Logger%Write( "-> potcl    : ", ti6, ic6, Fr="es15.8" )
!    call Logger%Write( "-> pot      : ", ti7, ic7, Fr="es15.8" )
!    call Logger%Write( "-> mmidv    : ", ti8, ic8, Fr="es15.8" )
!    call Logger%Write( "-> rzextrv  : ", ti9, ic9, Fr="es15.8" )
!    call Logger%Write( "-> SubTotal : ", ti1 + ti2 + ti3 + ti4 + ti5 + ti6 + ti7 + ti8 + ti9, Fr="es15.8" )
!    call Logger%Write( "-> Total    : ", t_total, Fr="es15.8" )
!  end if
!! ==============================================================================================================

  if (i_Debug_RT) call Logger%Write( "Normal termination" )


End Program
