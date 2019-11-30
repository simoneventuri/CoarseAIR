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

Module Trajectory_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Trajectory_Type

  Type      ::    Trajectory_Type
    integer                               ::    Idx
    integer                               ::    NEqtTot
    integer                               ::    NEqtVar

!   Constant parameters common to all trajectories
    integer                               ::    icoord                    ! Indicator of the coordinate system: jacobi coordinates => True, cartesiancoordinates => False. Old name: icoord=0 => cartesian coordinates, icoord=1 => jacobi coordinates

!   Constant parameter characterizing the trajectory
    real(rkp)                             ::    b

!   Variable different for each traj and characterizing the initial state: Set in setis
    integer                               ::    jqn                       ! Rotational quantum number of the inital state
    integer                               ::    vqn                       ! Vibrational quantum number of the inital state
    integer                               ::    iState                    ! Index of the selected initial state
    integer                               ::    itype                     ! Arrangement quantum number

    real(rkp)                             ::    rv                        ! Vibrational quantum number
    real(rkp)                             ::    rj                        ! Rotational quantum number
    real(rkp)                             ::    RN_VibPhase               ! Vibrational phase (random number)

!   Time-dependent variables
    real(rkp) ,dimension(:) ,allocatable  ::    PaQ                       ! Dim=(NEqtTot)
    real(rkp)                             ::    H                         ! Hamiltonian. Dim=(NEqtTot)
    real(rkp) ,dimension(:) ,allocatable  ::    erra                      ! (2*maxcd)
    real(rkp) ,dimension(:) ,allocatable  ::    Rpi
    real(rkp)                             ::    dt                        ! Time step
    real(rkp)                             ::    t                         ! Time
    real(rkp)                             ::    smax                      !
    real(rkp)                             ::    smin                      !
    integer                               ::    irej                      !
    integer                               ::    iste                      !
    integer                               ::    irec                      !

!   Input-related variables
    real(rkp)                               ::    tmax                    ! Maximum trajectory time above which the trajectory is stopped:      if t     > tmax     => Stop
    real(rkp)                               ::    RpiMax                  ! Inverse of the maximum distance between atom pair [1/bohr]: if Rpi <= RpiMax then the trajectories are stopped
    logical                                 ::    BwrdIntegrationFlg
    
    integer                                 ::    RestartUnit = 9
    character(:)  ,allocatable              ::    RestartFileName

  contains
    private
    procedure ,public   ::    Initialize              =>  Initialize_Trajectory
    procedure ,public   ::    IsConverged
  End Type

  logical   ,parameter  ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Trajectory( This, NEqtTot, NPairs, i_Debug )

  use Parameters_Module     ,only:  Zero

  class(Trajectory_Type)                    ,intent(out)    ::    This
  integer                                   ,intent(in)     ::    NEqtTot           ! Total number of equations (for both Coodinates/Momenta): 2 * NEqtVar
  integer                                   ,intent(in)     ::    NPairs            ! Number of connections between atoms (Old name was 'maxpr'): NAtoms*(NAtoms-1)/2
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Trajectory" )
  !i_Debug_Loc   =     Logger%On()
  

  This%NEqtTot    =     NEqtTot
  This%NEqtVar    =     This%NEqtTot / 2
  allocate( This%PaQ (NEqtTot) ); This%PaQ  =   Zero
  allocate( This%erra(NEqtTot) ); This%erra =   Zero
  allocate( This%Rpi (NPairs)  ); This%Rpi  =   Zero
  This%H     =   Zero
  This%dt    =   Zero
  This%t     =   Zero
  This%smax  =   Zero
  This%smin  =   Zero
  This%irej  =   0
  This%iste  =   0
  This%irec  =   0

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Function IsConverged( This ) result( Converged )
! This procedure assesses if we can stop the integration process.
! There are two criteria for assessing the convergence:
! * if at least one atom has exeeeded the maximum interatomic distance or,
! * if the trajectory time has exceeded the maximum time
! Note that via 'Hamiltonian', potcl returns the inverses of the atom-atom distances.
! For a given trajectory, the Converged variable is set to 0 if the traj is not converged and to -1 if it is
! the variable used to be called 'iDone'

  class(Trajectory_Type)                    ,intent(in)     ::    This
  
  logical                                                   ::    Converged
  
  Converged         =   .False.                   ! By default, set all trajectories as NOT converged
  if ( This%BwrdIntegrationFlg ) then    ! Case of backward integration
    Converged     =   abs(This%t) >= 1.0E-4_rkp
  else                                    ! Case of forward integration
    Converged     =   ( This%t >= This%tmax ) .or. any( This%Rpi(:) < This%RpiMax )
  end if
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
