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

Module Trajectories_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Trajectories_Type

  Type      ::    Trajectories_Type
    integer                                 ::    NTraj     =   0   !< Number of trajectories. Can be different that the dimension of the components
    integer                                 ::    NTrajInit =   0   !< Number of initialized trajectories so far. This value is used to set the index of trajectories Idx
    integer   ,dimension(:)   ,allocatable  ::    Idx               !< Global index of trajectories
    real(rkp) ,dimension(:)   ,allocatable  ::    b                 !< Impact parameters. Dim=(NTrajMax)
    real(rkp) ,dimension(:)   ,allocatable  ::    t                 !< Physical time. Dim=(NTrajMax)
    integer   ,dimension(:)   ,allocatable  ::    iPES              !< Index mapping Trajectory to PES to use. Dim=(NTrajMax)
    real(rkp) ,dimension(:)   ,allocatable  ::    H                 !< Hamiltonian at t=Time. Dim=(NTrajMax)
    real(rkp) ,dimension(:,:) ,allocatable  ::    PaQ               !< Vector of solution variables at t=Time. Dim=(NEqtTot,NTrajMax)
    real(rkp) ,dimension(:)   ,allocatable  ::    H0                !< Hamiltonian at initial time t=0. Dim=(NTrajMax)
    real(rkp) ,dimension(:,:) ,allocatable  ::    PaQ0              !< Vector of solution variables at initial time t=0. Dim=(NEqtTot,NTrajMax)
    real(rkp) ,dimension(:)   ,allocatable  ::    dt                !< Integration time step. Dim=(NTrajMax)
    real(rkp) ,dimension(:,:) ,allocatable  ::    erra              !< Integration tolerence. Dim=(NEqtTot,NTrajMax)
    real(rkp) ,dimension(:,:) ,allocatable  ::    Rpi               !< Inverse of the distances of atom-atom pairs [1/bohr]. Dim=(NPairs,NTrajMax)
    integer   ,dimension(:)   ,allocatable  ::    ib        
            
    integer   ,dimension(:)   ,allocatable  ::    iTrajPES          !< Index for the current Nb of Trajectories in PES iPES-th 
!   Integration related variables
    real(rkp) ,dimension(:)   ,allocatable  ::    Smax              !<
    real(rkp) ,dimension(:)   ,allocatable  ::    Smin              !<
    integer   ,dimension(:)   ,allocatable  ::    irej              !<
    integer   ,dimension(:)   ,allocatable  ::    iste              !<
!   Input-related variables
    real(rkp)                               ::    tmax              ! Maximum trajectory time above which the trajectory is stopped:      if t     > tmax     => Stop
    real(rkp)                               ::    RpiMax            ! Inverse of the maximum distance distances of atom-atom pairs [1/bohr]: if Rpi <= RpiMax then the trajectories are stopped
    logical                                 ::    BwrdIntegrationFlg
    integer                                 ::    RestartUnit = 9
    character(:)  ,allocatable              ::    RestartFileName
  contains
    private
    procedure ,public   ::    Initialize  => InitializeTrajectories
    procedure ,public   ::    ReSet       => ReSetTrajectories
    procedure ,public   ::    AssessConvergence
    procedure ,public   ::    Reorder
  End Type

  logical   ,parameter  ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeTrajectories( This, NTraj, NEqtTot, NPairs, NPESs )

  use Parameters_Module           ,only:  Zero
    
  class(Trajectories_Type)                  ,intent(out)    ::    This
  integer                                   ,intent(in)     ::    NTraj             ! Number of trajectories
  integer                                   ,intent(in)     ::    NEqtTot           ! Total number of equations (for both Coodinates/Momenta): 2 * NEqtVar
  integer                                   ,intent(in)     ::    NPairs            ! Number of atom-atom pairs (Old name was 'maxpr'): NAtoms*(NAtoms-1)/2
  integer                                   ,intent(in)     ::    NPESs             ! Number of PESs
  
  integer                                                   ::    i                 ! Index of trajectories
  
  This%NTraj      =   NTraj                                                         ! Setting the number of traj. to the input value
  This%NTrajInit  =   0                                                             ! Initializing to zero the number of initialized trajectories so far
  allocate( This%Idx  (        NTraj) ); This%Idx   =   [(i,i=1,NTraj)]
  allocate( This%b    (        NTraj) ); This%b     =   Zero
  allocate( This%t    (        NTraj) ); This%t     =   Zero
  allocate( This%iPES (        NTraj) ); This%iPES  =   1
  allocate( This%H    (        NTraj) ); This%H     =   Zero
  allocate( This%PaQ  (NEqtTot,NTraj) ); This%PaQ   =   Zero
  allocate( This%H0   (        NTraj) ); This%H0    =   Zero
  allocate( This%PaQ0 (NEqtTot,NTraj) ); This%PaQ0  =   Zero
  allocate( This%dt   (        NTraj) ); This%dt    =   Zero
  allocate( This%erra (NEqtTot,NTraj) ); This%erra  =   Zero
  allocate( This%Rpi  (NPairs, NTraj) ); This%Rpi   =   Zero
  allocate( This%smax (        NTraj) ); This%smax  =   Zero
  allocate( This%smin (        NTraj) ); This%smin  =   Zero
  allocate( This%irej (        NTraj) ); This%irej  =   0
  allocate( This%iste (        NTraj) ); This%iste  =   0
  allocate( This%ib   (        NTraj) ); This%ib    =   0
  allocate( This%iTrajPES (    NPESs) ); This%iTrajPES = 0
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReSetTrajectories( This, NTraj )

  use Parameters_Module           ,only:  Zero
    
  class(Trajectories_Type)                  ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    NTraj             ! Number of trajectories
  
  integer                                                   ::    i                 ! Index of trajectories
  
  This%NTraj      =   NTraj                                                         ! Setting the number of traj. to the input value
  This%NTrajInit  =   0                                                             ! Initializing to zero the number of initialized trajectories so far
  This%Idx        =   [(i,i=1,NTraj)]
  This%b          =   Zero
  This%ib         =   0
  This%t          =   Zero
  This%iPES       =   1
  This%H          =   Zero
  This%PaQ        =   Zero
  This%H0         =   Zero
  This%PaQ0       =   Zero
  This%dt         =   Zero
  This%ib         =   0
  This%iTrajPES   =   0
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Reorder( This, Converged, NTrajBatch_Remaining )
!               i      iTi                     iTf     i+1
!   -----------------------------------------------------------------
!   |       |*******|       |       |       |       |*******|       |                |*******|
!   |   1   |***2***|   3   |   4   |   5   |   6   |***7***|   8   |                |*******|  =>  Converged point
!   |       |*******|       |       |       |       |*******|       |                |*******|
!   -----------------------------------------------------------------
!
!              iTi                     iTf
!   -------------------------------------------------
!   |       |       |       |       |       |       |
!   |   1   |   3   |   4   |   5   |   6   |   8   |
!   |       |       |       |       |       |       |
!   -------------------------------------------------

  class(Trajectories_Type)                  ,intent(inout)  ::    This
  logical   ,dimension(:)                   ,intent(inout)  ::    Converged
  integer                                   ,intent(in)     ::    NTrajBatch_Remaining

  integer                                                   ::    i                   ! Index of converged trajectories
  integer                                                   ::    iT                  ! Index of trajectories
  integer                                                   ::    NCvg                ! Number of converged trajectories
  integer                                                   ::    iIniOld, iFinOld    ! Old intial and final trajectories index
  integer                                                   ::    iIniNew, iFinNew    ! New intial and final trajectories index
  integer   ,dimension( count(Converged(1:This%NTraj)) )    ::    IdxCvg              ! Vector containing the index of the converged trajectories

!   call Logger%Write( "-> This%NTraj = ", This%NTraj )
!   call Logger%Write( "-> count(Converged(1:This%NTraj)) = ", count(Converged(1:This%NTraj)) )

  i               =   0
  do iT = 1,This%NTraj
    if ( Converged(iT) ) then
      i           =   i + 1
      IdxCvg(i)   =   iT
    end if
    Converged(iT) =   .False.
  end do
  NCvg            =   i         ! Number of converged trajectories. Value is equal to size(IdxCvg)


!   call Logger%Write( "-> size(IdxCvg) = ", size(IdxCvg) )
!   call Logger%Write( "-> NCvg = ", NCvg )
!   call Logger%Write( "-> IdxCvg = ", IdxCvg )

  if ( NCvg == 1 ) then       ! Dealing with the frequent cases: only one trajectory is converged
    if ( IdxCvg(1) == This%NTraj ) return
    iIniOld       =   IdxCvg(1) + 1
    iFinOld       =   NTrajBatch_Remaining + 1
    iIniNew       =   iIniOld - 1
    iFinNew       =   iFinOld - 1
    call Move( This, iIniNew, iFinNew, iIniOld, iFinOld )
  else
    i = 0
    Loop: do
      i = i + 1
      !call Logger%Write( "-> i = ", i )
      
      if ( i > NCvg ) exit Loop                                                         ! If current index if equal or larger than the number of converged point, then exit the loop
      iIniOld     =   IdxCvg(i) + 1                                                     ! Setting the index of the initial old point just after the current converged point
      if ( i+1 <= NCvg ) then                                                           ! If there is at lease one other converged trajectory point, then ...
        iFinOld   =   IdxCvg(i+1) - 1                                                   ! ... setting the index of the final old point just before the next converged point
      else                                                                              ! If there is no other converged trajectory point, then ...
        if ( IdxCvg(i) == This%NTraj ) exit Loop                                        ! ... if the current point is the last trajectory point, then no re-ordering is required here 
        iFinOld   =   This%NTraj                                                        ! ... and so exiting the loop setting the index of the final old point to the last trajectory point
      end if
      
      iIniNew     =   iIniOld - i
      iFinNew     =   iFinOld - i
!       write(Logger%Unit,"(4x,'[Reorder]: Reordering trajectories from [',g0,':',g0,'] to [',g0,':',g0,']')") iIniOld,iFinOld,   iIniNew,iFinNew
      call Move( This, iIniNew, iFinNew, iIniOld, iFinOld )
    end do Loop
  end if
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine Move( Traj, iIniNew, iFinNew, iIniOld, iFinOld )
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  integer                                   ,intent(in)     ::    iIniOld, iFinOld    ! Old intial and final trajectories index
  integer                                   ,intent(in)     ::    iIniNew, iFinNew    ! New intial and final trajectories index
  Traj%Idx(   iIniNew:iFinNew)    =   Traj%Idx(   iIniOld:iFinOld)
  Traj%t(     iIniNew:iFinNew)    =   Traj%t(     iIniOld:iFinOld)
  Traj%dt(    iIniNew:iFinNew)    =   Traj%dt(    iIniOld:iFinOld)
  Traj%b(     iIniNew:iFinNew)    =   Traj%b(     iIniOld:iFinOld)
  Traj%ib(    iIniNew:iFinNew)    =   Traj%ib(    iIniOld:iFinOld)
  Traj%H(     iIniNew:iFinNew)    =   Traj%H(     iIniOld:iFinOld)
  Traj%H0(    iIniNew:iFinNew)    =   Traj%H0(    iIniOld:iFinOld)
  Traj%PaQ( :,iIniNew:iFinNew)    =   Traj%PaQ( :,iIniOld:iFinOld)
  Traj%PaQ0(:,iIniNew:iFinNew)    =   Traj%PaQ0(:,iIniOld:iFinOld)
  Traj%erra(:,iIniNew:iFinNew)    =   Traj%erra(:,iIniOld:iFinOld)
  Traj%iPES(  iIniNew:iFinNew)    =   Traj%IPES(  iIniOld:iFinOld)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function AssessConvergence( This ) result( Converged )
! This procedure assesses if we can stop the integration process.
! There are two criteria for assessing the convergence:
! * if at least one atom has exeeeded the maximum interatomic distance or,
! * if the trajectory time has exceeded the maximum time
! Note that via 'Hamiltonian', potcl returns the inverses of the atom-atom distances.
! For a given trajectory, the Converged variable is set to 0 if the traj is not converged and to -1 if it is
! the variable used to be called 'iDone'

  class(Trajectories_Type)                  ,intent(in)     ::    This
  
  logical   ,dimension( This%NTraj )                        ::    Converged
  integer                                                   ::    iT        ! Index of trajectories
  integer                                                   ::    iP        ! Index of pairs
  
  Converged           =   .False.               ! By default, set all trajectories as NOT converged
  if ( This%BwrdIntegrationFlg ) then          ! Case of backward integration
  
    do iT = 1,This%NTraj
      if ( abs(This%t(iT)) >= 1.0E-4_rkp ) cycle
      Converged(iT)   =   .True.
    end do
    
  else                                    ! Case of forward integration

    TrajLoop: do iT = 1,This%NTraj
      if ( This%t(iT) > This%tmax ) then
        Converged(iT)   =   .True.
        cycle TrajLoop
      end if
      do iP = 1,size(This%Rpi,1)
        if ( This%Rpi(iP,iT) < This%RpiMax ) then
          Converged(iT) =   .True.
          cycle TrajLoop
        end if
      end do
    end do TrajLoop
    
  end if
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
