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

Module ImpactParameter_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger

  implicit none

  private
  public    ::    ImpactParameter_Type

  Type      ::    ImpactParamRing_Type
    integer                                         ::    NTraj   =   0       ! Number of trajectories in current ring
    integer                                         ::    Weight  =   0
    real(rkp)                                       ::    b       =   Zero
    integer                                         ::    iTrajMax=   0
  End Type

  Type      ::    ImpactParameter_Type
    integer                                                 ::    NRings      ! Number of impact parameter rings
    type(ImpactParamRing_Type)  ,dimension(:) ,allocatable  ::    Rings       ! Vector of ImpactParamRing objects
  contains
    procedure ,public   ::  Initialize  =>    InitializeImpactParameter
    procedure ,public   ::  GetRing     =>    GetImpactParamRing
    procedure ,public   ::  GetValue    =>    GetImpactParamValue
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeImpactParameter( This, Input, NTrajTotTemp, i_Debug )
! The variables 'This%Erel' is the initial relative translational energy in hartree. It can tak the following values:
!  * < 0: this energy will be determined from a Boltzmann distribution at abs(This%Erel) degrees k.
!  * = 1234321: the trajectories will have fixed total energy. Then it is necessary to read in This%Etot.
! The variables 'bstart' is the initial impact parameter.
!  * >= 0:      all trajectories will have impact parameter bstart.
!  * <  0:      the impact parameters will be abs(bstart)*sqrt(eta), where eta is a random number on (0,1)
!  * = 1234321: stratified sampling will be used. This requires additional input - nring,ntring,bring
!               These parameters are:
!                * nring: number of stratified sampling ranges,
!                * ntring: number of trajectories in a particular range, (-1 means there is no limit)
!                * bring: like bstart (e.g. signs), except it applies to the range.
!  * = 2345432: Impact parameter chosen so that the trajectory has a particular value of the total angular momentum.
!               This requires additional input, namely, This%TotalAngularMomentum, the total angular momentum.

  use Input_Class           ,only:  Input_Type

  class(ImpactParameter_Type)               ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    NTrajTotTemp      ! Total number of trajectories
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iR          ! Index of ring
  integer                                                   ::    SumNTraj
  integer                                                   ::    SumAllWeight
  integer                                                   ::    SumNegWeight
  integer                                                   ::    iTrajMax
  character(:)  ,allocatable                                ::    String
  integer                                                   ::    NTrajTot
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeImpactParameter")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  NTrajTot = NTrajTotTemp !floor( real(NTrajTotTemp) / real(Input%NPESs) )
  if (i_Debug_Loc) call Logger%Write( "-> Number of Traj.s per PES:   NTrajTot = ", NTrajTot )

! ==============================================================================================================
!     TREATING THE CASE WHEN THERE IS NO PARAMETER TO INITIALIZE
! ==============================================================================================================
  if ( Input%NRing == 0) then  ! Case when we are Initializing from Analysis
    if (i_Debug_Loc) then
      call Logger%Write( "Input%NRing = ", Input%NRing )
      call Logger%Write( "-> Nothing to do => Exiting')")
    end if
    call Logger%Exiting()
    return
  end if

  This%NRings     =   Input%NRing
  allocate( This%Rings(This%NRings) )
  do iR = 1,This%NRings
    associate( Ring => This%Rings(iR) )
      Ring%b        =   Input%bring(iR)
      Ring%Weight   =   Input%ntring(iR)
      Ring%NTraj    =   Input%ntring(iR)
      Ring%iTrajMax =   0
    end associate
  end do

!   if (i_Debug_Loc) then
!     call Logger%Write( "-> Number of stratified samping for the impact parameter: This%NRings = ", This%NRings )
!     do iR = 1,This%NRings
!     associate( Ring => This%Rings(iR) )
!     call Logger%Write( "   * iR = ", iR, "b = ", Ring%b, "Weight = ", Ring%Weight, Fi="i3", F4="es15.8" )
!     end associate
!     end do
!   end if

  SumNegWeight    =   0
  do iR = 1,This%NRings
    associate( Ring => This%Rings(iR) )
      if ( Ring%Weight < 0 ) SumNegWeight = SumNegWeight + 1
    end associate
  end do
!   if (i_Debug_Loc) call Logger%Write( "SumNegWeight = ", SumNegWeight )

  SumAllWeight    =   sum( This%Rings(:)%Weight )
  
  if ( SumNegWeight == 0 ) then
  
    SumNTraj      =   0
    do iR = 1,This%NRings
      associate( Ring => This%Rings(iR) )
        Ring%NTraj  =   int( dfloat( Ring%Weight * NTrajTot ) / SumAllWeight )
        SumNTraj    =   SumNTraj + Ring%NTraj
!         if (i_Debug_Loc) call Logger%Write( "iR = ", iR, "Ring%NTraj = ", Ring%NTraj )
      end associate
    end do
    if ( SumNTraj /= NTrajTot) This%Rings(1)%NTraj = This%Rings(1)%NTraj + NTrajTot - SumNTraj

    iTrajMax        =   0
    do iR = 1,This%NRings
      associate( Ring => This%Rings(iR) )
        Ring%iTrajMax =   iTrajMax + Ring%NTraj
        iTrajMax      =   Ring%iTrajMax
      end associate
    end do
    
    !This%Rings(This%NRings)%NTraj    = This%Rings(This%NRings)%NTraj    !+ mod(NTrajTotTemp, Input%NPESs)
    !This%Rings(This%NRings)%iTrajMax = This%Rings(This%NRings)%iTrajMax !+ mod(NTrajTotTemp, Input%NPESs)

  end if

  if (i_Debug_Loc) then
    call Logger%Write( "-> Number of stratified samping for the impact parameter: This%NRings = ", This%NRings )
    do iR = 1,This%NRings
      associate( Ring => This%Rings(iR) )
        if ( allocated(String) ) deallocate(String)
        if ( Ring%b < Zero ) then;  allocate( String , source = 'bmax' )
        else;                       allocate( String , source = 'b' ); end if
        if ( Ring%NTraj > 0 ) then
          call Logger%Write( "   * iR = ", iR, "NTraj = ", Ring%NTraj, "iTrajMax = ", Ring%iTrajMax, String//" = ", abs(Ring%b), F2="i2", F4="i9", F6="i9", F8="es15.8" )
        else
          call Logger%Write( "   * iR = ", iR, "NTraj = ", 'unlimited', "iTrajMax = ", Ring%iTrajMax, String//" = ", abs(Ring%b), F2="i2", F4="a9", F6="i9", F8="es15.8" )
        end if
      end associate
    end do
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________! Original David and Bruno
!Pure Subroutine GetImpactParamValue( This, Idx, b )
!  
!  class(ImpactParameter_Type)               ,intent(in)     ::    This
!  integer                                   ,intent(in)     ::    Idx         ! Global index of current trajectory
!  real(rkp)                                 ,intent(out)    ::    b           ! Impact parameter
!  
!  integer                                                   ::    iR          ! Index of ring
!  
!  do iR = 1,This%NRings
!    if ( Idx > This%Rings(iR)%iTrajMax ) cycle
!    b = This%Rings(iR)%b
!    exit
!  end do

!End Subroutine
!!--------------------------------------------------------------------------------------------------------------------------------! Original David and Bruno


!________________________________________________________________________________________________________________________________!
Pure Subroutine GetImpactParamValue( This, Idx, bmin, bmax )
  
  class(ImpactParameter_Type)               ,intent(in)     ::    This
  integer                                   ,intent(inout)  ::    Idx         ! Global index of current trajectory
  real(rkp)                                 ,intent(out)    ::    bmin          
  real(rkp)                                 ,intent(out)    ::    bmax           
    
  integer                                                   ::    iR          ! Index of ring

  Idx = Idx + 1

  do iR = 1,This%NRings
    if ( Idx > This%Rings(iR)%iTrajMax ) cycle
    bmax = This%Rings(iR)%b
    bmin = Zero
    if (iR>1) bmin = This%Rings(iR-1)%b
    exit
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine GetImpactParamRing( This, Idx, iRing)

  class(ImpactParameter_Type)               ,intent(in)     ::    This
  integer                                   ,intent(inout)  ::    Idx       ! Index of trajectory
  integer                                   ,intent(out)    ::    iRing     ! Index of ring corresoinding to the input trajectory
  
  integer                                                   ::    iR        ! Index of ring

  iRing   =   0
  do iR = 1,This%NRings
    if ( Idx > This%Rings(iR)%iTrajMax ) cycle
    iRing   =   iR
    exit
  end do

  !Idx = Idx + 1
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
