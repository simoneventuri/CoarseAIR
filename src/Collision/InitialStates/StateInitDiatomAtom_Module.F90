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

Module StateInitDiatomAtom_Module

  use Parameters_Module         ,only:  rkp
  use Logger_Class              ,only:  Logger
  use Error_Class               ,only:  Error
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type
  use File_Class                ,only:  File_Type
  use Timing_Module

  implicit none

  private
  public    ::    InitializeLevels_DiatomAtom
  public    ::    SetInitialState_DiatomAtom
  public    ::    ComputeKineticEnergy
  public    ::    BondLength
  public    ::    Compute_Coordinates_Velocities

  public    ::    ParamsFile, Params
  
  type(File_Type)                         ::    ParamsFile
  real(rkp)       ,dimension(19)          ::    Params
  real(rkp)                               ::    xmui_Global
  real(rkp)                               ::    Vc_R2_Global
  class(DiatomicPotential_Type) ,pointer  ::    DiaPot_Global                            !< Local pointer for the intra-nuclear diatomic potential object used in procedure 'FuncRHS'
  logical   ,parameter                    ::    i_Debug_Global = .False.
  
  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeLevels_DiatomAtom( Input, Species, i_Debug )
! Read in info for target and projectile initial states.
! The target is ab and the projectile is c, i.e. diatom+atom collisions.
! This procedure was initially called preini3
! Subroutine InitializeLevels_DiatomAtom( Input, Species, InitialStateSetter, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Species_Class         ,only:  Species_Type

  type(Input_Type)                          ,intent(in)     ::    Input
  type(Species_Type)  ,dimension(:)         ,intent(inout)  ::    Species                ! Species objects which always have 2 elements [Target,Projectile]. Since we are in the case of a Diatom-Atom collision, 
                                                                                         !    only the 1st element is used, which corresponds to the Target species. Dim=(NSpecies)=(2).                                                                            
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  character(:) ,allocatable                                 ::    FileName
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeLevels_DiatomAtom")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  FileName = './levels_' // trim(adjustl(Input%Molecules_Name(1))) // '.inp'          

  if (i_Debug_Loc) call Logger%Write( "Initializing the set of levels for a 'Diatom-Atom' collision" )
  if (i_Debug_Loc) call Logger%Write( " Calling Species(1)%ListStates%Initialize with FileName = ", FileName )
  call Species(1)%ListStates%Initialize( Input, 1, FileName=FileName, i_Debug=i_Debug_Loc )                                                
  if (i_Debug_Loc) call Logger%Write( " Done initializing the Species(1)%ListStates object" )

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetInitialState_DiatomAtom( Species, Qijk, dQijk, iTraj, iPES, i_Debug, i_Debug_Deep )
! @TODO: Remove The Traj,iTraj argument since it is not used
! This procedure sets the initial states.
! The procedure was initially called 'setis'.
! The target species (associated to the object 'TargSpecies') will always be a diatomic species with 2 atoms.
! The projectile species is not used here
! Subroutine SetInitialState_DiatomAtom( q1, q1dot, iste, irej, smin, smax, i_Debug )

  use Species_Class         ,only:  Species_Type
  use Parameters_Module     ,only:  Zero, One, Two, Half
  use RandomVector_Module   ,only:  RanddwsVec
  use Level_Class           ,only:  Level_Type

  type(Species_Type)  ,dimension(:)         ,intent(in)     ::    Species       ! Species objects which always have 2 elements [Target,Projectile]. Since we are in the case of a Diatom-Atom collision, 
                                                                                !    only the 1st element is used, which corresponds to the Target species. Dim=(NSpecies)=(2)
  real(rkp)           ,dimension(:,:,:)     ,intent(out)    ::    Qijk          ! Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species 
                                                                                !    (dim-3: 1=target, 2:projectile). Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp)           ,dimension(:,:,:)     ,intent(out)    ::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species 
                                                                                !    (dim-3: 1=target, 2:projectile). Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  integer                                   ,intent(in)     ::    iTraj
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    i
  real(rkp)                                                 ::    RandNum        ! Random number
  real(rkp)                                                 ::    Omega         ! Angular velocity for rotation
  real(rkp)                                                 ::    Vc_R2         ! Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                                 ::    Ekin
  real(rkp)                                                 ::    BL            ! bond length
  real(rkp)                                                 ::    dBL           ! Derivative of the ound length
  type(Level_Type)                                          ::    State
  integer                                                   ::    Traj_iState, Traj_jqn, Traj_vqn

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetInitialState_DiatomAtom")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


! ======================================================================================================================
!     SELECTING THE INTERNAL STATE (Traj_iState, Traj_jqn and Traj_vqn) AND COMPUTING THE CENTRIFUGAL POTENTIAL (Vc_R2)
! ======================================================================================================================
  if ( .not. (    ( ( trim(adjustl(Species(1)%ListStates%vInMethod)) .eq. "Fixed"        ) .and. ( trim(adjustl(Species(1)%ListStates%jInMethod)) .eq. "Fixed" ) ) &
               .or. ( trim(adjustl(Species(1)%ListStates%vInMethod)) .eq. "MostProbable" )                                                                         ) ) then     ! Case of sampled initial internal states
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Using a Thermal Distribution for Choosing Initial State')")

    RandNum = RanddwsVec(iPES)
    if (i_Debug_Loc) call Logger%Write("RAND_DA: ", RandNum)
    do i = 1,Species(1)%ListStates%NStates
      if ( Species(1)%ListStates%States(i)%rlim >= RandNum ) exit
    end do
    Traj_iState   =   i
    Traj_jqn      =   Species(1)%ListStates%States(i)%jqn
  else     
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Either Fixed State or State that Maximizes Boltzmann Distribution')")
    Traj_iState   =   Species(1)%ListStates%StateIn
    Traj_jqn      =   Species(1)%ListStates%jIn
  end if
  Traj_vqn        =   Species(1)%ListStates%States(i)%vqn
  if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Species(1)%DiaPot%xmui2 = ',es15.8,'; Traj_jqn = ',g0)") Species(1)%DiaPot%xmui2, Traj_jqn
  Vc_R2           =   Species(1)%DiaPot%xmui2 * ( Traj_jqn + Half )**2
  associate( S => Species(1)%ListStates%States(Traj_iState) )
    State%rmin    =   S%rmin
    State%tau     =   S%tau
    State%Eint    =   S%Eint
    State%ri      =   S%ri
    State%ro      =   S%ro
    State%egam    =   S%egam
    State%rmax    =   S%rmax
    State%vmin    =   S%vmin
    State%vmax    =   S%vmax
    State%rlim    =   Zero      ! NOT USED
  end associate
  
  
! ! ==============================================================================================================
! !   RE-COMPUTING THE STATE PROPERTIES
! ! ==============================================================================================================
! !!! SIMONE added this subroutine. For some (~20) of the N2 Quasi-Bound levels from NASA Ames (e.g., 7489 from n2.leroy.Sort.dat), 
! ! the difference between the maximum values of the diatomic potential from the Levels list and the one computed in the "period" 
! ! subroutine was of the order of 1.e-5; this fact was making CoarseAIR not entering the if "(Ekin .ge. (State%vmax - Emin)) )" 
! ! statement in the "ComputeKineticEnergy" subroutine, and finally crashing with the error: "NaN found in scalar variable"
! ! in "* Procedure: Integrate_ODE_1d", "* Variable:  h before mmidv".

!   if ( State%egam > Zero ) then
!     if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Found a Quasi Bound Level. Calling RecomputeLevelProperties')")
!     call Species(1)%DiaPot%RecomputeLevelProperties( State%Eint, Vc_R2, State%rMin, State%VMin, State%rMax, State%VMax, State%ri, State%ro, i_Debug=i_Debug_Loc)
!   end if
  
!   if (i_Debug_Loc) then
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Internal state:                   Traj_iState = ',g0)")     Traj_iState
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Rotational quantum number:        Traj_jqn    = ',g0)")     Traj_jqn
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Vibrational quantum number:       Traj_vqn    = ',g0)")     Traj_vqn
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Centrifugal potential times R2:   Vc_R2       = ',es15.8)") Vc_R2
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%rmin  = ',es15.8)") State%rmin
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%tau   = ',es15.8)") State%tau
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%Eint  = ',es15.8)") State%Eint
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%ri    = ',es15.8)") State%ri
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%ro    = ',es15.8)") State%ro
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%egam  = ',es15.8)") State%egam
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%rmax  = ',es15.8)") State%rmax
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%vmin  = ',es15.8)") State%vmin
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%vmax  = ',es15.8)") State%vmax
!     write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: ->                                   State%rlim  = ',es15.8)") State%rlim
!   end if
! ! ==============================================================================================================


! ==============================================================================================================
!   COMPUTING THE KINETIC ENERGY AND UPDATING SOME COMPONENT IN THE STATE OBJECT
! ==============================================================================================================
  if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Calling ComputeKineticEnergy')")
  call ComputeKineticEnergy( Ekin, State, Vc_R2, Species(1)%DiaPot, iPES, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) then
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Kinetic energy:                   Ekin      = ',es15.8)") Ekin
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Estimate of outter turning point: State%ro  = ',es15.8)") State%ro
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Period:                           State%Tau = ',es15.8)") State%Tau
  end if
! ==============================================================================================================


! ==============================================================================================================
!   COMPUTING THE INITIAL BOND LENGTH AND ITS DERIVATIVES
! ==============================================================================================================
  if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Computing the bond length and its derivatives. Calling BondLength')")
  Ekin = Zero                                                                                                                     ! Setting the kinetic energy at the potential minimum
  call BondLength( BL, dBL, State%ro, State%Tau, Ekin, Vc_R2, Species(1)%DiaPot, iPES, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Deep )! Computng BL/bBL output, all the rest is input
  if (i_Debug_Loc) then
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Bond length:        BL       = ',es15.8)") BL
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Bond length Deriv.: dBL      = ',es15.8)") dBL
  end if
! ==============================================================================================================


! ==============================================================================================================
!   COMPUTING THE MOLECULE COORDINATES AND VELOCITIES
! ==============================================================================================================
! choose random variables for diatomic orientation; set some variables defining the angular momentum vector Omega
! ==============================================================================================================
  if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Computing the diatomic coordinates')")
  ! Shouldn't it be:    Omega = sqrt(Two * Vc_R2 * xmui2) ??? (SIMONE)
  Omega = Two * sqrt(Vc_R2 * Species(1)%DiaPot%xmui2) / BL**2                                                                                         ! Setting the angular velocity for rotation
  if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Angular velocity: Omega = ',d20.10)")  Omega
                                                                                        
  if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Calling Compute_Coordinates_Velocities')")
  call Compute_Coordinates_Velocities( Species(1)%GetAtomsMass(), BL, dBL, Omega, Qijk(:,:,1), dQijk(:,:,1), iTraj, iPES, i_Debug=i_Debug_Loc )! Computing the coordinates of the target
  Qijk( :,:,2) = Zero                                                                                                             ! Setting the coordinates of the projectile to zero
  dQijk(:,:,2) = Zero                                                                                                             ! Setting the coordinates of the projectile to zero
  if (i_Debug_Loc) then
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Positions of atoms in target:')")
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]:    -> Atom A: ',*(d20.10,3x))")           Qijk(:,1,1)
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]:    -> Atom B: ',*(d20.10,3x))")           Qijk(:,2,1)
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: -> Velocities of atoms in target:')")
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]:    -> Atom A: ',*(d20.10,3x))")           dQijk(:,1,1)
    write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]:    -> Atom B: ',*(d20.10,3x))")           dQijk(:,2,1)
  end if
! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ComputeKineticEnergy( Ekin, State, Vc_R2, DiaPot, iPES, i_Debug )
! This procedure sets the initial kinetic energy of diatom for BondLength.
! For bound states, this is just the difference between the potential minimum and the internal energy.
! For quasibond states, choose the internal energy distributed about the resonance energy.
! Use a normal distribution with lower limit e=0 and upper limit e=vmax.
! This procedure was initially called 'tener'
! On input: See below for the sub-intent of State
!  *   Eint is the internal energy for bound states and eres for q.b. states
!  *   egam is Zero for bound states and gamma/2 for q.b. states
!  *   vmin is the potential minimum (including centrifugal potential)
!  *   vmax is the potential local maximum (including centrifugal pot.)
!  *   tau is the vibrational period for e = Eint
!  *   ri is the inner turning point for e = Eint
!  *   ro is the outter turning point for e = Eint
!  *   Vc_R2 is the centrifugal potential times r**2
!  *   xmui is 1/reduced mass
!  *   rmax is the distance where the maximum of the centrifugal barrier occurs
!  *   rmin is the distance where the minimum of the effective potential occurs
! On return:
!  *   Ekin is the initial kinetic energy.
!
!   Component of State        Attribute
!   ------------------------------------
!     Eint                      in
!     egam                      in
!     rmin                      in
!     vmin                      in
!     vmax                      in
!     tau                       out
!     ri                        in
!     ro                        inout
!     rmax                      in
!     rlim                      <NOT USED>

  use Parameters_Module     ,only:  Zero, One, Two
  use RandomVector_Module   ,only:  RanddwsVec
  use Level_Class           ,only:  Level_Type

  real(rkp)                                 ,intent(out)    ::    Ekin
  type(Level_Type)                          ,intent(inout)  ::    State
  real(rkp)                                 ,intent(in)     ::    Vc_R2                   !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  class(DiatomicPotential_Type)             ,intent(in)     ::    DiaPot                  !< Intra-nuclear diatomic potential object
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    UpdatePeriod
  integer                                                   ::    Iter
  real(rkp)                                                 ::    Emin                    ! Minimum energy distance between the edges (0, vmax) which the q.b. state energy is allowed to be.
  real(rkp)                                                 ::    Eave                    ! Averaged energy
  real(rkp)                                                 ::    v1, v2, r, fac, g1, g2  ! Variables related to randum numbers
  real(rkp)                                                 ::    rms
  real(rkp)                                                 ::    RandNumA
  real(rkp)                                                 ::    RandNumB
  integer   ,parameter                                      ::    IterMax = 100000
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ComputeKineticEnergy")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  

  Emin          =     Zero                                                                ! Emin=4.0E-6 ! nonZero made sense for resonance transfer theory, but not for state-to-state
  Eave          =     Zero
  UpdatePeriod  =     .False.

  if ( State%egam == Zero ) then
    Ekin          =   State%Eint
  else
    UpdatePeriod  =     .True.
    Iter          =   0
    do
      Iter    =   Iter+1
      if ( Iter .ge. IterMax ) exit
      
      RandNumA = RanddwsVec(iPES)
      RandNumB = RanddwsVec(iPES)
      if (i_Debug_Loc) call Logger%Write("RANDA: ", RandNumA)
      if (i_Debug_Loc) call Logger%Write("RANDB: ", RandNumB)
      v1      =   Two * RandNumA - One
      v2      =   Two * RandNumB - One
      r       =   v1*v1 + v2*v2
      if ( r .ge. One) cycle
      fac     =   sqrt( -Two*log(r) / r )
      g1      =   v1 * fac                                                                !   g1 and g2 are from a normal distribution
      g2      =   v2 * fac
      Ekin    =   State%egam * g1 + State%Eint
      if (i_Debug_Loc) call Logger%Write( " -> Ekin        = ", Ekin, "; State%vmax = ", State%vmax )
      if   ( (Ekin < Emin) .or. (Ekin .ge. (State%vmax - Emin)) ) then
        Ekin  =   State%egam * g2 + State%Eint
        if ( (Ekin < Emin) .or. (Ekin .ge. (State%vmax - Emin)) ) cycle
      end if
      exit
    end do
!     if (Iter>=IterMax) call Error( "[ComputeKineticEnergy] Maximum iteration reached" )
    if (Iter>=IterMax) write(Logger%Unit,"(12x,'[ComputeKineticEnergy] MAXIMUM ITERATION REACHED')")
    Eave      =   Eave + Ekin
  end if

! For quasibond states, it is necessary to generate new periods. This procedure will only change the variables:
!  * 'State%ro'   :   estimate of outter turning point
!  * 'State%Tau'  :   the period in a.u.
  if ( UpdatePeriod )  call DiaPot%Period( Ekin, Vc_R2, .True., State%rMin, State%VMin, State%rMax, State%VMax, State%ri, State%ro, State%Tau, i_Debug=i_Debug_Loc )

  Eave      =   Eave
  rms       =   sqrt( (Eave - Ekin)**2 )
  Ekin      =   Ekin - State%vmin
  
  Params(7) = Ekin
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! !________________________________________________________________________________________________________________________________!
! Subroutine Period( State, Eint, Vc_R2, xmu, DiaPot, RecompPropFlg, i_Debug )
! ! This Subroutine calculates the vibrational period, performing a numerical derivative of the classical action.
! ! In input:
! !     NTraj - number of trajectories
! !     ri - estimate of inner turning point
! !     ro - estimate of outter turning point
! !     Eint - internal energy
! !     xmu - reduced mass of oscillator
! !     vmax - maximum of centrifugal barrier
! !     rm - location of centrifugal barrier maximum.
! !     rn - location of the minimum of the effective potential
! !     on return:
! !     tau - the period in a.u.
! !     ro - the outter turning point at Eint
! ! Subroutine period( tau,ri,ro,Eint,cent,xmu,vmax,diat,rm,rn)
! ! This procedure will only change the variables:
! !  * 'State%ro'   :   estimate of outter turning point
! !  * 'State%Tau'  :   the period in a.u.

!   use Level_Class           ,only:  Level_Type
!   use Parameters_Module     ,only:  Two, Half

!   type(Level_Type)                          ,intent(inout)  ::    State
!   real(rkp)                                 ,intent(in)     ::    Eint
!   real(rkp)                                 ,intent(in)     ::    Vc_R2                            !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
!   real(rkp)                                 ,intent(in)     ::    xmu                              !< Reduced mass
!   class(DiatomicPotential_Type)             ,intent(in)     ::    DiaPot                           !< Intra-nuclear diatomic potential object
!   logical                                   ,intent(in)     ::    RecompPropFlg
!   logical                         ,optional ,intent(in)     ::    i_Debug

!   logical                                                   ::    i_Debug_Loc
!   integer   ,parameter                                      ::    NBisect   =   10                 ! number of bisection steps to perform
!   integer   ,parameter                                      ::    NNewton   =   20                 ! number of newton-raphson steps to perform
!   integer   ,parameter                                      ::    NPtsQuad  =   40
!   real(rkp) ,parameter                                      ::    eps       =   1.0E-10_rkp
!   real(rkp)                                                 ::    h
!   real(rkp)                                                 ::    e
!   real(rkp)                                                 ::    r11, r12, r21, r22, r1, r2
!   real(rkp)                                                 ::    Ap, Am                           ! Action integrals
!   real(rkp) ,dimension(2)                                   ::    RrangeExtreme, Rrange0, Rrange1, Rrange2
!   real(rkp)                                                 ::    vmin, vmax
!   integer                                                   ::    ierr

!   i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
!   if (i_Debug_Loc) call Logger%Entering( "Period")  !, Active = i_Debug_Loc )
!   !i_Debug_Loc   =     Logger%On()

!   if (i_Debug_Loc)  call Logger%Write( " -> Eint                              = ", Eint)


!   if ( RecompPropFlg ) then  
    
!     ! Looking for Diat. Potential Maximum and its Spatial Coordinate
!     RrangeExtreme = [State%rmin*1.01d0, 20.d0] 
!     call DiaPot%FindMaximum( RrangeExtreme, Vc_R2, eps, r22, vmax, ierr )
!     if (i_Debug_Loc)  call Logger%Write( " -> vmax                              = ", vmax, "; r(vmax) = ", r22)

  
!     ! Looking for Diat. Potential Minimum and its Spatial Coordinate
!     r11       =   State%ri * Half
!     Rrange0   =   [r11,r22]
!     r21       =   State%rmin
!     call DiaPot%FindMinimum( Rrange0, Vc_R2, eps, r21, vmin, ierr )
!     if (i_Debug_Loc)  call Logger%Write( " -> vmin                              = ", vmin, "; r(vmin) = ", r21)

!     r12       =   r21

!     Rrange1   =   [r11,r21]
!     Rrange2   =   [r12,r22]

!     ! Looking for Inner Tourning Point
!     call DiaPot%TurningPoint( Rrange1, Vc_R2, Eint, r1, ierr, NBisection=NBisect , NNewton=NNewton )
!     State%ri  =   r1
!     if (i_Debug_Loc)  call Logger%Write( " -> State%ri                          = ", State%ri)

!     ! Looking for Outer Tourning Point
!     call DiaPot%TurningPoint( Rrange2, Vc_R2, Eint, r2, ierr, NBisection=NBisect , NNewton=NNewton )
!     State%ro  =   r2
!     if (i_Debug_Loc)  call Logger%Write( " -> State%ro                          = ", State%ro)
                                       

!   else

!     r11       =   State%ri * Half
!     r12       =   State%rmin
!     r21       =   State%rmin
!     r22       =   State%rmax

!     Rrange1   =   [r11,r21]
!     Rrange2   =   [r12,r22]

!   end if
  

!   h         =   eps * Eint
!   e         =   Eint + h
!   if ( e > vmax ) then
!     h       =   Half*(vmax-Eint)
!     e       =   Eint + h
!   end if
!   if (i_Debug_Loc)  call Logger%Write( " -> DeltaE                            = ", h )
!   if (i_Debug_Loc)  call Logger%Write( " -> E+ = E + DeltaE                   = ", e )
!   call DiaPot%TurningPoint( Rrange1, Vc_R2, e, r1, ierr, NBisection=NBisect , NNewton=NNewton )
!   call DiaPot%TurningPoint( Rrange2, Vc_R2, e, r2, ierr, NBisection=NBisect , NNewton=NNewton )
!   if (i_Debug_Loc)  call Logger%Write( " -> Inner Tourning Point for E+       = ", r1 )
!   if (i_Debug_Loc)  call Logger%Write( " -> Outer Tourning Point for E+       = ", r2 )
!   Ap = DiaPot%ActionIntegral( r1, r2, Vc_R2, e, NPtsQuad, Two * xmu, QuadratureType=1 )
!   if (i_Debug_Loc)  call Logger%Write( " -> Action Integral for E+            = ", Ap )


!   e  = Eint - h
!   if (i_Debug_Loc)  call Logger%Write( " -> E- = E - DeltaE                   = ", e )
!   call DiaPot%TurningPoint( Rrange1, Vc_R2, e, r1, ierr, NBisection=NBisect , NNewton=NNewton )
!   call DiaPot%TurningPoint( Rrange2, Vc_R2, e, r2, ierr, NBisection=NBisect , NNewton=NNewton )
!   if (i_Debug_Loc)  call Logger%Write( " -> Inner Tourning Point for E-       = ", r1 )
!   if (i_Debug_Loc)  call Logger%Write( " -> Outer Tourning Point for E-       = ", r2 )
!   Am = DiaPot%ActionIntegral( r1, r2, Vc_R2, e, NPtsQuad, Two * xmu, QuadratureType=1 )
!   if (i_Debug_Loc)  call Logger%Write( " -> Action Integral for E-            = ", Am )

  
!   State%tau = Half * ( Ap - Am ) / h
!   if (i_Debug_Loc)  call Logger%Write( " -> State%tau                         = ", State%tau)
  

!   if (i_Debug_Loc) call Logger%Exiting

! End Subroutine
! !--------------------------------------------------------------------------------------------------------------------------------!


! !________________________________________________________________________________________________________________________________!
! Subroutine ResonanceWidth( State, Eint, Vc_R2, xmu, DiaPot )
! ! This Subroutine calculates the vibrational period, performing a numerical derivative of the classical action.
! ! In input:
! !     NTraj - number of trajectories
! !     ri - estimate of inner turning point
! !     ro - estimate of outter turning point
! !     Eint - internal energy
! !     xmu - reduced mass of oscillator
! !     vmax - maximum of centrifugal barrier
! !     rm - location of centrifugal barrier maximum.
! !     rn - location of the minimum of the effective potential
! !     on return:
! !     tau - the period in a.u.
! !     ro - the outter turning point at Eint
! ! Subroutine period( tau,ri,ro,Eint,cent,xmu,vmax,diat,rm,rn)
! ! This procedure will only change the variables:
! !  * 'State%ro'   :   estimate of outter turning point
! !  * 'State%Tau'  :   the period in a.u.

!   use Level_Class           ,only:  Level_Type
!   use Parameters_Module     ,only:  Two, Half, hbarH

!   type(Level_Type)                          ,intent(inout)  ::    State
!   real(rkp)                                 ,intent(in)     ::    Eint
!   real(rkp)                                 ,intent(in)     ::    Vc_R2                            !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
!   real(rkp)                                 ,intent(in)     ::    xmu                              !< Reduced mass
!   class(DiatomicPotential_Type)             ,intent(in)     ::    DiaPot                           !< Intra-nuclear diatomic potential object

!   integer   ,parameter                                      ::    NBisect   =   10                 ! number of bisection steps to perform
!   integer   ,parameter                                      ::    NNewton   =   20                 ! number of newton-raphson steps to perform
!   integer   ,parameter                                      ::    NPtsQuad  =   40
!   real(rkp) ,parameter                                      ::    eps       =   1.0E-5_rkp
!   real(rkp)                                                 ::    h
!   real(rkp)                                                 ::    e
!   real(rkp)                                                 ::    r11, r12, r21, r22, r1, r2
!   real(rkp)                                                 ::    Ap, Am                           ! Action integrals
!   real(rkp) ,dimension(2)                                   ::    RrangeExtreme, Rrange0, Rrange1, Rrange2
!   integer                                                   ::    ierr

!   r11       =   State%rmin 
!   r21       =   State%rmax
!   r12       =   r21
!   r22       =   1.d2
  
!   Rrange0   =   [r11,r22]
!   Rrange1   =   [r11,r21]
!   Rrange2   =   [r12,r22]
 
!   call DiaPot%TurningPoint( Rrange1, Vc_R2, Eint, r1, ierr, NBisection=NBisect , NNewton=NNewton )
!   call DiaPot%TurningPoint( Rrange2, Vc_R2, Eint, r2, ierr, NBisection=NBisect , NNewton=NNewton ) 
!   Ap = DiaPot%ActionIntegralGamma( r1, r2, Vc_R2, Eint, NPtsQuad, Two * xmu, QuadratureType=3 )

!   State%egam = Ap / State%tau

! End Subroutine
! !--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine BondLength( r, rdot, rm, tau, tm, Vc_R2, DiaPot, iPES, i_Debug, i_Debug_Deep )
! This procedures determines the initial vibrational displacement and velocity.
! On input:
!  * phas:  vibrational phase
!  * rm:    potential minimum for this j
!  * tau:   vibrational period
!  * tm:    kinetic energy at the potential minimum
!  * xmui:  inverse of the reduced mass.
! On return:
!  * r:     bond length
!  * rdot:  time derivative of the bond length
! Subroutine BondLength(r,rdot,phas,rm,tau,tm,NTraj,xmui,derd,cent,smax,smin,irej,iste)     ! Initial
! @TODO: REmove The Traj,iTraj argument since it is not used

  use Parameters_Module     ,only:  Zero, Two, Half
  use ODE_Solver            ,only:  ODE_Solver_Type
  use RandomVector_Module   ,only:  RanddwsVec

  real(rkp)                                 ,intent(out)    ::    r
  real(rkp)                                 ,intent(out)    ::    rdot
  real(rkp)                                 ,intent(in)     ::    rm                                        !< Initial solution vector: r
  real(rkp)                                 ,intent(in)     ::    tau                                       !< vibrational period
  real(rkp)                                 ,intent(in)     ::    tm                                        !< kinetic energy at the potential minimum: Always 0 in input variable
  real(rkp)                                 ,intent(in)     ::    Vc_R2                                     !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  class(DiatomicPotential_Type)     ,target ,intent(in)     ::    DiaPot                                    !< Intra-nuclear diatomic potential object
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  logical                                                   ::    i_Debug_Loc
  real(rkp)                                                 ::    VibPhase                                  ! Vibrational phase (random number)
  integer   ,parameter                                      ::    IterMax =   10000                         ! Maximum number of iteration
  integer   ,parameter                                      ::    N       =   2                             ! Number of coupled ODE
  real(rkp) ,parameter                                      ::    eps     =   1.0E-08_rkp                   ! Global tolerence for the integration
  integer   ,parameter                                      ::    NSteps  =   4   !4                          ! Number of stepsizes to use in the burlisch-stoer integration
  real(rkp) ,parameter                                      ::    epse    =   1.0E-10_rkp                   ! Relative error control parameter in the burlisch-stoer integration
  real(rkp) ,parameter                                      ::    Relax   =   Half                          ! Relaxation parameter for increasing stepsize in the burlisch-stoer integration
  real(rkp) ,parameter                                      ::    fnum    =   10.0_rkp !10.0_rkp

  integer                                                   ::    Iter
  real(rkp) ,dimension(2)                                   ::    y
  real(rkp) ,dimension(2)                                   ::    yscal                                     ! Absolute error control vector for the burlisch-stoer integration
  real(rkp)                                                 ::    t_Final                                   ! Final integration time
  real(rkp)                                                 ::    t                                         ! Current integration time
  real(rkp)                                                 ::    V                                         ! Potential
  real(rkp)                                                 ::    d0
  real(rkp)                                                 ::    dt                                        ! Timestep. The initial stepsize is:    dt = phas*tau/fnum
  integer                                                   ::    iste
  integer                                                   ::    irej
  real(rkp)                                                 ::    smin
  real(rkp)                                                 ::    smax
  type(ODE_Solver_Type)                                     ::    ODE
  integer   ,dimension(2)                                   ::    rdum0
  real(rkp) ,dimension(2,1)                                 ::    rdum1, rdum2
  real(rkp)                                                 ::    StartTime, EndTime
  real(rkp)                                                 ::    EkinVerify

  call CPU_Time( StartTime )
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "BondLength")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  VibPhase  =   RanddwsVec(iPES)                                                                                                  ! Getting random number
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: -> Vibrational phase:   VibPhase = ',es15.8)") VibPhase

!   Not sure the initialization is required here
  smin    =   Zero
  smax    =   Zero
  irej    =   0
  iste    =   0


! ==============================================================================================================
!   SET-UP GLOBAL VARIABLES
! ==============================================================================================================
! We need to store the variable associated to the centrifugal potential Vc_R2 and to the reduced mass inverse xmui
! in the module scope in order for them to be accessbile from the 'FuncRHS' procedure which is passed to the
! integrator procedure 'bsstepv'. There might be a better (ie. more oriented-object) way of doing this.
! This need to be done before calling the ODE integration procedure 'bsstepv'.
! ==============================================================================================================
  xmui_Global     =     DiaPot%xmui
  Vc_R2_Global    =     Vc_R2
  DiaPot_Global   =>    DiaPot
! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING THE ODE SOLVER OBJECT
! ==============================================================================================================
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: Calling ODE%Initialize')")
  call ODE%Initialize( eps=epse, NSteps=NSteps, Relax=Relax, NanCheck=.True., EvaluateRHS=FuncRHS, i_Debug=i_Debug_Loc )
  if (.False.) call FuncRHS( rdum0, rdum1, rdum2 )      ! Juste to remove the call warning: ‘funcrhs’ defined but not used
! ==============================================================================================================


! ==============================================================================================================
!   SETTING THE INTEGRATION PARAMETERS
! ==============================================================================================================
  t_Final   =   VibPhase * tau
  t         =   Zero                                                                                                              ! Setting the initial time
  dt        =   t_Final / fnum                                                                                                    ! Setting the initial timestep
  y(1)      =   rm                                                                                                                ! Setting the initial solution vector
  y(2)      =   sqrt( Two * DiaPot%xmui * tm )                                                                                    ! Setting the initial solution vector
  yscal(1)  =   rm * epse
  yscal(2)  =   yscal(1)
  V         =   DiaPot%DiatomicPotential( y(1) )
  d0        =   V  + Vc_R2/y(1)**2 + ( Half *y(2)**2 / DiaPot%xmui )
! ==============================================================================================================


! ==============================================================================================================
!   INTEGRATING FROM T=0 TO T=VibPhase*TAU
! ==============================================================================================================
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: Integrating in time')")
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: -> Initial solution:   [r,rdot] = ',*(es15.8,3x))") y
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: -> Final time:         t_Final  = ',es15.8)") t_Final
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: Calling ODE%Integrate')")
  Iter      =   0
  do
    Iter = Iter + 1
    dt   = min( dt, t_Final-t )                                                                                                   ! Setting the timestep: Make sure we do not integrate too far
    call ODE%Integrate( 1, y, t, dt, yscal, i_Debug=i_Debug_Deep )                                                                   ! Integrating forward by dt
    if (i_Debug_Loc) write(Logger%Unit,"(12x,'[BondLength]: -> Iter = ',i3,3x,'t = ',es15.8,3x,'abs(t-t_Final) = ',es15.8,3x,'[r,rdot] = ',*(es15.8,3x))") Iter, t, abs(t-t_Final), y
    if ( abs(t-t_Final) < eps ) exit                                                                                              ! If the resudial is higher than the tolerence, then go to the next time step, otherwise, we are done
    if ( Iter > IterMax ) exit
  end do
  r         =   y(1)
  rdot      =   y(2)
  if (i_Debug_Loc) then
    write(Logger%Unit,"(12x,'[BondLength]: Done integrating')")
    write(Logger%Unit,"(12x,'[BondLength]: -> Number of iterations: Iter = ',g0)")         Iter
    write(Logger%Unit,"(12x,'[BondLength]: -> Final solution: [r,rdot] = ',*(es15.8,3x))") y
  end if
  if ( Iter > IterMax ) then
    write(Logger%Unit,"(12x,'[BondLength]: <<< Failure >>> Iter reached IterMax')") ; error stop
  end if
  nullify(DiaPot_Global) 
  Params(8:9) = y(1:2)                                                                                                            !  Tear-down global variables
! ==============================================================================================================


  EkinVerify = Vc_R2/y(1)**2 + ( Half *y(2)**2 / DiaPot%xmui )
  if (i_Debug_Loc) then
    write(Logger%Unit,"(12x,'[BondLength]: -> EkinVerify = ',es15.8)") EkinVerify
  end if


! ! ==============================================================================================================
! !      CHECKING ENERGY CONSERVATION
! ! ==============================================================================================================
!   call DiaPot%Compute_V( y(1), V )

!
!   imax      =   0
!   Emax      =   zero
!   Eave      =   zero
!   do nt = 1,ntraj
!     ef(nt)  =   V(nt) + ( Half * ( rdot(nt)**2 )/xmui)+Vc_R2(nt)/(r(nt)**2)
!     V(nt)   =   abs(ef(nt)-d0(nt))
!     Eave    =   Eave + ef(nt)
!     if ( V(nt)> Emax ) then
!       Emax  =   V(nt)
!       imax  =   nt
!     end if
!   end do
!   Eave      =   abs( Eave / float(ntraj) )
!
!
! !   ef        =   V + ( Half * rdot**2 / xmui ) + Vc_R2 / r**2
! !   V         =   abs( ef - d0 )
! !   Eave      =   abs( Eave + ef )
! !   Emax      =   V
!   if ( Emax > max(Eave*1.d-6,1d-9) ) then
!     write(Logger%Unit,"(12x,'[BondLength]: energy not well converged')")
!     write(Logger%Unit,"(12x,'[BondLength]: abs error =',1pe15.7,' average Eint =',e15.7)")Emax,Eave
!     stop
!   end if
! ! ==============================================================================================================

  call CPU_Time( EndTime )
  t_vphase  =   t_vphase + EndTime - StartTime
  i_vphase  =   i_vphase + 1
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_Coordinates_Velocities( Masses, r, dr, Omega, q, qdot, iTraj, iPES, i_Debug )
! This procedure computes the initial coordinates and velocities for a diatomic molecule.
! The coordinates are cartesian with the center of mass at the origin.
! This procedure was initially called 'inid'.

  use Parameters_Module     ,only:  One, Two, Half, TwoPi
  use RandomVector_Module   ,only:  RanddwsVec

  real(rkp) ,dimension(2)                   ,intent(in)     ::    Masses        !> Mass of atoms
  real(rkp)                                 ,intent(in)     ::    r             !> Bond length
  real(rkp)                                 ,intent(in)     ::    dr            !> Bond length derivative
  real(rkp)                                 ,intent(in)     ::    Omega         !> initial angular velocity for rotation
  real(rkp) ,dimension(3,2)                 ,intent(out)    ::    q             !> Coordinates stored in the order x,y,z for atom a, x,y,z for atom b. (DIM=dof,NAtoms)
  real(rkp) ,dimension(3,2)                 ,intent(out)    ::    qdot          !> Velocities stored in the order x,y,z for atom a, x,y,z for atom b.
  integer                                   ,intent(in)     ::    iTraj
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  real(rkp)                                                 ::    ma, mb        ! Mass of atom A and B
  real(rkp)                                                 ::    fa, fb
  real(rkp)                                                 ::    gam, alp, bet
  real(rkp)                                                 ::    EkinVerify
  real(rkp)                                                 ::    RandNumA, RandNumB, RandNumC

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Compute_Coordinates_Velocities")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  ma  =   Masses(1)
  mb  =   Masses(2)
  fa  =   mb / ( ma + mb )
  fb  = - ma / ( ma + mb )

! Setting the rotational angles: phase, azimuthal and inclination angles
  RandNumA = RanddwsVec(iPES)
  RandNumB = RanddwsVec(iPES)
  RandNumC = RanddwsVec(iPES)
  if (i_Debug_Loc) call Logger%Write("RANDA: ", RandNumA)
  if (i_Debug_Loc) call Logger%Write("RANDB: ", RandNumB)
  if (i_Debug_Loc) call Logger%Write("RANDC: ", RandNumC)

  gam =   RandNumA * TwoPi                        ! Setting the azimuthal angle of the vector perpendicular to the rotation: Random
  bet =   acos( Two * (RandNumB-Half) )           ! Setting the inclination angle of the vector perpendicular to the rotation: Random
  alp =   RandNumC * TwoPi                        ! Setting the rotational phase angle: Random
  
  if (i_Debug_Loc) then
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> Azimuthal Angle:    alpha = ',es15.8,' = ',es15.8)") alp, alp * 180.d0 / (TwoPi/Two)
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> Inclination Angle:  beta  = ',es15.8,' = ',es15.8)") bet, bet * 180.d0 / (TwoPi/Two)
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> Rotational Phase:   gamma = ',es15.8,' = ',es15.8)") gam, gam * 180.d0 / (TwoPi/Two)
  end if
  
  Params(1)   = real(iTraj,8)
  Params(3:5) = [alp * 180.d0 / (TwoPi/Two), bet * 180.d0 / (TwoPi/Two), gam * 180.d0 / (TwoPi/Two)]


!! DAVID'S VERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! construct coordinates of atoms
  q(1,1)    =     fa * r * (   cos(bet) * cos(alp) * cos(gam) - sin(alp) * sin(gam) )                                             
  q(2,1)    =     fa * r * ( - cos(bet) * cos(alp) * sin(gam) - sin(alp) * cos(gam) )                                             
  q(3,1)    =     fa * r * sin(bet) * cos(alp)
  q(1,2)    =     fb * r * (   cos(bet) * cos(alp) * cos(gam) - sin(alp) * sin(gam) )                                             
  q(2,2)    =     fb * r * ( - cos(bet) * cos(alp) * sin(gam) - sin(alp) * cos(gam) )                                             
  q(3,2)    =     fb * r * sin(bet) * cos(alp)
  
! construct velocities of atoms
  qdot(1,1) =     fa * ( dr * (   cos(bet) * cos(alp) * cos(gam) - sin(alp) * sin(gam) ) + r * Omega * (   cos(bet) * sin(alp) * cos(gam) + sin(gam) * cos(alp) ) )    
  qdot(2,1) =     fa * ( dr * ( - cos(bet) * sin(gam) * cos(alp) - sin(alp) * cos(gam) ) + r * Omega * ( - cos(bet) * sin(alp) * sin(gam) + cos(alp) * cos(gam) ) )    
  qdot(3,1) =     fa * ( dr * (   sin(bet) * cos(alp) )                                  + r * Omega * sin(bet) * sin(alp) )
  qdot(1,2) =     fb * ( dr * (   cos(bet) * cos(alp) * cos(gam) - sin(alp) * sin(gam) ) + r * Omega * (   cos(bet) * sin(alp) * cos(gam) + sin(gam) * cos(alp) ) )    
  qdot(2,2) =     fb * ( dr * ( - cos(bet) * sin(gam) * cos(alp) - sin(alp) * cos(gam) ) + r * Omega * ( - cos(bet) * sin(alp) * sin(gam) + cos(alp) * cos(gam) ) )    
  qdot(3,2) =     fb * ( dr * (   sin(bet) * cos(alp) )                                  + r * Omega * sin(bet) * sin(alp) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (i_Debug_Loc) then
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> fa = ',es15.8,'; fb = ',es15.8)") fa, fb
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> q(:,1) = ',es15.8,', ',es15.8,', ',es15.8)") q(:,1)
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> q(:,2) = ',es15.8,', ',es15.8,', ',es15.8)") q(:,2)
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> qdot(:,1) = ',es15.8,', ',es15.8,', ',es15.8)") qdot(:,1)
    write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> qdot(:,2) = ',es15.8,', ',es15.8,', ',es15.8)") qdot(:,2)
  end if  
  
  EkinVerify = Half * Masses(1) * (qdot(1,1)**2 + qdot(2,1)**2 + qdot(3,1)**2) + Half * Masses(2) * (qdot(1,2)**2 + qdot(2,2)**2 + qdot(3,2)**2)             
  if (i_Debug_Loc) write(Logger%Unit,"(12x,'[Compute_Coordinates_Velocities]: -> EkinVerify = ',es15.8)") EkinVerify
  
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

!________________________________________________________________________________________________________________________________!
Pure Subroutine FuncRHS( iPES, r, drdx )
! This procedure compute the derivatives of the position and velocity.
! We are integrating a single oscillator.
! The interface must conform the one used in the ODE integrator procedure 'bsstepv'
! As a result, several arguments are required even if they are not used.
! This procedure will alway be called with 'NTraj=1'.

  use Parameters_Module       ,only:  Two
  
  integer   ,dimension(:)             ,intent(in)     ::    iPES
  real(rkp) ,dimension(:,:)           ,intent(in)     ::    r         ! Particle position and velocity. Dim=(NVar,NTraj)=(2,1)
  real(rkp) ,dimension(:,:)           ,intent(out)    ::    drdx      ! Derivatives of the particle position and velocity wrt position. Dim=(NVar,NTraj)=(2,1)

  real(rkp) ,dimension( size(r,2) )                   ::    Pos       ! Position
  real(rkp) ,dimension( size(r,2) )                   ::    Vel       ! Velocity
  real(rkp) ,dimension( size(r,2) )                   ::    V         ! Potential (DIM=NTraj)
  real(rkp) ,dimension( size(r,2) )                   ::    dV        ! Potential derivative (DIM=NTraj)

  Pos     =   r(1,:)
  Vel     =   r(2,:)
  call DiaPot_Global%Compute_Vd_dVd( Pos, V, dV )                                                                              !
  drdx(1,:) =   Vel
  drdx(2,:) =   - xmui_Global * ( dV - Two * Vc_R2_Global / Pos**3 )
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
