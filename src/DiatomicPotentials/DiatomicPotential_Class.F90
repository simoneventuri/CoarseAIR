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

Module DiatomicPotential_Class

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two
  use Logger_Class          ,only:  Logger

  implicit none

  private
  public  ::    DiatomicPotential_Type

  Type    ,abstract                       ::    DiatomicPotential_Type
    logical                               ::    Initialized
    character(:)  ,allocatable            ::    Name
    character(:)  ,allocatable            ::    SpeciesName
    integer                               ::    iMol
    real(rkp)                             ::    RedMass
    real(rkp)                             ::    xmui
    real(rkp)                             ::    xmui2
  contains
    private
    
    procedure                                     ,public ::    Initialize    =>    Initialize_DiatPot
    procedure                                     ,public ::    Output        =>    Output_DiatPot
    procedure                                     ,public ::    Period
    procedure                                     ,public ::    ResonanceWidth
    procedure                                     ,public ::    CheckMaxAndMin
    procedure                                     ,public ::    CheckTurningPoints
    procedure                                     ,public ::    TurningPoint
    procedure                                     ,public ::    ActionIntegral
    procedure                                     ,public ::    ActionIntegralGamma
    procedure                                     ,public ::    FindMinimum
    procedure                                     ,public ::    FindMaximum

    procedure                                     ,public ::    Compute_Vd_dVd    => Compute_Vd_dVd_DiatPot
    procedure                                     ,public ::    DiatomicPotential => DiatomicPotential_DiatPot

    generic                                       ,public ::    CentrifualPotential => CenPot_0d, CenPot_1d       
    procedure                             ,nopass         ::    CenPot_0d
    procedure                             ,nopass         ::    CenPot_1d

    generic                                       ,public ::    EffectivePotential  => EffPot_From_R_Vc_0d,    EffPot_From_R_Vc_1d
    procedure                                             ::    EffPot_From_R_Vc_0d
    procedure                                             ::    EffPot_From_R_Vc_1d
    
  End Type

  integer   ,parameter                                    ::    NBisectGlobal      =   30                 ! number of bisection steps to perform
  integer   ,parameter                                    ::    NNewtonGlobal      =   50                 ! number of newton-raphson steps to perform
  integer   ,parameter                                    ::    NPtsQuadGlobal     =   40
  real(rkp) ,parameter                                    ::    EpsStatPtsGlobal   =   1.0E-10_rkp
  real(rkp) ,parameter                                    ::    RStepStatPtsGlobal =   0.005_rkp
  real(rkp) ,parameter                                    ::    rEpsResWidthGlobal =   0.005_rkp

  Interface           CentrifualPotential
  
    Module Procedure  CenPot_0d
    Module Procedure  CenPot_1d
    
  End interface

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_DiatPot( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class   ,only:  Input_Type

  class(DiatomicPotential_Type)             ,intent(out)    ::    This          !< Diatomic-Potential object
  type(Input_Type)                          ,intent(in)     ::    Input
  character(:)  ,allocatable                ,intent(in)     ::    SpeciesName
  integer                                   ,intent(in)     ::    iMol
  real(rkp)                                 ,intent(in)     ::    Mass1
  real(rkp)                                 ,intent(in)     ::    Mass2
  logical ,optional                         ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_DiatPot" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = '' )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =   iMol
  This%Initialized  =   .True.

  This%RedMass = Mass1
  if (Mass2 /= Zero) This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_DiatPot( This, Unit )

  class(DiatomicPotential_Type)             ,intent(in)     ::    This          !< Diatomic-Potential object
  integer                                   ,intent(in)     ::    Unit
  
  integer                                                   ::    idum
  logical                                                   ::    ldum
  
  idum  =   Unit              ! Just to remove warning message about unused argument
  ldum  =   This%Initialized  ! Just to remove warning message about unused argument
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_DiatPot( This, R, V, dV )
  use Parameters_Module    ,only:  rkp
  class(DiatomicPotential_Type) ,intent(in)     ::    This          !< Diatomic-Potential object
  real(rkp)                     ,intent(in)     ::    R             !< Internuclear distance [bohr]
  real(rkp)                     ,intent(out)    ::    V             !< Diatomic potential energy [hartree]
  real(rkp)                     ,intent(out)    ::    dV            !< First derivative of the diatomic potential energy [hartree/bohr]
  V  = Zero
  dV = Zero
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_DiatPot( This, R ) result( V )
  use Parameters_Module    ,only:  rkp
  class(DiatomicPotential_Type) ,intent(in)     ::    This          !< Diatomic-Potential object
  real(rkp)                     ,intent(in)     ::    R             !< Internuclear distance [bohr]
  real(rkp)                                     ::    V             !< Diatomic potential energy [hartree]
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


! ********************************************************************************************************************************
! ********************************************************************************************************************************
!                                        PROCEDURES FOR COMPUTING THE CENTRIFUAL POTENTIAL                                        
! ********************************************************************************************************************************
! ********************************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Pure Function CenPot_0d( R, Vc_R2 ) result(Vc)

  real(rkp)                                 ,intent(in)     ::    R             !< Internuclear distance [bohr]
  real(rkp)                                 ,intent(in)     ::    Vc_R2         !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  
  real(rkp)                                                 ::    Vc            ! Centrifual potential
  
  Vc        =   Vc_R2 / R**2
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function CenPot_1d( R, Vc_R2 ) result(Vc)

  real(rkp) ,dimension(:)                   ,intent(in)     ::    R             !< Internuclear distance [bohr]
  real(rkp)                                 ,intent(in)     ::    Vc_R2         !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  
  real(rkp) ,dimension( size(R) )                           ::    Vc            ! Centrifual potential
  
  Vc        =   Vc_R2 / R**2
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!



! ********************************************************************************************************************************
! ********************************************************************************************************************************
!                                        PROCEDURES FOR COMPUTING THE EFFECTIVE POTENTIAL
! ********************************************************************************************************************************
! ********************************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Pure Function EffPot_From_R_Vc_0d( This, R, Vc_R2 ) result(Ve)

  class(DiatomicPotential_Type)             ,intent(in)     ::    This          !< Diatomic-Potential object
  real(rkp)                                 ,intent(in)     ::    R             !< Internuclear distance [bohr]
  real(rkp)                                 ,intent(in)     ::    Vc_R2         !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  
  real(rkp)                                                 ::    Ve            !< Effective potential [hartree]
  
  Ve      =   This%DiatomicPotential(R) + CentrifualPotential(R,Vc_R2)
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function EffPot_From_R_Vc_1d( This, R, Vc_R2 ) result(Ve)

  class(DiatomicPotential_Type)             ,intent(in)     ::    This          !< Diatomic-Potential object
  real(rkp) ,dimension(:)                   ,intent(in)     ::    R             !< Internuclear distance [bohr]
  real(rkp)                                 ,intent(in)     ::    Vc_R2         !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  
  real(rkp) ,dimension( size(R) )                           ::    Ve            !< Effective potential [hartree]
  
  Ve      =   This%DiatomicPotential(R) + CentrifualPotential(R,Vc_R2)
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Period( This, Eint, Vc_R2, RecompPropFlg, rMin, VMin, rMax, VMax, rIn, rOut, Tau, i_Debug )
! This Subroutine calculates the vibrational period, performing a numerical derivative of the classical action.
! In input:
!     NTraj - number of trajectories
!     ri - estimate of inner turning point
!     ro - estimate of outter turning point
!     Eint - internal energy
!     vmax - maximum of centrifugal barrier
!     rm - location of centrifugal barrier maximum.
!     rn - location of the minimum of the effective potential
!     on return:
!     tau - the period in a.u.
!     ro - the outter turning point at Eint
! Subroutine period( tau,ri,ro,Eint,cent,xmu,vmax,diat,rm,rn)
! This procedure will only change the variables:
!  * 'State%ro'   :   estimate of outter turning point
!  * 'State%Tau'  :   the period in a.u.

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                             !< Intra-nuclear diatomic potential object
  real(rkp)                                 ,intent(in)     ::    Eint
  real(rkp)                                 ,intent(in)     ::    Vc_R2                            !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  logical                                   ,intent(in)     ::    RecompPropFlg
  real(rkp)                                 ,intent(inout)  ::    rMin
  real(rkp)                                 ,intent(inout)  ::    VMin
  real(rkp)                                 ,intent(inout)  ::    rMax
  real(rkp)                                 ,intent(inout)  ::    VMax
  real(rkp)                                 ,intent(inout)  ::    rIn
  real(rkp)                                 ,intent(inout)  ::    rOut
  real(rkp)                                 ,intent(out)    ::    Tau
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    h
  real(rkp)                                                 ::    Em, Ep
  real(rkp)                                                 ::    r11, r12, r21, r22, r1, r2
  real(rkp)                                                 ::    Ap, Am                           ! Action integrals
  real(rkp) ,dimension(2)                                   ::    RrangeExtreme, Rrange0, Rrange1, Rrange2
  integer                                                   ::    ierr

  logical                                                   ::    i_Debug_Loc

  ! i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  ! if (i_Debug_Loc) call Logger%Entering( "Period")  !, Active = i_Debug_Loc )
  ! !i_Debug_Loc   =     Logger%On()

  ! if (i_Debug_Loc)  call Logger%Write( " -> Eint                              = ", Eint)


  if ( RecompPropFlg ) then  
    
    ! Looking for Diat. Potential Maximum and its Spatial Coordinate
    RrangeExtreme = [rMin*1.01d0, 20.d0] 
    call This%FindMaximum( RrangeExtreme, Vc_R2, EpsStatPtsGlobal, r22, vmax, ierr )
    ! if (i_Debug_Loc)  call Logger%Write( " -> VMax                              = ", vmax, "; r(VMax) = ", r22)

  
    ! Looking for Diat. Potential Minimum and its Spatial Coordinate
    r11       =   rIn * Half
    Rrange0   =   [r11,r22]
    r21       =   rMin
    call This%FindMinimum( Rrange0, Vc_R2, EpsStatPtsGlobal, r21, vmin, ierr )
    ! if (i_Debug_Loc)  call Logger%Write( " -> VMin                              = ", vmin, "; r(VMin) = ", r21)

    r12       =   r21

    Rrange1   =   [r11,r21]
    Rrange2   =   [r12,r22]

    ! Looking for Inner Tourning Point
    call This%TurningPoint( Rrange1, Vc_R2, Eint, r1, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
    rIn  =   r1
    ! if (i_Debug_Loc)  call Logger%Write( " -> rIn                                = ", rIn)

    ! Looking for Outer Tourning Point
    call This%TurningPoint( Rrange2, Vc_R2, Eint, r2, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
    rOut  =   r2
    ! if (i_Debug_Loc)  call Logger%Write( " -> rOut                               = ", rOut)
                                       

  else

    r11       =   rIn * Half
    r12       =   rMin
    r21       =   rMin
    r22       =   rMax

    Rrange1   =   [r11,r21]
    Rrange2   =   [rOut,r22]

  end if
  

  h         =   EpsStatPtsGlobal
  Ep        =   Eint + h
  if ( Ep > VMax ) then
    h       =   Half*(VMax-Eint)
    Ep      =   Eint + h
  end if

  
  ! if (i_Debug_Loc)  call Logger%Write( " -> DeltaE                            = ", h )
  ! if (i_Debug_Loc)  call Logger%Write( " -> E+ = E + DeltaE                   = ", Ep )

  r1 = Zero
  r2 = Zero
  call This%TurningPoint( Rrange1, Vc_R2, Ep, r1, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
  call This%TurningPoint( Rrange2, Vc_R2, Ep, r2, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Inner Tourning Point for E+          = ", r1 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Outer Tourning Point for E+          = ", r2 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> ierr for Outer Tourning Point for E+ = ", ierr )

  Ap = This%ActionIntegral( r1, r2, Vc_R2, Ep, NPtsQuadGlobal, Two * This%RedMass, QuadratureType=1 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Action Integral for E+               = ", Ap )


  Em = Eint - h
  ! if (i_Debug_Loc)  call Logger%Write( " -> E- = E - DeltaE                      = ", Em )

  r1 = Zero
  r2 = Zero
  Rrange2   =   [rMin,rOut]
  call This%TurningPoint( Rrange1, Vc_R2, Em, r1, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
  call This%TurningPoint( Rrange2, Vc_R2, Em, r2, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Inner Tourning Point for E-          = ", r1 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Outer Tourning Point for E-          = ", r2 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> ierr for Outer Tourning Point for E- = ", ierr )

  Am = This%ActionIntegral( r1, r2, Vc_R2, Em, NPtsQuadGlobal, Two * This%RedMass, QuadratureType=1 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Action Integral for E-               = ", Am )

  
  Tau = Half * ( Ap - Am ) / h
  ! if (i_Debug_Loc)  call Logger%Write( " -> Tau                                  = ", Tau)
  

  ! if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ResonanceWidth( This, Eint, Vc_R2, rMin, rMax, Tau, Egam, i_Debug )
! This Subroutine calculates 

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                             !< Intra-nuclear diatomic potential object
  real(rkp)                                 ,intent(in)     ::    Eint
  real(rkp)                                 ,intent(in)     ::    Vc_R2                            !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    rMin                             !< Location of Diat Potential Minimum
  real(rkp)                                 ,intent(in)     ::    rMax                             !< Location of Diat Potential Maximum
  real(rkp)                                 ,intent(in)     ::    Tau                              !< Vibrational Period
  real(rkp)                                 ,intent(out)    ::    Egam                             
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    r11, r12, r21, r22, r1, r2
  real(rkp)                                                 ::    Ap                               ! Action integrals
  real(rkp) ,dimension(2)                                   ::    Rrange0, Rrange1, Rrange2
  integer                                                   ::    ierr
  real(rkp)                                                 ::    rTemp, VTemp

  logical                                                   ::    i_Debug_Loc

  ! i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  ! if (i_Debug_Loc) call Logger%Entering( "ResonanceWidth")  !, Active = i_Debug_Loc )
  ! !i_Debug_Loc   =     Logger%On()

  r11   = rMin 
  ! if (i_Debug_Loc)  call Logger%Write( " -> r11                               = ", r11)

  r21   = rMax
  r12   = rMax
  ! if (i_Debug_Loc)  call Logger%Write( " -> r12 = r21                         = ", r12)

  ! r22   = 1.e2 Original from David
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simone
  !!! For Diat. Pot.s in Polynomial Form, we might have Oscillations around the dissociation energy;
  !!! Through the next lines we exclude the contributions to the ActionIntegralGamma coming from 
  !!! Regions Below the Dissociation Value. 
  !!! 
  rTemp = rMax
  VTemp = This%EffectivePotential( rTemp, Vc_R2 )
  do while (VTemp >= Eint)
    rTemp = rTemp + rEpsResWidthGlobal
    VTemp = This%EffectivePotential( rTemp, Vc_R2 )
  end do
  r22 = rTemp
  ! if (i_Debug_Loc)  call Logger%Write( " -> r11                               = ", r11)
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simone

  Rrange1   =   [r11,r21]
  Rrange2   =   [r12,r22]
  ! if (i_Debug_Loc)  call Logger%Write( " -> Rrange1                           = ", Rrange1)
  ! if (i_Debug_Loc)  call Logger%Write( " -> Rrange2                           = ", Rrange2)

  call This%TurningPoint( Rrange1, Vc_R2, Eint, r1, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal )
  ! if (i_Debug_Loc)  call Logger%Write( " -> r1                                = ", r1)

  call This%TurningPoint( Rrange2, Vc_R2, Eint, r2, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal ) 
  ! if (i_Debug_Loc)  call Logger%Write( " -> r2                                = ", r2)

  Ap = This%ActionIntegralGamma( r1, r2, Vc_R2, Eint, NPtsQuadGlobal, Two * This%RedMass, QuadratureType=3 )
  ! if (i_Debug_Loc)  call Logger%Write( " -> Ap                                = ", Ap)

  Egam = Ap / Tau
  ! if (i_Debug_Loc)  call Logger%Write( " -> Egam                              = ", Egam)


  ! if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine CheckMaxAndMin( This, State, iState, i_Debug )
! This Subroutine 

  use Level_Class           ,only:  Level_Type

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                           !< Intra-nuclear diatomic potential object
  type(Level_Type)                          ,intent(inout)  ::    State 
  integer                                   ,intent(in)     ::    iState
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    h
  real(rkp)                                                 ::    e
  real(rkp) ,dimension(2)                                   ::    RrangeExtreme, Rrange0, Rrange1, Rrange2
  integer                                                   ::    ierr
  real(rkp)                                                 ::    rMaxErrAbs, VMaxErrAbs, rMinErrAbs, VMinErrAbs
  real(rkp)                                                 ::    rMaxErrRel, VMaxErrRel, rMinErrRel, VMinErrRel
  logical                                                   ::    i_Debug_Loc
  
  real(rkp) ,parameter                                      ::    ToleranceR = 1.d-3 !1.d-5
  real(rkp) ,parameter                                      ::    ToleranceV = 1.d-3 !1.d-8

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "CheckMaxAndMin")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc)  call Logger%Write( " -> State%jqn = ", State%jqn)
  
  
  RrangeExtreme = [State%rMax*0.8d0, State%rMax*1.2d0] 
  call This%FindMaximum( RrangeExtreme, State%Vc_R2, EpsStatPtsGlobal, State%RMaxNew, State%VMaxNew, ierr )
  if (i_Debug_Loc)  call Logger%Write( " -> VMax    = ", State%VMax,    "; r(VMax) = ", State%rMax)
  if (i_Debug_Loc)  call Logger%Write( " -> VMaxNew = ", State%VMaxNew, "; r(VMax) = ", State%rMaxNew)
  
  rMaxErrAbs = dabs( State%rMax - State%rMaxNew ) 
  rMaxErrRel = dabs( rMaxErrAbs / State%rMax )
  if ( rMaxErrRel > ToleranceR ) then
    if (i_Debug_Loc)  call Logger%Write( " WARNING! Found a difference between the r@VMax from Levels List and the one ricomputed. Absolute Error = ", rMaxErrAbs, "; Relative Error = ", rMaxErrRel*100.d0, "%" )
  end if
  
  VMaxErrAbs = dabs( State%VMax - State%VMaxNew )
  VMaxErrRel = dabs( VMaxErrAbs / State%VMax )
  if ( VMaxErrRel > ToleranceV ) then
    if (i_Debug_Loc)  call Logger%Write( " WARNING! Found a difference between the VMax from Levels List and the one ricomputed. Absolute Error = ", VMaxErrAbs, "; Relative Error = ", VMaxErrRel*100.d0, "%" )
  end if
  

  Rrange0 = [State%rMin*0.8d0, State%rMin*1.2d0]
  call This%FindMinimum( Rrange0, State%Vc_R2, EpsStatPtsGlobal, State%rMinNew, State%VMinNew, ierr )
  if (i_Debug_Loc)  call Logger%Write( " -> VMin    = ", State%VMin,    "; r(VMin) = ", State%rMin)
  if (i_Debug_Loc)  call Logger%Write( " -> VMinNew = ", State%VMinNew, "; r(VMin) = ", State%rMinNew)
  
  rMinErrAbs = dabs( State%rMin - State%rMinNew )
  rMinErrRel = dabs( rMinErrAbs / State%rMin )
  if ( rMinErrRel > ToleranceR ) then
    if (i_Debug_Loc)  call Logger%Write( " WARNING! Found a difference between the r@VMin from Levels List and the one ricomputed. Absolute Error = ", rMinErrAbs, "; Relative Error = ", rMinErrRel*100.d0, "%" )
  end if
  
  VMinErrAbs = dabs( State%VMin - State%VMinNew )
  VMinErrRel = dabs( VMinErrAbs / State%VMin )
  if ( VMinErrRel > ToleranceV ) then
    if (i_Debug_Loc)  call Logger%Write( " WARNING! Found a difference between the VMin from Levels List and the one ricomputed. Absolute Error = ", VMinErrAbs, "; Relative Error = ", VMinErrRel*100.d0, "%" )
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine CheckTurningPoints( This, State, iState, i_Debug )
! This Subroutine 

  use Level_Class           ,only:  Level_Type

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                           !< Intra-nuclear diatomic potential object
  type(Level_Type)                          ,intent(inout)  ::    State 
  integer                                   ,intent(in)     ::    iState
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    h
  real(rkp)                                                 ::    e
  real(rkp) ,dimension(2)                                   ::    RrangeExtreme, Rrange0, Rrange1, Rrange2
  integer                                                   ::    ierr
  logical                                                   ::    i_Debug_Loc
  
  real(rkp)                                                 ::    roErrRel, riErrRel
  real(rkp)                                                 ::    roErrAbs, riErrAbs
  real(rkp) ,parameter                                      ::    ToleranceR = 1.d-3 !1.d-5
  real(rkp) ,parameter                                      ::    ToleranceV = 1.d-3 !1.d-8

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "CheckTurningPoints")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc)  call Logger%Write( " -> Level ", iState, "; State%Eint = ", State%Eint, "; State%vqn = ", State%vqn, "; State%jqn = ", State%jqn)
  

  Rrange1 = [State%ri*0.8d0, State%rMin]
  call This%TurningPoint( Rrange1, State%Vc_R2, State%Eint, State%riNew, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal)
  if (i_Debug_Loc)  call Logger%Write( " -> ri    = ", State%ri)
  if (i_Debug_Loc)  call Logger%Write( " -> riNew = ", State%riNew)
  
  riErrAbs = dabs( State%ri - State%riNew )
  riErrRel = dabs( riErrAbs / State%ri )  
  if ( riErrRel > ToleranceR ) then
    if (i_Debug_Loc)  call Logger%Write( " WARNING! Found a difference between the ri from Levels List and the one ricomputed. Absolute Error = ", riErrAbs, "; Relative Error = ", riErrRel*100.d0, "%")
  end if
  
  
  Rrange2 = [State%rMin, State%rMax]
  call This%TurningPoint( Rrange2, State%Vc_R2, State%Eint, State%roNew, ierr, NBisection=NBisectGlobal, NNewton=NNewtonGlobal)
  if (i_Debug_Loc)  call Logger%Write( " -> ro    = ", State%ro)
  if (i_Debug_Loc)  call Logger%Write( " -> roNew = ", State%roNew)

  roErrAbs = dabs( State%ri - State%riNew )
  roErrRel = dabs( roErrAbs / State%ri ) 
  if ( roErrRel > ToleranceR ) then
    if (i_Debug_Loc)  call Logger%Write( " WARNING! Found a difference between the ro from Levels List and the one ricomputed. Absolute Error = ", roErrAbs, "; Relative Error = ", roErrRel*100.d0, "%")
  end if
  

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine TurningPoint( This, Rrange, Vc_R2, Eint, Rt, ier, NBisection, NNewton)
! This procedure searchs for turning point, that is, one of the spatial points in which the kinetics energy is equal to zero.
! r1 is starting place of search, and r2 signals the end of the search ( in which case the search fails).
! r2 should not signficantly be beyond the turning point.
! r1 can be > r2 or < r2.
! If the number of bisection steps nbis is 0, pass the initial guess for the turning point in rt.
! If the output variable Rt is zero, then an error has occured.

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                            !< Diatomic-Potential object
  real(rkp) ,dimension(2)                   ,intent(in)     ::    Rrange                          !< Lower/upper bound for the search range for the turning point [bohr]
  real(rkp)                                 ,intent(in)     ::    Vc_R2                           !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    Eint                            !< Internal energy level
  real(rkp)                                 ,intent(out)    ::    Rt                              !< Turning point [bohr]
  integer                                   ,intent(out)    ::    ier
  integer                         ,optional ,intent(in)     ::    NBisection                      !< Number of bisection steps
  integer                         ,optional ,intent(in)     ::    NNewton                         !< Number of newton-raphson steps
 
  integer                                                   ::    NBis                            ! Number of bisection steps
  integer                                                   ::    NNew                            ! Number of newton-raphson steps
  integer                                                   ::    Iter                            ! Iteration index of the bisection/newton steps
  real(rkp)                                                 ::    Ra, Rc                          ! Internuclear distance [bohr]
  real(rkp)                                                 ::    Va, Vb, Vc, V                   ! Potentials
  real(rkp)                                                 ::    V2                              ! Potential squared
  real(rkp)                                                 ::    dVb, dV
 
  ier   =   0

  NBis  = NBisectGlobal
  NNew  = NNewtonGlobal
  if ( present(NBisection) ) NBis  = NBisection
  if ( present(NNewton   ) ) NNew  = NNewton

! ==============================================================================================================
!     PERFORMING THE BISECTION STEPS
! ==============================================================================================================
  if ( NBis > 0 ) then
    Ra      =   Rrange(1)
    Rc      =   Rrange(2)
    Va      =   This%EffectivePotential( Ra, Vc_R2 ) - Eint
    Vc      =   This%EffectivePotential( Rc, Vc_R2 ) - Eint
    V2      =   Vc * Va
    if ( ( V2 > Zero ) .and. ( V2 > EpsStatPtsGlobal ) ) then
      Rt    =   Zero
      ier   =   1                                                                                 ! <- The turning point is not contained in the R-range that has been given as input
      return
    end if
    do Iter = 1,NBis
      Rt    =   Half * (Ra+Rc)
      Vb    =   This%EffectivePotential( Rt, Vc_R2 ) - Eint
      V2    =   Vb * Va
      if ( V2 > Zero ) then
        Ra  =   Rt
        Va  =   Vb
      else
        Rc  =   Rt
        Vc  =   Vb
      end if
    end do
  end if
! ==============================================================================================================


! ==============================================================================================================
!     PERFORMING THE NEWTON-RAPHSON STEPS
! ==============================================================================================================
  do Iter = 1,NNew
    Vc = CentrifualPotential( Rt, Vc_R2 )                         ! Compiler bug with gcc-5.3 with  This%CentrifualPotential
    call This%Compute_Vd_dVd( Rt, Vb, dVb )
    V  = Vb + Vc
    dV = dVb - Two * Vc / Rt
    Rt = Rt + ( Eint - V ) / dV                           ! N.B.: dVb_e (Effective Diatomic Force) = dVb (Diatomic Force) + dVb_c (Centrifugal force)
  end do
! ==============================================================================================================

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function ActionIntegral( This, Ri, Ro, Vc_R2, Eint, n, TwoMu, QuadratureType ) result(Actionn)
! This procedure computes the action integral using a n point gauss-legendre quadrature.

  use GaussQuadrature_Module ,only: gaussq

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                            !< Diatomic-Potential object
  real(rkp)                                 ,intent(in)     ::    Ri                              !< Inner turning points [bohr]
  real(rkp)                                 ,intent(in)     ::    Ro                              !< Outer turning points
  real(rkp)                                 ,intent(in)     ::    Vc_R2                           !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    Eint                            !< internal energy
  integer                                   ,intent(in)     ::    n                               !< number of points in quadrature. If n is negitive, the weights and nodes will not be calculated but will be assumed to be passed in the array's x and w.
  real(rkp)                                 ,intent(in)     ::    TwoMu                           !< Two times the reduced mass
  integer                                   ,intent(in)     ::    QuadratureType
  real(rkp)                                                 ::    Actionn                          !< output action

  integer                                                   ::    i
  integer                                                   ::    nneg, npos
  integer                                                   ::    nu
  real(rkp)                                                 ::    Rp, Rm
  real(rkp)                                                 ::    Ve                              ! Effective potential
  real(rkp)                                                 ::    Energy
  real(rkp) ,dimension(2)                                   ::    ep
  real(rkp) ,dimension( abs(n) )                            ::    b
  real(rkp) ,dimension( abs(n) )                            ::    x
  real(rkp) ,dimension( abs(n) )                            ::    w
  real(rkp) ,dimension(2,abs(n))                            ::    tmp


  if ( n > 0 ) then
    call gaussq(QuadratureType,n,Zero,Zero,0,ep,b,x,w)      ! Weighting function (w) and evaluation points (x)
    if ( QuadratureType == 3 ) w = w / sqrt( One - x**2 )   ! Weighting function => Chebyshev–Gauss
  end if

  Rm          =   ( Ro - Ri ) * Half
  Rp          =   ( Ro + Ri ) * Half + Rm * x(1) 
  
  Actionn      =   Zero
  nu          =   abs(n)
  nneg        =   0
  npos        =   0
  
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> nu = ',g0)") nu
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Rm = ',es15.8)") Rm
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Rp = ',es15.8)") Rp
  
  do i = 1,nu
!    write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Rp     = ',es15.8)") Rp
!    write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Vc_R2  = ',es15.8)") Vc_R2
    Ve        =   This%EffectivePotential( Rp, Vc_R2 )      ! Computing the effective potential
!    write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Eint   = ',es15.8)") Eint
!    write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Ve     = ',es15.8)") Ve
    Energy    =   Eint - Ve
!    write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Energy = ',es15.8)") Energy
    if ( QuadratureType == 3 ) then
      tmp(1,i)  =   Rp
      tmp(2,i)  =   Energy
      if ( Energy > Zero ) then
        npos    =   npos + 1
      else
        nneg    =   nneg + 1
        Energy  =   - Energy
      end if
    end if
    Actionn    =   Actionn + sqrt(Energy) * w(i)
!    write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Action = ',es15.8)") Action
    if ( i /= nu ) Rp = Rp + Rm * ( x(i+1) - x(i) )
  end do
  
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Action = ',es15.8)") Action
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Rm = ',es15.8)") Rm
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> TwoMu = ',g0)") TwoMu
  
  Actionn =  Two * sqrt(TwoMu) * Actionn * Rm
  
  if ( QuadratureType == 3 ) then
    if ( (nneg > 0) .and. (npos > 0) ) then
!      write(Logger%Unit,"(2x,'[Initialize_Analyzor]: Error')")
!      write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Number of negative sqrts = ',g0)") nneg
!      write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Number of positive sqrts = ',g0)") npos
!      write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Eint  = ',es15.8)") Eint
!      write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Vc_R2 = ',es15.8)") Vc_R2
      do i = 1,nu
!        write(*,*) 'Rp = ', tmp(1,i), '   Energy = ', tmp(2,i)
!        write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Rp = ',es15.8,3x,'Energy = ',es15.8)") tmp(1:2,i)
      end do
!      stop
    end if
  end if
  
!  write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Action = ',es15.8)") Action

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function ActionIntegralGamma( This, Ri, Ro, Vc_R2, Eint, n, TwoMu, QuadratureType ) result(Action)
! This procedure computes the action integral using a n point gauss-legendre quadrature.
  
  use GaussQuadrature_Module ,only: gaussq
  
  class(DiatomicPotential_Type)             ,intent(in)     ::    This                            !< Diatomic-Potential object
  real(rkp)                                 ,intent(in)     ::    Ri                              !< Inner turning points [bohr]
  real(rkp)                                 ,intent(in)     ::    Ro                              !< Outer turning points
  real(rkp)                                 ,intent(in)     ::    Vc_R2                           !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    Eint                            !< internal energy
  integer                                   ,intent(in)     ::    n                               !< number of points in quadrature. If n is negitive, the weights and nodes will not be calculated but will be assumed to be passed in the array's x and w.
  real(rkp)                                 ,intent(in)     ::    TwoMu                           !< Two times the reduced mass
  integer                                   ,intent(in)     ::    QuadratureType
  real(rkp)                                                 ::    Action                          !< output action

  integer                                                   ::    i
  integer                                                   ::    nneg, npos
  integer                                                   ::    nu
  real(rkp)                                                 ::    Rp, Rm
  real(rkp)                                                 ::    Ve                              ! Effective potential
  real(rkp)                                                 ::    Energy
  real(rkp) ,dimension(2)                                   ::    ep
  real(rkp) ,dimension( abs(n) )                            ::    b
  real(rkp) ,dimension( abs(n) )                            ::    x
  real(rkp) ,dimension( abs(n) )                            ::    w
  real(rkp) ,dimension(2,abs(n))                            ::    tmp

  if ( n > 0 ) then
    call gaussq(QuadratureType,n,Zero,Zero,0,ep,b,x,w)      ! Weighting function (w) and evaluation points (x)
    if ( QuadratureType == 3 ) w = w / sqrt( One - x**2 )   ! Weighting function => Chebyshev–Gauss
  end if

  Rm          =   ( Ro - Ri ) * Half
  Rp          =   ( Ro + Ri ) * Half + Rm * x(1)
  Action      =   Zero
  nu          =   abs(n)
  nneg        =   0
  npos        =   0
  
  do i = 1,nu

    Ve     = This%EffectivePotential( Rp, Vc_R2 )      ! Computing the effective potential
    Energy = Ve - Eint

    if ( QuadratureType == 3 ) then
      tmp(1,i)  =   Rp
      tmp(2,i)  =   Energy
      if ( Energy > Zero ) then
        npos    =   npos + 1
      else
        nneg    =   nneg + 1
        Energy  =   - Energy
      end if
    end if
    Action      =   Action + sqrt(Energy) * w(i)

    if ( i /= nu ) Rp = Rp + Rm * ( x(i+1) - x(i) )
  end do
  
  Action = exp( - Two * sqrt(TwoMu) * Action * Rm )
  
!   if ( QuadratureType == 3 ) then
!     if ( (nneg > 0) .and. (npos > 0) ) then
!       do i = 1,nu
! !        write(*,*) 'Rp = ', tmp(1,i), '   Energy = ', tmp(2,i)
! !        write(Logger%Unit,"(2x,'[Initialize_Analyzor]: -> Rp = ',es15.8,3x,'Energy = ',es15.8)") tmp(1:2,i)
!       end do
!       stop
!     end if
!   end if
  

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine FindMinimum( This, Rrange, Vc_R2, Tolerence, Rmin, Vmin, ierr )
  
  class(DiatomicPotential_Type)             ,intent(in)     ::    This                            !< Diatomic-Potential object
  real(rkp) ,dimension(2)                   ,intent(in)     ::    Rrange                          !< Lower/upper bound for the search range of the potential min/max location [bohr]
  real(rkp)                                 ,intent(in)     ::    Vc_R2                           !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    Tolerence                       !< Absolute convergence tolerence for the search
  real(rkp)                                 ,intent(out)    ::    Rmin
  real(rkp)                                 ,intent(out)    ::    Vmin
  integer                                   ,intent(out)    ::    ierr
  
  real(rkp) ,parameter                                      ::    Sig = + One
  
  call FindExtremum( This, Sig, Rrange, Vc_R2, Tolerence, Rmin, Vmin, ierr )
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine FindMaximum( This, Rrange, Vc_R2, Tolerence, Rmax, Vmax, ierr )
  
  class(DiatomicPotential_Type)             ,intent(in)     ::    This                            !< Diatomic-Potential object
  real(rkp) ,dimension(2)                   ,intent(in)     ::    Rrange                          !< Lower/upper bound for the search range of the potential min/max location
  real(rkp)                                 ,intent(in)     ::    Vc_R2                           !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    Tolerence                       !< Absolute convergence tolerence for the search
  real(rkp)                                 ,intent(out)    ::    Rmax
  real(rkp)                                 ,intent(out)    ::    Vmax
  integer                                   ,intent(out)    ::    ierr
  
  real(rkp) ,parameter                                      ::    Sig = - One
  
  call FindExtremum( This, Sig, Rrange, Vc_R2, Tolerence, Rmax, Vmax, ierr )
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine FindExtremum( This, Sig, Rrange, Vc_R2, Tolerence, Rout, Vout, ierr )
! This procedure looks for min/max of potential function.
!  * sig
!  * re
!  * Tolerence
!  * Rmin   ! Location of min/max
!  * Vmin   ! Value of potential at Rmin
!  * ierr   ! is 0 for normal return, otherwise it is 1.
!
!  Find location of potential minimum. First find crude location minimum, then switch over to a more clever method.

  class(DiatomicPotential_Type)             ,intent(in)     ::    This                            !< Diatomic-Potential object
  real(rkp)                                 ,intent(in)     ::    Sig                             !< Indicator of the type of extremum to search for: +1=>minimum and -1=>maximum.
  real(rkp) ,dimension(2)                   ,intent(in)     ::    Rrange                          !< Lower/upper bound for the search range of the potential min/max location
  real(rkp)                                 ,intent(in)     ::    Vc_R2                           !< Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                 ,intent(in)     ::    Tolerence                       !< Absolute convergence tolerence for the search
  real(rkp)                                 ,intent(out)    ::    Rout
  real(rkp)                                 ,intent(out)    ::    Vout
  integer                                   ,intent(out)    ::    ierr

  real(rkp)                                                 ::    Ra, Rb, Rc, Ri                  ! Internuclear distances
  real(rkp)                                                 ::    Va, Vb, Vc                      ! Potential
  real(rkp)                                                 ::    Rt1
  real(rkp)                                                 ::    Vlast, Vt1, Vtry
  real(rkp)                                                 ::    Re
  real(rkp)                                                 ::    Rlim

  real(rkp) ,parameter                                      ::    Factor  = 0.3819660115_rkp

  ierr        =   0

  Re          =   Rrange(1)
  Rlim        =   Rrange(2)

  Ra          =   Re - RStepStatPtsGlobal
  Va          =   This%EffectivePotential( Ra, Vc_R2 ) * Sig
  Vlast       =   This%EffectivePotential( Re, Vc_R2 ) * Sig

  Ri          =   Re
  do
    Ri        =   Ri + RStepStatPtsGlobal
    if ( Ri >= Rlim ) then
      Rout    =   Ri - RStepStatPtsGlobal
      Vout    =   Vlast * Sig
      ierr    =   1
      return
    end if
    Vtry      =   This%EffectivePotential( Ri, Vc_R2 ) * Sig
    if ( Vtry > Vlast ) exit
    Vlast     =   Vtry
  end do


  Rb          =   Ri - RStepStatPtsGlobal
  Vb          =   Vlast
  Rc          =   Ri
  Vc          =   Vtry

  Ra          =   Rb - RStepStatPtsGlobal
  do
    if ( Rb-Ra > Rc-Rb ) then
      rt1     =   (Ra-Rb) * Factor + Rb
      vt1     =   This%EffectivePotential( rt1, Vc_R2 ) * Sig
! ! ! ! ! ! ! !       if ( abs(rt1) < 1d-3) stop 5
      if ( Vb < vt1 ) then
        Ra    =   rt1
      else
        Rc    =   Rb
        Vc    =   Vb
        Rb    =   rt1
        Vb    =   vt1
      end if
    else
      rt1     =   (Rc-Rb) * Factor + Rb
      vt1     =   This%EffectivePotential( rt1, Vc_R2 ) * Sig
      if ( Vb < vt1 ) then
        Rc    =   rt1
        Vc    =   vt1
      else
        Ra    =   Rb
        Va    =   Vb
        Rb    =   rt1
        Vb    =   vt1
      end if
    end if

    if ( abs(Rc-Ra) <= Tolerence ) exit

  end do

  Rout      =     Rb
  Vout      =     Vb * Sig

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module