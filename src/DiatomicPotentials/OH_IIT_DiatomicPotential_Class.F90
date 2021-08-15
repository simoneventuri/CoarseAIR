! -*-F90-*-
!===============================================================================================================
!
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR)
!
! Copyright (C) 2018 Berkan Bolkan and Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign).
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
! This module computes the diatomic potential of ground state OH. The energy values were computed and
! reported on "Ab initio potential energy curves for the ground and low-lying excited states of OH and OH- and
! a study of rotational fine structure in photodetachment" by Srivastava, S. and Sathyamurthy, N. in the
! Journal of Physical Chemistry A 118, 2014. The potential did not have an analytical form but was fitted
! to the following expression with fixed equilibrium distance and dissociation energy reported in NIST:
! V(r) = De*( 1 - exp( -f*(r-re) ) )ˆ2
! f    = a + b*rd + c*rdˆ2 + d*rdˆ3 + e*rdˆ4 + f*rdˆ5 + g*rdˆ6
! rd   = (rˆ4 - reˆ4)/(rˆ4 + reˆ4)
! with re = 0.96966 Angstrom and
! De = 4.63 eV from Arnold, J.O.; Whiting, E.E.; Sharbaugh, L.F., J. Chem. Phys., 1976, 64, 3251
! The coefficient values obtained from the fit are:
! a =  1.294987373533374
! b = -0.12116470501146599
! c = -0.7814954684504813
! d =  0.4582199547136645
! e =  2.2001736402873724
! f = -0.36592711896172825
! g = -1.4691984743053073
! Module implemented by João Vargas (joao.francisco.vargas@gmail.com)
!===============================================================================================================

Module OH_IIT_DiatomicPotential_Class

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Eight, B_To_Ang, eV_To_Hartree
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  Private
  public  ::    OH_IIT_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    OH_IIT_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_OH_IIT_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_OH
    procedure         ::    DiatomicPotential =>    DiatomicPotential_OH
  End Type

  character(*)  ,parameter  ::    Name_DiatPot    = 'IIT'
  logical       ,parameter  ::    i_Debug_Global = .False.

  real(rkp) ,dimension(0:6)                 ,parameter      ::    cs   = [  1.294987373533374_rkp, -0.12116470501146599_rkp, -0.7814954684504813_rkp, 0.4582199547136645_rkp, 2.2001736402873724_rkp, -0.36592711896172825_rkp, -1.4691984743053073_rkp ]
  real(rkp)                                 ,parameter      ::    de   = 4.63_rkp
  real(rkp)                                 ,parameter      ::    red  = 0.96966_rkp
  real(rkp)                                 ,parameter      ::    red4 = red**(4.0_rkp)
  real(rkp)                                 ,parameter      ::    VRef = Zero

  contains
!--------------------------------------------------------------------------------------------------------------------------------!
Subroutine Initialize_OH_IIT_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(OH_IIT_DiatomicPotential_Type)  ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug

  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_OH_IIT_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()

  allocate( This%Name        ,source = trim(Name_DiatPot) )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =    iMol
  This%Initialized  =   .True.

  This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------------------------------------!
Elemental Function DiatomicPotential_OH( This, R ) result( V )

  class(OH_IIT_DiatomicPotential_Type)        ,intent(in)  :: This
  real(rkp)                                   ,intent(in)  :: R
  real(rkp)                                                :: V

  real(rkp)                                                :: RAng
  real(rkp)                                                :: VDiat

  RAng = R * B_To_Ang

  call Ev2gm2(RAng, VDiat)

  V = VDiat * eV_To_Hartree + VRef

End Function
!--------------------------------------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------------------------------------!
Elemental Subroutine Compute_Vd_dVd_OH( This, R, V, dV )

  class(OH_IIT_DiatomicPotential_Type)        ,intent(in)  :: This
  real(rkp)                                   ,intent(in)  :: R
  real(rkp)                                   ,intent(out) :: V
  real(rkp)                                   ,intent(out) :: dV

  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat
  real(rkp)                                    :: dVDiat

  RAng = R * B_To_Ang

  call Ev2gm2_Grad(RAng, VDiat, dVDiat)

  V  =  VDiat * eV_To_Hartree + VRef
  dV = dVDiat * eV_To_Hartree * B_To_Ang

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!
!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2( RAng, V)

  use Parameters_Module     ,only:  Zero, Two

  real(rkp)                                 ,intent(in)     ::    RAng                                                            ! Distances of atom-atom pairs [Angstrom]
  real(rkp)                                 ,intent(out)    ::    V

  real(rkp)                                                 ::    y, y2, y3, y4, y5, y6
  real(rkp)                                                 ::    RAng3, RAng4
  real(rkp)                                                 ::    TempSum, u, minu, fy, dfdy, dydr, dfdr

  RAng3   = RAng**3;
  RAng4   = RAng3*RAng;
  TempSum = (RAng4 + red4);
  y       = (RAng4 - red4) / TempSum;
  y2      = y**2;
  y3      = y2*y;
  y4      = y2*y2;
  y5      = y3*y2;
  y6      = y3*y3;

  fy      =   cs(0) + cs(1)*y + cs(2)*y2 + cs(3)*y3 + cs(4)*y4 + cs(5)*y5 + cs(6)*y6
  u       =   exp(-fy * (RAng-red))
  minu    =   One - u

  V       =   de * minu**Two - de

End Subroutine Ev2gm2
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2_Grad( RAng, V, dV)

  use Parameters_Module     ,only:  Zero, Two

  real(rkp)                                 ,intent(in)     ::    RAng              ! Distances of atom-atom pairs [Angstrom]
  real(rkp)                                 ,intent(out)    ::    V
  real(rkp)                                 ,intent(out)    ::    dV

  real(rkp)                                                 ::    y, y2, y3, y4, y5, y6
  real(rkp)                                                 ::    RAng3, RAng4
  real(rkp)                                                 ::    TempSum, u, minu, fy, dfdy, dydr, dfdr

  RAng3   = RAng**3;
  RAng4   = RAng3*RAng;
  TempSum = (RAng4 + red4);
  y       = (RAng4 - red4) / TempSum;
  y2      = y**2;
  y3      = y2*y;
  y4      = y2*y2;
  y5      = y3*y2;
  y6      = y3*y3;

  fy      =   cs(0) + cs(1)*y + cs(2)*y2 + cs(3)*y3 + cs(4)*y4 + cs(5)*y5 + cs(6)*y6
  u       =   exp(-fy * (RAng-red))
  minu    =   One - u

  dfdy    =   cs(1) + Two*cs(2)*y + Three*cs(3)*y2 + Four*cs(4)*y3 + Five*cs(5)*y4 + Six*cs(6)*y5

  dydr    =   Eight * RAng3 * red4 / TempSum**2
  dfdr    =   dfdy * dydr

  V       =   de * minu**Two - de
  dV      =   Two * de * minu * u * (dfdr * (RAng-red) + fy)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module




