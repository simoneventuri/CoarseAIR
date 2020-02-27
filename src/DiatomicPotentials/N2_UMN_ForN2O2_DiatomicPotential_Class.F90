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

Module N2_UMN_ForN2O2_DiatomicPotential_Class

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  Private
  public  ::    N2_UMN_ForN2O2_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    N2_UMN_ForN2O2_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_N2_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_N2
    procedure         ::    DiatomicPotential =>    DiatomicPotential_N2
  End Type

  character(*)  ,parameter  ::    Name_DiaPot    = 'UMN_ForN2O2'
  logical       ,parameter  ::    i_Debug_Global = .False.
  
  !!! For N2+O
  real(rkp) ,dimension(0:6)                 ,parameter      ::    cs   = [  2.71405774451e0_rkp, 1.32757649829e-1_rkp, 2.66756890408e-1_rkp, 1.95350725241e-1_rkp, -4.08663480982e-1_rkp, 3.92451705557e-1_rkp, 1.13006674877e0_rkp ]
  real(rkp)                                 ,parameter      ::    de   = 228.4_rkp 
  real(rkp)                                 ,parameter      ::    red  = 1.098_rkp
  real(rkp)                                 ,parameter      ::    red4 = 1.4534810048_rkp
  real(rkp)                                 ,parameter      ::    VRef = Zero!0.191619504727d0 !Reference energy of infinitely separated N2 + O in hartree (taken from UMN DSEC corrected calculations)

  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_N2_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(N2_UMN_ForN2O2_DiatomicPotential_Type)  ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N2_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = trim(Name_DiaPot) )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =    iMol
  This%Initialized  =   .True.

  This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_N2( This, R ) result( V )

  class(N2_UMN_ForN2O2_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                                     ,intent(in)  :: R
  real(rkp)                                                  :: V
    
  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat
  
  RAng = R * B_To_Ang                                                                                                                   

  call Ev2gm2(RAng, VDiat)
  
  V = VDiat * Kcm_To_Hartree + VRef
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!  

!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_N2( This, R, V, dV )

  class(N2_UMN_ForN2O2_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                                     ,intent(in)  :: R
  real(rkp)                                     ,intent(out) :: V
  real(rkp)                                     ,intent(out) :: dV
  
  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat
  real(rkp)                                    :: dVDiat
  
  RAng = R * B_To_Ang
  
  call Ev2gm2_Grad(RAng, VDiat, dVDiat)
  
  V  =  VDiat * Kcm_To_Hartree + VRef                                                                       
  dV = dVDiat * Kcm_To_Hartree * B_To_Ang
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2( RAng, V)

!**********************************************************************
!
! This subroutine evaulates the 2-body potential energy and gradient
! for given r with two kinds of potentials:
!
! A polinomial based on MEG variables and a Lennard-Jones term
! V(r) = 0 for r -> inifity
! u=exp(-(r-red)/alpha-(r-red)^2/beta))
! V(r)=-De*(cs1*u^1+cs2*u^2+cs3*u^3+cs4*u^4+cs5*u^5+cs6*u^6)
!
!**********************************************************************

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

  V       =   de * minu**2 - de
  
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
  minu    =   1.d0 - u

  dfdy    =   cs(1) + 2.0d0*cs(2)*y + 3.0d0*cs(3)*y2 + 4.0d0*cs(4)*y3 + 5.0d0*cs(5)*y4 + 6.0d0*cs(6)*y5

  dydr    =   8.0d0 * RAng3 * red4 / TempSum**2
  dfdr    =   dfdy * dydr

  V       =   de * minu**2 - de
  dV      =   2.0d0 * de * minu * u * (dfdr * (RAng-red) + fy)
   
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
