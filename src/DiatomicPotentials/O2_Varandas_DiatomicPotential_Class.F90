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

Module O2_Varandas_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    O2_Varandas_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    O2_Varandas_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_O2_Varandas_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_O2_Varandas
    procedure         ::    DiatomicPotential =>    DiatomicPotential_O2_Varandas
  End Type

  character(*)  ,parameter  ::    Name_DiaPot    = 'O2_Varandas'
  logical       ,parameter  ::    i_Debug_Global = .False.
  
  real(rkp)                            ,parameter              :: D     = -0.142912d0
  real(rkp)                            ,parameter              :: Gamma =  3.3522d0
  real(rkp)                            ,parameter              :: a1    =  3.644590d0
  real(rkp)                            ,parameter              :: a2    =  3.928120d0
  real(rkp)                            ,parameter              :: a3    =  2.0986665d0
  real(rkp)                            ,parameter              :: c6    =  15.4d0 
  real(rkp)                            ,parameter              :: c8    =  235.2d0
  real(rkp)                            ,parameter              :: c10   =  4066.d0
  real(rkp)                            ,parameter              :: Rm    =  2.2818d0
  real(rkp)                            ,parameter              :: R0    =  5.66169d0
  real(rkp)                            ,parameter              :: alpha0 = 25.952776d0
  real(rkp)                            ,parameter              :: alpha1 = 1.186793d0
  real(rkp)                            ,parameter              :: beta0  = 15.738100d0
  real(rkp)                            ,parameter              :: beta1  = 0.097291d0
  
  real(rkp)                                                    :: rho    = 8.2180125
  real(rkp) ,dimension(3)                                      :: A      = [3.095133326007682, 2.199900018659006, 1.688071424840357]
  real(rkp) ,dimension(3)                                      :: B      = [8.778789466267723, 7.226512279520406, 5.948710801953137]
  !A(1) = alpha0 *  6.0d0**(-alpha1)
  !A(2) = alpha0 *  8.0d0**(-alpha1)
  !A(3) = alpha0 * 10.0d0**(-alpha1)
  !B(1) = beta0 * exp(-beta1 *  6.0d0)
  !B(2) = beta0 * exp(-beta1 *  8.0d0)
  !B(3) = beta0 * exp(-beta1 * 10.0d0)
  !rho  = (Rm + 2.5d0 * R0) / Two
 
  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_O2_Varandas_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(O2_Varandas_DiatomicPotential_Type)  ,intent(out)  ::    This
  type(Input_Type)                           ,intent(in)   ::    Input
  character(:) ,allocatable                  ,intent(in)   ::    SpeciesName
  integer                                    ,intent(in)   ::    iMol
  real(rkp)                                  ,intent(in)   ::    Mass1
  real(rkp)                                  ,intent(in)   ::    Mass2
  logical ,optional                          ,intent(in)   ::    i_Debug
  
  logical                                                  ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O2_Varandas_DiatomicPotential" )
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
Elemental Function DiatomicPotential_O2_Varandas( This, R ) result( V )

  class(O2_Varandas_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                                  ,intent(in)  :: R
  real(rkp)                                               :: V
    
  real(rkp)                                    :: Rd
  real(rkp)                                    :: V_EHF, V_Corr
  real(rkp) ,dimension(3)                      :: Xn
  
  Rd     = R - Rm;
  
  V_EHF  = D * (One + a1 * Rd + a2 * Rd**2 + a3 * Rd**3) * exp(-Gamma * Rd)
  
  Xn(1) = (One - exp(- A(1) * R / rho - B(1) * (R / rho)**2) )**6
  Xn(2) = (One - exp(- A(2) * R / rho - B(2) * (R / rho)**2) )**8
  Xn(3) = (One - exp(- A(3) * R / rho - B(3) * (R / rho)**2) )**10
  
  V_Corr = - (C6 * Xn(1) * R**(-6) + C8 * Xn(2) * R**(-8) + C10 * Xn(3) * R**(-10))
  
  V = V_EHF + V_Corr
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!  


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_O2_Varandas( This, R, V, dV )

  class(O2_Varandas_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                                  ,intent(in)  :: R
  real(rkp)                                  ,intent(out) :: V
  real(rkp)                                  ,intent(out) :: dV
  
  real(rkp)                                    :: Rd
  real(rkp)                                    :: V_EHF, V_Corr
  real(rkp) ,dimension(3)                      :: Xn
  
  Rd     = R - Rm;
  
  V_EHF  = D * (One + a1 * Rd + a2 * Rd**2 + a3 * Rd**3) * exp(-Gamma * Rd)
  
  Xn(1) = (One - exp(- A(1) * R / rho - B(1) * (R / rho)**2) )**6
  Xn(2) = (One - exp(- A(2) * R / rho - B(2) * (R / rho)**2) )**8
  Xn(3) = (One - exp(- A(3) * R / rho - B(3) * (R / rho)**2) )**10
  
  V_Corr = - (C6 * Xn(1) * R**(-6) + C8 * Xn(2) * R**(-8) + C10 * Xn(3) * R**(-10))
  
  V = V_EHF + V_Corr
  
  dV = 0.0
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
