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

Module CH_UIUC_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Six
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    CH_UIUC_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    CH_UIUC_DiatomicPotential_Type
  Contains
    procedure         ::    Initialize        =>    Initialize_CH_UIUC_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_CH_UIUC
    procedure         ::    DiatomicPotential =>    DiatomicPotential_CH_UIUC
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiaPot = 'UIUC'

  integer                        ,parameter  :: PolOrd =  20 !13
  real(rkp)                      ,parameter  :: re     =  1.301179737758152_rkp !1.498759159484161_rkp
  real(rkp)                      ,parameter  :: Beta   =  1.044155249442320_rkp !1.765060790703594_rkp
  real(rkp) ,dimension(PolOrd+1) ,parameter  :: cPol   =  (/ 0.0_rkp, 0.00379971564237380_rkp, -0.953712712168401_rkp, -3.59328529208477_rkp, 37.8323517379414_rkp, -124.786301375258_rkp, 226.722421321307_rkp, -218.909452090161_rkp, 48.2886301126978_rkp, 107.010402227210_rkp, -68.7353658518817_rkp, -29.4476842865390_rkp, -9.83386396216202_rkp, 79.7356679254129_rkp, -3.82511534297961_rkp, -132.925148097350_rkp, 165.220696460728_rkp, -100.141985282533_rkp, 34.5002667646521_rkp, -6.50481292023681_rkp, 0.524699734678209 /)
                                                          !(/  0.000000000000000e+00_rkp, &
                                                            ! -0.004030637922008e+02_rkp, &
                                                            !  0.003370387345772e+02_rkp, &
                                                            !  0.031210072284260e+02_rkp, &
                                                            ! -0.193537545390564e+02_rkp, &
                                                            !  0.773889740488166e+02_rkp, &
                                                            ! -2.150805024065447e+02_rkp, &
                                                            !  4.174683339555139e+02_rkp, &
                                                            ! -5.702374626716050e+02_rkp, &
                                                            !  5.457748545456931e+02_rkp, &
                                                            ! -3.577278681230736e+02_rkp, &
                                                            !  1.527254028673570e+02_rkp, &
                                                            ! -0.382338403243996e+02_rkp, &
                                                            !  0.042584888475653e+02_rkp /)

  Contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_CH_UIUC_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(CH_UIUC_DiatomicPotential_Type) ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_CH_UIUC_DiatomicPotential" )
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
Elemental Subroutine Compute_Vd_dVd_CH_UIUC( This, R, V, dV )

  class(CH_UIUC_DiatomicPotential_Type) ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                             ,intent(out) :: V
  real(rkp)                             ,intent(out) :: dV
 
  integer                               :: i
  real(rkp)                             :: Xpnt

  Xpnt = exp( - Beta * (R - re) )

  V  = cPol(1)
  dV = cPol(2)
  do i=1,PolOrd-1
    V  = V  +         cPol(i+1) * Xpnt**i
    dV = dV + (i+1) * cPol(i+2) * Xpnt**i
  end do
  V  = V  + cPol(PolOrd+1) * Xpnt**PolOrd
  dV = dV * (- Beta * Xpnt)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_CH_UIUC( This, R ) result( V )

  class(CH_UIUC_DiatomicPotential_Type) ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                                          :: V
 
  integer                               :: i
  real(rkp)                             :: Xpnt

  Xpnt = exp( - Beta * (R - re) )

  V  = cPol(1)
  do i=1,PolOrd
    V  = V + cPol(i+1) * Xpnt**i
  end do
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
