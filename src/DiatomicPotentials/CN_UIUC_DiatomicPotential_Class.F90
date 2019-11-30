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
! You should have received a Copy of the GNU Lesser General Public License along with this library; 
! if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
! 
!---------------------------------------------------------------------------------------------------------------
!===============================================================================================================

Module CN_UIUC_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Six
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    CN_UIUC_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    CN_UIUC_DiatomicPotential_Type
  Contains
    procedure         ::    Initialize        =>    Initialize_CN_UIUC_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_CN_UIUC
    procedure         ::    DiatomicPotential =>    DiatomicPotential_CN_UIUC
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiaPot = 'CN_UIUC'

  integer                        ,parameter  :: PolOrd = 14
  real(rkp)                      ,parameter  :: re     = 2.022652218205421_rkp
  real(rkp)                      ,parameter  :: Beta   = 0.688711263655651_rkp
  real(rkp) ,dimension(PolOrd+1) ,parameter  :: cPol   = (/  0.000000000000000e+04_rkp, &
                                                             0.000003674932525e+04_rkp, &
                                                            -0.000191729648640e+04_rkp, &
                                                             0.002270200689903e+04_rkp, &
                                                            -0.008509563620165e+04_rkp, &
                                                            -0.011879551093198e+04_rkp, &
                                                             0.171199809102443e+04_rkp, &
                                                            -0.574183242420719e+04_rkp, &
                                                             1.076964416909275e+04_rkp, &
                                                            -1.297336272731918e+04_rkp, &
                                                             1.045262673602783e+04_rkp, &
                                                            -0.562525357308539e+04_rkp, &
                                                             0.194604438421548e+04_rkp, &
                                                            -0.039203917425613e+04_rkp, &
                                                             0.003499620526327e+04_rkp /)

  Contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_CN_UIUC_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(CN_UIUC_DiatomicPotential_Type) ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_CN_UIUC_DiatomicPotential" )
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
Elemental Subroutine Compute_Vd_dVd_CN_UIUC( This, R, V, dV )
  
  class(CN_UIUC_DiatomicPotential_Type) ,intent(in)  :: This
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
Elemental Function DiatomicPotential_CN_UIUC( This, R ) result( V )

  class(CN_UIUC_DiatomicPotential_Type) ,intent(in)  :: This
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
