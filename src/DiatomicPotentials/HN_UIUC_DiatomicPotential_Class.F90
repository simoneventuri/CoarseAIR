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

Module HN_UIUC_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Six
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    HN_UIUC_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    HN_UIUC_DiatomicPotential_Type
  Contains
    procedure ::    Initialize        =>    Initialize_HN_UIUC_DiatomicPotential
    procedure ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_HN_UIUC
    procedure ::    DiatomicPotential =>    DiatomicPotential_HN_UIUC
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiaPot = 'UIUC'

  integer                        ,parameter  :: PolOrd = 20 !13
  real(rkp)                      ,parameter  :: re     = 1.434108289721189_rkp !1.241140761868844_rkp
  real(rkp)                      ,parameter  :: Beta   = 0.874341879007108_rkp !0.712073604463818_rkp
  real(rkp) ,dimension(PolOrd+1) ,parameter  :: cPol   = (/ 0.0_rkp, 0.0148065083264408_rkp, -1.51689356091355_rkp, 2.58050537015949_rkp, -0.605095080945340_rkp, -1.21820066050174_rkp, 0.0227002981216103_rkp, 1.00822355421032_rkp, -0.127698709570522_rkp, 0.439980532910101_rkp, -0.823690344520123_rkp, -0.314537373409736_rkp, 0.678543670187740_rkp, -0.0986245407910460_rkp, -0.0690968172738827_rkp, 0.159938800081593_rkp, -0.359237573116952_rkp, 0.304260671802539_rkp, -0.115358008199976_rkp, 0.0190834162424563_rkp, -0.000935340646967855 /) 

                                                          ! (/  0.000000000000000e+00_rkp, &
                                                          !    0.000036749325248e+03_rkp, &
                                                          !   -0.001260798447708e+03_rkp, &
                                                          !    0.005933292934428e+03_rkp, &
                                                          !   -0.055299209746028e+03_rkp, &
                                                          !    0.328109229946940e+03_rkp, &
                                                          !   -1.195310344676442e+03_rkp, &
                                                          !    2.899290647273216e+03_rkp, &
                                                          !   -4.841694257967588e+03_rkp, &
                                                          !    5.604600957770337e+03_rkp, &
                                                          !   -4.425577141363808e+03_rkp, &
                                                          !    2.276637422748264e+03_rkp, &
                                                          !   -0.687960511117803e+03_rkp, &
                                                          !    0.092629844574817e+03_rkp /)

  Contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_HN_UIUC_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(HN_UIUC_DiatomicPotential_Type) ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_HN_UIUC_DiatomicPotential" )
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
Elemental Subroutine Compute_Vd_dVd_HN_UIUC( This, R, V, dV )

  class(HN_UIUC_DiatomicPotential_Type) ,intent(in)  :: This
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
Elemental Function DiatomicPotential_HN_UIUC( This, R ) result( V )

  class(HN_UIUC_DiatomicPotential_Type) ,intent(in)  :: This
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
