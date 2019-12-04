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

Module CO_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Six
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    CO_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    CO_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_CO_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_CO
    procedure         ::    DiatomicPotential =>    DiatomicPotential_CO
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiaPot = 'NASA'

  integer                   ,parameter  :: nf     = 12
  real(rkp) ,dimension(nf) ,parameter   :: coef   = (/  -4.5761074811943886e-2_rkp,    0.22111734303145603_rkp, &
                                                          -0.34188485707539029_rkp,    0.21975372852913833_rkp, &
                                                          -0.18415829369278044_rkp,    0.19162121009025507_rkp, &
                                                          -0.14681584359973315_rkp,  7.3117032359438006e-2_rkp, & 
                                                        -2.2588593079026137e-2_rkp,  4.1361679444378470e-3_rkp, &
                                                        -4.0747646620877898e-4_rkp,  1.6857235957160763e-5_rkp /)
  real(rkp)                 ,parameter  :: c6     =  284.57619366543872_rkp
  real(rkp)                 ,parameter  :: cinf   = -112.90087267991528_rkp
  real(rkp)                 ,parameter  :: alpha  = -2.2343545089620513_rkp
  real(rkp)                 ,parameter  :: damp6  =  418.89225320735858_rkp
  integer                   ,parameter  :: npow   =  9
  real(rkp)                 ,parameter  :: a      =  4.18604651162790730_rkp
  real(rkp)                 ,parameter  :: frel   =  0.12114229643464522_rkp
  real(rkp)                 ,parameter  :: R0     =  2.14999999999999999_rkp

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_CO_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(CO_DiatomicPotential_Type)      ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_CO_DiatomicPotential" )
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
Elemental Subroutine Compute_Vd_dVd_CO( This, R, V, dV )

  class(CO_DiatomicPotential_Type)      ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                             ,intent(out) :: V
  real(rkp)                             ,intent(out) :: dV
 
  integer                               :: nfm
  real(rkp)                             :: freli
  real(rkp)                             :: vsr, dvsr
  real(rkp)                             :: vlr, dvlr
  real(rkp)                             :: vint, dvint
  
  integer                               :: j
 
  nfm   = nf - 1
  freli = One / frel
  
  vsr  = 48.0_rkp * dexp(alpha*R) / R
  dvsr = (alpha - (One / R)) * vsr
  vlr  = -c6 / (R**6 + damp6)
  dvlr = Six * (R**5) * (vlr*vlr) / c6
  
  vint  = coef(nf)
  dvint = Zero
  do j = nfm,1,-1
    dvint = vint    + dvint * (R-R0)
    vint  = coef(j) + vint  * (R-R0)
  end do
  dvint = dvint * (R**npow) * dexp(-a*R) * freli
  vint  = vint  * (R**npow) * dexp(-a*R) * freli
  dvint = dvint + ((dfloat(npow)/R)-a) * vint
  
  V  = vsr  + vlr  + vint
  dV = dvsr + dvlr + dvint

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_CO( This, R ) result( V )

  class(CO_DiatomicPotential_Type)      ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                                          :: V
 
  integer                               :: nfm
  real(rkp)                             :: freli
  real(rkp)                             :: vsr
  real(rkp)                             :: vlr
  real(rkp)                             :: vint
  integer                               :: j
 
  nfm   = nf - 1
  freli = One / frel
  
  vsr  = 48.0_rkp * dexp(alpha*R) / R
  vlr  = -c6 / (R**6 + damp6)
  
  vint  = coef(nf)
  do j = nfm,1,-1
    vint  = coef(j) + vint  * (R-R0)
  end do
  vint  = vint  * (R**npow) * dexp(-a*R) * freli
  
  V  = vsr  + vlr  + vint
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
