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

Module O2_NASA_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    O2_NASA_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    O2_NASA_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_O2_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_O2
    procedure         ::    DiatomicPotential =>    DiatomicPotential_O2
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiaPot   = 'O2_NASA'

  integer                         ,parameter   :: nfitco  = 6
  real(rkp) ,dimension(nfitco+1)  ,parameter   :: coefo2  = (/  -187.888473490235_rkp, -47.6881525811959_rkp, &
                                                                 13.8451858338530_rkp, -42.2149670400162_rkp, &
                                                                 15.9948962333531_rkp, -35.3271903135257_rkp, &  
                                                                 32.0805760535816_rkp /)
  
  real(rkp)                       ,parameter   :: c6o2    = 97.3185530097047_rkp
  real(rkp)                       ,parameter   :: damp6co = 418.892_rkp
  real(rkp)                       ,parameter   :: ao2     = 5.68_rkp
  real(rkp)                       ,parameter   :: R0      = 2.15_rkp
  real(rkp)                       ,parameter   :: alphaco = -2.234_rkp
  integer                         ,parameter   :: npowco  = 9

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_O2_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(O2_NASA_DiatomicPotential_Type) ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O2_DiatomicPotential" )
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
Elemental Subroutine Compute_Vd_dVd_O2( This, R, V, dV )

  class(O2_NASA_DiatomicPotential_Type) ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                             ,intent(out) :: V
  real(rkp)                             ,intent(out) :: dV
  
  real(rkp)                                    :: vsr,  dvsr
  real(rkp)                                    :: vlr,  dvlr
  real(rkp)                                    :: fact, dfact
  real(rkp)                                    :: cfcn, dcfcn
  integer                                      :: k

  vsr   = 64.0_rkp * dexp(alphaco*R) / R
  dvsr  = vsr * (alphaco-(One/R))
  vlr   = -c6o2 / (R**6+damp6co)
  dvlr  = -6.0_rkp * (R**5) * vlr / (R**6+damp6co)
  fact  = (R**npowco) * dexp(-ao2*R)
  dfact = fact * ((dfloat(npowco)/R)-ao2)
  
  cfcn  = coefo2(nfitco+1)
  dcfcn = Zero
  do k = nfitco,2,-1
   dcfcn = cfcn      + dcfcn * (R-R0)
   cfcn  = coefo2(k) + cfcn  * (R-R0)
  end do
  dcfcn = dfact     * cfcn + fact * dcfcn
  cfcn  = coefo2(1) + fact * cfcn
  
  v  = cfcn  + vsr  + vlr  + 187.888473490235_rkp
  dV = dcfcn + dvlr + dvsr

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_O2( This, R ) result( V )

  class(O2_NASA_DiatomicPotential_Type) ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                                          :: V
  
  real(rkp)                                    :: vsr
  real(rkp)                                    :: vlr
  real(rkp)                                    :: fact
  real(rkp)                                    :: cfcn
  integer                                      :: k

  vsr   = 64.0_rkp * dexp(alphaco*R) / R
  vlr   = -c6o2 / (R**6+damp6co)
  fact  = (R**npowco) * dexp(-ao2*R)
  
  cfcn  = coefo2(nfitco+1)
  do k = nfitco,2,-1
   cfcn  = coefo2(k) + cfcn  * (R-R0)
  end do
  cfcn  = coefo2(1) + fact * cfcn
  
  v  = cfcn  + vsr  + vlr  + 187.888473490235_rkp
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
