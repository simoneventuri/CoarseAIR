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

Module N2_LeRoy_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Half, One
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    N2_LeRoy_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    N2_LeRoy_DiatomicPotential_Type
  contains
    procedure               ::    Initialize        =>    Initialize_N2_LeRoy_DiatomicPotential
    procedure               ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_N2_LeRoy
    procedure               ::    DiatomicPotential =>    DiatomicPotential_N2_LeRoy
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiatPot = 'LeRoy'

  integer   ,parameter                                      ::    nrep  =  12
  real(rkp) ,parameter                                      ::    co    =  49.0_rkp
  real(rkp) ,parameter                                      ::    ao    =  2.4_rkp
  real(rkp) ,parameter                                      ::    a1    =  5.38484562702061_rkp
  real(rkp) ,parameter                                      ::    c     =  23.9717901220746_rkp
  real(rkp) ,parameter                                      ::    cx20  =  c*20.0_rkp
  real(rkp) ,parameter                                      ::    cx500 =  c*500.0_rkp
  real(rkp) ,parameter                                      ::    d     =  3.0_rkp
  real(rkp) ,parameter                                      ::    d2    =  d*d
  real(rkp) ,parameter                                      ::    d4    =  d2*d2
  real(rkp) ,parameter                                      ::    d6    =  d4*d2
  real(rkp) ,parameter                                      ::    re    =  2.1_rkp
  real(rkp) ,parameter  ,dimension(nrep)                    ::    csr   = [-3.89327515161896e+2_rkp,-6.19410598346194e+2_rkp,-5.51461129947346e+2_rkp, &
                                                                           -3.54896837660797e+2_rkp,-1.08347448451266e+2_rkp,-6.59348244094835e+1_rkp, &
                                                                            1.30457802135760e+1_rkp,7.20557758084161e+1_rkp,-1.81062803146583e+1_rkp, &
                                                                           -2.84137057950277e+1_rkp,1.40509401686544e+1_rkp,-1.84651681798865_rkp]

  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_N2_LeRoy_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )
  
  use Input_Class   ,only:  Input_Type
  
  class(N2_LeRoy_DiatomicPotential_Type)      ,intent(out)  ::    This
  type(Input_Type)                            ,intent(in)   ::    Input
  character(:)                   ,allocatable ,intent(in)   ::    SpeciesName
  integer                                     ,intent(in)   ::    iMol
  real(rkp)                                   ,intent(in)   ::    Mass1
  real(rkp)                                   ,intent(in)   ::    Mass2
  logical ,optional                           ,intent(in)   ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N2_LeRoy_DiatomicPotential" )
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


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_N2_LeRoy( This, R, V, dV )

  use Parameters_Module   ,only:  Zero, One

  class(N2_LeRoy_DiatomicPotential_Type)    ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Distances between nuclear centers [bohr]
  real(rkp)                                 ,intent(out)    ::    V                         ! Potential energy [hartree]
  real(rkp)                                 ,intent(out)    ::    dV                        ! First derivative of the potential energy wrt the distance [hartree/bohr]

  integer                                                   ::    n
  real(rkp)                                                 ::    dr, r2, r3, r4, r5, r6, ri
  real(rkp)                                                 ::    Vdlr, Vrep, drep, suma, dsuma, ex, Vsr

! Computing r-related terms
  ri      =   One / r
  r2      =   r   * r
  r3      =   r2  * r
  r4      =   r3  * r
  r5      =   r4  * r
  r6      =   r5  * r

  Vrep    =   co * exp(-ao*r) * ri
  drep    =   - Vrep * ( ao + ri )

  suma    =   csr(nrep)
  dsuma   =   Zero
  dr      =   r - re

  do n = nrep-1,1,-1
    dsuma =   suma + dr * dsuma
    suma  =   csr(n) + suma * dr
  end do

  ex      =   exp( - a1 * r ) * r6
  Vsr     =   suma * ex

  drep    =   drep + ( (-a1+(6.0_rkp*ri)) * suma + dsuma ) * ex

  Vdlr    =        - c     / ( r6 + d6 )
  Vdlr    =   Vdlr - cx20  / ( r4 + d4 )**2
  Vdlr    =   Vdlr - cx500 / ( r2 + d2 )**5

  dV      =   drep + 6.0_rkp  * c     * r5 / ( r6 + d6 )**2
  dV      =   dV   + 8.0_rkp  * cx20  * r3 / ( r4 + d4 )**3
  dV      =   dV   + 10.0_rkp * cx500 * r  / ( r2 + d2 )**6

  V       =   Vrep + Vdlr + Vsr

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_N2_LeRoy( This, R ) result( V )

  class(N2_LeRoy_DiatomicPotential_Type)    ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Internuclear distance [bohr]
  real(rkp)                                                 ::    V                         ! Diatomic potential energy [hartree]

  integer                                                   ::    i
  real(rkp)                                                 ::    r2, r4, r6, dr
  real(rkp)                                                 ::    Vdlr, Vrep, suma, Vsr
  
  r2      =   r   * r
  r4      =   r2  * r2
  r6      =   r4  * r2
  suma    =   csr(nrep)
  dr      =   r - re
  
  do i = nrep-1,1,-1
    suma  =   csr(i) + suma * dr
  end do
  
  Vrep    =   co * exp(-ao*r) / r
  Vdlr    =   - c / ( r6 + d6 ) - cx20  / ( r4 + d4 )**2 - cx500 / ( r2 + d2 )**5
  Vsr     =   suma * exp( - a1 * r ) * r6
  V       =   Vrep + Vdlr + Vsr
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module