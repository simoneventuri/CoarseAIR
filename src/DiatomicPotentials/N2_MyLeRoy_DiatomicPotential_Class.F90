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
  character(*)  ,parameter  ::    Name_DiaPot = 'LeRoy'

  integer   ,parameter                                      ::    nrep  = 12
  real(rkp) ,parameter                                      ::    co    =  49.0E+00_rkp
  real(rkp) ,parameter                                      ::    ao    =  2.4E+00_rkp
  real(rkp) ,parameter                                      ::    a1    =  5.38484562702061E+00_rkp
  real(rkp) ,parameter                                      ::    c     =  23.9717901220746E+00_rkp
  real(rkp) ,parameter                                      ::    cx20  =  c * 20
  real(rkp) ,parameter                                      ::    cx500 =  c * 500
  real(rkp) ,parameter                                      ::    d     =  3E+00_rkp
  real(rkp) ,parameter                                      ::    d2    =  d**2
  real(rkp) ,parameter                                      ::    d4    =  d**4
  real(rkp) ,parameter                                      ::    d6    =  d**6
  real(rkp) ,parameter                                      ::    re    =  2.1E+00_rkp
  real(rkp) ,parameter  ,dimension(nrep)                    ::    csr   = [ -3.89327515161896D+02, -6.19410598346194D+02, -5.51461129947346D+02, &
                                                                            -3.54896837660797D+02, -1.08347448451266D+02, -6.59348244094835D+01, &
                                                                             1.30457802135760D+01,  7.20557758084161D+01, -1.81062803146583D+01, &
                                                                            -2.84137057950277D+01,  1.40509401686544D+01, -1.84651681798865D+00]

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_N2_LeRoy_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )
  
  use Input_Class   ,only:  Input_Type
  
  class(N2_LeRoy_DiatomicPotential_Type)      ,intent(out)  ::    This
  type(Input_Type)                            ,intent(in)   ::    Input
  character(:)  ,allocatable                  ,intent(in)   ::    SpeciesName
  integer                                     ,intent(in)   ::    iMol
  real(rkp)                                   ,intent(in)   ::    Mass1
  real(rkp)                                   ,intent(in)   ::    Mass2
  logical ,optional                           ,intent(in)   ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N2_LeRoy_DiatomicPotential" )
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
Pure Elemental Subroutine Compute_Vd_dVd_N2_LeRoy( This, R, V, dV )

  use Parameters_Module   ,only:  Zero, One

  class(N2_LeRoy_DiatomicPotential_Type)    ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Distances between nuclear centers [bohr]
  real(rkp)                                 ,intent(out)    ::    V                         ! Potential energy [hartree]
  real(rkp)                                 ,intent(out)    ::    dV                        ! First derivative of the potential energy wrt the distance [hartree/bohr]

  integer                                                   ::    n
  real(rkp)                                                 ::    r2, r3, r4, r5, r6, ri
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
  dV      =   r - re

!   write(Logger%Unit,"(6x,'[vdiat]: drep = ',es15.8)") drep

  do n = nrep-1,1,-1
    dsuma =   suma + dV * dsuma
    suma  =   csr(n) + suma * dV
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

!   write(Logger%Unit,"(6x,'[vdiat]: Vrep = ',es15.8)") Vrep
!   write(Logger%Unit,"(6x,'[vdiat]: Vdlr = ',es15.8)") Vdlr
!   write(Logger%Unit,"(6x,'[vdiat]: Vsr  = ',es15.8)") Vsr
!   write(Logger%Unit,"(6x,'[vdiat]: V    = ',es15.8)") V
!   write(Logger%Unit,"(6x,'[vdiat]: dV   = ',es15.8)") dV

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Elemental Function DiatomicPotential_N2_LeRoy( This, R ) result( V )

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
