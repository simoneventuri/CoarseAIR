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

Module Null_DiatomicPotential_Class

  use Parameters_Module         ,only:  rkp, Zero, Half, One
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    Null_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    Null_DiatomicPotential_Type
  contains
    procedure               ::    Initialize        =>    Initialize_Null_DiatomicPotential
    procedure               ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_Null
    procedure               ::    DiatomicPotential =>    DiatomicPotential_Null
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiatPot   =   '<NULL>'

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Null_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class   ,only:  Input_Type

  class(Null_DiatomicPotential_Type)    ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Null_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = trim(Name_DiatPot) )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =    iMol
  This%Initialized  =   .True.

  if( Mass2 == Zero ) then
    This%RedMass = One
  else
    This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  end if
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_Null( This, R, V, dV )

  use Parameters_Module   ,only:  Zero

  class(Null_DiatomicPotential_Type)        ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Distances between nuclear centers [bohr]
  real(rkp)                                 ,intent(out)    ::    V                         ! Potential energy [hartree]
  real(rkp)                                 ,intent(out)    ::    dV                        ! First derivative of the potential energy wrt the distance [hartree/bohr]
  
  V   =   Zero * R  ! Just to remove warning message about unused argument
  dV  =   Zero
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_Null( This, R ) result( V )

  use Parameters_Module   ,only:  Zero

  class(Null_DiatomicPotential_Type)        ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Internuclear distance [bohr]
  real(rkp)                                                 ::    V                         ! Diatomic potential energy [hartree]
  
  V   =   Zero * R
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
