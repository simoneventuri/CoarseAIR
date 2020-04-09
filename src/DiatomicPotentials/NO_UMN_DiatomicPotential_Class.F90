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

Module NO_UMN_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Zero, Half, One, Two, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type

  implicit none

  private
  public  ::    NO_UMN_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    NO_UMN_DiatomicPotential_Type
  contains
    procedure         ::    Initialize        =>    Initialize_NO_UMN_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_NO
    procedure         ::    DiatomicPotential =>    DiatomicPotential_NO
  End Type

  character(*)  ,parameter  ::    Name_DiatPot    = 'UMN'
  logical       ,parameter  ::    i_Debug_Global = .False.
  
  real(rkp) ,dimension(1:10)                ,parameter      ::    cs   = [ 0.322338e0_rkp, 5.878590e0_rkp, -12.790761e0_rkp, 13.320811e0_rkp, -7.516309e0_rkp, 1.875839e0_rkp, -0.052723e0_rkp, -0.037783e0_rkp, 0.48294e0_rkp, 1.98697e0_rkp] 
  real(rkp)                                 ,parameter      ::    red  = 1.1508_rkp ! Parameter for NO dissociation with SEC
  real(rkp)                                 ,parameter      ::    de   = 152.6_rkp  ! Parameter for NO dissociation with SEC
  !real(rkp)                                 ,parameter      ::    VRef = 0.191619504727d0 ! Same reference as UMN database
  real(rkp)                                 ,parameter      ::    VRef = Zero ! Same reference as David Schwenke

  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_NO_UMN_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(NO_UMN_DiatomicPotential_Type)  ,intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_NO_UMN_DiatomicPotential" )
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
Elemental Function DiatomicPotential_NO( This, R ) result( V )

  class(NO_UMN_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                                          :: V
    
  real(rkp)                                    :: RAng
  real(rkp)                                    :: VDiat
  
  !Reference energy of infinitely separated N2 + O in hartree (taken from DSEC corrected calculations) = 0.191619504727d0 Hartree
  
  RAng = R * B_To_Ang                                                                                                                   

  call Ev2gm2(RAng, VDiat)
  
  V = VDiat * Kcm_To_Hartree + VRef
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!  

!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_NO( This, R, V, dV )

  class(NO_UMN_DiatomicPotential_Type)  ,intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                             ,intent(out) :: V
  real(rkp)                             ,intent(out) :: dV
  
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

  real(rkp)                                 ,intent(in)     ::    RAng                                            ! Distances of atom-atom pairs [Angstrom]
  real(rkp)                                 ,intent(out)    ::    V   
  
  integer                                                   ::    i
  real(rkp)                                                 ::    u
  real(rkp), dimension(8)                                   ::    uvec

  ! MEG variable
  u    =   exp( -(RAng-red)/cs(9) - (RAng-red)**2.0d0 / (cs(10)) ) 

  ! Pairwise potential
  uvec(1) = u
  do i = 2,8
    uvec(i) = uvec(i - 1)*u
  enddo
  V = - de*dot_product(cs(1:8), uvec)
  !V  =  - de * (cs(1)*u + cs(2)*u**2.0d0 + cs(3)*u**3.0d0 + cs(4)*u**4.0d0 + cs(5)*u**5.0d0 + cs(6)*u**6.0d0 + cs(7)*u**7.0d0 + cs(8)*u**8.0d0)
  
End Subroutine Ev2gm2
!--------------------------------------------------------------------------------------------------------------------------------! 


!________________________________________________________________________________________________________________________________!
Pure Subroutine Ev2gm2_Grad( RAng, V, dV) 

  use Parameters_Module     ,only:  Zero, Two

  real(rkp)                                 ,intent(in)     ::    RAng              ! Distances of atom-atom pairs [Angstrom]
  real(rkp)                                 ,intent(out)    ::    V   
  real(rkp)                                 ,intent(out)    ::    dV 
  
  real(rkp)                                                 ::    u, dfdr

  ! MEG variable
  u    =   exp( -(RAng-red)/cs(9) - (RAng-red)**2.0d0 / (cs(10)) ) 
  dfdr =  (-2.0d0 * (RAng-red)/cs(10) - 1.0d0/cs(9))  

  ! Pairwise potential
  V  =  - de * (cs(1)*u + cs(2)*u**2.0d0 + cs(3)*u**3.0d0 + cs(4)*u**4.0d0 + cs(5)*u**5.0d0 + cs(6)*u**6.0d0 + cs(7)*u**7.0d0 + cs(8)*u**8.0d0)

  dV =  - de*(cs(1)      *dfdr*u            + &
              cs(2)*2.0d0*dfdr*u**2.0d0     + &
              cs(3)*3.0d0*dfdr*u**3.0d0     + &
              cs(4)*4.0d0*dfdr*u**4.0d0     + &
              cs(5)*5.0d0*dfdr*u**5.0d0     + &
              cs(6)*6.0d0*dfdr*u**6.0d0     + &
              cs(7)*7.0d0*dfdr*u**7.0d0     + &
              cs(8)*8.0d0*dfdr*u**8.0d0)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
