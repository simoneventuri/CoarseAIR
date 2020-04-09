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

Module MolecularStateParam_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger

  implicit none

  private
  public    ::    MolecularStateParam_Type


  Type      ::    MolecularStateParam_Type
    integer               ::    nquad       =   0
    real(rkp)             ::    rcent       =   Zero      ! Upper limit on position of centrifugal barriar maximum
    real(rkp)             ::    rsmal       =   Zero
    real(rkp)             ::    Vinf        =   Zero      ! @TODO: This should be different for different pairs
    real(rkp)             ::    rfact       =   Zero
    real(rkp)             ::    tthresau    =   Zero
  contains
    procedure ,public   ::    Initialize    =>    InitializeMolecularStateParam
  End Type

  logical   ,parameter    ::    i_Debug_Global = .True.!.False.

  integer                                         ::    iSpeTar = 1
  integer                                         ::    iSpePro = 2

  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeMolecularStateParam( This, Input, Collision, i_Debug )

  use Input_Class               ,only:  Input_Type
  use Collision_Class           ,only:  Collision_Type

  class(MolecularStateParam_Type)           ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  real(rkp)                                    ,parameter   ::    Rinf = 1.0E+03

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeMolecularStateParam")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()


! ==============================================================================================================
!     SETTING VARIABLES FROM THE INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Setting the parameters from the Input object")
  This%nquad      =   Input%nquad
  This%rcent      =   Input%rcent
  This%rsmal      =   Input%rsmal
  This%rfact      =   Input%rfact
  This%tthresau   =   Input%tthresau
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING THE DIATOMIC POTENITAL OF THE target DISTANCE FOR AN INFINITY BOND LENGTH
! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Computing the diatomic potenital of the target distance for an infinity bond length")
  if (i_Debug_Loc) call Logger%Write( "-> Calling Collision%Species(iSpeTar)%DiatPot%Compute_V")
  This%Vinf       =   Collision%Species(iSpeTar)%DiatPot%DiatomicPotential( Rinf )
  if (i_Debug_Loc) call Logger%Write( "-> Rinf = ", Rinf, "This%Vinf = ", This%Vinf, Fr="es15.8" )
! ==============================================================================================================


  if (i_Debug_Loc) then
    call Logger%Write( "-> Number of quadrature points:      This%nquad    = ", This%nquad                  )
    call Logger%Write( "-> Upper bound for Vmax position:    This%rcent    = ", This%rcent    , Fr="es15.8" )
    call Logger%Write( "-> Lower bound for Vmin position:    This%rsmal    = ", This%rsmal    , Fr="es15.8" )
    call Logger%Write( "-> Potenital of target at infinity:  This%Vinf     = ", This%Vinf     , Fr="es15.8" )
    call Logger%Write( "-> Factor for q.b. states:           This%rfact    = ", This%rfact    , Fr="es15.8" )
    call Logger%Write( "-> Time threshold for q.b. states:   This%tthresau = ", This%tthresau , Fr="es15.8" )
    call Logger%Exiting
  end if

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
