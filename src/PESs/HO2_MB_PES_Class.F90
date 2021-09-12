! -*-F90-*-
!===============================================================================================================
!
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR)
!
! Copyright (C) 2021 João Vargas (King Abdullah University of Science and Technology).
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
! Based on the HO2 MB PES found in POTLIB
!===============================================================================================================

Module HO2_MB_PES_Class

#include "../qct.inc"

!   use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Kcm_To_Hartree, B_To_Ang
!   use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
!   use Logger_Class          ,only:  Logger
!   use Error_Class           ,only:  Error
!   use Input_Class           ,only:  Input_Type


  implicit none

  private
  public    :: HO2_MB_PES_Type

  Type    ,extends(PES_Type)                :: HO2_MB_PES_Type
    real(rkp), dimension(3)                 :: re
    real(rkp), dimension(3)                 :: DeB
    real(rkp), dimension(3)                 :: DeA
    real(rkp), dimension(3)                 :: alpha
    real(rkp), dimension(3)                 :: beta
    real(rkp)                               :: S
    real(rkp)                               :: k
    real(rkp)                               :: A0
    real(rkp)                               :: A2
    real(rkp)                               :: A4
    character(10)                           :: model_output
    real(rkp)                               :: Rmin
    real(rkp)                               :: Vmin
    logical                                 :: ComputeDiatPotFlg = .False.
  contains
    procedure                               ::  Initialize        =>    Initialize_HO2_MB_PES
    procedure                               ::  Output            =>    Output_HO2_MB_PES
    procedure                               ::  Compute           =>    Compute_HO2_MB_PES_1d
    procedure                               ::  Potential         =>    HO2_MB_Potential_From_R
    procedure                               ::  TriatPotential    =>    HO2_MB_Potential_From_R_OnlyTriat
  End Type

  type(Input_Type)                          ::    Input

  logical                     ,parameter    :: i_Debug_Global = .False.

  contains
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_HO2_MB_PES()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_HO2_MB_PES()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function HO2_MB_Potential_From_R() result(V)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function HO2_MB_Potential_From_R_OnlyTriat() result(V)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_HO2_MB_PES_1d()

End Subroutine
!________________________________________________________________________________________________________________________________!

End Module



