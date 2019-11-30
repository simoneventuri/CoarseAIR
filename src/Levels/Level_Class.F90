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

Module Level_Class

  use Parameters_Module       ,only:  rkp, Zero

  implicit none

  private
  public    ::    Level_Type

  Type            ::    Level_Type
    integer       ::    jqn           =   0           ! the rotational q.n. of the i'th quantum state
    integer       ::    vqn           =   0           ! the vibrational q.n. of the i'th quantum state
    real(rkp)     ::    eint          =   Zero        ! internal energy of i'th quantum state in Hartree
    real(rkp)     ::    eint_scaled   =   Zero        ! internal energy of i'th quantum state in Hartree, with reference to the first level
    real(rkp)     ::    einteV        =   Zero        ! internal energy of i'th quantum state in eV
    real(rkp)     ::    einteV_scaled =   Zero        ! internal energy of i'th quantum state in eV, with reference to the first level
    real(rkp)     ::    eintK         =   Zero        ! internal energy of i'th quantum state in K
    real(rkp)     ::    egam          =   Zero        ! Half width of i'th quantum state
    real(rkp)     ::    rmin          =   Zero        ! the position of the potential minimum (included centrifugal potential) for i'th quantum state
    real(rkp)     ::    rmax          =   Zero        ! location of maximum in centrifugal barrier
    real(rkp)     ::    vmin          =   Zero        ! the value of the potential minimun (inc. cent. pot.)
    real(rkp)     ::    vmax          =   Zero        ! the value of the local potential maximum (inc. cent. pot.)
    real(rkp)     ::    tau           =   Zero        ! the vibrational period of the i'th quantum state
    real(rkp)     ::    ri            =   Zero        ! inner turning point
    real(rkp)     ::    ro            =   Zero        ! outter turning point
    real(rkp)     ::    rlim          =   Zero        !
    real(rkp)     ::    g             =   Zero        ! degeneracy of the i'th quantum state
    real(rkp)     ::    to_level      =   Zero        ! in case the Levels Container is a Container for either only Bound States or only Quasi-Bound States, to_level maps (B/QB)Levels -> (Molecule)Levels
    integer       ::    ToBin         =   0
    real(rkp)     ::    gInit         =   Zero        ! degeneracy of the i'th quantum state
    real(rkp)     ::    Q_CG          =   Zero
  End Type

End Module
