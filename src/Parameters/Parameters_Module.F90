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

Module Parameters_Module

  implicit none

  public

  integer ,parameter       ::    rkp = 8
  integer ,parameter       ::    LogUnit = 6

! Math constant
  real(rkp) ,parameter     ::     Onem              = - 1.0_rkp
  real(rkp) ,parameter     ::     Zero              = 0.0_rkp
  real(rkp) ,parameter     ::     Quart             = 0.25_rkp
  real(rkp) ,parameter     ::     Third             = 1.0_rkp / 3.0_rkp
  real(rkp) ,parameter     ::     Half              = 0.5_rkp
  real(rkp) ,parameter     ::     One               = 1.0_rkp
  real(rkp) ,parameter     ::     F4o3              = 4.0_rkp / 3.0_rkp
  real(rkp) ,parameter     ::     Two               = 2.0_rkp
  real(rkp) ,parameter     ::     Three             = 3.0_rkp
  real(rkp) ,parameter     ::     Four              = 4.0_rkp
  real(rkp) ,parameter     ::     Five              = 5.0_rkp
  real(rkp) ,parameter     ::     Six               = 6.0_rkp
  real(rkp) ,parameter     ::     Seven             = 7.0_rkp
  real(rkp) ,parameter     ::     Eight             = 8.0_rkp
  real(rkp) ,parameter     ::     Nine              = 9.0_rkp
  real(rkp) ,parameter     ::     Ten               = 10.0_rkp

  real(rkp) ,parameter     ::     Pi                = acos( - One )
  real(rkp) ,parameter     ::     TwoPi             = Pi * Two

  real(rkp) ,parameter     ::     epss              = 1.0E-4_rkp
  
  real(rkp) ,parameter     ::     UKb               = 1.380658e-23_rkp                    !< Boltzmann's constant [J/K]
  real(rkp) ,parameter     ::     Ue                = 1.602191e-19_rkp                    !< Electron charge in [C]
  real(rkp) ,parameter     ::     Una               = 6.0221367e23_rkp                    !< Avogadro's number [1/mol]
  real(rkp) ,parameter     ::     Rugc              = UKb * Una                           !< Universal gas constant [J/mol/K0]

  real(rkp) ,parameter     ::     K_c0              = 2.99792458e+08_rkp                  !< Speed of light in vacuum [m/s]
  real(rkp) ,parameter     ::     rydberg           = 10973731.568527_rkp
  real(rkp) ,parameter     ::     autime_to_sec     = 0.25_rkp / ( Pi * K_c0 * rydberg )
  
  real(rkp) ,parameter     ::     Hartree_To_Kelvin = 3.1577513E+05_rkp                   !< Conversion factor from [Hartree] to [K]
  real(rkp) ,parameter     ::     Kelvin_To_Hartree = One / Hartree_To_Kelvin             !< Conversion factor from [K] to [Hartree]
  
  real(rkp) ,parameter     ::     Hartree_To_eV     = 27.2113839712790_rkp                !< Conversion factor from [Hartree] to [eV]
  !real(rkp) ,parameter     ::     Hartree_To_eV     = 27.211399_rkp                      !< [Hartree] to [Electron Volt]
  real(rkp) ,parameter     ::     eV_To_Hartree     = One / Hartree_To_eV                 !< Conversion factor from [eV] to [Hartree]

  real(rkp) ,parameter     ::     J_To_eV           = 6.2415093433e+18_rkp                !< Conversion factor from [J] to [eV]
  real(rkp) ,parameter     ::     eV_To_J           = One / J_To_eV                       !< Conversion factor from [eV] to [J]

  real(rkp) ,parameter     ::     Jm_To_eV          = J_To_eV / Una                       !< Conversion factor from [J/mol] to [eV]
  real(rkp) ,parameter     ::     eV_To_Jm          = One / Jm_To_eV                      !< Conversion factor from [eV] to [J/mol]
  
  real(rkp) ,parameter     ::     Kcm_To_Hartree    = 0.159360144e-2_rkp                  !< [kcal/mol] to [hartree] 
  real(rkp) ,parameter     ::     Hartree_To_Kcm    = One / Kcm_To_Hartree
  
  real(rkp) ,parameter     ::     cmm1_To_Hartree   = 0.0000045563352812122_rkp           !< [cm^-1] to [hartree] 
  real(rkp) ,parameter     ::     Hartree_To_cmm1   = One / Kcm_To_Hartree

  real(rkp) ,parameter     ::     hbareV            = 6.582119514e-16_rkp
  real(rkp) ,parameter     ::     hbarH             = hbareV * eV_To_Hartree

  real(rkp) ,parameter     ::     B_To_Ang          = 0.52917721092_rkp                   !< [bohr] to [Angstrom]
  real(rkp) ,parameter     ::     KcmAng_To_HartB   = 0.843297564e-3_rkp                  !< [kcal/(mol*Angstrom)] to [hartree/bohr]

  real(rkp) ,parameter     ::     kgm_To_au         = 1822888.486209_rkp                  !< [kg/mol] to [a.u.] ! 31.9988d-3 : 5.82978912e4 = 1 : x ->>>> x = 5.82978912e4 / 31.9988d-3
  real(rkp) ,parameter     ::     au_To_kgm         = One / kgm_To_au                     !< [a.u.] to [kg/mol]



  integer ,dimension(3)    ,parameter ::     Atom_To_Species_3At  = [1,1,2]
  integer ,dimension(2)    ,parameter ::     NAtomsPerSpecies_3At = [2,1]
  integer ,dimension(2,2)  ,parameter ::     AtomsToSpecies_3At   = Reshape( [1,2, 3,0], [2,2] )
  integer ,dimension(3)    ,parameter ::     OppositePair_3At     = [3,2,1]
  integer ,dimension(3,3)  ,parameter ::     Atoms_To_Pair_3At    = Reshape([0,1,2, 1,0,3, 2,3,0], [3,3] )
  integer ,dimension(2,3)  ,parameter ::     Pair_To_Atoms_3At    = Reshape([1,2, 1,3, 2,3], [2,3] )

  integer ,dimension(4)    ,parameter ::     Atom_To_Species_4At  = [1,1,2,2]
  integer ,dimension(2)    ,parameter ::     NAtomsPerSpecies_4At = [2,2]
  integer ,dimension(2,2)  ,parameter ::     AtomsToSpecies_4At   = Reshape([1,2, 3,4], [2,2] )
  integer ,dimension(6)    ,parameter ::     OppositePair_4At     = [6,5,4,3,2,1]
  integer ,dimension(4,4)  ,parameter ::     Atoms_To_Pair_4At    = Reshape( [0,1,2,3, 1,0,4,5, 2,4,0,6, 3,5,6,0], [4,4] )
  integer ,dimension(2,6)  ,parameter ::     Pair_To_Atoms_4At    = Reshape( [1,2, 1,3, 1,4, 2,3, 2,4, 3,4], [2,6] )

End Module
