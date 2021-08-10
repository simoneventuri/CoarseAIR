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
! The following PES was published by Zuo et al in "Theoretical Investigations of Rate Coefficients for H+O3
! and HO2+O Reactions on a Full-Dimensional Potential Energy Surface" in the
! Journal of Physical Chemistry A 124, 2020. The PES was shared in the online database POTLIB
! (https://comp.chem.umn.edu/potlib/index.html) retrieved on the 2nd of August 2021.
! Re-implemented from the analytic form shared in POTLIB by J. Vargas (joao.francisco.vargas@gmail.com)
!===============================================================================================================
