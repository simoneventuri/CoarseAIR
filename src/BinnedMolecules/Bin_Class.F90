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

Module Bin_Class

  use Parameters_Module       ,only:  rkp, Zero

  implicit none

  private
  public    ::    Bin_Type

  Type   ::    Bin_Type
    real(rkp)                               ::    MineinteV       = Zero             !< Vector (Bin Nb)   -> Bin Min Energy [eV]: from user's input
    real(rkp)                               ::    MinLevEintev    = Zero             !< Vector (Bin Nb)   -> Bin Min Energy [eV]: from minimum energy in the bin levels
    integer                                 ::    vqnFirst        = 0                !< Vector (Bin Nb)   -> First vqn
    integer                                 ::    jqnFirst        = 0                !< Vector (Bin Nb)   -> First jqn
    real(rkp)                               ::    EeVFirst        = Zero
    real(rkp)                               ::    EFirst          = Zero
    integer                                 ::    jqnLast         = 0                !< Vector (Bin Nb)   -> Last jqn
    real(rkp)                               ::    vqnFirst_rkp    = Zero             !< Vector (Bin Nb)   -> First vqn + 0.5
    real(rkp)                               ::    jqnFirst_rkp    = Zero             !< Vector (Bin Nb)   -> First jqn + 0.5
    integer                                 ::    NLevels         = 0                !< Vector (Bin Nb)   -> Nb of Levels contained in the Bin
    real(rkp)                               ::    ToteinteV       = Zero             !< Vector (Bin Nb)   -> Bin Total Energy [eV]
    real(rkp)                               ::    Q               = Zero             !< Vector (Bin Nb)   -> Bin Part Function
    real(rkp)                               ::    QRatio          = Zero             !< Vector (Bin Nb)   -> Bin Part Function Ratio
    real(rkp)                               ::    ToteinteVInit   = Zero             !< Vector (Bin Nb)   -> Bin Total Energy [eV] at Initial Temperature
    real(rkp)                               ::    QInit           = Zero             !< Vector (Bin Nb)   -> Bin Part Function at Initial Temperature
    real(rkp)                               ::    QRatioInit      = Zero             !< Vector (Bin Nb)   -> Bin Part Function Ratio at Initial Temperature
    real(rkp) ,dimension(:,:) ,allocatable  ::    CrossSec  
    real(rkp) ,dimension(:)   ,allocatable  ::    CrossSec_Final 
    real(rkp) ,dimension(:)   ,allocatable  ::    CrossSec_Sigma2
    real(rkp) ,dimension(:,:) ,allocatable  ::    RateConst       
    real(rkp) ,dimension(:,:) ,allocatable  ::    RateConst_Norm  
    real(rkp) ,dimension(:)   ,allocatable  ::    RateConst_Final 
    real(rkp) ,dimension(:)   ,allocatable  ::    RateConst_Sigma2
    real(rkp) ,dimension(:,:) ,allocatable  ::    RateConst_Arr
    real(rkp) ,dimension(:,:) ,allocatable  ::    RateConst_ArrSigma2
    real(rkp) ,dimension(:,:) ,allocatable  ::    RateConst_Arr_Rnd
    real(rkp) ,dimension(:,:) ,allocatable  ::    CArr
    integer                                 ::    Hybrid          = 0
  End Type
  
End Module
