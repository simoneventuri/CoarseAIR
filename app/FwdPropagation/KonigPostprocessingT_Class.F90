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

Module KonigPostprocessingT_Class

  use Parameters_Module       ,only:  rkp, Zero

  implicit none

  private
  public    ::    KonigPostprocessingT_Type

  Type      ::    KonigPostprocessingT_Type
      
    real(rkp)    ,dimension(:)     ,allocatable             ::    Time
    integer                                                 ::    NTimeSteps
    real(rkp)    ,dimension(:)     ,allocatable             ::    TimeGrid
    
    real(rkp)    ,dimension(:)     ,allocatable             ::    TTranslational
    real(rkp)    ,dimension(:)     ,allocatable             ::    Rho
    real(rkp)    ,dimension(:)     ,allocatable             ::    P
    real(rkp)    ,dimension(:)     ,allocatable             ::    Nd
    real(rkp)    ,dimension(:)     ,allocatable             ::    Energy

    real(rkp)    ,dimension(:,:)   ,allocatable             ::    X
    real(rkp)    ,dimension(:,:,:) ,allocatable             ::    XAtNodes
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    XSum
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    XSqSum
    integer      ,dimension(:,:,:) ,allocatable             ::    XHist
    real(rkp)    ,dimension(:)     ,allocatable             ::    XBinsExtremes
    
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    Tint
    real(rkp)    ,dimension(:,:,:) ,allocatable             ::    TIntAtNodes
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    TIntSum
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    TIntSqSum
    integer      ,dimension(:,:,:) ,allocatable             ::    TIntHist
    real(rkp)    ,dimension(:)     ,allocatable             ::    TIntBinsExtremes
    
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    Pop
    real(rkp)    ,dimension(:,:,:) ,allocatable             ::    PopAtNodes
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    PopSum
    real(rkp)    ,dimension(:,:)   ,allocatable             ::    PopSqSum
    integer      ,dimension(:,:,:) ,allocatable             ::    PopHist
    real(rkp)    ,dimension(:)     ,allocatable             ::    PopBinsExtremes
    
  End Type
    
End Module
