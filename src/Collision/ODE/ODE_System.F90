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

Module ODE_System

#include "../../qct.inc"

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error, CheckVariable
  use Collision_Class       ,only:  Collision_Type

  implicit none

  private
  public    ::    HamiltonODESystem
  public    ::    Collision

  logical   ,parameter                ::    NanCheck    =   .False.
  type(Collision_Type)    ,pointer    ::    Collision => null()

  contains


!________________________________________________________________________________________________________________________________!
PURITY Subroutine HamiltonODESystem( iPES, y, yp )
! This procedure computes the derivatives of the hamiltonian for Hamilton's equations of motion.

  use Parameters_Module       ,only:  Zero

  integer   ,dimension(:)                   ,intent(in)     ::    iPES 
  real(rkp) ,dimension(:,:)                 ,intent(in)     ::    y                     ! Solution             Dim=(NEqtTot,NTraj)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    yp                    ! Solution derivatives Dim=(NEqtTot,NTraj)

  integer                                                   ::    NTraj
  integer                                                   ::    NEqtTot
  integer                                                   ::    NEqtVar
  integer                                                   ::    i, j, jmin, im, k
  real(rkp) ,dimension(                  size(y,2))         ::    V
  real(rkp) ,dimension(Collision%NEqtVar,size(y,2))         ::    Q                     ! Position. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension(Collision%NEqtVar,size(y,2))         ::    dVdQ                  ! potential derivatives. Dim=(NEqtVar,NTraj)
  real(rkp) ,dimension(Collision%NPairs, size(y,2))         ::    Rpi                   ! Inverse of the distances of atom-atom pairs [1/bohr]. Dim=(NPairs,NTraj)

  NEqtTot   =   size(y,1)
  NEqtVar   =   NEqtTot / 2
  NTraj     =   size(y,2)

  if ( .Not. Collision%NormalKineticEnergy ) then
    do k = 1,NTraj
      do i = 1,NEqtVar
        yp(NEqtVar+i,k)   =   Zero
      end do
      do i = 1,NEqtVar
        jmin    =   mod(i,3)
        if ( jmin == 0 ) jmin = 3
        do j = jmin,NEqtVar,3
          yp(NEqtVar+j,k) = yp(NEqtVar+j,k) + Collision%Transformation%Tpq(j,i) * y(i,k)
        end do
      end do
      Q(:,k)    =     y(NEqtVar+1:,k)                                                   ! Copying data in contiguous variable
    end do
  else
    do k = 1,NTraj
      do i = 1,NEqtVar
        im      =     i / 3
        if ( im*3 /= i ) im = im+1
        yp(NEqtVar+i,:)   =   y(i,:) * Collision%mMiMn(im)
      end do
      Q(:,k)    =     y(NEqtVar+1:,k)                                                   ! Copying data in contiguous variable
    end do
  end if

  call Collision%Compute_PES( iPES, Q, dVdQ, V, Rpi )                                   ! Computing the inter-molecular potential using only contiguous variables
  
  yp(1:NEqtVar,:)   =   dVdQ                                                            ! Copying back the data to the final variable

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
