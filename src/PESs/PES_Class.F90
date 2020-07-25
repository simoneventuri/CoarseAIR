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

Module PES_Class

#include "../qct.inc"

  use Parameters_Module         ,only:  rkp, Zero, One
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type
  use File_Class                ,only:  File_Type
  
  implicit none

  private
  public  ::    PES_Type
  public  ::    DiatPotContainer_Type
  public  ::    PESEvoFile, PESEvoFlg
  
  Type                                                      ::    DiatPotContainer_Type
    class(DiatomicPotential_Type) ,allocatable              ::    Vd
  End Type

  Type    ,abstract                                         ::    PES_Type
    logical                                                 ::    Initialized         !< Indicator whether the object is initialized
    integer                                                 ::    NPairs              !< Number of pairs
    character(:)                               ,allocatable ::    Name                !< Name of current PES
    character(:)                               ,allocatable ::    Model               !< Name of current PES
    real(rkp)                    ,dimension(:) ,allocatable ::    mMiMn               !< Opposite of the ratio of the mass of the N-1 first atoms over the mass of 
                                                                                      !< the last atoms: -Mi(1:N-1)/M(N). Dim=(NAtoms-1)
    type(DiatPotContainer_Type)  ,dimension(:) ,allocatable ::    Pairs               !< Vector of Pair objects. Each pair can have a different diatomic-potential object
    logical                                                 ::    CartCoordFlg

    integer   ,dimension(3)                                 ::    UseSurface
  contains
    private
    procedure              ,public                          ::    Initialize       =>    Initialize_PES
    procedure              ,public                          ::    Output           =>    Output_PES

    procedure              ,public                          ::    Potential        =>    Potential_0d
    procedure              ,public                          ::    TriatPotential   =>    TriatPotential_0d
    procedure              ,public                          ::    DiatPotential    =>    DiatPotential_0d
    
    procedure              ,public                          ::    Compute          =>    Compute_0d
    
    procedure              ,public                          ::    SampleParamPost  =>    SampleParamPost_PES
    procedure              ,public                          ::    ReadParamPost    =>    ReadParamPost_PES
    procedure              ,public                          ::    WriteParamSample =>    WriteParamSample_PES

    procedure              ,public                          ::    TransToCart_3Atoms
    procedure              ,public                          ::    TransToCart_4Atoms

  End Type
  
  type(File_Type) ,dimension(:) ,allocatable                ::    PESEvoFile
  real(rkp)       ,dimension(:) ,allocatable                ::    tPES
  logical                                                   ::    PESEvoFlg = .false.

  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                 ,only:  Input_Type
  use Atom_Class                  ,only:  Atom_Type
  
  class(PES_Type)                           ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_PES" )
  !i_Debug_Loc   =     Logger%On()
  
  
  This%Name         =   '<Unknown>'
  This%Model        =   '<Unknown>'
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   0
  
  if (i_Debug_Loc) write(Logger%Unit,"(4x,'[Initialize_PES]: Nothing to do here')")
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_PES( This, Unit )
  class(PES_Type)                           ,intent(in)     ::    This
  integer                                   ,intent(in)     ::    Unit
  integer                                                   ::    idum
  logical                                                   ::    ldum
  idum = Unit
  ldum = This%Initialized
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_0d( This, R, Q, V, dVdR, dVdQ )
  class(PES_Type)                           ,intent(in)     ::    This
  real(rkp) ,dimension(:)  CONTIGUOUS       ,intent(in)     ::    R            ! Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)  CONTIGUOUS       ,intent(in)     ::    Q            ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                 ,intent(out)    ::    V            ! Potential energy in [hartree]
  real(rkp) ,dimension(:)  CONTIGUOUS       ,intent(out)    ::    dVdR         ! Derivative of the Potential energy [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)  CONTIGUOUS       ,intent(out)    ::    dVdQ         ! Derivative of the Potential energy [hartree/bohr]. Dim=(NPairs)
  
  integer                                   ,parameter      ::    N = 1       ! Number of trajectories
  real(rkp) ,dimension(size(R),N)                           ::    R_          ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(size(Q),N)                           ::    Q_          ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(N)                                   ::    V_          ! Potential energy in [hartree]. Dim=(NTraj)
  real(rkp) ,dimension(size(R),N)                           ::    dVdR_       ! Derivative of the Potential energy [hartree/bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(size(Q),N)                           ::    dVdQ_       ! Derivative of the Potential energy [hartree/bohr]. Dim=(NPairs,NTraj)
  
  V    = Zero
  dVdR = Zero
  dVdQ = Zero
   
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function Potential_0d( This, R, Q ) result( V )
  class(PES_Type)                              ,intent(in)  ::    This
  real(rkp) ,dimension(:)         CONTIGUOUS   ,intent(in)  ::    R           ! Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)         CONTIGUOUS   ,intent(in)  ::    Q           ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  
  real(rkp)                                                 ::    V           ! Potential energy in [hartree].
  
  integer                                   ,parameter      ::    N = 1       ! Number of trajectories
  real(rkp) ,dimension(size(R),N)                           ::    R_          ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(size(Q),N)                           ::    Q_          ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3,NTraj) 
  real(rkp) ,dimension(N)                                   ::    V_          ! Potential energy in [hartree]. Dim=(NTraj)
  
  V  = Zero
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function DiatPotential_0d( This, R, Q ) result( V )
  class(PES_Type)                              ,intent(in)  ::    This
  real(rkp) ,dimension(:)         CONTIGUOUS   ,intent(in)  ::    R           ! Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)         CONTIGUOUS   ,intent(in)  ::    Q           ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  
  real(rkp)                                                 ::    V           ! Potential energy in [hartree].
  
  integer                                   ,parameter      ::    N = 1       ! Number of trajectories
  real(rkp) ,dimension(size(R),N)                           ::    R_          ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(size(Q),N)                           ::    Q_          ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3,NTraj) 
  real(rkp) ,dimension(N)                                   ::    V_          ! Potential energy in [hartree]. Dim=(NTraj)
  
  V  = Zero
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function TriatPotential_0d( This, R, Q ) result( V )
  class(PES_Type)                              ,intent(in)  ::    This
  real(rkp) ,dimension(:)         CONTIGUOUS   ,intent(in)  ::    R           ! Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)         CONTIGUOUS   ,intent(in)  ::    Q           ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  
  real(rkp)                                                 ::    V           ! Potential energy in [hartree].
  
  integer                                   ,parameter      ::    N = 1       ! Number of trajectories
  real(rkp) ,dimension(size(R),N)                           ::    R_          ! Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(size(Q),N)                           ::    Q_          ! Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3,NTraj) 
  real(rkp) ,dimension(N)                                   ::    V_          ! Potential energy in [hartree]. Dim=(NTraj)
  
  V  = Zero
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SampleParamPost_PES(This)
! This procedure determines the terms for an analytic fit for a molecule A3 using Jacobian coordinates
! Removed:
!  - Arguments: dr, drs, dep, dpl, dct
!  - Local:     dd

  class(PES_Type)                          ,intent(inout)  ::    This

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadParamPost_PES(This, PESiseed)
! This procedure determines the terms for an analytic fit for a molecule A3 using Jacobian coordinates
! Removed:
!  - Arguments: dr, drs, dep, dpl, dct
!  - Local:     dd

  class(PES_Type)                          ,intent(inout)  ::    This
  integer                                   ,intent(in)    ::    PESiseed

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteParamSample_PES(This, PESiseed, LevelOutputDir)
! This procedure determines the terms for an analytic fit for a molecule A3 using Jacobian coordinates
! Removed:
!  - Arguments: dr, drs, dep, dpl, dct
!  - Local:     dd

  class(PES_Type)                          ,intent(in)  ::    This
  integer                                  ,intent(in)  ::    PESiseed
  character(150)                           ,intent(in)  ::    LevelOutputDir

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



!________________________________________________________________________________________________________________________________!
Subroutine TransToCart_3Atoms(This, R, Q, dVdR, mdVdQ)

  use Transformation_Class, only: dR_To_dX

  class(PES_Type)                           ,intent(in)  ::    This
  real(rkp)       ,dimension(3)             ,intent(in)  ::    R
  real(rkp)       ,dimension(9)             ,intent(in)  ::    Q
  real(rkp)       ,dimension(3)             ,intent(in)  ::    dVdR
  real(rkp)       ,dimension(6)             ,intent(out) ::    mdVdQ

  real(rkp)       ,dimension(9)                          ::    dVdQ
  real(rkp)       ,dimension(3)                          ::    Rx, Ry, Rz, dVdRi_

  Rx        = [Q(1)-Q(4), Q(1)-Q(7), Q(4)-Q(7)]
  Ry        = [Q(2)-Q(5), Q(2)-Q(8), Q(5)-Q(8)]
  Rz        = [Q(3)-Q(6), Q(3)-Q(9), Q(6)-Q(9)]
  dVdRi_    = dVdR / R


  ! !!! OLD IMPLEMENTATION
  ! !
  ! mdVdQ(1)  =   -dVdRi_(1) * (Q(1)-Q(4)) - dVdRi_(2) * (Q(1)-Rxyz(1)) * (one-This%mMiMn(1)) + dVdRi_(3) * (Q(4)-Rxyz(1))*This%mMiMn(1)
  ! mdVdQ(2)  =   -dVdRi_(1) * (Q(2)-Q(5)) - dVdRi_(2) * (Q(2)-Rxyz(2)) * (one-This%mMiMn(1)) + dVdRi_(3) * (Q(5)-Rxyz(2))*This%mMiMn(1)
  ! mdVdQ(3)  =   -dVdRi_(1) * (Q(3)-Q(6)) - dVdRi_(2) * (Q(3)-Rxyz(3)) * (one-This%mMiMn(1)) + dVdRi_(3) * (Q(6)-Rxyz(3))*This%mMiMn(1)
  ! mdVdQ(4)  =    dVdRi_(1) * (Q(1)-Q(4)) + dVdRi_(2) * (Q(1)-Rxyz(1)) * This%mMiMn(2)       - dVdRi_(3) * (Q(4)-Rxyz(1))*(one-This%mMiMn(2))
  ! mdVdQ(5)  =    dVdRi_(1) * (Q(2)-Q(5)) + dVdRi_(2) * (Q(2)-Rxyz(2)) * This%mMiMn(2)       - dVdRi_(3) * (Q(5)-Rxyz(2))*(one-This%mMiMn(2))
  ! mdVdQ(6)  =    dVdRi_(1) * (Q(3)-Q(6)) + dVdRi_(2) * (Q(3)-Rxyz(3)) * This%mMiMn(2)       - dVdRi_(3) * (Q(6)-Rxyz(3))*(one-This%mMiMn(2))


  !!! IMPLEMENTATION 2
  !
  mdVdQ(1)  =  dVdRi_(1)*Rx(1) +(One-This%mMiMn(1))*dVdRi_(2)*Rx(2)       -This%mMiMn(1)*dVdRi_(3)*Rx(3)
  mdVdQ(2)  =  dVdRi_(1)*Ry(1) +(One-This%mMiMn(1))*dVdRi_(2)*Ry(2)       -This%mMiMn(1)*dVdRi_(3)*Ry(3)
  mdVdQ(3)  =  dVdRi_(1)*Rz(1) +(One-This%mMiMn(1))*dVdRi_(2)*Rz(2)       -This%mMiMn(1)*dVdRi_(3)*Rz(3)

  mdVdQ(4)  = -dVdRi_(1)*Rx(1)       -This%mMiMn(2)*dVdRi_(2)*Rx(2) +(One-This%mMiMn(2))*dVdRi_(3)*Rx(3)
  mdVdQ(5)  = -dVdRi_(1)*Ry(1)       -This%mMiMn(2)*dVdRi_(2)*Ry(2) +(One-This%mMiMn(2))*dVdRi_(3)*Ry(3)
  mdVdQ(6)  = -dVdRi_(1)*Rz(1)       -This%mMiMn(2)*dVdRi_(2)*Rz(2) +(One-This%mMiMn(2))*dVdRi_(3)*Rz(3)
  

  ! !!! IMPLEMENTATION 3
  ! !
  ! call dR_To_dX(R, Q, dVdR=dVdR, dVdX=dVdQ)

  ! mdVdQ(1) = dVdQ(1) + This%mMiMn(1)*dVdQ(7)
  ! mdVdQ(2) = dVdQ(2) + This%mMiMn(1)*dVdQ(8)
  ! mdVdQ(3) = dVdQ(3) + This%mMiMn(1)*dVdQ(9)

  ! mdVdQ(4) = dVdQ(4) + This%mMiMn(2)*dVdQ(7)
  ! mdVdQ(5) = dVdQ(5) + This%mMiMn(2)*dVdQ(8)
  ! mdVdQ(6) = dVdQ(6) + This%mMiMn(2)*dVdQ(9)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



!________________________________________________________________________________________________________________________________!
Subroutine TransToCart_4Atoms(This, R, Q, dVdR, mdVdQ)

  use Transformation_Class, only: dR_To_dX

  class(PES_Type)                           ,intent(in)  ::    This
  real(rkp)       ,dimension(6)             ,intent(in)  ::    R
  real(rkp)       ,dimension(12)            ,intent(in)  ::    Q
  real(rkp)       ,dimension(6)             ,intent(in)  ::    dVdR
  real(rkp)       ,dimension(9)             ,intent(out) ::    mdVdQ

  real(rkp)       ,dimension(12)                         ::    dVdQ
  real(rkp)       ,dimension(6)                          ::    Rx, Ry, Rz, dVdRi_


  !!! IMPLEMENTATION 2
  !
  Rx        = [Q(1)-Q(4), Q(1)-Q(7), Q(1)-Q(10), Q(4)-Q(7), Q(4)-Q(10), Q(7)-Q(10)]
  Ry        = [Q(2)-Q(5), Q(2)-Q(8), Q(2)-Q(11), Q(5)-Q(8), Q(5)-Q(11), Q(8)-Q(11)]
  Rz        = [Q(3)-Q(6), Q(3)-Q(9), Q(3)-Q(12), Q(6)-Q(9), Q(6)-Q(12), Q(9)-Q(12)]
  dVdRi_    = dVdR / R

  mdVdQ(1)  =  dVdRi_(1)*Rx(1) +dVdRi_(2)*Rx(2) +(One-This%mMiMn(1))*dVdRi_(3)*Rx(3)                        -This%mMiMn(1)*dVdRi_(5)*Rx(5)       -This%mMiMn(1)*dVdRi_(6)*Rx(6)
  mdVdQ(2)  =  dVdRi_(1)*Ry(1) +dVdRi_(2)*Ry(2) +(One-This%mMiMn(1))*dVdRi_(3)*Ry(3)                        -This%mMiMn(1)*dVdRi_(5)*Ry(5)       -This%mMiMn(1)*dVdRi_(6)*Ry(6)
  mdVdQ(3)  =  dVdRi_(1)*Rz(1) +dVdRi_(2)*Rz(2) +(One-This%mMiMn(1))*dVdRi_(3)*Rz(3)                        -This%mMiMn(1)*dVdRi_(5)*Rz(5)       -This%mMiMn(1)*dVdRi_(6)*Rz(6)

  mdVdQ(4)  = -dVdRi_(1)*Rx(1)                        -This%mMiMn(2)*dVdRi_(3)*Rx(3) +dVdRi_(4)*Rx(4) +(One-This%mMiMn(2))*dVdRi_(5)*Rx(5)       -This%mMiMn(2)*dVdRi_(6)*Rx(6)
  mdVdQ(5)  = -dVdRi_(1)*Ry(1)                        -This%mMiMn(2)*dVdRi_(3)*Ry(3) +dVdRi_(4)*Ry(4) +(One-This%mMiMn(2))*dVdRi_(5)*Ry(5)       -This%mMiMn(2)*dVdRi_(6)*Ry(6)
  mdVdQ(6)  = -dVdRi_(1)*Rz(1)                        -This%mMiMn(2)*dVdRi_(3)*Rz(3) +dVdRi_(4)*Rz(4) +(One-This%mMiMn(2))*dVdRi_(5)*Rz(5)       -This%mMiMn(2)*dVdRi_(6)*Rz(6)

  mdVdQ(7)  =                  -dVdRi_(2)*Rx(2)       -This%mMiMn(3)*dVdRi_(3)*Rx(3) -dVdRi_(4)*Rx(4)       -This%mMiMn(3)*dVdRi_(5)*Rx(5) +(One-This%mMiMn(3))*dVdRi_(6)*Rx(6)
  mdVdQ(8)  =                  -dVdRi_(2)*Ry(2)       -This%mMiMn(3)*dVdRi_(3)*Ry(3) -dVdRi_(4)*Ry(4)       -This%mMiMn(3)*dVdRi_(5)*Ry(5) +(One-This%mMiMn(3))*dVdRi_(6)*Ry(6)
  mdVdQ(9)  =                  -dVdRi_(2)*Rz(2)       -This%mMiMn(3)*dVdRi_(3)*Rz(3) -dVdRi_(4)*Rz(4)       -This%mMiMn(3)*dVdRi_(5)*Rz(5) +(One-This%mMiMn(3))*dVdRi_(6)*Rz(6)


  ! !!! IMPLEMENTATION 3
  ! !
  ! call dR_To_dX(R, Q, dVdR=dVdR, dVdX=dVdQ)

  ! mdVdQ(1) = dVdQ(1) + This%mMiMn(1)*dVdQ(10)
  ! mdVdQ(2) = dVdQ(2) + This%mMiMn(1)*dVdQ(11)
  ! mdVdQ(3) = dVdQ(3) + This%mMiMn(1)*dVdQ(12)

  ! mdVdQ(4) = dVdQ(4) + This%mMiMn(2)*dVdQ(10)
  ! mdVdQ(5) = dVdQ(5) + This%mMiMn(2)*dVdQ(11)
  ! mdVdQ(6) = dVdQ(6) + This%mMiMn(2)*dVdQ(12)

  ! mdVdQ(7) = dVdQ(7) + This%mMiMn(3)*dVdQ(10)
  ! mdVdQ(8) = dVdQ(8) + This%mMiMn(3)*dVdQ(11)
  ! mdVdQ(9) = dVdQ(9) + This%mMiMn(3)*dVdQ(12)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
