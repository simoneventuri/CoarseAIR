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

Module CHN_UIUC_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Kcm_To_Hartree, B_To_Ang
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use Input_Class           ,only:  Input_Type

  implicit none

  private
  public    :: CHN_UIUC_PES_Type 
  
  Type    ,extends(PES_Type)                :: CHN_UIUC_PES_Type
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
    procedure                               ::  Initialize        =>    Initialize_CHN_UIUC_PES
    procedure                               ::  Output            =>    Output_CHN_UIUC_PES
    procedure                               ::  Compute           =>    Compute_CHN_UIUC_PES_1d
    procedure                               ::  Potential         =>    CHN_UIUC_Potential_From_R
    procedure                               ::  TriatPotential    =>    CHN_UIUC_Potential_From_R_OnlyTriat
  End Type
  
  type(Input_Type)                          ::    Input
  
  logical                     ,parameter    :: i_Debug_Global = .False.

  contains


! **************************************************************************************************************
! **************************************************************************************************************
!                                       DEFERRED PROCEDURES for CHN UIUC PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_CHN_UIUC_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(CHN_UIUC_PES_Type)                  ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms  
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  type(DiatomicPotential_Factory_Type)                       ::    DiatPotFactory   
  integer                                                   ::    iP
  character(*)                    ,parameter                ::    Name_PES = 'CHN_UIUC'
  integer         ,dimension(3,2)                           ::    iA

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_CHN_UIUC_PES" )
  !i_Debug_Loc   =     Logger%On()


  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  ! allocate( This%mMiMn(3) )
  ! This%mMiMn(1:2) = - Atoms(1:2)%Mass / Atoms(3)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )
  
  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiatPotFactory%Construct" )
  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )
  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
  ! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_CHN_UIUC_PES( This, Unit )

  class(CHN_UIUC_PES_Type)                ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function CHN_UIUC_Potential_From_R( This, R, Q ) result( V )
   
  class(CHN_UIUC_PES_Type)                          ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat
  integer                                                    ::    iP 
  
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  !call ComputeTriatomic( A, dAdR, VTriat, dVTriat )
  VTriat = Zero

  !write(*,*) 'R     = ', R
  !write(*,*) 'VDiat = ', VDiat
  !write(*,*) ' '

  V = sum(VDiat) + VTriat 

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function CHN_UIUC_Potential_From_R_OnlyTriat( This, R, Q ) result( V )
   
  class(CHN_UIUC_PES_Type)                          ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat
  integer                                                    ::    iP 
  
  !call ComputeTriatomic( A, dAdR, VTriat, dVTriat )
  VTriat = Zero

  V = VTriat 

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_CHN_UIUC_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(CHN_UIUC_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat
  integer                                                    ::    iP 
  
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  !call ComputeTriatomic( A, dAdR, VTriat, dVTriat )
  VTriat  = Zero
  dVTriat = Zero


  V    = sum(VDiat) + VTriat
  dVdR = dVDiat     + dVTriat

  dVdQ = Zero
  call This%TransToCart_3Atoms( R, Q, dVdR, dVdQ)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
! **************************************************************************************************************
!                                         PRIVATE PROCEDURES for LEPS
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!     

End Module