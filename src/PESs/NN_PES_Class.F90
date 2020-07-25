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

Module NN_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, eV_To_Hartree, Zero, Two, Three, Atoms_To_Pair_3At
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    NN_PES_Type

  Type    ,extends(PES_Type)               :: NN_PES_Type
  
    integer   ,dimension(2)                :: NHL
    character(3)                           :: NHL_Char1
    character(3)                           :: NHL_Char2
    
    real(rkp) ,dimension(3)                :: Scaling_MEAN
    real(rkp) ,dimension(3)                :: Scaling_SD
    real(rkp)                              :: y_Scaling
    
    real(rkp) ,dimension(3)                :: Lambda
    real(rkp) ,dimension(3)                :: re
    real(rkp) ,dimension(:,:) ,allocatable :: W1
    real(rkp) ,dimension(:,:) ,allocatable :: W2
    real(rkp) ,dimension(:)   ,allocatable :: W3
    real(rkp) ,dimension(:)   ,allocatable :: b1
    real(rkp) ,dimension(:)   ,allocatable :: b2
    real(rkp)                              :: b3
    character(3) ,dimension(3)             :: AtomName
    integer      ,dimension(3)             :: AtomIdx   
    integer      ,dimension(3)             :: PairToIdx
    integer      ,dimension(3)             :: IdxToPair
    integer                                :: SystemType
    integer                                :: NG   
  contains
    procedure          ::  Initialize        =>    Initialize_NN_PES
    procedure          ::  Output            =>    Output_NN_PES
    procedure          ::  Compute           =>    Compute_NN_PES_1d
    procedure          ::  Potential         =>    NN_Potential_From_R
    procedure          ::  TriatPotential    =>    NN_Potential_From_R_OnlyTriat
  End Type

  logical                  ,parameter    ::    i_Debug_Global = .False.
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_NN_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class    ,only:  DiatomicPotential_Factory_Type

  class(NN_PES_Type)                        ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP, i
  character(:)       ,allocatable                           ::    NN_Weights_Folder
  character(:)       ,allocatable                           ::    NN_Weights_File
  character(*)                    ,parameter                ::    Name_PES = 'NN'
  integer                                                   ::    Status
  integer                                                   ::    Unit
  type(DiatomicPotential_Factory_Type)                      ::    DiatPotFactory
  integer         ,dimension(3,2)                           ::    iA
  integer                                                   ::    jA, kA
    
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_NN_PES" )
  !i_Debug_Loc   =     Logger%On()
    
  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair
  
  ! allocate( This%mMiMn(3) )
  ! This%mMiMn(1:2) = - Atoms(1:2)%Mass / Atoms(3)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )
  
  This%AtomName   = '---'
  This%AtomIdx    = 0           
  This%PairToIdx    = 0  
  This%SystemType = 0  
  This%NG         = 0   

  iA(1,:)         =   [1,2]
  iA(2,:)         =   [1,3]
  iA(3,:)         =   [2,3]


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


  This%NHL = Input%NHiddenLayerNeurons

  write(This%NHL_Char1, '(I3)') This%NHL(1)
  write(This%NHL_Char2, '(I3)') This%NHL(2)

  NN_Weights_Folder = trim(adjustl(Input%DtbPath))  // '/Systems/' // trim(adjustl(Input%System)) // '/PESs/NN/' // trim(adjustl(Input%PES_ParamsFldr(iPES))) // '/' // trim(adjustl(This%NHL_Char1)) // '_' // trim(adjustl(This%NHL_Char2)) // '/'
  if (i_Debug_Loc) call Logger%Write( "Reading NN PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening files in the folder: ", NN_Weights_Folder)
  

  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/Order.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*)
    read(Unit,*) This%AtomName(1)
    if (i_Debug_Loc) call Logger%Write( "Atom 1: ", This%AtomName(1) )
    read(Unit,*) This%AtomName(2)
    if (i_Debug_Loc) call Logger%Write( "Atom 2: ", This%AtomName(2) )
    read(Unit,*) This%AtomName(3)
    if (i_Debug_Loc) call Logger%Write( "Atom 3: ", This%AtomName(3) )
  close(Unit)

  if ( ( adjustl(trim(This%AtomName(1))) == adjustl(trim(This%AtomName(2))) ) .and. ( adjustl(trim(This%AtomName(1))) == adjustl(trim(This%AtomName(3))) ) ) then
    This%SystemType = 3
    This%NG         = 6
  elseif ( ( adjustl(trim(This%AtomName(1))) == adjustl(trim(This%AtomName(2))) ) .or. ( adjustl(trim(This%AtomName(1))) == adjustl(trim(This%AtomName(3))) ) .or. ( adjustl(trim(This%AtomName(3))) == adjustl(trim(This%AtomName(2))) ) ) then
    This%SystemType = 2
    This%NG         = 8
  else 
    This%SystemType = 1
    This%NG         = 8
  end if
  if (i_Debug_Loc) call Logger%Write( "This%SystemType = ", This%SystemType )

  do jA=1,3
    do kA=1,3
      if ( ( adjustl(trim(Atoms(jA)%Name)) == adjustl(trim(This%AtomName(kA))) ) .and. (This%AtomIdx(kA) == 0) ) then
        This%AtomIdx(kA) = jA
        exit
      end if
    end do
  end do

  do iP=1,3
    This%PairToIdx(iP)                 = Atoms_To_Pair_3At( This%AtomIdx( iA(iP,1) ), This%AtomIdx( iA(iP,2) ))
    This%IdxToPair(This%PairToIdx(iP)) = iP
  end do
  if (i_Debug_Loc) call Logger%Write( "Real Pair -> NN   Pair, PairToIdx = ", This%PairToIdx )
  if (i_Debug_Loc) call Logger%Write( "NN Pair   -> Real Pair, IdxToPair = ", This%IdxToPair )

  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/ScalingValues.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*)
    read(Unit,*) This%Scaling_MEAN
    if (i_Debug_Loc) call Logger%Write( "Scaling_MEAN: ", This%Scaling_MEAN )
    read(Unit,*) This%Scaling_SD
    if (i_Debug_Loc) call Logger%Write( "Scaling_SD: ",   This%Scaling_SD )
    read(Unit,*) This%y_Scaling
    if (i_Debug_Loc) call Logger%Write( "y_Scaling: ",    This%y_Scaling )
  close(Unit)



  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Lambda.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do iP=1,3
      read(Unit,*) This%Lambda(iP)
    end do
    if (i_Debug_Loc) call Logger%Write( "Lambda: ", This%Lambda )
  close(Unit)

  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/re.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do iP=1,3
      read(Unit,*) This%re(iP)
    end do
    if (i_Debug_Loc) call Logger%Write( "re: ", This%re )
  close(Unit)



  allocate(This%W1(This%NG,This%NHL(1)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%W1" )
  allocate(This%b1(This%NHL(1)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%b1" )
  allocate(This%W2(This%NHL(1),This%NHL(2)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%W2" )
  allocate(This%b2(This%NHL(2)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%b2" )
  allocate(This%W3(This%NHL(2)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%W3" )

  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W1.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NG
      read(Unit,*) This%W1(i,1:This%NHL(1))
      if (i_Debug_Loc) call Logger%Write( "W1, Line ", i, " = ", This%W1(i,1:This%NHL(1)))
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b1.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(1)
      read(Unit,*) This%b1(i)
    end do
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "b1 = ", This%b1)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W2.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(1)
      read(Unit,*) This%W2(i,1:This%NHL(2))
      if (i_Debug_Loc) call Logger%Write( "W2, Line ", i, " = ", This%W2(i,1:This%NHL(2)))
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b2.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(2)
      read(Unit,*) This%b2(i)
    end do
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "b2 = ", This%b2)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W3.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(2)
      read(Unit,*) This%W3(i)
    end do
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "W3 = ", This%W3)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b3.csv'
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*) This%b3
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "b3 = ", This%b3)

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine



Subroutine Output_NN_PES( This, Unit )

  class(NN_PES_Type)                  ,intent(in)  ::    This
  integer                             ,intent(in)  ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function NN_Potential_From_R( This, R, Q ) result( V )
  
  class(NN_PES_Type)                            ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(size(R))                 ::    p          
  real(rkp) ,dimension(size(R))                 ::    pSqr   
  real(rkp) ,dimension(size(R))                 ::    PSqrSum
  real(rkp) ,dimension(size(R))                 ::    PCube
  real(rkp) ,dimension(size(R))                 ::    PCubeSum
  real(rkp) ,dimension(2*size(R))               ::    GComp
  real(rkp) ,dimension(size(R))                 ::    GVecSum
  real(rkp) ,dimension(size(R))                 ::    GVecProd
  real(rkp) ,dimension(size(R))                 ::    GVecProdSum
  real(rkp)                                     ::    GProd
  real(rkp) ,dimension(This%NG)                 ::    G          
  real(rkp) ,dimension(This%NHL(1))             ::    z1
  real(rkp) ,dimension(This%NHL(1))             ::    y1
  real(rkp) ,dimension(This%NHL(2))             ::    z2
  real(rkp) ,dimension(This%NHL(2))             ::    y2
  real(rkp) ,dimension(size(R),2*size(R))       ::    dG
  integer                                       ::    i,j,iP
  integer                                       ::    NData
  real(rkp)                                     ::    t1, t2
  real(rkp) ,dimension(This%NHL(1),3)           ::    dz1dR
  real(rkp) ,dimension(This%NHL(1),This%NHL(1)) ::    dy1dz1
  real(rkp) ,dimension(This%NHL(1),3)           ::    dy1dR
  real(rkp) ,dimension(This%NHL(2),3)           ::    dz2dR
  real(rkp) ,dimension(This%NHL(2),This%NHL(2)) ::    dVdz2
  real(rkp) ,dimension(size(R))                 ::    GVec1
  real(rkp) ,dimension(2*size(R))               ::    GVec2
  real(rkp)                                     ::    VExp
  real(rkp) ,dimension(This%NHL(2),3)           ::    dVTemp
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
  real(rkp)                                     ::    G2p3, G2t3
      
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  p(1) = dexp( - This%Lambda(1) * (R(This%PairToIdx(1)) - This%re(1)) )
  p(2) = dexp( - This%Lambda(2) * (R(This%PairToIdx(2)) - This%re(2)) )
  p(3) = dexp( - This%Lambda(3) * (R(This%PairToIdx(3)) - This%re(3)) )

  if (This%SystemType == 3) then
    GVecProd    = [p(1)*p(2), p(2)*p(3), p(3)*p(1)]
    GVecSum     = [p(2)+p(3), p(3)+p(1), p(1)+p(2)]
    GComp       = [GVecProd(1)*p(1), GVecProd(1)*p(2), GVecProd(2)*p(2), GVecProd(2)*p(3), GVecProd(3)*p(1), GVecProd(3)*p(3)]
    GProd       = product(p)

    G(1) = sum(GVecProd)
    G(2) = GProd
    G(3) = sum(GComp)  
    G(4) = GComp(1)*p(1) + GComp(2)*p(2) + GComp(3)*p(2) + GComp(4)*p(3) + GComp(5)*p(1) + GComp(6)*p(3)
    G(5) = GComp(1)*p(3) + GComp(3)*p(1) + GComp(6)*p(2) 
    G(6) = GComp(1)*p(2) + GComp(3)*p(3) + GComp(5)*p(3)
  elseif (This%SystemType == 2) then
    G2p3  = p(2) + p(3)
    G2t3  = p(2) * p(3)

    G(1) = p(1)                * G2p3
    G(2) = G2t3
    G(3) = p(1)**2             * G2p3
    G(4) = G2t3                * G2p3
    G(5) = p(1)**3             * G2p3
    G(6) = p(1)**2             * G2p3**2
    G(7) = (p(2)**2 + p(3)**2) * G2t3
    G(8) = G2t3**2
  end if


!  z1   = This%b1
!  call dgemm("N","N",1,This%NHL(1),6,1.0,G,1,This%W1,6,1.0,z1,1)
!  y1   = tanh(z1)
!  z2   = This%b2
!  call dgemm("N","N",1,This%NHL(2),This%NHL(1),1.0,y1,1,This%W2,This%NHL(1),1.0,z2,1)
!  y2   = tanh(z2)
!  V    = sum(y2 * This%W3) + This%b3
  z1   = matmul(G, This%W1)  + This%b1
  y1   = tanh(z1)
  z2   = matmul(y1, This%W2) + This%b2
  y2   = tanh(z2)
  V    = sum(y2 * This%W3) + This%b3
  VExp = exp(V)
  
    
  V  = (VExp - This%y_Scaling ) * eV_To_Hartree + sum(VDiat) 
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function NN_Potential_From_R_OnlyTriat( This, R, Q ) result( V )
  
  class(NN_PES_Type)                            ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(size(R))                 ::    p          
  real(rkp) ,dimension(size(R))                 ::    pSqr   
  real(rkp) ,dimension(size(R))                 ::    PSqrSum
  real(rkp) ,dimension(size(R))                 ::    PCube
  real(rkp) ,dimension(size(R))                 ::    PCubeSum
  real(rkp) ,dimension(2*size(R))               ::    GComp
  real(rkp) ,dimension(size(R))                 ::    GVecSum
  real(rkp) ,dimension(size(R))                 ::    GVecProd
  real(rkp) ,dimension(size(R))                 ::    GVecProdSum
  real(rkp)                                     ::    GProd
  real(rkp) ,dimension(This%NG)                 ::    G          
  real(rkp) ,dimension(This%NHL(1))             ::    z1
  real(rkp) ,dimension(This%NHL(1))             ::    y1
  real(rkp) ,dimension(This%NHL(2))             ::    z2
  real(rkp) ,dimension(This%NHL(2))             ::    y2
  real(rkp) ,dimension(size(R),2*size(R))       ::    dG
  integer                                       ::    i,j,iP
  integer                                       ::    NData
  real(rkp)                                     ::    t1, t2
  real(rkp) ,dimension(This%NHL(1),3)           ::    dz1dR
  real(rkp) ,dimension(This%NHL(1),This%NHL(1)) ::    dy1dz1
  real(rkp) ,dimension(This%NHL(1),3)           ::    dy1dR
  real(rkp) ,dimension(This%NHL(2),3)           ::    dz2dR
  real(rkp) ,dimension(This%NHL(2),This%NHL(2)) ::    dVdz2
  real(rkp) ,dimension(size(R))                 ::    GVec1
  real(rkp) ,dimension(2*size(R))               ::    GVec2
  real(rkp)                                     ::    VExp
  real(rkp) ,dimension(This%NHL(2),3)           ::    dVTemp
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
  real(rkp)                                     ::    G2p3, G2t3
  real(rkp)                                     ::    p0, p1, p2
    
  p(1) = dexp( - This%Lambda(1) * (R(This%PairToIdx(1)) - This%re(1)) )
  p(2) = dexp( - This%Lambda(2) * (R(This%PairToIdx(2)) - This%re(2)) )
  p(3) = dexp( - This%Lambda(3) * (R(This%PairToIdx(3)) - This%re(3)) )

  if (This%SystemType == 3) then
    GVecProd    = [p(1)*p(2), p(2)*p(3), p(3)*p(1)]
    GVecSum     = [p(2)+p(3), p(3)+p(1), p(1)+p(2)]
    GComp       = [GVecProd(1)*p(1), GVecProd(1)*p(2), GVecProd(2)*p(2), GVecProd(2)*p(3), GVecProd(3)*p(1), GVecProd(3)*p(3)]
    GProd       = product(p)

    G(1) = sum(GVecProd)
    G(2) = GProd
    G(3) = sum(GComp)  
    G(4) = GComp(1)*p(1) + GComp(2)*p(2) + GComp(3)*p(2) + GComp(4)*p(3) + GComp(5)*p(1) + GComp(6)*p(3)
    G(5) = GComp(1)*p(3) + GComp(3)*p(1) + GComp(6)*p(2) 
    G(6) = GComp(1)*p(2) + GComp(3)*p(3) + GComp(5)*p(3)
  elseif (This%SystemType == 2) then
    G2p3  = p(2) + p(3)
    G2t3  = p(2) * p(3)

    G(1) = p(1)                * G2p3
    G(2) = G2t3
    G(3) = p(1)**2             * G2p3
    G(4) = G2t3                * G2p3
    G(5) = p(1)**3             * G2p3
    G(6) = p(1)**2             * G2p3**2
    G(7) = (p(2)**2 + p(3)**2) * G2t3
    G(8) = G2t3**2
  end if

!  z1   = This%b1
!  call dgemm("N","N",1,This%NHL(1),6,1.0,G,1,This%W1,6,1.0,z1,1)
!  y1   = tanh(z1)
!  z2   = This%b2
!  call dgemm("N","N",1,This%NHL(2),This%NHL(1),1.0,y1,1,This%W2,This%NHL(1),1.0,z2,1)
!  y2   = tanh(z2)
!  V    = sum(y2 * This%W3) + This%b3
  z1   = matmul(G, This%W1)  + This%b1
  y1   = tanh(z1)
  z2   = matmul(y1, This%W2) + This%b2
  y2   = tanh(z2)
  V    = sum(y2 * This%W3) + This%b3
  VExp = exp(V)
  
  V    = (VExp - This%y_Scaling) * eV_To_Hartree
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_NN_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(NN_PES_Type)                            ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp) ,dimension(size(R))                 ::    p          
  real(rkp) ,dimension(size(R))                 ::    pSqr   
  real(rkp) ,dimension(size(R))                 ::    PSqrSum
  real(rkp) ,dimension(size(R))                 ::    PCube
  real(rkp) ,dimension(size(R))                 ::    PCubeSum
  real(rkp) ,dimension(2*size(R))               ::    GComp
  real(rkp) ,dimension(size(R))                 ::    GVecSum
  real(rkp) ,dimension(size(R))                 ::    GVecProd
  real(rkp) ,dimension(size(R))                 ::    GVecProdSum
  real(rkp)                                     ::    GProd
  real(rkp) ,dimension(This%NG)                 ::    G          
  real(rkp) ,dimension(This%NHL(1))             ::    z1
  real(rkp) ,dimension(This%NHL(1))             ::    y1
  real(rkp) ,dimension(This%NHL(2))             ::    z2
  real(rkp) ,dimension(This%NHL(2))             ::    y2
  real(rkp) ,dimension(size(R),This%NG)         ::    dG
  integer                                       ::    i,j,iP
  integer                                       ::    NData
  real(rkp)                                     ::    t1, t2
  real(rkp) ,dimension(This%NHL(1),3)           ::    dz1dR
  real(rkp) ,dimension(This%NHL(1),This%NHL(1)) ::    dy1dz1
  real(rkp) ,dimension(This%NHL(1),3)           ::    dy1dR
  real(rkp) ,dimension(This%NHL(2),3)           ::    dz2dR
  real(rkp) ,dimension(This%NHL(2),This%NHL(2)) ::    dVdz2
  real(rkp) ,dimension(size(R))                 ::    GVec1
  real(rkp) ,dimension(2*size(R))               ::    GVec2
  real(rkp)                                     ::    VExp
  real(rkp) ,dimension(This%NHL(2),3)           ::    dVTemp
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat, dVTriat
  real(rkp)                                     ::    p1, p1s, p1c, p2, p2s, p2c, p3, p3s, p3c, G2p3, G2t3, G2t3s, p2p3s, p2sp3
    
  !call cpu_time ( t1 )

  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  p(1) = dexp( - This%Lambda(1) * (R(This%PairToIdx(1)) - This%re(1)) )
  p(2) = dexp( - This%Lambda(2) * (R(This%PairToIdx(2)) - This%re(2)) )
  p(3) = dexp( - This%Lambda(3) * (R(This%PairToIdx(3)) - This%re(3)) )

  if (This%SystemType == 3) then
    GVecProd    = [p(1)*p(2), p(2)*p(3), p(3)*p(1)]
    GVecSum     = [p(2)+p(3), p(3)+p(1), p(1)+p(2)]
    GComp       = [GVecProd(1)*p(1), GVecProd(1)*p(2), GVecProd(2)*p(2), GVecProd(2)*p(3), GVecProd(3)*p(1), GVecProd(3)*p(3)]
    GProd       = product(p)

    G(1) = sum(GVecProd)
    G(2) = GProd
    G(3) = sum(GComp)  
    G(4) = GComp(1)*p(1) + GComp(2)*p(2) + GComp(3)*p(2) + GComp(4)*p(3) + GComp(5)*p(1) + GComp(6)*p(3)
    G(5) = GComp(1)*p(3) + GComp(3)*p(1) + GComp(6)*p(2) 
    G(6) = GComp(1)*p(2) + GComp(3)*p(3) + GComp(5)*p(3)
    
    PSqr        = This%Lambda * p**2
    PSqrSum     = [PSqr(2) + PSqr(3), PSqr(3) + PSqr(1), PSqr(1) + PSqr(2)]
    PCube       = This%Lambda * p**3
    PCubeSum    = [PCube(2) + PCube(3), PCube(3) + PCube(1), PCube(1) + PCube(2)]
    GVecProdSum = - This%Lambda * [GVecProd(1) + GVecProd(3), GVecProd(2) + GVecProd(1), GVecProd(3) + GVecProd(2)]
    
    dG(:,1) =   GVecProdSum
    dG(:,2) = - This%Lambda * GProd
    dG(:,3) = - (2.0*GVecSum*PSqr  + p*PSqrSum)
    dG(:,4) = - (3.0*GVecSum*PCube + p*PCubeSum)
    dG(:,5) = - This%Lambda * GProd * (2.0*p + GVecSum)
    dG(:,6) = - 2.0*PSqr*PSqrSum / This%Lambda
  elseif (This%SystemType == 2) then
    p1    = p(1)
    p1s   = p1**2
    p1c   = p1s*p1
    p2    = p(2)
    p2s   = p2**2
    p2c   = p2s*p2
    p3    = p(3)
    p3s   = p3**2
    p3c   = p3s*p3
    G2p3  = p2 + p3
    G2t3  = p2 * p3
    G2t3s = G2t3**2
    p2p3s = p2  * p3s
    p2sp3 = p2s * p3

    G(1) = p1                  * G2p3
    G(2) = G2t3
    G(3) = p1s                 * G2p3
    G(4) = G2t3                * G2p3
    G(5) = p1c                 * G2p3
    G(6) = p1s                 * G2p3**2
    G(7) = (p2s + p3s) * G2t3
    G(8) = G2t3s

    dG(:,1) = [G2p3,                        p1,              p1]
    dG(:,2) = [Zero,                        p3,              p2]
    dG(:,3) = [Two*G(1),                   p1s,             p1s]
    dG(:,4) = [Zero,              Two*G2t3+p3s,    Two*G2t3+p2s]
    dG(:,5) = [Three*G(3),                 p1c,             p1c]
    dG(:,6) = [Two*G(1)**2/p1,        Two*G(3),        Two*G(3)]
    dG(:,7) = [Zero,           Three*p2sp3+p3c, Three*p2p3s+p2c]
    dG(:,8) = [Zero,                 Two*p2p3s,       Two*p2sp3]

    dG(1,:) = - This%Lambda(1) * dG(1,:) * p1
    dG(2,:) = - This%Lambda(2) * dG(2,:) * p2
    dG(3,:) = - This%Lambda(3) * dG(3,:) * p3
  end if
  

!  z1   = This%b1
!  call dgemm("N","N",1,This%NHL(1),6,1.0,G,1,This%W1,6,1.0,z1,1)
!  y1   = tanh(z1)
!  z2   = This%b2
!  call dgemm("N","N",1,This%NHL(2),This%NHL(1),1.0,y1,1,This%W2,This%NHL(1),1.0,z2,1)
!  y2   = tanh(z2)
!  V    = sum(y2 * This%W3) + This%b3
  z1   = matmul(G, This%W1)  + This%b1
  y1   = tanh(z1)
  z2   = matmul(y1, This%W2) + This%b2
  y2   = tanh(z2)
  V    = sum(y2 * This%W3) + This%b3
  VExp = exp(V)
  
  V    = (VExp - This%y_Scaling) * eV_To_Hartree + sum(VDiat) 
 
 
  !call dgemm("T","N",This%NHL(1),3,3,1.0,This%W1,3,dG,3,0.0,dz1dR,This%NHL(1))
  dz1dR = transpose(matmul(dG,This%W1))
  
  dy1dz1 = 0.d0
  do j=1,This%NHL(1)
    dy1dz1(j,j) =  (1.d0-y1(j)**2)
  end do
  
!  call dgemm("N","N",This%NHL(1),3,This%NHL(1),1.0,dy1dz1,This%NHL(1),dz1dR,This%NHL(1),0.0,dy1dR,This%NHL(1))
!  call dgemm("T","N",This%NHL(2),3,This%NHL(1),1.0,This%W2,This%NHL(1),dy1dR,This%NHL(1),0.0,dz2dR,This%NHL(2))
  dy1dR = matmul(dy1dz1, dz1dR)
  dz2dR = matmul(transpose(This%W2), dy1dR)
 
  dVdz2 = 0.d0
  do i=1,This%NHL(2)
    dVdz2(i,i) =  This%W3(i) * (1.d0-y2(i)**2)
  end do

!  call dgemm("N","N",This%NHL(2),3,This%NHL(2),1.0,dVdz2,This%NHL(2),dz2dR,This%NHL(2),0.0,dVTemp,This%NHL(2))
  dVTemp  = matmul(dVdz2, dz2dR)
  dVTriat = VExp * sum(dVTemp,1) * eV_To_Hartree
  
  dVdQ         = Zero
  dVdR(:)      = dVTriat(This%IdxToPair(:)) + dVDiat(:)
  
  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential Calculations = ', t2-t1
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
