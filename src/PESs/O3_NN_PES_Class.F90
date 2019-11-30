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

Module O3_NN_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, eV_To_Hartree, Zero, Two
  use PES_Class             ,only:  PES_Type, DiaPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    O3_NN_PES_Type

  Type    ,extends(PES_Type)               :: O3_NN_PES_Type
  
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
  
  contains
    procedure          ::  Initialize        =>    Initialize_O3_NN_PES
    procedure          ::  Output            =>    Output_O3_NN_PES
    procedure          ::  Compute           =>    Compute_O3_NN_PES_1d
    procedure          ::  Potential         =>    O3_NN_Potential_From_R
    procedure          ::  TriatPotential    =>    O3_NN_Potential_From_R_OnlyTriat
  End Type

  logical                  ,parameter    ::    i_Debug_Global = .False.
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                    DEFERRED PROCEDURES for O3 NN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_O3_NN_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type

  class(O3_NN_PES_Type)                     ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP, i
  character(:)       ,allocatable                           ::    O3_NN_Weights_Folder
  character(:)       ,allocatable                           ::    O3_NN_Weights_File
  real(rkp)                                                 ::    LambdaTemp, reTemp        
  character(*)                    ,parameter                ::    Name_PES = 'O3_NN'
  integer                                                   ::    Status
  integer                                                   ::    Unit
  type(DiatomicPotential_Factory_Type)                       ::    DiaPotFactory
  integer         ,dimension(3,2)                           ::    iA
    
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O3_NN_PES" )
  !i_Debug_Loc   =     Logger%On()
    
  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiaPotFactory%Construct" )
  do iP = 1,This%NPairs
    call DiaPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )
  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
 ! ==============================================================================================================

  This%NHL = Input%NHiddenLayerNeurons

  write(This%NHL_Char1, '(I3)') This%NHL(1)
  write(This%NHL_Char2, '(I3)') This%NHL(2)

  allocate(This%W1(6,This%NHL(1)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%W1" )
  allocate(This%b1(This%NHL(1)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%b1" )
  allocate(This%W2(This%NHL(1),This%NHL(2)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%W2" )
  allocate(This%b2(This%NHL(2)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%b2" )
  allocate(This%W3(This%NHL(2)), stat=Status)
  if (Status/=0) call Error( "Error allocating This%W3" )

  
  O3_NN_Weights_Folder = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/PESs/NN/' // trim(adjustl(Input%PES_Model(iPES))) // '/' // trim(adjustl(This%NHL_Char1)) // '_' // trim(adjustl(This%NHL_Char2)) // '/'
  if (i_Debug_Loc) call Logger%Write( "Reading O3_NN PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening files in the folder: ", O3_NN_Weights_Folder)
  
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/ScalingValues.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    read(Unit,*)
    if (i_Debug_Loc) call Logger%Write( "Scaling_MEAN: ", This%Scaling_MEAN )
    read(Unit,*) This%Scaling_MEAN
    if (i_Debug_Loc) call Logger%Write( "Scaling_SD: ",   This%Scaling_SD )
    read(Unit,*) This%Scaling_SD
    if (i_Debug_Loc) call Logger%Write( "y_Scaling: ",    This%y_Scaling )
    read(Unit,*) This%y_Scaling
  close(Unit)
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/Lambda.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    read(Unit,*) LambdaTemp
    This%Lambda = LambdaTemp
    if (i_Debug_Loc) call Logger%Write( "Lambda: ", This%Lambda )
  close(Unit)

  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/re.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    read(Unit,*) reTemp
    This%re = reTemp
    if (i_Debug_Loc) call Logger%Write( "re: ", This%re )
  close(Unit)


  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/W1.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    do i = 1,6
      read(Unit,*) This%W1(i,1:This%NHL(1))
      if (i_Debug_Loc) call Logger%Write( "W1, Line ", i, " = ", This%W1(i,1:This%NHL(1)))
    end do
  close(Unit)
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/b1.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    do i = 1,This%NHL(1)
      read(Unit,*) This%b1(i)
    end do
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "b1 = ", This%b1)
  
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/W2.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    do i = 1,This%NHL(1)
      read(Unit,*) This%W2(i,1:This%NHL(2))
      if (i_Debug_Loc) call Logger%Write( "W2, Line ", i, " = ", This%W2(i,1:This%NHL(2)))
    end do
  close(Unit)
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/b2.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    do i = 1,This%NHL(2)
      read(Unit,*) This%b2(i)
    end do
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "b2 = ", This%b2)
  
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/W3.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    do i = 1,This%NHL(2)
      read(Unit,*) This%W3(i)
    end do
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "W3 = ", This%W3)
  
  O3_NN_Weights_File  = trim(adjustl(O3_NN_Weights_Folder)) // '/CalibratedParams/b3.csv'
  open( File=O3_NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // O3_NN_Weights_File )
    read(Unit,*) This%b3
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "b3 = ", This%b3)

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine



Subroutine Output_O3_NN_PES( This, Unit )

  class(O3_NN_PES_Type)                  ,intent(in)  ::    This
  integer                             ,intent(in)  ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O3_NN_Potential_From_R( This, R, Q ) result( V )
  
  class(O3_NN_PES_Type)                         ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(size(R,1))               ::    p          
  real(rkp) ,dimension(size(R,1))               ::    G   
  real(rkp) ,dimension(size(R))                 ::    GVecSum
  real(rkp) ,dimension(size(R))                 ::    GVecProd        
  real(rkp) ,dimension(This%NHL(1))             ::    z1
  real(rkp) ,dimension(This%NHL(1))             ::    y1
  real(rkp) ,dimension(This%NHL(2))             ::    z2
  real(rkp) ,dimension(This%NHL(2))             ::    y2
  integer                                       ::    i,iP
  
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
      
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  p = dexp( - This%Lambda * (R - This%re) )

  GVecProd = [p(1)*p(2), p(2)*p(3), p(3)*p(1)]
  GVecSum  = [p(1)+p(2), p(2)+p(3), p(3)+p(1)]

  G(1) = sum(GVecProd)
  G(2) = product(p)
  G(3) = sum(GVecProd * GVecSum)  
  G = (G - This%Scaling_MEAN) / This%Scaling_SD

  z1 = This%b1
  call dgemm("N","N",1,This%NHL(1),3,1.0,G,1,This%W1,3,1.0,z1,1)
  y1 = tanh(z1)
  z2 = This%b2
  call dgemm("N","N",1,This%NHL(2),This%NHL(1),1.0,y1,1,This%W2,This%NHL(1),1.0,z2,1)
  y2 = tanh(z2)
  !V  = ddot(This%NHL(2),y2,1,This%W3,1) + This%b3
  V  = sum(y2 * This%W3)   + This%b3
  
!  z1 = matmul(G, This%W1)  + This%b1
!  y1 = tanh(z1)
!  z2 = matmul(y1, This%W2) + This%b2
!  y2 = tanh(z2)
!  V  = sum(y2 * This%W3)   + This%b3
  
  V  = (V - This%y_Scaling ) * eV_To_Hartree + sum(VDiat) 
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O3_NN_Potential_From_R_OnlyTriat( This, R, Q ) result( V )
  
  class(O3_NN_PES_Type)                         ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(size(R,1))               ::    p          
  real(rkp) ,dimension(size(R,1))               ::    G   
  real(rkp) ,dimension(size(R))                 ::    GVecSum
  real(rkp) ,dimension(size(R))                 ::    GVecProd       
  real(rkp) ,dimension(This%NHL(1))             ::    z1
  real(rkp) ,dimension(This%NHL(1))             ::    y1
  real(rkp) ,dimension(This%NHL(2))             ::    z2
  real(rkp) ,dimension(This%NHL(2))             ::    y2
  
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
    
  p = dexp( - This%Lambda * (R - This%re) )

  GVecProd = [p(1)*p(2), p(2)*p(3), p(3)*p(1)]
  GVecSum  = [p(1)+p(2), p(2)+p(3), p(3)+p(1)]

  G(1) = sum(GVecProd)
  G(2) = product(p)
  G(3) = sum(GVecProd * GVecSum)  
  G = (G - This%Scaling_MEAN) / This%Scaling_SD
  
  z1 = This%b1
  call dgemm("N","N",1,This%NHL(1),3,1.0,G,1,This%W1,3,1.0,z1,1)
  y1 = tanh(z1)
  z2 = This%b2
  call dgemm("N","N",1,This%NHL(2),This%NHL(1),1.0,y1,1,This%W2,This%NHL(1),1.0,z2,1)
  y2 = tanh(z2)
  !V  = ddot(This%NHL(2),y2,1,This%W3,1) + This%b3
  V  = sum(y2 * This%W3)   + This%b3
  
!  z1 = matmul(G, This%W1)  + This%b1
!  y1 = tanh(z1)
!  z2 = matmul(y1, This%W2) + This%b2
!  y2 = tanh(z2)
!  V  = sum(y2 * This%W3)   + This%b3
  
  V  = (V - This%y_Scaling) * eV_To_Hartree
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_O3_NN_PES_1d( This, R, Q, V, dVdR, dVdQ )


  class(O3_NN_PES_Type)                         ,intent(in)  ::    This
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
  real(rkp) ,dimension(2*size(R))               ::    G          
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
  
  real(rkp) ,dimension(This%NHL(2),3)           ::    dVTemp
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
    
  !call cpu_time ( t1 )
  
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  p = dexp( - This%Lambda * (R - This%re) )

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
  
  V    = (V - This%y_Scaling) * eV_To_Hartree + sum(VDiat) 
 
 
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
  dVTemp = matmul(dVdz2, dz2dR)
  

  dVdQ         = Zero
  dVdR         = sum(dVTemp,1) * eV_To_Hartree + dVDiat 
  
  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential Calculations = ', t2-t1
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
