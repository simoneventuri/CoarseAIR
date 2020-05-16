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

Module BNN_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, eV_To_Hartree
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    BNN_PES_Type


  Type    ,extends(PES_Type)               :: BNN_PES_Type
  
    integer                                :: iPES
    character(5)                           :: iPES_Char
  
    integer   ,dimension(2)                :: NHL
    real(rkp) ,dimension(3)                :: Lambda
    real(rkp) ,dimension(3)                :: re
    real(rkp) ,dimension(3)                :: Scaling_MEAN
    real(rkp) ,dimension(3)                :: Scaling_SD
    real(rkp)                              :: y_Scaling
    
    real(rkp) ,dimension(:,:) ,allocatable :: W1
    real(rkp) ,dimension(:,:) ,allocatable :: W2
    real(rkp) ,dimension(:)   ,allocatable :: W3
    real(rkp) ,dimension(:)   ,allocatable :: b1
    real(rkp) ,dimension(:)   ,allocatable :: b2
    real(rkp)                              :: b3
    real(rkp)                              :: Sigma
    real(rkp)                              :: Noise
    
    real(rkp)                              :: Lambda_MEAN
    real(rkp)                              :: re_MEAN
    real(rkp) ,dimension(:,:) ,allocatable :: W1_MEAN
    real(rkp) ,dimension(:,:) ,allocatable :: W2_MEAN
    real(rkp) ,dimension(:)   ,allocatable :: W3_MEAN
    real(rkp) ,dimension(:)   ,allocatable :: b1_MEAN
    real(rkp) ,dimension(:)   ,allocatable :: b2_MEAN
    real(rkp)                              :: b3_MEAN
    real(rkp)                              :: Sigma_MEAN
    
    real(rkp)                              :: Lambda_SD
    real(rkp)                              :: re_SD
    real(rkp) ,dimension(:,:) ,allocatable :: W1_SD
    real(rkp) ,dimension(:,:) ,allocatable :: W2_SD
    real(rkp) ,dimension(:)   ,allocatable :: W3_SD
    real(rkp) ,dimension(:)   ,allocatable :: b1_SD
    real(rkp) ,dimension(:)   ,allocatable :: b2_SD
    real(rkp)                              :: b3_SD
    real(rkp)                              :: Sigma_SD
    
    character(:)              ,allocatable :: NN_Weights_Folder
  
  contains
    procedure          ::  Initialize        =>    Initialize_BNN_PES
    procedure          ::  Output            =>    Output_BNN_PES
    procedure          ::  Compute           =>    Compute_BNN_PES_1d
    procedure          ::  Potential         =>    BNN_Potential_From_R
    procedure          ::  TriatPotential    =>    BNN_Potential_From_R_OnlyTriat
    procedure          ::  SampleParamPost   =>    SampleParamPost_BNN
    procedure          ::  ReadParamPost     =>    ReadParamPost_BNN
    procedure          ::  WriteParamSample  =>    WriteParamSample_BNN
  End Type

  logical                  ,parameter    ::    i_Debug_Global = .False.
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_BNN_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(BNN_PES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms 
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP, i
  character(:)       ,allocatable                           ::    NN_Weights_Folder
  character(:)       ,allocatable                           ::    NN_Weights_File
  character(*)                    ,parameter                ::    Name_PES = 'BNN'
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(5)                                              ::    iPES_Char
  type(DiatomicPotential_Factory_Type)                       ::    DiatPotFactory
  integer         ,dimension(3,2)                           ::    iA
    
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_BNN_PES" )
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
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiatPotFactory%Construct" )
  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )
  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
 ! ==============================================================================================================
  
  This%iPES = iPES
  write(iPES_Char, "(I5)") This%iPES 
  This%iPES_Char = iPES_Char
  
  This%NHL = Input%NHiddenLayerNeurons

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


  NN_Weights_Folder = trim(adjustl(Input%DtbPath))  // '/Systems/' // trim(adjustl(Input%System)) // '/PESs/BNN/' // trim(adjustl(Input%PES_ParamsFldr(1)))
  This%NN_Weights_Folder = NN_Weights_Folder
  
  if (i_Debug_Loc) call Logger%Write( "Reading BNN PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening files in the folder: ", NN_Weights_Folder)
  
  
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


  if (Input%SampleParamsStochPES) then
    
    allocate(This%W1_MEAN(6,This%NHL(1)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%W1_MEAN" )
    allocate(This%b1_MEAN(This%NHL(1)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%b1_MEAN" )
    allocate(This%W2_MEAN(This%NHL(1),This%NHL(2)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%W2_MEAN" )
    allocate(This%b2_MEAN(This%NHL(2)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%b2_MEAN" )
    allocate(This%W3_MEAN(This%NHL(2)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%W3_MEAN" )
    
    allocate(This%W1_SD(6,This%NHL(1)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%W1_SD" )
    allocate(This%b1_SD(This%NHL(1)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%b1_SD" )
    allocate(This%W2_SD(This%NHL(1),This%NHL(2)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%W2_SD" )
    allocate(This%b2_SD(This%NHL(2)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%b2_SD" )
    allocate(This%W3_SD(This%NHL(2)), stat=Status)
    if (Status/=0) call Error( "Error allocating This%W3_SD" )
  
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Lambda_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%Lambda_MEAN
      if (i_Debug_Loc) call Logger%Write( "Lambda_MEAN: ", This%Lambda_MEAN )
    close(Unit)

    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/re_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%re_MEAN
      if (i_Debug_Loc) call Logger%Write( "re_MEAN: ", This%re_MEAN )
    close(Unit)

    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W1_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,6
        read(Unit,*) This%W1_MEAN(i,1:This%NHL(1))
        if (i_Debug_Loc) call Logger%Write( "W1_MEAN, Line ", i, " = ", This%W1_MEAN(i,1:This%NHL(1)))
      end do
    close(Unit)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b1_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(1)
        read(Unit,*) This%b1_MEAN(i)
      end do
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "b1_MEAN = ", This%b1_MEAN)
    
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W2_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(1)
        read(Unit,*) This%W2_MEAN(i,1:This%NHL(2))
        if (i_Debug_Loc) call Logger%Write( "W2_MEAN, Line ", i, " = ", This%W2_MEAN(i,1:This%NHL(2)))
      end do
    close(Unit)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b2_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(2)
        read(Unit,*) This%b2_MEAN(i)
      end do
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "b2_MEAN = ", This%b2_MEAN)
    
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W3_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(2)
        read(Unit,*) This%W3_MEAN(i)
      end do
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "W3_MEAN = ", This%W3_MEAN)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b3_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%b3_MEAN
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "b3_MEAN = ", This%b3_MEAN)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Sigma_mu.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%Sigma_MEAN
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "Sigma_MEAN = ", This%Sigma_MEAN)
    
    
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Lambda_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%Lambda_SD
      if (i_Debug_Loc) call Logger%Write( "Lambda_SD: ", This%Lambda_SD )
    close(Unit)

    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/re_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%re_SD
      if (i_Debug_Loc) call Logger%Write( "re_SD: ", This%re_SD )
    close(Unit)
    
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W1_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,6
        read(Unit,*) This%W1_SD(i,1:This%NHL(1))
        if (i_Debug_Loc) call Logger%Write( "W1_SD, Line ", i, " = ", This%W1_SD(i,1:This%NHL(1)))
      end do
    close(Unit)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b1_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(1)
        read(Unit,*) This%b1_SD(i)
      end do
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "b1_SD = ", This%b1_SD)
    
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W2_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(1)
        read(Unit,*) This%W2_SD(i,1:This%NHL(2))
        if (i_Debug_Loc) call Logger%Write( "W2_SD, Line ", i, " = ", This%W2_SD(i,1:This%NHL(2)))
      end do
    close(Unit)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b2_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(2)
        read(Unit,*) This%b2_SD(i)
      end do
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "b2_SD = ", This%b2_SD)
    
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W3_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      do i = 1,This%NHL(2)
        read(Unit,*) This%W3_SD(i)
      end do
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "W3_SD = ", This%W3_SD)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b3_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%b3_SD
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "b3_SD = ", This%b3_SD)
    
    NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Sigma_sd.csv'
    open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
      read(Unit,*) This%Sigma_SD
    close(Unit)
    if (i_Debug_Loc) call Logger%Write( "Sigma_SD = ", This%Sigma_SD)

  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine



Subroutine Output_BNN_PES( This, Unit )

  class(BNN_PES_Type)                  ,intent(in)  ::    This
  integer                              ,intent(in)  ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function BNN_Potential_From_R( This, R, Q ) result( V )
  
  class(BNN_PES_Type)                           ,intent(in)  ::    This
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
  real(rkp)                                     ::    VExp
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

  VExp = exp(V + This%Noise)
  V    = (VExp - This%y_Scaling ) * eV_To_Hartree + sum(VDiat) 
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function BNN_Potential_From_R_OnlyTriat( This, R, Q ) result( V )
  
  class(BNN_PES_Type)                           ,intent(in)  ::    This
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
  real(rkp)                                     ::    VExp
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

  VExp = exp(V + This%Noise)
  V    = (VExp - This%y_Scaling) * eV_To_Hartree
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_BNN_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(BNN_PES_Type)                           ,intent(in)  ::    This
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
  real(rkp)                                     ::    VExp
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

  
  z1   = This%b1
  !call dgemm("N","N",1,This%NHL(1),3,1.0,G,1,This%W1,3,1.0,z1,1)
  call dgemm("N","N",1,This%NHL(1),6,1.0,G,1,This%W1,6,1.0,z1,1)
  y1   = tanh(z1)
  z2   = This%b2
  call dgemm("N","N",1,This%NHL(2),This%NHL(1),1.0,y1,1,This%W2,This%NHL(1),1.0,z2,1)
  y2   = tanh(z2)
  V    = sum(y2 * This%W3) + This%b3
!  z1   = matmul(G, This%W1)  + This%b1
!  y1   = tanh(z1)
!  z2   = matmul(y1, This%W2) + This%b2
!  y2   = tanh(z2)
!  V    = sum(y2 * This%W3) + This%b3
  VExp = exp(V + This%Noise)
  V    = (VExp - This%y_Scaling) * eV_To_Hartree + sum(VDiat) 
 
 
  !call dgemm("T","N",This%NHL(1),3,3,1.0,This%W1,3,dG,3,0.0,dz1dR,This%NHL(1))
  dz1dR = transpose(matmul(dG,This%W1))
  
  dy1dz1 = 0.d0
  do j=1,This%NHL(1)
    dy1dz1(j,j) =  (1.d0-y1(j)**2)
  end do
  
  call dgemm("N","N",This%NHL(1),3,This%NHL(1),1.0,dy1dz1,This%NHL(1),dz1dR,This%NHL(1),0.0,dy1dR,This%NHL(1))
  call dgemm("T","N",This%NHL(2),3,This%NHL(1),1.0,This%W2,This%NHL(1),dy1dR,This%NHL(1),0.0,dz2dR,This%NHL(2))
!  dy1dR = matmul(dy1dz1, dz1dR)
!  dz2dR = matmul(transpose(This%W2), dy1dR)
 
  dVdz2 = 0.d0
  do i=1,This%NHL(2)
    dVdz2(i,i) =  This%W3(i) * (1.d0-y2(i)**2)
  end do

  call dgemm("N","N",This%NHL(2),3,This%NHL(2),1.0,dVdz2,This%NHL(2),dz2dR,This%NHL(2),0.0,dVTemp,This%NHL(2))
!  dVTemp = matmul(dVdz2, dz2dR)
  

  dVdQ         = Zero
  dVdR         = VExp * sum(dVTemp,1) * eV_To_Hartree + dVDiat 
  
  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential Calculations = ', t2-t1
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SampleParamPost_BNN(This)

  use RandomVector_Module   ,only:  RandNorm
  
  class(BNN_PES_Type)                       ,intent(inout)  ::    This

  integer                                                   ::    i, j

  This%Lambda = RandNorm(This%Lambda_MEAN, This%Lambda_SD) 

  This%re     = RandNorm(This%re_MEAN, This%re_SD) 
  
  do j = 1,size(This%W1_MEAN,2)
    do i= 1,size(This%W1_MEAN,1)
      This%W1(i,j) =  RandNorm(This%W1_MEAN(i,j), This%W1_SD(i,j)) !This%W1_MEAN(i,j) + fgsl_ran_gaussian(r_n_fgsl,This%W1_SD(i,j))
    end do
  end do
  
  do j = 1,size(This%W2_MEAN,2)
    do i= 1,size(This%W2_MEAN,1)
      This%W2(i,j) = RandNorm(This%W2_MEAN(i,j), This%W2_SD(i,j)) !This%W2_MEAN(i,j) + fgsl_ran_gaussian(r_n_fgsl,This%W2_SD(i,j))
    end do
  end do
  
  do i = 1,size(This%W3_MEAN,1)
    This%W3(i) = RandNorm(This%W3_MEAN(i), This%W3_SD(i)) !This%W3_MEAN(i) + fgsl_ran_gaussian(r_n_fgsl,This%W3_SD(i))
  end do
  
  
  do i = 1,size(This%b1_MEAN,1)
    This%b1(i) = RandNorm(This%b1_MEAN(i), This%b1_SD(i)) !This%b1_MEAN(i) + fgsl_ran_gaussian(r_n_fgsl,This%b1_SD(i))
  end do
  
  do i = 1,size(This%b2_MEAN,1)
    This%b2(i) = RandNorm(This%b2_MEAN(i), This%b2_SD(i)) !This%b2_MEAN(i) + fgsl_ran_gaussian(r_n_fgsl,This%b2_SD(i))
  end do
  
  This%b3    = RandNorm(This%b3_MEAN, This%b3_SD) !This%b3_MEAN + fgsl_ran_gaussian(r_n_fgsl,This%b3_SD)
  
  This%Sigma = RandNorm(This%Sigma_MEAN, This%Sigma_SD) !This%Sigma_MEAN + fgsl_ran_gaussian(r_n_fgsl,This%Sigma_SD)
  !write(*,*) 'Sigma for this PES = ', This%Sigma

  This%Noise = RandNorm(0.0d0, This%Sigma) !This%Sigma_MEAN + fgsl_ran_gaussian(r_n_fgsl,This%Sigma_SD)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadParamPost_BNN(This, PESiseed)
  
  class(BNN_PES_Type)                       ,intent(inout)  ::    This
  integer                                   ,intent(in)     ::    PESiseed
  
  integer                                                   ::    i, j
  character(:)       ,allocatable                           ::    NN_Weights_File
  character(6)                                              ::    PESiseed_Char
  real(rkp)                                                 ::    LambdaTemp, reTemp        
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)       ,allocatable                           ::    NN_Weights_Folder
  
  NN_Weights_Folder = This%NN_Weights_Folder
  
  write(PESiseed_Char, "(I6)") PESiseed

  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Lambda.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*) LambdaTemp
    This%Lambda = LambdaTemp
  close(Unit)

  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/re.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*) reTemp
    This%re = reTemp
  close(Unit)


  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W1.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,6
      read(Unit,*) This%W1(i,1:This%NHL(1))
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b1.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(1)
      read(Unit,*) This%b1(i)
    end do
  close(Unit)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W2.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(1)
      read(Unit,*) This%W2(i,1:This%NHL(2))
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b2.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(2)
      read(Unit,*) This%b2(i)
    end do
  close(Unit)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/W3.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i = 1,This%NHL(2)
      read(Unit,*) This%W3(i)
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/b3.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*) This%b3
  close(Unit)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Sigma.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*) This%Sigma
  close(Unit)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/CalibratedParams/Noise.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    read(Unit,*) This%Noise
  close(Unit)
  !This%Noise = 0.0
  !write(*,*) This%Noise
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteParamSample_BNN(This, PESiseed, LevelOutputDir)
  
  class(BNN_PES_Type)                       ,intent(in)  ::    This
  integer                                   ,intent(in)  ::    PESiseed
  character(150)                            ,intent(in)  ::    LevelOutputDir
  
  integer                                                ::    i, j
  character(3)                                           ::    Temp
  character(100)                                         ::    Str
  character(:)   ,allocatable                            ::    NN_Weights_Folder
  character(:)   ,allocatable                            ::    NN_Weights_File
  integer                                                ::    Status
  integer                                                ::    Unit
  character(6)                                           ::    PESiseed_Char

    
  NN_Weights_Folder = trim(adjustl(LevelOutputDir))  // '/StochPESParams/'
  call system('mkdir -p ' // trim(adjustl(NN_Weights_Folder)) )

  write(PESiseed_Char, "(I6)") PESiseed

  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/Lambda.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    write(Unit,"(es15.8)") This%Lambda
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/re.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    write(Unit,"(es15.8)") This%re
  close(Unit) 
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/W1.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i=1,size(This%W1,1)
      write(Unit,"(SP,es15.8,*(',',SP,es15.8:))") This%W1(i,:)
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/W2.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i=1,size(This%W2,1)
      write(Unit,"(SP,es15.8,*(',',SP,es15.8:))") This%W2(i,:)
    end do
  close(Unit) 
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/W3.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i=1,size(This%W3,1)
      write(Unit,"(es15.8)") This%W3(i)
    end do
  close(Unit)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/b1.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i=1,size(This%b1,1)
      write(Unit,"(es15.8)") This%b1(i)
    end do
  close(Unit)
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/b2.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    do i=1,size(This%b2,1)
      write(Unit,"(es15.8)") This%b2(i)
    end do
  close(Unit) 
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/b3.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    write(Unit,"(es15.8)") This%b3
  close(Unit)
  
  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/Sigma.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    write(Unit,"(es15.8)") This%Sigma
  close(Unit)

  
  NN_Weights_File  = trim(adjustl(NN_Weights_Folder)) // '/Noise.csv.' // trim(adjustl(PESiseed_Char))
  open( File=NN_Weights_File, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // NN_Weights_File )
    write(Unit,"(es15.8)") This%Noise
  close(Unit)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
