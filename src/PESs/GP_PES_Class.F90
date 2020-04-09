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

Module GP_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, eV_To_Hartree
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    GP_PES_Type


  Type    ,extends(PES_Type)               :: GP_PES_Type
  
    real(rkp) ,dimension(3)                :: Scaling_MEAN
    real(rkp) ,dimension(3)                :: Scaling_SD  
  
    integer                                :: NIdxPnts
    integer                                :: NLengthsKernel
    real(rkp) ,dimension(:,:) ,allocatable :: ObsIdxPnts
    real(rkp) ,dimension(:,:) ,allocatable :: LChol
    real(rkp) ,dimension(:)   ,allocatable :: Alpha
    real(rkp) ,dimension(:)   ,allocatable :: LengthKernel
    real(rkp) ,dimension(:)   ,allocatable :: LengthKernelInvSqr
    real(rkp) ,dimension(:)   ,allocatable :: ModPip
    real(rkp)                              :: AmpKernel
    real(rkp)                              :: AmpKernelSqreV
    real(rkp)                              :: y_Scaling
    character(3)                           :: PIP

  contains
    procedure          ::  Initialize        =>    Initialize_GP_PES
    procedure          ::  Output            =>    Output_GP_PES
    procedure          ::  Compute           =>    Compute_GP_PES_1d
    procedure          ::  Potential         =>    GP_Potential_From_R
    procedure          ::  TriatPotential    =>    GP_Potential_From_R_OnlyTriat
  End Type

  logical                  ,parameter    ::    i_Debug_Global = .False.
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_GP_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(GP_PES_Type)                        ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms 
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP, i, iLines
  character(:)       ,allocatable                           ::    GP_Weights_Folder
  character(:)       ,allocatable                           ::    GP_File
  character(*)                    ,parameter                ::    Name_PES = 'GP'
  integer                                                   ::    Status
  integer                                                   ::    Unit
  type(DiatomicPotential_Factory_Type)                       ::    DiatPotFactory
  integer         ,dimension(3,2)                           ::    iA

  logical                                                   ::    i_Debug_Loc
    
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_GP_PES" )
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

 
  This%PIP = 'DEF'


  GP_Weights_Folder = trim(adjustl(Input%DtbPath))  // '/Systems/' // trim(adjustl(Input%System)) // '/PESs/GP/' // trim(adjustl(Input%PES_Model(iPES))) // '/' // trim(adjustl(Input%PES_ParamsFldr(1)))
  if (i_Debug_Loc) call Logger%Write( "Reading O3 GP PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening files in the folder: ", GP_Weights_Folder)
  
   
  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/ScalingValues.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // GP_File )
    read(Unit,*)
    read(Unit,*) This%Scaling_MEAN
    if (i_Debug_Loc) call Logger%Write( "Scaling_MEAN: ", This%Scaling_MEAN )
    read(Unit,*) This%Scaling_SD
    if (i_Debug_Loc) call Logger%Write( "Scaling_SD: ", This%Scaling_SD )
  close(Unit)


  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/ObsIdxPnts.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) then
      call Error( "Error opening file: " // GP_File )
    else
      !read(Unit,*,iostat=Status)
      iLines = 0
      do while (Status==0)
        read(Unit,*,iostat=Status)
        iLines = iLines + 1
      end do
      iLines = iLines - 1
    end if
    This%NIdxPnts = iLines
    if (i_Debug_Loc) call Logger%Write( "Found This%NIdxPnts = ", This%NIdxPnts, " Index Points" )
    
    allocate(This%ObsIdxPnts(3,This%NIdxPnts), stat=Status)
    if (Status/=0) call Error( "Error allocating This%ObsIdxPnts" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%ObsIdxPnts with dimension (3,", This%NIdxPnts,")" )
    
    allocate(This%LChol(This%NIdxPnts,This%NIdxPnts), stat=Status)
    if (Status/=0) call Error( "Error allocating This%LChol" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%LChol with dimension (", This%NIdxPnts, ",", This%NIdxPnts, ")" )
    
    allocate(This%ModPip(6), stat=Status)
    if (Status/=0) call Error( "Error allocating This%ModPip" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%ModPip with dimension (6)")

    allocate(This%Alpha(This%NIdxPnts), stat=Status)
    if (Status/=0) call Error( "Error allocating This%Alpha" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%Alpha with dimension (", This%NIdxPnts, ")" )
    
    write(*,*) 'Nidx = ', This%NIdxPnts

  rewind(Unit)
  
    !read(Unit,*,iostat=Status)
    do i = 1,This%NIdxPnts
      read(Unit,*) This%ObsIdxPnts(1:3,i)
    end do
    if (i_Debug_Loc) call Logger%Write( "This%ObsIdxPnts = ", This%ObsIdxPnts )
    
    !WRITE(*,*) "Is arg contiguous? ", is_contiguous( This%ObsIdxPnts )

    
  close(Unit)
  
  
  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/LChol.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // GP_File )
    do i = 1,This%NIdxPnts
      read(Unit,*) This%LChol(1:This%NIdxPnts,i)
    end do
  if (i_Debug_Loc) call Logger%Write("This%LChol", This%LChol)
  close(Unit)
  
  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/Alpha.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // GP_File )
    !read(Unit,*)
    do i = 1,This%NIdxPnts
      read(Unit,*) This%Alpha(i)
    end do
    if (i_Debug_Loc) call Logger%Write( "This%Alpha = ", This%Alpha)
  close(Unit)
  
    GP_File  = trim(adjustl(GP_Weights_Folder)) // '/ModPip.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // GP_File )
    !read(Unit,*)
    do i = 1,6
      read(Unit,*) This%ModPip(i)
    end do
    This%y_Scaling = This%ModPip(6)
    if (i_Debug_Loc) call Logger%Write( "This%ModPip = ", This%ModPip)
  close(Unit)
  
  
!  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/SigmaNoise.csv'
!  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
!    if (Status/=0) call Error( "Error opening file: " // GP_File )
!    read(Unit,*)
!    read(Unit,*) This%SigmaNoise
!    if (i_Debug_Loc) call Logger%Write( "This%SigmaNoise = ", This%SigmaNoise)
!  close(Unit)
  
  
  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/AmpKernel.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // GP_File )
    !read(Unit,*)
    read(Unit,*) This%AmpKernel
    if (i_Debug_Loc) call Logger%Write( "This%AmpKernel = ", This%AmpKernel)
  close(Unit)
  This%AmpKernelSqreV = This%AmpKernel**2
  
  
  GP_File  = trim(adjustl(GP_Weights_Folder)) // '/LengthKernel.csv'
  open( File=GP_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) then
      call Error( "Error opening file: " // GP_File )
    else
      !read(Unit,*,iostat=Status)
      iLines = 0
      do while (Status==0)
        read(Unit,*,iostat=Status)
        iLines = iLines + 1
      end do
    end if
    iLines = iLines - 1
    This%NLengthsKernel = iLines
    if (i_Debug_Loc) call Logger%Write( "Found This%NLengthsKernel = ", This%NLengthsKernel, " Kernel Length Scales" )
    
    allocate(This%LengthKernel(This%NLengthsKernel), stat=Status)
    if (Status/=0) call Error( "Error allocating This%LengthKernel" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%LengthKernel with dimension (", This%NLengthsKernel, ")" )
    
    allocate(This%LengthKernelInvSqr(This%NLengthsKernel), stat=Status)
    if (Status/=0) call Error( "Error allocating This%LengthKernelInvSqr" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%LengthKernelInvSqr with dimension (", This%NLengthsKernel, ")" )

  rewind(Unit)
  
    !read(Unit,*,iostat=Status)
    do i = 1,This%NLengthsKernel
      read(Unit,*) This%LengthKernel(i)
    end do
    write(*,*) This%LengthKernel
    if (i_Debug_Loc) call Logger%Write( "This%LengthKernel = ", This%LengthKernel )
    This%LengthKernelInvSqr = 5.d-1 / This%LengthKernel**2
    if (i_Debug_Loc) call Logger%Write( "This%LengthKernelInvSqr = ", This%LengthKernelInvSqr )
     
  close(Unit)
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine



Subroutine Output_GP_PES( This, Unit )

  class(GP_PES_Type)                   ,intent(in)  ::    This
  integer                                 ,intent(in)  ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  !write(Unit,"('O3 PES with cta bug fixed')")
  !write(Unit,"('acpf energies, 1s2s core and leroy n2 with hyper radius repulsion')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function GP_Potential_From_R( This, R, Q ) result( V )

  class(GP_PES_Type)                            ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(size(R,1))               ::    p          
  real(rkp) ,dimension(size(R,1))               ::    G   
  real(rkp)                                     ::    KSum
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
  real(rkp)                         ,parameter  ::    ebn2      =   0.3554625704d0
  integer                                       ::    iP
  
  !real(rkp)                         ,parameter  ::    RCut      =   9.90d0
  real(rkp)                                     ::    t1, t2
  
  !call cpu_time ( t1 )
  
  do iP=1,3
    VDiat = This%Pairs(iP)%Vd%DiatomicPotential(R(iP))
  end do 
  
  V  = 0.d0
    
  p = exp( - This%ModPip(1)* R  ** (This%ModPip(2)) )
  
  if (This%PIP ==      'OLD') then
    G(1) = ( p(1)        + p(2)        + p(3)                      ) ** (This%ModPip(3))
  elseif (This%PIP ==  'NEW') then
    G(1) = ( ( p(1) + p(2) + p(3) )*( p(1) * p(2) * p(3) )         ) ** (This%ModPip(3))
  elseif (This%PIP ==  'DEF') then
    G(1) = ( (p(1)*p(2) + p(2)*p(3) + p(3)*p(1))*(p(1)+p(2)+p(3))  ) ** (This%ModPip(3))
  end if
  G(2) = ( p(1) * p(2) + p(2) * p(3) + p(3) * p(1)                 ) ** (This%ModPip(4))
  G(3) = ( p(1) * p(2) * p(3)                                      ) ** (This%ModPip(5))
  G    = (G - This%Scaling_MEAN) / This%Scaling_SD

  call ComputeMean(G, This%NIdxPnts, This%LengthKernelInvSqr, This%Alpha, This%ObsIdxPnts, KSum)
  
  V = KSum * This%AmpKernelSqreV + sum(VDiat) + ebn2
  
  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential Calculations = ', t2-t1
    
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function GP_Potential_From_R_OnlyTriat( This, R, Q ) result( V )


  class(GP_PES_Type)                            ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(size(R,1))               ::    p          
  real(rkp) ,dimension(size(R,1))               ::    G   
  real(rkp) ,dimension(This%NIdxPnts)           ::    K
  real(rkp)                                     ::    KSum
  
  real(rkp)                                     ::    t1, t2
  
  !call cpu_time ( t1 )
  
  p = exp( - This%ModPip(1)* R  ** (This%ModPip(2)) )
  
  if (This%PIP ==      'OLD') then
    G(1) = ( p(1)        + p(2)        + p(3)                     ) ** (This%ModPip(3))
  elseif (This%PIP ==  'NEW') then
    G(1) = ( ( p(1) + p(2) + p(3) )*( p(1) * p(2) * p(3) )        ) ** (This%ModPip(3))
  elseif (This%PIP ==  'DEF') then
    G(1) = ( (p(1)*p(2) + p(2)*p(3) + p(3)*p(1))*(p(1)+p(2)+p(3)) ) ** (This%ModPip(3))
  end if
  G(2) = ( p(1) * p(2) + p(2) * p(3) + p(3) * p(1)                ) ** (This%ModPip(4))
  G(3) = ( p(1) * p(2) * p(3)                                     ) ** (This%ModPip(5))
  G    = (G - This%Scaling_MEAN) / This%Scaling_SD

  call ComputeMean(G, This%NIdxPnts, This%LengthKernelInvSqr, This%Alpha, This%ObsIdxPnts, KSum)

  V = KSum * This%AmpKernelSqreV
  
  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential Calculations = ', t2-t1
    
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_GP_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(GP_PES_Type)                            ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp) ,dimension(3)                       ::    dVFin
  real(rkp) ,dimension(size(R,1))               ::    p          
  real(rkp) ,dimension(size(R,1))               ::    G   
  real(rkp) ,dimension(3,3)                     ::    dG
  
  real(rkp)                                     ::    t1, t2, t3, t4
  real(rkp) ,dimension(3)                       ::    VDiat          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp) ,dimension(3)                       ::    dVDiat
  real(rkp)                                     ::    KSum
  real(rkp) ,dimension(3)                       ::    KDerSum
  real(rkp) ,dimension(3,3)                     ::    Var
  real(rkp) ,dimension(3)                       ::    R_el
  integer                                       ::    iP
  real(rkp)                                     ::    Temp1, Temp2, Temp3, GG1, GG2, GG3
  real(rkp) ,dimension(3)                       ::    TempVec1, TempVec2
  
  call cpu_time ( t1 )
    
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do 
  
  p = exp( - This%ModPip(1)* R  ** (This%ModPip(2)) )
  
  TempVec1 = [p(2)*p(3), p(3)*p(1), p(1)*p(2)]
  TempVec2 = [p(2)+p(3), p(3)+p(1), p(1)+p(2)]

  Temp1    = sum(TempVec2) / 2.d0
  Temp2    = p(1) * p(2) * p(3)
  Temp3    = sum(TempVec1)

  if (This%PIP ==      'OLD') then
    G(1)  = Temp1**This%ModPip(3)
  elseif (This%PIP ==  'NEW') then
    G(1)  = (Temp1*Temp2)**This%ModPip(3)
  elseif (This%PIP ==  'DEF') then  
    G(1)  = (Temp1*Temp3)**This%ModPip(3)
  end if
  G(2) = Temp3**This%ModPip(4)
  G(3) = Temp2**This%ModPip(5)
  

  R_el = - R(1:3) ** (This%ModPip(2)-1) * (This%ModPip(1) * This%ModPip(2))
  
  if (This%PIP ==      'OLD') then
    dG(1,1:3) = p                              * R_el * This%ModPip(3) * G(1) / Temp1         / This%Scaling_SD(1)
  elseif (This%PIP ==  'NEW') then
    dG(1,1:3) = (p + Temp1)                    * R_el * This%ModPip(3) * G(1) / Temp1         / This%Scaling_SD(1)
  elseif (This%PIP ==  'DEF') then
    dG(1,1:3) = (TempVec2 * Temp1 + Temp3) * p * R_el * This%ModPip(3) * G(1) / (Temp1*Temp3) / This%Scaling_SD(1)
  end if
  dG(2,1:3) = TempVec2 * p                     * R_el * This%ModPip(4) * G(2) / Temp3         / This%Scaling_SD(2)
  dG(3,1:3) =                                    R_el * This%ModPip(5) * G(3)                 / This%Scaling_SD(3)


  ! GG2  = p(1) * p(2) + p(2) * p(3) + p(3) * p(1)
  ! GG1  = p(1)        + p(2)        + p(3)

  ! dG(1,1:3) = p(1:3) * R_el(1:3) * [GG2 + GG1 * ( p(2) + p(3) ), GG2 + GG1 * ( p(1) + p(3) ), GG2 + GG1 * ( p(2) + p(1) )] * This%ModPip(3) * G(1) / ( GG1 * GG2)            /This%Scaling_SD(1)
  ! dG(2,1:3) = p(1:3) * R_el(1:3) * [p(2) + p(3), p(1) + p(3), p(2) + p(1)]                                                 * This%ModPip(4) * G(2) / ( GG2 )                 /This%Scaling_SD(2)
  ! dG(3,1:3) = R_el(1:3)                                                                                                    * This%ModPip(5) * G(3)                           /This%Scaling_SD(3)
  

  !Normalized later in order to use G in dG formulation
  G    = (G - This%Scaling_MEAN) / This%Scaling_SD

  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential CalculationsA = ', t2-t1
  
  call Compute_Mean_Var_dV(G, This%NIdxPnts, This%LengthKernelInvSqr, This%LChol, This%Alpha, This%ObsIdxPnts, This%AmpKernelSqreV, KSum, KDerSum, Var)
  
  !call cpu_time ( t3 )
  !write(*,*) 'Time for Potential CalculationsB = ', t3-t2
  
  !exp(Ksum) from Compute_Mean_Var_dV
  V     = (KSum - This%y_Scaling)*eV_To_Hartree + sum(VDiat)  
  

  dVdQ         = Zero
  dVdR         = dVDiat + matmul(KDerSum, dG)*KSum*eV_To_Hartree
  
  !call cpu_time ( t4 )
  !write(*,*) 'Time for Potential CalculationsC = ', t4-t3
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
! **************************************************************************************************************
!                                   PRIVATE PROCEDURES for GP PES
! **************************************************************************************************************
! **************************************************************************************************************
Pure Subroutine ComputeMean(G, NIdxPnts, LengthKernelInvSqr, Alpha, ObsIdxPnts, KSum)
  
  real(rkp) ,dimension(3)          ,intent(in)  ::    G   
  integer                          ,intent(in)  ::    NIdxPnts
  real(rkp) ,dimension(3)          ,intent(in)  ::    LengthKernelInvSqr
  real(rkp) ,dimension(NIdxPnts)   ,intent(in)  ::    Alpha
  real(rkp) ,dimension(3,NIdxPnts) ,intent(in)  ::    ObsIdxPnts
  real(rkp)                        ,intent(out) ::    KSum
  
  real(rkp)                                     ::    Ki
  real(rkp) ,dimension(3)                       ::    DeltaX
  real(rkp)                                     ::    Temp
  
  integer                                       ::    i

  KSum    = 0.d0
  do i = 1,NIdxPnts
    Temp = Alpha(i)
    
    DeltaX = G - ObsIdxPnts(:,i)
    
    Ki     = 1.d0 / (1.d0 + sum(DeltaX**2 * LengthKernelInvSqr) ) 
          
    KSum    = KSum    + Ki    *  Temp
    
  end do

end Subroutine


Subroutine Compute_Mean_Var_dV(G, NIdxPnts, LengthKernelInvSqr, LChol, Alpha, ObsIdxPnts, AmpKernelSqreV, KSum, KDerSum, Var)
  
  real(rkp) ,dimension(3)                        ,intent(in)  ::    G   
  integer                                        ,intent(in)  ::    NIdxPnts
  real(rkp) ,dimension(3)                        ,intent(in)  ::    LengthKernelInvSqr
  real(rkp) ,dimension(NIdxPnts,NIdxPnts)        ,intent(in)  ::    LChol
  real(rkp) ,dimension(NIdxPnts)                 ,intent(in)  ::    Alpha
  real(rkp) ,dimension(3,NIdxPnts)               ,intent(in)  ::    ObsIdxPnts
  real(rkp)                                      ,intent(in)  ::    AmpKernelSqreV
  real(rkp)                                      ,intent(out) ::    KSum
  real(rkp) ,dimension(3)                        ,intent(out) ::    KDerSum
  real(rkp) ,dimension(3,3)                      ,intent(out) ::    Var
  
  real(rkp)                                      ::    Ki
  real(rkp) ,dimension(NIdxPnts)                 ::    KK
  real(rkp) ,dimension(3,NIdxPnts)               ::    KDer
  real(rkp) ,dimension(NIdxPnts,3)               ::    KDer_1
  real(rkp) ,dimension(3)                        ::    DeltaX
  real(rkp)                                      ::    Temp
  
  integer                                        ::    i, p, j, k

  KSum    = 0.d0
  KDerSum = 0.d0
  
  do i = 1,NIdxPnts
    
    DeltaX     = G - ObsIdxPnts(:,i)
    
    Ki         = 1 / (1.d0 + sum(DeltaX**2 * LengthKernelInvSqr) )

    KDer(:,i)  = 2.d0 * LengthKernelInvSqr * DeltaX * Ki ** 2

    KK(i)      = Ki
    
  end do 
  
  KDerSum = - AmpKernelSqreV * matmul(KDer, Alpha)

  KSum    =   exp(AmpKernelSqreV * dot_product(Alpha, KK))
  
  Var     = 0.d0

  !do p=1,3

  !  Var(p,p) = AmpKernelSqreV*LengthKernelInvSqr(p)
  
  !end do

  !call dgemm("N","T",NIdxPnts,3,NIdxPnts,1.d0,LChol,NIdxPnts,KDer,3,0.d0,KDer_1,NIdxPnts)

  !KDer_1    = matmul(LChol,transpose(KDer))
  !call dgemm("T","N",3,3,NIdxPnts,1.d0,KDer1,NIdxPnts,KDer1,NIdxPnts,0.d0,Var,3)
  
  !Var = Var - AmpKernelSqreV ** 2 * matmul(transpose(KDer_1),KDer_1)        

end Subroutine



End Module
