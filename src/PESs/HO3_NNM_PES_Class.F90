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
Module HO3_NNM_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, eV_To_Hartree, B_To_Ang
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use Input_Class           ,only:  Input_Type


  implicit none

  private
  public    :: HO3_NNM_PES_Type

  Type    ,extends(PES_Type)                ::  HO3_NNM_PES_Type
    integer     ,dimension(3)               ::  indOH
    integer     ,dimension(3)               ::  indO2
  contains
    procedure                               ::  Initialize        =>    Initialize_HO3_NNM_PES
    procedure                               ::  Output            =>    Output_HO3_NNM_PES
    procedure                               ::  Compute           =>    Compute_HO3_NNM_PES_1d
    procedure                               ::  Potential         =>    HO3_NNM_Potential_From_R
    procedure                               ::  TriatPotential    =>    HO3_NNM_Potential_From_R_OnlyTriat
  End Type

  type(Input_Type)                          ::    Input

  real(rkp)                         ,parameter    :: vpescut        = 6.0_rkp*eV_To_Hartree   ! Superior limit of PES [Hartree]
  real(rkp)                         ,parameter    :: vinflim        = -1.5_rkp*eV_To_Hartree  ! Inferior limit of PES [Hartree]
  integer                           ,parameter    :: ninputs        = 22
  integer                           ,parameter    :: nnodes2        = 50
  integer                           ,parameter    :: nnodes3        = 80
  integer                           ,parameter    :: nodemax        = nnodes3
  integer                           ,parameter    :: noutputs       = 1
  integer                           ,parameter    :: nlayer         = 4
  integer                           ,parameter    :: nscale         = ninputs + noutputs
  integer       ,dimension(4)       ,parameter    :: nodes          = [ninputs,nnodes2,nnodes3,noutputs]
  logical                           ,parameter    :: i_Debug_Global = .True.
  integer                                         :: iH
  integer       ,dimension(3)                     :: iO 
  real(rkp)     ,dimension(nscale)                :: pdela, pdelb, pdelc
  real(rkp)     ,dimension(nscale)                :: pavga, pavgb, pavgc


  real(rkp)     ,dimension(nodemax,nodemax,2:nlayer)        :: weightsa, weightsb,weightsc
  real(rkp)     ,dimension(nodemax,2:nlayer)                :: biasa, biasb, biasc


  contains
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_HO3_NNM_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                         ,only:  Input_Type
  use Atom_Class                          ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type

  class(HO3_NNM_PES_Type)                   ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug

  character(*)                    ,parameter                ::    Name_PES = 'HO3_NNM'
  integer                                                   ::    iP
  integer                                                   ::    Unit
  integer                                                   ::    ilayer
  integer                                                   ::    inode1
  integer                                                   ::    inode2
  integer                                                   ::    Status
  integer                                                   ::    iO2, iOH
  character(150)                                            ::    InputsFile
  integer         ,dimension(6,2)                           ::    iA
  type(DiatomicPotential_Factory_Type)                      ::    DiatPotFactory

  logical                                                   ::    i_Debug_Loc
  
  real(rkp)       ,dimension(3)                             ::    QH, QO1, QO2, QO3
  real(rkp)       ,dimension(12)                            ::    Q
  real(rkp)       ,dimension(6)                             ::    R
  real(rkp)                                                 ::    V


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_HO3_NNM_PES" )


  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   6
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  iA(1,:) = [1, 2]
  iA(2,:) = [1, 3]
  iA(3,:) = [1, 4]
  iA(4,:) = [2, 3]
  iA(5,:) = [2, 4]
  iA(6,:) = [3, 4]

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiatPotFactory%Construct" )
  if (i_Debug_Loc) call Logger%Write( "-> Atom 1 Name: ", Atoms(1)%Name )
  if (i_Debug_Loc) call Logger%Write( "-> Atom 2 Name: ", Atoms(2)%Name )
  if (i_Debug_Loc) call Logger%Write( "-> Atom 3 Name: ", Atoms(3)%Name )
  if (i_Debug_Loc) call Logger%Write( "-> Atom 4 Name: ", Atoms(4)%Name )
  
  iO2           = 1
  iOH           = 1
  This%indOH    = -1
  This%indO2    = -1
  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )
    if( trim(adjustl(This%Pairs(iP)%Vd%SpeciesName)) == 'O2' ) then
      This%indO2(iO2)   = iP
      iO2               = iO2+1
      if (i_Debug_Loc) call Logger%Write( "This%inO2 pair = ", iP )
    else 
      This%indOH(iOH)   = iP
      iOH               = iOH+1
      if (i_Debug_Loc) call Logger%Write( "This%inOH pair = ", iP )
    end if

  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
  ! ==============================================================================================================
  !   DETERMINING THE INDEX OF THE H ATOM AND THE O ATOMS
  ! ==============================================================================================================
  if( This%indOH(1) == 1 ) then
    if( This%indOH(2) == 2 ) then
      iH  = 1
      iO  = [2,3,4]
    else
      iH  = 2
      iO  = [1,3,4]
    end if
  end if
  if( This%indO2(1) == 1 ) then
    if( This%indO2(2) == 2 ) then
      iH  = 4
      iO  = [1,2,3]
    else
      iH  = 3
      iO  = [1,2,4]
    end if
  end if
  ! ==============================================================================================================
  !   READ PARAMETERS
  ! ==============================================================================================================
  InputsFile  = trim(adjustl(Input%DtbPath))  // '/Systems/HO3/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading HO3 NNM PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", InputsFile)
  weightsa    = Zero
  weightsb    = Zero
  weightsc    = Zero
  pdela       = Zero
  pavga       = Zero
  pdelb       = Zero
  pavgb       = Zero
  pdelc       = Zero
  pavgc       = Zero
  biasa       = Zero
  biasb       = Zero
  biasc       = Zero
  open( File=InputsFile, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // InputsFile )

    read(Unit,*) (pdela(inode1),inode1=1,nscale)
    read(Unit,*) (pavga(inode1),inode1=1,nscale)
    read(Unit,*) (pdelb(inode1),inode1=1,nscale)
    read(Unit,*) (pavgb(inode1),inode1=1,nscale)
    read(Unit,*) (pdelc(inode1),inode1=1,nscale)
    read(Unit,*) (pavgc(inode1),inode1=1,nscale)

    do ilayer = 2,nlayer
      do inode1 = 1,nodes(ilayer)
        do inode2 = 1,nodes(ilayer-1)
          read(Unit,*) weightsa(inode2,inode1,ilayer), weightsb(inode2,inode1,ilayer), weightsc(inode2,inode1,ilayer)
        end do
      end do
    end do

    do ilayer = 2,nlayer
      do inode1 = 1,nodes(ilayer)
        read(Unit,*) biasa(inode1,ilayer), biasb(inode1,ilayer), biasc(inode1,ilayer)
      end do
    end do
  close(Unit)
  ! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_HO3_NNM_PES( This, Unit )

  class(HO3_NNM_PES_Type)                 ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit

  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('N.N. HO3 PES  from U. New Mexico and Nanjing U.')")

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!
!
!
!________________________________________________________________________________________________________________________________!
Function HO3_NNM_Potential_From_R( This, R, Q ) result(V)

  class(HO3_NNM_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3)
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(6)                                     ::    inputR     !< Atom-Atom distances (bohr). Dim=(NPairs)
  integer                                                     ::    iP
  real(rkp)                                                   ::    vpesa, vpesb, vpesc
  real(rkp)                                                   ::    Vtria
  real(rkp)                                                   ::    Vdiat
  logical                                                     ::    Done
  
  Done    = .False.
  
  VDiat   = Zero
  do iP=1,6
    Vdiat = Vdiat + This%Pairs(iP)%Vd%DiatomicPotential(R(iP))
  end do
  
  if( any( R(This%indOH) <= 0.65_rkp ) .or. any( R(This%indO2) <= 0.65_rkp ) ) then
    vpesa = vpescut
    vpesb = vpescut
    vpesc = vpescut
    Vtria = vpescut
    Done  = .True.
  end if
  
  ! this input order was defined by the authors of the article in their PES code
  inputR(1) = R(This%indO2(3))
  inputR(2) = R(This%indO2(2))
  inputR(3) = R(This%indOH(3))
  inputR(4) = R(This%indO2(1))
  inputR(5) = R(This%indOH(2))
  inputR(6) = R(This%indOH(1))
  
  if( .not. Done ) then
    call ho3NN( inputR, Vtria, vpesa, vpesb, vpesc )
    if( (Vtria < vinflim ) .or. (Vtria > vpescut) ) then
      vpesa = vpescut
      vpesb = vpescut
      vpesc = vpescut
      Vtria = vpescut
    end if
  end if
  
  V       = Vtria + Vdiat

End Function
!--------------------------------------------------------------------------------------------------------------------------------!
!
!
!________________________________________________________________________________________________________________________________!
Function HO3_NNM_Potential_From_R_OnlyTriat( This, R, Q ) result(V)

  class(HO3_NNM_PES_Type)                       ,intent(in)   ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)   ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)   ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3)
  real(rkp)                                                   ::    V           !< Potential energy in [hartree].
  
  real(rkp) ,dimension(6)                                     ::    inputR      !< Atom-Atom distances (bohr). Dim=(NPairs)
  integer                                                     ::    iP
  real(rkp)                                                   ::    Vtria       !< Potential energy in [hartree]
  real(rkp)                                                   ::    vpesa, vpesb, vpesc
  logical                                                     ::    Done
  
  Done    = .False.
  
  if( any( R(This%indOH) <= 0.65_rkp ) .or. any( R(This%indO2) <= 0.65_rkp ) ) then
    vpesa = vpescut
    vpesb = vpescut
    vpesc = vpescut
    Vtria = vpescut
    Done  = .True.
  end if
  
  ! this input order was defined by the authors of the article in their PES code
  inputR(1) = R(This%indO2(3))
  inputR(2) = R(This%indO2(2))
  inputR(3) = R(This%indOH(3))
  inputR(4) = R(This%indO2(1))
  inputR(5) = R(This%indOH(2))
  inputR(6) = R(This%indOH(1))
  
  if( .not. Done ) then
    call ho3NN( inputR, Vtria, vpesa, vpesb, vpesc )
    if( (Vtria < vinflim ) .or. (Vtria > vpescut) ) then
      vpesa = vpescut
      vpesb = vpescut
      vpesc = vpescut
      Vtria = vpescut
      Done  = .True.
    end if
  end if
  
  V         = Vtria

End Function
!--------------------------------------------------------------------------------------------------------------------------------!
!
!
!________________________________________________________________________________________________________________________________!
Subroutine Compute_HO3_NNM_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(HO3_NNM_PES_Type)                       ,intent(in)   ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)   ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)   ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out)  ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out)  ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out)  ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)
  
  integer                                                     ::    iP
  real(rkp)                                                   ::    h
  real(rkp)                                                   ::    Vtria
  real(rkp) ,dimension(6)                                     ::    Rp, Rm
  real(rkp)                                                   ::    Vdiat
  real(rkp)                                                   ::    dVdiatdR
  real(rkp)                                                   ::    C
  real(rkp)                                                   ::    D
  
  V           = Zero
  dVdR        = Zero
  do iP=1,6
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), Vdiat, dVdiatdR )
    V         = Vdiat     + V
    dVdR(iP)  = dVdiatdR
  end do
  Vtria       = This%TriatPotential( R, Q )
  V           = V         + Vtria
  
!   write(*,"(a,6es15.8)") "R     = ", R
  do iP=1,6
    h         = 1.0E-7 !R(iP)*1.0E-8_rkp
    Rp        = R
    Rp(iP)    = Rp(iP)    + h
    Rm        = R
    Rm(iP)    = Rm(iP)    - h
    C         = (This%TriatPotential( Rp, Q ) -  This%TriatPotential( Rm, Q ))/(Two*h) ! Two*Vtria +
    Rp(iP)    = Rp(iP)    + h
    Rm(iP)    = Rm(iP)    - h
    D         = (This%TriatPotential( Rp, Q ) -  This%TriatPotential( Rm, Q ))/(Four*h)
    dVdR(iP)  = dVdR(iP)  + (Four*C-D)/Three
  end do
  
  write(*,"(a,6es15.8)") "Vtria = ", dVdR
!   stop
  dVdQ        = Zero

  call This%TransToCart_4Atoms( R, Q, dVdR, dVdQ)
  
End Subroutine
!________________________________________________________________________________________________________________________________!

!--------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------!
! PRIVATE ROUTINES
!--------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------------------!

Subroutine ho3NN( R, V, va, vb, vc )
  real(rkp) ,dimension(6)                       ,intent(in)   ::    R               !< Atom-Atom distances (bohr). Dim=(NPairs)
  real(rkp)                                     ,intent(out)  ::    V               !< Potential [Hartree]
  real(rkp)                                     ,intent(out)  ::    va              !< Potential component [eV]
  real(rkp)                                     ,intent(out)  ::    vb              !< Potential component [eV]
  real(rkp)                                     ,intent(out)  ::    vc              !< Potential component [eV]
  real(rkp)                 ,parameter                        ::    alfa  = 2.5_rkp !< Bohrs
  real(rkp)                 ,parameter                        ::    alfaAng = 1.322943123_rkp
  real(rkp) ,dimension(6)                                     ::    xbond
  real(rkp) ,dimension(0:ninputs)                             ::    basis
  real(rkp) ,dimension(ninputs)                               ::    xinput
  
  xbond = exp(-R/(alfa*B_To_Ang))
  basis = Zero
  
  call pip( xbond, basis )
  
  xinput(1:ninputs) = basis(1:ninputs)
  
  call getpot( xinput, pavga, pdela, weightsa, biasa, va )
  call getpot( xinput, pavgb, pdelb, weightsb, biasb, vb )
  call getpot( xinput, pavgc, pdelc, weightsc, biasc, vc )
  
  V     = ((va+vb+vc)/Three)*eV_To_Hartree
  

End Subroutine 

Subroutine pip( x, p )
  real(rkp) ,dimension(1:6)                     ,intent(in)   ::    x
  real(rkp) ,dimension(0:22)                    ,intent(out)  ::    p
  real(rkp) ,dimension(0:29)                                  ::    m

  call evmono( x, m )
  
  call evpoly( m, p )
  
End Subroutine

Subroutine evmono( x, m )
  real(rkp) ,dimension(1:6)                     ,intent(in)   ::    x
  real(rkp) ,dimension(0:29)                    ,intent(out)  ::    m
  m(0)    = 1.0
  m(1)    = x(6)
  m(2)    = x(5)
  m(3)    = x(3)
  m(4)    = x(4)
  m(5)    = x(2)
  m(6)    = x(1)
  m(7)    = m(1)*m(2)
  m(8)    = m(1)*m(3)
  m(9)    = m(2)*m(3)
  m(10)   = m(3)*m(4)
  m(11)   = m(2)*m(5)
  m(12)   = m(1)*m(6)
  m(13)   = m(4)*m(5)
  m(14)   = m(4)*m(6)
  m(15)   = m(5)*m(6)
  m(16)   = m(1)*m(9)
  m(17)   = m(1)*m(10)
  m(18)   = m(2)*m(10)
  m(19)   = m(1)*m(11)
  m(20)   = m(3)*m(11)
  m(21)   = m(2)*m(12)
  m(22)   = m(3)*m(12)
  m(23)   = m(2)*m(13)
  m(24)   = m(3)*m(13)
  m(25)   = m(1)*m(14)
  m(26)   = m(3)*m(14)
  m(27)   = m(1)*m(15)
  m(28)   = m(2)*m(15)
  m(29)   = m(4)*m(15)
End Subroutine 

Subroutine evpoly( m, p )
  real(rkp) ,dimension(0:29)                    ,intent(in)   ::    m
  real(rkp) ,dimension(0:22)                    ,intent(out)  ::    p
  p(0)    = m(0)                                          
  p(1)    = m(1) + m(2) + m(3)
  p(2)    = m(4) + m(5) + m(6)
  p(3)    = m(7) + m(8) + m(9)
  p(4)    = m(10) + m(11) + m(12)
  p(5)    = p(1)*p(2) - p(4)
  p(6)    = m(13) + m(14) + m(15)
  p(7)    = p(1)*p(1) - p(3) - p(3)
  p(8)    = p(2)*p(2) - p(6) - p(6)
  p(9)    = m(16)
  p(10)   = m(17) + m(18) + m(19) + m(20) + m(21) + m(22)
  p(11)   = p(2)*p(3) - p(10)
  p(12)   = m(23) + m(24) + m(25) + m(26) + m(27) + m(28)
  p(13)   = m(29)
  p(14)   = p(1)*p(6) - p(12)
  p(15)   = p(1)*p(3) - p(9) - p(9) - p(9)
  p(16)   = p(1)*p(4) - p(10)
  p(17)   = p(2)*p(7) - p(16)
  p(18)   = p(2)*p(4) - p(12)
  p(19)   = p(1)*p(8) - p(18)
  p(20)   = p(2)*p(6) - p(13) - p(13) - p(13)
  p(21)   = p(1)*p(7) - p(15)
  p(22)   = p(2)*p(8) - p(20) 
End Subroutine

Subroutine getpot( x, pavg, pdel, weights, bias, v )
  real(rkp)     ,dimension(ninputs)                   ,intent(in)     ::  x
  real(rkp)     ,dimension(nscale)                    ,intent(in)     ::  pavg
  real(rkp)     ,dimension(nscale)                    ,intent(in)     ::  pdel
  real(rkp)     ,dimension(nodemax,nodemax,2:nlayer)  ,intent(in)     ::  weights
  real(rkp)     ,dimension(nodemax,2:nlayer)          ,intent(in)     ::  bias
  real(rkp)                                           ,intent(out)    ::  v
  real(rkp)     ,dimension(nodemax,nlayer)                            ::  y
  integer                                                             ::  inode1, inode2, ilayer1, ilayer2
  
  y   = Zero
  
  do inode1 = 1,ninputs
    y(inode1,1)   = (x(inode1) - pavg(inode1))/pdel(inode1)
  end do
  
  do ilayer1 = 2,nlayer-1
    ilayer2 = ilayer1-1
    do inode1 = 1,nodes(ilayer1)
      y(inode1,ilayer1) = bias(inode1,ilayer1)
      do inode2 = 1,nodes(ilayer2)
        y(inode1,ilayer1) = y(inode1,ilayer1) + y(inode2,ilayer2)*weights(inode2,inode1,ilayer1)
      end do
      y(inode1,ilayer1) = tanh(y(inode1,ilayer1))
    end do
  end do
  
  ilayer1 = nlayer
  ilayer2 = ilayer1-1
  do inode1 = 1,nodes(ilayer1)
    y(inode1,ilayer1) = bias(inode1,ilayer1)
    do inode2 = 1,nodes(ilayer2)
      y(inode1,ilayer1) = y(inode1,ilayer1) + y(inode2,ilayer2)*weights(inode2,inode1,ilayer1)
    end do
  end do
  
  v   = y(nodes(nlayer),nlayer)*pdel(nscale) + pavg(nscale)
  
End Subroutine

End Module































