! -*-F90-*
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Berkan Bolkan and Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
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

Module O3_UMN_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, One, Two, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    O3_UMN_PES_Type


  Type    ,extends(PES_Type)          ::    O3_UMN_PES_Type
    
    real(rkp)                                                 ::    a                            ! a(in Ang)
    real(rkp)                                                 ::    ab                           ! ab (in Ang**2)
    real(rkp)                                                 ::    ra                           ! ra (in Ang)
    real(rkp)                                                 ::    rb                           ! and rb (in Ang)
    real(rkp) ,dimension(56)                                  ::     C        
  
  contains
    procedure          ::  Initialize        =>    Initialize_O3_UMN_PES
    procedure          ::  Output            =>    Output_O3_UMN_PES
    procedure          ::  Compute           =>    Compute_O3_UMN_PES_1d
    procedure          ::  Potential         =>    O3_UMN_Potential_From_R
    procedure          ::  DiatPotential     =>    O3_UMN_Potential_From_R_OnlyDiat
    procedure          ::  TriatPotential    =>    O3_UMN_Potential_From_R_OnlyTriat
    
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.
  real(rkp) ,dimension(3,9)                     ::    dRdX
                                                    
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for UMN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_O3_UMN_PES( This, Input, Atoms, iPES, i_Debug )

  use O2_UMN_DiatomicPotential_Class  ,only:  O2_UMN_DiatomicPotential_Type
  use Input_Class                     ,only:  Input_Type
  use Atom_Class                      ,only:  Atom_Type
  
  class(O3_UMN_PES_Type)                    ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP
  integer                                                   ::    k
  integer                                                   ::    Status
  integer                                                   ::    Unit
  character(:)       ,allocatable                           ::    UMN_O3_file
  character(*)                    ,parameter                ::    Name_PES = 'O3_UMN'
  character(80)                                             ::    line
  integer         ,dimension(3,2)                           ::    iA

  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O3_UMN_PES" )
  !i_Debug_Loc   =     Logger%On()
  
  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs

  ! allocate( This%mMiMn(3) )
  ! This%mMiMn(1:2) = - Atoms(1:2)%Mass / Atoms(3)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )

  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]
  
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair
  do iP = 1,This%NPairs
    allocate( O2_UMN_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
  end do
  
  
  UMN_O3_file = trim(adjustl(Input%DtbPath))  // '/Systems/' // trim(adjustl(Input%System)) // '/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading UMN O3 PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", UMN_O3_file)
  open( File=UMN_O3_file, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // UMN_O3_file )
  
  read(Unit,"(A80)") line
  if (i_Debug_Loc) call Logger%Write( "Line: ", line )

  read(Unit,*) This%a, This%ab, This%ra, This%rb
  if (i_Debug_Loc) call Logger%Write( "a  = ", This%a)
  if (i_Debug_Loc) call Logger%Write( "ab = ", This%ab )
  if (i_Debug_Loc) call Logger%Write( "ra = ", This%ra )
  if (i_Debug_Loc) call Logger%Write( "rb = ", This%rb )


  read(Unit,*) ( This%C(k), k = 1,56 )
  if (i_Debug_Loc) call Logger%Write( "C = ", This%C )

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_O3_UMN_PES( This, Unit )

  class(O3_UMN_PES_Type)                  ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('O3 PES with cta bug fixed')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O3_UMN_Potential_From_R( This, R, Q ) result( V )

  use O2_UMN_DiatomicPotential_Class   ,only: O2_UMN_DiatomicPotential_Type

  class(O3_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp)                                               ::    VDiat
  real(rkp) ,dimension(3)                                 ::    RAng  
  real(rkp) ,dimension(3)                                 ::    rMs
  real(rkp) ,dimension(0:7)                               ::    rM
  real(rkp) ,dimension(0:66)                              ::    P
  real(rkp) ,dimension(56)                                ::    B       
  integer                                                 ::    i
  real(rkp) ,dimension(3)                                 ::    Vd
  type(O2_UMN_DiatomicPotential_Type)                     ::    DiatPot      ! Diatomic potential object
  
  ! The soubroutines receives in input distances in Bohr and gives derivatives in dR
  RAng = R * B_To_Ang

  do i=1,3
    rMs(i) = exp( - (RAng(i) - This%ra)/This%a - ((RAng(i) - This%rb)**Two)/This%ab )
  end do
  call Unified( This, RAng, rMs, rM, P, B )

  ! Evaluate 2-body interactions
  Vd    = DiatPot%DiatomicPotential( R )
  VDiat = sum(Vd)

  ! Evaluate V by taken the product of C and Basis function array
  V = Zero
  do i=1,56
    V = V + This%C(i)*B(i)
  end do
  V = (V * Kcm_To_Hartree) + VDiat

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O3_UMN_Potential_From_R_OnlyDiat( This, R, Q ) result( V )

  use O2_UMN_DiatomicPotential_Class   ,only: O2_UMN_DiatomicPotential_Type

  class(O3_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].
  
  real(rkp)                                               ::    VDiat
  real(rkp) ,dimension(3)                                 ::    RAng  
  real(rkp) ,dimension(3)                                 ::    rMs
  real(rkp) ,dimension(0:7)                               ::    rM
  real(rkp) ,dimension(0:66)                              ::    P
  real(rkp) ,dimension(56)                                ::    B       
  integer                                                 ::    i
  real(rkp) ,dimension(3)                                 ::    Vd
  type(O2_UMN_DiatomicPotential_Type)                     ::    DiatPot      ! Diatomic potential object

  ! Evaluate 2-body interactions
  Vd    = DiatPot%DiatomicPotential( R )
  V     = sum(Vd)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O3_UMN_Potential_From_R_OnlyTriat( This, R, Q ) result( V )


  class(O3_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].
  
  real(rkp) ,dimension(3)                                 ::    RAng  
  real(rkp) ,dimension(3)                                 ::    rMs
  real(rkp) ,dimension(0:7)                               ::    rM
  real(rkp) ,dimension(0:66)                              ::    P
  real(rkp) ,dimension(56)                                ::    B       
  integer                                                 ::    i
  
  ! The soubroutines receives in input distances in Bohr and gives derivatives in dR
  RAng = R * B_To_Ang

  do i=1,3
    rMs(i) = exp( - (RAng(i) - This%ra)/This%a - ((RAng(i) - This%rb)**Two)/This%ab )
  end do
  call Unified( This, RAng, rMs, rM, P, B )

  ! Evaluate V by taken the product of C and Basis function array
  V = Zero
  do i=1,56
    V = V + This%C(i)*B(i)
  end do
  V = (V * Kcm_To_Hartree)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_O3_UMN_PES_1d( This, R, Q, V, dVdR, dVdQ )

  use O2_UMN_DiatomicPotential_Class   ,only: O2_UMN_DiatomicPotential_Type

  class(O3_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp)                                             ::    VDiat
  real(rkp) ,dimension(3)                               ::    Vd
  real(rkp) ,dimension(3)                               ::    dVDiat
  real(rkp) ,dimension(3)                               ::    RAng
  real(rkp) ,dimension(3)                               ::    rMs
  real(rkp) ,dimension(0:7)                             ::    rM
  real(rkp) ,dimension(0:66)                            ::    P
  real(rkp) ,dimension(56)                              ::    B       
  real(rkp) ,dimension(3,3)                             ::    dMsdR
  real(rkp) ,dimension(3,56)                            ::    dBdR           
  integer                                               ::    i, j
  type(O2_UMN_DiatomicPotential_Type)                   ::    DiatPot      ! Diatomic potential object
  
  ! The soubroutines receives in input distances in Bohr and gives derivatives in dR
  RAng = R * B_To_Ang

  do i=1,3
    rMs(i) = exp( - (RAng(i) - This%ra)/This%a - ((RAng(i) - This%rb)**Two)/This%ab )
  end do
  call Unified( This, RAng, rMs, rM, P, B )

  ! Evaluate 2-body interactions
  VDiat = Zero 
  call DiatPot%Compute_Vd_dVd( R, Vd, dVDiat )
  VDiat = sum(Vd)

  ! Evaluate V by taken the product of C and Basis function array
  V = Zero
  do i=1,56
    V = V + This%C(i)*B(i)
  end do
  V = (V * Kcm_To_Hartree) + VDiat
  
  ! Calculate dMEG/dr(3,3) for giving R(3)
  ! Calculate the monomials for each point by using six MEG terms
  ! Calculate the polynomials by using monomials
  ! Remove 2-body interactions and unconnected terms from polynomials
  ! Initialize dMsdR
  dMsdR = Zero
  do i=1,3
    dMsdR(i,i) = (-Two * (RAng(i)-This%rb)/This%ab -One/This%a) * exp( - (RAng(i) - This%ra)/This%a - ((RAng(i) - This%rb)**Two)/This%ab )
  end do 
  call UnifiedBis( This, RAng, rM, P, dMsdR, dBdR )
  
  ! Evaluate dV(3) by taken the product of C(j) and dPdR(i,j) 
  dVdR = Zero
  do j=1,56                                                                                                                       
    do i=1,3
     dVdR(i) = dVdR(i) + This%C(j)*dBdR(i,j)
    end do
  end do
  dVdR = dVdR * KcmAng_To_HartB + dVDiat

  dVdQ = Zero
  call This%TransToCart_3Atoms( R, Q, dVdR, dVdQ)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
! **************************************************************************************************************
!                                   PRIVATE PROCEDURES for UMS PES
! **************************************************************************************************************
! **************************************************************************************************************      
!________________________________________________________________________________________________________________________________!
Subroutine Unified( This, R, rMs, rM, P, B )
! mixed exponential gaussian term ms = exp(-(r-ra)/a-(r-rb)**2/ab)
! ra:   reference bond length
! rb:   reference bond length
! a:    nonlinear paramter, unit Anstrom
! ab:   nonlinear paramter, unit Anstrom**2
!   +
!  The Subroutine reads six MEG variables(X) and calculate the monomials(M) that do not have usable decomposition. For A4 with max. degree 10, the number of monomials is nom.
!   +
!  The Subroutine reads monomials(m) and calculate the permutation invariant polynomials(. p)For A4 with max. degree 10, the number of polynomials is nob.
!   +
!  The Subroutine eliminate the 2-body terms in Bowman's approach
!
  
  class(O3_UMN_PES_Type)                    ,intent(in)     ::    This
  real(rkp) ,dimension(3)                   ,intent(in)     ::    R 
  real(rkp) ,dimension(3)                   ,intent(in)     ::    rMs
  real(rkp) ,dimension(0:7)                 ,intent(out)    ::    rM
  real(rkp) ,dimension(0:66)                ,intent(out)    ::    P 
  real(rkp) ,dimension(56)                  ,intent(out)    ::    B       
    
  real(rkp) ,dimension(67)                                  ::    b1
  integer                                                   ::    i
  
  rM(0) = One
  rM(1) = rMs(3)
  rM(2) = rMs(2)
  rM(3) = rMs(1)
  rM(4) = rM(1)*rM(2)
  rM(5) = rM(1)*rM(3)
  rM(6) = rM(2)*rM(3)
  rM(7) = rM(1)*rM(6)
 
  P( 0) = rM(0)
  P( 1) = rM(1) + rM(2) + rM(3)
  P( 2) = rM(4) + rM(5) + rM(6)
  P( 3) = P( 1)*P( 1) - P( 2) - P( 2)
  P( 4) = rM(7)
  P( 5) = P( 1)*P( 2) - P( 4) - P( 4) - P( 4)
  P( 6) = P( 1)*P( 3) - P( 5)
  P( 7) = P( 1)*P( 4)
  P( 8) = P( 2)*P( 2) - P( 7) - P( 7)
  P( 9) = P( 2)*P( 3) - P( 7)
  P(10) = P( 1)*P( 6) - P( 9)
  P(11) = P( 2)*P( 4)
  P(12) = P( 3)*P( 4)
  P(13) = P( 1)*P( 8) - P(11)
  P(14) = P( 2)*P( 6) - P(12)
  P(15) = P( 1)*P(10) - P(14)
  P(16) = P( 4)*P( 4)
  P(17) = P( 4)*P( 5)
  P(18) = P( 4)*P( 6)
  P(19) = P( 2)*P( 8) - P(17)
  P(20) = P( 1)*P(13) - P(17) - P(19) - P(19)
  P(21) = P( 2)*P(10) - P(18)
  P(22) = P( 1)*P(15) - P(21)
  P(23) = P( 1)*P(16)
  P(24) = P( 4)*P( 8)
  P(25) = P( 3)*P(11) - P(23)
  P(26) = P( 4)*P(10)
  P(27) = P( 1)*P(19) - P(24)
  P(28) = P( 6)*P( 8) - P(23)
  P(29) = P( 2)*P(15) - P(26)
  P(30) = P( 1)*P(22) - P(29)
  P(31) = P( 2)*P(16)
  P(32) = P( 3)*P(16)
  P(33) = P( 1)*P(24) - P(31)
  P(34) = P( 4)*P(14)
  P(35) = P( 4)*P(15)
  P(36) = P( 2)*P(19) - P(33)
  P(37) = P( 3)*P(19) - P(31)
  P(38) = P( 8)*P(10) - P(32)
  P(39) = P( 2)*P(22) - P(35)
  P(40) = P( 1)*P(30) - P(39)
  P(41) = P( 4)*P(16)
  P(42) = P( 4)*P(17)
  P(43) = P( 4)*P(19)
  P(44) = P( 6)*P(16)
  P(45) = P( 4)*P(20)
  P(46) = P( 4)*P(21)
  P(47) = P( 4)*P(22)
  P(48) = P( 1)*P(36) - P(43)
  P(49) = P( 1)*P(37) - P(45) - P(48)
  P(50) = P( 8)*P(15) - P(44)
  P(51) = P( 2)*P(30) - P(47)
  P(52) = P( 1)*P(40) - P(51)
  P(53) = P( 1)*P(41)
  P(54) = P( 4)*P(24)
  P(55) = P( 3)*P(31) - P(53)
  P(56) = P( 1)*P(43) - P(54)
  P(57) = P(10)*P(16)
  P(58) = P( 4)*P(28)
  P(59) = P( 4)*P(29)
  P(60) = P( 4)*P(30)
  P(61) = P( 2)*P(36) - P(56)
  P(62) = P( 3)*P(36) - P(54)
  P(63) = P(10)*P(19) - P(53)
  P(64) = P( 8)*P(22) - P(57)
  P(65) = P( 2)*P(40) - P(60)
  P(66) = P( 1)*P(52) - P(65)


  ! Pass P(0:66) to BM1(1:67)
  do i=1,67
    b1(i) = P(i-1)
  end do

  ! Remove unconnected terms and 2-body terms and pass to B(1:430)
  B(1)=b1(3)

  do i=2,3
    B(i)=b1(i+3)
  end do

  do i=4,6
    B(i)=b1(i+4)
  end do

  do i=7,10
    B(i)=b1(i+5)
  end do

  do i=11,16
    B(i)=b1(i+6)
  end do

  do i=17,23
    B(i)=b1(i+7)
  end do

  do i=24,32
    B(i)=b1(i+8)
  end do

  do i=33,43
    B(i)=b1(i+9)
  end do

  do i=44,56
    B(i)=b1(i+10)
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------! 


!________________________________________________________________________________________________________________________________! 
Pure Subroutine UnifiedBis( This, R, rM, P, dMsdR, dBdR )
! Subroutine to evalute the derivatives of MEG term X
! w.r.t. interatomic distance R(3)
! dmsdR:        Local variables, dirm(3,3)
! a:            Nonlinear pamameter(Angstrom)
! ab:           Nonlinear pamameter(Angstrom**2)
! ra:           reference bond length(Angstrom)
! rb:           reference bond length(Angstrom)
!    +
!  The Subroutine reads M(nom) and dMSdR(3,3) and calculate the dMdR(3,nom) that do not have usable decomposition. For A4 with max. degree 10, the number of monomials is nom.
!    +
!  The Subroutine reads monomials(m) and calculate the permutation invariant polynomials(p) For A4 with max. degree 10, the number of polynomials is nob.
!    +
!  The Subroutine elminate the 2-body terms in Bowman's approach            
!

  class(O3_UMN_PES_Type)                    ,intent(in)     ::    This
  real(rkp) ,dimension(3)                   ,intent(in)     ::    R         ! Distances of atom-atom pairs [Angstrom]
  real(rkp) ,dimension(0:7)                 ,intent(in)     ::    rM
  real(rkp) ,dimension(0:66)                ,intent(in)     ::    P
  real(rkp) ,dimension(3,3)                 ,intent(in)     ::    dMsdR
  real(rkp) ,dimension(3,56)                ,intent(out)    ::    dBdR           
  
  real(rkp) ,dimension(3,67)                                ::    db1dr
  integer                                                   ::    i ,j
  real(rkp) ,dimension(3,0:7)                               ::    dMdR 
  real(rkp) ,dimension(3,0:66)                              ::    dPdR
                                             
  do i=1,3
    dMdR(i,0) = Zero
    dMdR(i,1) = dMsdR(i,3)
    dMdR(i,2) = dMsdR(i,2)
    dMdR(i,3) = dMsdR(i,1)
    dMdR(i,4) = dMdR(i,1)*rM(2) + rM(1)*dMdR(i,2)
    dMdR(i,5) = dMdR(i,1)*rM(3) + rM(1)*dMdR(i,3)
    dMdR(i,6) = dMdR(i,2)*rM(3) + rM(2)*dMdR(i,3)
    dMdR(i,7) = dMdR(i,1)*rM(6) + rM(1)*dMdR(i,6)
  end do
  
  do i=1,3
    dPdR(i, 0) = dMdR(i,0)
    dPdR(i, 1) = dMdR(i,1) + dMdR(i,2) + dMdR(i,3)
    dPdR(i, 2) = dMdR(i,4) + dMdR(i,5) + dMdR(i,6)
    dPdR(i, 3) = dPdR(i, 1)*P( 1) + P( 1)*dPdR(i, 1) - dPdR(i, 2) - dPdR(i, 2)
    dPdR(i, 4) = dMdR(i,7)
    dPdR(i, 5) = dPdR(i, 1)*P( 2) + P( 1)*dPdR(i, 2) - dPdR(i, 4) - dPdR(i, 4) - dPdR(i, 4)
    dPdR(i, 6) = dPdR(i, 1)*P( 3) + P( 1)*dPdR(i, 3) - dPdR(i, 5)
    dPdR(i, 7) = dPdR(i, 1)*P( 4) + P( 1)*dPdR(i, 4)
    dPdR(i, 8) = dPdR(i, 2)*P( 2) + P( 2)*dPdR(i, 2) - dPdR(i, 7) - dPdR(i, 7)
    dPdR(i, 9) = dPdR(i, 2)*P( 3) + P( 2)*dPdR(i, 3) - dPdR(i, 7)
    dPdR(i,10) = dPdR(i, 1)*P( 6) + P( 1)*dPdR(i, 6) - dPdR(i, 9)
    dPdR(i,11) = dPdR(i, 2)*P( 4) + P( 2)*dPdR(i, 4)
    dPdR(i,12) = dPdR(i, 3)*P( 4) + P( 3)*dPdR(i, 4)
    dPdR(i,13) = dPdR(i, 1)*P( 8) + P( 1)*dPdR(i, 8) - dPdR(i,11)
    dPdR(i,14) = dPdR(i, 2)*P( 6) + P( 2)*dPdR(i, 6) - dPdR(i,12)
    dPdR(i,15) = dPdR(i, 1)*P(10) + P( 1)*dPdR(i,10) - dPdR(i,14)
    dPdR(i,16) = dPdR(i, 4)*P( 4) + P( 4)*dPdR(i, 4)
    dPdR(i,17) = dPdR(i, 4)*P( 5) + P( 4)*dPdR(i, 5)
    dPdR(i,18) = dPdR(i, 4)*P( 6) + P( 4)*dPdR(i, 6)
    dPdR(i,19) = dPdR(i, 2)*P( 8) + P( 2)*dPdR(i, 8) - dPdR(i,17)
    dPdR(i,20) = dPdR(i, 1)*P(13) + P( 1)*dPdR(i,13) - dPdR(i,17) - dPdR(i,19) - dPdR(i,19)
    dPdR(i,21) = dPdR(i, 2)*P(10) + P( 2)*dPdR(i,10) - dPdR(i,18)
    dPdR(i,22) = dPdR(i, 1)*P(15) + P( 1)*dPdR(i,15) - dPdR(i,21)
    dPdR(i,23) = dPdR(i, 1)*P(16) + P( 1)*dPdR(i,16)
    dPdR(i,24) = dPdR(i, 4)*P( 8) + P( 4)*dPdR(i, 8)
    dPdR(i,25) = dPdR(i, 3)*P(11) + P( 3)*dPdR(i,11) - dPdR(i,23)
    dPdR(i,26) = dPdR(i, 4)*P(10) + P( 4)*dPdR(i,10)
    dPdR(i,27) = dPdR(i, 1)*P(19) + P( 1)*dPdR(i,19) - dPdR(i,24)
    dPdR(i,28) = dPdR(i, 6)*P( 8) + P( 6)*dPdR(i, 8) - dPdR(i,23)
    dPdR(i,29) = dPdR(i, 2)*P(15) + P( 2)*dPdR(i,15) - dPdR(i,26)
    dPdR(i,30) = dPdR(i, 1)*P(22) + P( 1)*dPdR(i,22) - dPdR(i,29)
    dPdR(i,31) = dPdR(i, 2)*P(16) + P( 2)*dPdR(i,16)
    dPdR(i,32) = dPdR(i, 3)*P(16) + P( 3)*dPdR(i,16)
    dPdR(i,33) = dPdR(i, 1)*P(24) + P( 1)*dPdR(i,24) - dPdR(i,31)
    dPdR(i,34) = dPdR(i, 4)*P(14) + P( 4)*dPdR(i,14)
    dPdR(i,35) = dPdR(i, 4)*P(15) + P( 4)*dPdR(i,15)
    dPdR(i,36) = dPdR(i, 2)*P(19) + P( 2)*dPdR(i,19) - dPdR(i,33)
    dPdR(i,37) = dPdR(i, 3)*P(19) + P( 3)*dPdR(i,19) - dPdR(i,31)
    dPdR(i,38) = dPdR(i, 8)*P(10) + P( 8)*dPdR(i,10) - dPdR(i,32)
    dPdR(i,39) = dPdR(i, 2)*P(22) + P( 2)*dPdR(i,22) - dPdR(i,35)
    dPdR(i,40) = dPdR(i, 1)*P(30) + P( 1)*dPdR(i,30) - dPdR(i,39)
    dPdR(i,41) = dPdR(i, 4)*P(16) + P( 4)*dPdR(i,16)
    dPdR(i,42) = dPdR(i, 4)*P(17) + P( 4)*dPdR(i,17)
    dPdR(i,43) = dPdR(i, 4)*P(19) + P( 4)*dPdR(i,19)
    dPdR(i,44) = dPdR(i, 6)*P(16) + P( 6)*dPdR(i,16)
    dPdR(i,45) = dPdR(i, 4)*P(20) + P( 4)*dPdR(i,20)
    dPdR(i,46) = dPdR(i, 4)*P(21) + P( 4)*dPdR(i,21)
    dPdR(i,47) = dPdR(i, 4)*P(22) + P( 4)*dPdR(i,22)
    dPdR(i,48) = dPdR(i, 1)*P(36) + P( 1)*dPdR(i,36) - dPdR(i,43)
    dPdR(i,49) = dPdR(i, 1)*P(37) + P( 1)*dPdR(i,37) - dPdR(i,45) - dPdR(i,48)
    dPdR(i,50) = dPdR(i, 8)*P(15) + P( 8)*dPdR(i,15) - dPdR(i,44)
    dPdR(i,51) = dPdR(i, 2)*P(30) + P( 2)*dPdR(i,30) - dPdR(i,47)
    dPdR(i,52) = dPdR(i, 1)*P(40) + P( 1)*dPdR(i,40) - dPdR(i,51)
    dPdR(i,53) = dPdR(i, 1)*P(41) + P( 1)*dPdR(i,41)
    dPdR(i,54) = dPdR(i, 4)*P(24) + P( 4)*dPdR(i,24)
    dPdR(i,55) = dPdR(i, 3)*P(31) + P( 3)*dPdR(i,31) - dPdR(i,53)
    dPdR(i,56) = dPdR(i, 1)*P(43) + P( 1)*dPdR(i,43) - dPdR(i,54)
    dPdR(i,57) = dPdR(i,10)*P(16) + P(10)*dPdR(i,16)
    dPdR(i,58) = dPdR(i, 4)*P(28) + P( 4)*dPdR(i,28)
    dPdR(i,59) = dPdR(i, 4)*P(29) + P( 4)*dPdR(i,29)
    dPdR(i,60) = dPdR(i, 4)*P(30) + P( 4)*dPdR(i,30)
    dPdR(i,61) = dPdR(i, 2)*P(36) + P( 2)*dPdR(i,36) - dPdR(i,56)
    dPdR(i,62) = dPdR(i, 3)*P(36) + P( 3)*dPdR(i,36) - dPdR(i,54)
    dPdR(i,63) = dPdR(i,10)*P(19) + P(10)*dPdR(i,19) - dPdR(i,53)
    dPdR(i,64) = dPdR(i, 8)*P(22) + P( 8)*dPdR(i,22) - dPdR(i,57)
    dPdR(i,65) = dPdR(i, 2)*P(40) + P( 2)*dPdR(i,40) - dPdR(i,60)
    dPdR(i,66) = dPdR(i, 1)*P(52) + P( 1)*dPdR(i,52) - dPdR(i,65)
  end do
  
  
  do i=1,67                                                                                                                       
    do j=1,3
      db1dr(j,i) = dPdR(j,i-1)
    end do
  end do

  ! Remove unconnected terms and 2-body terms and pass to B(1:56)
  do j=1,3

    dBdR(j,1)=db1dr(j,3)

    do i=2,3
      dBdR(j,i)=db1dr(j,i+3)
    end do

    do i=4,6
      dBdR(j,i)=db1dr(j,i+4)
    end do

    do i=7,10
      dBdR(j,i)=db1dr(j,i+5)
    end do

    do i=11,16
      dBdR(j,i)=db1dr(j,i+6)
    end do

    do i=17,23
      dBdR(j,i)=db1dr(j,i+7)
    end do

    do i=24,32
      dBdR(j,i)=db1dr(j,i+8)
    end do

    do i=33,43
      dBdR(j,i)=db1dr(j,i+9)
    end do

    do i=44,56
      dBdR(j,i)=db1dr(j,i+10)
    end do

  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!      

End Module
