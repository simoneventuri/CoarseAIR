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

Module CO2_NASA_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Four, Six
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    CO2_NASA_PES_Type

  Type    ,extends(PES_Type)          ::    CO2_NASA_PES_Type
    integer                                                      :: idc = 128
    integer                                                      :: idx = 100
  
    integer                                                      :: iO2
    integer                                                      :: iCO
    integer                                                      :: jCO
    
    character(80)                                                :: line
    integer                                                      :: npowco
    integer                                                      :: nfitco
    integer                                                      :: nparh
    integer                                                      :: nfcnt
    integer                                                      :: kdum
    real(rkp)                                                    :: c6o2, c6co
    real(rkp)                                                    :: ao2
    real(rkp)                                                    :: aco
    real(rkp)                                                    :: r0co
    real(rkp)                                                    :: alphaco
    real(rkp)                                                    :: alphax1, alphax2, alphax3
    real(rkp)                                                    :: ac6b, ac6oo, ac6co
    real(rkp)                                                    :: gdamp, damp6co, damp6cox
    real(rkp)                                                    :: gaus, gaus1, gaus2, gausv
    real(rkp)                                                    :: goff, goff3
    real(rkp)                                                    :: alpc0
    real(rkp)                                                    :: c6o0, c6c0
    real(rkp)                                                    :: alpo0
    real(rkp)                                                    :: et0
    real(rkp)                                                    :: vcupt0, vcupt1
    real(rkp)                                                    :: ret0
    real(rkp)                                                    :: cte0
    real(rkp)                                                    :: a1, a2
    real(rkp)                                                    :: diatshift
    real(rkp)                                                    :: Vd_Inf
    
    integer        ,dimension(3,128)                             :: iqn                  !!!!!!!!! <<<<< Because idc = 128
    integer        ,dimension(2,100)                             :: ixpd                 !!!!!!!!! <<<<< Because idx = 100
    real(rkp)      ,dimension(15)                                :: coefco
    real(rkp)      ,dimension(15)                                :: coefo2
    real(rkp)      ,dimension(100*2)                             :: par                  !!!!!!!!! <<<<< Because idx = 100
    real(rkp)      ,dimension(128)                               :: qff
    real(rkp)      ,dimension(7)                                 :: xm 
    
  contains
    procedure          ::  Initialize        =>    Initialize_CO2_NASA_PES
    procedure          ::  Output            =>    Output_CO2_NASA_PES
    procedure          ::  Compute           =>    Compute_CO2_NASA_PES_1d
    procedure          ::  Potential         =>    CO2_NASA_Potential_From_R
    procedure ,private ::  Compute_CO2_NASA_PES_1d_NoTraj
    procedure ,private ::  CO2_NASA_Potential_From_R_NoTraj
  End Type
   
  real(rkp)      ,dimension(3,3)                               :: cart
  real(rkp)      ,dimension(3)                                 :: dcc
  real(rkp)      ,dimension(7,3)                               :: fmat 
  real(rkp)      ,dimension(7,3)                               :: dfmat
  real(rkp)      ,dimension(4,4,3)                             :: dham
  real(rkp)      ,dimension(4)                                 :: eig
  real(rkp)      ,dimension(4)                                 :: fv1
  real(rkp)      ,dimension(4)                                 :: fv2
  real(rkp)      ,dimension(4,4)                               :: ham
  real(rkp)      ,dimension(4,4)                               :: ham2
  real(rkp)      ,dimension(200)                               :: scr 
  
  logical                                        ,parameter    :: i_Debug_Global = .False.
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NASA PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_CO2_NASA_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                      ,only:  Input_Type
  use Atom_Class                       ,only:  Atom_Type
  use CO_DiatomicPotential_Class       ,only:  CO_DiatomicPotential_Type
  use O2_NASA_DiatomicPotential_Class  ,only:  O2_NASA_DiatomicPotential_Type
  
  class(CO2_NASA_PES_Type)                  ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms 
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP
  character(*)                    ,parameter                ::    Name_PES = 'CO2_NASA'
  character(:)       ,allocatable                           ::    NASA_CO2_file
  integer                                                   ::    Unit
  integer                                                   ::    Status
  integer                                                   ::    k
  integer                                                   ::    kp
  integer                                                   ::    ii
  integer            ,dimension(6)                          ::    MatrixTemp = (/ 1, 2, 1, 3, 2, 3 /)
  type(CO_DiatomicPotential_Type)                           ::    CO_DiatPot
  integer         ,dimension(3,2)                           ::    iA

  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize CO2_NASA PES")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
    
  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair
  ii = 1

  ! allocate( This%mMiMn(3) )
  ! This%mMiMn(1:2) = - Atoms(1:2)%Mass / Atoms(3)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )
  
  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]
  
  do iP = 1,This%NPairs
    if ( ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+1))) ) .ne. "C") .and. ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+2)))) .ne. "C" ) ) then
      allocate( O2_NASA_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      This%iO2 = iP
      if (i_Debug_Loc) call Logger%Write( 'O2 pair =',iP )
    else
      allocate( CO_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      if (ii == 1) then
        This%iCO=iP
        if (i_Debug_Loc) call Logger%Write( 'First CO pair =',iP )
      else
        This%jCO=iP
        if (i_Debug_Loc) call Logger%Write( 'Second CO pair =',iP )
      end if
      ii=ii+1
    end if
  end do
  
  NASA_CO2_file = trim(adjustl(Input%DtbPath))  // '/Systems/CO2/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading NASA CO2 PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", NASA_CO2_file)
  open( File=NASA_CO2_file, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // NASA_CO2_file )
  
  read(Unit,"(A80)") This%line
  if (i_Debug_Loc) call Logger%Write( "Line: ", This%line )

  This%xm(1)=16.0_rkp
  This%xm(2)=12.0_rkp
  This%xm(3)=16.0_rkp
  
  read(Unit,*) This%c6o2, This%c6co, This%damp6co, This%ao2, This%aco, This%r0co, This%alphaco, This%npowco, This%nfitco
  if (i_Debug_Loc) call Logger%Write( "c6o2 = ", This%c6o2)
  if (i_Debug_Loc) call Logger%Write( "c6co = ", This%c6co )
  if (i_Debug_Loc) call Logger%Write( "damp6co = ", This%damp6co )
  if (i_Debug_Loc) call Logger%Write( "ao2 = ", This%ao2 )
  if (i_Debug_Loc) call Logger%Write( "aco = ", This%aco )
  if (i_Debug_Loc) call Logger%Write( "r0co = ", This%r0co )
  if (i_Debug_Loc) call Logger%Write( "npowco = ", This%npowco )
  if (i_Debug_Loc) call Logger%Write( "nfitco = ", This%nfitco )

  read(Unit,*) ( This%coefo2(k), k = 1,This%nfitco+1 )
  if (i_Debug_Loc) call Logger%Write( "coefo2 = ", This%coefo2 )
  
  read(Unit,*) ( This%coefco(k), k = 1,This%nfitco+1 )
  if (i_Debug_Loc) call Logger%Write( "coefco = ", This%coefco )
  
  read(Unit,*) This%alphax1, This%alphax2, This%alphax3, This%ac6oo, This%ac6co, This%ac6b, This%damp6cox
  if (i_Debug_Loc) call Logger%Write( "alphax1 = ", This%alphax1 )
  if (i_Debug_Loc) call Logger%Write( "alphax2 = ", This%alphax2 )
  if (i_Debug_Loc) call Logger%Write( "alphax3 = ", This%alphax3 )
  if (i_Debug_Loc) call Logger%Write( "ac6oo = ", This%ac6oo )
  if (i_Debug_Loc) call Logger%Write( "ac6co = ", This%ac6co )
  if (i_Debug_Loc) call Logger%Write( "ac6b = ", This%ac6b )
  if (i_Debug_Loc) call Logger%Write( "damp6cox = ", This%damp6cox )   
  
  read(Unit,*) This%nparh
  if (i_Debug_Loc) call Logger%Write( "Modulating function parameters.     nparh = ", This%nparh )
 
  do k = 1,This%nparh
    
    kp = k + This%nparh
    read(Unit,*) This%kdum, This%par(k), This%par(kp), This%ixpd(1,k), This%ixpd(2,k)
    if (i_Debug_Loc) call Logger%Write( "k = ", k )
    if (i_Debug_Loc) call Logger%Write( "kp = ", kp )
    if (i_Debug_Loc) call Logger%Write( "kdum = ", This%kdum )
    if (i_Debug_Loc) call Logger%Write( "par(", k, ") = ", This%par(k) )
    if (i_Debug_Loc) call Logger%Write( "par(", kp, ") = ", This%par(kp ))
    if (i_Debug_Loc) call Logger%Write( "ixpd(1,", k, ") = ", This%ixpd(1,k) )
    if (i_Debug_Loc) call Logger%Write( "ixpd(2,", k, ") = ", This%ixpd(2,k) )

  end do
   
  read(Unit,*) This%gdamp, This%gaus, This%gaus1, This%gaus2, This%gausv
  if (i_Debug_Loc) call Logger%Write( "gdamp = ", This%gdamp )
  if (i_Debug_Loc) call Logger%Write( "gaus = ", This%gaus )
  if (i_Debug_Loc) call Logger%Write( "gaus1 = ", This%gaus1 )
  if (i_Debug_Loc) call Logger%Write( "gaus2 = ", This%gaus2 )
  if (i_Debug_Loc) call Logger%Write( "gausv = ", This%gausv )
 
  read(Unit,*) This%goff, This%goff3, This%xm(4), This%xm(6)
  if (i_Debug_Loc) call Logger%Write( "goff = ", This%goff )
  if (i_Debug_Loc) call Logger%Write( "goff3 = ", This%goff3 )
  if (i_Debug_Loc) call Logger%Write( "xm(4) = ", This%xm(4) )
  if (i_Debug_Loc) call Logger%Write( "xm(6) = ", This%xm(6) )
 
  read(Unit,*) This%alpc0, This%c6c0, This%alpo0, This%c6o0, This%damp6co, This%et0, This%xm(5), This%vcupt0, This%xm(7), This%vcupt1
  if (i_Debug_Loc) call Logger%Write( "alpc0 = ", This%alpc0 )
  if (i_Debug_Loc) call Logger%Write( "c6c0 = ", This%c6c0 )
  if (i_Debug_Loc) call Logger%Write( "alpo0 = ", This%alpo0 )
  if (i_Debug_Loc) call Logger%Write( "c6o0 = ", This%c6o0 )
  if (i_Debug_Loc) call Logger%Write( "damp6co = ", This%damp6co )
  if (i_Debug_Loc) call Logger%Write( "et0 = ", This%et0 )
  if (i_Debug_Loc) call Logger%Write( "xm(5) = ", This%xm(5) )
  if (i_Debug_Loc) call Logger%Write( "vcupt0 = ", This%vcupt0 )
  if (i_Debug_Loc) call Logger%Write( "xm(7) = ", This%xm(7) )
  if (i_Debug_Loc) call Logger%Write( "vcupt1 = ", This%vcupt1 )
 
  read(Unit,*) This%ret0, This%cte0, This%a1, This%a2
  if (i_Debug_Loc) call Logger%Write( "Global Min Force Field Parameters: ")
  if (i_Debug_Loc) call Logger%Write( "ret0 = ", This%ret0 )
  if (i_Debug_Loc) call Logger%Write( "cte0 = ", This%cte0 )
  if (i_Debug_Loc) call Logger%Write( "a1 = ", This%a1 )
  if (i_Debug_Loc) call Logger%Write( "a2 = ", This%a2 )
 
  read(Unit,*) This%nfcnt
  if (i_Debug_Loc) call Logger%Write( "nfcnt = ", This%nfcnt)
 
  This%Vd_Inf = CO_DiatPot%DiatomicPotential( 1.0e4_rkp )
  if (i_Debug_Loc) call Logger%Write( "Vd at C+O Limit:     This%Vd_Inf", This%Vd_Inf )
  if (i_Debug_Loc) call Logger%Write( "coefco(1) = ", This%coefco(1) )
  This%diatshift = This%coefco(1) - This%Vd_inf
  if (i_Debug_Loc) call Logger%Write( "diatshift = ", This%diatshift )
   
  do k = 1,This%nfcnt
 
    read(Unit,*) This%kdum, This%qff(k), This%iqn(1,k), This%iqn(2,k), This%iqn(3,k)
    if (i_Debug_Loc) call Logger%Write( "k = ", k )
    if (i_Debug_Loc) call Logger%Write( "kdum = ", This%kdum )
    if (i_Debug_Loc) call Logger%Write( "qff(", k, ") = ", This%qff(k) )
    if (i_Debug_Loc) call Logger%Write( "iqn(1,", k, ") = ", This%iqn(1,k) )
    if (i_Debug_Loc) call Logger%Write( "iqn(2,", k, ") = ", This%iqn(2,k) )
    if (i_Debug_Loc) call Logger%Write( "iqn(3,", k, ") = ", This%iqn(3,k) )
    
  end do
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_CO2_NASA_PES( This, Unit )

  class(CO2_NASA_PES_Type)                ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('CO2 PES with cta bug fixed')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function CO2_NASA_Potential_From_R( This, R, Q ) result( V )

  class(CO2_NASA_PES_Type)                      ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].
  
  call CO2_NASA_Potential_From_R_NoTraj( This, R(This%iO2), R(This%iCO), R(This%jCO), V )

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_CO2_NASA_PES_1d( This, R, Q, V, dVdR, dVdQ )

  use CO_DiatomicPotential_Class        ,only:  CO_DiatomicPotential_Type
  use O2_NASA_DiatomicPotential_Class   ,only:  O2_NASA_DiatomicPotential_Type
    
  class(CO2_NASA_PES_Type)                      ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  type(CO_DiatomicPotential_Type)                            :: CO_DiatPot
  type(O2_NASA_DiatomicPotential_Type)                       :: O2_NASA_DiatPot
  real(rkp)                                                  :: VTemp

  call Compute_CO2_NASA_PES_1d_NoTraj( This, R(This%iO2), R(This%iCO), R(This%jCO), V, dVdR(:) )

  dVdQ = Zero 
  call This%TransToCart_3Atoms( R, Q, dVdR, dVdQ)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_CO2_NASA_PES_1d_NoTraj( This, R_O2, R_CO_1, R_CO_2, V, dV, i_Debug )
  
  use CO_DiatomicPotential_Class        ,only:  CO_DiatomicPotential_Type
  use O2_NASA_DiatomicPotential_Class   ,only:  O2_NASA_DiatomicPotential_Type
  
  class(CO2_NASA_PES_Type)                  ,intent(in)     :: This
  real(rkp)                                 ,intent(in)     :: R_O2
  real(rkp)                                 ,intent(in)     :: R_CO_1
  real(rkp)                                 ,intent(in)     :: R_CO_2
  real(rkp)                                 ,intent(out)    :: V
  real(rkp)       ,dimension(3)             ,intent(out)    :: dV
  logical                         ,optional ,intent(in)     :: i_Debug
  
  type(CO_DiatomicPotential_Type)                           :: CO_DiatPot
  
  real(rkp)                                                 :: roo 
   
  integer                                                   :: i, j
  integer                                                   :: k, k1, k2, k3
  integer                                                   :: kp
  integer                                                   :: km, k1m, k2m, k3m
  integer                                                   :: ierr
  integer                                                   :: ndim 
  
  real(rkp)                                                 :: calf, dcalf3, dcalf3j, dcalfj, dcalft
  real(rkp)                                                 :: cc
  real(rkp)                                                 :: cfcn, dcfcn
  real(rkp)                                                 :: costheta
  real(rkp)                                                 :: cthet, sthet, cott
  real(rkp)                                                 :: d131, d132, d133, d231, d232, d233
  real(rkp)                                                 :: damp, damp2
  real(rkp)                                                 :: dcostj, dcost3j, dcost3 
  real(rkp)                                                 :: dct1, dct2, dct3
  real(rkp)                                                 :: dm1, dm2, dm3
  real(rkp)                                                 :: doff, ddoff, ddoff3, ddoff3j, ddoffj
  real(rkp)                                                 :: drbigj, drbig3j, drbig3
  real(rkp)                                                 :: drj1, drj2, drj3
  real(rkp)                                                 :: drh1, drh2, drh3
  real(rkp)                                                 :: dvdr
  real(rkp)                                                 :: dvoolr
  real(rkp)                                                 :: ex1, ex2, ex3, ex4, ex5, ex6, ex7, exg, exr
  real(rkp)                                                 :: fact, dfact
  real(rkp)                                                 :: fc, fo
  real(rkp)                                                 :: part1, part2, part3
  real(rkp)                                                 :: raw, rawi
  real(rkp)                                                 :: r15, r25, r35   
  real(rkp)                                                 :: rbig, rbig2, rbigi
  real(rkp)                                                 :: rh
  real(rkp)                                                 :: rjac
  real(rkp)                                                 :: ro
  real(rkp)                                                 :: rho, rhoi
  real(rkp)                                                 :: realf, drealf3, drealf3j, drealfj 
  real(rkp)                                                 :: sum1, sum2, dsum11, dsum12, dsum21, dsum22
  real(rkp)                                                 :: term, term0, term1, term2, term3
  real(rkp)                                                 :: V_O2
  real(rkp)                                                 :: v0o, dv0o1, dv0o2, dv0o3
  real(rkp)                                                 :: v1, v2, v3, v4
  real(rkp)                                                 :: vco, dvco, vcolr, dvcolr, voco, voco0
  real(rkp)                                                 :: vdiatm
  real(rkp)                                                 :: vlr, dvlr, vlr1, vlr2, vlr3
  real(rkp)                                                 :: voo, dvoo, voolr, dvool
  real(rkp)                                                 :: vpair
  real(rkp)                                                 :: vsr, dvsr
  real(rkp)                                                 :: xmv, xmr, xms
  real(rkp)                                                 :: xx
  integer                                                   :: iO2
  integer                                                   :: iCO
  integer                                                   :: jCO
  
  iO2=This%iO2
  iCO=This%iCO
  jCO=This%jCO
    
  roo    = R_O2
  fv1(1) = R_CO_1
  fv1(2) = R_CO_2
  

  fc  = This%xm(1) / (This%xm(1)+This%xm(2))
  fo  = This%xm(2) / (This%xm(1)+This%xm(2))
  xmv = This%xm(1) * This%xm(2) / (This%xm(1)+This%xm(2))
  xmr = This%xm(3) * (This%xm(1)+This%xm(2)) / (This%xm(1)+This%xm(2)+This%xm(3))
  xms = dsqrt(xmr/xmv)
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This is Simply:    O2_NASA_DiatPot%Compute_Vd_dVd_O2( R_O2, V_O2=ham(3,3), dham(3,3,3) )
  vsr   = 64.0_rkp * dexp(This%alphaco*roo) / roo
  dvsr  = vsr * (This%alphaco-(One/roo))
  vlr   = -This%c6o2 / (roo**6+This%damp6co)
  dvlr  = -Six * (roo**5) * vlr / (roo**6+This%damp6co)
  fact  = (roo**This%npowco) * dexp(-This%ao2*roo)
  dfact = fact * ((dfloat(This%npowco)/roo)-This%ao2)
  cfcn  = This%coefo2(This%nfitco+1)
  dcfcn = Zero
  do k=This%nfitco,2,-1
        
    dcfcn = cfcn      + dcfcn * (roo-This%r0co)
    cfcn  = This%coefo2(k) + cfcn  * (roo-This%r0co)
    
  end do    
  dcfcn   = dfact * cfcn + fact * dcfcn               
  cfcn    = This%coefo2(1)    + fact * cfcn
  
  V_O2        = cfcn + vsr + vlr !+ 187.888473490235_rkp !????
  !write(222,*) "V_O2", V_O2 
  ham(3,3)    = V_O2
  dham(3,3,3) = dcfcn + dvlr + dvsr
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  voo = 64.0_rkp * dexp(This%alphax1*roo) / roo
  dvoo = voo * (This%alphax1-(One/roo))
  voolr = -This%ac6oo / (roo**6+This%damp6cox)
  dvoolr = -Six * (roo**5) * voolr / (roo**6+This%damp6cox)
 


  cc     = - (roo**2 - fv1(1)**2 - fv1(2)**2) *Half / (fv1(1)*fv1(2))
  dcc(1) = (One/fv1(2)) - (cc/fv1(1))
  dcc(2) = (One/fv1(1)) - (cc/fv1(2))
  dcc(3) = -roo / (fv1(1)*fv1(2))
      
      
      
  do j = 1,2
  
    call CO_DiatPot%Compute_Vd_dVd( fv1(j), vdiatm, dvdr )   
   
    vdiatm      = vdiatm + This%diatshift
    dham(j,j,j) = dvdr
    
    vco=48.0_rkp*dexp(This%alphax3*fv1(3-j))/fv1(3-j)
    dvco=vco*(This%alphax3-(One/fv1(3-j)))
    !write(222,*) "vco 1", vco
    !write(222,*) "dvco 1", dvco
    
    vcolr=-This%ac6b/(fv1(3-j)**6+This%damp6cox)
    dvcolr=-Six*(fv1(3-j)**5)*vcolr/(fv1(3-j)**6+This%damp6cox)
    
    dham(3,3,3-j)=dvco+dvcolr
    ham(3,3)=ham(3,3)+vco+vcolr
    
    vco=48.0_rkp*dexp(This%alphax2*fv1(3-j))/fv1(3-j)
    dvco=vco*(This%alphax2-(One/fv1(3-j)))
    !write(222,*) "vco 2", vco
    !write(222,*) "dvco 2", dvco
    
    vcolr=-This%ac6co/(fv1(3-j)**6+This%damp6cox)
    dvcolr=-Six*(fv1(3-j)**5)*vcolr/(fv1(3-j)**6+This%damp6cox)
    
    dham(j,j,3)=dvoolr
    dham(j,j,3-j)=dvcolr
    
    vpair=voo+vco
    
    rh=fv1(j)*Half
    ro=fv1(3-j)
    rbig2=rh*rh+ro*ro-Two*rh*ro*cc
    rbig=dsqrt(rbig2)
    rbigi=One/rbig
    costheta=-(ro*ro-rh*rh-rbig2)*Half/(rbig*rh)
    drbigj=(Half*rh-ro*(Half*cc+rh*dcc(j)))*rbigi
    drbig3j=(ro-rh*(cc+ro*dcc(3-j)))*rbigi
    drbig3=-rh*ro*dcc(3)*rbigi
    dcostj=(Half/rbig)+(drbigj/rh)-costheta*(drbigj/rbig+Half/rh)
    dcost3j=(((-ro/rbig)+drbig3j)/rh)-costheta*drbig3j/rbig
    dcost3=drbig3/rh-costheta*drbig3/rbig
   
    fmat(1,1)  = One
    fmat(2,1)  = fv1(j)
    dfmat(1,1) = Zero
    dfmat(2,1) = One
    fmat(1,2)  = One
    fmat(2,2)  = costheta
    dfmat(1,2) = Zero
    dfmat(2,2) = One
       
    do k = 3,7
    
      km = k - 1
      fmat(k,1)  = fmat(2,1)  * fmat(km,1)
      dfmat(k,1) = dfmat(2,1) * fmat(km,1) + fmat(2,1) * dfmat(km,1)
      fmat(k,2)  = fmat(2,2)  * fmat(km,2)
      dfmat(k,2) = dfmat(2,2) * fmat(km,2) + fmat(2,2) * dfmat(km,2)
      
    end do
       
    sum1   = Zero
    sum2   = Zero
    dsum11 = Zero
    dsum12 = Zero
    dsum21 = Zero
    dsum22 = Zero
       
    do k = 1,This%nparh
    
      kp=k+This%nparh
      !write(222,*) 'kp', kp
      term=fmat(This%ixpd(1,k),1)*fmat(This%ixpd(2,k),2)
      !write(222,*) 'par(k)', This%par(k)
      !write(222,*) 'term', term
      sum1=sum1+This%par(k)*term
      !write(222,*) 'dfmat', dfmat(This%ixpd(1,k),1)
      !write(222,*) 'fmat',  fmat(This%ixpd(2,k),2)
      dsum11 = dsum11 + This%par(k) * dfmat(This%ixpd(1,k),1) * fmat(This%ixpd(2,k),2)
      dsum12 = dsum12 + This%par(k) * fmat(This%ixpd(1,k),1)  * dfmat(This%ixpd(2,k),2)
      sum2=sum2+This%par(kp)*term
      dsum21 = dsum21 + This%par(kp) * dfmat(This%ixpd(1,k),1) * fmat(This%ixpd(2,k),2)
      dsum22 = dsum22 + This%par(kp) * fmat(This%ixpd(1,k),1)  * dfmat(This%ixpd(2,k),2)
      !write(222,*) 'sum1', sum1
      !write(222,*) 'sum2', sum2
      !write(222,*) 'dsum11', dsum11
      !write(222,*) 'dsum12', dsum12
      !write(222,*) 'dsum21', dsum21
      !write(222,*) 'dsum22', dsum22
        
    end do
    
    ex1=dexp(-This%gaus*((rbig-This%gaus1)**2))
    ex2=dexp(-This%gaus*((rbig-This%gaus2)**2))
    exr=dexp(-This%gdamp*fv1(j))
    calf=exr*(ex1*sum1+ex2*sum2)
    dcalft=Two*This%gaus*exr*((rbig-This%gaus1)*ex1*sum1+(rbig-This%gaus2)*ex2*sum2)
    dcalfj=-drbigj*dcalft-This%gdamp*calf
    dcalf3j=-drbig3j*dcalft
    dcalf3=-drbig3*dcalft+exr*dcost3*(ex1*dsum12+ex2*dsum22)
    dcalfj=dcalfj+exr*(ex1*(dsum11+dsum12*dcostj)+ex2*(dsum21+dsum22*dcostj))
    dcalf3j=dcalf3j+exr*(ex1*dsum12*dcost3j+ex2*dsum22*dcost3j)
    exg=dexp(-This%gausv*calf)
    doff=One/(One+exg)
    ddoff=doff*doff*This%gausv*exg
    ddoffj=ddoff*dcalfj
    ddoff3j=ddoff*dcalf3j
    ddoff3=ddoff*dcalf3
    realf=One+calf*doff
    drealfj=dcalfj*doff+calf*ddoffj
    drealf3j=dcalf3j*doff+calf*ddoff3j
    drealf3=dcalf3*doff+calf*ddoff3
    ham(j,j)=vdiatm+vpair*realf+voolr+vcolr
    dham(j,j,j)=dham(j,j,j)+vpair*drealfj
    dham(j,j,3-j)=dham(j,j,3-j)+dvco*realf+vpair*drealf3j
    dham(j,j,3)=dham(j,j,3)+dvoo*realf+vpair*drealf3
       
  end do
      
  cart(1,1)=Zero
  cart(2,1)=fv1(2)*fo
  cart(1,2)=Zero
  cart(2,2)=-fv1(2)*fc
  
  cthet = -(roo**2 - fv1(1)**2 - fv1(2)**2) * Half / (fv1(1)*fv1(2))
  
  dct1  = (One/fv1(2)) - (cthet/fv1(1))
  dct2  = (One/fv1(1)) - (cthet/fv1(2))
  dct3  = -roo / (fv1(1)*fv1(2))
  
  sthet=dsqrt(dabs(One-cthet**2))
  cott=cthet/sthet
  
  cart(1,3)=fv1(1)*sthet
  cart(2,3)=fv1(1)*cthet+cart(2,2)
  
  d131=dsqrt(dabs(One-cthet**2))-fv1(1)*dct1*cott
  d132=-fv1(1)*dct2*cott
  d133=-fv1(1)*dct3*cott
  d231=cthet+fv1(1)*dct1
  d232=fv1(1)*dct2-fc
  d233=fv1(1)*dct3
  
  raw=dsqrt(cart(1,3)**2+cart(2,3)**2)
  rawi=xms/raw
  rjac=xms*raw
  
  drj1=(cart(1,3)*d131+cart(2,3)*d231)*rawi
  drj2=(cart(1,3)*d132+cart(2,3)*d232)*rawi
  drj3=(cart(1,3)*d133+cart(2,3)*d233)*rawi
  
  rho=dsqrt(fv1(2)**2+rjac**2)
  rhoi=One/rho
  
  drh1=drj1*rjac*rhoi
  drh2=(fv1(2)+drj2*rjac)*rhoi
  drh3=drj3*rjac*rhoi

  ex4=dexp(-This%xm(4)*rho) + 1.e-7_rkp
  ex5=dexp(-This%xm(5)*rho) + 1.e-7_rkp
  ex6=dexp(-This%xm(6)*rho) + 1.e-7_rkp
  ex7=dexp(-This%xm(7)*rho) + 1.e-7_rkp
  
  
  
  ham(1,2)=This%goff*ex4
  ham(2,1)=ham(1,2)
  ham(1,3)=This%goff3*ex6
  ham(2,3)=ham(1,3)
  ham(3,1)=ham(1,3)
  ham(3,2)=ham(1,3)
  
  xx=-This%xm(4)*(ham(1,2)-This%goff * 1.e-7_rkp)
  
  dham(1,2,1)=xx*drh1
  dham(1,2,2)=xx*drh2
  dham(1,2,3)=xx*drh3
  dham(2,1,1)=dham(1,2,1)
  dham(2,1,2)=dham(1,2,2)
  dham(2,1,3)=dham(1,2,3)
  
  xx=-This%xm(6)*(ham(1,3)-This%goff3 * 1.e-7_rkp)
  
  dham(1,3,1)=xx*drh1
  dham(1,3,2)=xx*drh2
  dham(1,3,3)=xx*drh3
  dham(2,3,1)=dham(1,3,1)
  dham(2,3,2)=dham(1,3,2)
  dham(2,3,3)=dham(1,3,3)
  dham(3,2,1)=dham(1,3,1)
  dham(3,2,2)=dham(1,3,2)
  dham(3,2,3)=dham(1,3,3)
  dham(3,1,1)=dham(1,3,1)
  dham(3,1,2)=dham(1,3,2)
  dham(3,1,3)=dham(1,3,3)
  ndim=3
  
  
      
  if (This%et0 .ne. Zero) then
      
    ndim=4
    ex1 = 48.0_rkp * dexp(-This%alpc0*fv1(1)) / fv1(1)
    r15 = fv1(1)**5
    vlr1=-One/(r15*fv1(1)+This%damp6co)
    part1=ex1+This%c6c0*vlr1
    ex2=48.0_rkp*dexp(-This%alpc0*fv1(2))/fv1(2)
    r25=fv1(2)**5
    vlr2=-One/(fv1(2)*r25+This%damp6co)
    part2=ex2+This%c6c0*vlr2
    ex3=64.0_rkp*dexp(-This%alpo0*roo)/roo
    r35=roo**5
    vlr3=-One/(roo*r35+This%damp6co)
    part3=ex3+This%c6o0*vlr3
    v0o=This%et0+part1+part2+part3
    xx=Six*This%c6c0
    dv0o1=-ex1*(This%alpc0+One/fv1(1))+xx*r15*vlr1*vlr1
    dv0o2=-ex2*(This%alpc0+One/fv1(2))+xx*r25*vlr2*vlr2
    dv0o3=-ex3*(This%alpo0+One/roo)+Six*This%c6o0*r35*vlr3*vlr3
       
    if (This%nfcnt .gt. 0) then
       
      scr(1)=One
      scr(1+5)=One
      scr(1+10)=One
      scr(2)=fv1(1)-This%ret0
      scr(2+5)=fv1(2)-This%ret0
      scr(2+10)=cc-This%cte0
      damp=dexp(-This%a1*(scr(2)**2+scr(7)**2)-This%a2*(scr(12)**2))
      damp2=damp*Two
      xx=-This%a2*scr(12)*damp2
      dm1=-This%a1*scr(2)*damp2+xx*dcc(1)
      dm2=-This%a1*scr(7)*damp2+xx*dcc(2)
      dm3=xx*dcc(3)
        
      do k = 3,5
      
        k1  = k
        k1m = k1 - 1 
        scr(k1) = scr(k1m) * scr(2)
        k2  = k + 5
        k2m = k2 - 1
        scr(k2) = scr(k2m) * scr(2+5)
        k3  = k + 10
        k3m = k3 - 1 
        scr(k3) = scr(k3m) * scr(2+10)
         
      end do
        
      voco=Zero
        
      dV = Zero
        
      do k = 1,This%nfcnt
        
        term0 = (scr(This%iqn(1,k)) * scr(This%iqn(2,k)+5) + scr(This%iqn(1,k)+5) * scr(This%iqn(2,k)))
        term  = scr(This%iqn(3,k)+10) * term0
        term3 = dfloat(This%iqn(3,k)-1) * scr(int(max(This%iqn(3,k)-One,One)+10.0_rkp)) * term0
        term1 = scr(This%iqn(3,k)+10) * (dfloat(This%iqn(1,k)-1)  * scr(int(max(This%iqn(1,k)-One,One))) * scr(This%iqn(2,k)+5) + &
                   dfloat(This%iqn(2,k)-1) * scr(This%iqn(1,k)+5) * scr(int(max(This%iqn(2,k)-1,1)  )))
        term2 = scr(This%iqn(3,k)+10) * (dfloat(This%iqn(1,k)-1)  * scr(int(max(This%iqn(1,k)-1,1))+5) * scr(This%iqn(2,k)  ) + &
                   dfloat(This%iqn(2,k)-1) * scr(This%iqn(1,k)  ) * scr(int(max(This%iqn(2,k)-1,1))+5))
   
        dV(1) = dV(1) +This%qff(k) * (term3*dcc(1)+term1)
        dV(2) = dV(2) +This%qff(k) * (term3*dcc(2)+term2)
        dV(3) = dV(3) +This%qff(k) * term3 * dcc(3)
       
        voco = voco + This%qff(k) * term
         
      end do
        
      dham(4,4,1) = dV(1) *damp + voco * dm1 + dv0o1
      dham(4,4,2) = dV(2) *damp + voco * dm2 + dv0o2
      dham(4,4,3) = dV(3) *damp + voco * dm3 + dv0o3
      voco0 = voco
      voco  = v0o + damp * voco
        
    else
        
      voco  = v0o
      voco0 = Zero
      dham(4,4,1) = dv0o1
      dham(4,4,2) = dv0o2
      dham(4,4,3) = dv0o3
         
    end if
        
    ham(4,4) = voco
    ham(1,4) = This%vcupt0 * ex5
    ham(2,4) = This%vcupt0 * ex5
    ham(4,2) = This%vcupt0 * ex5
    ham(4,1) = This%vcupt0 * ex5
    xx = -This%xm(5) * (ham(1,4)-This%vcupt0 * 1.e-7_rkp)
    
    dham(1,4,1)=xx * drh1
    dham(1,4,2)=xx * drh2
    dham(1,4,3)=xx * drh3
    
    dham(2,4,1) = dham(1,4,1)
    dham(2,4,2) = dham(1,4,2)
    dham(2,4,3) = dham(1,4,3)
    dham(4,2,1) = dham(1,4,1)
    dham(4,2,2) = dham(1,4,2)
    dham(4,2,3) = dham(1,4,3)
    dham(4,1,1) = dham(1,4,1)
    dham(4,1,2) = dham(1,4,2)
    dham(4,1,3) = dham(1,4,3)
    ham(3,4)    = This%vcupt1 * ex7
    
    xx = -This%xm(7) * (ham(3,4)-This%vcupt1 * 1.e-7_rkp)
    
    dham(3,4,1) = xx * drh1
    dham(3,4,2) = xx * drh2
    dham(3,4,3) = xx * drh3
    dham(4,3,1) = dham(3,4,1)
    dham(4,3,2) = dham(3,4,2)
    dham(4,3,3) = dham(3,4,3)
    ham(4,3)    = ham(3,4)
        
  end if
  
   
      
  do i = 1,ndim
    do j = 1,ndim
      ham2(j,i) = ham(j,i)
     end do
  end do
      
  call rs(4,ndim,ham,eig,1,ham,fv1,fv2,ierr)
      
  V = eig(1)
      
  v1 = ham(1,1)*dham(1,1,1)+ham(2,1)*dham(2,1,1)+ham(3,1)*dham(3,1,1)+ham(4,1)*dham(4,1,1)
  v2 = ham(1,1)*dham(1,2,1)+ham(2,1)*dham(2,2,1)+ham(3,1)*dham(3,2,1)+ham(4,1)*dham(4,2,1)
  v3 = ham(1,1)*dham(1,3,1)+ham(2,1)*dham(2,3,1)+ham(3,1)*dham(3,3,1)+ham(4,1)*dham(4,3,1)
  v4 = ham(1,1)*dham(1,4,1)+ham(2,1)*dham(2,4,1)+ham(3,1)*dham(3,4,1)+ham(4,1)*dham(4,4,1)
  dV(iCO) = ham(1,1)*v1+ham(2,1)*v2+ham(3,1)*v3+ham(4,1)*v4
  
  v1 = ham(1,1)*dham(1,1,2)+ham(2,1)*dham(2,1,2)+ham(3,1)*dham(3,1,2)+ham(4,1)*dham(4,1,2)
  v2 = ham(1,1)*dham(1,2,2)+ham(2,1)*dham(2,2,2)+ham(3,1)*dham(3,2,2)+ham(4,1)*dham(4,2,2)
  v3 = ham(1,1)*dham(1,3,2)+ham(2,1)*dham(2,3,2)+ham(3,1)*dham(3,3,2)+ham(4,1)*dham(4,3,2)
  v4 = ham(1,1)*dham(1,4,2)+ham(2,1)*dham(2,4,2)+ham(3,1)*dham(3,4,2)+ham(4,1)*dham(4,4,2)
  dV(jCO) = ham(1,1)*v1+ham(2,1)*v2+ham(3,1)*v3+ham(4,1)*v4
  
  v1 = ham(1,1)*dham(1,1,3)+ham(2,1)*dham(2,1,3)+ham(3,1)*dham(3,1,3)+ham(4,1)*dham(4,1,3)
  v2 = ham(1,1)*dham(1,2,3)+ham(2,1)*dham(2,2,3)+ham(3,1)*dham(3,2,3)+ham(4,1)*dham(4,2,3)
  v3 = ham(1,1)*dham(1,3,3)+ham(2,1)*dham(2,3,3)+ham(3,1)*dham(3,3,3)+ham(4,1)*dham(4,3,3)
  v4 = ham(1,1)*dham(1,4,3)+ham(2,1)*dham(2,4,3)+ham(3,1)*dham(3,4,3)+ham(4,1)*dham(4,4,3)
  dV(iO2) = ham(1,1)*v1+ham(2,1)*v2+ham(3,1)*v3+ham(4,1)*v4

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine CO2_NASA_Potential_From_R_NoTraj( This, R_O2, R_CO_1, R_CO_2, V, i_Debug )
  
  use CO_DiatomicPotential_Class   ,only:  CO_DiatomicPotential_Type
  use O2_NASA_DiatomicPotential_Class   ,only:  O2_NASA_DiatomicPotential_Type
  
  class(CO2_NASA_PES_Type)                  ,intent(in)     :: This
  real(rkp)                                 ,intent(in)     :: R_O2
  real(rkp)                                 ,intent(in)     :: R_CO_1
  real(rkp)                                 ,intent(in)     :: R_CO_2
  real(rkp)                                 ,intent(out)    :: V
  logical                         ,optional ,intent(in)     :: i_Debug
  
  type(CO_DiatomicPotential_Type)                           :: CO_DiatPot
  
  real(rkp)                                                 :: roo 
   
  integer                                                   :: i, j
  integer                                                   :: k, k1, k2, k3
  integer                                                   :: kp
  integer                                                   :: km, k1m, k2m, k3m
  integer                                                   :: ierr
  integer                                                   :: ndim 
  
  real(rkp)                                                 :: calf, dcalf3, dcalf3j, dcalfj, dcalft
  real(rkp)                                                 :: cc
  real(rkp)                                                 :: cfcn, dcfcn
  real(rkp)                                                 :: costheta
  real(rkp)                                                 :: cthet, sthet, cott
  real(rkp)                                                 :: d131, d132, d133, d231, d232, d233
  real(rkp)                                                 :: damp, damp2
  real(rkp)                                                 :: dcostj, dcost3j, dcost3 
  real(rkp)                                                 :: dct1, dct2, dct3
  real(rkp)                                                 :: dm1, dm2, dm3
  real(rkp)                                                 :: doff, ddoff, ddoff3, ddoff3j, ddoffj
  real(rkp)                                                 :: drbigj, drbig3j, drbig3
  real(rkp)                                                 :: drj1, drj2, drj3
  real(rkp)                                                 :: drh1, drh2, drh3
  real(rkp)                                                 :: dvoolr
  real(rkp)                                                 :: ex1, ex2, ex3, ex4, ex5, ex6, ex7, exg, exr
  real(rkp)                                                 :: fact, dfact
  real(rkp)                                                 :: fc, fo
  real(rkp)                                                 :: part1, part2, part3
  real(rkp)                                                 :: raw, rawi
  real(rkp)                                                 :: r15, r25, r35   
  real(rkp)                                                 :: rbig, rbig2, rbigi
  real(rkp)                                                 :: rh
  real(rkp)                                                 :: rjac
  real(rkp)                                                 :: ro
  real(rkp)                                                 :: rho, rhoi
  real(rkp)                                                 :: realf, drealf3, drealf3j, drealfj 
  real(rkp)                                                 :: sum1, sum2, dsum11, dsum12, dsum21, dsum22
  real(rkp)                                                 :: term, term0, term1, term2, term3
  real(rkp)                                                 :: V_O2
  real(rkp)                                                 :: v0o, dv0o1, dv0o2, dv0o3
  real(rkp)                                                 :: v1, v2, v3, v4
  real(rkp)                                                 :: vco, dvco, vcolr, dvcolr, voco, voco0
  real(rkp)                                                 :: vdiatm
  real(rkp)                                                 :: vlr, dvlr, vlr1, vlr2, vlr3
  real(rkp)                                                 :: voo, dvoo, voolr, dvool
  real(rkp)                                                 :: vpair
  real(rkp)                                                 :: vsr, dvsr
  real(rkp)                                                 :: xmv, xmr, xms
  real(rkp)                                                 :: xx
  
  roo    = R_O2
  fv1(1) = R_CO_1
  fv1(2) = R_CO_2
  

  fc  = This%xm(1) / (This%xm(1)+This%xm(2))
  fo  = This%xm(2) / (This%xm(1)+This%xm(2))
  xmv = This%xm(1) * This%xm(2) / (This%xm(1)+This%xm(2))
  xmr = This%xm(3) * (This%xm(1)+This%xm(2)) / (This%xm(1)+This%xm(2)+This%xm(3))
  xms = dsqrt(xmr/xmv)
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This is Simply:    O2_NASA_DiatPot%Compute_Vd_dVd_O2( R_O2, V_O2=ham(3,3), dham(3,3,3) )
  vsr   = 64.0_rkp * dexp(This%alphaco*roo) / roo
  dvsr  = vsr * (This%alphaco-(One/roo))
  vlr   = -This%c6o2 / (roo**6+This%damp6co)
  dvlr  = -Six * (roo**5) * vlr / (roo**6+This%damp6co)
  fact  = (roo**This%npowco) * dexp(-This%ao2*roo)
  dfact = fact * ((dfloat(This%npowco)/roo)-This%ao2)
  cfcn  = This%coefo2(This%nfitco+1)
  dcfcn = Zero
  do k=This%nfitco,2,-1
        
    dcfcn = cfcn      + dcfcn * (roo-This%r0co)
    cfcn  = This%coefo2(k) + cfcn  * (roo-This%r0co)
    
  end do    
  dcfcn   = dfact * cfcn + fact * dcfcn               
  cfcn    = This%coefo2(1)    + fact * cfcn
  
  V_O2        = cfcn + vsr + vlr !+ 187.888473490235_rkp !????
  ham(3,3)    = V_O2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  voo = 64.0_rkp * dexp(This%alphax1*roo) / roo
  dvoo = voo * (This%alphax1-(One/roo))
  voolr = -This%ac6oo / (roo**6+This%damp6cox)
  dvoolr = -Six * (roo**5) * voolr / (roo**6+This%damp6cox)
 


  cc     = - (roo**2 - fv1(1)**2 - fv1(2)**2) *Half / (fv1(1)*fv1(2))
  dcc(1) = (One/fv1(2)) - (cc/fv1(1))
  dcc(2) = (One/fv1(1)) - (cc/fv1(2))
  dcc(3) = -roo / (fv1(1)*fv1(2))
      
      
      
  do j = 1,2
  
    vdiatm = CO_DiatPot%DiatomicPotential( fv1(j) ) 
   
    vdiatm      = vdiatm + This%diatshift
    
    vco=48.0_rkp*dexp(This%alphax3*fv1(3-j))/fv1(3-j)
    dvco=vco*(This%alphax3-(One/fv1(3-j)))
    
    vcolr=-This%ac6b/(fv1(3-j)**6+This%damp6cox)
    dvcolr=-Six*(fv1(3-j)**5)*vcolr/(fv1(3-j)**6+This%damp6cox)
    
    ham(3,3)=ham(3,3)+vco+vcolr
    
    vco=48.0_rkp*dexp(This%alphax2*fv1(3-j))/fv1(3-j)
    dvco=vco*(This%alphax2-(One/fv1(3-j)))
    
    vcolr=-This%ac6co/(fv1(3-j)**6+This%damp6cox)
    dvcolr=-Six*(fv1(3-j)**5)*vcolr/(fv1(3-j)**6+This%damp6cox)
    
    vpair=voo+vco
    
    rh=fv1(j)*Half
    ro=fv1(3-j)
    rbig2=rh*rh+ro*ro-Two*rh*ro*cc
    rbig=dsqrt(rbig2)
    rbigi=One/rbig
    costheta=-(ro*ro-rh*rh-rbig2)*Half/(rbig*rh)
    drbigj=(Half*rh-ro*(Half*cc+rh*dcc(j)))*rbigi
    drbig3j=(ro-rh*(cc+ro*dcc(3-j)))*rbigi
    drbig3=-rh*ro*dcc(3)*rbigi
    dcostj=(Half/rbig)+(drbigj/rh)-costheta*(drbigj/rbig+Half/rh)
    dcost3j=(((-ro/rbig)+drbig3j)/rh)-costheta*drbig3j/rbig
    dcost3=drbig3/rh-costheta*drbig3/rbig
   
    fmat(1,1)  = One
    fmat(2,1)  = fv1(j)
    dfmat(1,1) = Zero
    dfmat(2,1) = One
    fmat(1,2)  = One
    fmat(2,2)  = costheta
    dfmat(1,2) = Zero
    dfmat(2,2) = One
       
    do k = 3,7
    
      km = k - 1
      fmat(k,1)  = fmat(2,1)  * fmat(km,1)
      dfmat(k,1) = dfmat(2,1) * fmat(km,1) + fmat(2,1) * dfmat(km,1)
      fmat(k,2)  = fmat(2,2)  * fmat(km,2)
      dfmat(k,2) = dfmat(2,2) * fmat(km,2) + fmat(2,2) * dfmat(km,2)
      
    end do
       
    sum1   = Zero
    sum2   = Zero
    dsum11 = Zero
    dsum12 = Zero
    dsum21 = Zero
    dsum22 = Zero
       
    do k = 1,This%nparh
    
      kp=k+This%nparh
      term=fmat(This%ixpd(1,k),1)*fmat(This%ixpd(2,k),2)
      sum1=sum1+This%par(k)*term
      dsum11=dsum11+This%par(k)*dfmat(This%ixpd(1,k),1)*fmat(This%ixpd(2,k),2)
      dsum12=dsum12+This%par(k)*fmat(This%ixpd(1,k),1)*dfmat(This%ixpd(2,k),2)
      sum2=sum2+This%par(kp)*term
      dsum21=dsum21+This%par(kp)*dfmat(This%ixpd(1,k),1)*fmat(This%ixpd(2,k),2)
      dsum22=dsum22+This%par(kp)*fmat(This%ixpd(1,k),1)*dfmat(This%ixpd(2,k),2)
        
    end do
       
    ex1=dexp(-This%gaus*((rbig-This%gaus1)**2))
    ex2=dexp(-This%gaus*((rbig-This%gaus2)**2))
    exr=dexp(-This%gdamp*fv1(j))
    calf=exr*(ex1*sum1+ex2*sum2)
    dcalft=Two*This%gaus*exr*((rbig-This%gaus1)*ex1*sum1+(rbig-This%gaus2)*ex2*sum2)
    dcalfj=-drbigj*dcalft-This%gdamp*calf
    dcalf3j=-drbig3j*dcalft
    dcalf3=-drbig3*dcalft+exr*dcost3*(ex1*dsum12+ex2*dsum22)
    dcalfj=dcalfj+exr*(ex1*(dsum11+dsum12*dcostj)+ex2*(dsum21+dsum22*dcostj))
    dcalf3j=dcalf3j+exr*(ex1*dsum12*dcost3j+ex2*dsum22*dcost3j)
    exg=dexp(-This%gausv*calf)
    doff=One/(One+exg)
    ddoff=doff*doff*This%gausv*exg
    ddoffj=ddoff*dcalfj
    ddoff3j=ddoff*dcalf3j
    ddoff3=ddoff*dcalf3
    realf=One+calf*doff
    drealfj=dcalfj*doff+calf*ddoffj
    drealf3j=dcalf3j*doff+calf*ddoff3j
    drealf3=dcalf3*doff+calf*ddoff3
    ham(j,j)=vdiatm+vpair*realf+voolr+vcolr
       
  end do
      
  cart(1,1)=Zero
  cart(2,1)=fv1(2)*fo
  cart(1,2)=Zero
  cart(2,2)=-fv1(2)*fc
  
  cthet = -(roo**2 - fv1(1)**2 - fv1(2)**2) * Half / (fv1(1)*fv1(2))
  
  dct1  = (One/fv1(2)) - (cthet/fv1(1))
  dct2  = (One/fv1(1)) - (cthet/fv1(2))
  dct3  = -roo / (fv1(1)*fv1(2))
  
  sthet=dsqrt(dabs(One-cthet**2))
  cott=cthet/sthet
  
  cart(1,3)=fv1(1)*sthet
  cart(2,3)=fv1(1)*cthet+cart(2,2)
  
  d131=dsqrt(dabs(One-cthet**2))-fv1(1)*dct1*cott
  d132=-fv1(1)*dct2*cott
  d133=-fv1(1)*dct3*cott
  d231=cthet+fv1(1)*dct1
  d232=fv1(1)*dct2-fc
  d233=fv1(1)*dct3
  
  raw=dsqrt(cart(1,3)**2+cart(2,3)**2)
  rawi=xms/raw
  rjac=xms*raw
  
  drj1=(cart(1,3)*d131+cart(2,3)*d231)*rawi
  drj2=(cart(1,3)*d132+cart(2,3)*d232)*rawi
  drj3=(cart(1,3)*d133+cart(2,3)*d233)*rawi
  
  rho=dsqrt(fv1(2)**2+rjac**2)
  rhoi=One/rho
  
  drh1=drj1*rjac*rhoi
  drh2=(fv1(2)+drj2*rjac)*rhoi
  drh3=drj3*rjac*rhoi

  ex4=dexp(-This%xm(4)*rho) + 1.e-7_rkp
  ex5=dexp(-This%xm(5)*rho) + 1.e-7_rkp
  ex6=dexp(-This%xm(6)*rho) + 1.e-7_rkp
  ex7=dexp(-This%xm(7)*rho) + 1.e-7_rkp
  
  
  
  ham(1,2)=This%goff*ex4
  ham(2,1)=ham(1,2)
  ham(1,3)=This%goff3*ex6
  ham(2,3)=ham(1,3)
  ham(3,1)=ham(1,3)
  ham(3,2)=ham(1,3)
  
  xx=-This%xm(4)*(ham(1,2)-This%goff * 1.e-7_rkp)
  
  xx=-This%xm(6)*(ham(1,3)-This%goff3 * 1.e-7_rkp)
  
  ndim=3
  
  
      
  if (This%et0 .ne. Zero) then
      
    ndim=4
    ex1 = 48.0_rkp * dexp(-This%alpc0*fv1(1)) / fv1(1)
    r15 = fv1(1)**5
    vlr1=-One/(r15*fv1(1)+This%damp6co)
    part1=ex1+This%c6c0*vlr1
    ex2=48.0_rkp*dexp(-This%alpc0*fv1(2))/fv1(2)
    r25=fv1(2)**5
    vlr2=-One/(fv1(2)*r25+This%damp6co)
    part2=ex2+This%c6c0*vlr2
    ex3=64.0_rkp*dexp(-This%alpo0*roo)/roo
    r35=roo**5
    vlr3=-One/(roo*r35+This%damp6co)
    part3=ex3+This%c6o0*vlr3
    v0o=This%et0+part1+part2+part3
    xx=Six*This%c6c0
    dv0o1=-ex1*(This%alpc0+One/fv1(1))+xx*r15*vlr1*vlr1
    dv0o2=-ex2*(This%alpc0+One/fv1(2))+xx*r25*vlr2*vlr2
    dv0o3=-ex3*(This%alpo0+One/roo)+Six*This%c6o0*r35*vlr3*vlr3
       
    if (This%nfcnt .gt. 0) then
       
      scr(1)=One
      scr(1+5)=One
      scr(1+10)=One
      scr(2)=fv1(1)-This%ret0
      scr(2+5)=fv1(2)-This%ret0
      scr(2+10)=cc-This%cte0
      damp=dexp(-This%a1*(scr(2)**2+scr(7)**2)-This%a2*(scr(12)**2))
      damp2=damp*Two
      xx=-This%a2*scr(12)*damp2
      dm1=-This%a1*scr(2)*damp2+xx*dcc(1)
      dm2=-This%a1*scr(7)*damp2+xx*dcc(2)
      dm3=xx*dcc(3)
        
      do k = 3,5
      
        k1  = k
        k1m = k1 - 1 
        scr(k1) = scr(k1m) * scr(2)
        k2  = k + 5
        k2m = k2 - 1
        scr(k2) = scr(k2m) * scr(2+5)
        k3  = k + 10
        k3m = k3 - 1 
        scr(k3) = scr(k3m) * scr(2+10)
         
      end do
        
      voco=Zero
        
      do k = 1,This%nfcnt
        
        term0 = (scr(This%iqn(1,k)) * scr(This%iqn(2,k)+5) + scr(This%iqn(1,k)+5) * scr(This%iqn(2,k)))
        term  = scr(This%iqn(3,k)+10) * term0
        term3 = dfloat(This%iqn(3,k)-1) * scr(int(max(This%iqn(3,k)-One,One)+10.0_rkp)) * term0
        term1 = scr(This%iqn(3,k)+10) * (dfloat(This%iqn(1,k)-1)*scr(int(max(This%iqn(1,k)-One,One))  ) * scr(This%iqn(2,k)+5) + &
                   dfloat(This%iqn(2,k)-1) * scr(This%iqn(1,k)+5) * scr(int(max(This%iqn(2,k)-One,One))  ))
        term2 = scr(This%iqn(3,k)+10) * (dfloat(This%iqn(1,k)-1)*scr(int(max(This%iqn(1,k)-One,One)+5)) * scr(This%iqn(2,k)  ) + &
                   dfloat(This%iqn(2,k)-1) * scr(This%iqn(1,k)  ) * scr(int(max(This%iqn(2,k)-One,One))+5))
       
        voco = voco + This%qff(k) * term
         
      end do
      
      voco0 = voco
      voco  = v0o + damp * voco
        
    else
        
      voco  = v0o
      voco0 = Zero
         
    end if
        
    ham(4,4) = voco
    ham(1,4) = This%vcupt0 * ex5
    ham(2,4) = This%vcupt0 * ex5
    ham(4,2) = This%vcupt0 * ex5
    ham(4,1) = This%vcupt0 * ex5
    xx = -This%xm(5) * (ham(1,4)-This%vcupt0 * 1.e-7_rkp)
  
    ham(3,4)    = This%vcupt1 * ex7
    
    xx = -This%xm(7) * (ham(3,4)-This%vcupt1 * 1.e-7_rkp)
    
    ham(4,3)    = ham(3,4)
        
  end if
  
   
      
  do i = 1,ndim
    do j = 1,ndim
      ham2(j,i) = ham(j,i)
     end do
  end do
      
  call rs(4,ndim,ham,eig,1,ham,fv1,fv2,ierr)
      
  V = eig(1)
      
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
