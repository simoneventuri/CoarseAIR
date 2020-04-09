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

Module COAr_NASA_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Four, Six, Eight
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    COAr_NASA_PES_Type

  Type    ,extends(PES_Type)          ::    COAr_NASA_PES_Type
  
    integer                                                      :: idp = 150
  
    integer                                                      :: iCO
    integer                                                      :: iCAr
    integer                                                      :: iOAr
  
    integer                                                      :: i1
    integer                                                      :: i2
    integer                                                      :: i3
    integer                                                      :: nx
    real(rkp)                                                    :: a
    real(rkp)                                                    :: b
    real(rkp)                                                    :: c6c
    real(rkp)                                                    :: c6o
    real(rkp)                                                    :: rvib0
    real(rkp)                                                    :: galp
    real(rkp)                                                    :: g0
    real(rkp)                                                    :: dgo
    real(rkp)                                                    :: damp6
    real(rkp)                                                    :: vpdamp
    real(rkp)                                                    :: govrl
    real(rkp)                                                    :: gdel
    real(rkp)                                                    :: alfo
    real(rkp)                                                    :: alfc
    
    integer        ,dimension(2,150)                             :: ixpd         !!!!!!!!! <<<<< Because idp = 150
    real(rkp)      ,dimension(150)                               :: par          !!!!!!!!! <<<<< Because idp = 150
    
  contains
    procedure          ::  Initialize        =>    Initialize_COAr_NASA_PES
    procedure          ::  Output            =>    Output_COAr_NASA_PES
    procedure          ::  Compute           =>    Compute_COAr_NASA_PES_1d
    procedure          ::  Potential         =>    COAr_NASA_Potential_From_R
    procedure ,private ::  Compute_COAr_NASA_PES_1d_NoTraj
    procedure ,private ::  COAr_NASA_Potential_From_R_NoTraj
  End Type
  
  real(rkp)      ,dimension(10,2)                              :: fmat 
  real(rkp)      ,dimension(10,2,3)                            :: dfmat
  real(rkp)      ,dimension(20)                                :: parts2
  
  logical                                        ,parameter    :: i_Debug_Global = .False.
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NASA PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_COAr_NASA_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                  ,only:  Input_Type
  use Atom_Class                   ,only:  Atom_Type
  use CO_DiatomicPotential_Class   ,only:  CO_DiatomicPotential_Type
  
  class(COAr_NASA_PES_Type)                 ,intent(out)    :: This
  type(Input_Type)                          ,intent(in)     :: Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     :: Atoms 
  integer                                   ,intent(in)     :: iPES
  logical                         ,optional ,intent(in)     :: i_Debug
  

  character(80)                                             :: line
  integer                                                   :: iP
  character(*)                    ,parameter                :: Name_PES = 'COAr_NASA'
  character(:)       ,allocatable                           :: NASA_COAr_file
  integer                                                   :: Unit
  integer                                                   :: Status
  integer                                                   :: ipp
  integer                                                   :: nf
  integer            ,dimension(6)                          :: MatrixTemp = (/ 1, 2, 1, 3, 2, 3 /)
  type(CO_DiatomicPotential_Type)                           :: CO_DiatPot
  integer         ,dimension(3,2)                           :: iA

  logical                                                   :: i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize COAr_NASA PES")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
    
  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  
  do iP = 1,This%NPairs
    if ( ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+1))) ) .ne. "C" ) .and. ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+2)))) .ne. "O" ) ) then
      allocate( CO_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      This%iCO = iP
      if (i_Debug_Loc) call Logger%Write( 'CO pair =',iP )
    elseif ( ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+1))) ) .ne. "C" ) .and. ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+2)))) .ne. "Ar" ) ) then
      !allocate( CAr_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      This%iCAr = iP
      if (i_Debug_Loc) call Logger%Write( 'CAr pair =',iP )
    elseif ( ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+1))) ) .ne. "O" ) .and. ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+2)))) .ne. "Ar" ) ) then
      !allocate( OAr_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      This%iOAr = iP
      if (i_Debug_Loc) call Logger%Write( 'OAr pair =',iP )
    end if
  end do
  
  
  NASA_COAr_file = trim(adjustl(Input%DtbPath))  // '/Systems/' // trim(adjustl(Input%System)) // '/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '.inp'
  if (i_Debug_Loc) call Logger%Write( "Reading NASA COAr PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", NASA_COAr_file)
  open( File=NASA_COAr_file, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // NASA_COAr_file )
  
  read(Unit,'(A80)') line
  if (i_Debug_Loc) call Logger%Write( "Line:    ", line )

  read(Unit,*) This%c6c
  if (i_Debug_Loc) call Logger%Write( "c6c    = ", This%c6c )

  read(Unit,*) This%c6o
  if (i_Debug_Loc) call Logger%Write( "c6o    = ", This%c6o )
  
  read(Unit,*) This%rvib0
  if (i_Debug_Loc) call Logger%Write( "rvib0  = ", This%rvib0 )
  
  read(Unit,*) This%galp, This%g0
  if (i_Debug_Loc) call Logger%Write( "galp   = ", This%galp )
  if (i_Debug_Loc) call Logger%Write( "g0     = ", This%g0 )
  
  read(Unit,*) This%dgo
  if (i_Debug_Loc) call Logger%Write( "dgo    = ", This%dgo )
  
  read(Unit,*) This%damp6
  if (i_Debug_Loc) call Logger%Write( "damp6  = ", This%damp6 )
  This%damp6 = This%damp6**6
  
  read(Unit,'(A80)') line                                                    
  if (line(2:9) .eq. 'vpdamp =') then                                  
    read(line(10:80),*) This%vpdamp 
    if (i_Debug_Loc) call Logger%Write( "vpdamp = ", This%vpdamp )                                      
  end if  
                                                           
  This%govrl = log(This%dgo)  
  !if (i_Debug_Loc) call Logger%Write( "govrl  = ", This%govrl )                                      
                                                                                                                        
  read(Unit,*) This%par(1), This%par(2), This%par(3)
  if (i_Debug_Loc) call Logger%Write( "par    = ", This%par )         
                               
  This%gdel = sqrt(- Two * This%govrl / This%galp)                                      
  !if (i_Debug_Loc) call Logger%Write( "gdel  = ", This%gdel )         
  
  This%alfo = This%par(2)**2
  This%alfc = This%par(3)**2
  
  ipp = 6
  nf  = 0
  Status = 0
  do while (Status == 0) 
    
    read(Unit,*,iostat=Status) This%i1, This%i2, This%i3, This%a, This%b
    if (i_Debug_Loc) call Logger%Write( "i1     = ", This%i1 )
    if (i_Debug_Loc) call Logger%Write( "i2     = ", This%i2 )
    if (i_Debug_Loc) call Logger%Write( "i3     = ", This%i3 )
    if (i_Debug_Loc) call Logger%Write( "a      = ", This%a )
    if (i_Debug_Loc) call Logger%Write( "b      = ", This%b )
    
    nf = nf + 1
    if (ipp+2 .gt. This%idp) stop 'idp'
    This%ixpd(1,nf)  = This%i2 + 1
    This%ixpd(2,nf)  = This%i3 + 1
    This%par(ipp)    = This%a
    This%par(ipp+1)  = This%b
    ipp              = ipp+2
    This%nx          = This%i1
    
  end do 
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_COAr_NASA_PES( This, Unit )

  class(COAr_NASA_PES_Type)               ,intent(in)     :: This
  integer                                 ,intent(in)     :: Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('COAr PES')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function COAr_NASA_Potential_From_R( This, R, Q ) result( V )

  class(COAr_NASA_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  integer                                                    :: iCO
  integer                                                    :: iCAr
  integer                                                    :: iOAr
  
  iCO  = This%iCO
  iCAr = This%iCAr
  iOAr = This%iOAr
      
  call COAr_NASA_Potential_From_R_NoTraj( This, R(iCO), R(iCAr), R(iOAr), V )

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________!
!Function COAr_NASA_Potential_From_R( This, R ) result( V )

!  class(COAr_NASA_PES_Type)               ,intent(in)     :: This
!  real(rkp) ,dimension(:,:)  CONTIGUOUS   ,intent(in)     :: R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs,NTraj)
!  real(rkp) ,dimension( size(R,2) )                       :: V           !< Potential energy in [hartree]. Dim=(NTraj)

!  integer                                                 :: iTraj
!  integer                                                 :: iCO
!  integer                                                 :: iCAr
!  integer                                                 :: iOAr
!  
!  iCO  = This%iCO
!  iCAr = This%iCAr
!  iOAr = This%iOAr
!  
!  do iTraj = 1,size(V)
!      
!    call COAr_NASA_Potential_From_R_NoTraj( This, R(iCO,iTraj), R(iCAr,iTraj), R(iOAr,iTraj), V(iTraj) )
!    
!  end do

!End Function
!!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_COAr_NASA_PES_1d( This, R, Q, V, dVdR, dVdQ )
    
  class(COAr_NASA_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)
      
  integer                                                      :: iCO
  integer                                                      :: iCAr
  integer                                                      :: iOAr
  
  iCO  = This%iCO
  iCAr = This%iCAr
  iOAr = This%iOAr

  dVdQ         = Zero
  call Compute_COAr_NASA_PES_1d_NoTraj( This, R(iCO), R(iCAr), R(iOAr), V, dVdR(:) )
 
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_COAr_NASA_PES_1d_NoTraj( This, R_CO, R_CAr, R_OAr, V, dV, i_Debug )
  
  use CO_DiatomicPotential_Class   ,only:  CO_DiatomicPotential_Type
  
  class(COAr_NASA_PES_Type)                 ,intent(in)     :: This
  real(rkp)                                 ,intent(in)     :: R_CO
  real(rkp)                                 ,intent(in)     :: R_CAr
  real(rkp)                                 ,intent(in)     :: R_OAr
  real(rkp)                                 ,intent(out)    :: V
  real(rkp)       ,dimension(3)             ,intent(out)    :: dV
  logical                         ,optional ,intent(in)     :: i_Debug
  
  type(CO_DiatomicPotential_Type)                           :: CO_DiatPot
  
  integer                                                   :: iCO
  integer                                                   :: iCAr
  integer                                                   :: iOAr
  
  real(rkp)                                                 :: rca
  real(rkp)                                                 :: roa
  real(rkp)                                                 :: rvib
  real(rkp)                                                 :: rvibh
  real(rkp)                                                 :: cosa 
  real(rkp)                                                 :: rbig2
  real(rkp)                                                 :: rbig
  real(rkp)                                                 :: dot
  real(rkp)                                                 :: varo
  real(rkp)                                                 :: dvaro
  real(rkp)                                                 :: varc
  real(rkp)                                                 :: dvarc
  real(rkp)                                                 :: da1
  real(rkp)                                                 :: da2
  real(rkp)                                                 :: da3
  real(rkp)                                                 :: rbigi
  real(rkp)                                                 :: dr1
  real(rkp)                                                 :: dr2
  real(rkp)                                                 :: dr3
  real(rkp)                                                 :: ex1
  real(rkp)                                                 :: dex1
  real(rkp)                                                 :: dex11
  real(rkp)                                                 :: dex12
  real(rkp)                                                 :: dex13
  real(rkp)                                                 :: ex2
  real(rkp)                                                 :: dex2
  real(rkp)                                                 :: dex21
  real(rkp)                                                 :: dex22
  real(rkp)                                                 :: dex23
  real(rkp)                                                 :: exr
  real(rkp)                                                 :: dexr
  real(rkp)                                                 :: ddot1
  real(rkp)                                                 :: ddot2
  real(rkp)                                                 :: ddot3
  real(rkp)                                                 :: term1
  real(rkp)                                                 :: term2           
  real(rkp)                                                 :: dterm11
  real(rkp)                                                 :: dterm12
  real(rkp)                                                 :: dterm13
  real(rkp)                                                 :: dterm21
  real(rkp)                                                 :: dterm22
  real(rkp)                                                 :: dterm23
  real(rkp)                                                 :: vlr
  real(rkp)                                                 :: dvlr1
  real(rkp)                                                 :: dvlr2
  real(rkp)                                                 :: dvlr3
  real(rkp)                                                 :: farg
  real(rkp)                                                 :: dfarg1
  real(rkp)                                                 :: dfarg2
  real(rkp)                                                 :: dfarg3
  real(rkp)                                                 :: xxe
  real(rkp)                                                 :: xx
  real(rkp)                                                 :: xx2
  real(rkp)                                                 :: dxx1
  real(rkp)                                                 :: dxx2
  real(rkp)                                                 :: dxx3
  real(rkp)                                                 :: fix
  real(rkp)                                                 :: dfix1
  real(rkp)                                                 :: dfix2
  real(rkp)                                                 :: dfix3
   
  integer                                                   :: j, jm, ip
  
  iCO  = This%iCO
  iCAr = This%iCAr
  iOAr = This%iOAr
  
  dvlr1 = Zero
                     
  call CO_DiatPot%Compute_Vd_dVd( R_CO, V, dV(iCO) )

  rca   = R_CAr
  roa   = R_OAr
  rvib  = R_CO
  rvibh = Half * rvib 
  cosa  = Half * (rca**2 + rvib**2 - roa**2) / (rca*rvib)
  rbig2 = rca**2 + rvibh**2 - rca * rvib * cosa
  rbig  = sqrt(rbig2)
  dot   = (rvibh**2 + rbig2 - rca**2) / (rvib * rbig)

  dV(iCAr) = Zero
  dV(iOAr) = Zero

  varo  = 18.0_rkp * Eight * exp(- This%alfo * roa) / roa
  dvaro = - varo * (This%alfo + (One / roa))

  varc  = 18.0_rkp * Six * exp(- This%alfc * rca) / rca
  dvarc = - varc * (This%alfc + (One / rca))

  da1 = (One /rca)  - (cosa / rvib)
  da2 = (One /rvib) - (cosa / rca)
  da3 = - roa / (rca * rvib)

  rbigi = Half / rbig
  dr1   = rbigi * (Half * rvib - rca * (cosa + rvib * da1))
  dr2   = rbigi * (Two * rca  - rvib * (cosa + rca * da2))
  dr3   = rbigi * (- rca * rvib * da3)

  ex1   = exp(- This%galp * ((rbig - This%g0)**2))
  dex1  = - This%galp * Two * (rbig - This%g0) *ex1
  dex11 = dex1 * dr1
  dex12 = dex1 * dr2
  dex13 = dex1 * dr3
  ex2   = exp(- This%galp * ((rbig - This%g0 - This%gdel)**2))
  dex2  = - This%galp * Two * (rbig - This%g0 - This%gdel) * ex2
  dex21 = dex2 * dr1
  dex22 = dex2 * dr2
  dex23 = dex2 * dr3
  exr   =  exp( - This%par(1) * ((rvib - This%rvib0)**2))
  dexr  = - Two * This%par(1) *  (rvib - This%rvib0) * exr

  term1 = Zero
  term2 = Zero

  ddot1 =  (Two * dr1 / rvib) - (dot * dr1 / rbig) + (Half / rbig) - (dot / rvib)
  ddot2 = -(Two * dr2 / rvib) - (One / rbig) * dot * dr2
  ddot3 =  (Two * dr3 / rvib) - (dot * dr3 / rbig)

  fmat(1,1)    = One 
  fmat(1,2)    = One 
  dfmat(1,1,1) = Zero
  dfmat(1,1,2) = Zero
  dfmat(1,1,3) = Zero
  dfmat(1,2,1) = Zero
  dfmat(1,2,2) = Zero
  dfmat(1,2,3) = Zero
  fmat(2,1)    = rvib
  dfmat(2,1,1) = One 
  dfmat(2,1,2) = Zero
  dfmat(2,1,3) = Zero
  fmat(2,2)    = dot
  dfmat(2,2,1) = ddot1
  dfmat(2,2,2) = ddot2
  dfmat(2,2,3) = ddot3
  
  do j = 3,9
    jm           = j - 1
    
    fmat(j,1)    = fmat(jm,1)    * fmat(2,1)
    dfmat(j,1,1) = dfmat(jm,1,1) * fmat(2,1) + fmat(jm,1) * dfmat(2,1,1)
    fmat(j,2)    = fmat(jm,2)    * fmat(2,2)
    dfmat(j,2,1) = dfmat(jm,2,1) * fmat(2,2) + fmat(jm,2) * dfmat(2,2,1)
    dfmat(j,2,2) = dfmat(jm,2,2) * fmat(2,2) + fmat(jm,2) * dfmat(2,2,2)
    dfmat(j,2,3) = dfmat(jm,2,3) * fmat(2,2) + fmat(jm,2) * dfmat(2,2,3)
    
  end do
  
  ip = 6
  dterm11 = Zero
  dterm12 = Zero
  dterm13 = Zero
  dterm21 = Zero
  dterm22 = Zero
  dterm23 = Zero
  do j = 1,This%nx
  
    term1   = term1   + This%par(ip) *  fmat(This%ixpd(1,j),1) *  fmat(This%ixpd(2,j),2)
    dterm11 = dterm11 + This%par(ip) * (fmat(This%ixpd(2,j),2) * dfmat(This%ixpd(1,j),1,1) + fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,1))
    dterm12 = dterm12 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,2)
    dterm13 = dterm13 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,3)

    ip      = ip + 1
    term2   = term2   + This%par(ip) *  fmat(This%ixpd(1,j),1) *  fmat(This%ixpd(2,j),2)
    dterm21 = dterm21 + This%par(ip) * (fmat(This%ixpd(2,j),2) * dfmat(This%ixpd(1,j),1,1) + fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,1))
    dterm22 = dterm22 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,2)
    dterm23 = dterm23 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,3)

    ip = ip + 1
  end do

  vlr   = (This%c6c / (rca**6 + This%damp6)) + (This%c6o / (roa**6 + This%damp6))
  dvlr2 = -Six * This%c6c * (rca**5) / ((rca**6 + This%damp6)**2)
  dvlr3 = -Six * This%c6o * (roa**5) / ((roa**6 + This%damp6)**2)

  parts2(1) = exr
  parts2(2) = ex1
  parts2(3) = term1
  parts2(4) = ex2
  parts2(5) = term2
  farg      = (ex1 * term1 + ex2 * term2) * exr                                   
  xxe       = exp(- This%vpdamp * farg)                                            
  xx        = One / (One + xxe)                                                 
  xx2       = xx * xx * This%vpdamp * xxe                                             
  fix       = One +xx * farg                                                  

  dfarg1 = exr  * (dex11 * term1 + ex1 * dterm11 + dex21 * term2 + ex2 * dterm21) + dexr * (ex1 * term1 + ex2 * term2)
  dfarg2 = exr  * (dex12 * term1 + ex1 * dterm12 + dex22 * term2 + ex2 * dterm22)     
  dfarg3 = exr  * (dex13 * term1 + ex1 * dterm13 + dex23 * term2 + ex2 * dterm23)     
  dxx1   = xx2  * dfarg1                                                  
  dxx2   = xx2  * dfarg2                                                  
  dxx3   = xx2  * dfarg3                                                  
  dfix1  = dxx1 * farg + xx * dfarg1                                        
  dfix2  = dxx2 * farg + xx * dfarg2                                        
  dfix3  = dxx3 * farg + xx * dfarg3                                        

  V        = V           + (varo + varc) * fix   + vlr
  
  dV(iCO)  = dV(iCO)     + (varo + varc) * dfix1 + dvlr1
  dV(iCAr) = dvarc * fix + (varo + varc) * dfix2 + dvlr2
  dV(iOAr) = dvaro * fix + (varo + varc) * dfix3 + dvlr3
    
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine COAr_NASA_Potential_From_R_NoTraj( This, R_CO, R_CAr, R_OAr, V, i_Debug )
  
  use CO_DiatomicPotential_Class   ,only:  CO_DiatomicPotential_Type
  
  class(COAr_NASA_PES_Type)                 ,intent(in)     :: This
  real(rkp)                                 ,intent(in)     :: R_CO
  real(rkp)                                 ,intent(in)     :: R_CAr
  real(rkp)                                 ,intent(in)     :: R_OAr
  real(rkp)                                 ,intent(out)    :: V
  logical                         ,optional ,intent(in)     :: i_Debug
  
  type(CO_DiatomicPotential_Type)                           :: CO_DiatPot
  
  integer                                                   :: iCO
  integer                                                   :: iCAr
  integer                                                   :: iOAr
  
  real(rkp)                                                 :: rca
  real(rkp)                                                 :: roa
  real(rkp)                                                 :: rvib
  real(rkp)                                                 :: rvibh
  real(rkp)                                                 :: cosa 
  real(rkp)                                                 :: rbig2
  real(rkp)                                                 :: rbig
  real(rkp)                                                 :: dot
  real(rkp)                                                 :: varo
  real(rkp)                                                 :: dvaro
  real(rkp)                                                 :: varc
  real(rkp)                                                 :: dvarc
  real(rkp)                                                 :: da1
  real(rkp)                                                 :: da2
  real(rkp)                                                 :: da3
  real(rkp)                                                 :: rbigi
  real(rkp)                                                 :: dr1
  real(rkp)                                                 :: dr2
  real(rkp)                                                 :: dr3
  real(rkp)                                                 :: ex1
  real(rkp)                                                 :: dex1
  real(rkp)                                                 :: dex11
  real(rkp)                                                 :: dex12
  real(rkp)                                                 :: dex13
  real(rkp)                                                 :: ex2
  real(rkp)                                                 :: dex2
  real(rkp)                                                 :: dex21
  real(rkp)                                                 :: dex22
  real(rkp)                                                 :: dex23
  real(rkp)                                                 :: exr
  real(rkp)                                                 :: dexr
  real(rkp)                                                 :: ddot1
  real(rkp)                                                 :: ddot2
  real(rkp)                                                 :: ddot3
  real(rkp)                                                 :: term1
  real(rkp)                                                 :: term2           
  real(rkp)                                                 :: dterm11
  real(rkp)                                                 :: dterm12
  real(rkp)                                                 :: dterm13
  real(rkp)                                                 :: dterm21
  real(rkp)                                                 :: dterm22
  real(rkp)                                                 :: dterm23
  real(rkp)                                                 :: vlr
  real(rkp)                                                 :: dvlr1
  real(rkp)                                                 :: dvlr2
  real(rkp)                                                 :: dvlr3
  real(rkp)                                                 :: farg
  real(rkp)                                                 :: dfarg1
  real(rkp)                                                 :: dfarg2
  real(rkp)                                                 :: dfarg3
  real(rkp)                                                 :: xxe
  real(rkp)                                                 :: xx
  real(rkp)                                                 :: xx2
  real(rkp)                                                 :: dxx1
  real(rkp)                                                 :: dxx2
  real(rkp)                                                 :: dxx3
  real(rkp)                                                 :: fix
  real(rkp)                                                 :: dfix1
  real(rkp)                                                 :: dfix2
  real(rkp)                                                 :: dfix3
  real(rkp)                                                 :: dVTemp
  
  integer                                                   :: j, jm, ip
  
  iCO  = This%iCO
  iCAr = This%iCAr
  iOAr = This%iOAr
  
  dvlr1 = Zero                     

  call CO_DiatPot%Compute_Vd_dVd( R_CO, V, dVTemp )

  rca   = R_CAr
  roa   = R_OAr
  rvib  = R_CO
  rvibh = Half * rvib 
  cosa  = Half * (rca**2 + rvib**2 - roa**2) / (rca*rvib)
  rbig2 = rca**2 + rvibh**2 - rca * rvib * cosa
  rbig  = sqrt(rbig2)
  dot   = (rvibh**2 + rbig2 - rca**2) / (rvib * rbig)

  varo  = 18.0_rkp * Eight * exp(- This%alfo * roa) / roa
  dvaro = - varo * (This%alfo + (One / roa))

  varc  = 18.0_rkp * Six * exp(- This%alfc * rca) / rca
  dvarc = - varc * (This%alfc + (One / rca))

  da1 = (One /rca)  - (cosa / rvib)
  da2 = (One /rvib) - (cosa / rca)
  da3 = - roa / (rca * rvib)

  rbigi = Half / rbig
  dr1   = rbigi * (Half * rvib - rca * (cosa + rvib * da1))
  dr2   = rbigi * (Two * rca  - rvib * (cosa + rca * da2))
  dr3   = rbigi * (- rca * rvib * da3)

  ex1   = exp(- This%galp * ((rbig - This%g0)**2))
  dex1  = - This%galp * Two * (rbig - This%g0) *ex1
  dex11 = dex1 * dr1
  dex12 = dex1 * dr2
  dex13 = dex1 * dr3
  ex2   = exp(- This%galp * ((rbig - This%g0 - This%gdel)**2))
  dex2  = - This%galp * Two * (rbig - This%g0 - This%gdel) * ex2
  dex21 = dex2 * dr1
  dex22 = dex2 * dr2
  dex23 = dex2 * dr3
  exr   =  exp( - This%par(1) * ((rvib - This%rvib0)**2))
  dexr  = - Two * This%par(1) *  (rvib - This%rvib0) * exr

  term1 = Zero
  term2 = Zero

  ddot1 =  (Two * dr1 / rvib) - (dot * dr1 / rbig) + (Half / rbig) - (dot / rvib)
  ddot2 = -(Two * dr2 / rvib) - (One / rbig) * dot * dr2
  ddot3 =  (Two * dr3 / rvib) - (dot * dr3 / rbig)

  fmat(1,1)    = One 
  fmat(1,2)    = One 
  dfmat(1,1,1) = Zero
  dfmat(1,1,2) = Zero
  dfmat(1,1,3) = Zero
  dfmat(1,2,1) = Zero
  dfmat(1,2,2) = Zero
  dfmat(1,2,3) = Zero
  fmat(2,1)    = rvib
  dfmat(2,1,1) = One 
  dfmat(2,1,2) = Zero
  dfmat(2,1,3) = Zero
  fmat(2,2)    = dot
  dfmat(2,2,1) = ddot1
  dfmat(2,2,2) = ddot2
  dfmat(2,2,3) = ddot3
  
  do j = 3,9
    jm           = j - 1
    
    fmat(j,1)    = fmat(jm,1)    * fmat(2,1)
    dfmat(j,1,1) = dfmat(jm,1,1) * fmat(2,1) + fmat(jm,1) * dfmat(2,1,1)
    fmat(j,2)    = fmat(jm,2)    * fmat(2,2)
    dfmat(j,2,1) = dfmat(jm,2,1) * fmat(2,2) + fmat(jm,2) * dfmat(2,2,1)
    dfmat(j,2,2) = dfmat(jm,2,2) * fmat(2,2) + fmat(jm,2) * dfmat(2,2,2)
    dfmat(j,2,3) = dfmat(jm,2,3) * fmat(2,2) + fmat(jm,2) * dfmat(2,2,3)
    
  end do
  
  ip = 6
  dterm11 = Zero
  dterm12 = Zero
  dterm13 = Zero
  dterm21 = Zero
  dterm22 = Zero
  dterm23 = Zero
  do j = 1,This%nx
  
    term1   = term1   + This%par(ip) *  fmat(This%ixpd(1,j),1) *  fmat(This%ixpd(2,j),2)
    dterm11 = dterm11 + This%par(ip) * (fmat(This%ixpd(2,j),2) * dfmat(This%ixpd(1,j),1,1) + fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,1))
    dterm12 = dterm12 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,2)
    dterm13 = dterm13 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,3)

    ip      = ip + 1
    term2   = term2   + This%par(ip) *  fmat(This%ixpd(1,j),1) *  fmat(This%ixpd(2,j),2)
    dterm21 = dterm21 + This%par(ip) * (fmat(This%ixpd(2,j),2) * dfmat(This%ixpd(1,j),1,1) + fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,1))
    dterm22 = dterm22 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,2)
    dterm23 = dterm23 + This%par(ip) *  fmat(This%ixpd(1,j),1) * dfmat(This%ixpd(2,j),2,3)

    ip = ip + 1
  end do

  vlr   = (This%c6c / (rca**6 + This%damp6)) + (This%c6o / (roa**6 + This%damp6))
  dvlr2 = -Six * This%c6c * (rca**5) / ((rca**6 + This%damp6)**2)
  dvlr3 = -Six * This%c6o * (roa**5) / ((roa**6 + This%damp6)**2)

  parts2(1) = exr
  parts2(2) = ex1
  parts2(3) = term1
  parts2(4) = ex2
  parts2(5) = term2
  farg      = (ex1 * term1 + ex2 * term2) * exr                                   
  xxe       = exp(- This%vpdamp * farg)                                            
  xx        = One / (One + xxe)                                                 
  xx2       = xx * xx * This%vpdamp * xxe                                             
  fix       = One +xx * farg                                                  

  dfarg1 = exr  * (dex11 * term1 + ex1 * dterm11 + dex21 * term2 + ex2 * dterm21) + dexr * (ex1 * term1 + ex2 * term2)
  dfarg2 = exr  * (dex12 * term1 + ex1 * dterm12 + dex22 * term2 + ex2 * dterm22)     
  dfarg3 = exr  * (dex13 * term1 + ex1 * dterm13 + dex23 * term2 + ex2 * dterm23)     
  dxx1   = xx2  * dfarg1                                                  
  dxx2   = xx2  * dfarg2                                                  
  dxx3   = xx2  * dfarg3                                                  
  dfix1  = dxx1 * farg + xx * dfarg1                                        
  dfix2  = dxx2 * farg + xx * dfarg2                                        
  dfix3  = dxx3 * farg + xx * dfarg3                                        

  V = V + (varo + varc) * fix + vlr

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
