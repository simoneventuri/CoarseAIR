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

Module O4_UMN_PES_Class
!**********************************************************************
!   System:                     O4
!   Functional form:            permutation-invariant polynomials
!   Common name:                O4(adiabatic ground state)
!   Number of derivatives:      1
!   Number of bodies:           4
!   Number of electronic surfaces: 1
!   Interface: Section-2
!
!   References:: Bender, Valentini, Nompelis, Paukku, Varga, Truhlar,
!                Schwartzentruber, and Candler, Journal of Chemical 
!                Physics, submitted. Manuscript no. A15.02.0237
!
!   Notes:    PES of O4 with special emphasize for
!             N2 + N2 --> N2 + N + N
!            - New fit based on extended dataset
!             
!            - Instead of Morse variables
!              mixed-exponential- gaussian (MEG)
!              variable is applied to describe the 
!              long-range interactions in a better way  
!             (named O4pes-gpip-meg)
!
!     N1--N2
!
!     N3--O4
!
!
!
!   Input: X(4),Y(4),Z(4)               in units of bohr
!   Output: E                           in units of hartree
!   Output: dEdX(4),dEdY(4),dEdZ(4)     hartree/bohr
!**********************************************************************

!**********************************************************************
! This code is based on O4pes-gpip-meg.f (FORTRAN77 version) where some  
! keywords were repleced by more modern ones.
! The surface itself was not changed.
!**********************************************************************

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten, Pi, B_To_Ang, Kcm_To_Hartree
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    O4_UMN_PES_Type


  Type, extends(PES_Type) :: O4_UMN_PES_Type
    real(rkp) ::  shiftof0
    ! Nonlinear parameters:
    ! a(in Ang)
    ! ab (in Ang^2)
    ! ra (in Ang)
    ! and rb (in Ang)
    real(rkp) :: a
    real(rkp) :: ab
    real(rkp) :: ra
    real(rkp) :: rb
    ! Reference energy of infinitely separated O2 + O2 in hartree (taken from calculations)
    real(rkp) :: Eref
    ! For O2 + O2 framework total diss. energy is 2*120.243 kcal/mol 
    real(rkp) :: totdiss
    ! Linear parameters optimized by the weighted-least square fitting
    real(rkp) ,dimension(430) :: C
  contains
    procedure ::  Initialize     =>    Initialize_O4_UMN_PES
    procedure ::  Output         =>    Output_O4_UMN_PES
    procedure ::  Compute        =>    Compute_O4_UMN_PES_1d
    procedure ::  Potential      =>    O4_UMN_Potential_From_R
    procedure ::  TriatPotential =>    O4_UMN_Potential_From_R_OnlyTriat
  End Type

  logical                 ,parameter      :: i_Debug_Global = .False.

  real(rkp) ,dimension(0:111)   :: rM       ! Array to store monomials
  real(rkp) ,dimension(0:465)   :: P        ! Array to store polynomials
  real(rkp) ,dimension(430)     :: B        ! Array to store basis functions
  real(rkp) ,dimension(6,0:111) :: dMdR     ! The derivative of monomials w.r.t. R
  real(rkp) ,dimension(6,0:465) :: dPdR     ! The derivative of basis functions w.r.t. R
  real(rkp) ,dimension(6,430)   :: dBdR     ! The derivative of B w.r.t. R 

  contains

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for UMN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_O4_UMN_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(O4_UMN_PES_Type)                    ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms  
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP
  character(*)                    ,parameter                ::    Name_PES = 'O4_UMN'
  real(rkp)                                                 ::    Temp
  real(rkp)                                                 ::    veccsd
  real(rkp)                                                 ::    veuse
  integer         ,dimension(6,2)                           ::    iA
  type(DiatomicPotential_Factory_Type)                      ::    DiatPotFactory
  character(150)                                            ::    Weights_File
  integer                                                   ::    Unit, Status, i
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O4_UMN_PES" )
  !i_Debug_Loc   =     Logger%On()
    

  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   6               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  iA(1,:) = [1, 2]
  iA(2,:) = [1, 3]
  iA(3,:) = [1, 4]
  iA(4,:) = [2, 3]
  iA(5,:) = [2, 4]
  iA(6,:) = [3, 4]  

  ! allocate( This%mMiMn(4) )
  ! This%mMiMn(1:3) = - Atoms(1:3)%Mass / Atoms(4)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )


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


  ! ==============================================================================================================
  !   READ PARAMETERS
  ! ============================================================================================================== 
  Weights_File = trim(adjustl(Input%DtbPath))  // '/Systems/O4/PESs/UMN/' // trim(adjustl(Input%PES_Model(iPES))) // '.dat'
  if (i_Debug_Loc) call Logger%Write( "Reading O4 UMN PES Parameters" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", Weights_File)
  open( File=Weights_File, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // Weights_File )

    read(Unit,*) This%a 
    read(Unit,*) This%ab
    read(Unit,*) This%ra
    read(Unit,*) This%rb
    read(Unit,*) This%Eref
    read(Unit,*) This%totdiss
    if (i_Debug_Loc) call Logger%Write( "This%a       = ", This%a  )
    if (i_Debug_Loc) call Logger%Write( "This%ab      = ", This%ab  )
    if (i_Debug_Loc) call Logger%Write( "This%ra      = ", This%ra  )
    if (i_Debug_Loc) call Logger%Write( "This%rb      = ", This%rb  )
    if (i_Debug_Loc) call Logger%Write( "This%Eref    = ", This%Eref  )
    if (i_Debug_Loc) call Logger%Write( "This%totdiss = ", This%totdiss  )
    do i=1,143
      read(Unit,*) This%C((i-1)*3+1), This%C((i-1)*3+2), This%C((i-1)*3+3)
    end do
    read(Unit,*) This%C(430)
    if (i_Debug_Loc) call Logger%Write( "This%C = ", This%C )

  close(Unit)
  ! ==============================================================================================================



  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_O4_UMN_PES( This, Unit )

  class(O4_UMN_PES_Type)                  ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('O4 PES  from UMN')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O4_UMN_Potential_From_R( This, R, Q ) result( V )

  class(O4_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp)                                                  ::    VTemp, VDiat, VTriat
  real(rkp) ,dimension(6)                                    ::    dVdRTriat
  integer                                                    ::    iP

  VDiat = Zero
  do iP=1,6
    VTemp = This%Pairs(iP)%Vd%DiatomicPotential(R(iP))
    VDiat = VDiat + VTemp
  end do

  call O4pes(This%a, This%ab, This%ra, This%rb, This%Eref, This%totdiss, This%C, R, VTriat, dVdRTriat, 0)
  V = VTriat + VDiat

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function O4_UMN_Potential_From_R_OnlyTriat( This, R, Q ) result( V )

  class(O4_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(6)                                    ::    dVdR
  
  call O4pes(This%a, This%ab, This%ra, This%rb, This%Eref, This%totdiss, This%C, R, V, dVdR, 0)
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_O4_UMN_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(O4_UMN_PES_Type)                        ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp)                                                  ::    VTemp, VDiat, VTriat
  real(rkp) ,dimension(6)                                    ::    dVdRDiat, dVdRTriat
  integer                                                    ::    iP

  dVdR = Zero

  VDiat = Zero
  do iP=1,6
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VTemp, dVdRDiat(iP) )
    VDiat = VDiat + VTemp
  end do

  call O4pes(This%a, This%ab, This%ra, This%rb, This%Eref, This%totdiss, This%C, R, VTriat, dVdRTriat, 1)

  V    = VTriat    + VDiat
  dVdR = dVdRTriat + dVdRDiat

  dVdQ = Zero
  call This%TransToCart_4Atoms( R, Q, dVdR, dVdQ)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



! ! **************************************************************************************************************
! ! **************************************************************************************************************
! !                                   PRIVATE PROCEDURES for UMN PES
! ! **************************************************************************************************************
! ! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!
!________________________________________________________________________________________________________________________________!
Subroutine O4pes(a, ab, ra, rb, Eref, totdiss, C, R, V, dVdR, igrad)
!Subroutine O4pes(X, V, dVdX, igrad)
!**********************************************************************
! Subroutine to calculate the potential energy V and gradient dVdX
! for given Cartesian coordinate X(12)  
! R:    Interatomic bond distance (6)
! V:    Calculated potential energy
! dVdX:   The derivative of V w.r.t. X, dim(12)
! dVdR:   The derivative of V w.r.t. R, dim(6) 
! dPdR:   The derivative of basis functions w.r.t. R
!   dim(6*306)
! dMdR:   The derivative of monomials w.r.t. R
!   dim(6*112)
! dRdX:   The derivative of R w.r.t. X, dim(6*12)
!**********************************************************************
  real(rkp)               ,intent(in)  :: a
  real(rkp)               ,intent(in)  :: ab
  real(rkp)               ,intent(in)  :: ra
  real(rkp)               ,intent(in)  :: rb
  real(rkp)               ,intent(in)  :: Eref
  real(rkp)               ,intent(in)  :: totdiss
  real(rkp) ,dimension(430),intent(in) :: C
  real(rkp) ,dimension(6) ,intent(in)  :: R
  !real(rkp) ,dimension(12) ,intent(in)  :: X
  integer                 ,intent(in)  :: igrad
  real(rkp)               ,intent(out) :: V
  real(rkp) ,dimension(6) ,intent(out) :: dVdR
  !real(rkp) ,dimension(12) ,intent(out)  :: dVdX

  integer                              :: i, j, k
  real(rkp) ,dimension(6)              :: RAng
  integer                              :: nob

  ! The soubroutines receives in input distances in Bohr and gives derivatives in dR
  RAng = R * B_To_Ang

  ! Read cartesian coordinate from input file
  !call coord_convt(X)

  if (igrad .le. 1) then
  ! Call subroutine EvV to evaluate potential energy V
    call EvV(a, ab, ra, rb, C, RAng, V)

    if (igrad .eq. 1) then
      ! Call EvdVdR to evaluate dVdR(6)
      Call EvdVdR(a, ab, ra, rb, C, RAng, dVdR)
      !! Call EvdVdX to evaluate the derivatives of V w.r.t. X
      !call evdvdx(X,dVdX)
    endif

  else
    write (*,*) 'Only igrad = 0, 1 is allowed!'
  endif

  ! Initialized v to be totdiss
  
  V    = V * Kcm_To_Hartree !+ Eref
  !V    = (V + totdiss) * Kcm_To_Hartree !+ Eref
  dVdR = dVdR * Kcm_To_Hartree * B_To_Ang

end Subroutine 
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine EvV(a, ab, ra, rb, C, R, V)
!**********************************************************************
! Subroutine to evaluate V for giving R 
! V(R) = C*P
! C:    Coefficients, stored in 'dim.inc' 
! P:    Basis functions evaluated for giving R
! rMs:    rMs(6), six mixed exponential gaussian terms (MEG)
! a:    Nonlinear parameters in Morse terms(Angstrom)
! ab:   Nonlinear parameters in Gauss terms(Angstrom^2)
! re:   Equilibrium bond length(Angstrom)
! nop:    number of points
! nom:    number of monomials
! nob:    number of basis functions(polynomials)
! rM(0:111):  Array to store monomials
! P(0:305): Array to store polynomials
! B(1:276):     Array to store basis functions
!**********************************************************************
  real(rkp)               ,intent(in)  :: a
  real(rkp)               ,intent(in)  :: ab
  real(rkp)               ,intent(in)  :: ra
  real(rkp)               ,intent(in)  :: rb
  real(rkp) ,dimension(430),intent(in) :: C
  real(rkp) ,dimension(6) ,intent(in)  :: R
  real(rkp)               ,intent(out) :: V

  integer                              :: i, j, k
  real(rkp)                            :: dist
  real(rkp)                            :: dV2dR
  real(rkp)                            :: V2
  real(rkp) ,dimension(6)              :: rms


  ! Evaluate 2-body interactions
  ! do i=1,6
  !   dist = R(i)
  !   call ev2gm2(dist, V2, dV2dR, 1, 0)
  !   V    = V + V2
  ! enddo


  ! Calculate the six MEG terms for each point
  call evmorse(a, ab, ra, rb, R, rms)

  ! Calculate the monomials for each point by using six MEG terms
  call evmono(rms)

  ! Calculate the polynomials (basis functions) by using monomials
  call evpoly 

  ! Calculate the basis functions by removing unconnected and 2-body terms
  call evbas

  ! Evaluate V by taken the product of C and Basis function array
  do i=1,430
    V = V + C(i)*b(i)
  enddo

  !      Write(*,9999) V 
  ! 9999 Format('The potential energy is ',F20.14,' kcal/mol')

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine EvdVdR(a, ab, ra, rb, C, R, dVdR)
!**********************************************************************
! Subroutine to evaluate dVdR for giving R 
! dVdR = dV2dR + C*dBdR
! C:    Coefficients, stored in 'dim.inc' 
! P:    Basis functions evaluated for giving R
! M:    Monomials evaluated for giving R
! dV2dR:        Gradient of 2-body interactions
! dMsdR:  dMsdR(6,6), 6 MEG terms w.r.t. 6 bond length
! dMdR:   dMdR(6,nom), nom monomials w.r.t.6 bond length
! dPdR:   dPdR(6,nob), nop polynomial basis functions 
!   w.r.t. 6 bond length
! nom:    number of monomials
! nob:    number of basis functions(polynomials)
! M(nom): Array to store monomials
! P(nob): Array to store polynomials
!**********************************************************************
  real(rkp)               ,intent(in)  :: a
  real(rkp)               ,intent(in)  :: ab
  real(rkp)               ,intent(in)  :: ra
  real(rkp)               ,intent(in)  :: rb
  real(rkp) ,dimension(430),intent(in) :: C
  real(rkp) ,dimension(6) ,intent(in)  :: R
  real(rkp) ,dimension(6) ,intent(out) :: dVdR
  
  integer                              :: i, j
  real(rkp)                            :: dist
  real(rkp)                            :: V2
  real(rkp)                            :: dV2dR
  real(rkp) ,dimension(6,6)            :: dmsdr

  ! Initialize dVdR(6)
  do i=1,6
    dVdR(i) = Zero
  enddo


  ! Add dV2dR(i) to dVdR
  ! do i=1,6
  !   dist = R(i)
  !   call ev2gm2(dist, V2, dV2dR, 1, 1)
  !   dVdR(i) = dV2dR
  ! enddo


  ! Calculate dMEG/dr(6,6) for giving R(6)
  call evdmsdr(a, ab, ra, rb, R, dmsdr)

  ! Calculate the monomials for each point by using six MEG terms
  call evdmdr(dmsdr)

  ! Calculate the polynomials by using monomials
  call evdpdr 

  ! Remove 2-body interactions and unconnected terms from polynomials
  call evdbdr

  ! Evaluate dVdR(6) by taken the product of C(j) and dPdR(i,j)
  do i=1,6      
    do j=1,430
     dVdR(i) = dVdR(i) + C(j)*dBdR(i,j)
    enddo
  enddo

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine EvMorse(a, ab, ra, rb, R, rms)
!**********************************************************************
! mixed exponential gaussian term rms = exp(-(r-re)/a-(r-re)^2/ab)
! re: equlibrium bond length
! a:  nonlinear parameter, unit Ang
! ab: nonlinear parameter, unit Ang^2   
!**********************************************************************
  
  real(rkp)               ,intent(in)  :: a
  real(rkp)               ,intent(in)  :: ab
  real(rkp)               ,intent(in)  :: ra
  real(rkp)               ,intent(in)  :: rb
  real(rkp) ,dimension(6) ,intent(in)  :: R
  real(rkp) ,dimension(6) ,intent(out) :: rms

  integer                               :: i

  do i=1,6
    rms(i) = exp( - (R(i) - ra)/a - ((R(i)-rb)**Two)/ab )
  enddo

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine EvdMsdR(a, ab, ra, rb, R, dmsdr)
!**********************************************************************
! Subroutine to evaluate the derivatives of MEG term X
! w.r.t. interatomic distance R(6)
! dmsdR:  Local variables, dirm(6,6)
! a:    Nonlinear parameter(Angstrom)
! ab:       Nonlinear parameter(Angstrom^2)
! re:   equilibrium bond length(Angstrom)
!**********************************************************************
  
  real(rkp)                 ,intent(in)   :: a
  real(rkp)                 ,intent(in)   :: ab
  real(rkp)                 ,intent(in)   :: ra
  real(rkp)                 ,intent(in)   :: rb
  real(rkp) ,dimension(6)   ,intent(in)   :: R
  real(rkp) ,dimension(6,6) ,intent(out)  :: dmsdr

  integer :: i,j

  ! Initialize dmsdr
  dmsdr = Zero

  ! MEG term dmsdr = exp(-(r-re)/a-(r-re)^2/ab)
  ! dmsdr(i,j)=0  i!=j

  do i=1,6
    dmsdr(i,i) = ( - Two*(R(i) - rb)/ab - One/a) * exp(-(r(i)-ra)/a - ((r(i)-rb)**Two)/ab)
  end do 


End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine EvMono(rms)
!**********************************************************************
!  The subroutine reads six MEG variables(X) and calculates the
!  monomials(M) that do not have usable decomposition.
!  For A4 with max. degree 9, the number of monomials is nom.
!**********************************************************************

  real(rkp) ,dimension(6) ,intent(in) :: rms

  rm(0) = One
  rm(1) = rms(6)
  rm(2) = rms(5)
  rm(3) = rms(4)
  rm(4) = rms(3)
  rm(5) = rms(2)
  rm(6) = rms(1)
  rm(7) = rm(3)*rm(4)
  rm(8) = rm(2)*rm(5)
  rm(9) = rm(1)*rm(6)
  rm(10) = rm(1)*rm(2)
  rm(11) = rm(1)*rm(3)
  rm(12) = rm(2)*rm(3)
  rm(13) = rm(1)*rm(4)
  rm(14) = rm(2)*rm(4)
  rm(15) = rm(1)*rm(5)
  rm(16) = rm(3)*rm(5)
  rm(17) = rm(4)*rm(5)
  rm(18) = rm(2)*rm(6)
  rm(19) = rm(3)*rm(6)
  rm(20) = rm(4)*rm(6)
  rm(21) = rm(5)*rm(6)
  rm(22) = rm(1)*rm(7)
  rm(23) = rm(2)*rm(7)
  rm(24) = rm(1)*rm(8)
  rm(25) = rm(2)*rm(16)
  rm(26) = rm(2)*rm(17)
  rm(27) = rm(3)*rm(17)
  rm(28) = rm(1)*rm(18)
  rm(29) = rm(1)*rm(19)
  rm(30) = rm(1)*rm(20)
  rm(31) = rm(3)*rm(20)
  rm(32) = rm(1)*rm(21)
  rm(33) = rm(2)*rm(21)
  rm(34) = rm(1)*rm(12)
  rm(35) = rm(1)*rm(17)
  rm(36) = rm(2)*rm(20)
  rm(37) = rm(3)*rm(21)
  rm(38) = rm(1)*rm(14)
  rm(39) = rm(1)*rm(16)
  rm(40) = rm(2)*rm(19)
  rm(41) = rm(4)*rm(21)
  rm(42) = rm(2)*rm(27)
  rm(43) = rm(1)*rm(31)
  rm(44) = rm(1)*rm(33)
  rm(45) = rm(1)*rm(23)
  rm(46) = rm(1)*rm(25)
  rm(47) = rm(1)*rm(26)
  rm(48) = rm(1)*rm(27)
  rm(49) = rm(1)*rm(40)
  rm(50) = rm(1)*rm(36)
  rm(51) = rm(2)*rm(31)
  rm(52) = rm(1)*rm(37)
  rm(53) = rm(2)*rm(37)
  rm(54) = rm(1)*rm(41)
  rm(55) = rm(2)*rm(41)
  rm(56) = rm(3)*rm(41)
  rm(57) = rm(1)*rm(42)
  rm(58) = rm(1)*rm(51)
  rm(59) = rm(1)*rm(53)
  rm(60) = rm(1)*rm(55)
  rm(61) = rm(1)*rm(56)
  rm(62) = rm(2)*rm(56)
  rm(63) = rm(1)*rm(62)
  rm(64) = rm(2)*rm(57)
  rm(65) = rm(3)*rm(57)
  rm(66) = rm(4)*rm(57)
  rm(67) = rm(5)*rm(57)
  rm(68) = rm(1)*rm(58)
  rm(69) = rm(3)*rm(58)
  rm(70) = rm(4)*rm(58)
  rm(71) = rm(1)*rm(59)
  rm(72) = rm(2)*rm(59)
  rm(73) = rm(1)*rm(60)
  rm(74) = rm(2)*rm(60)
  rm(75) = rm(1)*rm(61)
  rm(76) = rm(2)*rm(62)
  rm(77) = rm(3)*rm(61)
  rm(78) = rm(3)*rm(62)
  rm(79) = rm(4)*rm(61)
  rm(80) = rm(4)*rm(62)
  rm(81) = rm(5)*rm(59)
  rm(82) = rm(5)*rm(60)
  rm(83) = rm(5)*rm(62)
  rm(84) = rm(6)*rm(58)
  rm(85) = rm(6)*rm(59)
  rm(86) = rm(6)*rm(60)
  rm(87) = rm(6)*rm(61)
  rm(88) = rm(2)*rm(64)
  rm(89) = rm(3)*rm(65)
  rm(90) = rm(4)*rm(66)
  rm(91) = rm(5)*rm(67)
  rm(92) = rm(1)*rm(68)
  rm(93) = rm(3)*rm(69)
  rm(94) = rm(4)*rm(70)
  rm(95) = rm(1)*rm(71)
  rm(96) = rm(2)*rm(72)
  rm(97) = rm(1)*rm(73)
  rm(98) = rm(2)*rm(74)
  rm(99) = rm(1)*rm(75)
  rm(100) = rm(2)*rm(76)
  rm(101) = rm(3)*rm(77)
  rm(102) = rm(3)*rm(78)
  rm(103) = rm(4)*rm(79)
  rm(104) = rm(4)*rm(80)
  rm(105) = rm(5)*rm(81)
  rm(106) = rm(5)*rm(82)
  rm(107) = rm(5)*rm(83)
  rm(108) = rm(6)*rm(84)
  rm(109) = rm(6)*rm(85)
  rm(110) = rm(6)*rm(86)
  rm(111) = rm(6)*rm(87)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine EvdMdR(dmsdr)
!**********************************************************************
!  The subroutine reads M(nom) and dMSdR(6,6) and calculates the
!  dMdR(6,nom) that do not have usable decomposition.
!  For A4 with max. degree 9, the number of monomials is nom.
!**********************************************************************

  real(rkp) ,dimension(6,6) ,intent(in) :: dmsdr

  integer :: i

  do i=1,6
    dmdr(i,0) = Zero
    dmdr(i,1) = dmsdr(i,6)
    dmdr(i,2) = dmsdr(i,5)
    dmdr(i,3) = dmsdr(i,4)
    dmdr(i,4) = dmsdr(i,3)
    dmdr(i,5) = dmsdr(i,2)
    dmdr(i,6) = dmsdr(i,1)
    dmdr(i,7) = dmdr(i,3)*rm(4) + rm(3)*dmdr(i,4)
    dmdr(i,8) = dmdr(i,2)*rm(5) + rm(2)*dmdr(i,5)
    dmdr(i,9) = dmdr(i,1)*rm(6) + rm(1)*dmdr(i,6)
    dmdr(i,10) = dmdr(i,1)*rm(2) + rm(1)*dmdr(i,2)
    dmdr(i,11) = dmdr(i,1)*rm(3) + rm(1)*dmdr(i,3)
    dmdr(i,12) = dmdr(i,2)*rm(3) + rm(2)*dmdr(i,3)
    dmdr(i,13) = dmdr(i,1)*rm(4) + rm(1)*dmdr(i,4)
    dmdr(i,14) = dmdr(i,2)*rm(4) + rm(2)*dmdr(i,4)
    dmdr(i,15) = dmdr(i,1)*rm(5) + rm(1)*dmdr(i,5)
    dmdr(i,16) = dmdr(i,3)*rm(5) + rm(3)*dmdr(i,5)
    dmdr(i,17) = dmdr(i,4)*rm(5) + rm(4)*dmdr(i,5)
    dmdr(i,18) = dmdr(i,2)*rm(6) + rm(2)*dmdr(i,6)
    dmdr(i,19) = dmdr(i,3)*rm(6) + rm(3)*dmdr(i,6)
    dmdr(i,20) = dmdr(i,4)*rm(6) + rm(4)*dmdr(i,6)
    dmdr(i,21) = dmdr(i,5)*rm(6) + rm(5)*dmdr(i,6)
    dmdr(i,22) = dmdr(i,1)*rm(7) + rm(1)*dmdr(i,7)
    dmdr(i,23) = dmdr(i,2)*rm(7) + rm(2)*dmdr(i,7)
    dmdr(i,24) = dmdr(i,1)*rm(8) + rm(1)*dmdr(i,8)
    dmdr(i,25) = dmdr(i,2)*rm(16) + rm(2)*dmdr(i,16)
    dmdr(i,26) = dmdr(i,2)*rm(17) + rm(2)*dmdr(i,17)
    dmdr(i,27) = dmdr(i,3)*rm(17) + rm(3)*dmdr(i,17)
    dmdr(i,28) = dmdr(i,1)*rm(18) + rm(1)*dmdr(i,18)
    dmdr(i,29) = dmdr(i,1)*rm(19) + rm(1)*dmdr(i,19)
    dmdr(i,30) = dmdr(i,1)*rm(20) + rm(1)*dmdr(i,20)
    dmdr(i,31) = dmdr(i,3)*rm(20) + rm(3)*dmdr(i,20)
    dmdr(i,32) = dmdr(i,1)*rm(21) + rm(1)*dmdr(i,21)
    dmdr(i,33) = dmdr(i,2)*rm(21) + rm(2)*dmdr(i,21)
    dmdr(i,34) = dmdr(i,1)*rm(12) + rm(1)*dmdr(i,12)
    dmdr(i,35) = dmdr(i,1)*rm(17) + rm(1)*dmdr(i,17)
    dmdr(i,36) = dmdr(i,2)*rm(20) + rm(2)*dmdr(i,20)
    dmdr(i,37) = dmdr(i,3)*rm(21) + rm(3)*dmdr(i,21)
    dmdr(i,38) = dmdr(i,1)*rm(14) + rm(1)*dmdr(i,14)
    dmdr(i,39) = dmdr(i,1)*rm(16) + rm(1)*dmdr(i,16)
    dmdr(i,40) = dmdr(i,2)*rm(19) + rm(2)*dmdr(i,19)
    dmdr(i,41) = dmdr(i,4)*rm(21) + rm(4)*dmdr(i,21)
    dmdr(i,42) = dmdr(i,2)*rm(27) + rm(2)*dmdr(i,27)
    dmdr(i,43) = dmdr(i,1)*rm(31) + rm(1)*dmdr(i,31)
    dmdr(i,44) = dmdr(i,1)*rm(33) + rm(1)*dmdr(i,33)
    dmdr(i,45) = dmdr(i,1)*rm(23) + rm(1)*dmdr(i,23)
    dmdr(i,46) = dmdr(i,1)*rm(25) + rm(1)*dmdr(i,25)
    dmdr(i,47) = dmdr(i,1)*rm(26) + rm(1)*dmdr(i,26)
    dmdr(i,48) = dmdr(i,1)*rm(27) + rm(1)*dmdr(i,27)
    dmdr(i,49) = dmdr(i,1)*rm(40) + rm(1)*dmdr(i,40)
    dmdr(i,50) = dmdr(i,1)*rm(36) + rm(1)*dmdr(i,36)
    dmdr(i,51) = dmdr(i,2)*rm(31) + rm(2)*dmdr(i,31)
    dmdr(i,52) = dmdr(i,1)*rm(37) + rm(1)*dmdr(i,37)
    dmdr(i,53) = dmdr(i,2)*rm(37) + rm(2)*dmdr(i,37)
    dmdr(i,54) = dmdr(i,1)*rm(41) + rm(1)*dmdr(i,41)
    dmdr(i,55) = dmdr(i,2)*rm(41) + rm(2)*dmdr(i,41)
    dmdr(i,56) = dmdr(i,3)*rm(41) + rm(3)*dmdr(i,41)
    dmdr(i,57) = dmdr(i,1)*rm(42) + rm(1)*dmdr(i,42)
    dmdr(i,58) = dmdr(i,1)*rm(51) + rm(1)*dmdr(i,51)
    dmdr(i,59) = dmdr(i,1)*rm(53) + rm(1)*dmdr(i,53)
    dmdr(i,60) = dmdr(i,1)*rm(55) + rm(1)*dmdr(i,55)
    dmdr(i,61) = dmdr(i,1)*rm(56) + rm(1)*dmdr(i,56)
    dmdr(i,62) = dmdr(i,2)*rm(56) + rm(2)*dmdr(i,56)
    dmdr(i,63) = dmdr(i,1)*rm(62) + rm(1)*dmdr(i,62)
    dmdr(i,64) = dmdr(i,2)*rm(57) + rm(2)*dmdr(i,57)
    dmdr(i,65) = dmdr(i,3)*rm(57) + rm(3)*dmdr(i,57)
    dmdr(i,66) = dmdr(i,4)*rm(57) + rm(4)*dmdr(i,57)
    dmdr(i,67) = dmdr(i,5)*rm(57) + rm(5)*dmdr(i,57)
    dmdr(i,68) = dmdr(i,1)*rm(58) + rm(1)*dmdr(i,58)
    dmdr(i,69) = dmdr(i,3)*rm(58) + rm(3)*dmdr(i,58)
    dmdr(i,70) = dmdr(i,4)*rm(58) + rm(4)*dmdr(i,58)
    dmdr(i,71) = dmdr(i,1)*rm(59) + rm(1)*dmdr(i,59)
    dmdr(i,72) = dmdr(i,2)*rm(59) + rm(2)*dmdr(i,59)
    dmdr(i,73) = dmdr(i,1)*rm(60) + rm(1)*dmdr(i,60)
    dmdr(i,74) = dmdr(i,2)*rm(60) + rm(2)*dmdr(i,60)
    dmdr(i,75) = dmdr(i,1)*rm(61) + rm(1)*dmdr(i,61)
    dmdr(i,76) = dmdr(i,2)*rm(62) + rm(2)*dmdr(i,62)
    dmdr(i,77) = dmdr(i,3)*rm(61) + rm(3)*dmdr(i,61)
    dmdr(i,78) = dmdr(i,3)*rm(62) + rm(3)*dmdr(i,62)
    dmdr(i,79) = dmdr(i,4)*rm(61) + rm(4)*dmdr(i,61)
    dmdr(i,80) = dmdr(i,4)*rm(62) + rm(4)*dmdr(i,62)
    dmdr(i,81) = dmdr(i,5)*rm(59) + rm(5)*dmdr(i,59)
    dmdr(i,82) = dmdr(i,5)*rm(60) + rm(5)*dmdr(i,60)
    dmdr(i,83) = dmdr(i,5)*rm(62) + rm(5)*dmdr(i,62)
    dmdr(i,84) = dmdr(i,6)*rm(58) + rm(6)*dmdr(i,58)
    dmdr(i,85) = dmdr(i,6)*rm(59) + rm(6)*dmdr(i,59)
    dmdr(i,86) = dmdr(i,6)*rm(60) + rm(6)*dmdr(i,60)
    dmdr(i,87) = dmdr(i,6)*rm(61) + rm(6)*dmdr(i,61)
    dmdr(i,88) = dmdr(i,2)*rm(64) + rm(2)*dmdr(i,64)
    dmdr(i,89) = dmdr(i,3)*rm(65) + rm(3)*dmdr(i,65)
    dmdr(i,90) = dmdr(i,4)*rm(66) + rm(4)*dmdr(i,66)
    dmdr(i,91) = dmdr(i,5)*rm(67) + rm(5)*dmdr(i,67)
    dmdr(i,92) = dmdr(i,1)*rm(68) + rm(1)*dmdr(i,68)
    dmdr(i,93) = dmdr(i,3)*rm(69) + rm(3)*dmdr(i,69)
    dmdr(i,94) = dmdr(i,4)*rm(70) + rm(4)*dmdr(i,70)
    dmdr(i,95) = dmdr(i,1)*rm(71) + rm(1)*dmdr(i,71)
    dmdr(i,96) = dmdr(i,2)*rm(72) + rm(2)*dmdr(i,72)
    dmdr(i,97) = dmdr(i,1)*rm(73) + rm(1)*dmdr(i,73)
    dmdr(i,98) = dmdr(i,2)*rm(74) + rm(2)*dmdr(i,74)
    dmdr(i,99) = dmdr(i,1)*rm(75) + rm(1)*dmdr(i,75)
    dmdr(i,100) = dmdr(i,2)*rm(76) + rm(2)*dmdr(i,76)
    dmdr(i,101) = dmdr(i,3)*rm(77) + rm(3)*dmdr(i,77)
    dmdr(i,102) = dmdr(i,3)*rm(78) + rm(3)*dmdr(i,78)
    dmdr(i,103) = dmdr(i,4)*rm(79) + rm(4)*dmdr(i,79)
    dmdr(i,104) = dmdr(i,4)*rm(80) + rm(4)*dmdr(i,80)
    dmdr(i,105) = dmdr(i,5)*rm(81) + rm(5)*dmdr(i,81)
    dmdr(i,106) = dmdr(i,5)*rm(82) + rm(5)*dmdr(i,82)
    dmdr(i,107) = dmdr(i,5)*rm(83) + rm(5)*dmdr(i,83)
    dmdr(i,108) = dmdr(i,6)*rm(84) + rm(6)*dmdr(i,84)
    dmdr(i,109) = dmdr(i,6)*rm(85) + rm(6)*dmdr(i,85)
    dmdr(i,110) = dmdr(i,6)*rm(86) + rm(6)*dmdr(i,86)
    dmdr(i,111) = dmdr(i,6)*rm(87) + rm(6)*dmdr(i,87)
  enddo

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine EvPoly
!**********************************************************************
!  The subroutine reads monomials(m) and calculates the
!  permutation-invariant polynomials(p)
!  For A4 with max. degree 9, the number of polynomials is nob.
!**********************************************************************

  p(0) = rm(0)
  p(1) = rm(1) + rm(2) + rm(3) + rm(4) + rm(5) + rm(6)
  p(2) = rm(7) + rm(8) + rm(9)
  p(3) = rm(10) + rm(11) + rm(12) + rm(13) + rm(14) + rm(15) + rm(16) + rm(17) + rm(18) + rm(19) + rm(20) + rm(21)
  p(4) = p(1)*p(1) - p(3) - p(2) - p(3) - p(2)
  p(5) = rm(22) + rm(23) + rm(24) + rm(25) + rm(26) + rm(27) + rm(28) + rm(29) + rm(30) + rm(31) + rm(32) + rm(33)
  p(6) = rm(34) + rm(35) + rm(36) + rm(37)
  p(7) = rm(38) + rm(39) + rm(40) + rm(41)
  p(8) = p(1)*p(2) - p(5)
  p(9) = p(1)*p(3) - p(6) - p(7) - p(5) - p(6) - p(7) - p(5) - p(6) - p(7)
  p(10) = p(1)*p(4) - p(9) - p(8)
  p(11) = rm(42) + rm(43) + rm(44)
  p(12) = rm(45) + rm(46) + rm(47) + rm(48) + rm(49) + rm(50) + rm(51) + rm(52) + rm(53) + rm(54) + rm(55) + rm(56)
  p(13) = p(2)*p(3) - p(12)
  p(14) = p(1)*p(5) - p(12) - p(11) - p(13) - p(12) - p(11) - p(11) - p(11)
  p(15) = p(1)*p(6) - p(12)
  p(16) = p(1)*p(7) - p(12)
  p(17) = p(2)*p(2) - p(11) - p(11)
  p(18) = p(3)*p(3) - p(12) - p(11) - p(15) - p(16) - p(14) - p(12) - p(11) - p(15) - p(16) - p(14) - p(12) - p(11) - p(12) - p(11)
  p(19) = p(2)*p(4) - p(14)
  p(20) = p(3)*p(4) - p(15) - p(16) - p(13)
  p(21) = p(1)*p(10) - p(20) - p(19)
  p(22) = rm(57) + rm(58) + rm(59) + rm(60) + rm(61) + rm(62)
  p(23) = p(1)*p(11) - p(22)
  p(24) = p(2)*p(6)
  p(25) = p(2)*p(7)
  p(26) = p(1)*p(12) - p(22) - p(24) - p(25) - p(22) - p(22) - p(22)
  p(27) = p(2)*p(5) - p(22) - p(23) - p(22)
  p(28) = p(3)*p(5) - p(22) - p(26) - p(24) - p(25) - p(23) - p(22) - p(24) - p(25) - p(23) - p(22) - p(22)
  p(29) = p(3)*p(6) - p(22) - p(26) - p(22)
  p(30) = p(3)*p(7) - p(22) - p(26) - p(22)
  p(31) = p(2)*p(9) - p(26) - p(28)
  p(32) = p(1)*p(14) - p(26) - p(23) - p(28)
  p(33) = p(4)*p(6) - p(25)
  p(34) = p(4)*p(7) - p(24)
  p(35) = p(1)*p(17) - p(27)
  p(36) = p(1)*p(18) - p(29) - p(30) - p(28)
  p(37) = p(2)*p(10) - p(32)
  p(38) = p(3)*p(10) - p(33) - p(34) - p(31)
  p(39) = p(1)*p(21) - p(38) - p(37)
  p(40) = rm(63)
  p(41) = rm(64) + rm(65) + rm(66) + rm(67) + rm(68) + rm(69) + rm(70) + rm(71) + rm(72) + rm(73) + rm(74) + rm(75) + rm(76) + rm(77) + rm(78) + rm(79) + rm(80) + rm(81) + rm(82) + rm(83) + rm(84) + rm(85) + rm(86) + rm(87)
  p(42) = p(1)*p(22) - p(40) - p(41) - p(40) - p(40) - p(40) - p(40) - p(40)
  p(43) = p(2)*p(11) - p(40) - p(40) - p(40)
  p(44) = p(2)*p(12) - p(41)
  p(45) = p(3)*p(11) - p(41)
  p(46) = p(5)*p(6) - p(41)
  p(47) = p(5)*p(7) - p(41)
  p(48) = p(6)*p(7) - p(40) - p(40) - p(40) - p(40)
  p(49) = p(4)*p(11) - p(42)
  p(50) = p(2)*p(15) - p(46)
  p(51) = p(2)*p(16) - p(47)
  p(52) = p(4)*p(12) - p(41) - p(50) - p(51)
  p(53) = p(2)*p(14) - p(42) - p(49) - p(42)
  p(54) = p(6)*p(6) - p(42) - p(42)
  p(55) = p(7)*p(7) - p(42) - p(42)
  p(56) = p(3)*p(17) - p(44)
  p(57) = p(2)*p(18) - p(48)
  p(58) = p(3)*p(14) - p(41) - p(52) - p(46) - p(47) - p(45) - p(45)
  p(59) = p(6)*p(9) - p(41) - p(52) - p(47)
  p(60) = p(7)*p(9) - p(41) - p(52) - p(46)
  p(61) = p(2)*p(20) - p(52) - p(58)
  p(62) = p(1)*p(32) - p(52) - p(49) - p(58)
  p(63) = p(6)*p(10) - p(51)
  p(64) = p(7)*p(10) - p(50)
  p(65) = p(2)*p(17) - p(43)
  p(66) = p(3)*p(18) - p(46) - p(47) - p(45) - p(59) - p(60) - p(58)
  p(67) = p(2)*p(19) - p(49)
  p(68) = p(1)*p(36) - p(59) - p(60) - p(58) - p(57) - p(66) - p(66)
  p(69) = p(2)*p(21) - p(62)
  p(70) = p(3)*p(21) - p(63) - p(64) - p(61)
  p(71) = p(1)*p(39) - p(70) - p(69)
  p(72) = p(40)*p(1)
  p(73) = p(2)*p(22) - p(72)
  p(74) = p(6)*p(11)
  p(75) = p(7)*p(11)
  p(76) = p(3)*p(22) - p(72) - p(74) - p(75) - p(72) - p(72) - p(72)
  p(77) = rm(88) + rm(89) + rm(90) + rm(91) + rm(92) + rm(93) + rm(94) + rm(95) + rm(96) + rm(97) + rm(98) + rm(99) + rm(100) + rm(101) + rm(102) + rm(103) + rm(104) + rm(105) + rm(106) + rm(107) + rm(108) + rm(109) + rm(110) + rm(111)
  p(78) = p(1)*p(42) - p(72) - p(76)
  p(79) = p(5)*p(11) - p(72) - p(73) - p(72)
  p(80) = p(2)*p(26) - p(76) - p(77)
  p(81) = p(6)*p(12) - p(72) - p(76) - p(72)
  p(82) = p(7)*p(12) - p(72) - p(76) - p(72)
  p(83) = p(8)*p(11) - p(72)
  p(84) = p(6)*p(17)
  p(85) = p(7)*p(17)
  p(86) = p(9)*p(11) - p(76) - p(77)
  p(87) = p(2)*p(29) - p(81)
  p(88) = p(6)*p(14) - p(75) - p(75)
  p(89) = p(2)*p(30) - p(82)
  p(90) = p(7)*p(14) - p(74) - p(74)
  p(91) = p(1)*p(48) - p(76) - p(81) - p(82)
  p(92) = p(10)*p(11) - p(78)
  p(93) = p(2)*p(33) - p(88)
  p(94) = p(2)*p(34) - p(90)
  p(95) = p(10)*p(12) - p(77) - p(93) - p(94)
  p(96) = p(2)*p(27) - p(73) - p(79)
  p(97) = p(2)*p(28) - p(76) - p(86)
  p(98) = p(1)*p(53) - p(80) - p(79) - p(97)
  p(99) = p(1)*p(54) - p(81)
  p(100) = p(1)*p(55) - p(82)
  p(101) = p(5)*p(18) - p(76) - p(91) - p(87) - p(89) - p(86)
  p(102) = p(6)*p(18) - p(74) - p(90)
  p(103) = p(7)*p(18) - p(75) - p(88)
  p(104) = p(2)*p(31) - p(77) - p(86)
  p(105) = p(2)*p(36) - p(91) - p(101)
  p(106) = p(3)*p(32) - p(77) - p(95) - p(88) - p(90) - p(86)
  p(107) = p(4)*p(29) - p(82) - p(80) - p(99)
  p(108) = p(4)*p(30) - p(81) - p(80) - p(100)
  p(109) = p(2)*p(38) - p(95) - p(106)
  p(110) = p(1)*p(62) - p(95) - p(92) - p(106)
  p(111) = p(6)*p(21) - p(94)
  p(112) = p(7)*p(21) - p(93)
  p(113) = p(1)*p(65) - p(96)
  p(114) = p(1)*p(66) - p(102) - p(103) - p(101)
  p(115) = p(2)*p(37) - p(92)
  p(116) = p(10)*p(18) - p(99) - p(100) - p(97)
  p(117) = p(2)*p(39) - p(110)
  p(118) = p(3)*p(39) - p(111) - p(112) - p(109)
  p(119) = p(1)*p(71) - p(118) - p(117)
  p(120) = p(40)*p(2)
  p(121) = p(40)*p(3)
  p(122) = p(40)*p(4)
  p(123) = p(11)*p(12) - p(121)
  p(124) = p(2)*p(42) - p(122)
  p(125) = p(6)*p(22) - p(121)
  p(126) = p(7)*p(22) - p(121)
  p(127) = p(2)*p(41) - p(121) - p(123) - p(121)
  p(128) = p(6)*p(23) - p(123)
  p(129) = p(7)*p(23) - p(123)
  p(130) = p(3)*p(41) - p(122) - p(121) - p(120) - p(125) - p(128) - p(126) - p(129) - p(124) - p(123) - p(122) - p(121) - p(120) - p(125) - p(126) - p(124) - p(123) - p(122) - p(121) - p(120) - p(122) - p(121) - p(120) - p(120) - p(120) - p(120) - p(120)
  p(131) = p(3)*p(42) - p(121) - p(125) - p(126) - p(121)
  p(132) = p(4)*p(41) - p(121) - p(131) - p(128) - p(129) - p(127) - p(121)
  p(133) = p(1)*p(78) - p(122) - p(131)
  p(134) = p(11)*p(11) - p(120) - p(120)
  p(135) = p(2)*p(48) - p(130)
  p(136) = p(11)*p(17) - p(120)
  p(137) = p(2)*p(44) - p(123)
  p(138) = p(2)*p(45) - p(121)
  p(139) = p(11)*p(14) - p(122) - p(124) - p(122)
  p(140) = p(6)*p(27) - p(127)
  p(141) = p(2)*p(54)
  p(142) = p(7)*p(27) - p(127)
  p(143) = p(2)*p(52) - p(131) - p(132)
  p(144) = p(1)*p(81) - p(125) - p(141) - p(135) - p(125)
  p(145) = p(2)*p(55)
  p(146) = p(1)*p(82) - p(126) - p(135) - p(145) - p(126)
  p(147) = p(11)*p(18) - p(130)
  p(148) = p(6)*p(28) - p(129) - p(123) - p(143)
  p(149) = p(7)*p(28) - p(128) - p(123) - p(143)
  p(150) = p(6)*p(30) - p(121) - p(146)
  p(151) = p(11)*p(19) - p(122)
  p(152) = p(2)*p(50) - p(128)
  p(153) = p(2)*p(51) - p(129)
  p(154) = p(11)*p(20) - p(131) - p(132)
  p(155) = p(2)*p(59) - p(144) - p(148)
  p(156) = p(6)*p(32) - p(129)
  p(157) = p(2)*p(60) - p(146) - p(149)
  p(158) = p(7)*p(32) - p(128)
  p(159) = p(6)*p(34) - p(122) - p(145) - p(122)
  p(160) = p(11)*p(21) - p(133)
  p(161) = p(2)*p(63) - p(156)
  p(162) = p(2)*p(64) - p(158)
  p(163) = p(12)*p(21) - p(132) - p(161) - p(162)
  p(164) = p(2)*p(53) - p(124) - p(139)
  p(165) = p(2)*p(58) - p(131) - p(154)
  p(166) = p(3)*p(54) - p(125) - p(144)
  p(167) = p(3)*p(55) - p(126) - p(146)
  p(168) = p(3)*p(65) - p(137)
  p(169) = p(17)*p(18) - p(135)
  p(170) = p(1)*p(98) - p(143) - p(139) - p(165)
  p(171) = p(4)*p(54) - p(135)
  p(172) = p(4)*p(55) - p(135)
  p(173) = p(2)*p(66) - p(150)
  p(174) = p(1)*p(101) - p(148) - p(149) - p(147) - p(173) - p(165) - p(147)
  p(175) = p(1)*p(102) - p(150) - p(148) - p(166)
  p(176) = p(1)*p(103) - p(150) - p(149) - p(167)
  p(177) = p(2)*p(61) - p(132) - p(154)
  p(178) = p(2)*p(68) - p(159) - p(174)
  p(179) = p(3)*p(62) - p(132) - p(163) - p(156) - p(158) - p(154)
  p(180) = p(6)*p(38) - p(132) - p(163) - p(157)
  p(181) = p(7)*p(38) - p(132) - p(163) - p(155)
  p(182) = p(2)*p(70) - p(163) - p(179)
  p(183) = p(1)*p(110) - p(163) - p(160) - p(179)
  p(184) = p(6)*p(39) - p(162)
  p(185) = p(7)*p(39) - p(161)
  p(186) = p(2)*p(65) - p(136)
  p(187) = p(3)*p(66) - p(148) - p(149) - p(147) - p(175) - p(176) - p(174)
  p(188) = p(2)*p(67) - p(151)
  p(189) = p(4)*p(66) - p(166) - p(167) - p(165)
  p(190) = p(2)*p(69) - p(160)
  p(191) = p(18)*p(21) - p(171) - p(172) - p(169)
  p(192) = p(2)*p(71) - p(183)
  p(193) = p(3)*p(71) - p(184) - p(185) - p(182)
  p(194) = p(1)*p(119) - p(193) - p(192)
  p(195) = p(40)*p(5)
  p(196) = p(40)*p(6)
  p(197) = p(40)*p(7)
  p(198) = p(40)*p(8)
  p(199) = p(40)*p(9)
  p(200) = p(40)*p(10)
  p(201) = p(11)*p(22) - p(195)
  p(202) = p(12)*p(22) - p(196) - p(197) - p(195) - p(196) - p(197) - p(195) - p(196) - p(197)
  p(203) = p(17)*p(22) - p(198)
  p(204) = p(6)*p(43)
  p(205) = p(7)*p(43)
  p(206) = p(11)*p(26) - p(199) - p(202)
  p(207) = p(2)*p(76) - p(199) - p(202)
  p(208) = p(2)*p(78) - p(200)
  p(209) = p(6)*p(41) - p(199) - p(195) - p(202) - p(195)
  p(210) = p(6)*p(42) - p(197) - p(197) - p(197)
  p(211) = p(7)*p(41) - p(199) - p(195) - p(202) - p(195)
  p(212) = p(7)*p(42) - p(196) - p(196) - p(196)
  p(213) = p(11)*p(29) - p(209)
  p(214) = p(11)*p(30) - p(211)
  p(215) = p(18)*p(22) - p(199) - p(213) - p(214)
  p(216) = p(2)*p(77) - p(199) - p(206)
  p(217) = p(6)*p(49) - p(205)
  p(218) = p(7)*p(49) - p(204)
  p(219) = p(3)*p(77) - p(200) - p(199) - p(198) - p(209) - p(217) - p(211) - p(218) - p(207) - p(204) - p(205) - p(200) - p(199) - p(198) - p(200) - p(198) - p(200) - p(198)
  p(220) = p(3)*p(78) - p(199) - p(210) - p(212)
  p(221) = p(10)*p(41) - p(199) - p(220) - p(217) - p(218) - p(216)
  p(222) = p(1)*p(133) - p(200) - p(220)
  p(223) = p(11)*p(27) - p(195) - p(203)
  p(224) = p(1)*p(134) - p(201)
  p(225) = p(2)*p(81) - p(209)
  p(226) = p(2)*p(80) - p(202) - p(206)
  p(227) = p(2)*p(82) - p(211)
  p(228) = p(2)*p(91) - p(215) - p(219)
  p(229) = p(11)*p(28) - p(199) - p(207)
  p(230) = p(6)*p(53) - p(205)
  p(231) = p(5)*p(54) - p(209)
  p(232) = p(7)*p(53) - p(204)
  p(233) = p(7)*p(54) - p(196)
  p(234) = p(5)*p(55) - p(211)
  p(235) = p(6)*p(55) - p(197)
  p(236) = p(11)*p(35) - p(198)
  p(237) = p(6)*p(65)
  p(238) = p(7)*p(65)
  p(239) = p(2)*p(86) - p(199) - p(229)
  p(240) = p(1)*p(139) - p(206) - p(229) - p(224)
  p(241) = p(17)*p(29) - p(225)
  p(242) = p(2)*p(99) - p(231)
  p(243) = p(17)*p(30) - p(227)
  p(244) = p(2)*p(95) - p(220) - p(221)
  p(245) = p(4)*p(81) - p(202) - p(242) - p(227)
  p(246) = p(2)*p(100) - p(234)
  p(247) = p(4)*p(82) - p(202) - p(225) - p(246)
  p(248) = p(11)*p(36) - p(215) - p(219)
  p(249) = p(2)*p(102) - p(233)
  p(250) = p(6)*p(58) - p(206) - p(214) - p(244) - p(214)
  p(251) = p(2)*p(103) - p(235)
  p(252) = p(7)*p(58) - p(213) - p(206) - p(244) - p(213)
  p(253) = p(1)*p(150) - p(215) - p(233) - p(235)
  p(254) = p(11)*p(37) - p(200)
  p(255) = p(2)*p(93) - p(217)
  p(256) = p(2)*p(94) - p(218)
  p(257) = p(11)*p(38) - p(220) - p(221)
  p(258) = p(2)*p(107) - p(245) - p(250)
  p(259) = p(6)*p(62) - p(218)
  p(260) = p(2)*p(108) - p(247) - p(252)
  p(261) = p(7)*p(62) - p(217)
  p(262) = p(6)*p(64) - p(200) - p(246) - p(200)
  p(263) = p(11)*p(39) - p(222)
  p(264) = p(2)*p(111) - p(259)
  p(265) = p(2)*p(112) - p(261)
  p(266) = p(12)*p(39) - p(221) - p(264) - p(265)
  p(267) = p(2)*p(98) - p(208) - p(240)
  p(268) = p(6)*p(54) - p(210)
  p(269) = p(7)*p(55) - p(212)
  p(270) = p(2)*p(96) - p(203) - p(223)
  p(271) = p(2)*p(97) - p(207) - p(229)
  p(272) = p(2)*p(101) - p(215) - p(248)
  p(273) = p(2)*p(106) - p(220) - p(257)
  p(274) = p(6)*p(59) - p(220) - p(215) - p(209)
  p(275) = p(7)*p(60) - p(220) - p(215) - p(211)
  p(276) = p(5)*p(66) - p(215) - p(253) - p(249) - p(251) - p(248)
  p(277) = p(6)*p(66) - p(213) - p(252)
  p(278) = p(7)*p(66) - p(214) - p(250)
  p(279) = p(2)*p(104) - p(216) - p(239)
  p(280) = p(2)*p(105) - p(219) - p(248)
  p(281) = p(1)*p(170) - p(244) - p(240) - p(273)
  p(282) = p(10)*p(54) - p(227)
  p(283) = p(10)*p(55) - p(225)
  p(284) = p(2)*p(114) - p(253) - p(276)
  p(285) = p(1)*p(174) - p(250) - p(252) - p(248) - p(276) - p(273)
  p(286) = p(6)*p(68) - p(217) - p(261) - p(251)
  p(287) = p(7)*p(68) - p(218) - p(259) - p(249)
  p(288) = p(2)*p(109) - p(221) - p(257)
  p(289) = p(2)*p(116) - p(262) - p(285)
  p(290) = p(3)*p(110) - p(221) - p(266) - p(259) - p(261) - p(257)
  p(291) = p(6)*p(70) - p(221) - p(266) - p(260)
  p(292) = p(7)*p(70) - p(221) - p(266) - p(258)
  p(293) = p(2)*p(118) - p(266) - p(290)
  p(294) = p(1)*p(183) - p(266) - p(263) - p(290)
  p(295) = p(6)*p(71) - p(265)
  p(296) = p(7)*p(71) - p(264)
  p(297) = p(1)*p(186) - p(270)
  p(298) = p(1)*p(187) - p(277) - p(278) - p(276)
  p(299) = p(2)*p(115) - p(254)
  p(300) = p(1)*p(189) - p(286) - p(287) - p(285) - p(284) - p(298)
  p(301) = p(2)*p(117) - p(263)
  p(302) = p(18)*p(39) - p(282) - p(283) - p(280)
  p(303) = p(2)*p(119) - p(294)
  p(304) = p(3)*p(119) - p(295) - p(296) - p(293)
  p(305) = p(1)*p(194) - p(304) - p(303)
  p(306) = p(40)*p(11)
  p(307) = p(40)*p(12)
  p(308) = p(40)*p(17)
  p(309) = p(40)*p(13)
  p(310) = p(40)*p(14)
  p(311) = p(40)*p(15)
  p(312) = p(40)*p(16)
  p(313) = p(40)*p(18)
  p(314) = p(40)*p(19)
  p(315) = p(40)*p(20)
  p(316) = p(40)*p(21)
  p(317) = p(11)*p(42) - p(310)
  p(318) = p(11)*p(44) - p(307)
  p(319) = p(12)*p(43) - p(309) - p(318)
  p(320) = p(17)*p(42) - p(314)
  p(321) = p(2)*p(125) - p(311)
  p(322) = p(2)*p(126) - p(312)
  p(323) = p(11)*p(48) - p(313)
  p(324) = p(12)*p(42) - p(311) - p(312) - p(307) - p(307)
  p(325) = p(6)*p(79) - p(319)
  p(326) = p(11)*p(54)
  p(327) = p(7)*p(79) - p(319)
  p(328) = p(2)*p(131) - p(315) - p(324)
  p(329) = p(22)*p(29) - p(313) - p(311) - p(326) - p(311)
  p(330) = p(11)*p(55)
  p(331) = p(22)*p(30) - p(313) - p(312) - p(330) - p(312)
  p(332) = p(2)*p(127) - p(309) - p(319)
  p(333) = p(6)*p(83) - p(318)
  p(334) = p(7)*p(83) - p(318)
  p(335) = p(11)*p(52) - p(315) - p(324)
  p(336) = p(2)*p(130) - p(313) - p(323) - p(313)
  p(337) = p(2)*p(133) - p(316)
  p(338) = p(6)*p(77) - p(315) - p(309) - p(322)
  p(339) = p(6)*p(78) - p(312)
  p(340) = p(7)*p(77) - p(315) - p(309) - p(321)
  p(341) = p(7)*p(78) - p(311)
  p(342) = p(11)*p(59) - p(329) - p(338)
  p(343) = p(11)*p(60) - p(331) - p(340)
  p(344) = p(3)*p(130) - p(315) - p(311) - p(312) - p(309) - p(329) - p(338) - p(331) - p(340) - p(328) - p(321) - p(322) - p(311) - p(312)
  p(345) = p(18)*p(42) - p(310) - p(326) - p(330) - p(310)
  p(346) = p(2)*p(132) - p(315) - p(335)
  p(347) = p(6)*p(92) - p(334)
  p(348) = p(7)*p(92) - p(333)
  p(349) = p(3)*p(132) - p(316) - p(315) - p(314) - p(338) - p(347) - p(340) - p(348) - p(336) - p(333) - p(334) - p(316) - p(315) - p(314) - p(316) - p(314) - p(316) - p(314)
  p(350) = p(3)*p(133) - p(315) - p(339) - p(341)
  p(351) = p(4)*p(132) - p(315) - p(344) - p(342) - p(343) - p(332)
  p(352) = p(1)*p(222) - p(316) - p(350)
  p(353) = p(2)*p(134) - p(306)
  p(354) = p(2)*p(135) - p(323)
  p(355) = p(3)*p(134) - p(319)
  p(356) = p(11)*p(53) - p(310) - p(320)
  p(357) = p(2)*p(144) - p(329) - p(338)
  p(358) = p(2)*p(143) - p(324) - p(335)
  p(359) = p(6)*p(81) - p(311) - p(324)
  p(360) = p(2)*p(146) - p(331) - p(340)
  p(361) = p(2)*p(150) - p(344)
  p(362) = p(7)*p(82) - p(312) - p(324)
  p(363) = p(11)*p(65) - p(308)
  p(364) = p(2)*p(137) - p(318)
  p(365) = p(17)*p(45) - p(307)
  p(366) = p(4)*p(134) - p(317)
  p(367) = p(6)*p(96) - p(332)
  p(368) = p(17)*p(54)
  p(369) = p(7)*p(96) - p(332)
  p(370) = p(17)*p(55)
  p(371) = p(2)*p(159) - p(345) - p(349)
  p(372) = p(2)*p(147) - p(313)
  p(373) = p(11)*p(58) - p(315) - p(328)
  p(374) = p(2)*p(148) - p(329) - p(342)
  p(375) = p(6)*p(98) - p(327)
  p(376) = p(2)*p(166) - p(359)
  p(377) = p(14)*p(54) - p(323)
  p(378) = p(2)*p(149) - p(331) - p(343)
  p(379) = p(7)*p(98) - p(325)
  p(380) = p(7)*p(99) - p(311) - p(359)
  p(381) = p(2)*p(167) - p(362)
  p(382) = p(14)*p(55) - p(323)
  p(383) = p(6)*p(100) - p(312) - p(362)
  p(384) = p(11)*p(66) - p(344)
  p(385) = p(6)*p(101) - p(325) - p(343) - p(379)
  p(386) = p(7)*p(101) - p(342) - p(327) - p(375)
  p(387) = p(6)*p(103) - p(313) - p(382)
  p(388) = p(11)*p(67) - p(314)
  p(389) = p(2)*p(152) - p(333)
  p(390) = p(2)*p(153) - p(334)
  p(391) = p(2)*p(154) - p(315) - p(373)
  p(392) = p(1)*p(240) - p(335) - p(373) - p(366)
  p(393) = p(2)*p(155) - p(338) - p(342)
  p(394) = p(2)*p(171) - p(377)
  p(395) = p(2)*p(157) - p(340) - p(343)
  p(396) = p(2)*p(163) - p(350) - p(351)
  p(397) = p(6)*p(95) - p(315) - p(350) - p(340)
  p(398) = p(2)*p(172) - p(382)
  p(399) = p(7)*p(95) - p(315) - p(350) - p(338)
  p(400) = p(11)*p(68) - p(345) - p(349)
  p(401) = p(2)*p(175) - p(380) - p(385)
  p(402) = p(6)*p(106) - p(335) - p(343) - p(396)
  p(403) = p(2)*p(176) - p(383) - p(386)
  p(404) = p(7)*p(106) - p(342) - p(335) - p(396)
  p(405) = p(4)*p(150) - p(328) - p(359) - p(362)
  p(406) = p(11)*p(69) - p(316)
  p(407) = p(2)*p(161) - p(347)
  p(408) = p(2)*p(162) - p(348)
  p(409) = p(11)*p(70) - p(350) - p(351)
  p(410) = p(2)*p(180) - p(397) - p(402)
  p(411) = p(6)*p(110) - p(348)
  p(412) = p(2)*p(181) - p(399) - p(404)
  p(413) = p(7)*p(110) - p(347)
  p(414) = p(6)*p(112) - p(316) - p(398) - p(316)
  p(415) = p(11)*p(71) - p(352)
  p(416) = p(2)*p(184) - p(411)
  p(417) = p(2)*p(185) - p(413)
  p(418) = p(12)*p(71) - p(351) - p(416) - p(417)
  p(419) = p(2)*p(164) - p(320) - p(356)
  p(420) = p(2)*p(165) - p(328) - p(373)
  p(421) = p(2)*p(170) - p(337) - p(392)
  p(422) = p(1)*p(268) - p(359)
  p(423) = p(1)*p(269) - p(362)
  p(424) = p(2)*p(174) - p(345) - p(400)
  p(425) = p(6)*p(102) - p(345) - p(326)
  p(426) = p(7)*p(103) - p(345) - p(330)
  p(427) = p(3)*p(186) - p(364)
  p(428) = p(18)*p(65) - p(354)
  p(429) = p(17)*p(66) - p(361)
  p(430) = p(2)*p(179) - p(350) - p(409)
  p(431) = p(4)*p(166) - p(361) - p(357) - p(422)
  p(432) = p(4)*p(167) - p(361) - p(360) - p(423)
  p(433) = p(2)*p(187) - p(387)
  p(434) = p(14)*p(66) - p(328) - p(405) - p(376) - p(381) - p(373)
  p(435) = p(1)*p(277) - p(387) - p(385) - p(425)
  p(436) = p(1)*p(278) - p(387) - p(386) - p(426)
  p(437) = p(2)*p(177) - p(346) - p(391)
  p(438) = p(2)*p(178) - p(349) - p(400)
  p(439) = p(1)*p(281) - p(396) - p(392) - p(430)
  p(440) = p(21)*p(54) - p(370)
  p(441) = p(21)*p(55) - p(368)
  p(442) = p(2)*p(189) - p(405) - p(434)
  p(443) = p(1)*p(285) - p(402) - p(404) - p(400) - p(434) - p(430)
  p(444) = p(6)*p(116) - p(347) - p(413) - p(403)
  p(445) = p(7)*p(116) - p(348) - p(411) - p(401)
  p(446) = p(2)*p(182) - p(351) - p(409)
  p(447) = p(2)*p(191) - p(414) - p(443)
  p(448) = p(3)*p(183) - p(351) - p(418) - p(411) - p(413) - p(409)
  p(449) = p(6)*p(118) - p(351) - p(418) - p(412)
  p(450) = p(7)*p(118) - p(351) - p(418) - p(410)
  p(451) = p(2)*p(193) - p(418) - p(448)
  p(452) = p(1)*p(294) - p(418) - p(415) - p(448)
  p(453) = p(6)*p(119) - p(417)
  p(454) = p(7)*p(119) - p(416)
  p(455) = p(2)*p(186) - p(363)
  p(456) = p(3)*p(187) - p(385) - p(386) - p(384) - p(435) - p(436) - p(434)
  p(457) = p(2)*p(188) - p(388)
  p(458) = p(4)*p(187) - p(425) - p(426) - p(424)
  p(459) = p(2)*p(190) - p(406)
  p(460) = p(21)*p(66) - p(422) - p(423) - p(420)
  p(461) = p(2)*p(192) - p(415)
  p(462) = p(18)*p(71) - p(440) - p(441) - p(438)
  p(463) = p(2)*p(194) - p(452)
  p(464) = p(3)*p(194) - p(453) - p(454) - p(451)
  p(465) = p(1)*p(305) - p(464) - p(463)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine evbas
!**********************************************************************
!  The subroutine eliminates the 2-body terms in Bowman's approach
!**********************************************************************

  integer                   :: i
  real(rkp) ,dimension(466) :: b1 

  ! Pass P(0:465) to BM1(1:466)
  do i=1,466
    b1(i)=p(i-1)
  enddo

  ! Remove unconnected terms and 2-body terms and pass to B(1:430)
  b(1)=b1(4)

  do i=2,4
    b(i)=b1(i+4)
  enddo

  b(5)=b1(10)

  do i=6,11
    b(i)=b1(i+6)
  enddo

  b(12)=b1(19)
  b(13)=b1(21)

  do i=14,26
    b(i)=b1(i+9)
  enddo

  b(27)=b1(37)
  b(28)=b1(39)

  do i=29,53
    b(i)=b1(i+12)
  enddo

  b(54)=b1(67)
  b(55)=b1(69)
  b(56)=b1(71)

  do i=57,97
    b(i)=b1(i+16)
  enddo

  b(98)=b1(115)
  b(99)=b1(117)
  b(100)=b1(119)

  do i=101,166
    b(i)=b1(i+20)
  enddo

  b(167)=b1(188)
  b(168)=b1(190)
  b(169)=b1(192)
  b(170)=b1(194)

  do i=171,272
    b(i)=b1(i+25)
  enddo

  b(273)=b1(299)
  b(274)=b1(301)
  b(275)=b1(303)
  b(276)=b1(305)

  do i=277,425
    b(i)=b1(i+30)
  enddo

  b(426)=b1(457)
  b(427)=b1(459)
  b(428)=b1(461)
  b(429)=b1(463)
  b(430)=b1(465)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine EvdPdr
!**********************************************************************
!  The subroutine reads monomials(m) and calculates the
!  permutation-invariant polynomials(p)
!  For A4 with max. degree 9, the number of polynomials is nob.
!**********************************************************************

  integer i

  do i=1,6
    dpdr(i,0) = dmdr(i,0)
    dpdr(i,1) = dmdr(i,1) + dmdr(i,2) + dmdr(i,3) + dmdr(i,4) + dmdr(i,5) + dmdr(i,6)
    dpdr(i,2) = dmdr(i,7) + dmdr(i,8) + dmdr(i,9)
    dpdr(i,3) = dmdr(i,10) + dmdr(i,11) + dmdr(i,12) + dmdr(i,13) + dmdr(i,14) + dmdr(i,15) + dmdr(i,16) + dmdr(i,17) + dmdr(i,18) + dmdr(i,19) + dmdr(i,20) + dmdr(i,21)
    dpdr(i,4) = dpdr(i,1)*p(1) + p(1)*dpdr(i,1) - dpdr(i,3) - dpdr(i,2) - dpdr(i,3) - dpdr(i,2)
    dpdr(i,5) = dmdr(i,22) + dmdr(i,23) + dmdr(i,24) + dmdr(i,25) + dmdr(i,26) + dmdr(i,27) + dmdr(i,28) + dmdr(i,29) + dmdr(i,30) + dmdr(i,31) + dmdr(i,32) + dmdr(i,33)
    dpdr(i,6) = dmdr(i,34) + dmdr(i,35) + dmdr(i,36) + dmdr(i,37)
    dpdr(i,7) = dmdr(i,38) + dmdr(i,39) + dmdr(i,40) + dmdr(i,41)
    dpdr(i,8) = dpdr(i,1)*p(2) + p(1)*dpdr(i,2) - dpdr(i,5)
    dpdr(i,9) = dpdr(i,1)*p(3) + p(1)*dpdr(i,3) - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7) - dpdr(i,5) - dpdr(i,6) - dpdr(i,7)
    dpdr(i,10) = dpdr(i,1)*p(4) + p(1)*dpdr(i,4) - dpdr(i,9) - dpdr(i,8)
    dpdr(i,11) = dmdr(i,42) + dmdr(i,43) + dmdr(i,44)
    dpdr(i,12) = dmdr(i,45) + dmdr(i,46) + dmdr(i,47) + dmdr(i,48) + dmdr(i,49) + dmdr(i,50) + dmdr(i,51) + dmdr(i,52) + dmdr(i,53) + dmdr(i,54) + dmdr(i,55) + dmdr(i,56)
    dpdr(i,13) = dpdr(i,2)*p(3) + p(2)*dpdr(i,3) - dpdr(i,12)
    dpdr(i,14) = dpdr(i,1)*p(5) + p(1)*dpdr(i,5) - dpdr(i,12) - dpdr(i,11) - dpdr(i,13) - dpdr(i,12) - dpdr(i,11) - dpdr(i,11) - dpdr(i,11)
    dpdr(i,15) = dpdr(i,1)*p(6) + p(1)*dpdr(i,6) - dpdr(i,12)
    dpdr(i,16) = dpdr(i,1)*p(7) + p(1)*dpdr(i,7) - dpdr(i,12)
    dpdr(i,17) = dpdr(i,2)*p(2) + p(2)*dpdr(i,2) - dpdr(i,11) - dpdr(i,11)
    dpdr(i,18) = dpdr(i,3)*p(3) + p(3)*dpdr(i,3) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,15) - dpdr(i,16) - dpdr(i,14) - dpdr(i,12) - dpdr(i,11) - dpdr(i,12) - dpdr(i,11)
    dpdr(i,19) = dpdr(i,2)*p(4) + p(2)*dpdr(i,4) - dpdr(i,14)
    dpdr(i,20) = dpdr(i,3)*p(4) + p(3)*dpdr(i,4) - dpdr(i,15) - dpdr(i,16) - dpdr(i,13)
    dpdr(i,21) = dpdr(i,1)*p(10) + p(1)*dpdr(i,10) - dpdr(i,20) - dpdr(i,19)
    dpdr(i,22) = dmdr(i,57) + dmdr(i,58) + dmdr(i,59) + dmdr(i,60) + dmdr(i,61) + dmdr(i,62)
    dpdr(i,23) = dpdr(i,1)*p(11) + p(1)*dpdr(i,11) - dpdr(i,22)
    dpdr(i,24) = dpdr(i,2)*p(6) + p(2)*dpdr(i,6)
    dpdr(i,25) = dpdr(i,2)*p(7) + p(2)*dpdr(i,7)
    dpdr(i,26) = dpdr(i,1)*p(12) + p(1)*dpdr(i,12) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,22) - dpdr(i,22) - dpdr(i,22)
    dpdr(i,27) = dpdr(i,2)*p(5) + p(2)*dpdr(i,5) - dpdr(i,22) - dpdr(i,23) - dpdr(i,22)
    dpdr(i,28) = dpdr(i,3)*p(5) + p(3)*dpdr(i,5) - dpdr(i,22) - dpdr(i,26) - dpdr(i,24) - dpdr(i,25) - dpdr(i,23) - dpdr(i,22) - dpdr(i,24) - dpdr(i,25) - dpdr(i,23) - dpdr(i,22) - dpdr(i,22)
    dpdr(i,29) = dpdr(i,3)*p(6) + p(3)*dpdr(i,6) - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
    dpdr(i,30) = dpdr(i,3)*p(7) + p(3)*dpdr(i,7) - dpdr(i,22) - dpdr(i,26) - dpdr(i,22)
    dpdr(i,31) = dpdr(i,2)*p(9) + p(2)*dpdr(i,9) - dpdr(i,26) - dpdr(i,28)
    dpdr(i,32) = dpdr(i,1)*p(14) + p(1)*dpdr(i,14) - dpdr(i,26) - dpdr(i,23) - dpdr(i,28)
    dpdr(i,33) = dpdr(i,4)*p(6) + p(4)*dpdr(i,6) - dpdr(i,25)
    dpdr(i,34) = dpdr(i,4)*p(7) + p(4)*dpdr(i,7) - dpdr(i,24)
    dpdr(i,35) = dpdr(i,1)*p(17) + p(1)*dpdr(i,17) - dpdr(i,27)
    dpdr(i,36) = dpdr(i,1)*p(18) + p(1)*dpdr(i,18) - dpdr(i,29) - dpdr(i,30) - dpdr(i,28)
    dpdr(i,37) = dpdr(i,2)*p(10) + p(2)*dpdr(i,10) - dpdr(i,32)
    dpdr(i,38) = dpdr(i,3)*p(10) + p(3)*dpdr(i,10) - dpdr(i,33) - dpdr(i,34) - dpdr(i,31)
    dpdr(i,39) = dpdr(i,1)*p(21) + p(1)*dpdr(i,21) - dpdr(i,38) - dpdr(i,37)
    dpdr(i,40) = dmdr(i,63)
    dpdr(i,41) = dmdr(i,64) + dmdr(i,65) + dmdr(i,66) + dmdr(i,67) + dmdr(i,68) + dmdr(i,69) + dmdr(i,70) + dmdr(i,71) + dmdr(i,72) + dmdr(i,73) + dmdr(i,74) + dmdr(i,75) + dmdr(i,76) + dmdr(i,77) + dmdr(i,78) + dmdr(i,79) + dmdr(i,80) + dmdr(i,81) + dmdr(i,82) + dmdr(i,83) + dmdr(i,84) + dmdr(i,85) + dmdr(i,86) + dmdr(i,87)
    dpdr(i,42) = dpdr(i,1)*p(22) + p(1)*dpdr(i,22) - dpdr(i,40) - dpdr(i,41) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
    dpdr(i,43) = dpdr(i,2)*p(11) + p(2)*dpdr(i,11) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
    dpdr(i,44) = dpdr(i,2)*p(12) + p(2)*dpdr(i,12) - dpdr(i,41)
    dpdr(i,45) = dpdr(i,3)*p(11) + p(3)*dpdr(i,11) - dpdr(i,41)
    dpdr(i,46) = dpdr(i,5)*p(6) + p(5)*dpdr(i,6) - dpdr(i,41)
    dpdr(i,47) = dpdr(i,5)*p(7) + p(5)*dpdr(i,7) - dpdr(i,41)
    dpdr(i,48) = dpdr(i,6)*p(7) + p(6)*dpdr(i,7) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40) - dpdr(i,40)
    dpdr(i,49) = dpdr(i,4)*p(11) + p(4)*dpdr(i,11) - dpdr(i,42)
    dpdr(i,50) = dpdr(i,2)*p(15) + p(2)*dpdr(i,15) - dpdr(i,46)
    dpdr(i,51) = dpdr(i,2)*p(16) + p(2)*dpdr(i,16) - dpdr(i,47)
    dpdr(i,52) = dpdr(i,4)*p(12) + p(4)*dpdr(i,12) - dpdr(i,41) - dpdr(i,50) - dpdr(i,51)
    dpdr(i,53) = dpdr(i,2)*p(14) + p(2)*dpdr(i,14) - dpdr(i,42) - dpdr(i,49) - dpdr(i,42)
    dpdr(i,54) = dpdr(i,6)*p(6) + p(6)*dpdr(i,6) - dpdr(i,42) - dpdr(i,42)
    dpdr(i,55) = dpdr(i,7)*p(7) + p(7)*dpdr(i,7) - dpdr(i,42) - dpdr(i,42)
    dpdr(i,56) = dpdr(i,3)*p(17) + p(3)*dpdr(i,17) - dpdr(i,44)
    dpdr(i,57) = dpdr(i,2)*p(18) + p(2)*dpdr(i,18) - dpdr(i,48)
    dpdr(i,58) = dpdr(i,3)*p(14) + p(3)*dpdr(i,14) - dpdr(i,41) - dpdr(i,52) - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,45)
    dpdr(i,59) = dpdr(i,6)*p(9) + p(6)*dpdr(i,9) - dpdr(i,41) - dpdr(i,52) - dpdr(i,47)
    dpdr(i,60) = dpdr(i,7)*p(9) + p(7)*dpdr(i,9) - dpdr(i,41) - dpdr(i,52) - dpdr(i,46)
    dpdr(i,61) = dpdr(i,2)*p(20) + p(2)*dpdr(i,20) - dpdr(i,52) - dpdr(i,58)
    dpdr(i,62) = dpdr(i,1)*p(32) + p(1)*dpdr(i,32) - dpdr(i,52) - dpdr(i,49) - dpdr(i,58)
    dpdr(i,63) = dpdr(i,6)*p(10) + p(6)*dpdr(i,10) - dpdr(i,51)
    dpdr(i,64) = dpdr(i,7)*p(10) + p(7)*dpdr(i,10) - dpdr(i,50)
    dpdr(i,65) = dpdr(i,2)*p(17) + p(2)*dpdr(i,17) - dpdr(i,43)
    dpdr(i,66) = dpdr(i,3)*p(18) + p(3)*dpdr(i,18) - dpdr(i,46) - dpdr(i,47) - dpdr(i,45) - dpdr(i,59) - dpdr(i,60) - dpdr(i,58)
    dpdr(i,67) = dpdr(i,2)*p(19) + p(2)*dpdr(i,19) - dpdr(i,49)
    dpdr(i,68) = dpdr(i,1)*p(36) + p(1)*dpdr(i,36) - dpdr(i,59) - dpdr(i,60) - dpdr(i,58) - dpdr(i,57) - dpdr(i,66) - dpdr(i,66)
    dpdr(i,69) = dpdr(i,2)*p(21) + p(2)*dpdr(i,21) - dpdr(i,62)
    dpdr(i,70) = dpdr(i,3)*p(21) + p(3)*dpdr(i,21) - dpdr(i,63) - dpdr(i,64) - dpdr(i,61)
    dpdr(i,71) = dpdr(i,1)*p(39) + p(1)*dpdr(i,39) - dpdr(i,70) - dpdr(i,69)
    dpdr(i,72) = dpdr(i,40)*p(1) + p(40)*dpdr(i,1)
    dpdr(i,73) = dpdr(i,2)*p(22) + p(2)*dpdr(i,22) - dpdr(i,72)
    dpdr(i,74) = dpdr(i,6)*p(11) + p(6)*dpdr(i,11)
    dpdr(i,75) = dpdr(i,7)*p(11) + p(7)*dpdr(i,11)
    dpdr(i,76) = dpdr(i,3)*p(22) + p(3)*dpdr(i,22) - dpdr(i,72) - dpdr(i,74) - dpdr(i,75) - dpdr(i,72) - dpdr(i,72) - dpdr(i,72)
    dpdr(i,77) = dmdr(i,88) + dmdr(i,89) + dmdr(i,90) + dmdr(i,91) + dmdr(i,92) + dmdr(i,93) + dmdr(i,94) + dmdr(i,95) + dmdr(i,96) + dmdr(i,97) + dmdr(i,98) + dmdr(i,99) + dmdr(i,100) + dmdr(i,101) + dmdr(i,102) + dmdr(i,103) + dmdr(i,104) + dmdr(i,105) + dmdr(i,106) + dmdr(i,107) + dmdr(i,108) + dmdr(i,109) + dmdr(i,110) + dmdr(i,111)
    dpdr(i,78) = dpdr(i,1)*p(42) + p(1)*dpdr(i,42) - dpdr(i,72) - dpdr(i,76)
    dpdr(i,79) = dpdr(i,5)*p(11) + p(5)*dpdr(i,11) - dpdr(i,72) - dpdr(i,73) - dpdr(i,72)
    dpdr(i,80) = dpdr(i,2)*p(26) + p(2)*dpdr(i,26) - dpdr(i,76) - dpdr(i,77)
    dpdr(i,81) = dpdr(i,6)*p(12) + p(6)*dpdr(i,12) - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
    dpdr(i,82) = dpdr(i,7)*p(12) + p(7)*dpdr(i,12) - dpdr(i,72) - dpdr(i,76) - dpdr(i,72)
    dpdr(i,83) = dpdr(i,8)*p(11) + p(8)*dpdr(i,11) - dpdr(i,72)
    dpdr(i,84) = dpdr(i,6)*p(17) + p(6)*dpdr(i,17)
    dpdr(i,85) = dpdr(i,7)*p(17) + p(7)*dpdr(i,17)
    dpdr(i,86) = dpdr(i,9)*p(11) + p(9)*dpdr(i,11) - dpdr(i,76) - dpdr(i,77)
    dpdr(i,87) = dpdr(i,2)*p(29) + p(2)*dpdr(i,29) - dpdr(i,81)
    dpdr(i,88) = dpdr(i,6)*p(14) + p(6)*dpdr(i,14) - dpdr(i,75) - dpdr(i,75)
    dpdr(i,89) = dpdr(i,2)*p(30) + p(2)*dpdr(i,30) - dpdr(i,82)
    dpdr(i,90) = dpdr(i,7)*p(14) + p(7)*dpdr(i,14) - dpdr(i,74) - dpdr(i,74)
    dpdr(i,91) = dpdr(i,1)*p(48) + p(1)*dpdr(i,48) - dpdr(i,76) - dpdr(i,81) - dpdr(i,82)
    dpdr(i,92) = dpdr(i,10)*p(11) + p(10)*dpdr(i,11) - dpdr(i,78)
    dpdr(i,93) = dpdr(i,2)*p(33) + p(2)*dpdr(i,33) - dpdr(i,88)
    dpdr(i,94) = dpdr(i,2)*p(34) + p(2)*dpdr(i,34) - dpdr(i,90)
    dpdr(i,95) = dpdr(i,10)*p(12) + p(10)*dpdr(i,12) - dpdr(i,77) - dpdr(i,93) - dpdr(i,94)
    dpdr(i,96) = dpdr(i,2)*p(27) + p(2)*dpdr(i,27) - dpdr(i,73) - dpdr(i,79)
    dpdr(i,97) = dpdr(i,2)*p(28) + p(2)*dpdr(i,28) - dpdr(i,76) - dpdr(i,86)
    dpdr(i,98) = dpdr(i,1)*p(53) + p(1)*dpdr(i,53) - dpdr(i,80) - dpdr(i,79) - dpdr(i,97)
    dpdr(i,99) = dpdr(i,1)*p(54) + p(1)*dpdr(i,54) - dpdr(i,81)
    dpdr(i,100) = dpdr(i,1)*p(55) + p(1)*dpdr(i,55) - dpdr(i,82)
    dpdr(i,101) = dpdr(i,5)*p(18) + p(5)*dpdr(i,18) - dpdr(i,76) - dpdr(i,91) - dpdr(i,87) - dpdr(i,89) - dpdr(i,86)
    dpdr(i,102) = dpdr(i,6)*p(18) + p(6)*dpdr(i,18) - dpdr(i,74) - dpdr(i,90)
    dpdr(i,103) = dpdr(i,7)*p(18) + p(7)*dpdr(i,18) - dpdr(i,75) - dpdr(i,88)
    dpdr(i,104) = dpdr(i,2)*p(31) + p(2)*dpdr(i,31) - dpdr(i,77) - dpdr(i,86)
    dpdr(i,105) = dpdr(i,2)*p(36) + p(2)*dpdr(i,36) - dpdr(i,91) - dpdr(i,101)
    dpdr(i,106) = dpdr(i,3)*p(32) + p(3)*dpdr(i,32) - dpdr(i,77) - dpdr(i,95) - dpdr(i,88) - dpdr(i,90) - dpdr(i,86)
    dpdr(i,107) = dpdr(i,4)*p(29) + p(4)*dpdr(i,29) - dpdr(i,82) - dpdr(i,80) - dpdr(i,99)
    dpdr(i,108) = dpdr(i,4)*p(30) + p(4)*dpdr(i,30) - dpdr(i,81) - dpdr(i,80) - dpdr(i,100)
    dpdr(i,109) = dpdr(i,2)*p(38) + p(2)*dpdr(i,38) - dpdr(i,95) - dpdr(i,106)
    dpdr(i,110) = dpdr(i,1)*p(62) + p(1)*dpdr(i,62) - dpdr(i,95) - dpdr(i,92) - dpdr(i,106)
    dpdr(i,111) = dpdr(i,6)*p(21) + p(6)*dpdr(i,21) - dpdr(i,94)
    dpdr(i,112) = dpdr(i,7)*p(21) + p(7)*dpdr(i,21) - dpdr(i,93)
    dpdr(i,113) = dpdr(i,1)*p(65) + p(1)*dpdr(i,65) - dpdr(i,96)
    dpdr(i,114) = dpdr(i,1)*p(66) + p(1)*dpdr(i,66) - dpdr(i,102) - dpdr(i,103) - dpdr(i,101)
    dpdr(i,115) = dpdr(i,2)*p(37) + p(2)*dpdr(i,37) - dpdr(i,92)
    dpdr(i,116) = dpdr(i,10)*p(18) + p(10)*dpdr(i,18) - dpdr(i,99) - dpdr(i,100) - dpdr(i,97)
    dpdr(i,117) = dpdr(i,2)*p(39) + p(2)*dpdr(i,39) - dpdr(i,110)
    dpdr(i,118) = dpdr(i,3)*p(39) + p(3)*dpdr(i,39) - dpdr(i,111) - dpdr(i,112) - dpdr(i,109)
    dpdr(i,119) = dpdr(i,1)*p(71) + p(1)*dpdr(i,71) - dpdr(i,118) - dpdr(i,117)
    dpdr(i,120) = dpdr(i,40)*p(2) + p(40)*dpdr(i,2)
    dpdr(i,121) = dpdr(i,40)*p(3) + p(40)*dpdr(i,3)
    dpdr(i,122) = dpdr(i,40)*p(4) + p(40)*dpdr(i,4)
    dpdr(i,123) = dpdr(i,11)*p(12) + p(11)*dpdr(i,12) - dpdr(i,121)
    dpdr(i,124) = dpdr(i,2)*p(42) + p(2)*dpdr(i,42) - dpdr(i,122)
    dpdr(i,125) = dpdr(i,6)*p(22) + p(6)*dpdr(i,22) - dpdr(i,121)
    dpdr(i,126) = dpdr(i,7)*p(22) + p(7)*dpdr(i,22) - dpdr(i,121)
    dpdr(i,127) = dpdr(i,2)*p(41) + p(2)*dpdr(i,41) - dpdr(i,121) - dpdr(i,123) - dpdr(i,121)
    dpdr(i,128) = dpdr(i,6)*p(23) + p(6)*dpdr(i,23) - dpdr(i,123)
    dpdr(i,129) = dpdr(i,7)*p(23) + p(7)*dpdr(i,23) - dpdr(i,123)
    dpdr(i,130) = dpdr(i,3)*p(41) + p(3)*dpdr(i,41) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125) - dpdr(i,128) - dpdr(i,126) - dpdr(i,129) - dpdr(i,124) - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,125) - dpdr(i,126) - dpdr(i,124) - dpdr(i,123) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,122) - dpdr(i,121) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120) - dpdr(i,120)
    dpdr(i,131) = dpdr(i,3)*p(42) + p(3)*dpdr(i,42) - dpdr(i,121) - dpdr(i,125) - dpdr(i,126) - dpdr(i,121)
    dpdr(i,132) = dpdr(i,4)*p(41) + p(4)*dpdr(i,41) - dpdr(i,121) - dpdr(i,131) - dpdr(i,128) - dpdr(i,129) - dpdr(i,127) - dpdr(i,121)
    dpdr(i,133) = dpdr(i,1)*p(78) + p(1)*dpdr(i,78) - dpdr(i,122) - dpdr(i,131)
    dpdr(i,134) = dpdr(i,11)*p(11) + p(11)*dpdr(i,11) - dpdr(i,120) - dpdr(i,120)
    dpdr(i,135) = dpdr(i,2)*p(48) + p(2)*dpdr(i,48) - dpdr(i,130)
    dpdr(i,136) = dpdr(i,11)*p(17) + p(11)*dpdr(i,17) - dpdr(i,120)
    dpdr(i,137) = dpdr(i,2)*p(44) + p(2)*dpdr(i,44) - dpdr(i,123)
    dpdr(i,138) = dpdr(i,2)*p(45) + p(2)*dpdr(i,45) - dpdr(i,121)
    dpdr(i,139) = dpdr(i,11)*p(14) + p(11)*dpdr(i,14) - dpdr(i,122) - dpdr(i,124) - dpdr(i,122)
    dpdr(i,140) = dpdr(i,6)*p(27) + p(6)*dpdr(i,27) - dpdr(i,127)
    dpdr(i,141) = dpdr(i,2)*p(54) + p(2)*dpdr(i,54)
    dpdr(i,142) = dpdr(i,7)*p(27) + p(7)*dpdr(i,27) - dpdr(i,127)
    dpdr(i,143) = dpdr(i,2)*p(52) + p(2)*dpdr(i,52) - dpdr(i,131) - dpdr(i,132)
    dpdr(i,144) = dpdr(i,1)*p(81) + p(1)*dpdr(i,81) - dpdr(i,125) - dpdr(i,141) - dpdr(i,135) - dpdr(i,125)
    dpdr(i,145) = dpdr(i,2)*p(55) + p(2)*dpdr(i,55)
    dpdr(i,146) = dpdr(i,1)*p(82) + p(1)*dpdr(i,82) - dpdr(i,126) - dpdr(i,135) - dpdr(i,145) - dpdr(i,126)
    dpdr(i,147) = dpdr(i,11)*p(18) + p(11)*dpdr(i,18) - dpdr(i,130)
    dpdr(i,148) = dpdr(i,6)*p(28) + p(6)*dpdr(i,28) - dpdr(i,129) - dpdr(i,123) - dpdr(i,143)
    dpdr(i,149) = dpdr(i,7)*p(28) + p(7)*dpdr(i,28) - dpdr(i,128) - dpdr(i,123) - dpdr(i,143)
    dpdr(i,150) = dpdr(i,6)*p(30) + p(6)*dpdr(i,30) - dpdr(i,121) - dpdr(i,146)
    dpdr(i,151) = dpdr(i,11)*p(19) + p(11)*dpdr(i,19) - dpdr(i,122)
    dpdr(i,152) = dpdr(i,2)*p(50) + p(2)*dpdr(i,50) - dpdr(i,128)
    dpdr(i,153) = dpdr(i,2)*p(51) + p(2)*dpdr(i,51) - dpdr(i,129)
    dpdr(i,154) = dpdr(i,11)*p(20) + p(11)*dpdr(i,20) - dpdr(i,131) - dpdr(i,132)
    dpdr(i,155) = dpdr(i,2)*p(59) + p(2)*dpdr(i,59) - dpdr(i,144) - dpdr(i,148)
    dpdr(i,156) = dpdr(i,6)*p(32) + p(6)*dpdr(i,32) - dpdr(i,129)
    dpdr(i,157) = dpdr(i,2)*p(60) + p(2)*dpdr(i,60) - dpdr(i,146) - dpdr(i,149)
    dpdr(i,158) = dpdr(i,7)*p(32) + p(7)*dpdr(i,32) - dpdr(i,128)
    dpdr(i,159) = dpdr(i,6)*p(34) + p(6)*dpdr(i,34) - dpdr(i,122) - dpdr(i,145) - dpdr(i,122)
    dpdr(i,160) = dpdr(i,11)*p(21) + p(11)*dpdr(i,21) - dpdr(i,133)
    dpdr(i,161) = dpdr(i,2)*p(63) + p(2)*dpdr(i,63) - dpdr(i,156)
    dpdr(i,162) = dpdr(i,2)*p(64) + p(2)*dpdr(i,64) - dpdr(i,158)
    dpdr(i,163) = dpdr(i,12)*p(21) + p(12)*dpdr(i,21) - dpdr(i,132) - dpdr(i,161) - dpdr(i,162)
    dpdr(i,164) = dpdr(i,2)*p(53) + p(2)*dpdr(i,53) - dpdr(i,124) - dpdr(i,139)
    dpdr(i,165) = dpdr(i,2)*p(58) + p(2)*dpdr(i,58) - dpdr(i,131) - dpdr(i,154)
    dpdr(i,166) = dpdr(i,3)*p(54) + p(3)*dpdr(i,54) - dpdr(i,125) - dpdr(i,144)
    dpdr(i,167) = dpdr(i,3)*p(55) + p(3)*dpdr(i,55) - dpdr(i,126) - dpdr(i,146)
    dpdr(i,168) = dpdr(i,3)*p(65) + p(3)*dpdr(i,65) - dpdr(i,137)
    dpdr(i,169) = dpdr(i,17)*p(18) + p(17)*dpdr(i,18) - dpdr(i,135)
    dpdr(i,170) = dpdr(i,1)*p(98) + p(1)*dpdr(i,98) - dpdr(i,143) - dpdr(i,139) - dpdr(i,165)
    dpdr(i,171) = dpdr(i,4)*p(54) + p(4)*dpdr(i,54) - dpdr(i,135)
    dpdr(i,172) = dpdr(i,4)*p(55) + p(4)*dpdr(i,55) - dpdr(i,135)
    dpdr(i,173) = dpdr(i,2)*p(66) + p(2)*dpdr(i,66) - dpdr(i,150)
    dpdr(i,174) = dpdr(i,1)*p(101) + p(1)*dpdr(i,101) - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,173) - dpdr(i,165) - dpdr(i,147)
    dpdr(i,175) = dpdr(i,1)*p(102) + p(1)*dpdr(i,102) - dpdr(i,150) - dpdr(i,148) - dpdr(i,166)
    dpdr(i,176) = dpdr(i,1)*p(103) + p(1)*dpdr(i,103) - dpdr(i,150) - dpdr(i,149) - dpdr(i,167)
    dpdr(i,177) = dpdr(i,2)*p(61) + p(2)*dpdr(i,61) - dpdr(i,132) - dpdr(i,154)
    dpdr(i,178) = dpdr(i,2)*p(68) + p(2)*dpdr(i,68) - dpdr(i,159) - dpdr(i,174)
    dpdr(i,179) = dpdr(i,3)*p(62) + p(3)*dpdr(i,62) - dpdr(i,132) - dpdr(i,163) - dpdr(i,156) - dpdr(i,158) - dpdr(i,154)
    dpdr(i,180) = dpdr(i,6)*p(38) + p(6)*dpdr(i,38) - dpdr(i,132) - dpdr(i,163) - dpdr(i,157)
    dpdr(i,181) = dpdr(i,7)*p(38) + p(7)*dpdr(i,38) - dpdr(i,132) - dpdr(i,163) - dpdr(i,155)
    dpdr(i,182) = dpdr(i,2)*p(70) + p(2)*dpdr(i,70) - dpdr(i,163) - dpdr(i,179)
    dpdr(i,183) = dpdr(i,1)*p(110) + p(1)*dpdr(i,110) - dpdr(i,163) - dpdr(i,160) - dpdr(i,179)
    dpdr(i,184) = dpdr(i,6)*p(39) + p(6)*dpdr(i,39) - dpdr(i,162)
    dpdr(i,185) = dpdr(i,7)*p(39) + p(7)*dpdr(i,39) - dpdr(i,161)
    dpdr(i,186) = dpdr(i,2)*p(65) + p(2)*dpdr(i,65) - dpdr(i,136)
    dpdr(i,187) = dpdr(i,3)*p(66) + p(3)*dpdr(i,66) - dpdr(i,148) - dpdr(i,149) - dpdr(i,147) - dpdr(i,175) - dpdr(i,176) - dpdr(i,174)
    dpdr(i,188) = dpdr(i,2)*p(67) + p(2)*dpdr(i,67) - dpdr(i,151)
    dpdr(i,189) = dpdr(i,4)*p(66) + p(4)*dpdr(i,66) - dpdr(i,166) - dpdr(i,167) - dpdr(i,165)
    dpdr(i,190) = dpdr(i,2)*p(69) + p(2)*dpdr(i,69) - dpdr(i,160)
    dpdr(i,191) = dpdr(i,18)*p(21) + p(18)*dpdr(i,21) - dpdr(i,171) - dpdr(i,172) - dpdr(i,169)
    dpdr(i,192) = dpdr(i,2)*p(71) + p(2)*dpdr(i,71) - dpdr(i,183)
    dpdr(i,193) = dpdr(i,3)*p(71) + p(3)*dpdr(i,71) - dpdr(i,184) - dpdr(i,185) - dpdr(i,182)
    dpdr(i,194) = dpdr(i,1)*p(119) + p(1)*dpdr(i,119) - dpdr(i,193) - dpdr(i,192)
    dpdr(i,195) = dpdr(i,40)*p(5) + p(40)*dpdr(i,5)
    dpdr(i,196) = dpdr(i,40)*p(6) + p(40)*dpdr(i,6)
    dpdr(i,197) = dpdr(i,40)*p(7) + p(40)*dpdr(i,7)
    dpdr(i,198) = dpdr(i,40)*p(8) + p(40)*dpdr(i,8)
    dpdr(i,199) = dpdr(i,40)*p(9) + p(40)*dpdr(i,9)
    dpdr(i,200) = dpdr(i,40)*p(10) + p(40)*dpdr(i,10)
    dpdr(i,201) = dpdr(i,11)*p(22) + p(11)*dpdr(i,22) - dpdr(i,195)
    dpdr(i,202) = dpdr(i,12)*p(22) + p(12)*dpdr(i,22) - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197) - dpdr(i,195) - dpdr(i,196) - dpdr(i,197)
    dpdr(i,203) = dpdr(i,17)*p(22) + p(17)*dpdr(i,22) - dpdr(i,198)
    dpdr(i,204) = dpdr(i,6)*p(43) + p(6)*dpdr(i,43)
    dpdr(i,205) = dpdr(i,7)*p(43) + p(7)*dpdr(i,43)
    dpdr(i,206) = dpdr(i,11)*p(26) + p(11)*dpdr(i,26) - dpdr(i,199) - dpdr(i,202)
    dpdr(i,207) = dpdr(i,2)*p(76) + p(2)*dpdr(i,76) - dpdr(i,199) - dpdr(i,202)
    dpdr(i,208) = dpdr(i,2)*p(78) + p(2)*dpdr(i,78) - dpdr(i,200)
    dpdr(i,209) = dpdr(i,6)*p(41) + p(6)*dpdr(i,41) - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
    dpdr(i,210) = dpdr(i,6)*p(42) + p(6)*dpdr(i,42) - dpdr(i,197) - dpdr(i,197) - dpdr(i,197)
    dpdr(i,211) = dpdr(i,7)*p(41) + p(7)*dpdr(i,41) - dpdr(i,199) - dpdr(i,195) - dpdr(i,202) - dpdr(i,195)
    dpdr(i,212) = dpdr(i,7)*p(42) + p(7)*dpdr(i,42) - dpdr(i,196) - dpdr(i,196) - dpdr(i,196)
    dpdr(i,213) = dpdr(i,11)*p(29) + p(11)*dpdr(i,29) - dpdr(i,209)
    dpdr(i,214) = dpdr(i,11)*p(30) + p(11)*dpdr(i,30) - dpdr(i,211)
    dpdr(i,215) = dpdr(i,18)*p(22) + p(18)*dpdr(i,22) - dpdr(i,199) - dpdr(i,213) - dpdr(i,214)
    dpdr(i,216) = dpdr(i,2)*p(77) + p(2)*dpdr(i,77) - dpdr(i,199) - dpdr(i,206)
    dpdr(i,217) = dpdr(i,6)*p(49) + p(6)*dpdr(i,49) - dpdr(i,205)
    dpdr(i,218) = dpdr(i,7)*p(49) + p(7)*dpdr(i,49) - dpdr(i,204)
    dpdr(i,219) = dpdr(i,3)*p(77) + p(3)*dpdr(i,77) - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,209) - dpdr(i,217) - dpdr(i,211) - dpdr(i,218) - dpdr(i,207) - dpdr(i,204) - dpdr(i,205) - dpdr(i,200) - dpdr(i,199) - dpdr(i,198) - dpdr(i,200) - dpdr(i,198) - dpdr(i,200) - dpdr(i,198)
    dpdr(i,220) = dpdr(i,3)*p(78) + p(3)*dpdr(i,78) - dpdr(i,199) - dpdr(i,210) - dpdr(i,212)
    dpdr(i,221) = dpdr(i,10)*p(41) + p(10)*dpdr(i,41) - dpdr(i,199) - dpdr(i,220) - dpdr(i,217) - dpdr(i,218) - dpdr(i,216)
    dpdr(i,222) = dpdr(i,1)*p(133) + p(1)*dpdr(i,133) - dpdr(i,200) - dpdr(i,220)
    dpdr(i,223) = dpdr(i,11)*p(27) + p(11)*dpdr(i,27) - dpdr(i,195) - dpdr(i,203)
    dpdr(i,224) = dpdr(i,1)*p(134) + p(1)*dpdr(i,134) - dpdr(i,201)
    dpdr(i,225) = dpdr(i,2)*p(81) + p(2)*dpdr(i,81) - dpdr(i,209)
    dpdr(i,226) = dpdr(i,2)*p(80) + p(2)*dpdr(i,80) - dpdr(i,202) - dpdr(i,206)
    dpdr(i,227) = dpdr(i,2)*p(82) + p(2)*dpdr(i,82) - dpdr(i,211)
    dpdr(i,228) = dpdr(i,2)*p(91) + p(2)*dpdr(i,91) - dpdr(i,215) - dpdr(i,219)
    dpdr(i,229) = dpdr(i,11)*p(28) + p(11)*dpdr(i,28) - dpdr(i,199) - dpdr(i,207)
    dpdr(i,230) = dpdr(i,6)*p(53) + p(6)*dpdr(i,53) - dpdr(i,205)
    dpdr(i,231) = dpdr(i,5)*p(54) + p(5)*dpdr(i,54) - dpdr(i,209)
    dpdr(i,232) = dpdr(i,7)*p(53) + p(7)*dpdr(i,53) - dpdr(i,204)
    dpdr(i,233) = dpdr(i,7)*p(54) + p(7)*dpdr(i,54) - dpdr(i,196)
    dpdr(i,234) = dpdr(i,5)*p(55) + p(5)*dpdr(i,55) - dpdr(i,211)
    dpdr(i,235) = dpdr(i,6)*p(55) + p(6)*dpdr(i,55) - dpdr(i,197)
    dpdr(i,236) = dpdr(i,11)*p(35) + p(11)*dpdr(i,35) - dpdr(i,198)
    dpdr(i,237) = dpdr(i,6)*p(65) + p(6)*dpdr(i,65)
    dpdr(i,238) = dpdr(i,7)*p(65) + p(7)*dpdr(i,65)
    dpdr(i,239) = dpdr(i,2)*p(86) + p(2)*dpdr(i,86) - dpdr(i,199) - dpdr(i,229)
    dpdr(i,240) = dpdr(i,1)*p(139) + p(1)*dpdr(i,139) - dpdr(i,206) - dpdr(i,229) - dpdr(i,224)
    dpdr(i,241) = dpdr(i,17)*p(29) + p(17)*dpdr(i,29) - dpdr(i,225)
    dpdr(i,242) = dpdr(i,2)*p(99) + p(2)*dpdr(i,99) - dpdr(i,231)
    dpdr(i,243) = dpdr(i,17)*p(30) + p(17)*dpdr(i,30) - dpdr(i,227)
    dpdr(i,244) = dpdr(i,2)*p(95) + p(2)*dpdr(i,95) - dpdr(i,220) - dpdr(i,221)
    dpdr(i,245) = dpdr(i,4)*p(81) + p(4)*dpdr(i,81) - dpdr(i,202) - dpdr(i,242) - dpdr(i,227)
    dpdr(i,246) = dpdr(i,2)*p(100) + p(2)*dpdr(i,100) - dpdr(i,234)
    dpdr(i,247) = dpdr(i,4)*p(82) + p(4)*dpdr(i,82) - dpdr(i,202) - dpdr(i,225) - dpdr(i,246)
    dpdr(i,248) = dpdr(i,11)*p(36) + p(11)*dpdr(i,36) - dpdr(i,215) - dpdr(i,219)
    dpdr(i,249) = dpdr(i,2)*p(102) + p(2)*dpdr(i,102) - dpdr(i,233)
    dpdr(i,250) = dpdr(i,6)*p(58) + p(6)*dpdr(i,58) - dpdr(i,206) - dpdr(i,214) - dpdr(i,244) - dpdr(i,214)
    dpdr(i,251) = dpdr(i,2)*p(103) + p(2)*dpdr(i,103) - dpdr(i,235)
    dpdr(i,252) = dpdr(i,7)*p(58) + p(7)*dpdr(i,58) - dpdr(i,213) - dpdr(i,206) - dpdr(i,244) - dpdr(i,213)
    dpdr(i,253) = dpdr(i,1)*p(150) + p(1)*dpdr(i,150) - dpdr(i,215) - dpdr(i,233) - dpdr(i,235)
    dpdr(i,254) = dpdr(i,11)*p(37) + p(11)*dpdr(i,37) - dpdr(i,200)
    dpdr(i,255) = dpdr(i,2)*p(93) + p(2)*dpdr(i,93) - dpdr(i,217)
    dpdr(i,256) = dpdr(i,2)*p(94) + p(2)*dpdr(i,94) - dpdr(i,218)
    dpdr(i,257) = dpdr(i,11)*p(38) + p(11)*dpdr(i,38) - dpdr(i,220) - dpdr(i,221)
    dpdr(i,258) = dpdr(i,2)*p(107) + p(2)*dpdr(i,107) - dpdr(i,245) - dpdr(i,250)
    dpdr(i,259) = dpdr(i,6)*p(62) + p(6)*dpdr(i,62) - dpdr(i,218)
    dpdr(i,260) = dpdr(i,2)*p(108) + p(2)*dpdr(i,108) - dpdr(i,247) - dpdr(i,252)
    dpdr(i,261) = dpdr(i,7)*p(62) + p(7)*dpdr(i,62) - dpdr(i,217)
    dpdr(i,262) = dpdr(i,6)*p(64) + p(6)*dpdr(i,64) - dpdr(i,200) - dpdr(i,246) - dpdr(i,200)
    dpdr(i,263) = dpdr(i,11)*p(39) + p(11)*dpdr(i,39) - dpdr(i,222)
    dpdr(i,264) = dpdr(i,2)*p(111) + p(2)*dpdr(i,111) - dpdr(i,259)
    dpdr(i,265) = dpdr(i,2)*p(112) + p(2)*dpdr(i,112) - dpdr(i,261)
    dpdr(i,266) = dpdr(i,12)*p(39) + p(12)*dpdr(i,39) - dpdr(i,221) - dpdr(i,264) - dpdr(i,265)
    dpdr(i,267) = dpdr(i,2)*p(98) + p(2)*dpdr(i,98) - dpdr(i,208) - dpdr(i,240)
    dpdr(i,268) = dpdr(i,6)*p(54) + p(6)*dpdr(i,54) - dpdr(i,210)
    dpdr(i,269) = dpdr(i,7)*p(55) + p(7)*dpdr(i,55) - dpdr(i,212)
    dpdr(i,270) = dpdr(i,2)*p(96) + p(2)*dpdr(i,96) - dpdr(i,203) - dpdr(i,223)
    dpdr(i,271) = dpdr(i,2)*p(97) + p(2)*dpdr(i,97) - dpdr(i,207) - dpdr(i,229)
    dpdr(i,272) = dpdr(i,2)*p(101) + p(2)*dpdr(i,101) - dpdr(i,215) - dpdr(i,248)
    dpdr(i,273) = dpdr(i,2)*p(106) + p(2)*dpdr(i,106) - dpdr(i,220) - dpdr(i,257)
    dpdr(i,274) = dpdr(i,6)*p(59) + p(6)*dpdr(i,59) - dpdr(i,220) - dpdr(i,215) - dpdr(i,209)
    dpdr(i,275) = dpdr(i,7)*p(60) + p(7)*dpdr(i,60) - dpdr(i,220) - dpdr(i,215) - dpdr(i,211)
    dpdr(i,276) = dpdr(i,5)*p(66) + p(5)*dpdr(i,66) - dpdr(i,215) - dpdr(i,253) - dpdr(i,249) - dpdr(i,251) - dpdr(i,248)
    dpdr(i,277) = dpdr(i,6)*p(66) + p(6)*dpdr(i,66) - dpdr(i,213) - dpdr(i,252)
    dpdr(i,278) = dpdr(i,7)*p(66) + p(7)*dpdr(i,66) - dpdr(i,214) - dpdr(i,250)
    dpdr(i,279) = dpdr(i,2)*p(104) + p(2)*dpdr(i,104) - dpdr(i,216) - dpdr(i,239)
    dpdr(i,280) = dpdr(i,2)*p(105) + p(2)*dpdr(i,105) - dpdr(i,219) - dpdr(i,248)
    dpdr(i,281) = dpdr(i,1)*p(170) + p(1)*dpdr(i,170) - dpdr(i,244) - dpdr(i,240) - dpdr(i,273)
    dpdr(i,282) = dpdr(i,10)*p(54) + p(10)*dpdr(i,54) - dpdr(i,227)
    dpdr(i,283) = dpdr(i,10)*p(55) + p(10)*dpdr(i,55) - dpdr(i,225)
    dpdr(i,284) = dpdr(i,2)*p(114) + p(2)*dpdr(i,114) - dpdr(i,253) - dpdr(i,276)
    dpdr(i,285) = dpdr(i,1)*p(174) + p(1)*dpdr(i,174) - dpdr(i,250) - dpdr(i,252) - dpdr(i,248) - dpdr(i,276) - dpdr(i,273)
    dpdr(i,286) = dpdr(i,6)*p(68) + p(6)*dpdr(i,68) - dpdr(i,217) - dpdr(i,261) - dpdr(i,251)
    dpdr(i,287) = dpdr(i,7)*p(68) + p(7)*dpdr(i,68) - dpdr(i,218) - dpdr(i,259) - dpdr(i,249)
    dpdr(i,288) = dpdr(i,2)*p(109) + p(2)*dpdr(i,109) - dpdr(i,221) - dpdr(i,257)
    dpdr(i,289) = dpdr(i,2)*p(116) + p(2)*dpdr(i,116) - dpdr(i,262) - dpdr(i,285)
    dpdr(i,290) = dpdr(i,3)*p(110) + p(3)*dpdr(i,110) - dpdr(i,221) - dpdr(i,266) - dpdr(i,259) - dpdr(i,261) - dpdr(i,257)
    dpdr(i,291) = dpdr(i,6)*p(70) + p(6)*dpdr(i,70) - dpdr(i,221) - dpdr(i,266) - dpdr(i,260)
    dpdr(i,292) = dpdr(i,7)*p(70) + p(7)*dpdr(i,70) - dpdr(i,221) - dpdr(i,266) - dpdr(i,258)
    dpdr(i,293) = dpdr(i,2)*p(118) + p(2)*dpdr(i,118) - dpdr(i,266) - dpdr(i,290)
    dpdr(i,294) = dpdr(i,1)*p(183) + p(1)*dpdr(i,183) - dpdr(i,266) - dpdr(i,263) - dpdr(i,290)
    dpdr(i,295) = dpdr(i,6)*p(71) + p(6)*dpdr(i,71) - dpdr(i,265)
    dpdr(i,296) = dpdr(i,7)*p(71) + p(7)*dpdr(i,71) - dpdr(i,264)
    dpdr(i,297) = dpdr(i,1)*p(186) + p(1)*dpdr(i,186) - dpdr(i,270)
    dpdr(i,298) = dpdr(i,1)*p(187) + p(1)*dpdr(i,187) - dpdr(i,277) - dpdr(i,278) - dpdr(i,276)
    dpdr(i,299) = dpdr(i,2)*p(115) + p(2)*dpdr(i,115) - dpdr(i,254)
    dpdr(i,300) = dpdr(i,1)*p(189) + p(1)*dpdr(i,189) - dpdr(i,286) - dpdr(i,287) - dpdr(i,285) - dpdr(i,284) - dpdr(i,298)
    dpdr(i,301) = dpdr(i,2)*p(117) + p(2)*dpdr(i,117) - dpdr(i,263)
    dpdr(i,302) = dpdr(i,18)*p(39) + p(18)*dpdr(i,39) - dpdr(i,282) - dpdr(i,283) - dpdr(i,280)
    dpdr(i,303) = dpdr(i,2)*p(119) + p(2)*dpdr(i,119) - dpdr(i,294)
    dpdr(i,304) = dpdr(i,3)*p(119) + p(3)*dpdr(i,119) - dpdr(i,295) - dpdr(i,296) - dpdr(i,293)
    dpdr(i,305) = dpdr(i,1)*p(194) + p(1)*dpdr(i,194) - dpdr(i,304) - dpdr(i,303)
    dpdr(i,306) = dpdr(i,40)*p(11) + p(40)*dpdr(i,11)  
    dpdr(i,307) = dpdr(i,40)*p(12) + p(40)*dpdr(i,12)  
    dpdr(i,308) = dpdr(i,40)*p(17) + p(40)*dpdr(i,17)  
    dpdr(i,309) = dpdr(i,40)*p(13) + p(40)*dpdr(i,13)  
    dpdr(i,310) = dpdr(i,40)*p(14) + p(40)*dpdr(i,14)  
    dpdr(i,311) = dpdr(i,40)*p(15) + p(40)*dpdr(i,15)  
    dpdr(i,312) = dpdr(i,40)*p(16) + p(40)*dpdr(i,16)  
    dpdr(i,313) = dpdr(i,40)*p(18) + p(40)*dpdr(i,18)  
    dpdr(i,314) = dpdr(i,40)*p(19) + p(40)*dpdr(i,19)  
    dpdr(i,315) = dpdr(i,40)*p(20) + p(40)*dpdr(i,20)  
    dpdr(i,316) = dpdr(i,40)*p(21) + p(40)*dpdr(i,21)  
    dpdr(i,317) = dpdr(i,11)*p(42) + p(11)*dpdr(i,42) - dpdr(i,310)
    dpdr(i,318) = dpdr(i,11)*p(44) + p(11)*dpdr(i,44) - dpdr(i,307)
    dpdr(i,319) = dpdr(i,12)*p(43) + p(12)*dpdr(i,43) - dpdr(i,309) -  dpdr(i,318)
    dpdr(i,320) = dpdr(i,17)*p(42) + p(17)*dpdr(i,42) - dpdr(i,314)
    dpdr(i,321) = dpdr(i,2)*p(125) + p(2)*dpdr(i,125) - dpdr(i,311)
    dpdr(i,322) = dpdr(i,2)*p(126) + p(2)*dpdr(i,126) - dpdr(i,312)
    dpdr(i,323) = dpdr(i,11)*p(48) + p(11)*dpdr(i,48) - dpdr(i,313)
    dpdr(i,324) = dpdr(i,12)*p(42) + p(12)*dpdr(i,42) - dpdr(i,311) -  dpdr(i,312) - dpdr(i,307) - dpdr(i,307)
    dpdr(i,325) = dpdr(i,6)*p(79)  + p(6)*dpdr(i,79) - dpdr(i,319)
    dpdr(i,326) = dpdr(i,11)*p(54) + p(11)*dpdr(i,54)  
    dpdr(i,327) = dpdr(i,7)*p(79)  + p(7)*dpdr(i,79) - dpdr(i,319)
    dpdr(i,328) = dpdr(i,2)*p(131) + p(2)*dpdr(i,131) - dpdr(i,315) -  dpdr(i,324)
    dpdr(i,329) = dpdr(i,22)*p(29) + p(22)*dpdr(i,29) - dpdr(i,313) -  dpdr(i,311) - dpdr(i,326) - dpdr(i,311)
    dpdr(i,330) = dpdr(i,11)*p(55) + p(11)*dpdr(i,55)  
    dpdr(i,331) = dpdr(i,22)*p(30) + p(22)*dpdr(i,30) - dpdr(i,313) -  dpdr(i,312) - dpdr(i,330) - dpdr(i,312)
    dpdr(i,332) = dpdr(i,2)*p(127) + p(2)*dpdr(i,127) - dpdr(i,309) -  dpdr(i,319)
    dpdr(i,333) = dpdr(i,6)*p(83)  + p(6)*dpdr(i,83) - dpdr(i,318)
    dpdr(i,334) = dpdr(i,7)*p(83)  + p(7)*dpdr(i,83) - dpdr(i,318)
    dpdr(i,335) = dpdr(i,11)*p(52) + p(11)*dpdr(i,52) - dpdr(i,315) -  dpdr(i,324)
    dpdr(i,336) = dpdr(i,2)*p(130) + p(2)*dpdr(i,130) - dpdr(i,313) -  dpdr(i,323) - dpdr(i,313)
    dpdr(i,337) = dpdr(i,2)*p(133) + p(2)*dpdr(i,133) - dpdr(i,316)
    dpdr(i,338) = dpdr(i,6)*p(77)  + p(6)*dpdr(i,77)  - dpdr(i,315) -  dpdr(i,309) - dpdr(i,322)
    dpdr(i,339) = dpdr(i,6)*p(78)  + p(6)*dpdr(i,78) - dpdr(i,312)
    dpdr(i,340) = dpdr(i,7)*p(77)  + p(7)*dpdr(i,77) - dpdr(i,315) -  dpdr(i,309) - dpdr(i,321)
    dpdr(i,341) = dpdr(i,7)*p(78)  + p(7)*dpdr(i,78) - dpdr(i,311)
    dpdr(i,342) = dpdr(i,11)*p(59) + p(11)*dpdr(i,59) - dpdr(i,329) -  dpdr(i,338)
    dpdr(i,343) = dpdr(i,11)*p(60) + p(11)*dpdr(i,60) - dpdr(i,331) -  dpdr(i,340)
    dpdr(i,344) = dpdr(i,3)*p(130) + p(3)*dpdr(i,130) - dpdr(i,315) -  dpdr(i,311) - dpdr(i,312) - dpdr(i,309) - dpdr(i,329) -  dpdr(i,338) - dpdr(i,331) - dpdr(i,340) - dpdr(i,328) -  dpdr(i,321) - dpdr(i,322) - dpdr(i,311) - dpdr(i,312)
    dpdr(i,345) = dpdr(i,18)*p(42) + p(18)*dpdr(i,42) - dpdr(i,310) -  dpdr(i,326) - dpdr(i,330) - dpdr(i,310)
    dpdr(i,346) = dpdr(i,2)*p(132) + p(2)*dpdr(i,132) - dpdr(i,315) -  dpdr(i,335)
    dpdr(i,347) = dpdr(i,6)*p(92)  + p(6)*dpdr(i,92) - dpdr(i,334)
    dpdr(i,348) = dpdr(i,7)*p(92)  + p(7)*dpdr(i,92) - dpdr(i,333)
    dpdr(i,349) = dpdr(i,3)*p(132) + p(3)*dpdr(i,132) - dpdr(i,316) -  dpdr(i,315) - dpdr(i,314) - dpdr(i,338) - dpdr(i,347) -  dpdr(i,340) - dpdr(i,348) - dpdr(i,336) - dpdr(i,333) -  dpdr(i,334) - dpdr(i,316) - dpdr(i,315) - dpdr(i,314) -  dpdr(i,316) - dpdr(i,314) - dpdr(i,316) - dpdr(i,314)
    dpdr(i,350) = dpdr(i,3)*p(133) + p(3)*dpdr(i,133) - dpdr(i,315) -  dpdr(i,339) - dpdr(i,341)
    dpdr(i,351) = dpdr(i,4)*p(132) + p(4)*dpdr(i,132) - dpdr(i,315) -  dpdr(i,344) - dpdr(i,342) - dpdr(i,343) - dpdr(i,332)
    dpdr(i,352) = dpdr(i,1)*p(222) + p(1)*dpdr(i,222) - dpdr(i,316) -  dpdr(i,350)
    dpdr(i,353) = dpdr(i,2)*p(134) + p(2)*dpdr(i,134) - dpdr(i,306)
    dpdr(i,354) = dpdr(i,2)*p(135) + p(2)*dpdr(i,135) - dpdr(i,323)
    dpdr(i,355) = dpdr(i,3)*p(134) + p(3)*dpdr(i,134) - dpdr(i,319)
    dpdr(i,356) = dpdr(i,11)*p(53) + p(11)*dpdr(i,53) - dpdr(i,310) -  dpdr(i,320)
    dpdr(i,357) = dpdr(i,2)*p(144) + p(2)*dpdr(i,144) - dpdr(i,329) -  dpdr(i,338)
    dpdr(i,358) = dpdr(i,2)*p(143) + p(2)*dpdr(i,143) - dpdr(i,324) -  dpdr(i,335)
    dpdr(i,359) = dpdr(i,6)*p(81)  + p(6)*dpdr(i,81) - dpdr(i,311) -  dpdr(i,324)
    dpdr(i,360) = dpdr(i,2)*p(146) + p(2)*dpdr(i,146) - dpdr(i,331) -  dpdr(i,340)
    dpdr(i,361) = dpdr(i,2)*p(150) + p(2)*dpdr(i,150) - dpdr(i,344)
    dpdr(i,362) = dpdr(i,7)*p(82)  + p(7)*dpdr(i,82) - dpdr(i,312) -  dpdr(i,324)
    dpdr(i,363) = dpdr(i,11)*p(65) + p(11)*dpdr(i,65) - dpdr(i,308)
    dpdr(i,364) = dpdr(i,2)*p(137) + p(2)*dpdr(i,137) - dpdr(i,318)
    dpdr(i,365) = dpdr(i,17)*p(45) + p(17)*dpdr(i,45) - dpdr(i,307)
    dpdr(i,366) = dpdr(i,4)*p(134) + p(4)*dpdr(i,134) - dpdr(i,317)
    dpdr(i,367) = dpdr(i,6)*p(96)  + p(6)*dpdr(i,96) - dpdr(i,332)
    dpdr(i,368) = dpdr(i,17)*p(54) + p(17)*dpdr(i,54)  
    dpdr(i,369) = dpdr(i,7)*p(96)  + p(7)*dpdr(i,96) - dpdr(i,332)
    dpdr(i,370) = dpdr(i,17)*p(55) + p(17)*dpdr(i,55)  
    dpdr(i,371) = dpdr(i,2)*p(159) + p(2)*dpdr(i,159) - dpdr(i,345) -  dpdr(i,349)
    dpdr(i,372) = dpdr(i,2)*p(147) + p(2)*dpdr(i,147) - dpdr(i,313)
    dpdr(i,373) = dpdr(i,11)*p(58) + p(11)*dpdr(i,58) - dpdr(i,315) -  dpdr(i,328)
    dpdr(i,374) = dpdr(i,2)*p(148) + p(2)*dpdr(i,148) - dpdr(i,329) -  dpdr(i,342)
    dpdr(i,375) = dpdr(i,6)*p(98)  + p(6)*dpdr(i,98) - dpdr(i,327)
    dpdr(i,376) = dpdr(i,2)*p(166) + p(2)*dpdr(i,166) - dpdr(i,359)
    dpdr(i,377) = dpdr(i,14)*p(54) + p(14)*dpdr(i,54) - dpdr(i,323)
    dpdr(i,378) = dpdr(i,2)*p(149) + p(2)*dpdr(i,149) - dpdr(i,331) -  dpdr(i,343)
    dpdr(i,379) = dpdr(i,7)*p(98)  + p(7)*dpdr(i,98) - dpdr(i,325)
    dpdr(i,380) = dpdr(i,7)*p(99)  + p(7)*dpdr(i,99) - dpdr(i,311) -  dpdr(i,359)
    dpdr(i,381) = dpdr(i,2)*p(167) + p(2)*dpdr(i,167) - dpdr(i,362)
    dpdr(i,382) = dpdr(i,14)*p(55) + p(14)*dpdr(i,55) - dpdr(i,323)
    dpdr(i,383) = dpdr(i,6)*p(100) + p(6)*dpdr(i,100) - dpdr(i,312) -  dpdr(i,362)
    dpdr(i,384) = dpdr(i,11)*p(66) + p(11)*dpdr(i,66) - dpdr(i,344)
    dpdr(i,385) = dpdr(i,6)*p(101) + p(6)*dpdr(i,101) - dpdr(i,325) -  dpdr(i,343) - dpdr(i,379)
    dpdr(i,386) = dpdr(i,7)*p(101) + p(7)*dpdr(i,101) - dpdr(i,342) -  dpdr(i,327) - dpdr(i,375)
    dpdr(i,387) = dpdr(i,6)*p(103) + p(6)*dpdr(i,103) - dpdr(i,313) -  dpdr(i,382)
    dpdr(i,388) = dpdr(i,11)*p(67) + p(11)*dpdr(i,67) - dpdr(i,314)
    dpdr(i,389) = dpdr(i,2)*p(152) + p(2)*dpdr(i,152) - dpdr(i,333)
    dpdr(i,390) = dpdr(i,2)*p(153) + p(2)*dpdr(i,153) - dpdr(i,334)
    dpdr(i,391) = dpdr(i,2)*p(154) + p(2)*dpdr(i,154) - dpdr(i,315) -  dpdr(i,373)
    dpdr(i,392) = dpdr(i,1)*p(240) + p(1)*dpdr(i,240) - dpdr(i,335) -  dpdr(i,373) - dpdr(i,366)
    dpdr(i,393) = dpdr(i,2)*p(155) + p(2)*dpdr(i,155) - dpdr(i,338) -  dpdr(i,342)
    dpdr(i,394) = dpdr(i,2)*p(171) + p(2)*dpdr(i,171) - dpdr(i,377)
    dpdr(i,395) = dpdr(i,2)*p(157) + p(2)*dpdr(i,157) - dpdr(i,340) -  dpdr(i,343)
    dpdr(i,396) = dpdr(i,2)*p(163) + p(2)*dpdr(i,163) - dpdr(i,350) -  dpdr(i,351)
    dpdr(i,397) = dpdr(i,6)*p(95)  + p(6)*dpdr(i,95) - dpdr(i,315) -  dpdr(i,350) - dpdr(i,340)
    dpdr(i,398) = dpdr(i,2)*p(172) + p(2)*dpdr(i,172) - dpdr(i,382)
    dpdr(i,399) = dpdr(i,7)*p(95)  + p(7)*dpdr(i,95) - dpdr(i,315) -  dpdr(i,350) - dpdr(i,338)
    dpdr(i,400) = dpdr(i,11)*p(68) + p(11)*dpdr(i,68) - dpdr(i,345) -  dpdr(i,349)
    dpdr(i,401) = dpdr(i,2)*p(175) + p(2)*dpdr(i,175) - dpdr(i,380) -  dpdr(i,385)
    dpdr(i,402) = dpdr(i,6)*p(106) + p(6)*dpdr(i,106) - dpdr(i,335) -  dpdr(i,343) - dpdr(i,396)
    dpdr(i,403) = dpdr(i,2)*p(176) + p(2)*dpdr(i,176) - dpdr(i,383) -  dpdr(i,386)
    dpdr(i,404) = dpdr(i,7)*p(106) + p(7)*dpdr(i,106) - dpdr(i,342) -  dpdr(i,335) - dpdr(i,396)
    dpdr(i,405) = dpdr(i,4)*p(150) + p(4)*dpdr(i,150) - dpdr(i,328) -  dpdr(i,359) - dpdr(i,362)
    dpdr(i,406) = dpdr(i,11)*p(69) + p(11)*dpdr(i,69) - dpdr(i,316)
    dpdr(i,407) = dpdr(i,2)*p(161) + p(2)*dpdr(i,161) - dpdr(i,347)
    dpdr(i,408) = dpdr(i,2)*p(162) + p(2)*dpdr(i,162) - dpdr(i,348)
    dpdr(i,409) = dpdr(i,11)*p(70) + p(11)*dpdr(i,70) - dpdr(i,350) -  dpdr(i,351)
    dpdr(i,410) = dpdr(i,2)*p(180) + p(2)*dpdr(i,180) - dpdr(i,397) -  dpdr(i,402)
    dpdr(i,411) = dpdr(i,6)*p(110) + p(6)*dpdr(i,110) - dpdr(i,348)
    dpdr(i,412) = dpdr(i,2)*p(181) + p(2)*dpdr(i,181) - dpdr(i,399) -  dpdr(i,404)
    dpdr(i,413) = dpdr(i,7)*p(110) + p(7)*dpdr(i,110) - dpdr(i,347)
    dpdr(i,414) = dpdr(i,6)*p(112) + p(6)*dpdr(i,112) - dpdr(i,316) -  dpdr(i,398) - dpdr(i,316)
    dpdr(i,415) = dpdr(i,11)*p(71) + p(11)*dpdr(i,71) - dpdr(i,352)
    dpdr(i,416) = dpdr(i,2)*p(184) + p(2)*dpdr(i,184) - dpdr(i,411)
    dpdr(i,417) = dpdr(i,2)*p(185) + p(2)*dpdr(i,185) - dpdr(i,413)
    dpdr(i,418) = dpdr(i,12)*p(71) + p(12)*dpdr(i,71) - dpdr(i,351) -  dpdr(i,416) - dpdr(i,417)
    dpdr(i,419) = dpdr(i,2)*p(164) + p(2)*dpdr(i,164) - dpdr(i,320) -  dpdr(i,356)
    dpdr(i,420) = dpdr(i,2)*p(165) + p(2)*dpdr(i,165) - dpdr(i,328) -  dpdr(i,373)
    dpdr(i,421) = dpdr(i,2)*p(170) + p(2)*dpdr(i,170) - dpdr(i,337) -  dpdr(i,392)
    dpdr(i,422) = dpdr(i,1)*p(268) + p(1)*dpdr(i,268) - dpdr(i,359)
    dpdr(i,423) = dpdr(i,1)*p(269) + p(1)*dpdr(i,269) - dpdr(i,362)
    dpdr(i,424) = dpdr(i,2)*p(174) + p(2)*dpdr(i,174) - dpdr(i,345) -  dpdr(i,400)
    dpdr(i,425) = dpdr(i,6)*p(102) + p(6)*dpdr(i,102) - dpdr(i,345) -  dpdr(i,326)
    dpdr(i,426) = dpdr(i,7)*p(103) + p(7)*dpdr(i,103) - dpdr(i,345) -  dpdr(i,330)
    dpdr(i,427) = dpdr(i,3)*p(186) + p(3)*dpdr(i,186) - dpdr(i,364)
    dpdr(i,428) = dpdr(i,18)*p(65) + p(18)*dpdr(i,65) - dpdr(i,354)
    dpdr(i,429) = dpdr(i,17)*p(66) + p(17)*dpdr(i,66) - dpdr(i,361)
    dpdr(i,430) = dpdr(i,2)*p(179) + p(2)*dpdr(i,179) - dpdr(i,350) -  dpdr(i,409)
    dpdr(i,431) = dpdr(i,4)*p(166) + p(4)*dpdr(i,166) - dpdr(i,361) -  dpdr(i,357) - dpdr(i,422)
    dpdr(i,432) = dpdr(i,4)*p(167) + p(4)*dpdr(i,167) - dpdr(i,361) -  dpdr(i,360) - dpdr(i,423)
    dpdr(i,433) = dpdr(i,2)*p(187) + p(2)*dpdr(i,187) - dpdr(i,387)
    dpdr(i,434) = dpdr(i,14)*p(66) + p(14)*dpdr(i,66) - dpdr(i,328) -  dpdr(i,405) - dpdr(i,376) - dpdr(i,381) - dpdr(i,373)
    dpdr(i,435) = dpdr(i,1)*p(277) + p(1)*dpdr(i,277) - dpdr(i,387) -  dpdr(i,385) - dpdr(i,425)
    dpdr(i,436) = dpdr(i,1)*p(278) + p(1)*dpdr(i,278) - dpdr(i,387) -  dpdr(i,386) - dpdr(i,426)
    dpdr(i,437) = dpdr(i,2)*p(177) + p(2)*dpdr(i,177) - dpdr(i,346) -  dpdr(i,391)
    dpdr(i,438) = dpdr(i,2)*p(178) + p(2)*dpdr(i,178) - dpdr(i,349) -  dpdr(i,400)
    dpdr(i,439) = dpdr(i,1)*p(281) + p(1)*dpdr(i,281) - dpdr(i,396) -  dpdr(i,392) - dpdr(i,430)
    dpdr(i,440) = dpdr(i,21)*p(54) + p(21)*dpdr(i,54) - dpdr(i,370)
    dpdr(i,441) = dpdr(i,21)*p(55) + p(21)*dpdr(i,55) - dpdr(i,368)
    dpdr(i,442) = dpdr(i,2)*p(189) + p(2)*dpdr(i,189) - dpdr(i,405) -  dpdr(i,434)
    dpdr(i,443) = dpdr(i,1)*p(285) + p(1)*dpdr(i,285) - dpdr(i,402) -  dpdr(i,404) - dpdr(i,400) - dpdr(i,434)  - dpdr(i,430)
    dpdr(i,444) = dpdr(i,6)*p(116) + p(6)*dpdr(i,116) - dpdr(i,347) -  dpdr(i,413) - dpdr(i,403)
    dpdr(i,445) = dpdr(i,7)*p(116) + p(7)*dpdr(i,116) - dpdr(i,348) -  dpdr(i,411) - dpdr(i,401)
    dpdr(i,446) = dpdr(i,2)*p(182) + p(2)*dpdr(i,182) - dpdr(i,351) -  dpdr(i,409)
    dpdr(i,447) = dpdr(i,2)*p(191) + p(2)*dpdr(i,191) - dpdr(i,414) -  dpdr(i,443)
    dpdr(i,448) = dpdr(i,3)*p(183) + p(3)*dpdr(i,183) - dpdr(i,351) -  dpdr(i,418) - dpdr(i,411) - dpdr(i,413) - dpdr(i,409)
    dpdr(i,449) = dpdr(i,6)*p(118) + p(6)*dpdr(i,118) - dpdr(i,351) -  dpdr(i,418) - dpdr(i,412)
    dpdr(i,450) = dpdr(i,7)*p(118) + p(7)*dpdr(i,118) - dpdr(i,351) -  dpdr(i,418) - dpdr(i,410)
    dpdr(i,451) = dpdr(i,2)*p(193) + p(2)*dpdr(i,193) - dpdr(i,418) -  dpdr(i,448)
    dpdr(i,452) = dpdr(i,1)*p(294) + p(1)*dpdr(i,294) - dpdr(i,418) -  dpdr(i,415) - dpdr(i,448)
    dpdr(i,453) = dpdr(i,6)*p(119) + p(6)*dpdr(i,119) - dpdr(i,417)
    dpdr(i,454) = dpdr(i,7)*p(119) + p(7)*dpdr(i,119) - dpdr(i,416)
    dpdr(i,455) = dpdr(i,2)*p(186) + p(2)*dpdr(i,186) - dpdr(i,363)
    dpdr(i,456) = dpdr(i,3)*p(187) + p(3)*dpdr(i,187) - dpdr(i,385) -  dpdr(i,386) - dpdr(i,384) - dpdr(i,435) - dpdr(i,436) -  dpdr(i,434)
    dpdr(i,457) = dpdr(i,2)*p(188) + p(2)*dpdr(i,188) - dpdr(i,388)
    dpdr(i,458) = dpdr(i,4)*p(187) + p(4)*dpdr(i,187) - dpdr(i,425) -  dpdr(i,426) - dpdr(i,424)
    dpdr(i,459) = dpdr(i,2)*p(190) + p(2)*dpdr(i,190) - dpdr(i,406)
    dpdr(i,460) = dpdr(i,21)*p(66) + p(21)*dpdr(i,66) - dpdr(i,422) -  dpdr(i,423) - dpdr(i,420)
    dpdr(i,461) = dpdr(i,2)*p(192) + p(2)*dpdr(i,192) - dpdr(i,415)
    dpdr(i,462) = dpdr(i,18)*p(71) + p(18)*dpdr(i,71) - dpdr(i,440) -  dpdr(i,441) - dpdr(i,438)
    dpdr(i,463) = dpdr(i,2)*p(194) + p(2)*dpdr(i,194) - dpdr(i,452)
    dpdr(i,464) = dpdr(i,3)*p(194) + p(3)*dpdr(i,194) - dpdr(i,453) -  dpdr(i,454) - dpdr(i,451)
    dpdr(i,465) = dpdr(i,1)*p(305) + p(1)*dpdr(i,305) - dpdr(i,464) -  dpdr(i,463)
  enddo

end subroutine EvdPdR
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine evdbdr
!**********************************************************************
!  The subroutine elminates the 2-body terms in Bowman's approach
!**********************************************************************

  integer                     :: i,j
  real(rkp) ,dimension(6,466) :: db1dr

  ! Pass P(0:465) to BM1(1:466)
  do j=1,6
    do i=1,466
      db1dr(j,i) = dpdr(j,i-1)
    enddo
  enddo

  ! Remove unconnected terms and 2-body terms and pass to B(1:430)
  do j=1,6
    dbdr(j,1)=db1dr(j,4)

    do i=2,4
      dbdr(j,i)=db1dr(j,i+4)
    enddo

    dbdr(j,5)=db1dr(j,10)

    do i=6,11
      dbdr(j,i)=db1dr(j,i+6)
    enddo

    dbdr(j,12)=db1dr(j,19)
    dbdr(j,13)=db1dr(j,21)

    do i=14,26
      dbdr(j,i)=db1dr(j,i+9)
    enddo

    dbdr(j,27)=db1dr(j,37)
    dbdr(j,28)=db1dr(j,39)

    do i=29,53
      dbdr(j,i)=db1dr(j,i+12)
    enddo

    dbdr(j,54)=db1dr(j,67)
    dbdr(j,55)=db1dr(j,69)
    dbdr(j,56)=db1dr(j,71)

    do i=57,97
      dbdr(j,i)=db1dr(j,i+16)
    enddo

    dbdr(j,98)=db1dr(j,115)
    dbdr(j,99)=db1dr(j,117)
    dbdr(j,100)=db1dr(j,119)
    
    do i=101,166
      dbdr(j,i)=db1dr(j,i+20)
    enddo

    dbdr(j,167)=db1dr(j,188)
    dbdr(j,168)=db1dr(j,190)
    dbdr(j,169)=db1dr(j,192)
    dbdr(j,170)=db1dr(j,194)

    do i=171,272
      dbdr(j,i)=db1dr(j,i+25)
    enddo

    dbdr(j,273)=db1dr(j,299)
    dbdr(j,274)=db1dr(j,301)
    dbdr(j,275)=db1dr(j,303)
    dbdr(j,276)=db1dr(j,305)

    do i=277,425
      dbdr(j,i)=db1dr(j,i+30)
    enddo

    dbdr(j,426)=db1dr(j,457)
    dbdr(j,427)=db1dr(j,459)
    dbdr(j,428)=db1dr(j,461)
    dbdr(j,429)=db1dr(j,463)
    dbdr(j,430)=db1dr(j,465) 
  enddo

end subroutine      
!--------------------------------------------------------------------------------------------------------------------------------!

End Module