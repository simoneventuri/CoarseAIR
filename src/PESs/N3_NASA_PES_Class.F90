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

Module N3_NASA_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    N3_NASA_PES_Type


  Type    ,extends(PES_Type)          ::    N3_NASA_PES_Type
  contains
    procedure         ::  Initialize     =>    Initialize_N3_NASA_PES
    procedure         ::  Output         =>    Output_N3_NASA_PES
    procedure         ::  Compute        =>    Compute_N3_NASA_PES_1d
    procedure         ::  Potential      =>    N3_NASA_Potential_From_R
    procedure         ::  TriatPotential =>    N3_NASA_Potential_From_R_OnlyTriat
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.

  integer                         ,parameter    ::    mp    =   7
  integer   ,dimension(14)        ,parameter    ::    np    =   [1,5,11,21,34,52,74,102,135,175,221,275,336,406]
  real(rkp)                       ,parameter    ::    b1    =   0.203751689792889D+00
  real(rkp)                       ,parameter    ::    b2    =   0.275065598827013D+00
  real(rkp) ,dimension(52)        ,parameter    ::    &
    csr   =   [ -4.88634974677321D+01,-1.76620388625638D+01, 1.18417470793576D+02,-1.86954134323098D+02,&
                 3.14771520925086D+02, 3.81831837270977D+01,-7.78124866957063D+01, 1.74537233441170D+02,&
                -4.33954446347525D+01,-5.04137033787021D+00,-3.30536046284311D+02,-1.58750076501267D+01,&
                 1.02810708092583D+01,-6.09335430973599D+01, 2.38071318440977D+01,-9.65143583479881D+01,&
                 1.62705073558967D+01, 9.70619283318971D+01, 1.34123205900499D+01, 1.25149282644730D+02,&
                -2.01680569695858D+01, 2.81313764400686D+00,-2.80603430286581D+00, 9.35730949310909D+00,&
                -3.21488975988290D+00, 3.11387276395095D+01, 1.96468734967068D+00,-3.15886921235925D+00,&
                -6.76716567398226D+00,-4.70199700309722D+00,-2.95076511748160D+01, 1.51820045120461D+01,&
                -2.08883534727594D+01, 3.61031308372215D+00,-1.91167601539005D-01, 2.74109937832430D-01,&
                -5.27297920356713D-01,-1.70481123134168D-01,-2.93916599319425D+00, 4.74527566713238D-01,&
                 1.02694890552781D+01, 8.27033415839729D-01,-4.10215809838363D-01,-9.91773175105535D+00,&
                -1.96256961047575D+00, 4.34270329548511D-01, 1.84504824582918D+00,-3.43477508526355D+00,&
                 5.56384160845247D-02, 1.45753977547667D+00, 1.35243114380008D-01,-3.16034746986889D-01 ]
  integer   ,dimension(52)          ,parameter  ::    indx1     =   [1,2,1,1,0,3,2,2,1,1,0,4,3,3,2,2,1,1,1,0,0,5,4,4,3,3,2,2,2,1,1,1,0,0,6,5,5,4,4,3,3,3,2,2,2,1,1,1,1,0,0,0]
  integer   ,dimension(52)          ,parameter  ::    indx2     =   [1,1,2,0,1,1,2,0,3,1,2,1,2,0,3,1,4,2,0,3,1,1,2,0,3,1,4,2,0,5,3,1,4,2,1,2,0,3,1,4,2,0,5,3,1,6,4,2,0,5,3,1]
  integer   ,dimension(52)          ,parameter  ::    indx3     =   [0,0,0,2,2,0,0,2,0,2,2,0,0,2,0,2,0,2,4,2,4,0,0,2,0,2,0,2,4,0,2,4,2,4,0,0,2,0,2,0,2,4,0,2,4,0,2,4,6,2,4,6]
  real(rkp)                         ,parameter  ::    adisp     =   2.87E0_rkp
  real(rkp)                         ,parameter  ::    bdisp     =   1.40E0_rkp
  real(rkp)                         ,parameter  ::    rdisp     =   3.15E0_rkp
  real(rkp)                         ,parameter  ::    c62       =   5.41E0_rkp
  real(rkp)                         ,parameter  ::    dc        =   5.5E0_rkp
  real(rkp)                         ,parameter  ::    dc6       =   dc**6
  real(rkp)                         ,parameter  ::    Tolerence =   1.0E-14_rkp
  real(rkp)                         ,parameter  ::    RaMax     =   -132.0_rkp
  real(rkp)                         ,parameter  ::    ebn2      =   0.3554625704d0
  
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NASA PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_N3_NASA_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(N3_NASA_PES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms  
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP
  character(*)                    ,parameter                ::    Name_PES = 'N3_NASA'
  integer         ,dimension(3,2)                           ::    iA
  type(DiatomicPotential_Factory_Type)                       ::    DiatPotFactory
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N3_NASA_PES" )
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



Subroutine Output_N3_NASA_PES( This, Unit )

  class(N3_NASA_PES_Type)                 ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('N3 PES with cta bug fixed')")
  write(Unit,"('acpf energies, 1s2s core and leroy n2 with hyper radius repulsion')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N3_NASA_Potential_From_R( This, R, Q ) result( V )

  class(N3_NASA_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  integer                                                    ::    nrep
  integer                                                    ::    i, k, iP
  real(rkp)                                                  ::    vrp
  real(rkp) ,dimension(406)                                  ::    tp
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp) ,dimension(3)                                    ::    dVDiat
  real(rkp)                                                  ::    t1, t2
  
  ! Computing the diatomic potential energies associated to the 3 internuclear distances
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do                  

  nrep        =   np(mp-1)

  call tapsrd_NoGrad( mp, R(1), R(2), R(3), tp, vrp )
  do i = 1,nrep
    vrp = vrp + csr(i) * tp(i)
  end do
  V = sum(VDiat(:)) + vrp + ebn2                               ! Terms 2 and 3 are the many-body interaction potential
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N3_NASA_Potential_From_R_OnlyTriat( This, R, Q ) result( V )

  class(N3_NASA_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  integer                                                 ::    nrep
  integer                                                 ::    i, k
  real(rkp)                                               ::    vrp
  real(rkp) ,dimension(406)                               ::    tp
  real(rkp) ,dimension(size(R))                           ::    Vd          ! Diatomic potential energies associated to distances of pairs. Dim=(NPairs,NTraj)
  real(rkp)                                               ::    t1, t2
  
  nrep        =   np(mp-1)

  call tapsrd_NoGrad( mp, R(1), R(2), R(3), tp, vrp )
  do i = 1,nrep
    vrp = vrp + csr(i) * tp(i)
  end do
  V = vrp 
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_N3_NASA_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(N3_NASA_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  integer                                                 ::    nrep
  integer                                                 ::    i, iP
  real(rkp)                                               ::    vrp
  real(rkp) ,dimension(3)                                 ::    dlr
  real(rkp) ,dimension(406)                               ::    tp
  real(rkp) ,dimension(3,52)                              ::    tpd
  real(rkp) ,dimension(3)                                 ::    Vd          ! Diatomic potential energies associated to the pair distances. Dim=(NPairs,NTraj)
  real(rkp)                                               ::    t1, t2

  !call cpu_time ( t1 )

  dVdQ         = Zero

  ! Computing the diatomic potential energies associated to the 3 internuclear distances
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), Vd(iP), dVdR(iP) )
  end do     

  nrep        =   np(mp-1)
  
  !dV = 0.d0
  call tapsrd( mp, R(1), R(2), R(3), tp, tpd, vrp, dlr )
  dVdR(:)   =   dVdR(:) + dlr
  do i = 1,nrep
    vrp   =   vrp + csr(i) * tp(i)
    dVdR(:) =   dVdR(:)  + csr(i) * tpd(:,i)
  end do
  V = vrp + sum(Vd(:)) + ebn2
                                                                           ! Terms 2 and 3 are the many-body interaction potential
  !call cpu_time ( t2 )
  !write(*,*) 'Time for Potential Calculations = ', t2-t1
                                                                           
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
! **************************************************************************************************************
!                                   PRIVATE PROCEDURES for NASA PES
! **************************************************************************************************************
! **************************************************************************************************************
Pure Subroutine tapsrd( mpmax, Ra, Rb, Rc, tp, tpd, vlr, dlr )
! This procedure determine the terms for an analytic fit to the potential for a molecule a3 from
! the internuclear separation distances Ra,Rb,Rc.
! The input variables are:
!  * kder = 0: fit a polynomial to powers of R,R,cos(\theta )
!         = 1:  ''      ''      ''   ''   '' rho, ...
!  * mpmax: maximum power of M < or = 14

  use Parameters_Module       ,only:  Zero, One, Half, Two, Quart

  integer                       ,intent(in)     ::    mpmax   ! maximum power of M < or = 14
  real(rkp)                     ,intent(in)     ::    Ra      ! Distance between nuclear centers [bohr]
  real(rkp)                     ,intent(in)     ::    Rb      ! Distance between nuclear centers [bohr]
  real(rkp)                     ,intent(in)     ::    Rc      ! Distance between nuclear centers [bohr]
  real(rkp) ,dimension(406)     ,intent(out)    ::    tp
  real(rkp) ,dimension(3,52)    ,intent(out)    ::    tpd
  real(rkp)                     ,intent(out)    ::    vlr
  real(rkp) ,dimension(3)       ,intent(out)    ::    dlr

  logical                                       ::    Compute
  integer                                       ::    n
  real(rkp)                                     ::    Ra2, Rb2, Rc2
  real(rkp)                                     ::    delta, rho
  real(rkp)                                     ::    cab, cbc, cca
  real(rkp)                                     ::    rsa, rsb, rsc
  real(rkp)                                     ::    rsa2, rsb2, rsc2
  real(rkp)                                     ::    cta, ctb, ctc
  real(rkp)                                     ::    bot, term0, term, damp, plm
  real(rkp) ,dimension(8,3)                     ::    dpl
  real(rkp) ,dimension(3,3)                     ::    dep
  real(rkp) ,dimension(3)                       ::    R
  real(rkp) ,dimension(3)                       ::    SumR2
  real(rkp) ,dimension(3)                       ::    SumRs2
  real(rkp) ,dimension(3)                       ::    dppl
  real(rkp) ,dimension(51)                      ::    pla, plb, plc
  real(rkp) ,dimension(3)                       ::    rs, ep, pl
  real(rkp) ,dimension(3)                       ::    dsa, dsb, dsc
  real(rkp) ,dimension(3,3)                     ::    drs, dct

  Compute     =   .True.
  R(1:3)      =   [Ra,  Rb,  Rc  ]
  Ra2         =   Ra**2
  Rb2         =   Rb**2
  Rc2         =   Rc**2
  delta       =   ( Ra + Rb - Rc ) / Rc

!     right triangle
  if ( abs(delta) < Tolerence ) then
!     if ( Ra > RaMax ) call Error( "[tapsrd]: Ra > RaMax" )
    cab       =   -One
    cbc       =    One
    cca       =    One
    rsa       =   Ra * Half + Rb
    rsb       =   Rb * Half + Ra
    rsc       =   abs(Rb-Ra) * Half
    dsa(1:3)  =   [ Half , One  , Zero ]
    dsb(1:3)  =   [ One  , Half , Zero ]
    if ( Rb > Ra ) then; dsc(1:3) = [-Half , Half , Zero ]
    else;                 dsc(1:3) = [ Half ,-Half , Zero ]
    end if
    rsa2      =   rsa * rsa
    rsb2      =   rsb * rsb
    rsc2      =   rsc * rsc
    cta       =    One
    ctb       =   -One
    ctc       =    One
    Compute   =   .False.
  end if

!     other triangle
  if (Compute) then

    bot       =   Half / (Ra*Rb)
    cab       =   ( Ra2 + Rb2 - Rc2 ) * bot

    bot       =   Half / (Rb*Rc)
    cbc       =   ( Rb2 + Rc2 - Ra2 ) * bot

    bot       =   Half / (Rc*Ra)
    cca       =   ( Rc2 + Ra2 - Rb2 ) * bot

    rsa2      =   Ra2 * Quart + Rb2 - Ra * Rb * cab
    rsb2      =   Rb2 * Quart + Rc2 - Rb * Rc * cbc
    rsc2      =   Rc2 * Quart + Ra2 - Rc * Ra * cca

    rsa       =   sqrt(rsa2)
    bot       =   Half / rsa
    drs(1,1)  =   -Half*Ra*bot
    drs(1,2)  =   Rb*bot
    drs(1,3)  =   Rc*bot

    rsb       =   sqrt(rsb2)
    bot       =   Half / rsb
    drs(2,1)  =   Ra*bot
    drs(2,2)  =   -Half*Rb*bot
    drs(2,3)  =   Rc*bot

    rsc       =   sqrt(rsc2)
    bot       =   Half / rsc
    drs(3,1)  =   Ra*bot
    drs(3,2)  =   Rb*bot
    drs(3,3)  =   -Half*Rc*bot

    bot       =   One/(Ra*rsa)
    cta       =   (Ra2/4.0d0+rsa2-Rb2)*bot
    dct(1,1)  =   (Half*Ra+Two*rsa*drs(1,1)-cta*(rsa+Ra*drs(1,1)))*bot
    dct(1,2)  =   (-Two*Rb+Two*rsa*drs(1,2)-cta*(Ra*drs(1,2)))*bot
    dct(1,3)  =   (Two*rsa*drs(1,3)-cta*(Ra*drs(1,3)))*bot

    bot       =   One/(Rb*rsb)
    ctb       =   (Rb2/4.0d0+rsb2-Rc2)*bot
    dct(2,1)  =   (Two*rsb*drs(2,1)-ctb*(Rb*drs(2,1)))*bot
    dct(2,2)  =   (Half*Rb+Two*rsb*drs(2,2)-ctb*(rsb+Rb*drs(2,2)))*bot
    dct(2,3)  =   (-Two*Rc+Two*rsb*drs(2,3)-ctb*(Rb*drs(2,3)))*bot

    bot       =   One/(Rc*rsc)
    ctc       =   (Rc2/4.0d0+rsc2-Ra2)*bot
    dct(3,1)  =   (-Two*Ra+Two*rsc*drs(3,1)-ctc*(Rc*drs(3,1)))*bot
    dct(3,2)  =   (Two*rsc*drs(3,2)-ctc*(Rc*drs(3,2)))*bot
    dct(3,3)  =   (Half*Rc+Two*rsc*drs(3,3)-ctc*(rsc+Rc*drs(3,3)))*bot

  end if

!   write(Logger%Unit,"(6x,'[tapsrd]: Calling plnd')")
  call plnd( mpmax, cta, pla, dpl(1,1) )
  call plnd( mpmax, ctb, plb, dpl(1,2) )
  call plnd( mpmax, ctc, plc, dpl(1,3) )
!   write(Logger%Unit,"(6x,'[tapsrd]: done plnd')")


  rs(1:3)     =   [rsa, rsb, rsc ]
  SumR2(1:3)  =   [Ra2, Rb2, Rc2 ]
  SumRs2(1:3) =   [rsa2,rsb2,rsc2]

  rho         =   sqrt( (Two*rsa2/3d0) + Half*Ra2 )
  vlr         =   exp( - 6d0 * (rho-1.85d0) )

  dlr         =   Zero

  do n = 1,3
    ep(n)     =   exp(-b1*SumR2(n)-b2*SumRs2(n))
    damp      =   adisp*exp(-bdisp*((R(n)-rdisp)**2))
    if      ( n == 1 ) then;  plm = pla(3)
    else if ( n == 2 ) then;  plm = plb(3)
    else;                     plm = plc(3)
    end if
    term0     =   - c62 * plm * damp / ( rs(n)**6 + dc6 )
    vlr       =   vlr + term0
    term      =   - dpl(3,n) * c62 * damp / ( rs(n)**6 + dc6 )
    dlr(1:3)  =   dlr(1:3) + term*dct(n,1:3)
    dlr(n)    =   dlr(n) - Two * term0 * bdisp * ( R(n) - rdisp )
    term      =   - term0 * 6d0 * rs(n)**5 / ( rs(n)**6 + dc6 )
    dlr(1:3)  =   dlr(1:3) + term * drs(n,1:3)
  end do

!   write(Logger%Unit,"(6x,'[tapsrd]: flag 1')")
  dep(1,1)    =   ep(1) * (-b1*Two*Ra-b2*Two*rsa*drs(1,1))
  dep(1,2)    =   ep(1) * (-b2*Two*rsa*drs(1,2))
  dep(1,3)    =   ep(1) * (-b2*Two*rsa*drs(1,3))
  dep(2,1)    =   ep(2) * (-b2*Two*rsb*drs(2,1))
  dep(2,2)    =   ep(2) * (-b1*Two*Rb-b2*Two*rsb*drs(2,2))
  dep(2,3)    =   ep(2) * (-b2*Two*rsb*drs(2,3))
  dep(3,1)    =   ep(3) * (-b2*Two*rsc*drs(3,1))
  dep(3,2)    =   ep(3) * (-b2*Two*rsc*drs(3,2))
  dep(3,3)    =   ep(3) * (-b1*Two*Rc-b2*Two*rsc*drs(3,3))

!   write(Logger%Unit,"(6x,'[tapsrd]: flag 2')")
  do n = 1,52
    pl(1)     =   pla( indx3(n)+1 )
    pl(2)     =   plb( indx3(n)+1 )
    pl(3)     =   plc( indx3(n)+1 )
    dppl(1)   =   dpl( indx3(n)+1,1)
    dppl(2)   =   dpl( indx3(n)+1,2)
    dppl(3)   =   dpl( indx3(n)+1,3)
!     write(Logger%Unit,"(6x,'[tapsrd]: calling ptud for n= ',g0)") n
    call ptud( indx1(n), indx2(n), R, rs, ep, pl, tp(n), dep, drs, tpd(1,n), dppl, dct )
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine ptud( ns, ms, R, rs, ep, pl, tps, dep, drs, dr, dpl, dct )
! This procedure determines the terms for an analytic fit for a molecule A3 using Jacobian coordinates

  use Parameters_Module       ,only:  Zero

  integer                                   ,intent(in)     ::    ns
  integer                                   ,intent(in)     ::    ms
  real(rkp) ,dimension(3)                   ,intent(in)     ::    R
  real(rkp) ,dimension(3)                   ,intent(in)     ::    rs
  real(rkp) ,dimension(3)                   ,intent(in)     ::    ep
  real(rkp) ,dimension(3)                   ,intent(in)     ::    pl
  real(rkp)                                 ,intent(out)    ::    tps
  real(rkp) ,dimension(3,3)                 ,intent(in)     ::    dep
  real(rkp) ,dimension(3,3)                 ,intent(in)     ::    drs
  real(rkp) ,dimension(3)                   ,intent(out)    ::    dr
  real(rkp) ,dimension(3)                   ,intent(in)     ::    dpl
  real(rkp) ,dimension(3,3)                 ,intent(in)     ::    dct

  integer                                                   ::    n
  real(rkp)                                                 ::    dd, dtp
  real(rkp) ,dimension(3)                                   ::    tp

  tps       =   Zero
  dr        =   Zero

  do n = 1,3

    tp(n)   =   R(n)**ns * rs(n)**ms * pl(n)
    dtp     =   tp(n) * ep(n)
    dr      =   dr + tp(n) * dep(n,:)
    tps     =   tps + dtp

    if ( ns /= 0 ) dr(n) = dr(n) + real(ns,8) * R(n)**(ns-1) * rs(n)**ms * pl(n) * ep(n)

    if ( ms /= 0 ) then
      dd    =   real(ms,8) * R(n)**ns * rs(n)**(ms-1) * pl(n) * ep(n)
      dr    =   dr + dd * drs(n,:)
    end if
    dd      =   R(n)**ns * rs(n)**ms * ep(n) * dpl(n)
    dr      =   dr + dd * dct(n,:)

  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine plnd( nf, x, p, dp )

  use Parameters_Module       ,only:  Zero, One, Two

  integer                                   ,intent(in)     ::    nf
  real(rkp)                                 ,intent(in)     ::    x
  real(rkp) ,dimension(51)                  ,intent(inout)  ::    p
  real(rkp) ,dimension(8)                   ,intent(inout)  ::    dp

  integer                                                   ::    i
  real(rkp)                                                 ::    fni

  p(1:2)    =   [ One , x  ]
  dp(1:2)   =   [ Zero, One]

  do i = 2,nf
    fni     =   One / real(i,rkp)
    p(i+1)  =   Two * x * p(i) - p(i-1)    - ( (x*p(i)-p(i-1)) * fni )
    dp(i+1) =   Two * ( p(i) + x * dp(i) ) - dp(i-1) - ( (p(i) + x * dp(i) - dp(i-1)) * fni )
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine tapsrd_NoGrad( mpmax, Ra, Rb, Rc, tp, vlr )

  use Parameters_Module       ,only:  Zero, One, Half, Two, Quart

  integer                       ,intent(in)     ::    mpmax   ! maximum power of M < or = 14
  real(rkp)                     ,intent(in)     ::    Ra      ! Distance between nuclear centers [bohr]
  real(rkp)                     ,intent(in)     ::    Rb      ! Distance between nuclear centers [bohr]
  real(rkp)                     ,intent(in)     ::    Rc      ! Distance between nuclear centers [bohr]
  real(rkp) ,dimension(406)     ,intent(out)    ::    tp
  real(rkp)                     ,intent(out)    ::    vlr

  logical                                       ::    Compute
  integer                                       ::    n
  real(rkp)                                     ::    Ra2, Rb2, Rc2
  real(rkp)                                     ::    delta, rho
  real(rkp)                                     ::    cab, cbc, cca
  real(rkp)                                     ::    rsa, rsb, rsc
  real(rkp)                                     ::    rsa2, rsb2, rsc2
  real(rkp)                                     ::    cta, ctb, ctc
  real(rkp)                                     ::    bot, term0, damp, plm
  real(rkp) ,dimension(3)                       ::    R
!   real(rkp)                                     ::    Ra      ! Distances between nuclear centers [bohr]
!   real(rkp)                                     ::    Rb      ! Distances between nuclear centers [bohr]
!   real(rkp)                                     ::    Rc      ! Distances between nuclear centers [bohr]
  real(rkp) ,dimension(8,3)                     ::    dpl
  real(rkp) ,dimension(3)                       ::    SumR2
  real(rkp) ,dimension(3)                       ::    SumRs2
  real(rkp) ,dimension(51)                      ::    pla, plb, plc
  real(rkp) ,dimension(3)                       ::    rs, ep, pl
  real(rkp) ,dimension(3)                       ::    dsa, dsb, dsc

  Compute     =   .True.
!   Ra          =     R(1)
!   Rb          =     R(2)
!   Rc          =     R(3)
  R(1:3)      =   [Ra,  Rb,  Rc  ]
  Ra2         =   Ra**2
  Rb2         =   Rb**2
  Rc2         =   Rc**2
  delta       =   ( Ra + Rb - Rc ) / Rc

!     right triangle
  if ( abs(delta) < Tolerence ) then
!     if ( Ra > RaMax ) call Error( "[tapsrd]: Ra > RaMax" )
    cab       =   -One
    cbc       =    One
    cca       =    One
    rsa       =   Ra * Half + Rb
    rsb       =   Rb * Half + Ra
    rsc       =   abs(Rb-Ra) * Half
    dsa(1:3)  =   [ Half , One  , Zero ]
    dsb(1:3)  =   [ One  , Half , Zero ]
    if ( Rb > Ra ) then; dsc(1:3) = [-Half , Half , Zero ]
    else;                 dsc(1:3) = [ Half ,-Half , Zero ]
    end if
    rsa2      =   rsa * rsa
    rsb2      =   rsb * rsb
    rsc2      =   rsc * rsc
    cta       =    One
    ctb       =   -One
    ctc       =    One
    Compute   =   .False.
  end if

!     other triangle
  if (Compute) then

    bot       =   Half / (Ra*Rb)
    cab       =   ( Ra2 + Rb2 - Rc2 ) * bot

    bot       =   Half / (Rb*Rc)
    cbc       =   ( Rb2 + Rc2 - Ra2 ) * bot

    bot       =   Half / (Rc*Ra)
    cca       =   ( Rc2 + Ra2 - Rb2 ) * bot

    rsa2      =   Ra2 * Quart + Rb2 - Ra * Rb * cab
    rsb2      =   Rb2 * Quart + Rc2 - Rb * Rc * cbc
    rsc2      =   Rc2 * Quart + Ra2 - Rc * Ra * cca

    rsa       =   sqrt(rsa2)
    bot       =   Half / rsa

    rsb       =   sqrt(rsb2)
    bot       =   Half / rsb

    rsc       =   sqrt(rsc2)
    bot       =   Half / rsc

    bot       =   One/(Ra*rsa)
    cta       =   (Ra2/4.0d0+rsa2-Rb2)*bot

    bot       =   One/(Rb*rsb)
    ctb       =   (Rb2/4.0d0+rsb2-Rc2)*bot

    bot       =   One/(Rc*rsc)
    ctc       =   (Rc2/4.0d0+rsc2-Ra2)*bot
  end if

!   write(Logger%Unit,"(6x,'[tapsrd]: Calling plnd')")
  call plnd( mpmax, cta, pla, dpl(1,1) )
  call plnd( mpmax, ctb, plb, dpl(1,2) )
  call plnd( mpmax, ctc, plc, dpl(1,3) )
!   write(Logger%Unit,"(6x,'[tapsrd]: done plnd')")


!   R(1:3)      =   [Ra,  Rb,  Rc  ]
  rs(1:3)     =   [rsa, rsb, rsc ]
  SumR2(1:3)  =   [Ra2, Rb2, Rc2 ]
  SumRs2(1:3) =   [rsa2,rsb2,rsc2]

  rho         =   sqrt( (Two*rsa2/3d0) + Half*Ra2 )
  vlr         =   exp( - 6d0 * (rho-1.85d0) )


  do n = 1,3
    ep(n)     =   exp(-b1*SumR2(n)-b2*SumRs2(n))
    damp      =   adisp*exp(-bdisp*((R(n)-rdisp)**2))
    if      ( n == 1 ) then;  plm = pla(3)
    else if ( n == 2 ) then;  plm = plb(3)
    else;                     plm = plc(3)
    end if
    term0     =   - c62 * plm * damp / ( rs(n)**6 + dc6 )
    vlr       =   vlr + term0
  end do

  do n = 1,52
    pl(1)     =   pla( indx3(n)+1 )
    pl(2)     =   plb( indx3(n)+1 )
    pl(3)     =   plc( indx3(n)+1 )
    call ptud_NoGrad( indx1(n), indx2(n), R, rs, ep, pl, tp(n) )
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine ptud_NoGrad( ns, ms, R, rs, ep, pl, tps )
! This procedure determines the terms for an analytic fit for a molecule A3 using Jacobian coordinates
! Removed:
!  - Arguments: dr, drs, dep, dpl, dct
!  - Local:     dd

  use Parameters_Module       ,only:  Zero

  integer                                   ,intent(in)     ::    ns
  integer                                   ,intent(in)     ::    ms
  real(rkp) ,dimension(3)                   ,intent(in)     ::    R
  real(rkp) ,dimension(3)                   ,intent(in)     ::    rs
  real(rkp) ,dimension(3)                   ,intent(in)     ::    ep
  real(rkp) ,dimension(3)                   ,intent(in)     ::    pl
  real(rkp)                                 ,intent(out)    ::    tps

  integer                                                   ::    n
  real(rkp)                                                 ::    dtp
  real(rkp) ,dimension(3)                                   ::    tp

  tps       =   Zero
  do n = 1,3
    tp(n)   =   R(n)**ns * rs(n)**ms * pl(n)
    dtp     =   tp(n) * ep(n)
    tps     =   tps + dtp
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
