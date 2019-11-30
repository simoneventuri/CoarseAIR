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

Module MolecularState_Class

#include "../../qct.inc"

  use Parameters_Module         ,only:  rkp, Zero
  use Logger_Class              ,only:  Logger
  use MolecularStateParam_Class ,only:  MolecularStateParam_Type

  implicit none

  private
  public    ::    MolecularState_Type

  logical   ,parameter    ::    i_Debug_Global = .True.!.False.
  integer   ,parameter    ::    NSpace = 3

  Type      ::    MolecularState_Type
    integer                               ::    iPair   =   0
    integer   ,dimension(2)               ::    iAtoms  =   0
    integer                               ::    i       =   0           ! Index of internal level: i(v,j)
    integer                               ::    v       =   0           ! Vibrational quantum number
    integer                               ::    j       =   0           ! Rotational quantum number
    integer                               ::    itype                   ! Arrangement quantum number
    real(rkp)                             ::    Rpi     =   Zero        ! Inverse of the distances of atom-atom pairs [1/bohr]
    real(rkp)                             ::    Eint    =   Zero        ! Internal energy
    real(rkp)                             ::    AngMom  =   Zero        ! Anglular momentum
    real(rkp)                             ::    viba
    type(MolecularStateParam_Type)        ::    Param
    real(rkp) ,dimension( NSpace,2 )      ::    x                       ! Positions of the two atoms making up current molecule. Dim=(NSpace,NAtoms)
    real(rkp) ,dimension( NSpace,2 )      ::    xdot                    ! Velocities of the two atoms making up current molecule. Dim=(NSpace,NAtoms)
  contains
    procedure ,public   ::    FindState
  End Type

  contains
  
  
!________________________________________________________________________________________________________________________________!
Subroutine FindState( This, Collision, ierr, i_Debug )
! This procedure determines the final state of the diatomic molecule.
!   xa contains the coordinates of atom a
!   xdota contains the time derivatives of the coordinates of atom a
!   xma is the mass of a
!   xb etc ditto except for b.
!   rcent upper limit on position of centrifugal barriar maximum
!   rsmal lower limit on effective potential minimum.
!   nquad is the order of quadrature to use for action integral
!   vibact,ti final vibrational action and rotational angular momentum
!   ibound is 1 for bound state, 2 for quasi-bound state, and 3 for not bound.
!   subtract 1 from ibound to be consistent with rest of code
!   ierr is Zero on normal return, otherwise it is equal to 1.

! This procedure identifies the internal state associated to a given molecule. In the passed-object argument, the following variables are used as inputs:
! * iPair:  the index of the atom-atom pair associated to current molecule
! * x:      the positions of each atoms within current state
! * xdot:   the velocities of each atoms within current state
! and, the following variables are used as outputs:
! * Eint:   the internal energy
! * itype:  the arrangement quantum number
! * AngMom: the angular momentum (rot. quantum number)
! * viba:   the real vibrational quantum number

  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero, One, Half, Two, onem, Pi

  class(MolecularState_Type)                ,intent(inout)  ::    This
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(out)    ::    ierr
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  real(rkp) ,parameter                                      ::    epss = 1.0E-04_rkp

  integer                                                   ::    iatom1, iatom2
  integer                                                   ::    i, imin, ivmx, itry, ipass
  integer                                                   ::    iA
  real(rkp)                                                 ::    xma
  real(rkp)                                                 ::    xmb
  real(rkp)                                                 ::    Fa, Fb
  real(rkp)                                                 ::    Rbond
  real(rkp)                                                 ::    tx, ty, tz, tis
  real(rkp)                                                 ::    rms, twomu, eps
  real(rkp)                                                 ::    Vc_R2         ! Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                                 ::    Rmin, Vmin
  real(rkp)                                                 ::    Rmax, Vmax
  real(rkp)                                                 ::    Rb, vib, r1, r00, ri, ro, rga, rgb

  real(rkp)                                                 ::    dEint, dAdE
  real(rkp)                                                 ::    Eint_p, Ri_p, Ro_p, Action_p
  real(rkp)                                                 ::    Eint_m, Ri_m, Ro_m, Action_m

  real(rkp)                                                 ::    gam, Action, r3, prob
  real(rkp)                                                 ::    xlife, rho2
  real(rkp)                                                 ::    Mass

  real(rkp) ,dimension(NSpace)                              ::    xcm
  real(rkp) ,dimension(NSpace)                              ::    xdotcm
  real(rkp) ,dimension(NSpace,2)                            ::    x
  real(rkp) ,dimension(NSpace,2)                            ::    xdot
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "FindState")
  !i_Debug_Loc   =     Logger%On()


  ierr      =   0
  
  
  if (trim(adjustl(Collision%Pairs(This%iPair)%DiaPot%Name)) == '<NULL>') then
    ! The final pair is not going to be treated as a molecule 
    
    This%itype     =   3
    This%viba      =   Half
    
  else
  
  ! Getting the properties of the two atoms wihting the pair This%iPair
    iatom1    =   Collision%Pairs(This%iPair)%To_Atoms(1)
    iatom2    =   Collision%Pairs(This%iPair)%To_Atoms(2)
    xma       =   Collision%Atoms(iatom1)%Mass
    xmb       =   Collision%Atoms(iatom2)%Mass


  ! ==============================================================================================================
  !     SHIFT ATOM COORDINATES AND VELOCITY SO CENTER-OF-MASS IS STATIONARY AT ORIGIN
  ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Shift atom coordinates and velocity so center-of-mass is stationary at origin')")
    fa                =   xma / ( xma + xmb )
    fb                =   xmb / ( xma + xmb )
    xcm               =   fa * This%x(:,1)    + fb * This%x(:,2)
    xdotcm            =   fa * This%xdot(:,1) + fb * This%xdot(:,2)
    This%itype        =   1
    This%Eint         =   Zero
    Rbond             =   Zero
    do i = 1,NSpace
      x(i,:)          =   This%x(i,:)    - xcm(i)
      xdot(i,:)       =   This%xdot(i,:) - xdotcm(i)
      Rbond           =   Rbond + ( x(i,1) - x(i,2) )**2
      This%Eint       =   This%Eint + Half * ( xma * xdot(i,1)**2 + xmb * xdot(i,2)**2 )
    end do
    Rbond             =   sqrt(Rbond)
    if (i_Debug_Loc) then
      write(Logger%Unit,"(10x,'[FindState]: -> x(:,1)    = ',*(es15.8,3x))") x(:,1)
      write(Logger%Unit,"(10x,'[FindState]: -> x(:,2)    = ',*(es15.8,3x))") x(:,2)
      write(Logger%Unit,"(10x,'[FindState]: -> xdot(:,1) = ',*(es15.8,3x))") xdot(:,1)
      write(Logger%Unit,"(10x,'[FindState]: -> xdot(:,2) = ',*(es15.8,3x))") xdot(:,2)
      write(Logger%Unit,"(10x,'[FindState]: -> Rbond     = ',es15.8)") Rbond
      write(Logger%Unit,"(10x,'[FindState]: -> This%Eint = ',es15.8)") This%Eint
    end if
  ! ==============================================================================================================


  ! ==============================================================================================================
  !     COMPUTING ANGULAR MOMENTUM
  ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Computing angular momentum')")
    tx            =   Zero
    ty            =   Zero
    tz            =   Zero
    do i = 1,2
      iA          =   Collision%Pairs(This%iPair)%To_Atoms(i)
      Mass        =   Collision%Atoms(iA)%Mass
      tx          =   tx + (x(2,i) * xdot(3,i) - x(3,i) * xdot(2,i)) * Mass
      ty          =   ty + (x(3,i) * xdot(1,i) - x(1,i) * xdot(3,i)) * Mass
      tz          =   tz + (x(1,i) * xdot(2,i) - x(2,i) * xdot(1,i)) * Mass
    end do
    tis           =   tx**2 + ty**2 + tz**2
    This%AngMom   =   sqrt(tis)
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%AngMom = ',es15.8)") This%AngMom
  ! ==============================================================================================================


  ! ==============================================================================================================
  !     ESTIMATING THE POSITION OF THE POTENTIAL MINIMUM AND CENTRIFUGAL BARRIER
  ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Estimating the position of the potential minimum and centrifugal barrier')")
    rms           =   xma * xmb / ( xma + xmb )     ! Reduced Mass
    twomu         =   rms * two
    Vc_R2         =   half * tis / rms
    eps           =   epss

    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling DiaPot%FindMinimum')")
    call Collision%Pairs(This%iPair)%DiaPot%FindMinimum( [This%Param%rsmal,This%Param%rcent], Vc_R2, eps, Rmin, Vmin, ierr )
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> Rmin = ',es15.8,3x,'Vmin = ',es15.8)") Rmin, Vmin

    imin          =   1
    if ( ierr /= 0 ) then
      Rmin        =   ( This%Param%rsmal + This%Param%rcent ) * Half
      imin        =   0
    end if
    

    if ( imin == 1 ) then
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling DiaPot%FindMaximum')")
      call Collision%Pairs(This%iPair)%DiaPot%FindMaximum( [Rmin,This%Param%rcent], Vc_R2, eps, Rmax, Vmax, ierr )
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> Rmax = ',es15.8,3x,'Vmax = ',es15.8)") Rmax, Vmax
      ivmx        =   1
      if ( ierr /= 0 ) then
        Rmax      =   This%Param%rcent
        ivmx      =   0
      end if
    end if
    

    Rb = Rbond


    if ( Rb > Rmax .or. imin == 0 ) then
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Case (Rb > Rmax .or. imin == 0)')")
      This%itype     =   3
      This%viba      =   Half
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%viba  = ',es15.8)") This%viba
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%itype = ',g0)") This%itype

    else
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Case (else)')")

      vib = Collision%Pairs(This%iPair)%DiaPot%DiatomicPotential( Rb )
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> Rbound = ',es15.8,3x,'VibEn = ',es15.8)") Rb, vib
      This%Eint  = This%Eint + vib
      if ( This%Eint < This%Param%Vinf ) then
        This%itype   =   1
      else
        This%itype   =   2
      end if
      r1  =   Rmin
    
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%Eint vs Vmax: ',es15.8,es15.8)") This%Eint, Vmax
      if ( This%Eint > Vmax ) then
        This%viba    =   Half
        This%itype   =   3
      else
      
        itry              =   1
        r00               =   This%Param%rsmal
        do
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling TurningPoint')")
          call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [r00,r1], Vc_R2, This%Eint, ri, ierr, NBisection=4, NNewton=4 )
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> ri  = ',es15.8)") ri

          if ( ri /= Zero ) exit
          itry            =   itry + 1
          r00             =   r00 * Half
          if ( itry <= 10) cycle
          write(Logger%Unit,"(10x,'[FindState]: Error locating inner turning point')")
          write(Logger%Unit,"(10x,'[FindState]: -> This%Eint = ',es15.8)") This%Eint
          write(Logger%Unit,"(10x,'[FindState]: -> Vc_R2           = ',es15.8)") Vc_R2
          write(Logger%Unit,"(10x,'[FindState]: -> r00            = ',es15.8)") r00
          write(Logger%Unit,"(10x,'[FindState]: -> r1             = ',es15.8)") r1
          write(Logger%Unit,"(10x,'[FindState]: -> Rb             = ',es15.8)") Rb
          write(Logger%Unit,"(10x,'[FindState]: -> imin           = ',g0)")     imin
          write(Logger%Unit,"(10x,'[FindState]: -> Rmin           = ',es15.8)") Rmin
          write(Logger%Unit,"(10x,'[FindState]: -> Vmin           = ',es15.8)") Vmin
          stop
        end do

        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling TurningPoint')")
        call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [r1,Rmax], Vc_R2, This%Eint, ro, ierr, NBisection=10, NNewton=4 )
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> ro  = ',es15.8)") ro

        if ( r1 == Zero ) then
          write(Logger%Unit,"(10x,'[FindState]: Error locating outer turning point')")
          write(Logger%Unit,"(10x,'[FindState]: -> This%Eint = ',es15.8)") This%Eint
          write(Logger%Unit,"(10x,'[FindState]: -> Vc_R2           = ',es15.8)") Vc_R2
          write(Logger%Unit,"(10x,'[FindState]: -> ri              = ',es15.8)") ri
          write(Logger%Unit,"(10x,'[FindState]: -> Rb              = ',es15.8)") Rb
          write(Logger%Unit,"(10x,'[FindState]: -> ivmx            = ',g0)")     ivmx
          write(Logger%Unit,"(10x,'[FindState]: -> Rmax            = ',es15.8)") Rmax
          stop
        end if
        

        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling ActionIntegral')")
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: ri, ro, Vc_R2 = ',3es15.8)") ri, ro, Vc_R2
        
        This%viba    =   Collision%Pairs(This%iPair)%DiaPot%ActionIntegral( ri, ro, Vc_R2, This%Eint, This%Param%nquad, TwoMu, QuadratureType=3 )
        
        
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%viba  = ',es15.8)") This%viba

        gam         =   Zero
        
    
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%Eint, This%Param%Vinf = ',2es15.8)") This%Eint, This%Param%Vinf
        if ( This%Eint > This%Param%Vinf ) then
          rga       =   Rmax
          rgb       =   Rmax*Two
          ipass     =   0
    

          do
            ipass   =   ipass+1
            if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling turn; ipass = ',g0)") ipass
            
            call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [rga,rgb], Vc_R2, This%Eint, r3, ierr, NBisection=10, NNewton=4 )
            if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> r3 = ',es15.8)") r3
            if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> ierr = ',g0)") ierr

            if ( ierr==0)exit
            rga     =   rgb
            rgb     =   rgb * Two
            if ( ipass < 10)cycle
            write(Logger%Unit,"(10x,'[FindState]: error locating outer tunneling turning point')")
            write(Logger%Unit,"(10x,'[FindState]: -> This%Eint = ',es15.8)") This%Eint
            write(Logger%Unit,"(10x,'[FindState]: -> Vc_R2     = ',es15.8)") Vc_R2
            stop
          end do
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> ipass    = ',g0)") ipass
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> rga, rgb = ',2es15.8)") rga, rgb


  ! ==============================================================================================================
  !        EVALUATE DERIVATIVE OF THE ACTION WRT THE ENERGY
  ! ==============================================================================================================
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Evaluate derivative wrt This%Eint')")
          Eint_p    =   This%Eint + epss
          if ( Eint_p > Vmax ) Eint_p = This%Eint + half * ( Vmax - This%Eint )
          dEint     =   Eint_p - This%Eint
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%Eint = ',es15.8,3x,'dEint = ',es15.8)") This%Eint, dEint
          call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [r00,r1] , Vc_R2, Eint_p, Ri_p, ierr, NBisection=4 , NNewton=4 )
          call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [r1,Rmax], Vc_R2, Eint_p, Ro_p, ierr, NBisection=10, NNewton=4 )
          Action_p  =   Collision%Pairs(This%iPair)%DiaPot%ActionIntegral( Ri_p, Ro_p, Vc_R2, Eint_p, This%Param%nquad, TwoMu, QuadratureType=3 )
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: ->  Eint_p = ',es15.8,3x,'Ri_p = ',es15.8,3x,'Ro_p = ',es15.8,3x,'Action_p = ',es15.8)") Eint_p, Ri_p, Ro_p, Action_p

          Eint_m    =   This%Eint - dEint
          call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [r00,r1] , Vc_R2, Eint_m, Ri_m, ierr, NBisection=4, NNewton=4 )
          call Collision%Pairs(This%iPair)%DiaPot%TurningPoint( [r1,Rmax], Vc_R2, Eint_m, Ro_m, ierr, NBisection=10, NNewton=4 )
          Action_m  =   Collision%Pairs(This%iPair)%DiaPot%ActionIntegral( Ri_m, Ro_m, Vc_R2, Eint_m, This%Param%nquad, TwoMu, QuadratureType=3 )
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: ->  Eint_m = ',es15.8,3x,'Ri_m = ',es15.8,3x,'Ro_m = ',es15.8,3x,'Action_m = ',es15.8)") Eint_m, Ri_m, Ro_m, Action_m

          dAdE      =   half * ( Action_p - Action_m ) / dEint
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%dAdE = ',es15.8)") dAdE
  ! ==============================================================================================================

          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Calling ActionIntegral for tunneling')")
          Action    =   Collision%Pairs(This%iPair)%DiaPot%ActionIntegral( ro, r3, Vc_R2, This%Eint, This%Param%nquad, TwoMu, QuadratureType=3 )
          prob      =   exp( - Action )
          gam       =   prob / dAdE
          xlife     =   One / gam
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> xlife = ',es15.8,3x,'This%Param%tthresau = ',es15.8,3x)") xlife, This%Param%tthresau
          if ( xlife < This%Param%tthresau ) then 
            if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: Case xlife < This%Param%tthresau; This%itype set to 4')")
            This%itype = 4
          end if
          
          if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: dup action')")
          Action    =   Collision%Pairs(This%iPair)%DiaPot%ActionIntegral( ro, r3, Vc_R2, This%Eint, This%Param%nquad, TwoMu, QuadratureType=3 ) * Half / Pi
          rho2      =   GammaFct(Action) - Action * ( log( abs(Action) ) - One )
          This%viba =   This%viba + rho2 * This%Param%rfact
       
        end if

        This%viba     =   This%viba * Half / Pi
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[FindState]: -> This%viba = ',es15.8)") This%viba

      end if

    end if
  
  end if

  This%itype = This%itype - 1
  if ( This%itype == 2) This%AngMom = Half
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Function GammaFct( y ) result( SumTerm )
! This procedure computes arg gamma(1/5+y*i)

  use Parameters_Module   ,only:  rkp, Zero, Half

  real(rkp)     ,intent(in)     ::    y
  real(rkp)                     ::    SumTerm

  integer   ,parameter          ::    IterMax =   1000000
  real(rkp) ,parameter          ::    g1      =  -1.963510026021423d0
  real(rkp) ,parameter          ::    TolAbs  =   1.0E-06_rkp
  integer                       ::    Iter
  real(rkp)                     ::    Term, T1, T2

  SumTerm   =   g1 * y
  if ( y == Zero )return

  Iter      =   0
  do
    T1      =   y / ( Half + Iter )
    T2      =   atan(t1)
    Term    =   T1 - T2
    SumTerm =   SumTerm + Term
    if ( abs(Term/SumTerm) < TolAbs ) return
    Iter    =   Iter + 1
    if ( Iter > IterMax ) exit
  end do

!   write(*,"('[argg]: Failure to compute gamma')")
!   write(*,"('[argg]: y       = ',es15.8)") y
!   write(*,"('[argg]: TolAbs  = ',es15.8)") TolAbs
!   write(*,"('[argg]: Iter    = ',g0)") Iter
!   write(*,"('[argg]: SumTerm = ',es15.8)") SumTerm
!   write(*,"('[argg]: Term    = ',es15.8)") Term
   stop

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
