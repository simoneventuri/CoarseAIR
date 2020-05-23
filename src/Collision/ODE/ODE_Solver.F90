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

Module ODE_Solver

#include "../../qct.inc"

  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  CheckVariable
  use Parameters_Module     ,only:  rkp
  use Parameters_Module     ,only:  Half

  implicit none

  private
  public    ::    ODE_Solver_Type

  Type      ::    ODE_Solver_Type
    logical                                             ::    Initialized =   .False.
    logical                                             ::    IncStpSzFlg =   .True.
    real(rkp)                                           ::    eps         =   1.0E-10_rkp       !< Relative error control parameter (NOT USED)
    integer                                             ::    NSteps      =   4                 !< Number of stepsizes to use.
    real(rkp)                                           ::    relax       =   Half              !< Relaxation parameter for stepsize increase: if = 0: stepsize increase, if > 1, no stepsize increase
    logical                                             ::    NanCheck    =   .False.
    procedure(RHS_Evaluation),pointer  ,nopass ,private ::    EvaluateRHS =>  null()
  contains
    procedure     ,public   ::    Initialize  =>    Initialize_ODE_Solver
    generic       ,public   ::    Integrate   =>    Integrate_ODE_0d, Integrate_ODE_1d, Integrate_Trajectories
    procedure     ,private  ::    Integrate_ODE_0d
    procedure     ,private  ::    Integrate_ODE_1d
    procedure     ,private  ::    Integrate_Trajectories
  End Type

  Abstract Interface
    PURITY Subroutine RHS_Evaluation( iPES, r, drdx )
      import rkp
      integer   ,dimension(:)         ,intent(in)     ::    iPES
      real(rkp) ,dimension(:,:)       ,intent(in)     ::    R
      real(rkp) ,dimension(:,:)       ,intent(out)    ::    dRdx
    End Subroutine
  End interface

  logical   ,parameter    ::    i_Debug_Global  =   .True.!.False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_ODE_Solver( This, eps, NSteps, IncStpSzFlg, Relax, NanCheck, EvaluateRHS, i_Debug )
  class(ODE_Solver_Type)                ,intent(out)    ::    This
  real(rkp)                   ,optional ,intent(in)     ::    eps                   !< Relative error control parameter
  integer                     ,optional ,intent(in)     ::    NSteps                !< Number of stepsizes to use.
  logical                     ,optional ,intent(in)     ::    IncStpSzFlg
  real(rkp)                   ,optional ,intent(in)     ::    relax                 !< Relaxation parameter for stepsize increase : if = 0: stepsize increase, if > 1, no stepsize increase
  logical                     ,optional ,intent(in)     ::    NanCheck              !< NaN check indicator (To be used with debug run)
  procedure(RHS_Evaluation)   ,optional                 ::    EvaluateRHS
  logical                     ,optional ,intent(in)     ::    i_Debug               !< Debugging indicator

  logical                                               ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_ODE_Solver" )
  !i_Debug_Loc   =     Logger%On()

  nullify( This%EvaluateRHS )
  if ( present(eps     ) )    This%eps            =    eps
  if ( present(NSteps  ) )    This%NSteps         =    NSteps
  if ( present(IncStpSzFlg) ) This%IncStpSzFlg    =    IncStpSzFlg
  if ( present(relax   ) )    This%relax          =    relax
  if ( present(NanCheck) )    This%NanCheck       =    NanCheck
  if ( present(EvaluateRHS) ) This%EvaluateRHS    =>   EvaluateRHS

  This%Initialized    =     .True.

  if (i_Debug_Loc) then
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: Solver parameters:')")
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: -> Relative error control parameter:           This%eps         = ',es15.8)") This%eps
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: -> Increase Step Size Falg:                    This%IncStpSzFlg = ',l)")      This%IncStpSzFlg
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: -> Number of stepsizes to use:                 This%NSteps      = ',g0)")     This%NSteps
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: -> Relaxation parameter for stepsize increase: This%relax       = ',es15.8)") This%relax
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: -> NaN check indicator:                        This%NanCheck    = ',g0)")     This%NanCheck
    write(Logger%Unit,"(6x,'[Initialize_ODE_Solver]: Exiting')")
  end if
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Integrate_Trajectories( This, Traj, i_Debug )
! This procedure performs a fixed order bulirsch-stoer integration of coupled ode's with error control.
! It uses the 'EvalRHS' procedures (passed as an input argument) to compute the derivatives given y, x and NEqtTot.
! The dimension of the solution vector is 'y(NTraj,NEqtTot)' with:
! * NTraj = number of trajectories
! * NEqtTot  = number of coupled ode's

  use Trajectories_Class     ,only:  Trajectories_Type
  
  class(ODE_Solver_Type)                    ,intent(in)     ::    This
  type(Trajectories_Type)                   ,intent(inout)  ::    Traj
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    N , i
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Integrate_Trajectories" )
  !i_Debug_Loc   =     Logger%On()
  
  N = Traj%NTraj
  
  call This%Integrate( Traj%iPES(1:N), Traj%PaQ(:,1:N), Traj%t(1:N), Traj%dt(1:N), Traj%erra(:,1:N), Traj%smax(1:N), Traj%smin(1:N), Traj%irej(1:N), Traj%iste(1:N), i_Debug=.false. )
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Integrate_ODE_0d( This, iPES, y, x, h, yscal, smax, smin, irej, iste, i_Debug )
! This procedure performs a fixed order bulirsch-stoer integration of coupled ode's with error control.
! It uses the 'EvalRHS' procedures (passed as an input argument) to compute the derivatives given y, x and NEqtTot.
! The dimension of the solution vector is 'y(NTraj,NEqtTot)' with:
! * NTraj = number of trajectories
! * NEqtTot  = number of coupled ode's

  class(ODE_Solver_Type)                  ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    iPES 
  real(rkp) ,dimension(:)                 ,intent(inout)  ::    y                     ! Input: Solution at the initial time. Output: Solution at the final time. Dim=(NEqtTot)
  real(rkp)                               ,intent(inout)  ::    x                     ! Input: initial "time". Output: final "time"
  real(rkp)                               ,intent(inout)  ::    h                     ! Input: initial "time" stepsize. Output: estimate of optimum stepsize for next "time" step
  real(rkp) ,dimension(:)                 ,intent(in)     ::    yscal                 ! absolute error control vector. Dim=(NEqtTot)
  real(rkp)                     ,optional ,intent(inout)  ::    smax                  ! maximum stepsize so far
  real(rkp)                     ,optional ,intent(inout)  ::    smin                  ! minimum stepsize so far
  integer                       ,optional ,intent(inout)  ::    irej                  ! number of rejected steps so far
  integer                       ,optional ,intent(inout)  ::    iste                  ! number of good steps so far
  logical                       ,optional ,intent(in)     ::    i_Debug

  integer   ,parameter                                    ::    NTraj = 1
  integer   ,dimension( NTraj )                           ::    iPES_
  real(rkp) ,dimension( size(y), NTraj )                  ::    y_
  real(rkp) ,dimension( NTraj )                           ::    x_
  real(rkp) ,dimension( NTraj )                           ::    h_
  real(rkp) ,dimension( size(y), NTraj )                  ::    yscal_
  real(rkp) ,dimension( NTraj )                   ,target ::    smax_1d
  real(rkp) ,dimension( NTraj )                   ,target ::    smin_1d
  integer   ,dimension( NTraj )                   ,target ::    irej_1d
  integer   ,dimension( NTraj )                   ,target ::    iste_1d
  real(rkp) ,dimension(:) ,pointer                        ::    smax_
  real(rkp) ,dimension(:) ,pointer                        ::    smin_
  integer   ,dimension(:) ,pointer                        ::    irej_
  integer   ,dimension(:) ,pointer                        ::    iste_

  y_(:,1)     =     y
  x_(1)       =     x
  h_(1)       =     h
  yscal_(:,1) =     yscal
  smin_       =>    null()
  irej_       =>    null()
  iste_       =>    null()


  smax_       =>    null()
  if ( present(smax) ) then
    smax_1d   =     smax
    smax_     =>    smax_1d
  end if

  smin_       =>    null()
  if ( present(smin) ) then
    smin_1d   =     smin
    smin_     =>    smin_1d
  end if

  irej_       =>    null()
  if ( present(irej) ) then
    irej_1d   =     irej
    irej_     =>    irej_1d
  end if

  iste_       =>    null()
  if ( present(iste) ) then
    iste_1d   =     iste
    iste_     =>    iste_1d
  end if


  call This%Integrate( iPES_, y_, x_, h_, yscal_, smax_, smin_, irej_, iste_, i_Debug )

  y           =     y_(:,1)
  x           =     x_(1)
  h           =     h_(1)
  if ( present(smax) ) smax = smax_(1)
  if ( present(smin) ) smin = smin_(1)
  if ( present(irej) ) irej = irej_(1)
  if ( present(iste) ) iste = iste_(1)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Integrate_ODE_1d( This, iPES, y, x, h, yscal, smax, smin, irej, iste, i_Debug )
! This procedure performs a fixed order bulirsch-stoer integration of coupled ode's with error control.
! It uses the 'EvalRHS' procedures (passed as an input argument) to compute the derivatives given y, x and NEqtTot.
! The dimension of the solution vector is 'y(NTraj,NEqtTot)' with:
! * NTraj    = number of trajectories
! * NEqtTot  = number of coupled ode's

  use Parameters_Module   ,only:  Zero, One
  use, intrinsic :: IEEE_ARITHMETIC

  class(ODE_Solver_Type)                  ,intent(in)     ::    This
  integer   ,dimension(:)                 ,intent(in)     ::    iPES 
  real(rkp) ,dimension(:,:)               ,intent(inout)  ::    y         ! Input: Solution at the initial time. Output: Solution at the final time. Dim=(NEqtTot,NTraj)
  real(rkp) ,dimension(:)                 ,intent(inout)  ::    x         ! Input: initial "time". Output: final "time". Dim=(NTraj)
  real(rkp) ,dimension(:)                 ,intent(inout)  ::    h         ! Input: initial "time" stepsize. Output: estimate of optimum stepsize for next "time" step. Dim=(NTraj)
  real(rkp) ,dimension(:,:)               ,intent(in)     ::    yscal     ! absolute error control vector. Dim=(NEqtTot,NTraj)
  real(rkp) ,dimension(:)       ,optional ,intent(inout)  ::    smax      ! maximum stepsize so far. Dim=(NTraj)
  real(rkp) ,dimension(:)       ,optional ,intent(inout)  ::    smin      ! minimum stepsize so far. Dim=(NTraj)
  integer   ,dimension(:)       ,optional ,intent(inout)  ::    irej      ! number of rejected steps so far. Dim=(NTraj)
  integer   ,dimension(:)       ,optional ,intent(inout)  ::    iste      ! number of good steps so far. Dim=(NTraj)
  logical                       ,optional ,intent(in)     ::    i_Debug

  logical                                                 ::    i_Debug_Loc

  character(*)                            ,parameter      ::    ProcName = 'Integrate_ODE_1d'
  integer                                 ,parameter      ::    imax = 11
  integer   ,dimension(imax)              ,parameter      ::    nseq = [2,4,6,8,12,16,24,32,48,64,96]

  integer                                                 ::    NTraj     ! number of trajectories
  integer                                                 ::    NEqtTot   ! number of coupled ode's
  integer                                                 ::    iSteps, iTraj, iEqtTot, j, iconv
  real(rkp)                                               ::    htemp
  real(rkp) ,dimension(             size(y,2) )           ::    fv1
  real(rkp) ,dimension( size(y,1) , size(y,2) )           ::    dydx
  real(rkp) ,dimension( size(y,1) , size(y,2) )           ::    yseq
  real(rkp) ,dimension(             size(y,2) , imax )    ::    x2
  real(rkp) ,dimension( size(y,1) , size(y,2) , imax )    ::    ytry
  real(rkp) ,dimension( size(y,1) , size(y,2) )           ::    yerr
!   NOT DONE
  real(rkp) ,dimension(             size(y,2) , imax )    ::    errmax  ! NTraj
  real(rkp) ,dimension( size(y,1) , size(y,2) , 7  )      ::    d

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Integrate_ODE_1d" )
  !i_Debug_Loc   =     Logger%On()
  
  NEqtTot   =     size(y,1)             ! Getting the number of coupled ode's
  NTraj     =     size(y,2)             ! Getting the number of trajectories


  if ( .Not. associated(This%EvaluateRHS) ) then
    write(*,*) ' Integrate_ODE_1d Error 1'
    write(Logger%Unit,"(6x,'[Integrate_ODE_1d]: This%EvaluateRHS is not associated')") ; error stop
  end if

  call This%EvaluateRHS( iPES, y, dydx )
  
  d   =   Zero
  x2  =   Zero

  do iSteps = 1,This%NSteps
    
    do iTraj=1,NTraj
      if (This%NanCheck) call CheckVariable( h(iTraj), ProcName=ProcName, VarName="h before mmidv" )
    end do
    call mmidv( iPES, y, dydx, h, nseq(iSteps), yseq, This%EvaluateRHS )    ! yseq:out, all other inputs
    do iTraj=1,NTraj
      if (This%NanCheck) call CheckVariable( h(iTraj), ProcName=ProcName, VarName="h after mmidv" )
    end do
    
    htemp       =   One / dfloat(nseq(iSteps))
    fv1(:)      =   ( h(:) * htemp )**2    ! Computing the current value of independent variable
    call RatPolExtrap( This%NSteps, iSteps, fv1(:), yseq, ytry(:,:,iSteps), imax, d, x2 )
    if (This%NanCheck) call CheckVariable( ytry(:,:,iSteps), ProcName=ProcName, VarName="ytry" )
    
  end do
  
  do iSteps=1,This%NSteps-1
    do iTraj = 1,NTraj
      errmax(iTraj,iSteps) = Zero
    end do
    do iEqtTot=1,NEqtTot
      do iTraj=1,NTraj
      
       yerr(iEqtTot,iTraj) = dabs( (ytry(iEqtTot,iTraj,iSteps)-ytry(iEqtTot,iTraj,This%NSteps)) / (yscal(iEqtTot,iTraj)) )

       if (IEEE_IS_NAN(yerr(iEqtTot,iTraj))) then
         write(Logger%Unit,"(4x,'[Initialize_ODE_Solver]: ERROR: Got NaN for yerr')")
         write(*,*) ' Integrate_ODE_1d Error 2'
         stop
       end if
         
      end do
    end do
    do iEqtTot=1,NEqtTot
      do iTraj = 1,NTraj
        errmax(iTraj,iSteps) = max(errmax(iTraj,iSteps), yerr(iEqtTot,iTraj))
      end do
    end do
  end do

  do iTraj = 1,NTraj
  
    do iSteps = 2,This%NSteps
      j         =   This%NSteps+1-iSteps
      if ( errmax(iTraj,j) >= One ) exit
      iconv     =   j
    end do

    if ( errmax(iTraj,This%NSteps-1) < One ) then                                                                                 ! Step converged
      x(iTraj)     =   x(iTraj) + h(iTraj)
      y(:,iTraj)   =   ytry(:,iTraj,This%NSteps)
      if ( present(smin) ) smin(iTraj)  =   min( smin(iTraj) , h(iTraj) )
      if ( present(smax) ) smax(iTraj)  =   max( smax(iTraj) , h(iTraj) )
      if ( present(iste) ) iste(iTraj)  =   iste(iTraj) + 1
      if ( ( iconv < This%NSteps-1 ) .and. (This%IncStpSzFlg) ) then                                                                                           ! Can we increase stepsize?: yes because after This%NSteps-1 stepsizes we were also converged
        h(iTraj)   =   h(iTraj) * ( dfloat(nseq(This%NSteps)) / dfloat(nseq(iconv+1)) + This%Relax ) / ( One + This%Relax )
        if (This%NanCheck) call CheckVariable( h(iTraj), ProcName=ProcName, VarName="h before mmidv" )
      end if
    else                                                                                                                          ! Step is not converged. decrease wanted step.
      h(iTraj)     =   h(iTraj) / dfloat( nseq(max(1,This%NSteps/2)) )
      if (This%NanCheck) call CheckVariable( h(iTraj), ProcName=ProcName, VarName="h before mmidv" )
      if ( present(irej) ) irej(iTraj)  =   irej(iTraj) + 1
    end if
  end do
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine mmidv( iPES, y, dydx, htot, nstep, yout, EvaluateRHS )
! This procedure performs a single modified midpoint integration step for coupled ode's.
! On input:
!  * y - initial conditions
!  * dydx - initial derivatives
!  * NEqtTot - number of coupled ode's
!  * htot - total "time" stepsize
!  * nstep - number of substeps to take
!  * NTraj - number of trajectories
! On output:
!  * yout - final conditions

  use Parameters_Module   ,only:  rkp, One, Two, Half

  integer   ,dimension(:)               ,intent(in)     ::    iPES 
  real(rkp) ,dimension(:,:)             ,intent(in)     ::    y             !< initial conditions. Dim=(NEqtTot,NTraj)
  real(rkp) ,dimension(:,:)             ,intent(in)     ::    dydx          !< initial derivatives. Dim=(NEqtTot,NTraj)
  real(rkp) ,dimension(:)               ,intent(in)     ::    htot          !< total "time" stepsize. Dim=(NTraj)
  integer                               ,intent(in)     ::    nstep         !< number of substeps to take
  real(rkp) ,dimension(:,:)             ,intent(out)    ::    yout          !< final conditions. Dim=(NEqtTot,NTraj)
  procedure(RHS_Evaluation) ,pointer                    ::    EvaluateRHS

  character(*)                              ,parameter  ::    ProcName = 'mmidv'
  integer                                               ::    NTraj         !< Number of trajectories
  integer                                               ::    NEqtTot          !< Number of coupled ode's
  integer                                               ::    n, k
  integer                                               ::    i, itemp
  integer                                               ::    iym, iyn, iys
  real(rkp)                                             ::    htemp
  real(rkp) ,dimension( size(y,1), size(y,2), 3 )       ::    yn            !< scratch arrays. Dim=(NEqtTot,NTraj,3)
  real(rkp) ,dimension( size(y,2) )                     ::    h, h2         !< Dim=(NTraj)

  NEqtTot   =   size(y,1)
  NTraj     =   size(y,2)
  htemp     =   One / dfloat(nstep)
  h         =   htot * htemp
  h2        =   Two * h
!   if (This%NanCheck) call CheckVariable( h, ProcName=ProcName, VarName="h" )
!   if (This%NanCheck) call CheckVariable( dydx, ProcName=ProcName, VarName="dydx" )
  
  !write(*,*)
  !write(*,*) 'y(:,1)    =', y(:,1)
  !write(*,*) 'h         =', h
  !write(*,*) 'dydx(:,1) =', dydx(:,1)
  
  do k = 1,NTraj
  do i = 1,NEqtTot,2
    yn(i  ,k,1)   =   y(i  ,k)
    yn(i  ,k,2)   =   y(i  ,k) + h(k) * dydx(i,k)
    yn(i+1,k,1)   =   y(i+1,k)
    yn(i+1,k,2)   =   y(i+1,k) + h(k) * dydx(i+1,k)
  end do
  end do

  iym             =     1
  iyn             =     2
  iys             =     3
  
  
  !write(*,*) 'yn(1,:,iyn) =', yn(1,:,iyn)
  !write(*,*) 'yn(2,:,iyn) =', yn(2,:,iyn)
  
  call EvaluateRHS( iPES, yn(:,:,iyn), yout(:,:) )
  
  !write(*,*) 'yout(1,:) =',yout(1,:)
  !write(*,*) 'yout(2,:) =',yout(2,:)
  !write(*,*)
  
  do n = 2,nstep
    do k = 1,NTraj
    do i = 1,NEqtTot,2
      yn(i  ,k,iys)   =   yn(i,  k,iym) + h2(k) * yout(i  ,k)
      yn(i+1,k,iys)   =   yn(i+1,k,iym) + h2(k) * yout(i+1,k)
!       if (This%NanCheck) call CheckVariable( yn(:,i,iys),   ProcName=ProcName, VarName="yn(:,i,iys)"   )
!       if (This%NanCheck) call CheckVariable( yn(:,i+1,iys), ProcName=ProcName, VarName="yn(:,i+1,iys)" )
    end do
    end do
    itemp   =   iym
    iym     =   iyn
    iyn     =   iys
    iys     =   itemp
    
    !write(*,*) 'n           =', n
    !write(*,*) 'h2          =', h2
    !write(*,*) 'yn(1,:,iyn) =', yn(1,:,iyn)
    !write(*,*) 'yn(1,:,iyn) =', yn(2,:,iyn)
    
    call EvaluateRHS( iPES, yn(:,:,iyn), yout(:,:) )
    
    !write(*,*) 'yout(1,:) =',yout(1,:)
    !write(*,*) 'yout(2,:) =',yout(2,:)
    !write(*,*)
    
  end do

  do k = 1,NTraj
    yout(:,k)   =   Half * ( yn(:,k,iym) + yn(:,k,iyn) + h(k) * yout(:,k) )
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
PURITY Subroutine RatPolExtrap( NSteps, iest, xest, yest, yz, imax, d, x )
! This procedure performs a rational polynomial extrapolation of a vector.
! Input variables:
!  * iest - number of this call in a sequence of calls
!  * xest - current value of independent variable
!  * yest - function at xest
!  * nv - number elements in vectors
!  * NSteps - number of previous vectors to use
!  * NTraj - number of trajectories
! Output variables:
!     yz - extrapolated result
!     nmax ge nv, ncol ge NSteps, imax ge iest
! Subroutine rzextrv( iest, xest, yest, yz, nv, NSteps, NTraj, imax )

  integer                               ,intent(in)     ::    NSteps
  integer                               ,intent(in)     ::    iest    ! number of this call in a sequence of calls
  real(rkp) ,dimension(:)               ,intent(in)     ::    xest    ! current value of independent variable. Dim=(NTraj)
  real(rkp) ,dimension(:,:)             ,intent(in)     ::    yest    ! function at xest. Dim=(NEqtTot,NTraj)
  real(rkp) ,dimension(:,:)             ,intent(out)    ::    yz      ! extrapolated result. Dim=(NEqtTot,NTraj)
  integer                               ,intent(in)     ::    imax
  real(rkp) ,dimension(:,:,:)           ,intent(inout)  ::    d       ! Dim=(NEqtTot,NTraj,7)
  real(rkp) ,dimension(:,:)             ,intent(inout)  ::    x       ! Dim=(NTraj,imax)=(NTraj,11)

  integer   ,parameter                                  ::    ncol  = 7
  real(rkp) ,parameter                                  ::    eps   = 1.0E-40_rkp

  integer                                               ::    NEqtTot
  integer                                               ::    m1, k, j
  real(rkp) ,dimension(size(xest) ,ncol)                ::    fx      ! Dim=(NTraj,7)
  real(rkp) ,dimension(size(xest) )                     ::    yy      ! Dim=(NTraj)
  real(rkp) ,dimension(size(xest) )                     ::    v       ! Dim=(NTraj)
  real(rkp) ,dimension(size(xest) )                     ::    c       ! Dim=(NTraj)
  real(rkp) ,dimension(size(xest) )                     ::    b1      ! Dim=(NTraj)
  real(rkp) ,dimension(size(xest) )                     ::    b       ! Dim=(NTraj)
  real(rkp) ,dimension(size(xest) )                     ::    ddy     ! Dim=(NTraj)

  NEqtTot  = size(yest,1)

!   if ( iest > imax ) then
!     write(Logger%Unit,"(8x,'[rzextr]: iest=',i15,' exceeds dimensions =',i5)") iest, imax
!     stop
!   end if
!   if ( This%NSteps > ncol ) then
!     write(Logger%Unit,"(8x,'[rzextr]: This%NSteps=',i5,' exceeds dimensions =',i5)") This%NSteps, ncol
!     stop
!   end if

  x(:,iest)       =     xest(:)
  if ( iest == 1 ) then
    yz(:,:)       =   yest(:,:)
    d(:,:,1)      =   yest(:,:)
  else
    m1            =   min(iest,NSteps)
    do k = 1,m1-1
      fx(:,k+1)   =   x(:,iest-k) / xest(:)
    end do
    do j = 1,size(yest,1)   ! Loop on equations
      yy(:)       =   yest(j,:)
      v(:)        =   d(j,:,1)
      c(:)        =   yy(:)
      d(j,:,1)    =   yy(:)
      do k = 2,m1
        b1(:)     =   fx(:,k) * v(:)
        b(:)      =   ( c(:) - v(:) ) / ( (b1(:) - c(:) ) + eps )
        ddy(:)    =   c(:) * b(:)
        c(:)      =   b1(:) * b(:)
        v(:)      =   d(j,:,k)
        d(j,:,k)  =   ddy(:)
        yy(:)     =   yy(:) + ddy(:)
      yz(j,:)     =   yy(:)
      end do
    end do
  end if

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
