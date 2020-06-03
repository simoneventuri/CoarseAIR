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

Module Transformation_Class

  use Parameters_Module     ,only:  rkp, Zero, One, Two
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Transformation_Type, R_to_X, X_to_R, dX_To_dR, dR_To_dX
  public    ::    Transformation

  Type      ::    Transformation_Type
    logical                                     ::    Initialized   =   .False.   ! Indicator whether or not the object was initialized
    integer                                     ::    N                           ! Dimension of the matrix: N = 3*(NAtoms-1)
    real(rkp) ,dimension(:,:)     ,allocatable  ::    Tqp                         ! Transformation matrix : the order is x,y,z of atom 1, x,y,z of atom 2, etc. (DIM=nh,nh)
    real(rkp) ,dimension(:,:)     ,allocatable  ::    Tpq                         ! Inverse transformation matrix (DIM=nh,nh)
  contains
    private
    procedure     ,public   ::    Initialize             =>  Initialize_Transformation
  End Type

  type(Transformation_Type)               ::    Transformation
  
  logical   ,parameter    ::    i_Debug_Global = .True.!.False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Transformation( This, Mi, i_Debug )
! This procedures generates the transformation from the time derivative of the generalized coordinates to the conjugate momenta.
! The generalized coordinates are the cartesian coordinates of the first NAtom-1 atoms.
! (i.e.: From Velocities & Positions (PaQ) ---> Momenta & Positions)
! In case of 3 Atoms:
!
!       |                                                                                                            |
!       | (m3+m1)*m1/m3           0                0            m1*m2/m3             0                0              |
!       |                                                                                                            |
!       |        0         (m3+m1)*m1/m3           0                0            m1*m2/m3             0              |
!       |                                                                                                            |
!       |        0                0         (m3+m1)*m1/m3           0                0            m1*m2/m3           |
! Tqp = |                                                                                                            |
!       |    m1*m2/m3             0                0         (m3+m2)*m2/m3           0                0              |
!       |                                                                                                            |
!       |        0            m1*m2/m3             0                0         (m3+m2)*m2/m3           0              |
!       |                                                                                                            |
!       |        0                0            m1*m2/m3             0                0         (m3+m2)*m2/m3         |
!       |                                                                                                            |
!

  use Parameters_Module         ,only:  Zero, One
  use Matrix_Inversion_Module   ,only:  Inverse_Matrix

  class(Transformation_Type)              ,intent(out)    ::    This
  real(rkp) ,dimension(:)                 ,intent(in)     ::    Mi                      ! Masses of atoms
  logical                       ,optional ,intent(in)     ::    i_Debug

  logical                                                 ::    i_Debug_Loc

  integer                                                 ::    i, j, ic, ip
  integer                                                 ::    NAtoms                  ! Number of atoms
  integer                                                 ::    N                       ! Number of dof for the generalized coordinates or momenta
  real(rkp)                                               ::    xmi, temp, temp1, temp2 ! Inverse of the atom mass of last atom

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Transformation" )
  !i_Debug_Loc   =     Logger%On()
  

  NAtoms    =     size(Mi)                                                              ! Getting the number of atoms
  N         =     3 * (NAtoms-1)                                                        ! Getting the number of dof for the generalized coordinates or momenta: 6/9 for 3/4 atoms
  This%N    =     N

  allocate( This%Tqp(N,N) )
  allocate( This%Tpq(N,N) )

  This%Tqp    =   Zero

  xmi         =   One / Mi(NAtoms)
  ic          =   0

  do i = 1,NAtoms-1
    temp1     =   Mi(i) * xmi
    temp      =   ( Mi(NAtoms) + Mi(i) ) * temp1
    do j = 1,3
      ic            =   ic+1
      This%Tqp(ic,ic)    =   temp
      if ( i+1 > NAtoms-1 ) cycle
      do ip = i+1,NAtoms-1
        temp2                    =   temp1 * Mi(ip)
        This%Tqp(ic,(ip-1)*3+j)  =   temp2
        This%Tqp((ip-1)*3+j,ic)  =   temp2
      end do
    end do
  end do
  This%Tpq    =   Inverse_Matrix( This%Tqp )
  
  if (i_Debug_Loc) write(Logger%Unit,"(2x,'[Transformation]: Tqp = ',*(3x,es15.8))") This%Tqp
  if (i_Debug_Loc) write(Logger%Unit,"(2x,'[Transformation]: Tpq = ',*(3x,es15.8))") This%Tpq

  This%Initialized  =   .True.

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine X_to_R(X, R)
!**********************************************************************
!  Program to calculate the six interatomic distance 
!  by reading XYZ coordinate
!**********************************************************************
  
  real(rkp) ,dimension(:) ,intent(in)  :: X
  real(rkp) ,dimension(:) ,intent(out) :: R
  
  integer                              :: NAtoms, NPairs
  
  NAtoms = size(X,1) / 3
  NPairs = size(R,1) !NAtoms * (NAtoms - 1) / 2

  if (NAtoms == 3) then
  !**********************************************************************
  !  Now, calculate the inter-atomic distance
  !  r1 = r(AB)    
  !  r2 = r(AC)    
  !  r3 = r(BC)    
  !**********************************************************************

    R(1) = Sqrt( (  X(4) -  X(1) )**2 + (  X(5) -  X(2) )**2 + (  X(6) -  X(3) )**2)
    R(2) = Sqrt( (  X(7) -  X(1) )**2 + (  X(8) -  X(2) )**2 + (  X(9) -  X(3) )**2)
    R(3) = Sqrt( ( X(10) -  X(4) )**2 + ( X(11) -  X(5) )**2 + ( X(12) -  X(6) )**2)
 
  elseif (NAtoms == 4) then
  !**********************************************************************
  !  Now, calculate the inter-atomic distance
  !  r1 = r(AB)    r2 = r(AC)
  !  r3 = r(AD)    r4 = r(BC)
  !  r5 = r(BC)    r6 = r(CD)       
  !**********************************************************************

    R(1) = Sqrt( (  X(4) -  X(1) )**2 + (  X(5) -  X(2) )**2 + (  X(6) -  X(3) )**2)
    R(2) = Sqrt( (  X(7) -  X(1) )**2 + (  X(8) -  X(2) )**2 + (  X(9) -  X(3) )**2)
    R(3) = Sqrt( ( X(10) -  X(1) )**2 + ( X(11) -  X(2) )**2 + ( X(12) -  X(3) )**2)
    R(4) = Sqrt( (  X(4) -  X(7) )**2 + (  X(5) -  X(8) )**2 + (  X(6) -  X(9) )**2)
    R(5) = Sqrt( (  X(4) - X(10) )**2 + (  X(5) - X(11) )**2 + (  X(6) - X(12) )**2)
    R(6) = Sqrt( (  X(7) - X(10) )**2 + (  X(8) - X(11) )**2 + (  X(9) - X(12) )**2)
 
  end if

End Subroutine 
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine R_to_X(R, X, Theta)
!**********************************************************************
!  Program to calculate the interatomic distance by reading XYZ coordinate
!**********************************************************************
  
  real(rkp) ,dimension(:) ,intent(in)            :: R
  real(rkp) ,dimension(:) ,intent(out)           :: X
  real(rkp)               ,intent(in)  ,optional :: Theta

  integer                                        :: NAtoms, NPairs
  real(rkp)                                      :: Theta2

  NAtoms = size(X,1) / 3
  NPairs = size(R,1) !NAtoms * (NAtoms - 1) / 2

  X = Zero

  if (NAtoms == 3) then

    if (present(Theta)) then
      Theta2 = Theta
    else
      Theta2 = cos( (R(1)**2 + R(3)**2 - R(2)**2 ) / (Two*R(1)*R(3)) )
    end if

    X(1) = R(1)
    X(7) = R(3) * cos(Theta2)
    X(8) = R(3) * sin(Theta2)

  elseif (NAtoms == 4) then

    stop

  end if

End Subroutine 
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine dR_To_dX(R, X, dVdR, dRdX, dVdX)
!**********************************************************************
! Subroutine to evaluate dRdX for giving R and X 
! R:    R(:), 3/6 bond lengths
! X:    X(:), 6/12 Cartesian coordinates
! 
!**********************************************************************

  real(rkp) ,dimension(:)   ,intent(in)              :: R
  real(rkp) ,dimension(:)   ,intent(in)              :: X
  real(rkp) ,dimension(:)   ,intent(in)    ,optional :: dVdR
  real(rkp) ,dimension(:,:) ,intent(out)   ,optional :: dRdX
  real(rkp) ,dimension(:)   ,intent(out)   ,optional :: dVdX

  integer                                :: i, j
  integer                                :: NAtoms, NPairs
  logical                                :: ComputeDerFlg  = .false.

  if (present(dVdX)) then
    ComputeDerFlg = .true.
  end if

  NAtoms = size(X,1) / 3
  NPairs = size(R,1) !NAtoms * (NAtoms - 1) / 2

  ! Initialize dRdX
  dRdX = Zero

  if (NAtoms == 3) then
    ! Start to calculate the non-zero dRdX
    ! dr1dx
    dRdX(1,1)=(x(1)-x(4))/r(1)
    dRdX(1,2)=(x(2)-x(5))/r(1)
    dRdX(1,3)=(x(3)-x(6))/r(1)
    dRdX(1,4)=-dRdX(1,1)
    dRdX(1,5)=-dRdX(1,2)
    dRdX(1,6)=-dRdX(1,3)

    ! dr2dx
    dRdX(2,1)=(x(1)-x(7))/r(2)
    dRdX(2,2)=(x(2)-x(8))/r(2)
    dRdX(2,3)=(x(3)-x(9))/r(2)
    dRdX(2,7)=-dRdX(2,1)
    dRdX(2,8)=-dRdX(2,2)
    dRdX(2,9)=-dRdX(2,3)

    ! dr3dx
    dRdX(3,7)=(x(7)-x(4))/r(3)
    dRdX(3,8)=(x(8)-x(5))/r(3)
    dRdX(3,9)=(x(9)-x(6))/r(3)
    dRdX(3,4)=-dRdX(3,7)
    dRdX(3,5)=-dRdX(3,8)
    dRdX(3,6)=-dRdX(3,9)

  elseif (NAtoms == 4) then
    ! Start to calculate the non-zero dRdX
    ! dr1dx
    dRdX(1,1)=(x(1)-x(4))/r(1)
    dRdX(1,2)=(x(2)-x(5))/r(1)
    dRdX(1,3)=(x(3)-x(6))/r(1)
    dRdX(1,4)=-dRdX(1,1)
    dRdX(1,5)=-dRdX(1,2)
    dRdX(1,6)=-dRdX(1,3)

    ! dr2dx
    dRdX(2,1)=(x(1)-x(7))/r(2)
    dRdX(2,2)=(x(2)-x(8))/r(2)
    dRdX(2,3)=(x(3)-x(9))/r(2)
    dRdX(2,7)=-dRdX(2,1)
    dRdX(2,8)=-dRdX(2,2)
    dRdX(2,9)=-dRdX(2,3)

    ! dr3dx
    dRdX(3,1)=(x(1)-x(10))/r(3)
    dRdX(3,2)=(x(2)-x(11))/r(3)
    dRdX(3,3)=(x(3)-x(12))/r(3)
    dRdX(3,10)=-dRdX(3,1)
    dRdX(3,11)=-dRdX(3,2)
    dRdX(3,12)=-dRdX(3,3)
  
    ! dr4dx
    dRdX(4,4)=(x(4)-x(7))/r(4)
    dRdX(4,5)=(x(5)-x(8))/r(4)
    dRdX(4,6)=(x(6)-x(9))/r(4)
    dRdX(4,7)=-dRdX(4,4)
    dRdX(4,8)=-dRdX(4,5)
    dRdX(4,9)=-dRdX(4,6)

    ! dr5dx
    dRdX(5,4)=(x(4)-x(10))/r(5)
    dRdX(5,5)=(x(5)-x(11))/r(5)
    dRdX(5,6)=(x(6)-x(12))/r(5)
    dRdX(5,10)=-dRdX(5,4)
    dRdX(5,11)=-dRdX(5,5)
    dRdX(5,12)=-dRdX(5,6)

    ! dr6dx
    dRdX(6,7)=(x(7)-x(10))/r(6)
    dRdX(6,8)=(x(8)-x(11))/r(6)
    dRdX(6,9)=(x(9)-x(12))/r(6)
    dRdX(6,10)=-dRdX(6,7)
    dRdX(6,11)=-dRdX(6,8)
    dRdX(6,12)=-dRdX(6,9)
  end if
  ! Finish the calculation of non-zero dRdX

  if (ComputeDerFlg) then
    do i=1,NAtoms*3
      do j=1,NPairs
        dVdX(i) = dVdX(i) + dRdX(j,i) * dVdR(j)
      enddo
    enddo
  end if

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine dX_To_dR(R, X, dVdX, dXdR, dVdR)
!**********************************************************************
! Subroutine to evaluate dRdX for giving R and X 
! R:    R(:), 3/6 bond lengths
! X:    X(:), 6/12 Cartesian coordinates
! 
!**********************************************************************

  real(rkp) ,dimension(:)   ,intent(in)              :: R
  real(rkp) ,dimension(:)   ,intent(in)              :: X
  real(rkp) ,dimension(:)   ,intent(in)    ,optional :: dVdX
  real(rkp) ,dimension(:,:) ,intent(out)   ,optional :: dXdR
  real(rkp) ,dimension(:)   ,intent(out)   ,optional :: dVdR

  integer                                :: i, j
  integer                                :: NAtoms, NPairs
  logical                                :: ComputeDerFlg  = .false.

  if (present(dVdX)) then
    ComputeDerFlg = .true.
  end if

  NAtoms = size(X,1) / 3
  NPairs = size(R,1) !NAtoms * (NAtoms - 1) / 2

  ! Initialize dRdX
  dXdR = Zero

  if (NAtoms == 3) then
    ! Start to calculate the non-zero dRdX
    ! dr1dx
    dXdR(1,1) = R(1) / (x(1) - x(4))
    dXdR(2,1) = R(1) / (x(2) - x(5))
    dXdR(3,1) = R(1) / (x(3) - x(6))
    dXdR(4,1) = - dXdR(1,1)
    dXdR(5,1) = - dXdR(2,1)
    dXdR(6,1) = - dXdR(3,1)

    dXdR(1,2) = R(2) / (x(1) - x(7))
    dXdR(2,2) = R(2) / (x(2) - x(8))
    dXdR(3,2) = R(2) / (x(3) - x(9))
    dXdR(7,2) = - dXdR(1,2)
    dXdR(8,2) = - dXdR(2,2)
    dXdR(9,2) = - dXdR(3,2)

    dXdR(7,3) = R(3) / (x(7) - x(4))
    dXdR(8,3) = R(3) / (x(8) - x(5))
    dXdR(9,3) = R(3) / (x(9) - x(6))
    dXdR(4,3) = - dXdR(7,3)
    dXdR(5,3) = - dXdR(8,3)
    dXdR(6,3) = - dXdR(9,3)

  elseif (NAtoms == 4) then
    ! Start to calculate the non-zero dRdX
    ! dr1dx
    dXdR(1,1) = R(1) / (x(1) -  x(4))
    dXdR(2,1) = R(1) / (x(2) -  x(5))
    dXdR(3,1) = R(1) / (x(3) -  x(6))
    dXdR(1,2) = R(2) / (x(1) -  x(7))
    dXdR(2,2) = R(2) / (x(2) -  x(8))
    dXdR(3,2) = R(2) / (x(3) -  x(9))
    dXdR(1,3) = R(3) / (x(1) - x(10))
    dXdR(2,3) = R(3) / (x(2) - x(11))
    dXdR(3,3) = R(3) / (x(3) - x(12))

    dXdR(4,1) = - dXdR(1,1)
    dXdR(5,1) = - dXdR(2,1)
    dXdR(6,1) = - dXdR(3,1)
    dXdR(4,4) = R(4) / (x(4) -  x(7))
    dXdR(5,4) = R(4) / (x(5) -  x(8))
    dXdR(6,4) = R(4) / (x(6) -  x(9))
    dXdR(4,5) = R(5) / (x(4) - x(10))
    dXdR(5,5) = R(5) / (x(5) - x(11))
    dXdR(6,5) = R(5) / (x(6) - x(12))

    dXdR(7,2) = - dXdR(1,2)
    dXdR(8,2) = - dXdR(2,2)
    dXdR(9,2) = - dXdR(3,2)
    dXdR(7,4) = - dXdR(4,4)
    dXdR(8,4) = - dXdR(5,4)
    dXdR(9,4) = - dXdR(6,4)
    dXdR(7,6) = R(6) / (x(7) - x(10))
    dXdR(8,6) = R(6) / (x(8) - x(11))
    dXdR(9,6) = R(6) / (x(9) - x(12))

    dXdR(10,3) = - dXdR(1,3)
    dXdR(11,3) = - dXdR(2,3)
    dXdR(12,3) = - dXdR(3,3)
    dXdR(10,5) = - dXdR(4,5)
    dXdR(11,5) = - dXdR(5,5)
    dXdR(12,5) = - dXdR(6,5)
    dXdR(10,6) = - dXdR(7,6)
    dXdR(11,6) = - dXdR(8,6)
    dXdR(12,6) = - dXdR(9,6)
  end if
  ! Finish the calculation of non-zero dRdX

  if (ComputeDerFlg) then
    do j=1,NPairs
     do i=1,NAtoms*3
        dVdR(j) = dVdR(j) + dXdR(i,j) * dVdX(i)
      enddo
    enddo
  end if

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module