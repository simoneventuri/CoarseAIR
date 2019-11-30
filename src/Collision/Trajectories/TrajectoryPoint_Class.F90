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

Module TrajectoryPoint_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use MolecularState_Class  ,only:  MolecularState_Type

  implicit none

  private
  public    ::    TrajectoryPoint_Type
  public    ::    ComputeCoordAndVeloc

  Type      ::    TrajectoryPoint_Type
    integer                                 ::    NEqtTot     =   0     !< Number of equations. Size of the component 'PaQ'
    integer                                 ::    NPairs      =   0     !< Number of atom-atom pairs. Size of the components 'Rbs' and iRbs'
    integer                                 ::    NMolecules  =   0     !< Effective number of molecules. For NAtoms=3, NMolecules<=1 and for NAtoms=4, NMolecules<=2. Thus, maximum number of molecules is 2.
    real(rkp)                               ::    t           =   Zero  !< Time of trajectory point
    real(rkp)                               ::    H           =   Zero  !< Hamiltonian of trajectory point
    real(rkp) ,dimension(:) ,allocatable    ::    PaQ                   !< Solution vector of trajectory point. Dim=(NEqtTot)
    real(rkp) ,dimension(:) ,allocatable    ::    Rbs                   !< Distances of each atom-atom pairs sorted in increasing values. Dim=(NPairs)
    integer   ,dimension(:) ,allocatable    ::    iRbs                  !< Index mapping from the sorted atom-atom distances to associated elements in the vector of Pair object (stored in the Collision object). Dim=(NPairs)
    real(rkp)                               ::    b           =   Zero  !< Impact parameter
    real(rkp)                               ::    xkin        =   Zero  !< Energy
    real(rkp) ,dimension(3)                 ::    Vrel        =   Zero  !< Relative velocity components
    type(MolecularState_Type) ,dimension(2) ::    Molecules             !< Array of MolecularState objects. Note: Harded-codded dimensions to avoid allocation for each trajectory. Number of effective molecuels is given by component NMolecules.
  contains
    private
    procedure ,public   ::    Initialize  =>  InitializeTrajectoryPoint
    procedure ,public   ::    Analyze     =>  AnalyzeTrajectoryPoint
    procedure ,public   ::    SetData     =>  SetDataTrajectoryPoint
  End Type

  logical   ,parameter    ::    i_Debug_Global = .True.!.False.
  integer   ,parameter    ::    NSpace = 3

  contains


!________________________________________________________________________________________________________________________________!
Subroutine InitializeTrajectoryPoint( This, Input, Collision, i_Debug )

  use Input_Class               ,only:  Input_Type
  use Collision_Class           ,only:  Collision_Type
  use Parameters_Module         ,only:  Zero
  use MolecularStateParam_Class ,only:  MolecularStateParam_Type

  class(TrajectoryPoint_Type)               ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    iM
  type(MolecularStateParam_Type)                            ::    Param

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeTrajectoryPoint")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  This%NEqtTot    =     Collision%NEqtTot
  This%NPairs     =     Collision%NPairs
  allocate( This%PaQ  (This%NEqtTot) ); This%PaQ   =   Zero
  allocate( This%Rbs  (This%NPairs)  ); This%Rbs   =   Zero
  allocate( This%iRbs (This%NPairs)  ); This%iRbs  =   0


! ==============================================================================================================
!     SETTING VARIABLES FROM THE INPUT OBJECT
! ==============================================================================================================
  if (i_Debug_Loc)  call Logger%Write( "Setting the parameters of the MolecularState objects")
  if (i_Debug_Loc)  call Logger%Write( "-> Calling Param%Initialize")
  call Param%Initialize( Input, Collision, i_Debug )
  do iM = 1,size(This%Molecules)
    This%Molecules(iM)%Param = Param
  end do
  if (i_Debug_Loc)  call Logger%Write( "-> Done with Param%Initialize")
! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetDataTrajectoryPoint( This, t, H, PaQ )

  class(TrajectoryPoint_Type)               ,intent(inout)  ::    This
  real(rkp)                                 ,intent(in)     ::    t           !< Time of trajectory point
  real(rkp)                                 ,intent(in)     ::    H           !< Hamiltonian of trajectory point
  real(rkp) ,dimension(:)                   ,intent(in)     ::    PaQ         !< Solution vector of trajectory point. Dim=(NEqtTot)
  
  This%t    =   t
  This%H    =   H
  This%PaQ  =   PaQ
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine AnalyzeTrajectoryPoint( This, Collision, i_Debug )
! This procedure analyzes the coordinates of atoms
! From 6*n-6 coordinates and velocities with cm stationary at the origin, calculate 6*n coordinates and velocities.
! Also compute impact parameter and load coordinates for analysis.
! Finaly, compute initial relative translational energy.
! The following variables are set in the 'This' object:
! * b
! * xkin:
! * Vrel:

  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero, One, Half
  use Sorting_Module        ,only:  hpsort

  class(TrajectoryPoint_Type)               ,intent(inout)  ::    This
  type(Collision_Type)                      ,intent(in)     ::    Collision
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    j, i
  integer                                                   ::    iA          ! Index of atoms: From 1 to NAtoms
  integer                                                   ::    iS          ! Index of space: From 1 to NSpace
  integer                                                   ::    iP          ! Index of Pairs
  integer                                                   ::    ierr
  integer                                                   ::    iM          ! index of molecules
  real(rkp)                                                 ::    xt, xp, xti, xpi, xnum, v, t
  real(rkp)                                                 ::    xjx, xjy, xjz, xjj
  real(rkp)                                                 ::    dist
  real(rkp) ,dimension(NSpace,Collision%NAtoms)             ::    xx          ! Dim=(NSpace,NAtoms)
  real(rkp) ,dimension(NSpace,Collision%NAtoms)             ::    xxdot       ! Dim=(NSpace,NAtoms)
  real(rkp) ,dimension(NSpace,2)                            ::    cm          ! Dim=(NSpace,2)
  real(rkp) ,dimension(NSpace,2)                            ::    cmdot       ! Dim=(NSpace,2)
  real(rkp) ,dimension(NSpace)                              ::    xrel        ! Dim=(NSpace)

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AnalyzeTrajectoryPoint" )
  !i_Debug_Loc   =     Logger%On()
  

!  if (i_Debug_Loc)  then
!    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Inputs')")
!    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> This%t   = ',es15.8)") This%t
!    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> This%H   = ',es15.8)") This%H
!    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> This%PaQ = ',*(es15.8,3x))") This%PaQ
!  end if


! ==============================================================================================================
!     COMPUTING THE 3 COORDINATES AND 3 VELOCITIES ASSOCIATED TO ALL ATOMS
! ==============================================================================================================
  if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Computing the 3 coordinates and 3 velocities associated to all atoms. Calling ComputeCoordAndVeloc')")
  call ComputeCoordAndVeloc( This, Collision, xx, xxdot )
  if (i_Debug_Loc)  then
    do iA = 1,Collision%NAtoms
      call Logger%Write( "xx(:,", iA, ")    = ", xx(:,iA), Fr="*(es15.8,3x)") 
      call Logger%Write( "xxdot(:,", iA, ") = ", xxdot(:,iA), Fr="*(es15.8,3x)")
    end do
  end if
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING BOND LENGTH AND SORTING THEM IN INCREASING ORDER
! ==============================================================================================================
  if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Computing bond length and sorting them in increasing order')")

  select case(Collision%NAtoms)

    ! 3 Atoms
    case(3)
      iP = 0
      do i = 1,Collision%NAtoms-1
         do j = i+1,Collision%NAtoms
            iP           = iP + 1
            This%Rbs(iP) = ( xx(1,i) - xx(1,j) )**2 + ( xx(2,i) - xx(2,j) )**2 + ( xx(3,i) - xx(3,j) )**2
         enddo
      enddo

    ! >=4 Atoms
    case default
      do iP = 1,Collision%NPairs
         i = Collision%Pairs(iP)%To_Atoms(1)
         j = Collision%Pairs(iP)%To_Atoms(2)
         This%Rbs(iP) = ( xx(1,i) - xx(1,j) )**2 + ( xx(2,i) - xx(2,j) )**2 + ( xx(3,i) - xx(3,j) )**2
     enddo 

  end select

  if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Calling hpsort')")
  call hpsort( This%Rbs, This%iRbs )
  if (i_Debug_Loc)  then
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Sorted bond lengths:')")
    do iP = 1,Collision%NPairs
      write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> iP = ',i3,3x,'This%Rbs =',es15.8,3x,'This%iRbs = ',i3)") iP, This%Rbs(iP), This%iRbs(iP)
    end do
  end if
! ==============================================================================================================


! ==============================================================================================================
!     SELECTING THE PAIR(S) AND ITS 2/4 ATOMS
! ==============================================================================================================
  if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Assembling pairs of bound atoms. Calling SetPairs')")
  call SetPairs( This, Collision, i_Debug=i_Debug_Loc )
  if (i_Debug_Loc)  then
    do iS = 1,2
      associate( Molecule => This%Molecules(iS) )
        write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> iS = ',i3,3x,'Molecule%iPair = ',i2,3x,'Molecule%iAtoms =',2(i3,3x))") iS, Molecule%iPair, Molecule%iAtoms(:)
      end associate
    end do
  end if
! ==============================================================================================================


! ==============================================================================================================
!     SELECTING THE POSITIONS/VELOCITIES OF THE TWO ATOMS ASSOCIATED TO EACH MOLECULE
! ==============================================================================================================
  if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Selecting the positions/velocities of the two atoms associated to each molecule')")
  do iS = 1,This%NMolecules    
    associate( Molecule => This%Molecules(iS) )
      if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> Molecule%iPair = ',g0)") Molecule%iPair
      do iA = 1,2
        Molecule%x(:,iA)    = xx(:,Molecule%iAtoms(iA))
        Molecule%xdot(:,iA) = xxdot(:,Molecule%iAtoms(iA))
      end do
    end associate
  end do
! ==============================================================================================================


! ==============================================================================================================
!     COMPUTING THE POSITIONS AND VELOCITIES OF THE CENTER OF MASS
! ==============================================================================================================
  cm    = Zero
  cmdot = Zero
  do i = 1,NSpace
    do j = 1,Collision%NAtoms
      if ( j == This%Molecules(1)%iAtoms(1) .or. j == This%Molecules(1)%iAtoms(2) ) then
        cm(i,1)    = cm(i,1)    + Collision%Atoms(j)%Mass * xx(i,j)
        cmdot(i,1) = cmdot(i,1) + Collision%Atoms(j)%Mass * xxdot(i,j)
      else
        cm(i,2)    = cm(i,2)    + Collision%Atoms(j)%Mass * xx(i,j)
        cmdot(i,2) = cmdot(i,2) + Collision%Atoms(j)%Mass * xxdot(i,j)
      end if
    end do
  end do
! ==============================================================================================================

  xt     =   Zero
  xp     =   Zero
  do i = 1,Collision%NAtoms
    if ( (i == This%Molecules(1)%iAtoms(1)) .or. (i == This%Molecules(1)%iAtoms(2)) ) then                                        ! Shouldn't it be xp->xt and xt->xp??? (Simone)
      xp = xp + Collision%Atoms(i)%Mass
    else
      xt = xt + Collision%Atoms(i)%Mass
    end if
  end do
  xpi    = One / xp
  xti    = One / xt

  do i = 1,NSpace
    cm(i,1)    = cm(i,1)    * xpi
    cm(i,2)    = cm(i,2)    * xti
    cmdot(i,1) = cmdot(i,1) * xpi
    cmdot(i,2) = cmdot(i,2) * xti
  end do
  if (i_Debug_Loc)  then
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> cm(:,1)    = ',3(d20.10,3x))") cm(:,1)
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> cm(:,2)    = ',3(d20.10,3x))") cm(:,2)
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> cmdot(:,1) = ',3(d20.10,3x))") cmdot(:,1)
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> cmdot(:,2) = ',3(d20.10,3x))") cmdot(:,2)
  end if

  do i = 1,NSpace
    xrel(i)      = cm(i,1)    - cm(i,2)
    This%Vrel(i) = cmdot(i,1) - cmdot(i,2)
  end do
  
  if (i_Debug_Loc)  then
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> xrel = ',3(d20.10,3x))") xrel(:)
    write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> Vrel = ',3(d20.10,3x))") This%Vrel(:)
  end if

  xnum = xrel(1)*This%Vrel(1) + xrel(2)*This%Vrel(2) + xrel(3)*This%Vrel(3)
  v    = This%Vrel(1)**2 + This%Vrel(2)**2 + This%Vrel(3)**2
  t    = - xnum / v

  This%b   =   sqrt( ( xrel(1) + t * This%Vrel(1) )**2 + ( xrel(2) + t * This%Vrel(2) )**2 + ( xrel(3) + t * This%Vrel(3) )**2  )

  This%xkin   = Zero
  do i = 1,NSpace
    This%xkin = This%xkin + Half * xp * cmdot(i,1)**2 + Half * xt * cmdot(i,2)**2
  end do

! #################### note cross product signs may not be correct ####################
  xjx   =   ( cm(2,1) * cmdot(3,1) - cm(3,1) * cmdot(2,1) ) * xp + ( cm(2,2) * cmdot(3,2) - cm(3,2) * cmdot(2,2) ) * xp           ! Shouldn't it be xp and xt??? (Simone)
  xjy   =   ( cm(1,1) * cmdot(3,1) - cm(3,1) * cmdot(1,1) ) * xp + ( cm(1,2) * cmdot(3,2) - cm(3,2) * cmdot(1,2) ) * xp
  xjz   =   ( cm(1,1) * cmdot(2,1) - cm(2,1) * cmdot(1,1) ) * xp + ( cm(1,2) * cmdot(2,2) - cm(2,2) * cmdot(1,2) ) * xp
! #################### note cross product signs may not be correct ####################

  xjj   =   sqrt( xjx**2 + xjy**2 + xjz**2 )
  dist  =   sqrt( xrel(1)**2 + xrel(2)**2 + xrel(3)**2 )
  v     =   sqrt(v)


  if (i_Debug_Loc)  then
    do iM = 1,2
      associate( Molecule => This%Molecules(iM) )
        if (i_Debug_Loc) write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: iM = ',i3,3x,'Molecule%iPair = ',i2,3x,'Molecule%iAtoms =',2(i3,3x))") iM, Molecule%iPair, Molecule%iAtoms(:)
      end associate
    end do
  end if

  if ( This%Molecules(1)%iPair < 0 ) then
    do iM = 1,2
      associate( Molecule => This%Molecules(iM) )
        Molecule%itype    =   10 * Molecule%iPair                                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: CHECK 10 or 12
        Molecule%viba     =   -One
        Molecule%AngMom   =   -One
      end associate
    end do
    return
  end if


! ==============================================================================================================
!     IDENTIFYING THE INTERNAL STATES ASSOCIATED TO EACH MOLECULE
! ==============================================================================================================
  if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: Identifying the internal states associated to each molecule. Number of molecules: This%NMolecules = ',g0)") This%NMolecules
  do iM = 1,This%NMolecules
    if (i_Debug_Loc)  write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: iM = ',i3,3x,'Calling This%Molecules(iM)%FindState')") iM
    call This%Molecules(iM)%FindState( Collision, ierr, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc)  then
      associate( Molecule => This%Molecules(iM) )
        write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> Molecule%viba   = ',f6.2)")   Molecule%viba
        write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> Molecule%AngMom = ',f6.2)")   Molecule%AngMom
        write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> Molecule%Eint   = ',es15.8)") Molecule%Eint
        write(Logger%Unit,"(8x,'[AnalyzeTrajectoryPoint]: -> Molecule%itype  = ',g0)")     Molecule%itype
      end associate
    end if
  end do
! ==============================================================================================================

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine ComputeCoordAndVeloc( This, Collision, xx, xxdot )
! From 6*n-6 coordinates and velocities with cm stationary at the origin, calculate 6*n coordinates and velocities.

  use Collision_Class       ,only:  Collision_Type
  use Parameters_Module     ,only:  Zero

  type(TrajectoryPoint_Type)                ,intent(in)     ::    This
  type(Collision_Type)                      ,intent(in)     ::    Collision
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    xx          ! Dim=(NSpace,NAtoms)
  real(rkp) ,dimension(:,:)                 ,intent(out)    ::    xxdot       ! Dim=(NSpace,NAtoms)

  integer                                                   ::    N           ! Number of atoms
  integer                                                   ::    iA          ! Index of atoms: From 1 to NAtoms
  integer                                                   ::    iS          ! Index of space: From 1 to NSpace

  N                 =   Collision%NAtoms

  do iA = 1,N-1
    do iS = 1,NSpace
      xxdot(iS,iA)  =   This%PaQ(iS+(iA-1)*3)
      xx(iS,iA)     =   This%PaQ(iS+(iA-1)*3+Collision%NEqtVar)
    end do
  end do

  do iS = 1,NSpace
    xxdot(iS,N)     =   Zero
    xx(iS,N)        =   Zero
    do iA = 1,N-1
      xxdot(iS,N)   =   xxdot(iS,N) + Collision%mMiMn(iA) * This%PaQ(iS+(iA-1)*3)
      xx(iS,N)      =   xx(iS,N)    + Collision%mMiMn(iA) * This%PaQ(iS+(iA-1)*3+Collision%NEqtVar)
    end do
  end do

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine SetPairs( This, Collision, i_Debug )

  use Collision_Class       ,only:  Collision_Type

  type(TrajectoryPoint_Type)                ,intent(inout)  ::    This
  type(Collision_Type)                      ,intent(in)     ::    Collision
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    iP          ! Index of Pairs
  integer                                                   ::    iS          ! Index of species

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "SetPairs" )
  !i_Debug_Loc   =     Logger%On()


  This%Molecules(:)%iPair        =   -1

  if ( Collision%NAtoms == 3 ) then

    iS                            =   1
    iP                            =   This%iRbs(1)
    This%Molecules(iS)%iPair      =   iP
    This%Molecules(iS)%iAtoms(:)  =   Collision%Pairs(iP)%To_Atoms(:)

    iS                            =   2
    This%Molecules(iS)%iPair      =   0
    This%Molecules(iS)%iAtoms(:)  =   0

  else if ( Collision%NAtoms == 4 ) then

    if ( This%Rbs(2) > 100.0_rkp ) then
      select case ( This%iRbs(1) )
        case (:3); iP  =   This%iRbs(1)
        case ( 4); iP  =   3
        case ( 5); iP  =   2
        case ( 6); iP  =   1
      end select
    else if ( This%Rbs(3) < 30.0_rkp ) then      ! If the atom-atom distance between 2 atoms from different species is smaller than a predifined threshold, then a cluster has been formed
      write(Logger%Unit,"(10x,'[SetPairs]: Cluster formed')")
      write(Logger%Unit,"(10x,'[SetPairs]: -> Sorted bond lengths:')")
      do iP = 1,Collision%NPairs
        write(Logger%Unit,"(10x,'[SetPairs]: iP = ',i3,3x,'This%Rbs =',es15.8,3x,'This%iRbs = ',i3)") iP, This%Rbs(iP), This%iRbs(iP)
      end do
      return
    else
      iP           = min( This%iRbs(1), This%iRbs(2) )
      if (iP>3) iP = 7 - iP  !MARCO
    end if

    iS                           = 1
    This%Molecules(iS)%iPair     = iP
    This%Molecules(iS)%iAtoms(:) = Collision%Pairs(iP)%To_Atoms(:)

    iS                           = 2
    iP                           = Collision%Pairs(iP)%Opposite
    This%Molecules(iS)%iPair     = iP
    This%Molecules(iS)%iAtoms(:) = Collision%Pairs(iP)%To_Atoms(:)

  end if

  This%NMolecules    =   0
  do iS = 1,size(This%Molecules)
    if ( This%Molecules(iS)%iPair /= 0 ) This%NMolecules = This%NMolecules + 1
  end do
 
  if (i_Debug_Loc)  then
    write(Logger%Unit,"(10x,'[SetPairs]: Number of molecules: This%NMolecules = ',g0)") This%NMolecules
    do iS = 1,2
    associate( State => This%Molecules(iS) )
      write(Logger%Unit,"(10x,'[SetPairs]: iS = ',i3,3x,'State%iPair = ',i2,3x,'State%iAtoms =',2(i3,3x))") iS, State%iPair, State%iAtoms(:)
    end associate
    end do
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
