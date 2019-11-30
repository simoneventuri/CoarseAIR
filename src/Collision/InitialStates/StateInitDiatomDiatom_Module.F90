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
Module StateInitDiatomDiatom_Module

  use Parameters_Module             ,only:  rkp
  use Logger_Class                  ,only:  Logger
  use Error_Class                   ,only:  Error
  use StateInitDiatomAtom_Module    ,only:  ComputeKineticEnergy, BondLength, Compute_Coordinates_Velocities

  implicit none

  private
  public    ::    InitializeLevels_DiatomDiatom
  public    ::    SetInitialState_DiatomDiatom

  logical   ,parameter    ::    i_Debug_Global = .True.!.False.

  contains

! Read in info for target and projectile initial states.
! The target is ab and the projectile is cd, i.e. diatom+atom collisions.
! This procedure was initially called preini4
Subroutine InitializeLevels_DiatomDiatom( Input, Species, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Species_Class         ,only:  Species_Type

  type(Input_Type)                          ,intent(in)     ::    Input
  type(Species_Type)  ,dimension(:)         ,intent(inout)  ::    Species       ! Species objects which always have 2 elements [Target,Projectile]. Since we are in the case of a Diatom-Diatom collision, both elements are used, which corresponds to the Target species. Dim=(NSpecies)=(2)
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iMol 
  logical                                                   ::    i_Debug_Loc
  character(2)                                              ::    iMol_char
  character(:) ,allocatable                                 ::    FileName
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeLevels_DiatomDiatom")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  do iMol = 1,size(Species)

     Write(iMol_char,'(I2)')iMol
     iMol_char = adjustl(iMol_char)

     FileName = './levels_' // trim(adjustl(Input%Molecules_Name(iMol)))// trim(iMol_char) // '.inp'

     if (i_Debug_Loc) write(Logger%Unit,"(6x,'[InitializeLevels_DiatomDiatom]: Calling Species(',a,')%ListStates%Initialize with FileName = ',a)")trim(iMol_char),FileName
     call Species(iMol)%ListStates%Initialize( Input, iMol, FileName=FileName, i_Debug=i_Debug_Loc )                                                  
     if (i_Debug_Loc) write(Logger%Unit,"(6x,'[InitializeLevels_DiatomDiatom]: -> Done initializing the Species(',a,')%ListStates object')")trim(iMol_char)

  enddo

  if (i_Debug_Loc) call Logger%Exiting
 
End Subroutine

!________________________________________________________________________________________________________________________________!
Subroutine SetInitialState_DiatomDiatom( Species, Qijk, dQijk, iTraj, iPES, i_Debug, i_Debug_Deep )

  use Species_Class          ,only:  Species_Type
  use Parameters_Module      ,only:  Zero, One, Two, Half
  use RandomVector_Module    ,only:  RanddwsVec
  use RandomVectorBis_Module ,only:  RanddwsVecBis
  use Level_Class            ,only:  Level_Type

  type(Species_Type)  ,dimension(:)         ,intent(in)     ::    Species       ! Species objects which always have 2 elements [Target,Projectile]. Since we are in the case of a Diatom-Diatom collision, both elements are used. Dim=(NSpecies)=(2)
  real(rkp)           ,dimension(:,:,:)     ,intent(out)    ::    Qijk          ! Coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  real(rkp)           ,dimension(:,:,:)     ,intent(out)    ::    dQijk         ! Time derivatives of the coordinates for each spatial direction (dim-1), for each atom in the species (dim-2) and for each species (dim-3: 1=target, 2:projectile). Dim=(NSpace,NAtoMaxSpe,NSpecies)=(3,2,2)
  integer                                   ,intent(in)     ::    iTraj
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    i, iMol
  real(rkp)                                                 ::    RandNum        ! Random number
  real(rkp)                                                 ::    Omega         ! Angular velocity for rotation
  real(rkp)                                                 ::    Vc_R2         ! Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                                 ::    Ekin
  real(rkp)                                                 ::    BL            ! bond length
  real(rkp)                                                 ::    dBL           ! Derivative of the ound length
  type(Level_Type)                                          ::    State
  integer                                                   ::    Traj_iState, Traj_jqn, Traj_vqn

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug  ! This line was missing
 ! ======================================================================================================================
 !     LOOP OVER MOLECULES
 ! ======================================================================================================================
  do iMol = 1,Size(Species)
  

! ======================================================================================================================
!     SELECTING THE INTERNAL STATE (Traj_iState, Traj_jqn and Traj_vqn) AND COMPUTING THE CENTRIFUGAL POTENTIAL (Vc_R2)
! ======================================================================================================================
     if ( .not. (    ( ( trim(adjustl(Species(iMol)%ListStates%vInMethod)) .eq. "Fixed"        ) .and. ( trim(adjustl(Species(iMol)%ListStates%jInMethod)) .eq. "Fixed" ) ) &
               .or. ( trim(adjustl(Species(iMol)%ListStates%vInMethod)) .eq. "MostProbable" ) ) ) then     ! Case of sampled initial internal states
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Using a Thermal Distribution for Choosing Initial State')")

        RandNum = RanddwsVec(iPES)
        do i = 1,Species(iMol)%ListStates%NStates
           if ( Species(iMol)%ListStates%States(i)%rlim >= RandNum ) exit
        end do
        Traj_iState   =   i
        Traj_jqn      =   Species(iMol)%ListStates%States(i)%jqn
     else     
        if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Either Fixed State or State that Maximizes Boltzmann Distribution')")
        Traj_iState   =   Species(iMol)%ListStates%StateIn
        Traj_jqn      =   Species(iMol)%ListStates%jIn
     end if
     Traj_vqn        =   Species(iMol)%ListStates%States(i)%vqn
     if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Species(iMol)%DiaPot%xmui2 = ',es15.8,'; Traj_jqn = ',g0)") Species(iMol)%DiaPot%xmui2, Traj_jqn
     Vc_R2           =   Species(iMol)%DiaPot%xmui2 * ( Traj_jqn + Half )**2
     associate( S => Species(iMol)%ListStates%States(Traj_iState) )
       State%rmin    =   S%rmin
       State%tau     =   S%tau
       State%Eint    =   S%Eint
       State%ri      =   S%ri
       State%ro      =   S%ro
       State%egam    =   S%egam
       State%rmax    =   S%rmax
       State%vmin    =   S%vmin
       State%vmax    =   S%vmax
       State%rlim    =   Zero      ! NOT USED
     end associate
   
    ! ! ==============================================================================================================
    ! !   RE-COMPUTING THE STATE PROPERTIES
    ! ! ==============================================================================================================
    ! !!! SIMONE added this subroutine. For some (~20) of the N2 Quasi-Bound levels from NASA Ames (e.g., 7489 from n2.leroy.Sort.dat), 
    ! ! the difference between the maximum values of the diatomic potential from the Levels list and the one computed in the "period" 
    ! ! subroutine was of the order of 1.e-5; this fact was making CoarseAIR not entering the if "(Ekin .ge. (State%vmax - Emin)) )" 
    ! ! statement in the "ComputeKineticEnergy" subroutine, and finally crashing with the error: "NaN found in scalar variable"
    ! ! in "* Procedure: Integrate_ODE_1d", "* Variable:  h before mmidv".

    ! if ( State%egam > Zero ) then
    !    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Found a Quasi Bound Level. Calling RecomputeLevelProperties')")
    !    call Species(iMol)%DiaPot%RecomputeLevelProperties( State%Eint, Vc_R2, State%rMin, State%VMin, State%rMax, State%VMax, State%ri, State%ro, i_Debug=i_Debug_Loc )
    ! end if
  
    ! if (i_Debug_Loc) then
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Internal state:                   Traj_iState = ',g0)")     Traj_iState
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Rotational quantum number:        Traj_jqn    = ',g0)")     Traj_jqn
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Vibrational quantum number:       Traj_vqn    = ',g0)")     Traj_vqn
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Centrifugal potential times R2:   Vc_R2       = ',es15.8)") Vc_R2
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%rmin  = ',es15.8)") State%rmin
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%tau   = ',es15.8)") State%tau
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%Eint  = ',es15.8)") State%Eint
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%ri    = ',es15.8)") State%ri
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%ro    = ',es15.8)") State%ro
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%egam  = ',es15.8)") State%egam
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%rmax  = ',es15.8)") State%rmax
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%vmin  = ',es15.8)") State%vmin
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%vmax  = ',es15.8)") State%vmax
    !    write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: ->                                   State%rlim  = ',es15.8)") State%rlim
    ! end if
    ! ! ==============================================================================================================


    ! ==============================================================================================================
    !   COMPUTING THE KINETIC ENERGY AND UPDATING SOME COMPONENT IN THE STATE OBJECT
    ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Calling ComputeKineticEnergy')")
    call ComputeKineticEnergy( Ekin, State, Vc_R2, Species(iMol)%DiaPot, iPES, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) then
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Kinetic energy:                   Ekin      = ',es15.8)") Ekin
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Estimate of outter turning point: State%ro  = ',es15.8)") State%ro
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Period:                           State%Tau = ',es15.8)") State%Tau
    end if
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   COMPUTING THE INITIAL BOND LENGTH AND ITS DERIVATIVES
    ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Computing the bond length and its derivatives. Calling BondLength')")
    Ekin = Zero                                                                                                                     ! Setting the kinetic energy at the potential minimum
    call BondLength( BL, dBL, State%ro, State%tau, Ekin, Vc_R2, Species(iMol)%DiaPot, iPES, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Deep )! Computng BL/bBL output, all the rest is input
    if (i_Debug_Loc) then
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Bond length:        BL       = ',es15.8)") BL
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Bond length Deriv.: dBL      = ',es15.8)") dBL
    end if
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   COMPUTING THE MOLECULE COORDINATES AND VELOCITIES
    ! ==============================================================================================================
    ! choose random variables for diatomic orientation; set some variables defining the angular momentum vector Omega
    ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Computing the diatomic coordinates')")
    ! Shouldn't it be:    Omega = sqrt(Two * Vc_R2 * xmui2) ??? (SIMONE)
    Omega = Two * sqrt( Vc_R2 * Species(iMol)%DiaPot%xmui2 ) / BL**2                                                                                         ! Setting the angular velocity for rotation
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Angular velocity: Omega = ',d20.10)")  Omega
                                                                                        
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Calling Compute_Coordinates_Velocities')")
    call Compute_Coordinates_Velocities( Species(iMol)%GetAtomsMass(), BL, dBL, Omega, Qijk(:,:,iMol), dQijk(:,:,iMol), iTraj, iPES, i_Debug=i_Debug_Loc )! Computing the coordinates of the target
    if (i_Debug_Loc) then
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Positions of atoms in target:')")
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]:    -> Atom A: ',*(d20.10,3x))")           Qijk(:,1,iMol)
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]:    -> Atom B: ',*(d20.10,3x))")           Qijk(:,2,iMol)
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: -> Velocities of atoms in target:')")
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]:    -> Atom A: ',*(d20.10,3x))")           dQijk(:,1,iMol)
       write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]:    -> Atom B: ',*(d20.10,3x))")           dQijk(:,2,iMol)
    end if
! ==============================================================================================================

  enddo

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine

End Module
