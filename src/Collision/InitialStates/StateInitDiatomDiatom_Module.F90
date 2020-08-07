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
  use DiatomicPotential_Class       ,only:  DiatomicPotential_Type
  use File_Class                    ,only:  File_Type
  use Timing_Module
  use StateInitDiatomAtom_Module    ,only:  ComputeKineticEnergy, BondLength, Compute_Coordinates_Velocities

  implicit none

  private
  public    ::    InitializeLevels_DiatomDiatom
  public    ::    SetInitialState_DiatomDiatom

  type(File_Type)                         ::    ParamsFile
  real(rkp)       ,dimension(19)          ::    Params
  real(rkp)                               ::    xmui_Global
  real(rkp)                               ::    Vc_R2_Global
  class(DiatomicPotential_Type) ,pointer  ::    DiatPot_Global                            !< Local pointer for the intra-nuclear diatomic potential object used in procedure 'FuncRHS'
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

  integer                                                   ::    iMol, iSpec
  logical                                                   ::    i_Debug_Loc
  character(2)                                              ::    iMol_char
  character(:) ,allocatable                                 ::    FileName
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "InitializeLevels_DiatomDiatom")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  do iSpec = 1,size(Species)
    iMol = Species(iSpec)%To_Molecule
    write(iMol_char,'(I2)') iMol
    iMol_char = adjustl(iMol_char)

    FileName = './levels_' // trim(adjustl(Input%Molecules_Name(iMol)))// trim(iMol_char) // '.inp'

    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[InitializeLevels_DiatomDiatom]: Calling Species(',a,')%ListStates%Initialize with FileName = ',a)")trim(iMol_char),FileName
    call Species(iSpec)%ListStates%Initialize( Input, Species(iSpec)%DiatPot, iMol, iSpec, FileName=FileName, i_Debug=i_Debug_Loc )                                                  
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[InitializeLevels_DiatomDiatom]: -> Done initializing the Species(',a,')%ListStates object')") trim(iMol_char)

  enddo

  if (i_Debug_Loc) call Logger%Exiting
 
End Subroutine

!________________________________________________________________________________________________________________________________!
Subroutine SetInitialState_DiatomDiatom( Species, Qijk, dQijk, iTraj, iPES, NTrajOverall, iProc, i_Debug, i_Debug_Deep )

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
  integer                                   ,intent(in)     ::    NTrajOverall
  integer                                   ,intent(in)     ::    iProc
  logical                         ,optional ,intent(in)     ::    i_Debug
  logical                         ,optional ,intent(in)     ::    i_Debug_Deep

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    i, iMol, iSpec
  real(rkp)                                                 ::    RandNum        ! Random number
  real(rkp)                                                 ::    Omega         ! Angular velocity for rotation
  real(rkp)                                                 ::    Vc_R2         ! Centrifual potential multiplied by r**2 [hartree.bohr^2]
  real(rkp)                                                 ::    Ekin
  real(rkp)                                                 ::    BL            ! bond length
  real(rkp)                                                 ::    dBL           ! Derivative of the ound length
  type(Level_Type)                                          ::    State
  integer                                                   ::    Traj_iState, Traj_jqn, Traj_vqn

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug  ! This line was missing
  if (i_Debug_Loc) call Logger%Entering( "SetInitialState_DiatomDiatom" )
 ! ======================================================================================================================
 !     LOOP OVER MOLECULES
 ! ======================================================================================================================
  do iSpec = 1,2
    iMol = Species(iSpec)%To_Molecule

! ======================================================================================================================
!     SELECTING THE INTERNAL STATE (Traj_iState, Traj_jqn and Traj_vqn) AND COMPUTING THE CENTRIFUGAL POTENTIAL (Vc_R2)
! ======================================================================================================================
    if ( .not. ( trim(adjustl(Species(iSpec)%BSortMethod)) .eq. "State-Specific" ) ) then
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Selected a NON-State-Specific Initial Target')")
      if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomAtom]: Sampling the Initial internal Level of the Target')")

      RandNum = RanddwsVec(iPES)
      do i = 1,Species(iSpec)%ListStates%NStates
        if ( Species(iSpec)%ListStates%States(i)%rlim >= RandNum ) exit
      end do
      Traj_iState   =   i

    else     
      if (i_Debug_Loc) call Logger%Write( "Either Fixed State or State that Maximizes Boltzmann Distribution" )

      Traj_iState   =   1

    end if

    Traj_jqn        =   Species(iSpec)%ListStates%States(Traj_iState)%jqn
    Traj_vqn        =   Species(iSpec)%ListStates%States(Traj_iState)%vqn
    if (i_Debug_Loc) call Logger%Write( "Indx of the Selected Level: ", Traj_iState )
    if (i_Debug_Loc) call Logger%Write( "jqn  of the Selected Level: ", Traj_jqn    )
    if (i_Debug_Loc) call Logger%Write( "vqn  of the Selected Level: ", Traj_vqn    )

    if (i_Debug_Loc) call Logger%Write( "Species(iSpec)%DiatPot%xmui2 = ", Species(iSpec)%DiatPot%xmui2, "; Traj_jqn = ", Traj_jqn )
    Vc_R2           =   Species(iSpec)%DiatPot%xmui2 * ( Traj_jqn + Half )**2
  
    associate( S => Species(iSpec)%ListStates%States(Traj_iState) )
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


    ! ==============================================================================================================
    !   COMPUTING THE KINETIC ENERGY AND UPDATING SOME COMPONENT IN THE STATE OBJECT
    ! ==============================================================================================================
    if (i_Debug_Loc) write(Logger%Unit,"(10x,'[SetInitialState_DiatomDiatom]: Calling ComputeKineticEnergy')")
    call ComputeKineticEnergy( Ekin, State, Vc_R2, Species(iSpec)%DiatPot, iPES, i_Debug=i_Debug_Loc )
    if (i_Debug_Loc) then
      if (i_Debug_Loc) call Logger%Write( "-> Kinetic energy:                   Ekin      = ", Ekin      )
      if (i_Debug_Loc) call Logger%Write( "-> Estimate of outter turning point: State%ro  = ", State%ro  )
      if (i_Debug_Loc) call Logger%Write( "-> Period:                           State%Tau = ", State%Tau )
    end if
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   COMPUTING THE INITIAL BOND LENGTH AND ITS DERIVATIVES
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Computing the bond length and its derivatives. Calling BondLength" )
    Ekin = Zero                                                                                                                     ! Setting the kinetic energy at the potential minimum
    call BondLength( BL, dBL, State%ro, State%tau, Ekin, Vc_R2, Species(iSpec)%DiatPot, iPES, i_Debug=i_Debug_Loc, i_Debug_Deep=i_Debug_Deep )! Computng BL/bBL output, all the rest is input
    if (i_Debug_Loc) then
      if (i_Debug_Loc) call Logger%Write( "-> Bond length:        BL       = ", BL  )
      if (i_Debug_Loc) call Logger%Write( "-> Bond length Deriv.: dBL      = ", dBL )
    end if
    ! ==============================================================================================================


    ! ==============================================================================================================
    !   COMPUTING THE MOLECULE COORDINATES AND VELOCITIES
    ! ==============================================================================================================
    ! choose random variables for diatomic orientation; set some variables defining the angular momentum vector Omega
    ! ==============================================================================================================
    if (i_Debug_Loc) call Logger%Write( "Computing the diatomic coordinates" )

    ! Shouldn't it be:    Omega = sqrt(Two * Vc_R2 * xmui2) ??? (SIMONE)
    Omega = Two * sqrt( Vc_R2 * Species(iSpec)%DiatPot%xmui2 ) / BL**2                                                                                         ! Setting the angular velocity for rotation
    if (i_Debug_Loc) call Logger%Write( "-> Angular velocity: Omega = ", Omega )
    if (i_Debug_Loc) call Logger%Write( "Calling Compute_Coordinates_Velocities" )
    
    call Compute_Coordinates_Velocities( Species(iSpec)%GetAtomsMass(), BL, dBL, Omega, Qijk(:,:,iMol), dQijk(:,:,iMol), iTraj, iPES, NTrajOverall, iProc, i_Debug=i_Debug_Loc )! Computing the coordinates of the target
    if (i_Debug_Loc) then
       call Logger%Write( "-> Positions of atoms in target:" )
       call Logger%Write( "  -> Atom A: ", Qijk(:,1,iMol) )
       call Logger%Write( "  -> Atom B: ", Qijk(:,2,iMol) )
       call Logger%Write( "-> Velocities of atoms in target:" )
       call Logger%Write( "  -> Atom A:", dQijk(:,1,iMol) )
       call Logger%Write( "  -> Atom B:", dQijk(:,2,iMol) )
    end if
! ==============================================================================================================

  enddo

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


End Module