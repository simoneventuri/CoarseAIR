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

Module Nb4_System_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Six
  use System_Class          ,only:  System_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Nb4_System_Type


  Type    ,extends(System_Type)                 ::    Nb4_System_Type
  contains
    procedure         ::  Initialize                 =>    Initialize_Nb4_System
    procedure         ::  AssignPairsArrangements    =>    AssignPairsArrangements_Nb4_System
    procedure         ::  AdjustArrangements         =>    AdjustArrangements_Nb4_System
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.
  
  contains
  


! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for 3Atoms System
! **************************************************************************************************************
! **************************************************************************************************************

Subroutine Initialize_Nb4_System( This, Input, i_Debug )

  use Input_Class                  ,only:  Input_Type
  
  class(Nb4_System_Type)                    ,intent(out)    ::    This
  Type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
    
  This%Name         =   Input%System
  This%Initialized  =   .True.
  
End Subroutine


Subroutine AssignPairsArrangements_Nb4_System( This, Collision, Input, i_Debug )

  use Collision_Class             ,only: Collision_Type
  use Input_Class                 ,only:  Input_Type

  class(Nb4_System_Type)                    ,intent(in)     ::    This
  Type(Collision_Type)                      ,intent(inout)  ::    Collision
  Type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iA, iP, jP, iS, iSMax
  integer                                                   ::    Status
  integer                                                   ::    pos_m, pos_e
  logical                                                   ::    i_Debug_Loc
  logical                                                   ::    found
  character(len=10)                                         ::    AtA, AtB
  character(len=10)                                         ::    MolNameA, MolNameB, MolName

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AssignPairsArrangements_Nb4_System" )
  !i_Debug_Loc   =     Logger%On()

  
  ! Collision%DiffPairs is a vector with dimension NPairs that identifies Pairs occupied by the same Molecules ( e.g.: System=CO2 -> DiffPairs=[1,1,2] )
  allocate(Collision%DiffPairs(Collision%NPairs), Stat=Status )
  if (Status/=0) call Error( "Error allocating Collision%DiffPairs" )
  if (i_Debug_Loc) call Logger%Write( "Collision%DiffPairs, Different Arrangements Vector, allocated " ) 

  if (i_Debug_Loc) call Logger%Write( "Nb of Pairs:     Collision%NPairs = ", Collision%NPairs ) 

  !-------------------------!
  ! Loop over pairs to find molecules
  iSMax = Input%NMolecules 
  do iP = 1,Collision%NPairs

     ! Form molecule name. Here one must be account for the fact the NO is the same as ON. 
     ! NOTE: the code below is valid ONLY for diatomic molecules
     ! Atom names 
     pos_e = len_trim(Collision%Pairs(iP)%Name)
     pos_m = index(Collision%Pairs(iP)%Name, '-')

     AtA = Collision%Pairs(iP)%Name(1:pos_m - 1)
     AtB = Collision%Pairs(iP)%Name(pos_m + 1:pos_e)

     ! Homonuclear molecule (e.g., N2, O2)
     if (AtA.eq.AtB) then 
        MolNameA = trim(AtA)//'2'
        MolNameB = MolNameA
     ! Heteronucler molecule (e.g., NO, CO)
     else
        MolNameA = trim(AtA)//trim(AtB)
        MolNameB = trim(AtB)//trim(AtA)
     endif
 
     !----------------------!
     ! Search for molecule
     found = .false.;  iS = 0
     search_mol : do while(.not.found)

       iS      = iS + 1
       MolName = trim(Input%Molecules_Name(iS))
       
       if ((trim(MolName).eq.trim(MolNameA)).or.(trim(MolName).eq.trim(MolNameB))) then 
          found = .true.
          Collision%DiffPairs(iP) = iS
          exit search_mol
       endif

       ! Safety check
       if (iS.eq.ISMax) call Error( "Error while filling Collision%DiffPairs array, iS>ISMax...")

     end do search_mol

  enddo

  !-------------------------!
  ! 1. All the Pre-Dissociated and Dissociated Arrangements will be identified as "2", no matter what are the Closest Atoms in the Final Arrangements
  ! 2. In obtaining Rates, do we want to distinguish for Bound/Quasi-Bound states and/or for same Molecules occupying different Pairs (i.e.: for Exchanges)?
  Collision%Arrangements = 2
  if ((trim(adjustl(Input%ConsiderQB)) .eq. "no") .or. (trim(adjustl(Input%ConsiderQB)) .eq. "NO")) then
    ! Changing Initial and Final Arrangements for avoiding Distinguishing Bound State from Quasi-Bound States
    if ((trim(adjustl(Input%ConsiderExc)) .eq. "no") .or. (trim(adjustl(Input%ConsiderExc)) .eq. "NO")) then
      ! Changing Initial and Final Arrangements for avoiding Exchange
      do iP = 1,Collision%NPairs
        do iA = 1,Collision%NAtoms - 2
           Collision%Arrangements(16*iP +4*(iA - 1):16*iP +4*(iA -1) + 1) = int(16.0_rkp*Collision%DiffPairs(iP))
        end do
      end do
    else
      do iP = 1,Collision%NPairs
         do iA = 1,Collision%NAtoms - 2
            Collision%Arrangements(16*iP + 4*(iA - 1):16*iP + 4*(iA - 1) + 1) = int(16.0_rkp*iP)
         enddo
      end do
    end if  
  else
    if ((trim(adjustl(Input%ConsiderExc)) .eq. "no") .or. (trim(adjustl(Input%ConsiderExc)) .eq. "NO")) then
      ! Changing Initial and Final Arrangements for avoiding Exchange
      do iP = 1,Collision%NPairs
         do iA = 1,Collision%NAtoms - 2
            Collision%Arrangements(16*iP + 4*(iA -1))    = int(16.0_rkp*Collision%DiffPairs(iP))
            Collision%Arrangements(16*iP + 4*(iA -1) +1) = int(16.0_rkp*Collision%DiffPairs(iP) + 1)
         enddo
      end do
    else
      do iP = 1,Collision%NPairs
         do iA = 1,Collision%NAtoms - 2
            Collision%Arrangements(16*iP + 4*(iA - 1))    = int(16.0_rkp*iP)
            Collision%Arrangements(16*iP + 4*(iA - 1) +1) = int(16.0_rkp*iP + 1)
         enddo
      end do
    end if 
  end if
  if (i_Debug_Loc) call Logger%Write( "Possible Arrangements, taking into account equal pairs:    Collision%Arrangements = ", Collision%Arrangements) 
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


Subroutine AdjustArrangements_Nb4_System( This, Collision, Input, i_Debug )

  use Collision_Class             ,only: Collision_Type
  use Input_Class                 ,only:  Input_Type

  class(Nb4_System_Type)                    ,intent(in)     ::    This
  Type(Collision_Type)                      ,intent(inout)  ::    Collision
  Type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iMol
  character(:) ,allocatable                                 ::    FileName
  integer                                                   ::    Unit
  integer                                                   ::    Status
  integer                                                   ::    AI
  integer                                                   ::    AF
  integer                                                   ::    Traj
  integer                                                   ::    iPES
  real(rkp)                                                 ::    bmax
  real(rkp)                                                 ::    b
  real(rkp)                                                 ::    jqnI
  real(rkp)                                                 ::    vqnI
  real(rkp)                                                 ::    ArrI
  real(rkp)                                                 ::    jqnF
  real(rkp)                                                 ::    vqnF
  real(rkp)                                                 ::    ArrF
  logical                                                   ::    flag
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AdjustArrangements_Nb4_System" )
  !i_Debug_Loc   =     Logger%On()

  print*,'HERE (4 atoms)!' 
  stop

  ! Reading /trajectories.out
  if (i_Debug_Loc) call Logger%Write( "Reading the output from trajectories.out and adjusting it." )   
  FileName = trim(adjustl(Input%LevelOutputDir)) // '/trajectories.out'
  open( File=FileName, Unit=10, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Opening file: " // FileName )
  read(10,*)
  
  
  ! Writing the Modified Version of trajectories.out in the Temporary File /trajectories-temp.out
  if (i_Debug_Loc) call Logger%Write( "Writing the modified output trajectories.out" )  
  FileName = trim(adjustl(Input%LevelOutputDir)) // '/trajectories-temp1.out'
  open( File=FileName, Unit=101, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )     
  write(101,*) '#    iTraj        bmax              b_i               j_i               v_i              arr_i              j_f               v_f              arr_f'                 
  
  
  do
    read(10,*,iostat=Status) Traj, iPES, bmax, b, jqnI, vqnI, ArrI, jqnF, vqnF, ArrF
    if (Status .ne. 0) exit

    
    ! Finding the new Initial Arrangements based on the Vector Collision%Arrangements
    AI    = nint(ArrI-Half)
    AI    = Collision%Arrangements(AI)

    ! Changing Initial qns for taking into account species that are not treated as molcules
    if (AI /= 2) then
      flag=.TRUE.
      do iMol = 1,Input%NMolecules
        if (Collision%Pairs(int(AI/16))%To_Molecule  == iMol) then
          flag=.FALSE.
        end if
      end do
      if (flag) then
        ! AI   = Two
        jqnI = Half
        vqnI = Half
      end if
    end if
    

    ! Finding the new Final Arrangements based on the Vector Collision%Arrangements
    AF    = nint(ArrF-Half)
    AF    = Collision%Arrangements(AF)
    
    ! Changing Final qns for taking into species that are not treated as molcules
    if (AF /= 2) then
      flag=.TRUE.
      do iMol = 1,Input%NMolecules
        if (Collision%Pairs(int(AF/16))%To_Molecule  == iMol) then
          flag=.FALSE.
        end if
      end do
      if (flag) then
        ! AF   = Two
        jqnF = Half
        vqnF = Half
      end if
    end if
    
    
    ! Changing Initial and Final Arrangements 
    ArrI = dble(AI) + Half
    ArrF = dble(AF) + Half
    
    if (AI /= 2) then
      write(101,"(i9,3x,*(es15.8,3x))") Traj, iPES, bmax, b, jqnI, vqnI, ArrI, jqnF, vqnF, ArrF
    else
      if (i_Debug_Loc) call Logger%Write( "The following trajectory has been rejected, because it started from a pre-dissociated condiion:" )   
      if (i_Debug_Loc) call Logger%Write( "Trajectory   = ", Traj, iPES, bmax, b, jqnI, vqnI, ArrI, jqnF, vqnF, ArrF )   
    end if
    
  end do
  
  
  close(10) 
  
  
  close(101)

  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine


End Module