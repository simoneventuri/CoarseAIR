! -*-F90-*
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Berkan Bolkan and Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
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
Module StS_BinsContainer_Class

#include "../../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, One, Two
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    StS_BinsContainer_Type


  Type    ,extends(BinsContainer_Type)          ::    StS_BinsContainer_Type
    
  contains
    procedure          ::  Initialize     =>    Initialize_StS_BinsContainer    
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.                                                    
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                               DEFERRED PROCEDURES for StS BinsContainer
! **************************************************************************************************************
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Subroutine Initialize_StS_BinsContainer( This, Input, LevelsContainer, iMol, WriteFlg, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(StS_BinsContainer_Type)               ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    WriteFlg
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  integer                                                     ::    iLevels
  integer                                                     ::    iBins
  integer                                                     ::    iBins_temp
  integer                                                     ::    ijqn
  integer                                                     ::    iwrite
  character(:)                                   ,allocatable ::    FileName
  character(:)                                   ,allocatable ::    FileName2
  character(5)                                                ::    iBins_char
  integer                                                     ::    Status
  integer                                                     ::    Unit

  logical                                                     ::    i_Debug_Loc
  logical                                                     ::    i_Debug_Deepest
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_StS_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()


  This%NBins = Input%NBins(iMol)


  ! Allocate Bins
  if (.not. allocated(This%Bin)) then 
    allocate(This%Bin(Input%NBins(iMol)), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%Bin" )
    if (i_Debug_Loc)  call Logger%Write( "Allocated ", This%NBins, " Bins for the Molecule Nb", iMol )
  end if


  ! Assigning the correspondence level -> bin through the vector Lvl_to_Bin
  if (i_Debug_Loc)  call Logger%Write( "Assigning the correspondence level -> bin through the vector This%Lvl_to_Bin")
  if (.not. allocated(This%Lvl_to_Bin)) then
    allocate(This%Lvl_to_Bin(LevelsContainer%NStates), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%Lvl_to_Bin" )
    This%Lvl_to_Bin=0
  end if
  do iLevels = 1,LevelsContainer%NStates
    This%Lvl_to_Bin(iLevels) = iLevels
  end do
  
  LevelsContainer%States%einteV_scaled = LevelsContainer%States%einteV - LevelsContainer%MineinteV
  
  if (WriteFlg) then
  
    FileName = This%PathToBinMolFldr // '/qns_Bin.dat'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName ) 
    write(Unit,*) "#      vqn       jqn     En [eV]"
  
  end if

  ! Writing all the bin vqns and jqns in the related Bin_*.dat file
  This%Bin%NLevels = 1
  do iBins = 1,This%NBins
    write(iBins_char,'(I5)') iBins
    
    if (WriteFlg) then
    
      FileName2 = This%PathToBinMolFldr // '/levels_Bin' // trim(adjustl(iBins_char)) // '.inp'
      open( File=FileName2, Unit=101, status='REPLACE', iostat=Status )
        if (Status/=0) call Error( "Error opening file: " // FileName2 ) 
        write(101,'(A150)') ('######################################################################################################################################################')
        write(101,'(A150)') ('# vqn   : the vibrational q.n. of the i''th quantum state                                                                                             ')
        write(101,'(A150)') ('# jqn   : the rotational q.n. of the i''th quantum state                                                                                              ')
        write(101,'(A150)') ('# eint  : internal energy of i''th quantum state                                                                                                      ')
        write(101,'(A150)') ('# egam  : Half width of i''th quantum state                                                                                                           ')
        write(101,'(A150)') ('# rmin  : the position of the potential minimum (included centrifugal potential) for i''th quantum state                                              ')
        write(101,'(A150)') ('# vmin  : the value of the potential minimun (inc. cent. pot.)                                                                                        ')
        write(101,'(A150)') ('# vmax  : the value of the local potential maximum (inc. cent. pot.)                                                                                  ')
        write(101,'(A150)') ('# tau   : the vibrational period of the i''th quantum state                                                                                           ')
        write(101,'(A150)') ('# ri    : inner turning point                                                                                                                         ')
        write(101,'(A150)') ('# ro    : outter turning point                                                                                                                        ')
        write(101,'(A150)') ('# rmax  : location of maximum in centrifugal barrier                                                                                                  ')
        write(101,'(A150)') ('######################################################################################################################################################')
        write(101,'(A150)') ('#    vqn  jqn   eint            egam            rmin          rmax          vmin            vmax            tau           ri             ro           ')
        write(101,'(A150)') ('######################################################################################################################################################')
        write(101,2) int(LevelsContainer%States(iBins)%vqn),  int(LevelsContainer%States(iBins)%jqn),  LevelsContainer%States(iBins)%eint, &
                        LevelsContainer%States(iBins)%egam,      LevelsContainer%States(iBins)%rmin,  LevelsContainer%States(iBins)%rmax, &
                        LevelsContainer%States(iBins)%vmin,      LevelsContainer%States(iBins)%vmax,  LevelsContainer%States(iBins)%tau,  &
                        LevelsContainer%States(iBins)%ri,        LevelsContainer%States(iBins)%ro
      
      close(101) 
      
    end if

    This%Bin(iBins)%vqnFirst     = LevelsContainer%States(iBins)%vqn
    This%Bin(iBins)%jqnFirst     = LevelsContainer%States(iBins)%jqn

    if (WriteFlg) then
      write(Unit,1) LevelsContainer%States(iBins)%vqn, LevelsContainer%States(iBins)%jqn, LevelsContainer%States(iBins)%einteV

      This%qns_to_Bin( LevelsContainer%States(iBins)%vqn, LevelsContainer%States(iBins)%jqn ) = iBins
    end if
    
  end do

  if (WriteFlg) then
    close(Unit)

    FileName = This%PathToBinMolFldr // '/../NLevels.inp'
    open( File=FileName, Unit=101, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName2 ) 
      write(101,'(I10)') LevelsContainer%NStates
    close(101) 
  
    call This%WriteBinsData( Input, LevelsContainer, iMol )
  
  end if

  
  1 format (3X  2I10, d20.10)
  2 format (1X, 2I5, 9es15.7)
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module