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
Module FromFile_BinsContainer_Class

#include "../../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, One, Two
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    FromFile_BinsContainer_Type


  Type    ,extends(BinsContainer_Type)          ::    FromFile_BinsContainer_Type
    
  contains
    procedure          ::  Initialize     =>    Initialize_FromFile_BinsContainer    
  End Type

  logical                         ,parameter    ::    i_Debug_Global = .False.                                                    
  contains
  

! **************************************************************************************************************
! **************************************************************************************************************
!                               DEFERRED PROCEDURES for FromFile BinsContainer
! **************************************************************************************************************
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!
Subroutine Initialize_FromFile_BinsContainer( This, Input, LevelsContainer, iMol, WriteFlg, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(FromFile_BinsContainer_Type)          ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    WriteFlg
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  integer                                                     ::    iLevels, jLevels
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
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_FromFile_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()

  
  ! Allocate Bins
  if (.not. allocated(This%Bin)) then 
    allocate(This%Bin(Input%NBins(iMol)), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%Bin" )
    if (i_Debug_Loc)  call Logger%Write( "Allocated ", Input%NBins(iMol), " Bins for the Molecule Nb", iMol )
  end if


  ! Assigning the correspondence level -> bin through the vector Lvl_to_Bin
  if (i_Debug_Loc)  call Logger%Write( "Assigning the correspondence level -> bin through the vector This%Lvl_to_Bin")
  if (.not. allocated(This%Lvl_to_Bin)) then 
    allocate(This%Lvl_to_Bin(LevelsContainer%NStates), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%Lvl_to_Bin" )
    This%Lvl_to_Bin=0
  end if
  
  
  ! Reading from "/E_low.dat" the Bins Energy Minima
  FileName = trim(adjustl(Input%DtbPath)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
             trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/' // trim(adjustl(Input%LevToBinFile(iMol)))
  if (i_Debug_Loc)  call Logger%Write( "Reading the File for the Levels-Bins Mapping" )
  if (i_Debug_Loc)  call Logger%Write( "-> Opening file: ", FileName )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )
  read(Unit,*)
  do iLevels = 1,LevelsContainer%NStates
    read(Unit,*) jLevels, This%Lvl_to_Bin(jLevels)
  end do
  close(Unit)


  ! Finding the first (v,J)s for each of the bins (vectors vqnFirst and jqnFirst)
  ! Assigning the correspondence (v,J)_level -> bin through the matrix qns_to_Bin
  ! Counting the number of levels per bin (vector NLevels)
  ! Writing all the bin vqns and jqns in the related Bin_*.dat file
  This%Bin%NLevels = 0
  if (i_Debug_Loc)  call Logger%Write( "Finding the first (v,J)s for each of the bins (vectors vqnFirst and jqnFirst) ... ")
  do iBins = 1,Input%NBins(iMol)
    write(iBins_char,'(I5)') iBins
    
    FileName = This%PathToBinMolFldr // '/qns_Bin' // trim(adjustl(iBins_char)) // '.dat'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName ) 
    write(Unit,*) "          vqn         jqn"
    
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
    
    iwrite = 1
    do iLevels = 1,LevelsContainer%NStates 
         
      if (This%Lvl_to_Bin(iLevels) == iBins) then
         
        if (iwrite == 1) then
          if (LevelsContainer%States(iLevels)%einteV < This%Bin(iBins)%MinLevEintev) then
            This%Bin(iBins)%MinLevEintev = LevelsContainer%States(iLevels)%einteV
          end if
          This%Bin(iBins)%vqnFirst     = LevelsContainer%States(iLevels)%vqn
          This%Bin(iBins)%jqnFirst     = LevelsContainer%States(iLevels)%jqn
          iwrite = 2
        end if
         
       This%qns_to_Bin( LevelsContainer%States(iLevels)%vqn, LevelsContainer%States(iLevels)%jqn ) = iBins
       This%Bin(iBins)%NLevels = This%Bin(iBins)%NLevels + 1
       
       write(Unit,1) LevelsContainer%States(iLevels)%vqn, LevelsContainer%States(iLevels)%jqn
       
       write(101,2) int(LevelsContainer%States(iLevels)%vqn),  int(LevelsContainer%States(iLevels)%jqn),  LevelsContainer%States(iLevels)%eint, &
                        LevelsContainer%States(iLevels)%egam,      LevelsContainer%States(iLevels)%rmin,  LevelsContainer%States(iLevels)%rmax, &
                        LevelsContainer%States(iLevels)%vmin,      LevelsContainer%States(iLevels)%vmax,  LevelsContainer%States(iLevels)%tau,  &
                        LevelsContainer%States(iLevels)%ri,        LevelsContainer%States(iLevels)%ro
      
      end if
       
    end do
   
    close(Unit) 
    
    close(101) 
   
  end do
  if (i_Debug_Loc)  then
    call Logger%Write( "        Bin    Min En [eV]    Nb Levels  Levels %")
    do iBins = 1,Input%NBins(iMol)
      call Logger%Write( "         ", iBins, This%Bin(iBins)%MinLevEinteV,  This%Bin(iBins)%NLevels, real(This%Bin(iBins)%NLevels) * 1.e2_rkp / real(LevelsContainer%NStates) )
    end do
  end if 
  
  
  call This%WriteBinsData( Input, LevelsContainer, iMol )
  
  1 format (3X  2I10)
  2 format (1X, 2I5, 9es15.7)
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module