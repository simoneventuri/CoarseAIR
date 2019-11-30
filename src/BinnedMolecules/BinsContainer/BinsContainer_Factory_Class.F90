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

Module BinsContainer_Factory_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public  ::    BinsContainer_Factory_Type

  Type      ::    BinsContainer_Factory_Type
  contains
    private
    procedure ,nopass ,public ::  Construct => Construct_BinsContainer
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Construct_BinsContainer( Input, LevelsContainer, iMol, FldrPath, BinsContainer, WriteFlg, i_Debug )

  use Input_Class                    ,only:    Input_Type
  use LevelsContainer_Class          ,only:    LevelsContainer_Type
  use BinsContainer_Class            ,only:    BinsContainer_Type
  use StS_BinsContainer_Class        ,only:    StS_BinsContainer_Type
  use VibSpec_BinsContainer_Class    ,only:    VibSpec_BinsContainer_Type
  use Energy_BinsContainer_Class     ,only:    Energy_BinsContainer_Type
  use FromFile_BinsContainer_Class   ,only:    FromFile_BinsContainer_Type
  use Hybrid_BinsContainer_Class     ,only:    Hybrid_BinsContainer_Type

  Type(Input_Type)                            ,intent(inout)  ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  character(*)                                ,intent(in)     ::    FldrPath
  class(BinsContainer_Type)    ,allocatable   ,intent(out)    ::    BinsContainer
  logical                           ,optional ,intent(in)     ::    WriteFlg
  logical                           ,optional ,intent(in)     ::    i_Debug

  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Construct_BinsContainer")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Write( "Sorting Method for Molecule Nb ", iMol, ": Input%BSortMethod(iMol) = ", Input%BSortMethod(iMol) )
  
  select case (Input%BSortMethod(iMol))    
    case('State-Specific')
      Input%NBins(iMol) = LevelsContainer%NStates
      if (i_Debug) call Logger%Write( "Sorting Method for Molecule Nb", iMol, " is State-Specific. Setting NBins = ", Input%NBins )
      if (i_Debug_Loc) call Logger%Write( "Treating the Levels State-Specific" )
      allocate( StS_BinsContainer_Type :: BinsContainer )

    case('Vib-Specific')
      Input%NBins(iMol) = LevelsContainer%maxvqn + 1
      if (i_Debug) call Logger%Write( "Sorting Method for Molecule Nb", iMol, " is Vib-Specific. Setting NBins = ", Input%NBins )
      if (i_Debug_Loc) call Logger%Write( "Sorting the Levels Vibrationally Specific" )
      allocate( VibSpec_BinsContainer_Type :: BinsContainer )

    case('RoVib-CG')
      if (i_Debug_Loc) call Logger%Write( "Sorting the Levels by increasing Ro-Vibrational Energy" )
      allocate( Energy_BinsContainer_Type :: BinsContainer )

    case('From-File')
      if (i_Debug_Loc) call Logger%Write( "Reading the Mapping Level-To-Bin from File" )
      allocate( FromFile_BinsContainer_Type :: BinsContainer )

    case('Hybrid')
      if (i_Debug_Loc) call Logger%Write( "Sorting the Levels Based on a Hybrid Strategy" )
      allocate( Hybrid_BinsContainer_Type :: BinsContainer )
      
    case default
      allocate( StS_BinsContainer_Type :: BinsContainer )
      if (i_Debug_Loc) call Logger%Write( "Input%BSortMethod not specified for Molecule Nb ", iMol, "; Levels are going to be sorted StS" )

  end select

  write(Input%NBins_char(iMol), '(I6)') Input%NBins(iMol)
  if (i_Debug) call Logger%Write( "Input%NBins_char(iMol) for Molecule Nb", iMol, " = ", Input%NBins_char(iMol) )


  allocate( BinsContainer%PathToBinMolFldr, source = trim(adjustl( trim(adjustl(FldrPath)) // adjustl(trim( Input%Molecules_Name(iMol) )) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/' )) )
  if (i_Debug_Loc)  call Logger%Write( "Creating Folder ", BinsContainer%PathToBinMolFldr )
  call system('mkdir -p ' // BinsContainer%PathToBinMolFldr )


  if (i_Debug_Loc) call Logger%Write( "Calling BinsContainer%Initialize" )
  call BinsContainer%Initialize( Input, LevelsContainer, iMol, WriteFlg, i_Debug=i_Debug_Loc )

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module