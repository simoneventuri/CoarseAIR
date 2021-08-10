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

Module Degeneracy_Factory_Class

  use Parameters_Module     ,only:  rkp
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public    ::    Degeneracy_Factory_Type

  Type      ::    Degeneracy_Factory_Type
  contains
    private
    procedure ,nopass ,public ::  Define_Degeneracy
  End Type

  logical   ,parameter    ::    i_Debug_Global = .False.

  contains

Subroutine Define_Degeneracy( Input, Degeneracy, iMol, i_Debug )

  use Input_Class          ,only:    Input_Type
  use Degeneracy_Class     ,only:    Degeneracy_Type
  use N2_Degeneracy_Class  ,only:    N2_Degeneracy_Type                                                                      
  use NO_Degeneracy_Class  ,only:    NO_Degeneracy_Type                                                                      
  use CO_Degeneracy_Class  ,only:    CO_Degeneracy_Type
  use O2_Degeneracy_Class  ,only:    O2_Degeneracy_Type
  use CH_Degeneracy_Class  ,only:    CH_Degeneracy_Type
  use CN_Degeneracy_Class  ,only:    CN_Degeneracy_Type
  use HN_Degeneracy_Class  ,only:    HN_Degeneracy_Type
  use G2_Degeneracy_Class  ,only:    G2_Degeneracy_Type
  use OH_Degeneracy_Class  ,only:    OH_Degeneracy_Type

  type(Input_Type)                                  ,intent(in)     ::    Input
  class(Degeneracy_Type)               ,allocatable ,intent(out)    ::    Degeneracy
  integer                                           ,intent(in)     ::    iMol
  logical                                 ,optional ,intent(in)     ::    i_Debug

  logical                                                           ::    i_Debug_Loc
  integer                                                           ::    Status

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Define_Degeneracy")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Write( "Molecule for the Degeneracy to be Defined: Input%Molecules_Name(i) = ", Input%Molecules_Name(iMol) )
  select case ( trim(adjustl(Input%Molecules_Name(iMol))) )    
    case('N2', 'NN', 'NaNb','NcNd','NbNc','NaNd','NaNc','NbNd','NbNa','NdNc','NcNb','NdNa','NcNa','NdNb')
      if (i_Debug_Loc) call Logger%Write( "Defining a N2_Degeneracy_Type object" )
      allocate( N2_Degeneracy_Type :: Degeneracy )
     case('NO', 'ON')
      if (i_Debug_Loc) call Logger%Write( "Defining a NO_Degeneracy_Type object" )
      allocate( NO_Degeneracy_Type :: Degeneracy )  
    case('CO', 'OC')
      if (i_Debug_Loc) call Logger%Write( "Defining a CO_Degeneracy_Type object" )
      allocate( CO_Degeneracy_Type :: Degeneracy )
    case('O2', 'OO', 'OaOb','OcOd','ObOc','OaOd','OaOc','ObOd','ObOa','OdOc','OcOb','OdOa','OcOa','OdOb','O16O16','O18O18','O16O18','O18O16')
      if (i_Debug_Loc) call Logger%Write( "Defining a O2_Degeneracy_Type object" )
      allocate( O2_Degeneracy_Type :: Degeneracy )
    case('CN', 'NC')
      if (i_Debug_Loc) call Logger%Write( "Defining a CN_Degeneracy_Type object" )
      allocate( CN_Degeneracy_Type :: Degeneracy )
    case('CH', 'HC')
      if (i_Debug_Loc) call Logger%Write( "Defining a CH_Degeneracy_Type object" )
      allocate( CH_Degeneracy_Type :: Degeneracy )
    case('HN', 'NH')
      if (i_Debug_Loc) call Logger%Write( "Defining a HN_Degeneracy_Type object" )
      allocate( HN_Degeneracy_Type :: Degeneracy )
    case('GG')
      if (i_Debug_Loc) call Logger%Write( "Defining a G2_Degeneracy_Type object" )
      allocate( G2_Degeneracy_Type :: Degeneracy )
    case('OH', 'HO')
      if (i_Debug_Loc) call Logger%Write( "Defining a OH_Degeneracy_Type object" )
      allocate( OH_Degeneracy_Type :: Degeneracy )
    case default
      call Error( "Molecule not supported yet for Degeneracy Definition: Input%Molecules_Name(i) = " // Input%Molecules_Name(iMol) )
  end select

  if (i_Debug_Loc) call Logger%Write( "Calling Degeneracy%Initialize" )
  call Degeneracy%Initialize(Input)

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine

End Module
