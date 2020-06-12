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

Module DiatomicPotential_Factory_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  private
  public  ::    DiatomicPotential_Factory_Type

  Type                          ::    DiatomicPotential_Factory_Type
  contains
    private
    procedure ,nopass ,public   ::    Construct =>  Construct_DiatomicPotential
  End Type

  logical   ,parameter          ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Construct_DiatomicPotential( Atoms, iA, Input, DiatPot, i_Debug )

  use Atom_Class                             ,only:  Atom_Type
  use Input_Class                            ,only:  Input_Type
  use DiatomicPotential_Class                ,only:  DiatomicPotential_Type
  use Null_DiatomicPotential_Class           ,only:  Null_DiatomicPotential_Type
  use Morse_DiatomicPotential_Class          ,only:  Morse_DiatomicPotential_Type
  use Modified_Morse_DiatomicPotential_Class ,only:  Modified_Morse_DiatomicPotential_Type
  use N2_LeRoy_DiatomicPotential_Class       ,only:  N2_LeRoy_DiatomicPotential_Type
  use O2_UMN_DiatomicPotential_Class         ,only:  O2_UMN_DiatomicPotential_Type
  use O2_Basel_DiatomicPotential_Class       ,only:  O2_Basel_DiatomicPotential_Type
  use O2_Varandas_DiatomicPotential_Class    ,only:  O2_Varandas_DiatomicPotential_Type
  use O2_NASA_DiatomicPotential_Class        ,only:  O2_NASA_DiatomicPotential_Type
  use CO_DiatomicPotential_Class             ,only:  CO_DiatomicPotential_Type
  use N2_UMN_ForN4_DiatomicPotential_Class   ,only:  N2_UMN_ForN4_DiatomicPotential_Type
  use N2_UMN_ForN2O2_DiatomicPotential_Class ,only:  N2_UMN_ForN2O2_DiatomicPotential_Type
  use NO_UMN_DiatomicPotential_Class         ,only:  NO_UMN_DiatomicPotential_Type
  use NO_Basel_DiatomicPotential_Class       ,only:  NO_Basel_DiatomicPotential_Type
  use CN_UIUC_DiatomicPotential_Class        ,only:  CN_UIUC_DiatomicPotential_Type
  use HN_UIUC_DiatomicPotential_Class        ,only:  HN_UIUC_DiatomicPotential_Type
  use CH_UIUC_DiatomicPotential_Class        ,only:  CH_UIUC_DiatomicPotential_Type

  type(Atom_Type) ,dimension(:)                 ,intent(in)     ::    Atoms
  integer         ,dimension(:)                 ,intent(in)     ::    iA
  type(Input_Type)                              ,intent(in)     ::    Input
  class(DiatomicPotential_Type) ,allocatable    ,intent(out)    ::    DiatPot
  logical                         ,optional     ,intent(in)     ::    i_Debug

  logical                                                       ::    i_Debug_Loc
  integer                                                       ::    iMol, jMol
  character(:)  ,allocatable                                    ::    SpeciesName
  real(rkp)                                                     ::    Mass1       
  real(rkp)                                                     ::    Mass2       
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Construct_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  select case ( size(iA) )

    case (2:)
      if (i_Debug_Loc) call Logger%Write( "Constructing diatomic potential for atoms: " ,Atoms(iA(1))%Name, Atoms(iA(2))%Name )
      allocate( SpeciesName , source = trim( GetSpeciesName(Atoms(iA(1)), Atoms(iA(2)) ) ) )
      Mass1 = Atoms(iA(1))%Mass
      Mass2 = Atoms(iA(2))%Mass

      iMol = 1
      do jMol = 1,Input%NMolecules
        if ( trim(SpeciesName) == trim(Input%Molecules_Name(jMol)) ) then
          iMol = jMol
          exit
        end if
      end do

      if (i_Debug_Loc) call Logger%Write( "Name of species: ", SpeciesName )

      if (adjustl(trim(Input%Diatomic_Model(iMol))) == 'Morse') then
        if (i_Debug_Loc) call Logger%Write( "Constructing a Morse_DiatomicPotential_Type object" )
        allocate( Morse_DiatomicPotential_Type :: DiatPot )

      elseif (adjustl(trim(Input%Diatomic_Model(iMol))) == 'Modified_Morse') then
        if (i_Debug_Loc) call Logger%Write( "Constructing a Modified_Morse_DiatomicPotential_Type object" )
        allocate( Modified_Morse_DiatomicPotential_Type :: DiatPot )
      
      else

        select case ( SpeciesName )

          case ('N2', 'NNa', 'NNb', 'NNc', 'NNd', 'NaNb','NcNd','NbNc','NaNd','NaNc','NbNd' )
            if (Input%Diatomic_Model(iMol) == 'LeRoy') then
              allocate( N2_LeRoy_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'UMN_ForN4') then
              allocate( N2_UMN_ForN4_DiatomicPotential_Type :: DiatPot )   
            elseif (Input%Diatomic_Model(iMol) == 'UMN_ForN2O2') then
              allocate( N2_UMN_ForN2O2_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if
            
          case ('CO')
            if (Input%Diatomic_Model(iMol) == 'NASA') then
              allocate( CO_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if

          case ('CH')
            if (Input%Diatomic_Model(iMol) == 'UIUC') then
              allocate( CH_UIUC_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if

          case ('CN')
            if (Input%Diatomic_Model(iMol) == 'UIUC') then
              allocate( CN_UIUC_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if

          case ('HN')
            if (Input%Diatomic_Model(iMol) == 'UIUC') then
              allocate( HN_UIUC_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if
            
          case ('O2', 'OOa', 'OOb', 'OOc', 'OOd', 'OaOb','OaOc','OaOd','ObOc','ObOd','OcOd','O16O16','O18O18','O16O18','O18O16')
            if (Input%Diatomic_Model(iMol) == 'NASA') then
              allocate( O2_NASA_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'UMN') then
              allocate( O2_UMN_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'Varandas') then
              allocate( O2_Varandas_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'Basel') then
              allocate( O2_Basel_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if
            
          case ('ArC');  
            if (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if
          
          case ('ArO');  
            if (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if
            
          case ('NO', 'ON')
            if (Input%Diatomic_Model(iMol) == 'UMN') then
              allocate( NO_UMN_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'Basel') then
              allocate( NO_Basel_DiatomicPotential_Type :: DiatPot )
            elseif (Input%Diatomic_Model(iMol) == 'NONE') then
              allocate( Null_DiatomicPotential_Type :: DiatPot )
            else
              call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
            end if
            
          case ('C2');  allocate( Null_DiatomicPotential_Type :: DiatPot )
            call Error( "Diatomic Potential Model for " // SpeciesName // " Molecule not implemented yet." )
          case ('H2');  allocate( Null_DiatomicPotential_Type :: DiatPot )
            call Error( "Diatomic Potential Model for " // SpeciesName // " Molecule not implemented yet." )
          case ('OH');  allocate( Null_DiatomicPotential_Type :: DiatPot )
            call Error( "Diatomic Potential Model for " // SpeciesName // " Molecule not implemented yet." )
        
          case ('NULL');  allocate( Null_DiatomicPotential_Type :: DiatPot )
            call Error( "Unrecognized Molcule." )
          
          case default; call Error( "Diatomic Potential Model not supported: Species Name = " // SpeciesName // '; Input%Diatomic_Model(iMol) = ' // Input%Diatomic_Model(iMol) )
          
        end select

      end if


    case (1)   
      if (i_Debug_Loc) call Logger%Write( "Constructing diatomic potential for atom: ", Atoms(1)%Name )
      allocate( SpeciesName , source = 'NULL' )
      Mass1 = Atoms(1)%Mass
      Mass2 = Zero

      allocate( Null_DiatomicPotential_Type :: DiatPot )


    case (:0)
      call Error( "Error in Construct_DiatomicPotential: size(Atoms) = 0" )

  end select


  if (i_Debug_Loc) call Logger%Write( "Calling DiatPot%Initialize" )
  call DiatPot%Initialize( Input, SpeciesName, iMol, Mass1, Mass2, i_Debug=i_Debug_Loc )

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function GetSpeciesName( Atom1, Atom2 ) result(SpeciesName)
  use Atom_Class                  ,only:  Atom_Type
  type(Atom_Type)                               ,intent(in)     ::    Atom1
  type(Atom_Type)                               ,intent(in)     ::    Atom2
  character(6)                                                  ::    SpeciesName
  character(6)                                                  ::    SpeciesName1
  character(6)                                                  ::    SpeciesName2
  character(6)                                                  ::    Name
  type(Atom_Type) ,dimension(2)                                 ::    Atoms

!   write(Logger%Unit,"(10x,'[GetSpeciesName]: Atoms(1)%Name = ',g0)") Atoms(1)%Name
!   write(Logger%Unit,"(10x,'[GetSpeciesName]: Atoms(2)%Name = ',g0)") Atoms(2)%Name

  Atoms(1)  = Atom1
  Atoms(2)  = Atom2

  SpeciesName     =     ''

  if ( (len(Atoms(1)%Name) == 1) .and. ( trim(Atoms(1)%Name) == trim(Atoms(2)%Name) ) ) then
    SpeciesName1 = trim(Atoms(1)%Name) // '2'
    SpeciesName2 = SpeciesName1
  else
    SpeciesName1 = trim(Atoms(1)%Name) // trim(Atoms(2)%Name)
    SpeciesName2 = trim(Atoms(2)%Name) // trim(Atoms(1)%Name)
    SpeciesName  = SpeciesName1
  end if

!   write(Logger%Unit,"(10x,'[GetSpeciesName]: SpeciesName1 = ',g0)") SpeciesName1
!   write(Logger%Unit,"(10x,'[GetSpeciesName]: SpeciesName2 = ',g0)") SpeciesName2

  Name = 'N2'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'O2'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'C2'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'H2'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'NO'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'CO'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'CN'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'HN'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'HO'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)
  Name = 'CH'; if ( (SpeciesName1==Name) .or. (SpeciesName2==Name) ) SpeciesName = trim(Name)

!   write(Logger%Unit,"(10x,'[GetSpeciesName]: SpeciesName = ',g0)") SpeciesName

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
