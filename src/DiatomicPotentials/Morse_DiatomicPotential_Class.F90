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

Module Morse_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Half, One, Two, Kcm_To_Hartree, B_To_Ang
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type
  use Error_Class               ,only:  Error

  implicit none

  private
  public  ::    Morse_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)    :: Morse_DiatomicPotential_Type
    real(rkp)                               :: re
    real(rkp)                               :: De
    real(rkp)                               :: Beta
  contains
    procedure         ::    Initialize        =>    Initialize_Morse_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_Morse
    procedure         ::    DiatomicPotential =>    DiatomicPotential_Morse
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiatPot = 'Morse'
  
  real(rkp)                               :: re
  real(rkp)                               :: De
  real(rkp)                               :: Beta

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Morse_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(Morse_DiatomicPotential_Type)    ,intent(out)       ::    This
  type(Input_Type)                       ,intent(in)        ::    Input
  character(:) ,allocatable              ,intent(in)        ::    SpeciesName
  integer                                ,intent(in)        ::    iMol
  real(rkp)                              ,intent(in)        ::    Mass1
  real(rkp)                              ,intent(in)        ::    Mass2
  logical ,optional                      ,intent(in)        ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  character(:)                    ,allocatable              ::    Morse_file
  integer                                                   ::    Status
  integer                                                   ::    Unit, idum
  
  integer                                                   ::    iostat_SPES_input
  character(:)                    ,allocatable              ::    i_case
  character(150)                                            ::    line_input
  integer                                                   ::    i_char
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_Morse_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = trim(Name_DiatPot) )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =    iMol
  This%Initialized  =   .True.

  This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui
  
  ! ==============================================================================================================
  !     READING MORSE INPUT FILE
  ! ==============================================================================================================
  if (trim(adjustl(Input%DiatPot_ParamsFile(iMol))) == 'NONE') then
    Morse_file = trim(adjustl(Input%DtbPath))  // '/Molecules/' // trim(adjustl(This%SpeciesName)) // '/Morse/Morse.dat'
  elseif (trim(adjustl(Input%DiatPot_ParamsFile(iMol))) == 'Local') then
    Morse_file = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/LEPS.dat'
  else
    Morse_file = trim(adjustl(Input%DtbPath))  // '/Molecules/' // trim(adjustl(This%SpeciesName)) // '/Morse/' // trim(adjustl(Input%DiatPot_ParamsFile(iMol)))
  end if
  if (i_Debug_Loc) call Logger%Write( "Reading the Morse Parameters file" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", Morse_file )
  open( File=Morse_file, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // Morse_file )
  
  do 
    
    read(Unit,'(A)',iostat=Status) line_input
         
    if (Status /= 0) then
      exit
    else
    
      if ( (line_input(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
      
        i_char = 1
        do
          if (line_input(i_char:i_char) == '=') exit
          i_char = i_char + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_char-2))))

        select case (adjustl((TRIM(i_case))))

          case("re")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%re
            call Logger%Write( "re:      this%re  = ", this%re )
         
          case("De")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%De
            call Logger%Write( "De:      this%De  = ", this%De )

          case("Beta")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Beta 
            call Logger%Write( "Beta:   this%Beta  = ", this%Beta )
            
        end select
      
      end if
   
    end if
   
  end do

  close(Unit)

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_Morse( This, R, V, dV )

  use Parameters_Module   ,only:  Zero, One

  class(Morse_DiatomicPotential_Type)       ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Distances between nuclear centers [bohr]
  real(rkp)                                 ,intent(out)    ::    V                         ! Potential energy [hartree]
  real(rkp)                                 ,intent(out)    ::    dV                        ! First derivative of the potential energy wrt the distance [hartree/bohr]
  
  real(rkp)                                                 ::    RTemp
  real(rkp)                                                 ::    Xpnt
  
  RTemp = R * B_To_Ang

  Xpnt = - This%Beta * (RTemp - This%re)
  V    = This%De * ( exp(Two*Xpnt) - Two*exp(Xpnt) ) 
  dV   = This%De * Two * This%Beta * exp(Xpnt) * (One -  exp(Xpnt))

  V  = V  * Kcm_To_Hartree
  dV = dV * Kcm_To_Hartree * B_To_Ang

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_Morse( This, R ) result( V )

  class(Morse_DiatomicPotential_Type)       ,intent(in)     ::    This
  real(rkp)                                 ,intent(in)     ::    R                         ! Internuclear distance [bohr]
  real(rkp)                                                 ::    V                         ! Diatomic potential energy [hartree]
  
  real(rkp)                                                 ::    RTemp
  real(rkp)                                                 ::    Xpnt
 
  RTemp = R * B_To_Ang
 
  Xpnt = - This%Beta * (RTemp - This%re)
  V    = This%De * ( exp(Two*Xpnt) - Two*exp(Xpnt) ) 

  V = V * Kcm_To_Hartree
 
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
