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

Module Modified_Morse_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Half, One, Two, Three, Four, Five, Six, Seven, &
                                        Eight, Nine, Ten, B_To_Ang, eV_To_Hartree
  use Logger_Class              ,only:  Logger
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type
  use Error_Class               ,only:  Error

  implicit none

  private
  public  ::    Modified_Morse_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)    :: Modified_Morse_DiatomicPotential_Type
    real(rkp)                               :: re
    real(rkp)                               :: De
    integer                                 :: PolyOrder
    real(rkp), allocatable, dimension(:)    :: cPol
  contains
    procedure         ::    Initialize        =>    Initialize_Modified_Morse_DiatomicPotential
    procedure         ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_Modified_Morse
    procedure         ::    DiatomicPotential =>    DiatomicPotential_Modified_Morse
  End Type

  logical       ,parameter  ::    i_Debug_Global = .False.
  character(*)  ,parameter  ::    Name_DiatPot = 'Modified_Morse'
  
  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_Modified_Morse_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(Modified_Morse_DiatomicPotential_Type)    ,intent(out)       ::    This
  type(Input_Type)                       ,intent(in)        ::    Input
  character(:) ,allocatable              ,intent(in)        ::    SpeciesName
  integer                                ,intent(in)        ::    iMol
  real(rkp)                              ,intent(in)        ::    Mass1
  real(rkp)                              ,intent(in)        ::    Mass2
  logical ,optional                      ,intent(in)        ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  character(:)                    ,allocatable              ::    Modified_Morse_file
  integer                                                   ::    Status
  integer                                                   ::    Unit, idum
  
  integer                                                   ::    iostat_SPES_input
  character(:)                    ,allocatable              ::    i_case
  character(150)                                            ::    line_input
  integer                                                   ::    i_char
  integer                                                   ::    iOrd
  
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
  !     READING MODIFIED MORSE INPUT FILE
  ! ==============================================================================================================
  if (trim(adjustl(Input%DiatPot_ParamsFile(iMol))) == 'NONE') then
    Modified_Morse_file = trim(adjustl(Input%DtbPath))  // '/Molecules/' // trim(adjustl(This%SpeciesName)) // '/Modified_Morse/Modified_Morse.dat'
  elseif (trim(adjustl(Input%DiatPot_ParamsFile(iMol))) == 'Local') then
    Modified_Morse_file = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/LEPS.dat'
  else
    Modified_Morse_file = trim(adjustl(Input%DtbPath))  // '/Molecules/' // trim(adjustl(This%SpeciesName)) // '/Modified_Morse/' // trim(adjustl(Input%DiatPot_ParamsFile(iMol)))
  end if

  if (i_Debug_Loc) call Logger%Write( "Reading the Modified Morse Parameters file" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", Modified_Morse_file )
  open( File=Modified_Morse_file, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // Modified_Morse_file )
      
    read(Unit,'(A)',iostat=Status) line_input

    READ(Unit, '(d20.10)') This%re
    call Logger%Write( "re:      This%re  = ", This%re )
         
    READ(Unit, '(d20.4)') This%De
    call Logger%Write( "De:      This%De  = ", This%De )

    READ(Unit, '(i3)') This%PolyOrder
    call Logger%Write( "PolyOrder:      This%PolyOrder  = ", This%PolyOrder )
    allocate(This%cPol(0:This%PolyOrder+1))
    
    do iOrd = 1,This%PolyOrder+1
      READ(Unit, *) This%cPol(iOrd-1)
      call Logger%Write( "Coeffecient Nb ", iOrd, ":   This%cPol(iOrd)  = ", This%cPol(iOrd-1) )
    end do            

  close(Unit)

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_Modified_Morse( This, R, V, dV )

  class(Modified_Morse_DiatomicPotential_Type) ,intent(in)     ::    This
  real(rkp)                                    ,intent(in)     ::    R                         ! Distances between nuclear centers [bohr]
  real(rkp)                                    ,intent(out)    ::    V                         ! Potential energy [hartree]
  real(rkp)                                    ,intent(out)    ::    dV                        ! First derivative of the potential energy wrt the distance [hartree/bohr]
  
  real(rkp)                                                    ::    RTemp, RTemp4, re4
  real(rkp)                                                    ::    poly, dpoly, dydR
  real(rkp)                                                    ::    y, yTemp
  integer                                                      ::    iOrd

  RTemp4 = RTemp**4
  re4    = This%re**4
  y      = (RTemp4-re4)/(RTemp4+re4)
   
  poly  = This%cPol(0)
  dpoly = This%cPol(1)
  yTemp = One
  do iOrd = 1,This%PolyOrder-1
    yTemp =     y * yTemp
    poly  =  poly +            This%cPol(iOrd)   * yTemp
    dpoly = dpoly + (iOrd+1) * This%cPol(iOrd+1) * yTemp
  end do
  yTemp   =     y * yTemp
  poly    =  poly +    This%cPol(This%PolyOrder) * yTemp 

  dpoly   = dpoly * Eight*re4*rTemp**3 / (RTemp4+re4)**2
  dydR    = Four * RTemp**3 * ((RTemp4+re4)-(RTemp4-re4)) / (RTemp4+re4)**2

  V       = This%De * (One - exp(-poly*(RTemp-This%re)))**2  - This%De + This%cPol(This%PolyOrder+1)
  V       =  V * eV_To_Hartree

  dV      = Two * This%De * (One - exp(-poly*(RTemp-This%re))) * (-exp(-poly*(RTemp-This%re))) * (-dpoly*dydR*(RTemp-This%re) - poly)
  dV      = dV * eV_To_Hartree 

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_Modified_Morse( This, R ) result( V )

  class(Modified_Morse_DiatomicPotential_Type) ,intent(in)     ::    This
  real(rkp)                                    ,intent(in)     ::    R                         ! Distances between nuclear centers [bohr]
  real(rkp)                                                    ::    V                         ! Potential energy [hartree]
  
  real(rkp)                                                    ::    RTemp, RTemp4, re4
  real(rkp)                                                    ::    poly
  real(rkp)                                                    ::    y, yTemp
  integer                                                      ::    iOrd

 
!  RTemp  = R * B_To_Ang
  RTemp4 = RTemp**4
  re4    = This%re**4
  y      = (RTemp4-re4)/(RTemp4+re4)

  poly  = This%cPol(0)
  yTemp = One
  do iOrd = 1,This%PolyOrder
    yTemp =     y * yTemp
    poly  =  poly + This%cPol(iOrd) * yTemp
  end do

  V      = This%De * (One - exp(-poly*(RTemp-This%re)))**2  - This%De + This%cPol(This%PolyOrder+1)
  V      = V * eV_To_Hartree

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
