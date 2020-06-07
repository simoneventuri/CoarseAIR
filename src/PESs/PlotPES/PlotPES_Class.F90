! -*-F90-*-
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

Module PlotPES_Class
  
#include "../../qct.inc"

  use Parameters_Module           ,only:  rkp, Zero, Half, One, Two, Pi, Kelvin_To_Hartree, Kcm_To_Hartree, Hartree_To_eV, B_To_Ang, KcmAng_To_HartB
  use Logger_Class                ,only:  Logger, LogLevel_INFO
  use Error_Class                 ,only:  Error

  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Transformation_Class        ,only:  R_to_X, X_to_R, dX_To_dR, dR_To_dX

  implicit none

  private
  public  ::    PlotPES_Type
  
  Type    ,abstract     :: PlotPES_Type
    logical                 ::    Initialized         !< Indicator whether the object is initialized
  contains
    private
    procedure              ,public                          ::    Initialize             =>    PlotPES_Initialize
    procedure              ,public                          ::    DiatPot                =>    PlotPES_DiatPot
    procedure              ,public                          ::    IsoTri                 =>    PlotPES_IsoTri
    procedure              ,public                          ::    Rot3rd                 =>    PlotPES_Rot3rd
    procedure              ,public                          ::    ComputeCuts            =>    PlotPES_ComputeCuts
    procedure              ,public                          ::    EvaluatePoints         =>    PlotPES_EvaluatePoints
    procedure              ,public                          ::    PlotsVargasPaper       =>    PlotPES_PlotsVargasPaper
    procedure              ,public                          ::    StochPESStats          =>    PlotPES_StochPESStats
    procedure              ,public                          ::    GridForScatter         =>    PlotPES_GridForScatter
    procedure              ,public                          ::    ReadPoints             =>    PlotPES_ReadPoints
    procedure              ,public                          ::    TripleGrid             =>    PlotPES_TripleGrid
    procedure              ,public                          ::    DoubleGrid             =>    PlotPES_DoubleGrid
    procedure              ,public                          ::    GridForStochPES        =>    PlotPES_GridForStochPES
    procedure              ,public                          ::    Grid                   =>    PlotPES_Grid
  End Type

  integer   ,parameter    ::    NSpace         = 3
  logical   ,parameter    ::    Formatted      = .True.
  integer                 ::    iSpeTar        = 1
  integer                 ::    iSpePro        = 2
  real(rkp)               ::    MostProbEr                            
  real(rkp)               ::    rVMin_Min      = 2.0d0
  real(rkp)               ::    rVMin_Max      = 2.5d0

  real(rkp) ,save         ::    RConverter     = One
  real(rkp) ,save         ::    VConverter     = One
  real(rkp) ,save         ::    dVConverter    = One
  
  logical   ,parameter    ::    i_Debug_Global = .False.
  
  contains


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_Initialize( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_Initialize" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_DiatPot( This, Input, Collision, NPairs, NAtoms,  i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    ang
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)                                                 ::    dV
  real(rkp)                                                 ::    Angle
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(2)                                ::    hGrid
  real(rkp)                                                 ::    Temp
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA, iP, iMol
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit, Unit1, Unit2
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    t2, t1
  character(4)                                              ::    iAngChar
  character(4)                                              ::    iPChar
  real(rkp)                                                 ::    iAng, RpLong
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_DiatPot" )


  RpLong = 1000.d0


  if (trim(adjustl(Input%UnitDist)) .eq. 'Angstrom') then           
    RConverter  = One         / B_To_Ang
    dVConverter = dVConverter / B_To_Ang
  end if
  
  if (trim(adjustl(Input%UnitPot)) .eq. 'KcalMol') then                                                                           
    VConverter  = One              / Kcm_To_Hartree
    dVConverter = dVConverter / Kcm_To_Hartree
  elseif (trim(adjustl(Input%UnitPot)) .eq. 'ElectronVolt') then                                                                
    VConverter  = One              * Hartree_To_eV
    dVConverter = dVConverter * Hartree_To_eV
  end if
  if (i_Debug_Loc) call Logger%Write( "RConverter  = ", RConverter )
  if (i_Debug_Loc) call Logger%Write( "VConverter  = ", VConverter )
  if (i_Debug_Loc) call Logger%Write( "dVConverter = ", dVConverter )
  
  call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES' )
  if (i_Debug_Loc) call Logger%Write( "Created PlotPES Output Folder" )

  VRef  = 0.0d0
  write(*,*) 'VRef = ', VRef*VConverter, ' eV'


  !if (i_Debug_Loc) call Logger%Write( "Write the Velocity Scaling Factor" )  
  FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/Info.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName )  
    write(Unit,'(A)')        "#       Nb Points R1        Nb Points R2"
    write(Unit,'(2I20)')     Input%NGridPlot(:)
    write(Unit,'(A)')        "#             Min R1              Min R2"
    write(Unit,'(2d20.10)')  Input%MinGridPlot(:)
    write(Unit,'(A)')        "#             Max R1              Max R2"
    write(Unit,'(2d20.10)') Input%MaxGridPlot(:)
  close(Unit)
  
  hGrid = (Input%MaxGridPlot - Input%MinGridPlot) / (Input%NGridPlot - 1)
  write(*,*) Input%NGridPlot


  do iMol = 1,Input%NMolecules

    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/VDiat_From_VDiat_' // trim(adjustl(Collision%MoleculesContainer(iMol)%Molecule%Name)) // '.csv'
    open( File=FileName, NewUnit=Unit1, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName )
    write(Unit1,'(A)') 'Variables = "r1", "V"'

    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/dVDiat_From_VDiat_' // trim(adjustl(Collision%MoleculesContainer(iMol)%Molecule%Name)) // '.csv'
    open( File=FileName, NewUnit=Unit2, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName )
    write(Unit2,'(A)') 'Variables = "r1", "V", "dV"'
      
        
      Rp = RpLong
      do i = 1,Input%NGridPlot(1)
        
        Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 
      
        V = Collision%MoleculesContainer(iMol)%Molecule%DiatPot%DiatomicPotential( Rp(1) * RConverter )  
        write(Unit1,'(f15.6,(A,f15.6))') Rp(1), ',', (V - VRef) * VConverter   

        call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%Compute_Vd_dVd( Rp(1) * RConverter, V, dV )                                          
        write(Unit2,'(f15.6,2(A,f15.6))') Rp(1), ',',  (V - VRef) * VConverter, ',', (dV) * dVConverter                                           
            
      end do

      
    close(Unit1)   

    close(Unit2)   

  end do


  do iP = 1,3!6
    write(iPChar,'(I1)') iP

    do iPES = 1,Input%NPESs
      write(iPESChar,'(I4)') iPES


      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/dVDiat_From_PES' // trim(adjustl(iPESChar)) // '.csv.' // trim(adjustl(iPChar))
      open( File=FileName, NewUnit=Unit1, status='REPLACE', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // FileName )
      write(Unit1,'(A)') 'Variables = "r1", "V", "dVdR"'

        Rp     = RpLong
        do i = 1,Input%NGridPlot(1)
          V    = Zero
          
          Rp(iP) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 
          
          if (Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
            select case(iP)
              case(1)
                Qp = [Zero, Zero, Zero, Rp(iP), Zero, Zero, RpLong, RpLong, Zero, Zero, RpLong, Zero] * RConverter
              case(2)
                Qp = [Zero, Zero, Zero, RpLong, Zero, Zero, Rp(iP)/sqrt(Two), Rp(iP)/sqrt(Two), Zero, Zero, RpLong, Zero] * RConverter
              case(3)
                Qp = [Zero, Zero, Zero, RpLong, Zero, Zero, RpLong, RpLong, Zero, Zero, Rp(iP), Zero] * RConverter
              case(4)
                Qp = [Zero, Zero, Zero, RpLong, Zero, Zero, RpLong, Rp(iP), Zero, Zero, RpLong, Zero] * RConverter
              case(5)
                Qp = [Zero, Zero, Zero, RpLong, Zero, Zero, RpLong-Rp(iP)/sqrt(Two), Rp(iP)/sqrt(Two), Zero, Zero, RpLong, Zero] * RConverter
              case(6)
                Qp = [Zero, Zero, Zero, RpLong, Zero, Zero, Rp(iP), RpLong, Zero, Zero, RpLong, Zero] * RConverter
            end select
          end if


          call Collision%PESsContainer(iPES)%PES%Compute( Rp, Qp, V, dVdR, dVdQ )    
          write(Unit1,'(f15.6,*(A,f15.6))') Rp(iP), ',', (V - VRef) * VConverter, ',', (dVdR(iP)) * dVConverter     
              
        end do

        
      close(Unit1)   


    end do

  end do


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_Grid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_Grid" )
  !i_Debug_Loc   =     Logger%On()

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_GridForStochPES( This, Input, Collision, NPairs, NAtoms,  i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_GridForStochPES" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_DoubleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_DoubleGrid" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_TripleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_TripleGrid" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_ReadPoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_ReadPoints" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_GridForScatter( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_GridForScatter" )
  !i_Debug_Loc   =     Logger%On()
        
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_StochPESStats( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_StochPESStats" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_PlotsVargasPaper( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_PlotsVargasPaper" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_EvaluatePoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_EvaluatePoints" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_ComputeCuts( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_ComputeCuts" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_Rot3rd( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_Rot3rd" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine PlotPES_IsoTri( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(PlotPES_Type)                       ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "PlotPES_IsoTri" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module