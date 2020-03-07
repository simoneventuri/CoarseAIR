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

Module Nb3_PlotPES_Class
  
#include "../../qct.inc"

  use Parameters_Module           ,only:  rkp, Zero, Half, One, Two, Pi, Kelvin_To_Hartree, Kcm_To_Hartree, Hartree_To_eV, B_To_Ang, KcmAng_To_HartB
  use Logger_Class                ,only:  Logger, LogLevel_INFO
  use Error_Class                 ,only:  Error

  use PlotPES_Class               ,only:  PlotPES_Type
  use Input_Class                 ,only:  Input_Type
  use Collision_Class             ,only:  Collision_Type
  use Transformation_Class        ,only:  R_to_X, X_to_R, dX_To_dR, dR_To_dX

  implicit none

  private
  public  ::    Nb3_PlotPES_Type
  
  Type    ,extends(PlotPES_Type)                            ::    Nb3_PlotPES_Type
  contains
    private
    procedure              ,public                          ::    Initialize             =>    Nb3_PlotPES_Initialize
    procedure              ,public                          ::    IsoTri                 =>    Nb3_PlotPES_IsoTri
    procedure              ,public                          ::    Rot3rd                 =>    Nb3_PlotPES_Rot3rd
    procedure              ,public                          ::    ComputeCuts            =>    Nb3_PlotPES_ComputeCuts
    procedure              ,public                          ::    EvaluatePoints         =>    Nb3_PlotPES_EvaluatePoints
    procedure              ,public                          ::    PlotsVargasPaper       =>    Nb3_PlotPES_PlotsVargasPaper
    procedure              ,public                          ::    StochPESStats          =>    Nb3_PlotPES_StochPESStats
    procedure              ,public                          ::    GridForScatter         =>    Nb3_PlotPES_GridForScatter
    procedure              ,public                          ::    ReadPoints             =>    Nb3_PlotPES_ReadPoints
    procedure              ,public                          ::    TripleGrid             =>    Nb3_PlotPES_TripleGrid
    procedure              ,public                          ::    DoubleGrid             =>    Nb3_PlotPES_DoubleGrid
    procedure              ,public                          ::    GridForStochPES        =>    Nb3_PlotPES_GridForStochPES
    procedure              ,public                          ::    Grid                   =>    Nb3_PlotPES_Grid
  End Type

  integer   ,parameter    ::    NSpace         = 3
  logical   ,parameter    ::    Formatted      = .True.
  integer                 ::    iSpeTar        = 1
  integer                 ::    iSpePro        = 2
  real(rkp)               ::    MostProbEr                            
  real(rkp)               ::    rVMin_Min      = 2.0d0
  real(rkp)               ::    rVMin_Max      = 2.5d0
  
  logical   ,parameter    ::    i_Debug_Global = .False.

  real(rkp) ,save         ::    RConverter     = One
  real(rkp) ,save         ::    VConverter     = One
  real(rkp) ,save         ::    dVConverter    = One
  
  contains

!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_Initialize( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_Initialize" )
  !i_Debug_Loc   =     Logger%On()


  if (trim(adjustl(Input%UnitDist)) .eq. 'Angstrom') then           
    RConverter  = One              / B_To_Ang
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

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_Grid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    ang
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    Theta2
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)                                                 ::    dV
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NAtoms*3)                         ::    QpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(2)                                ::    hGrid
  real(rkp)                                                 ::    Temp
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    t2, t1
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_Grid" )
  !i_Debug_Loc   =     Logger%On()

  
  if (trim(adjustl(Input%PESOrDiatFlg)) .eq. 'PES') then
  
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
     
  
    do iPES = 1,Input%NPESs
      write(iPESChar,'(I4)') iPES
      FolderName = 'PES_' // trim(adjustl(iPESChar))
      call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))
  
  
      if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
        RpInf = 1000.d0
        if (Input%PlotPES_OnlyTriatFlg) then 
          call R_to_X(RpInf, QpInf)
          VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, QpInf )
        else
          VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, QpInf )
        end if
      end if
      !VInf  = Zero
      write(*,*) 'VInf = ', VInf*VConverter, ' eV'
      
      
      if (Input%PESZeroRefIntFlg == 0) then
        VRef  = 0.0d0
      else
        if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
          Rp(1)   = rVMin_Min
          Rp(3)   = 1000.0
          VRef    = 1.d10
          do while(Rp(1) < rVMin_Max)
            Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
            call R_to_X(Rp, Qp, Theta=120.d0)
            if (Input%PlotPES_OnlyTriatFlg) then 
              VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Qp ), VRef)
            else
              VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Qp ), VRef)
            end if
            Rp(1) = Rp(1) + 1.d-4
          end do
          VRef  = VRef !- VInf
        end if
      end if
      write(*,*) 'VRef = ', VRef*VConverter, ' eV'
  
  
      if (trim(adjustl(Input%yAxisVar)) .eq. 'Distance') then 

        do iA = 1,size(Input%AnglesPlot,1)
        
          call cpu_time ( t1 )
            
          FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESFromGrid.csv.' // trim(adjustl(Input%AnglesPlotChar(iA)))
          open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
              write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
            else
              write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
            end if
            if (Status/=0) call Error( "Error opening file: " // FileName )
            
            
            iTot = 1
            do i = 1,Input%NGridPlot(1)
              
              Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 


              do j = 1,Input%NGridPlot(2)
                
                Rp(3)   = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) )     
                Theta2  = Input%AnglesPlot(iA) / 180.d0 * Pi
                Rp(2)   = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Theta2) ) 
                call R_to_X(Rp, Qp, Theta=Theta2)
                if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

                  if (Input%PlotPES_OnlyTriatFlg) then 
                    !V = Collision%PESsContainer(iPES)%PES%DiatPotential( Rp * RConverter )
                    V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Qp )
                  else
                    V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Qp )
                  end if
                  !Temp = (V - VRef) * VConverter / abs((V - VRef) * VConverter)
                  if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                    write(Unit,'(es17.6E3,3(A,es17.6E3))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter!Temp*max(abs((V - VRef) * VConverter), 1.d-90 )
                  end if
                  
                elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
                  
                  call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )     
                  if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                           
                    write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
                  end if
                  
                end if  
              
                iTot = iTot + 1
                  
              end do
              
            end do
            
          close(Unit)   
          
          call cpu_time ( t2 )
          write(*,*) 'Time for PES Calculations = ', t2-t1
          
        end do      
          
        
      elseif (trim(adjustl(Input%yAxisVar)) .eq. 'Angle') then  
      
        FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESDistVsAngles.csv'
        open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
            write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "Angle", "E"'
          else
            write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "Angle", "E", "dEdR1", "dEdR2", "dEdR3"'
          end if
          if (Status/=0) call Error( "Error opening file: " // FileName )
          
          iTot = 1
          do i = 1,Input%NGridPlot(1)!radiuses (constrained to be equal)
            !this case r1=r2 and angle is varying

            Rp(1) = 1.4d0!( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 
            Rp(3) = 1.7d0!Rp(1)
          
            do j = 1,Input%NGridPlot(2)!angle
            
              Angle = Input%MinGridPlot(2) + hGrid(2) * (j-1)
              Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Angle/ 180.d0 * Pi) ) 
                
              call R_to_X(Rp, Qp, Theta=Angle)
              if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

                if (Input%PlotPES_OnlyTriatFlg) then 
                  V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Qp )
                else
                  V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Qp )
                end if
                !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                  write(Unit,'(es15.6,4(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Angle, ',', (V - VRef) * VConverter
                !end if
                
              elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
                
                call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Qp, V, dVdR, dVdQ )    
                !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                       
                  write(Unit,'(es15.6,7(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Angle, ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
                !end if
              
              end if  
             
              iTot = iTot + 1
            end do
          
          end do
        
        close(Unit)   
        
          
      end if ! Distance/Angle
    
    end do


  elseif (trim(adjustl(Input%PESOrDiatFlg)) .eq. 'Diatomic') then
    
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

  
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/DiatomicPotential.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // FileName )
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then    
        write(Unit,'(A)') 'Variables = "r1", "V"'
      elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then 
        write(Unit,'(A)') 'Variables = "r1", "V", "dV"'
      end if
    
      
      iTot = 1
      do i = 1,Input%NGridPlot(1)
        
        Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 

        do j = 1,1!,2                          !  Input%NGridPlot(2) is not needed
          
          Rp(3) = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) )
          
          Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * 0.0000001_rkp ) 
          
          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then    
               
            V = Collision%Species(1)%Diapot%DiatomicPotential( Rp(1) * RConverter )  
            write(Unit,'(f15.6,(A,f15.6))') Rp(1), ',', (V - VRef) * VConverter   
       
          elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then 
          
            call Collision%Species(1)%Diapot%Compute_Vd_dVd( Rp(1) * RConverter, V, dV )                                          
            write(Unit,'(f15.6,2(A,f15.6))') Rp(1), ',',  (V - VRef) * VConverter, ',', (dV) * dVConverter                                           
                         
          end if     
              
          iTot = iTot + 1
        end do
      
      end do
      
    close(Unit)   
    
    
  end if ! PES/DIAT

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_GridForStochPES( This, Input, Collision, NPairs, NAtoms,  i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
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
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    t2, t1
  character(4)                                              ::    iAngChar
  real(rkp)                                                 ::    iAng
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_GridForStochPES" )
  !i_Debug_Loc   =     Logger%On()
  

  do iA = 1,size(Input%AnglesPlot,1)
    iAng = Input%AnglesPlot(iA)
    write(iAngChar,'(I4)') floor(iAng)
    FolderName = 'Ang' // trim(adjustl(iAngChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) )
  end do
  
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
    

  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES


    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf  = Zero
    write(*,*) 'VInf = ', VInf*VConverter, ' eV'
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      Rp(1)   = rVMin_Min
      Rp(3)   = 1000.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    end if
    write(*,*) 'VRef = ', VRef*VConverter, ' eV'


    do iA = 1,size(Input%AnglesPlot,1)
      iAng = Input%AnglesPlot(iA)
      write(iAngChar,'(I4)') floor(iAng)
        
      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/Ang' // trim(adjustl(iAngChar)) // '/PESFromGrid.csv.' // trim(adjustl(iPESChar)) 
      open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
        if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
          write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
        else
          write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
        end if
        if (Status/=0) call Error( "Error opening file: " // FileName )
        
        
        iTot = 1
        do i = 1,Input%NGridPlot(1)
          
          Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 


          do j = 1,Input%NGridPlot(2)
            
            Rp(3) = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) )             
            Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Input%AnglesPlot(iA) / 180.d0 * Pi) ) 
            
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

              if (Input%PlotPES_OnlyTriatFlg) then 
                !V = Collision%PESsContainer(iPES)%PES%DiatPotential( Rp * RConverter )
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              !Temp = (V - VRef) * VConverter / abs((V - VRef) * VConverter)
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es17.6E3,3(A,es17.6E3))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter!Temp*max(abs((V - VRef) * VConverter), 1.d-90 )
              end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
              call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )     
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                           
                write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
              end if
              
            end if  
          
            iTot = iTot + 1
              
          end do
          
        end do
        
      close(Unit)   
      
      call cpu_time ( t2 )
      
    end do      
      
  end do
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_DoubleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    R1Min    
  real(rkp)                                                 ::    hMin    
  real(rkp)                                                 ::    R1Max  
  real(rkp)                                                 ::    hMax, h  
  real(rkp)                                                 ::    R1Delta
  real(rkp)                                                 ::    hDelta 
  real(rkp)                                                 ::    AngMin
  real(rkp)                                                 ::    AngMax
  real(rkp)                                                 ::    AngDelta
  real(rkp)                                                 ::    Ang  
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf, VRef
  real(rkp)    ,dimension(2)                                ::    hGrid
  integer                                                   ::    i, j, iAngle, iR, jR
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  logical                                                   ::    i_Debug_Loc
  
  real(rkp)    ,dimension(12)                               ::    AngA
  real(rkp)    ,dimension(13)                               ::    RA1
  real(rkp)    ,dimension(13)                               ::    RA2
      
  real(rkp)    ,dimension(18)                               ::    AngB
  real(rkp)    ,dimension(12)                               ::    RB1
  real(rkp)    ,dimension(19)                               ::    RB2
  real(rkp)    ,dimension(3)                                ::    RpTemp
  real(rkp)                                                 ::    VTemp
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_DoubleGrid" )
  !i_Debug_Loc   =     Logger%On()
  

  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))
    
    
    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf = 0.d0
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      Rp(1)   = rVMin_Min
      Rp(3)   = 1000.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    end if
    

    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/PESFromDoubleGrid.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
      else
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
      end if
      if (Status/=0) call Error( "Error opening file: " // FileName )
      
      
      AngA = [60.d0, 70.d0, 80.d0, 90.d0, 100.d0, 110.d0, 120.d0, 130.d0, 140.d0, 150.d0, 160.d0, 170.d0]
      RA1  = [3.d0, 3.25d0, 3.5d0, 3.75d0, 4.d0, 4.5d0, 5.d0, 5.5d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0]
      RA2  = [3.d0, 3.25d0, 3.5d0, 3.75d0, 4.d0, 4.5d0, 5.d0, 5.5d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0]
       
      do iAngle = 1,size(AngA)
        Ang = AngA(iAngle)
        
        do iR = 1,size(RA1)
          Rp(1) = RA1(iR)
        
          do jR = 1,size(RA2)
            Rp(3) = RA2(jR)
          
             if (Rp(3) >= Rp(1)) then
          
              Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( Ang / 180.d0 * pi ) )
            
              if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
      
                if (Input%PlotPES_OnlyTriatFlg) then 
                  V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
                else
                  V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
                end if
                if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                  write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
                end if
                
              elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
                
                call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )  
                if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                         
                  write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter
                end if
                
              end if  
              
            end if
          
          end do
          
        end do
        
      end do
      
      
      AngB = [50.d0, 60.d0, 70.d0, 80.d0, 90.d0, 95.d0, 100.d0, 105.d0, 110.d0, 115.d0, 120.d0, 130.d0, 140.d0, 150.d0, 160.d0, 165.d0, 170.d0, 175.d0]
      RB1  = [1.6d0, 1.7d0, 1.8d0, 1.9d0, 2.d0, 2.15d0, 2.3d0, 2.5d0, 2.8d0, 3.d0, 3.5d0, 4.d0]
      RB2  = [1.6d0, 1.7d0, 1.8d0, 1.9d0, 2.d0, 2.15d0, 2.3d0, 2.5d0, 2.8d0, 3.d0, 3.5d0, 4.d0, 4.5d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0]
      
      do iAngle = 1,size(AngB)
        Ang = AngB(iAngle)
        
        do iR = 1,size(RB1)
          Rp(1) = RB1(iR)
        
          do jR = 1,size(RB2)
            Rp(3) = RB2(jR)
            
            if (Rp(3) >=Rp(1)) then
          
              Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( Ang / 180.d0 * pi ) )
            
              if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
      
                if (Input%PlotPES_OnlyTriatFlg) then 
                  V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
                else
                  V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
                end if
                if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                  write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
                end if
                
              elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           

                call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )    
                if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                       
                  write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter
                end if
                
              end if  
              
            end if
          
          end do
          
        end do
        
      end do
      
      
    close(Unit)
    
  end do

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_TripleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    R1Min    
  real(rkp)                                                 ::    R3Min    
  real(rkp)                                                 ::    R1Max  
  real(rkp)                                                 ::    R3Max    
  real(rkp)                                                 ::    R1Delta
  real(rkp)                                                 ::    R3Delta 
  real(rkp)                                                 ::    AngMin
  real(rkp)                                                 ::    AngMax
  real(rkp)                                                 ::    AngDelta
  real(rkp)                                                 ::    Ang  
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)                                                 ::    dV
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf, VRef
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)    ,dimension(2)                                ::    hGrid
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_TripleGrid" )
  !i_Debug_Loc   =     Logger%On()

  
  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))
    
    
    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf = 0.d0
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      Rp(1)   = rVMin_Min
      Rp(3)   = 1000.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    end if
    
    
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/PESFromTripleGrid.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
      else
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
      end if
      if (Status/=0) call Error( "Error opening file: " // FileName )
    
        
      R1Min    =   1.6d0
      R3Min    =   6.5d0
      R1Max    =   5.5d0
      R3Max    =  18.0d0
      R1Delta  =   0.3d0
      R3Delta  =   0.3d0
      AngMin   =  30.0d0
      AngMax   = 180.0d0
      AngDelta =   5.0d0
      
      Rp(1) = R1Min
      do while (Rp(1) <= R1Max)
        
        Rp(3) = R3Min
        do while (Rp(3) <= R3Max)

          Ang = AngMin
          do while (Ang < AngMax)
        
            Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( Ang / 180.d0 * pi ) )
              
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
    
              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
              end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
              call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )             
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                              
                write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter
              end if
              
            end if  
            
            Ang = Ang + AngDelta
          end do
        
          Rp(3) = Rp(3) + R3Delta
        end do
        
        Rp(1) = Rp(1) + R1Delta
      end do
      
      
      R1Min    =   1.6d0
      R3Min    =   1.6d0
      R1Max    =   7.5d0
      R3Max    =   7.5d0
      R1Delta  =   0.2d0
      R3Delta  =   0.2d0
      AngMin   =  30.0d0
      AngMax   = 180.0d0
      AngDelta =   5.0d0
      
      Rp(1) = R1Min
      do while (Rp(1) <= R1Max)
        
        Rp(3) = R3Min
        do while (Rp(3) <= R3Max)
        
          Ang = AngMin
          do while (Ang < AngMax)
            Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( Ang / 180.d0 * pi ) )
              
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
    
              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
              end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           

              call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )      
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                     
                write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter
              end if
              
            end if  
            
            Ang = Ang + AngDelta
          end do
        
          Rp(3) = Rp(3) + R3Delta
        end do
        
        Rp(1) = Rp(1) + R1Delta
      end do
      
      
      R1Min    =   1.6d0
      R3Min    =   1.6d0
      R1Max    =  15.0d0
      R3Max    =  15.0d0
      R1Delta  =   1.2d0
      R3Delta  =   1.2d0
      AngMin   =  30.0d0
      AngMax   = 180.0d0
      AngDelta =   5.0d0
      
      Rp(1) = R1Min
      do while (Rp(1) <= R1Max)
        
        Rp(3) = R3Min
        do while (Rp(3) <= R3Max)
        
          Ang = AngMin
          do while (Ang < AngMax)
            Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( Ang / 180.d0 * pi ) )
              
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
    
              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
              end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then 
                      
              call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )      
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                     
                write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter
              end if
              
            end if  
            
            Ang = Ang + AngDelta
          end do
        
          Rp(3) = Rp(3) + R3Delta
        end do
        
        Rp(1) = Rp(1) + R1Delta
      end do
      
    
    close(Unit) 


  end do


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_ReadPoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    R1Min    
  real(rkp)                                                 ::    hMin    
  real(rkp)                                                 ::    R1Max  
  real(rkp)                                                 ::    hMax, h  
  real(rkp)                                                 ::    R1Delta
  real(rkp)                                                 ::    hDelta 
  real(rkp)                                                 ::    AngMin
  real(rkp)                                                 ::    AngMax
  real(rkp)                                                 ::    AngDelta
  real(rkp)                                                 ::    Ang  
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf, VRef
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NAtoms*3)                         ::    QpInf
  real(rkp)    ,dimension(NAtoms*3)                         ::    QpDiat
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    hDiat
  real(rkp)                                                 ::    hDiatMax
  real(rkp)    ,dimension(2)                                ::    hGrid
  integer                                                   ::    i, j, iAngle, iR, jR
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName, FileNameE
  integer                                                   ::    Unit, UnitRead, UnitReadE
  integer                                                   ::    Status, StatusRead, StatusReadE
  character(:)                 ,allocatable                 ::    FolderName
  logical                                                   ::    i_Debug_Loc
  
  real(rkp)    ,dimension(12)                               ::    AngA
  real(rkp)    ,dimension(13)                               ::    RA1
  real(rkp)    ,dimension(13)                               ::    RA2
      
  real(rkp)    ,dimension(18)                               ::    AngB
  real(rkp)    ,dimension(12)                               ::    RB1
  real(rkp)    ,dimension(19)                               ::    RB2
  real(rkp)    ,dimension(3)                                ::    RpTemp
  real(rkp)                                                 ::    VTemp
  real(rkp)                                                 ::    EPoints, Temp, RNew
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_ReadPoints" )
  !i_Debug_Loc   =     Logger%On()

  
  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! IF PES WORKS WITH RELATIVE DISTANCES ...
    if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
    
      RpInf = 2.d0
      if (Input%PlotPES_OnlyTriatFlg) then 
        VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf)
      else
        VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
      end if
      write(*,*) 'VInf = ', VInf*VConverter, ' eV'
      !VInf = 0.d0
      
      
      if (Input%PESZeroRefIntFlg == 0) then
        VRef  = 0.0d0
      else
        Rp(1)   = rVMin_Min
        Rp(3)   = 100.0
        VRef    = 1.d10
        do while(Rp(1) < rVMin_Max)
          Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
          if (Input%PlotPES_OnlyTriatFlg) then 
            VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
          else
            VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
          end if
          Rp(1) = Rp(1) + 1.d-4
        end do
        VRef  = VRef !- VInf
      end if
      write(*,*) 'VRef = ', VRef, ' eV'


      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/PESFromReadPoints.csv'
      open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
        !write(Unit,'(A)')        "#                 R1                   R2                   R3                    V"
        if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
          write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
        else
          write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
        end if
        if (Status/=0) call Error( "Error opening file: " // FileName )   
        
        if ((Input%NPESs > 1) .and. (.not. Input%StochPESFlg)) then
          FileName = trim(adjustl(Input%OutputDir)) // '/RSampled.csv.' // trim(adjustl(iPESChar))
        else
          FileName = trim(adjustl(Input%OutputDir)) // '/RSampled.csv.1'
        end if
        write(*,*) 'Reading Points from the File: ', FileName
        open( File=FileName, NewUnit=UnitRead, status='OLD', iostat=StatusRead ) 
        if (StatusRead /= 0) then 
           write(*,'(A)')"File "//trim(FileName)//" NOT found!"
           call Error('ReadPoints: File ./RSampled.csv not Found!')
        endif
        
  !      FileNameE = trim(adjustl(Input%OutputDir)) // '/ESampled.csv.' // trim(adjustl(iPESChar))
  !      write(*,*) 'Reading E Points from the File: ', FileNameE
  !      open( File=FileNameE, NewUnit=UnitReadE, status='OLD', iostat=StatusReadE )
          
          do while (StatusRead==0)
            read(UnitRead,*,iostat=StatusRead) Rp(1), Rp(2), Rp(3)          
            
            if (StatusRead/=0) exit
            
  !          EPoints = 0.d0
  !          read(UnitReadE,'(e12.6,A,e12.6)',iostat=StatusReadE) EPoints, Temp
            
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
      
              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
              !end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
              call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )   
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then  
                write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter                                  
              end if
              
            end if  
            
          end do

        close(UnitRead)
        
      close(Unit)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! IF PES WORKS WITH CARTESIAN COORDINATES ...
    else

      !QpInf = [1.d3,1.d3,0.d0, -1.d3,1.d3,0.d0, -1.d3,-1.d3,0.d0, 1.d3,-1.d3,0.d0] 
      !if (Input%PlotPES_OnlyTriatFlg) then 
      !  VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( Zero*QpInf, QpInf)
      !else
      !  VInf  = Collision%PESsContainer(iPES)%PES%Potential( Zero*QpInf, QpInf )
      !end if
      !write(*,*) 'VInf = ', VInf*VConverter, ' eV'
      !VInf = 0.d0
      
      
      if (Input%PESZeroRefIntFlg == 0) then
        VRef  = 0.0d0
      else
        hDiat    = 1.d-4 
        hDiatMax = 4.d0
        VRef     = 1.d10
        do while(hDiat < hDiatMax)
          QpDiat = [1.d3,1.d3,0.d0, 1.d3-hDiat,1.d3,0.d0, -1.d3,-1.d3,0.d0, 1.d3,-1.d3,0.d0] 
          if (Input%PlotPES_OnlyTriatFlg) then 
            VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Zero*QpDiat, QpDiat ), VRef)
          else
            VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Zero*QpDiat, QpDiat ), VRef)
          end if
          hDiat = hDiat + 1.d-4
        end do
        VRef  = VRef !- VInf
      end if
      write(*,*) 'VRef = ', VRef, ' eV'

      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/PESFromReadPoints.csv'
      open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
        !write(Unit,'(A)')        "#                 R1                   R2                   R3                    V"
        if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
          write(Unit,'(A)') 'Variables = "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "E"'
        else
          write(Unit,'(A)') 'Variables = "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "E", "dEdx1", "dEdy1", "dEdz1", "dEdx2", "dEdy2", "dEdz2", "dEdx3", "dEdy3", "dEdz3", "dEdx4", "dEdy4", "dEdz4"'
        end if
        if (Status/=0) call Error( "Error opening file: " // FileName )   
        
        if ((Input%NPESs > 1) .and. (.not. Input%StochPESFlg)) then
          FileName = trim(adjustl(Input%OutputDir)) // '/QSampled.csv.' // trim(adjustl(iPESChar))
        else
          FileName = trim(adjustl(Input%OutputDir)) // '/QSampled.csv.1'
        end if
        write(*,*) 'Reading Points from the File: ', FileName
        open( File=FileName, NewUnit=UnitRead, status='OLD', iostat=StatusRead ) 
        
        
  !      FileNameE = trim(adjustl(Input%OutputDir)) // '/ESampled.csv.' // trim(adjustl(iPESChar))
  !      write(*,*) 'Reading E Points from the File: ', FileNameE
  !      open( File=FileNameE, NewUnit=UnitReadE, status='OLD', iostat=StatusReadE )
          
          do while (StatusRead==0)
            read(UnitRead,*,iostat=StatusRead) Qp
            
            if (StatusRead/=0) exit
            
  !          EPoints = 0.d0
  !          read(UnitReadE,'(e12.6,A,e12.6)',iostat=StatusReadE) EPoints, Temp
            
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
      
              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Zero*Qp, Qp * RConverter )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Zero*Qp, Qp * RConverter )
              end if
              !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es15.6,*(",",es15.6))') Qp(:), (V - VRef) * VConverter
              !end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
              !call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ ) 
              call Collision%PESsContainer(iPES)%PES%Compute( Zero*Rp, Qp * RConverter, V, dVdR, dVdQ ) 
              
              if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then  
                !write(Unit,'(es15.6,*(",",es15.6))') Qp(:), (V - VRef) * VConverter, dVdR(:) * dVConverter    
                write(Unit,'(es15.6,*(",",es15.6))') Qp(:), (V - VRef) * VConverter, dVdQ(1:9) * dVConverter    

              end if
              
            end if  
            
          end do

        close(UnitRead)
        
      close(Unit)

    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  end do

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_GridForScatter( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  real(rkp)                                                 ::    R1Min    
  real(rkp)                                                 ::    R3Min    
  real(rkp)                                                 ::    R1Max  
  real(rkp)                                                 ::    R3Max    
  real(rkp)                                                 ::    R1Delta
  real(rkp)                                                 ::    R3Delta 
  real(rkp)                                                 ::    Ang  
  character(3)                                              ::    AngChar
  real(rkp)    ,dimension(9)                                ::    AngVec   
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    RR
  real(rkp)                                                 ::    V
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)    ,dimension(2)                                ::    hGrid
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iAng
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  character(:)                 ,allocatable                 ::    FolderName
  integer                                                   ::    Unit
  integer                                                   ::    Status
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_GridForScatter" )
  !i_Debug_Loc   =     Logger%On()

  
  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    
    
    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf = 0.d0
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      Rp(1)   = rVMin_Min
      Rp(3)   = 1000.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    end if
    
    
!    R1Min    =   1.6d0
!    R3Min    =   1.6d0
!    R1Max    =   18.d0
!    R3Max    =   18.d0
!    R1Delta  =   0.1d0
!    R3Delta  =   0.1d0
!    AngVec   =   0.0d0
    R1Min    =   1.5d0
    R3Min    =   1.5d0
    R1Max    =   10.d0
    R3Max    =   10.d0
    R1Delta  =   0.1d0
    R3Delta  =   0.1d0
    AngVec   =   0.0d0
    AngVec   = [60.0, 80.0, 100.0, 110.0, 120.0, 130.0, 140.0, 160.0, 180.0]
    
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESForScatterPlot.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
      else
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
      end if
      if (Status/=0) call Error( "Error opening file: " // FileName )
  
      do iAng = 1,size(AngVec,1)
        Ang      = AngVec(iAng)
        write(AngChar,'(I3)') int(Ang)
          
        Rp(1) = R1Min
        do while (Rp(1) <= R1Max)
          
          Rp(3) = R3Min
          do while (Rp(3) <= R3Max)
            
            Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( Ang / 180.d0 * pi ) )
            RR    = dsqrt(  Rp(1)**Two + Rp(2)**Two + Rp(3)**Two )
              
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
              !end if
              
            elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
              call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )      
              !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                     
                write(Unit,'(es15.6,6(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
              !end if
            end if  
          
            Rp(3) = Rp(3) + R3Delta
          end do
          
          Rp(1) = Rp(1) + R1Delta
        end do
          
      end do
        
    close(Unit)
    
        
  end do
        
        
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_StochPESStats( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
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
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)    ,dimension(2)                                ::    hGrid
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    VRefAvg
  real(rkp)                                                 ::    VSum
  real(rkp)                                                 ::    VSqrSum
  real(rkp)                                                 ::    VMean
  real(rkp)                                                 ::    VSD
  real(rkp)                                                 ::    VMinus
  real(rkp)                                                 ::    VPlus
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_StochPESStats" )
  !i_Debug_Loc   =     Logger%On()
  

  hGrid = (Input%MaxGridPlot - Input%MinGridPlot) / (Input%NGridPlot - 1)  

  VRefAvg = 0.d0
  do iPES = 1,Input%NPESs


    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf  = Zero


    if (Input%PESZeroRefIntFlg == 0) then
      Rp(1)   = rVMin_Min
      Rp(3)   = 1000.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    else
      VRef  = VInf
    end if
   
    
  end do
    
    
  do iA = 1,size(Input%AnglesPlot,1)   
  
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/StochPESStatistics.csv.' // trim(adjustl(Input%AnglesPlotChar(iA)))
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // FileName )
      !write(Unit,'(A)') '#     Nb PES Samples'
      !write(Unit,'(I20)') Input%NPESs 
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'R1, R2, R3, EMean, ESD, EMinus, EPlus'
      !else
      !  write(Unit,'(A)') '# Variables = "dEdR1", "dEdR2", "dEdR3"'
      end if
  
          
      do i = 1,Input%NGridPlot(1)
        Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 

        do j = 1,i
          Rp(3) = ( Input%MinGridPlot(1) + hGrid(1) * (j-1) )             
          Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Input%AnglesPlot(iA) / 180.d0 * Pi) ) 
          
          VSum    = 0.d0
          VSqrSum = 0.d0
          do iPES = 1,Input%NPESs
          
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

              if (Input%PlotPES_OnlyTriatFlg) then 
                V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
              else
                V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
              end if
              !if ( ((V - Vinf/VConverter) <= Input%EnergyCutOff) ) then
                !write(Unit,'(es12.6)') max((V - VRef) * VConverter, 1.d-99)
              !end if
              
              VSum    = VSum    +  (V - VRef) * VConverter
              VSqrSum = VSqrSum + ((V - VRef) * VConverter)**2.d0 
              
            !elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
            !  call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, V, dVdR )                                          
            
            end if  
          
          end do
          
          VMean  =  VSum    / Input%NPESs
          VSD    = (VSqrSum / Input%NPESs - VMean**2)**(5.d-1)
          VPlus  = VMean + 3.d0*VSD
          VMinus = VMean - 3.d0*VSD
          write(Unit,'(es15.4e3,6(A,es15.4e3))') Rp(1), ',', Rp(2), ',', Rp(3), ',', VMean, ',', VSD, ',', VMinus, ',', VPlus  
            
        end do
        
      end do
      
      
    close(Unit)    
  end do
  
  
  FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/GridForStochPESInfo.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName )  
    write(Unit,'(A)')        "#     Nb PES Samples"
    write(Unit,'(I20)')      Input%NPESs 
    write(Unit,'(A)')        "#             Min R1"
    write(Unit,'(e20.10)')   Input%MinGridPlot(1)
    write(Unit,'(A)')        "#             Max R1"
    write(Unit,'(e20.10)')   Input%MaxGridPlot(1)
    write(Unit,'(A)')        "#       Nb Points R1"
    write(Unit,'(I20)')      Input%NGridPlot(1)
    write(Unit,'(A)')        "#         Spacing R1"
    write(Unit,'(e20.10)')   hGrid(1)
    write(Unit,'(A)')        "#          Nb Angles"
    write(Unit,'(I20)')      size(Input%AnglesPlot,1)
    write(Unit,'(A)')        "#             Angles"
    do iA = 1,size(Input%AnglesPlot,1)   
      write(Unit,'(e20.10)') Input%AnglesPlot(iA)
    end do
  close(Unit) 
   
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! !________________________________________________________________________________________________________________________________!
! Subroutine GridForStochPES( This, Input, Collision, i_Debug )

!   class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
!  type(Input_Type)                          ,intent(in)     ::    Input
!   type(Collision_Type)                      ,intent(in)     ::    Collision
!   logical                         ,optional ,intent(in)     ::    i_Debug

!   real(rkp)                                                 ::    ang
!   real(rkp)                                                 ::    beta
!   real(rkp)                                                 ::    dist
!   real(rkp)    ,dimension(3)                                ::    Rp
!   real(rkp)                                                 ::    V
!   real(rkp)                                                 ::    dV
!   real(rkp)    ,dimension(3)                                ::    dVdR
!   real(rkp)                                                 ::    Angle
!   real(rkp)    ,dimension(3)                                ::    RpInf
!   real(rkp)                                                 ::    VInf
!   real(rkp)                                                 ::    VRef
!   real(rkp)    ,dimension(2)                                ::    hGrid
!   integer                                                   ::    i, j
!   integer                                                   ::    iTot
!   integer                                                   ::    iA
!   integer                                                   ::    iPES
!   character(4)                                              ::    iPESChar
!   character(:)                 ,allocatable                 ::    FileName  
!   integer                                                   ::    Unit
!   integer                                                   ::    Status
!   character(:)                 ,allocatable                 ::    FolderName
!   real(rkp)                                                 ::    VRefAvg

!   logical                                                   ::    i_Debug_Loc
  
!   i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
!   if (i_Debug_Loc) call Logger%Entering( "GridForStochPES" )
!   !i_Debug_Loc   =     Logger%On()
  
      
!   if (trim(adjustl(Input%UnitDist)) .eq. 'Angstrom') then           
!     RConverter  = One         / B_To_Ang
!     dVConverter = dVConverter / B_To_Ang
!   end if
  
!   if (trim(adjustl(Input%UnitPot)) .eq. 'KcalMol') then                                                                           
!     VConverter  = One         / Kcm_To_Hartree
!     dVConverter = dVConverter / Kcm_To_Hartree
!   elseif (trim(adjustl(Input%UnitPot)) .eq. 'ElectronVolt') then                                                                
!     VConverter  = One         * Hartree_To_eV
!     dVConverter = dVConverter * Hartree_To_eV
!   end if
  
!   call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES' )
    
!   hGrid = (Input%MaxGridPlot - Input%MinGridPlot) / (Input%NGridPlot - 1)  

!   VRefAvg = 0.d0
!   do iPES = 1,Input%NPESs


!     RpInf = 1000.d0
!     if (Input%PlotPES_OnlyTriatFlg) then 
!       VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf )
!     else
!       VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf )
!     end if
!     !VInf = 0.d0
    
    
!     if (Input%PESZeroRefIntFlg == 1) then
!       Rp(1)   = rVMin_Min
!       Rp(3)   = 1000.0
!       VRef    = 1.d10
!       do while(Rp(1) < rVMin_Max)
!         Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
!         if (Input%PlotPES_OnlyTriatFlg) then 
!           VRef  = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp )
!         else
!           VRef  =Collision%PESsContainer(iPES)%PES%Potential( Rp )
!         end if
!         Rp(1) = Rp(1) + 1.d-4
!       end do
!       VRefAvg = VRefAvg + (VRef )!- VInf)
!     end if
    
!   end do
!   VRefAvg = VRefAvg / Input%NPESs
    
    
!   do iA = 1,size(Input%AnglesPlot,1)   
      
!     FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/StochPESFromGrid.csv.' // trim(adjustl(Input%AnglesPlotChar(iA)))
!     open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
!       write(Unit,'(A)') '#     Nb PES Samples'
!       write(Unit,'(I20)') Input%NPESs 
!       if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
!         write(Unit,'(A)') '# Variable = "E"'
!       else
!         write(Unit,'(A)') '# Variables = "dEdR1", "dEdR2", "dEdR3"'
!       end if
!       if (Status/=0) call Error( "Error opening file: " // FileName )
  
          
!       do i = 1,Input%NGridPlot(1)
!         Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 

!         do j = 1,i
!           Rp(3) = ( Input%MinGridPlot(1) + hGrid(1) * (j-1) )             
!           Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Input%AnglesPlot(iA) / 180.d0 * Pi) ) 
          
!           do iPES = 1,Input%NPESs
          
!             if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

!               if (Input%PlotPES_OnlyTriatFlg) then 
!                 V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter )
!               else
!                 V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter )
!               end if
!               !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
!                 write(Unit,*) (V - VRef) * VConverter
!               !end if
              
!             elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
!               call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, V, dVdR )                                          
!               write(Unit,*) (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
            
!             end if  
          
!           end do
!           write(Unit,'(A,es10.4,3(A,es10.4))') '# ', Rp(1), ',', Rp(2), ',', Rp(3), ',', Input%AnglesPlot(iA)
          
!         end do
        
!       end do
      
      
!     close(Unit)    
!   end do
  
  
!   FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/GridForStochPESInfo.dat'
!   open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
!     if (Status/=0) call Error( "Error opening file: " // FileName )  
!     write(Unit,'(A)')        "#     Nb PES Samples"
!     write(Unit,'(I20)')      Input%NPESs 
!     write(Unit,'(A)')        "#             Min R1"
!     write(Unit,'(e20.10)')   Input%MinGridPlot(1)
!     write(Unit,'(A)')        "#             Max R1"
!     write(Unit,'(e20.10)')   Input%MaxGridPlot(1)
!     write(Unit,'(A)')        "#       Nb Points R1"
!     write(Unit,'(I20)')      Input%NGridPlot(1)
!     write(Unit,'(A)')        "#         Spacing R1"
!     write(Unit,'(e20.10)')   hGrid(1)
!     write(Unit,'(A)')        "#          Nb Angles"
!     write(Unit,'(I20)')      size(Input%AnglesPlot,1)
!     write(Unit,'(A)')        "#             Angles"
!     do iA = 1,size(Input%AnglesPlot,1)   
!       write(Unit,'(e20.10)') Input%AnglesPlot(iA)
!     end do
!   close(Unit) 
   
  
!   if (i_Debug_Loc) call Logger%Exiting

! End Subroutine
! !--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_PlotsVargasPaper( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
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
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    hGrid
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(4)                                              ::    AngleChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    VRefAvg
  real(rkp)    ,dimension(8)                                ::    AnglesVec
  real(rkp)                                                 ::    FixedValue, FixedValue1, FixedValue2
  integer                                                   ::    NGrid
  real(rkp)                                                 ::    RStart, REnd
  real(rkp)                                                 ::    AlphaStart, AlphaEnd
  real(rkp)    ,dimension(:)   ,allocatable                 ::    SumVec
  real(rkp)    ,dimension(:)   ,allocatable                 ::    SumSqrVec
  real(rkp)    ,dimension(:)   ,allocatable                 ::    MeanVec
  real(rkp)    ,dimension(:)   ,allocatable                 ::    SDVec
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_PlotsVargasPaper" )
  !i_Debug_Loc   =     Logger%On()
  

  VRefAvg = 0.d0
  do iPES = 1,Input%NPESs


    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf = 0.d0
    
    
    if (Input%PESZeroRefIntFlg == 1) then
      Rp(1)   = rVMin_Min
      Rp(3)   = 1000.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRefAvg = VRefAvg + (VRef )!- VInf)
    end if
    
  end do
  VRefAvg = VRefAvg / Input%NPESs
  
  
  AnglesVec  = [110.d0, 120.d0, 130.d0, 140.d0, 150.d0, 160.d0, 170.d0, 180.d0]
  FixedValue = 1.09768
  NGrid      = 301
  RStart     = 0.8
  REnd       = 3.2
  hGrid      = (REnd - RStart) / (NGrid-1)
  if (Input%NPESs > 1) then
    allocate(SumVec(NGrid))
    allocate(SumSqrVec(NGrid))
    allocate(MeanVec(NGrid))
    allocate(SDVec(NGrid))
  end if
   
  do iA = 1,size(AnglesVec,1)   
    write(AngleChar,'(I4)') int(AnglesVec(iA))
      
    do iPES = 1,Input%NPESs  
      write(iPESChar,'(I4)') iPES
      
      if (Input%NPESs > 1) then
        FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/FixedR1_' // trim(adjustl(AngleChar)) // '.csv.' // trim(adjustl(iPESChar))
      else
        FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/FixedR1_' // trim(adjustl(AngleChar)) // '.csv'
      end if
      open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
        write(Unit,'(A)') 'R1, E'
            
        do i = 1,NGrid
          Rp(1) = ( RStart + hGrid * (i-1) ) 
          Rp(3) = FixedValue       
          Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(AnglesVec(iA) / 180.d0 * Pi) ) 
        
          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

            if (Input%PlotPES_OnlyTriatFlg) then 
              V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * 1.d0/B_To_Ang, Zero*Rp )
            else
              V = Collision%PESsContainer(iPES)%PES%Potential( Rp * 1.d0/B_To_Ang, Zero*Rp )
            end if
            !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
              write(Unit,'(es15.6E3,1(A,es15.6E3))') Rp(1), ',', (V - VRef) * VConverter
            !end if
            
            if (Input%NPESs > 1) then
              SumVec(i)    = SumVec(i)    +   (V - VRef) * VConverter
              SumSqrVec(i) = SumSqrVec(i) + ( (V - VRef) * VConverter )**2
            end if
            
          elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
            
            call Collision%PESsContainer(iPES)%PES%Compute( Rp * 1.d0/B_To_Ang, Zero*Rp, V, dVdR, dVdQ )                                          
            write(Unit,'(es15.6E3,3(A,es15.6E3))') Rp(1), ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
          
          end if  
        
        end do
        
      close(Unit)    
        
    end do
    
    
    if (Input%NPESs > 1) then
    
      MeanVec = (SumVec    / Input%NPESs)
      SDVec   = (SumSqrVec / Input%NPESs - MeanVec**2)**(5.d-1)
    
      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/FixedR1_Stats_' // trim(adjustl(AngleChar)) // '.csv' 
      open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
        write(Unit,'(A)') 'R1, EMean, ESD, EMinus, EPlus'
            
        do i = 1,NGrid
          Rp(1) = ( RStart + hGrid * (i-1) ) 
          Rp(3) = FixedValue       
          Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(AnglesVec(iA) / 180.d0 * Pi) ) 
        
          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

            !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
              write(Unit,'(es15.6E3,4(A,es15.6E3))') Rp(1), ',', MeanVec(i), ',', SDVec(i), ',', MeanVec(i) - 3.d0*SDVec(i), ',', MeanVec(i) + 3.d0*SDVec(i)
            !end if
          
          end if  
        
        end do
        
      close(Unit)    
    
    end if
    
  end do
  
  

  FixedValue1 = 1.09768
  FixedValue2 = 1.09768
  NGrid       = 301
  AlphaStart  = 35.d0
  AlphaEnd    = 175.d0
  hGrid      = (AlphaEnd - AlphaStart) / (NGrid-1)
  if (Input%NPESs > 1) then
    deallocate(SumVec)
    deallocate(SumSqrVec)
    deallocate(MeanVec)
    deallocate(SDVec)
    
    allocate(SumVec(NGrid))
    allocate(SumSqrVec(NGrid))
    allocate(MeanVec(NGrid))
    allocate(SDVec(NGrid))
  end if

      
  do iPES = 1,Input%NPESs  
    write(iPESChar,'(I4)') iPES
  
    if (Input%NPESs > 1) then
      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/FixedR1R2_1.csv.' // trim(adjustl(iPESChar))
    else
      FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/FixedR1R2_1.csv'
    end if
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Angle, E'
          
      do i = 1,NGrid
        
        Angle = ( AlphaStart + hGrid * (i-1) ) 
        Rp(1) = FixedValue1
        Rp(3) = FixedValue2
        Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Angle / 180.d0 * Pi) ) 
      
        if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

          if (Input%PlotPES_OnlyTriatFlg) then 
            V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * 1.d0/B_To_Ang, Zero*Rp )
          else
            V = Collision%PESsContainer(iPES)%PES%Potential( Rp * 1.d0/B_To_Ang, Zero*Rp )
          end if
          !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
            write(Unit,'(es15.6E3,1(A,es15.6E3))') Angle, ',', (V - VRef) * VConverter
          !end if
          
          if (Input%NPESs > 1) then
            SumVec(i)    = SumVec(i)    +   (V - VRef) * VConverter
            SumSqrVec(i) = SumSqrVec(i) + ( (V - VRef) * VConverter )**2
          end if
          
        elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
          
          call Collision%PESsContainer(iPES)%PES%Compute( Rp * 1.d0/B_To_Ang, Zero*Rp, V, dVdR, dVdQ )                                          
          write(Unit,'(es15.6E3,3(A,es15.6E3))') Angle, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
        
        end if  
      
      end do
      
    close(Unit)    
      
  end do
  
  
  if (Input%NPESs > 1) then
  
    MeanVec = (SumVec    / Input%NPESs)
    SDVec   = (SumSqrVec / Input%NPESs - MeanVec**2)**(5.d-1)
  
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/FixedR1R2_Stats_1.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Angle, EMean, ESD, EMinus, EPlus'
          
      do i = 1,NGrid
         
        Angle = ( AlphaStart + hGrid * (i-1) ) 
        Rp(1) = FixedValue1
        Rp(3) = FixedValue2
        Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Angle / 180.d0 * Pi) )  
         
        if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

          !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
            write(Unit,'(es15.6E3,4(A,es15.6E3))') Angle, ',', MeanVec(i), ',', SDVec(i), ',', MeanVec(i) - 3.d0*SDVec(i), ',', MeanVec(i) + 3.d0*SDVec(i)
          !end if
        
        end if  
      
      end do
      
    close(Unit)    
  
  end if
  
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_EvaluatePoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    R1Min    
  real(rkp)                                                 ::    hMin    
  real(rkp)                                                 ::    R1Max  
  real(rkp)                                                 ::    hMax, h  
  real(rkp)                                                 ::    R1Delta
  real(rkp)                                                 ::    hDelta 
  real(rkp)                                                 ::    AngMin
  real(rkp)                                                 ::    AngMax
  real(rkp)                                                 ::    AngDelta
  real(rkp)                                                 ::    Ang  
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf, VRef
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)    ,dimension(2)                                ::    hGrid
  integer                                                   ::    i, j, iAngle, iR, jR
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName, FileNameE
  integer                                                   ::    Unit, UnitRead, UnitReadE
  integer                                                   ::    Status, StatusRead, StatusReadE
  character(:)                 ,allocatable                 ::    FolderName
  logical                                                   ::    i_Debug_Loc
  
  real(rkp)    ,dimension(12)                               ::    AngA
  real(rkp)    ,dimension(13)                               ::    RA1
  real(rkp)    ,dimension(13)                               ::    RA2
      
  real(rkp)    ,dimension(18)                               ::    AngB
  real(rkp)    ,dimension(12)                               ::    RB1
  real(rkp)    ,dimension(19)                               ::    RB2
  real(rkp)    ,dimension(3)                                ::    RpTemp
  real(rkp)                                                 ::    VTemp
  real(rkp)                                                 ::    EPoints, Temp, RNew
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_EvaluatePoints" )
  !i_Debug_Loc   =     Logger%On()


  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))
    
    
    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*Rp )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*Rp )
    end if
    !VInf = 0.d0
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      Rp(1)   = rVMin_Min
      Rp(3)   = 100.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    end if


    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/PESFromReadPoints.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      !write(Unit,'(A)')        "#                 R1                   R2                   R3                    V"
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E"'
      else
        write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
      end if
      if (Status/=0) call Error( "Error opening file: " // FileName )   
      
      if ((Input%NPESs > 1) .and. (.not. Input%StochPESFlg)) then
        FileName = trim(adjustl(Input%OutputDir)) // '/RSampled.csv.' // trim(adjustl(iPESChar))
      else
        FileName = trim(adjustl(Input%OutputDir)) // '/RSampled.csv'
      end if
      write(*,*) 'Reading Points from the File: ', FileName
      open( File=FileName, NewUnit=UnitRead, status='OLD', iostat=StatusRead ) 
      
      
      FileNameE = trim(adjustl(Input%OutputDir)) // '/ESampled.csv.' // trim(adjustl(iPESChar))
      write(*,*) 'Reading E Points from the File: ', FileNameE
      open( File=FileNameE, NewUnit=UnitReadE, status='OLD', iostat=StatusReadE )
      
        do while (StatusRead==0)
          read(UnitRead,*,iostat=StatusRead) Rp(1), Rp(2), Rp(3)
          if (StatusRead/=0) exit
          
          EPoints = 0.d0
          read(UnitReadE,*,iostat=StatusReadE) EPoints, Temp

          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
    
            !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
              !write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter
              !VRef  = Collision%PESsContainer(iPES)%PES%Pairs(1)%Vd%DiatomicPotential( 100.0 )
              write(Unit,'(es15.6,3(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', EPoints
            !end if
        
          end if  
          
        end do

      close(UnitRead)
      close(UnitReadE)
      
    close(Unit)
    
    
  end do

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_ComputeCuts( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    R1Min    
  real(rkp)                                                 ::    hMin    
  real(rkp)                                                 ::    R1Max  
  real(rkp)                                                 ::    hMax, h  
  real(rkp)                                                 ::    R1Delta
  real(rkp)                                                 ::    hDelta 
  real(rkp)                                                 ::    AngMin
  real(rkp)                                                 ::    AngMax
  real(rkp)                                                 ::    AngDelta
  real(rkp)                                                 ::    Ang  
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    dist
  real(rkp)                                                 ::    V
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf, VRef
  real(rkp)                                                 ::    hGrid
  integer                                                   ::    i, j, iAngle, iR, jR
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName, FileNameE
  integer                                                   ::    Unit, UnitRead, UnitReadE
  integer                                                   ::    Status, StatusRead, StatusReadE
  character(:)                 ,allocatable                 ::    FolderName
  logical                                                   ::    i_Debug_Loc
  
  real(rkp)                                                 ::    RStart
  real(rkp)                                                 ::    REnd
  integer                                                   ::    NPoints, iPoints
  real(rkp)                                                 ::    alpha
    
  real(rkp)    ,dimension(12)                               ::    AngA
  real(rkp)    ,dimension(13)                               ::    RA1
  real(rkp)    ,dimension(13)                               ::    RA2
      
  real(rkp)    ,dimension(18)                               ::    AngB
  real(rkp)    ,dimension(12)                               ::    RB1
  real(rkp)    ,dimension(19)                               ::    RB2
  real(rkp)    ,dimension(3)                                ::    RpTemp
  real(rkp)                                                 ::    VTemp
  real(rkp)                                                 ::    EPoints, Temp, RNew
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_ComputeCuts" )
  !i_Debug_Loc   =     Logger%On()

  
  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))
    
    
    RpInf = 1000.d0
    if (Input%PlotPES_OnlyTriatFlg) then 
      VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
    else
      VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
    end if
    !VInf = 0.d0
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      Rp(1)   = rVMin_Min
      Rp(3)   = 100.0
      VRef    = 1.d10
      do while(Rp(1) < rVMin_Max)
        Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
        if (Input%PlotPES_OnlyTriatFlg) then 
          VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
        else
          VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
        end if
        Rp(1) = Rp(1) + 1.d-4
      end do
      VRef  = VRef !- VInf
    end if
    
    
    RStart  = 1.5d0
    REnd    = 10.d0
    NPoints = 1000
    hGrid   = (REnd - RStart) / (NPoints-1)
  
    
    alpha = 110.0
    Rp(1) = 2.26767
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/Cut1.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Variables = "R", "E"'
      if (Status/=0) call Error( "Error opening file: " // FileName )   
      
      Rp(3) = RStart
      do iPoints=1,NPoints
        
        Rp(2) = dsqrt(  Rp(1)**2 + Rp(3)**2 - 2.0d0 * Rp(1) * Rp(3) * dcos( alpha/180.0d0*pi ) )
        V  = Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp )
        write(Unit,*) Rp(3), ',', (V - VRef) * VConverter
        Rp(3) = Rp(3) + hGrid
      
      end do

    close(Unit)
    
    
    
    alpha = 170.0
    Rp(1) = 2.26767
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/Cut2.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Variables = "R", "E"'
      if (Status/=0) call Error( "Error opening file: " // FileName )   
      
      Rp(3) = RStart
      do iPoints=1,NPoints
        
        Rp(2) = dsqrt(  Rp(1)**2 + Rp(3)**2 - 2.0d0 * Rp(1) * Rp(3) * dcos( alpha/180.0d0*pi ) )
        V  = Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp )
        write(Unit,*) Rp(3), ',', (V - VRef) * VConverter
        Rp(3) = Rp(3) + hGrid
        
      end do

    close(Unit)
    
    
    
    alpha = 60.0
    Rp(1) = 2.64562
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/Cut3.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Variables = "R", "E"'
      if (Status/=0) call Error( "Error opening file: " // FileName )   
      
      Rp(3) = RStart
      do iPoints=1,NPoints
        
        Rp(2) = dsqrt(  Rp(1)**2 + Rp(3)**2 - 2.0d0 * Rp(1) * Rp(3) * dcos( alpha/180.0d0*pi ) )
        V  = Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp )
        write(Unit,*) Rp(3), ',', (V - VRef) * VConverter
        Rp(3) = Rp(3) + hGrid
        
      end do

    close(Unit)
    
    
    alpha = 116.75
    Rp(1) = 2.28203327
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName))  // '/Cut4.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Variables = "R", "E"'
      if (Status/=0) call Error( "Error opening file: " // FileName )   
      
      Rp(3) = RStart
      do iPoints=1,NPoints
        
        Rp(2) = dsqrt(  Rp(1)**2 + Rp(3)**2 - 2.0d0 * Rp(1) * Rp(3) * dcos( alpha/180.0d0*pi ) )
        V  = Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp )
        write(Unit,*) Rp(3), ',', (V - VRef) * VConverter
        Rp(3) = Rp(3) + hGrid
        
      end do

    close(Unit)
    
    
  end do

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_Rot3rd( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
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
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(2)                                ::    hGrid
  real(rkp)                                                 ::    Temp
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    t2, t1
  real(rkp)                                                 ::    x, y
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_Rot3rd" )
  !i_Debug_Loc   =     Logger%On()
  

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
  
  hGrid(1) = (Two * Input%MaxGridPlot(1))                  / (Input%NGridPlot(1) - 1)
  hGrid(2) = (Input%MaxGridPlot(2) - Input%MinGridPlot(2)) / (Input%NGridPlot(2) - 1)  
   

  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))


    if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
      RpInf = 1000.d0
      if (Input%PlotPES_OnlyTriatFlg) then 
        VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
      else
        VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
      end if
    end if
    !VInf  = Zero
    write(*,*) 'VInf = ', VInf*VConverter, ' eV'
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
        Rp(1)   = rVMin_Min
        Rp(3)   = 1000.0
        VRef    = 1.d10
        do while(Rp(1) < rVMin_Max)
          Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
          if (Input%PlotPES_OnlyTriatFlg) then 
            VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
          else
            VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
          end if
          Rp(1) = Rp(1) + 1.d-4
        end do
        VRef  = VRef !- VInf
      end if
    end if
    write(*,*) 'VRef = ', VRef*VConverter, ' eV'


    
    call cpu_time ( t1 )
      
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESFrom3dRot.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'Variables = "x", "y", "r1", "r2", "r3", "E"'
      else
        write(Unit,'(A)') 'Variables = "x", "y", "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
      end if
      if (Status/=0) call Error( "Error opening file: " // FileName )
      
      
      iTot = 1
      do i = 1,Input%NGridPlot(1)
        
        x = ( -Input%MaxGridPlot(1) + hGrid(1) * (i-1) ) 


        do j = 1,Input%NGridPlot(2)
          
          y = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) ) 
          write(*,*) 'x = ', x
          write(*,*) 'y = ', y  
          Rp(1) = Input%MinGridPlot(1) 
          Rp(2) = sqrt( (x - Rp(1)/Two)**2 + (y)**2 )
          Rp(3) = sqrt( (x + Rp(1)/Two)**2 + (y)**2 )
          
          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

            if (Input%PlotPES_OnlyTriatFlg) then 
              !V = Collision%PESsContainer(iPES)%PES%DiatPotential( Rp * RConverter )
              V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
            else
              V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
            end if
            !Temp = (V - VRef) * VConverter / abs((V - VRef) * VConverter)
            if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
              write(Unit,'(es17.6E3,5(A,es17.6E3))') x, ',', y, ',', Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter!Temp*max(abs((V - VRef) * VConverter), 1.d-90 )
            end if
            
          elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
            
            call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )     
            if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                           
              write(Unit,'(es15.6,8(A,es15.6))') x, ',', y, ',', Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
            end if
            
          end if  
        
          iTot = iTot + 1
            
        end do
        
      end do
      
    close(Unit)   
    
    call cpu_time ( t2 )
    write(*,*) 'Time for PES Calculations = ', t2-t1
      

  end do
  
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb3_PlotPES_IsoTri( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb3_PlotPES_Type)                   ,intent(out)    ::    This
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
  real(rkp)    ,dimension(NPairs)                           ::    Rp
  real(rkp)    ,dimension(NAtoms*3)                         ::    Qp
  real(rkp)    ,dimension(NPairs)                           ::    RpInf
  real(rkp)    ,dimension(NPairs)                           ::    dVdR
  real(rkp)    ,dimension(NAtoms*3)                         ::    dVdQ
  real(rkp)                                                 ::    Angle
  real(rkp)                                                 ::    VInf
  real(rkp)                                                 ::    VRef
  real(rkp)    ,dimension(2)                                ::    hGrid
  real(rkp)                                                 ::    Temp
  integer                                                   ::    i, j
  integer                                                   ::    iTot
  integer                                                   ::    iA
  integer                                                   ::    iPES
  character(4)                                              ::    iPESChar
  character(:)                 ,allocatable                 ::    FileName  
  integer                                                   ::    Unit
  integer                                                   ::    Status
  character(:)                 ,allocatable                 ::    FolderName
  real(rkp)                                                 ::    t2, t1
  real(rkp)                                                 ::    x, y
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb3_PlotPES_IsoTri" )
  !i_Debug_Loc   =     Logger%On()
  

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
  
  hGrid(1) = (Two * Input%MaxGridPlot(1))                  / (Input%NGridPlot(1) - 1)
  hGrid(2) = (Input%MaxGridPlot(2) - Input%MinGridPlot(2)) / (Input%NGridPlot(2) - 1)  
   

  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))


    if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
      RpInf = 1000.d0
      if (Input%PlotPES_OnlyTriatFlg) then 
        VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, Zero*RpInf )
      else
        VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, Zero*RpInf )
      end if
    end if
    !VInf  = Zero
    write(*,*) 'VInf = ', VInf*VConverter, ' eV'
    
    
    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
        Rp(1)   = rVMin_Min
        Rp(3)   = 1000.0
        VRef    = 1.d10
        do while(Rp(1) < rVMin_Max)
          Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
          if (Input%PlotPES_OnlyTriatFlg) then 
            VRef  = min(Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Zero*Rp ), VRef)
          else
            VRef  = min(Collision%PESsContainer(iPES)%PES%Potential( Rp, Zero*Rp ), VRef)
          end if
          Rp(1) = Rp(1) + 1.d-4
        end do
        VRef  = VRef !- VInf
      end if
    end if
    write(*,*) 'VRef = ', VRef*VConverter, ' eV'


    
    call cpu_time ( t1 )
      
    FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESIsoTri.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
      if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
        write(Unit,'(A)') 'Variables = "x", "y", "r1", "r2", "r3", "E"'
      else
        write(Unit,'(A)') 'Variables = "x", "y", "r1", "r2", "r3", "E", "dEdR1", "dEdR2", "dEdR3"'
      end if
      if (Status/=0) call Error( "Error opening file: " // FileName )
      
      
      iTot = 1
      do i = 1,Input%NGridPlot(1)
        
        x = ( -Input%MaxGridPlot(1) + hGrid(1) * (i-1) ) 

        do j = 1,Input%NGridPlot(2)
          
          y = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) ) 

          Rp(1) = x
          Rp(2) = sqrt((x*0.5d0)**2+y**2)
          Rp(3) = sqrt((x*0.5d0)**2+y**2)
          
          if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

            if (Input%PlotPES_OnlyTriatFlg) then 
              V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Zero*Rp )
            else
              V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Zero*Rp )
            end if
            !Temp = (V - VRef) * VConverter / abs((V - VRef) * VConverter)
            if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
              write(Unit,'(es17.6E3,5(A,es17.6E3))') x, ',', y, ',', Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter!Temp*max(abs((V - VRef) * VConverter), 1.d-90 )
            end if
            
          elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
            
            call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )     
            if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                           
              write(Unit,'(es15.6,8(A,es15.6))') x, ',', y, ',', Rp(1), ',', Rp(2), ',', Rp(3), ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
            end if
            
          end if  
        
          iTot = iTot + 1
            
        end do
        
      end do
      
    close(Unit)   
    
    call cpu_time ( t2 )
    write(*,*) 'Time for PES Calculations = ', t2-t1
     

  end do
  
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module