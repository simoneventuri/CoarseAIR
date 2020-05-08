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

Module Nb4_PlotPES_Class
  
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
  public  ::    Nb4_PlotPES_Type
  
  Type    ,extends(PlotPES_Type)                            ::    Nb4_PlotPES_Type
  contains
    private
    procedure              ,public                          ::    Initialize             =>    Nb4_PlotPES_Initialize
    procedure              ,public                          ::    IsoTri                 =>    Nb4_PlotPES_IsoTri
    procedure              ,public                          ::    Rot3rd                 =>    Nb4_PlotPES_Rot3rd
    procedure              ,public                          ::    ComputeCuts            =>    Nb4_PlotPES_ComputeCuts
    procedure              ,public                          ::    EvaluatePoints         =>    Nb4_PlotPES_EvaluatePoints
    procedure              ,public                          ::    PlotsVargasPaper       =>    Nb4_PlotPES_PlotsVargasPaper
    procedure              ,public                          ::    StochPESStats          =>    Nb4_PlotPES_StochPESStats
    procedure              ,public                          ::    GridForScatter         =>    Nb4_PlotPES_GridForScatter
    procedure              ,public                          ::    ReadPoints             =>    Nb4_PlotPES_ReadPoints
    procedure              ,public                          ::    TripleGrid             =>    Nb4_PlotPES_TripleGrid
    procedure              ,public                          ::    DoubleGrid             =>    Nb4_PlotPES_DoubleGrid
    procedure              ,public                          ::    GridForStochPES        =>    Nb4_PlotPES_GridForStochPES
    procedure              ,public                          ::    Grid                   =>    Nb4_PlotPES_Grid
  End Type

  integer   ,parameter    ::    NSpace         = 3
  logical   ,parameter    ::    Formatted      = .True.
  integer                 ::    iSpeTar        = 1
  integer                 ::    iSpePro        = 2
  real(rkp)               ::    MostProbEr                            
  real(rkp)               ::    rVMin_Min      = 2.0d0
  real(rkp)               ::    rVMin_Max      = 2.5d0
  
  real(rkp)    ,save      ::    RConverter     = One
  real(rkp)    ,save      ::    VConverter     = One
  real(rkp)    ,save      ::    dVConverter    = One

  logical   ,parameter    ::    i_Debug_Global = .False.
  
  contains

!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_Initialize( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_Initialize" )
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

      write(*,*) RConverter
      write(*,*) dVConverter
      write(*,*) dVConverter

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_Grid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    ang
  real(rkp)                                                 ::    beta
  real(rkp)                                                 ::    Theta2, Theta4
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
  real(rkp)                                                 ::    RInf
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
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_Grid" )
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
  
  hGrid = (Input%MaxGridPlot - Input%MinGridPlot) / (Input%NGridPlot - 1)  
   

  do iPES = 1,Input%NPESs
    write(iPESChar,'(I4)') iPES
    FolderName = 'PES_' // trim(adjustl(iPESChar))
    call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)))

    RInf  = 1000.d0
    RpInf = RInf

    if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
      if (Input%PlotPES_OnlyTriatFlg) then 
        if (Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
          call R_to_X(RpInf, QpInf)
        else
          QpInf = Zero
        end if
        VInf  = Collision%PESsContainer(iPES)%PES%TriatPotential( RpInf, QpInf )
      else
        VInf  = Collision%PESsContainer(iPES)%PES%Potential( RpInf, QpInf )
      end if
    end if
    !VInf  = Zero
    write(*,*) 'VInf = ', VInf*VConverter, ' eV'
    
    write(*,*) RConverter
    write(*,*) dVConverter
    write(*,*) dVConverter

    if (Input%PESZeroRefIntFlg == 0) then
      VRef  = 0.0d0
    else
      if (.not. Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
        Rp      = RpInf
        Rp(1)   = rVMin_Min
        VRef    = 1.d10
        do while(Rp(1) < rVMin_Max)
          Rp(2) = dsqrt(  Rp(1)**Two + Rp(3)**Two - Two * Rp(1) * Rp(3) * dcos( 120.d0 / 180.d0 * pi ) )
          if (Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
            call R_to_X(Rp, Qp, Theta=120.d0)
          else
            Qp = Zero
          end if
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


      if (Collision%PESsContainer(iPES)%PES%CartCoordFlg) then


        do iA = 1,size(Input%AnglesPlot,1)
        
          call cpu_time ( t1 )


          FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESFromGrid.csv.' // trim(adjustl(Input%AnglesPlotChar(iA)))
          open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
              !write(Unit,'(A)') 'Variables = "A1_x", "A1_y", "A1_z", "A2_x", "A2_y", "A2_z", "A3_x", "A3_y", "A3_z", "A4_x", "A4_y", "A4_z", "E"'
              write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "r4", "r5", "r6", "E"'
            else
              !write(Unit,'(A)') 'Variables = "A1_x", "A1_y", "A1_z", "A2_x", "A2_y", "A2_z", "A3_x", "A3_y", "A3_z", "A4_x", "A4_y", "A4_z", "E", "dA1_x", "dA1_y", "dA1_z", "dA2_x", "dA2_y", "dA2_z", "dA3_x", "dA3_y", "dA3_z", "dA4_x", "dA4_y", "dA4_z"'
              write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "r4", "r5", "r6", "E", "dEdR1", "dEdR2", "dEdR3", "dEdR4", "dEdR5", "dEdR6"'
            end if
            if (Status/=0) call Error( "Error opening file: " // FileName )


            iTot = 1
            Rp     = RpInf
            do i = 1,Input%NGridPlot(1)
              Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) * RConverter

              do j = 1,Input%NGridPlot(2)
                Rp(2)  = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) ) * RConverter     
                
                Theta4 = Input%AnglesPlot(iA) / 180.d0 * Pi
                Rp(4)  = sqrt( Rp(1)**2 + Rp(2)**2 - 2.d0 * Rp(1) * Rp(2) * dcos(Theta4) ) * RConverter 

                Qp =[              0.0,               0.0, 0.0, &
                                 Rp(1),               0.0, 0.0, &
                     Rp(2)*cos(Theta4), Rp(2)*sin(Theta4), 0.0, &
                                   0.0,              RInf, 0.0] 


                if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

                  if (Input%PlotPES_OnlyTriatFlg) then 
                    !V = Collision%PESsContainer(iPES)%PES%DiatPotential( Rp * RConverter )
                    V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp, Qp )
                  else
                    V = Collision%PESsContainer(iPES)%PES%Potential( Rp, Qp )
                  end if
                    
                  call X_to_R(Qp, Rp)

                  !Temp = (V - VRef) * VConverter / abs((V - VRef) * VConverter)
                  if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                    ! write(Unit,'(es17.6E3,12(A,es17.6E3))')  Qp(1), ',',  Qp(2), ',',  Qp(3), ',', &
                    !                                          Qp(4), ',',  Qp(5), ',',  Qp(6), ',', &
                    !                                          Qp(7), ',',  Qp(8), ',',  Qp(9), ',', &
                    !                                         Qp(10), ',', Qp(11), ',', Qp(12), ',', &
                    !                                         (V - VRef) * VConverter
                    write(Unit,'(es17.6E3,6(A,es17.6E3))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Rp(4), ',', Rp(5), ',', Rp(6), ',', &
                                                           (V - VRef) * VConverter!Temp*max(abs((V - VRef) * VConverter), 1.d-90 )
                  end if
                  
                elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
                  
                  call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )  

                  call X_to_R(Qp, Rp)
                  dVdQ = dVdQ * dVConverter
                  
                  if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                           
                    ! write(Unit,'(es15.6,24(A,es15.6))') Qp(1), ',',  Qp(2), ',',  Qp(3), ',',       &
                    !                                     Qp(4), ',',  Qp(5), ',',  Qp(6), ',',       &
                    !                                     Qp(7), ',',  Qp(8), ',',  Qp(9), ',',       &
                    !                                    Qp(10), ',', Qp(11), ',', Qp(12), ',',       &
                    !                                    (V - VRef) * VConverter, ',',                &
                    !                                     dVdQ(1), ',',  dVdQ(2), ',',  dVdQ(3), ',', &
                    !                                     dVdQ(4), ',',  dVdQ(5), ',',  dVdQ(6), ',', &
                    !                                     dVdQ(7), ',',  dVdQ(8), ',',  dVdQ(9), ',', &
                    !                                    dVdQ(10), ',', dVdQ(11), ',', dVdQ(12)
                    write(Unit,'(es15.6,12(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Rp(4), ',', Rp(5), ',', Rp(6), ',', &
                                                        (V - VRef) * VConverter, ',',                                           &
                                                        (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter, ',', (dVdR(4)) * dVConverter, ',', (dVdR(5)) * dVConverter, ',', (dVdR(6)) * dVConverter 
                  end if
                  
                end if  
              
                iTot = iTot + 1
                  
              end do
              
            end do
            
          close(Unit)   
          
          call cpu_time ( t2 )
          write(*,*) 'Time for PES Calculations = ', t2-t1
          
        end do      


      else


        do iA = 1,size(Input%AnglesPlot,1)
        
          call cpu_time ( t1 )
            
          FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESFromGrid.csv.' // trim(adjustl(Input%AnglesPlotChar(iA)))
          open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
            if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
              write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "r4", "r5", "r6", "E"'
            else
              write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "r4", "r5", "r6", "E", "dEdR1", "dEdR2", "dEdR3", "dEdR4", "dEdR5", "dEdR6"'
            end if
            if (Status/=0) call Error( "Error opening file: " // FileName )
            
            iTot = 1
            Rp   = RpInf
            do i = 1,Input%NGridPlot(1)
              Rp(1) = ( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 

              do j = 1,Input%NGridPlot(2)
                Rp(4)  = ( Input%MinGridPlot(2) + hGrid(2) * (j-1) )     
                
                Theta2 = Input%AnglesPlot(iA) / 180.d0 * Pi
                Rp(2)  = sqrt( Rp(1)**2 + Rp(4)**2 - 2.d0 * Rp(1) * Rp(4) * dcos(Theta2) ) 

                if (Collision%PESsContainer(iPES)%PES%CartCoordFlg) then
                  call R_to_X(Rp, Qp, Theta=Theta2)
                else
                  Qp = Zero
                end if
                if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

                  if (Input%PlotPES_OnlyTriatFlg) then 
                    !V = Collision%PESsContainer(iPES)%PES%DiatPotential( Rp * RConverter )
                    V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Qp )
                  else
                    V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Qp )
                  end if
                  !Temp = (V - VRef) * VConverter / abs((V - VRef) * VConverter)
                  if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
                    write(Unit,'(es17.6E3,6(A,es17.6E3))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Rp(4), ',', Rp(5), ',', Rp(6), ',', &
                                                           (V - VRef) * VConverter!Temp*max(abs((V - VRef) * VConverter), 1.d-90 )
                  end if
                  
                elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
                  
                  call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Zero*Rp, V, dVdR, dVdQ )     
                  if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                           
                    write(Unit,'(es15.6,12(A,es15.6))')    Rp(1), ',', Rp(2), ',', Rp(3), ',', Rp(4), ',', Rp(5), ',', Rp(6), ',', &
                                                           (V - VRef) * VConverter, ',',                                      &
                                                           (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter, ',', (dVdR(4)) * dVConverter, ',', (dVdR(5)) * dVConverter, ',', (dVdR(6)) * dVConverter 
                  end if
                  
                end if  
              
                iTot = iTot + 1
                  
              end do
              
            end do
            
          close(Unit)   
          
          call cpu_time ( t2 )
          write(*,*) 'Time for PES Calculations = ', t2-t1
          
        end do   


      end if   
        

    end if
    ! elseif (trim(adjustl(Input%yAxisVar)) .eq. 'Angle') then  
    
    !   FileName = trim(adjustl(Input%OutputDir)) // '/PlotPES/' // trim(adjustl(FolderName)) // '/PESDistVsAngles.csv'
    !   open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    !     if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then
    !       write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "Angle", "E"'
    !     else
    !       write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "Angle", "E", "dEdR1", "dEdR2", "dEdR3"'
    !     end if
    !     if (Status/=0) call Error( "Error opening file: " // FileName )
        
    !     iTot = 1
    !     do i = 1,Input%NGridPlot(1)!radiuses (constrained to be equal)
    !       !this case r1=r2 and angle is varying

    !       Rp(1) = 1.4d0!( Input%MinGridPlot(1) + hGrid(1) * (i-1) ) 
    !       Rp(3) = 1.7d0!Rp(1)
        
    !       do j = 1,Input%NGridPlot(2)!angle
          
    !         Angle = Input%MinGridPlot(2) + hGrid(2) * (j-1)
    !         Rp(2) = sqrt( Rp(1)**2 + Rp(3)**2 - 2.d0 * Rp(1) * Rp(3) * dcos(Angle/ 180.d0 * Pi) ) 
              
    !         call R_to_X(Rp, Qp, Theta=Angle)
    !         if (trim(adjustl(Input%POTorFR)) .eq. 'Potential') then

    !           if (Input%PlotPES_OnlyTriatFlg) then 
    !             V = Collision%PESsContainer(iPES)%PES%TriatPotential( Rp * RConverter, Qp )
    !           else
    !             V = Collision%PESsContainer(iPES)%PES%Potential( Rp * RConverter, Qp )
    !           end if
    !           !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then 
    !             write(Unit,'(es15.6,4(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Angle, ',', (V - VRef) * VConverter
    !           !end if
              
    !         elseif (trim(adjustl(Input%POTorFR)) .eq. 'Force') then           
              
    !           call Collision%PESsContainer(iPES)%PES%Compute( Rp * RConverter, Qp, V, dVdR, dVdQ )    
    !           !if ( (V - Vinf)*VConverter <= Input%EnergyCutOff ) then                                       
    !             write(Unit,'(es15.6,7(A,es15.6))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Angle, ',', (V - VRef) * VConverter, ',', (dVdR(1)) * dVConverter, ',', (dVdR(2)) * dVConverter, ',', (dVdR(3)) * dVConverter 
    !           !end if
            
    !         end if  
           
    !         iTot = iTot + 1
    !       end do
        
    !     end do
      
    !   close(Unit)   
      
        
    ! end if ! Distance/Angle
  
  end do

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_GridForStochPES( This, Input, Collision, NPairs, NAtoms,  i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_GridForStochPES" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_DoubleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_DoubleGrid" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_TripleGrid( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_TripleGrid" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_ReadPoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_ReadPoints" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_GridForScatter( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_GridForScatter" )
  !i_Debug_Loc   =     Logger%On()
        
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_StochPESStats( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_StochPESStats" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_PlotsVargasPaper( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_PlotsVargasPaper" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_EvaluatePoints( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_EvaluatePoints" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_ComputeCuts( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_ComputeCuts" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_Rot3rd( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_Rot3rd" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Nb4_PlotPES_IsoTri( This, Input, Collision, NPairs, NAtoms, i_Debug )

  class(Nb4_PlotPES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  integer                                   ,intent(in)     ::    NPairs
  integer                                   ,intent(in)     ::    NAtoms
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Nb4_PlotPES_IsoTri" )
  !i_Debug_Loc   =     Logger%On()
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module