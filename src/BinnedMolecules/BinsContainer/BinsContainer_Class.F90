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

Module BinsContainer_Class

#include "../../qct.inc"

  use Parameters_Module       ,only:  rkp, Zero, Half, Two, UKb, Ue, Hartree_To_eV
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error
  use Degeneracy_Class        ,only:  Degeneracy_Type
  use BinsContainerT_Class    ,only:  BinsContainerT_Type
  use Bin_Class               ,only:  Bin_Type

  implicit none

  private
  public  ::    BinsContainer_Type

  Type    ,abstract                                       ::    BinsContainer_Type
    logical                                               ::    Initialized         !< Indicator whether the object is initialized
    integer                                               ::    NBins
    integer        ,dimension(-1:1000,-1:1000)            ::    qns_to_Bin = -1          !< Matrix (vqn,jqn)  -> Bin Nb
    integer        ,dimension(:)             ,allocatable ::    Lvl_to_Bin               !< Vector (Level Nb) -> Bin Nb
    real(rkp)      ,dimension(:,:)           ,allocatable ::    OverallRate
    real(rkp)      ,dimension(:,:)           ,allocatable ::    OverallRateSigma2
    class(Degeneracy_Type)                   ,allocatable ::    Degeneracy
    type(Bin_Type) ,dimension(:)             ,allocatable ::    Bin                      !< Array of Molecule Bins 
    integer                                               ::    To_Pair
    integer                                               ::    To_Molecule
    real(rkp)                                             ::    Q
    real(rkp)                                             ::    QInit
    type(BinsContainerT_Type) ,dimension(:)  ,allocatable ::    Ttra
    character(:)                             ,allocatable ::    PathToBinMolFldr
  contains
    private
    procedure              ,public                        ::    Initialize            =>    Initialize_BinsContainer
    procedure              ,public                        ::    WriteBinsData         =>    WriteBinsData_BinsContainer
    procedure              ,public                        ::    ComputePartFunEnergy  =>    ComputePartFunEnergy_BinsContainer
    procedure              ,public                        ::    ReadPartFunEnergy     =>    ReadPartFunEnergy_BinsContainer
    procedure              ,public                        ::    ReadLevelToBin        =>    ReadLevelToBin_BinsContainer
  End Type
  
  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine Initialize_BinsContainer( This, Input, LevelsContainer, iMol, WriteFlg, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    WriteFlg
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  logical                                                     ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()
  
  
  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteBinsData_BinsContainer( This, Input, LevelsContainer, iMol, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)     ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(in)     ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
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


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteBinsData_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()  


  ! Write the Mapping Quantum Numbers -> Bin
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping Quantum Numbers -> Bin " )
  FileName = This%PathToBinMolFldr // '/qns_to_Bin.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )  
  do ijqn = 0,LevelsContainer%maxjqn
    write(Unit,*) This%qns_to_Bin(:,ijqn)
  end do   
  close(Unit) 
  
  
  ! Write the Mapping Bins -> First Quantum Number
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/qnsFirst.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )  
  write(Unit,*) ("#    1st vqn   1st jqn   N Levels % Levels ")
  do iBins = 1,Input%NBins(iMol)
    write(Unit,1) This%Bin(iBins)%vqnFirst, This%Bin(iBins)%jqnFirst, This%Bin(iBins)%NLevels, real(This%Bin(iBins)%NLevels) * 1.e2_rkp / real(LevelsContainer%NStates)
  end do   
  close(Unit) 


  ! Write the Mapping qns and Energy [eV]-> Bin
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping qns and Energy [eV]-> Bin" )  
  FileName = This%PathToBinMolFldr // '/qnsEnBin.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )  
  write(Unit,*) ("#  Level   vqn    jqn        E [Eh]      Degeneracy  Bin ")
  do iLevels = 1,LevelsContainer%NStates 
    write(Unit,2) iLevels, int(LevelsContainer%States(iLevels)%vqn),  int(LevelsContainer%States(iLevels)%jqn), LevelsContainer%States(iLevels)%einteV, LevelsContainer%States(iLevels)%g, This%Lvl_to_Bin(iLevels)
  end do     
  close(Unit) 
  
  
  1 format (3X 3I10, es10.2)
  2 format (1X, 3I7, 2es15.7, I10)
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ComputePartFunEnergy_BinsContainer( This, Input, LevelsContainer, iMol, i_Debug )

  use Input_Class              ,only:  Input_Type
  use Degeneracy_Factory_Class ,only:  Degeneracy_Factory_Type
  use LevelsContainer_Class    ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  Type(Degeneracy_Factory_Type)                               ::    Degeneracy_Factory
  integer                                                     ::    iBins
  integer                                                     ::    iLevels
  real(rkp) ,dimension(:) ,allocatable                        ::    StateQ
  real(rkp) ,dimension(:) ,allocatable                        ::    StateQInit
  integer                                                     ::    Status
  integer                                                     ::    Unit
  character(:)            ,allocatable                        ::    FileName
  character(10)                                               ::    T_char
  logical                                                     ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ComputePartFunEnergy_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
      
  ! Computing Level Degeneracies
  call Degeneracy_Factory%Define_Degeneracy( Input, This%Degeneracy, iMol, i_Debug=i_Debug_Loc )   
  do iLevels = 1,LevelsContainer%NStates
    LevelsContainer%States(iLevels)%g = This%Degeneracy%Compute_Degeneracy_State( LevelsContainer%States(iLevels)%jqn )
  end do
  
  ! Computing Bins Effective Energy Minima (the Minima specified by the User could differ from actual Molecule Levels Energies)
  LevelsContainer%States%einteV_scaled = LevelsContainer%States%einteV - LevelsContainer%MineinteV
  This%Bin%MineinteV = 1000.d0
  do iLevels = 1,LevelsContainer%NStates
    if (LevelsContainer%States(iLevels)%einteV < This%Bin(This%Lvl_to_Bin(iLevels))%MineinteV) then 
      This%Bin(This%Lvl_to_Bin(iLevels))%MineinteV = LevelsContainer%States(iLevels)%einteV
    end if
  end do
  
  ! Computing Bins Partition Functions and Energies
  if (.not. allocated(StateQ)) then 
    allocate(StateQ(LevelsContainer%NStates), Stat=Status )
    if (Status/=0) call Error( "Error allocating StateQ" )
    if (i_Debug_Loc)  call Logger%Write( "StateQ, Levels Partition Function Vector, allocated " ) 
  end if
  
  if (.not. allocated(StateQInit)) then 
    allocate(StateQInit(LevelsContainer%NStates), Stat=Status )
    if (Status/=0) call Error( "Error allocating StateQInit" )
    if (i_Debug_Loc)  call Logger%Write( "StateQInit, Levels Partition Function Vector at KONIG Initial Temperature, allocated " ) 
  end if
  
  This%Bin(:)%ToteinteV     = Zero
  This%Bin(:)%Q             = Zero
  This%Bin(:)%ToteinteVInit = Zero
  This%Bin(:)%QInit         = Zero
  
  do iLevels = 1,LevelsContainer%NStates
    StateQInit(iLevels)                              = LevelsContainer%States(iLevels)%g * dexp( - (LevelsContainer%States(iLevels)%einteV_scaled ) * Ue / (UKb * Input%TInit) )   
    This%Bin(This%Lvl_to_Bin(iLevels))%ToteinteVInit = This%Bin(This%Lvl_to_Bin(iLevels))%ToteinteVInit + StateQInit(iLevels) * LevelsContainer%States(iLevels)%einteV_scaled
    This%Bin(This%Lvl_to_Bin(iLevels))%QInit         = This%Bin(This%Lvl_to_Bin(iLevels))%QInit         + StateQInit(iLevels)
  end do
  This%Bin%ToteinteVInit = This%Bin%ToteinteVInit / This%Bin%QInit
   
   
  if (( Input%vInMethod .eq. "Uniform" ) .or. (trim(Input%BSortMethod(iMol)) .eq. "State-Specific" )) then
    if (i_Debug_Loc)  call Logger%Write( "Uniform Method for Bins Internal States Populations Reconstruction" )   
   
    do iLevels = 1,LevelsContainer%NStates
      StateQ(iLevels)                              = LevelsContainer%States(iLevels)%g
      This%Bin(This%Lvl_to_Bin(iLevels))%ToteinteV = This%Bin(This%Lvl_to_Bin(iLevels))%ToteinteV + StateQ(iLevels) * LevelsContainer%States(iLevels)%einteV_scaled
      This%Bin(This%Lvl_to_Bin(iLevels))%Q         = This%Bin(This%Lvl_to_Bin(iLevels))%Q         + StateQ(iLevels)
    end do
    This%Bin%ToteinteV = This%Bin%ToteinteV / This%Bin%Q
         
  elseif ( Input%vInMethod .eq. "Boltzmann" ) then
    if (i_Debug_Loc)  call Logger%Write( "Boltzmann Method for Bins Internal States Populations Reconstruction" )   
   
    do iLevels = 1,LevelsContainer%NStates    
      StateQ(iLevels)     = LevelsContainer%States(iLevels)%g * dexp( - (LevelsContainer%States(iLevels)%einteV_scaled ) * Ue / (UKb * Input%Tint)  )          
      This%Bin(This%Lvl_to_Bin(iLevels))%ToteinteV = This%Bin(This%Lvl_to_Bin(iLevels))%ToteinteV + StateQ(iLevels) * LevelsContainer%States(iLevels)%einteV_scaled
      This%Bin(This%Lvl_to_Bin(iLevels))%Q         = This%Bin(This%Lvl_to_Bin(iLevels))%Q         + StateQ(iLevels)
    end do
    This%Bin%ToteinteV = This%Bin%ToteinteV / This%Bin%Q
 
  else
    if (i_Debug_Loc)  call Logger%Write( "ERROR: Bins Internal States Populations Reconstruction Method NOT SPECIFIED" )   
    write(*,*) ' BinsContainer_Class Error 1 '
    stop ( "ERROR: Bins Internal States Populations Reconstruction Method NOT SPECIFIED" )
  end if
  
  
  ! Computing the Bins Ratios Partition Functions
  This%Q = Zero
  do iBins = 1,Input%NBins(iMol)
    This%Q = This%Q + This%Bin(iBins)%Q 
  end do
  do iBins = 1,Input%NBins(iMol)
    This%Bin(iBins)%QRatio = This%Bin(iBins)%Q / This%Q
  end do
  
  
  ! Computing the Bins Ratios Partition Functions
  This%QInit = Zero
  do iBins = 1,Input%NBins(iMol)
    This%QInit = This%QInit + This%Bin(iBins)%QInit 
  end do
  do iBins = 1,Input%NBins(iMol)
    This%Bin(iBins)%QRatioInit = This%Bin(iBins)%QInit / This%QInit
  end do
  
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%Tint)
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,'(A60)') "#          Bin Ratio           Bin PartF          Bin Energy "
  if (Status/=0) call Error( "Error opening file: " // FileName )                      
    do iBins = 1,Input%NBins(iMol)
      write(Unit,'(3es20.10E3)') This%Bin(iBins)%QRatio, This%Bin(iBins)%Q, This%Bin(iBins)%ToteinteV
    end do   
  close(Unit) 
  
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%TInit)
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.dat'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,'(A60)') "#          Bin Ratio           Bin PartF          Bin Energy "
  if (Status/=0) call Error( "Error opening file: " // FileName )                      
    do iBins = 1,Input%NBins(iMol)
      write(Unit,'(3es20.10E3)') This%Bin(iBins)%QRatioInit, This%Bin(iBins)%QInit, This%Bin(iBins)%ToteinteVInit
    end do   
  close(Unit) 
  

  deallocate(StateQ, Stat=Status )
  if (Status/=0) call Error( "Error deallocating StateQ" )
  if (i_Debug_Loc)  call Logger%Write( "StateQ, Levels Partition Function Vector, dellocated" )  
  
  deallocate(StateQInit, Stat=Status )
  if (Status/=0) call Error( "Error deallocating StateQInit" )
  if (i_Debug_Loc)  call Logger%Write( "StateQInit, Initial Levels Partition Function Vector, dellocated" )  
  
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadPartFunEnergy_BinsContainer( This, Input, LevelsContainer, iMol, i_Debug )

  use Input_Class              ,only:  Input_Type
  use LevelsContainer_Class    ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    iBins
  integer                                                     ::    iLevels
  real(rkp) ,dimension(:) ,allocatable                        ::    StateQ
  real(rkp) ,dimension(:) ,allocatable                        ::    StateQInit
  integer                                                     ::    Status
  integer                                                     ::    Unit
  character(:)            ,allocatable                        ::    FileName
  character(10)                                               ::    T_char

  logical                                                     ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadPartFunEnergy_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
      
  
  This%Bin(:)%ToteinteV     = Zero
  This%Bin(:)%Q             = Zero
  This%Bin(:)%ToteinteVInit = Zero
  This%Bin(:)%QInit         = Zero
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%Tint)
  if (i_Debug_Loc)  call Logger%Write( "Reading the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.dat'
  if (i_Debug_Loc)  call Logger%Write( "Open File ", FileName ) 
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error Opening file: " // FileName )                      
    read(Unit,*)
    do iBins = 1,Input%NBins(iMol)
      read(Unit,'(3es20.10E3)', iostat=Status) This%Bin(iBins)%QRatio, This%Bin(iBins)%Q, This%Bin(iBins)%ToteinteV
      if (Status/=0) call Error( "Error Reading file: " // FileName )  
    end do   
  close(Unit) 
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%TInit)
  if (i_Debug_Loc)  call Logger%Write( "Reading the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.dat'
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error Opening file: " // FileName )  
    read(Unit,*)       
    do iBins = 1,Input%NBins(iMol)
      read(Unit,'(3es20.10E3)', iostat=Status) This%Bin(iBins)%QRatioInit, This%Bin(iBins)%QInit, This%Bin(iBins)%ToteinteVInit
      if (Status/=0) call Error( "Error Reading file: " // FileName )  
    end do   
  close(Unit) 
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadLevelToBin_BinsContainer( This, Input, LevelsContainer, iMol, i_Debug )

  use Input_Class              ,only:  Input_Type
  use Degeneracy_Factory_Class ,only:  Degeneracy_Factory_Type
  use LevelsContainer_Class    ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    iLevels
  integer                                                     ::    iTemp1, iTemp2, iTemp3
  real(rkp)                                                   ::    xTemp1, xTemp2
  character                                                   ::    charTemp
  integer                                                     ::    Status
  integer                                                     ::    Unit
  character(:)            ,allocatable                        ::    FileName
  character(10)                                               ::    T_char

  logical                                                     ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadLevelToBin_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  ! Write the Mapping qns and Energy [eV]-> Bin
  if (i_Debug_Loc)  call Logger%Write( "Reading the Mapping Level -> Bin" )  
  FileName = trim(adjustl(Input%SystemCGPath)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
             trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBinsCG_char(iMol))) // '/qnsEnBin.dat'
  if (i_Debug_Loc)  call Logger%Write( "Opening the File: ", trim(adjustl(FileName)) )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error Opening File: " // FileName )  
  read(Unit,*)
  do iLevels = 1,LevelsContainer%NStates 
    !read(Unit,'(I8, 2I7, 2es15.7, I10)')  iTemp1, iTemp2, iTemp3, xTemp1, xTemp2, LevelsContainer%States(iLevels)%ToBin
    read(Unit,'(I8, 2I7, 2es15.7, I10)', iostat=Status)  iTemp1, iTemp2, iTemp3, xTemp1, xTemp2, LevelsContainer%States(iLevels)%ToBin
    if (Status/=0) call Error( "Error Reading File: " // FileName )
  end do     
  close(Unit) 

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module