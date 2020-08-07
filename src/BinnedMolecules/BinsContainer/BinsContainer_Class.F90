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
    character(:)                             ,allocatable ::    NBins_char
    integer        ,dimension(-1:1000,-1:1000)            ::    qns_to_Bin = -1          !< Matrix (vqn,jqn)  -> Bin Nb
    real(rkp)      ,dimension(:,:)           ,allocatable ::    OverallRate
    real(rkp)      ,dimension(:,:)           ,allocatable ::    OverallRateSigma2
    class(Degeneracy_Type)                   ,allocatable ::    Degeneracy
    type(Bin_Type) ,dimension(:)             ,allocatable ::    Bin                      !< Array of Molecule Bins 
    integer                                               ::    To_Pair
    integer                                               ::    To_Molecule
    real(rkp)                                             ::    Q
    real(rkp)                                             ::    QInit
    type(BinsContainerT_Type) ,dimension(:)  ,allocatable ::    Ttra
    character(:)                             ,allocatable ::    PathToBinMolDtbFldr
    character(:)                             ,allocatable ::    PathToBinMolFldr
    character(150)            ,dimension(2)               ::    BinDataFile = ' '
  contains
    private
    procedure              ,public                        ::    MainInitialize        =>    MainInitialize_BinsContainer
    procedure              ,public                        ::    Initialize            =>    Initialize_BinsContainer
    procedure              ,public                        ::    Compute_Level_To_Bin  =>    Compute_Level_To_Bin_BinsContainer
    procedure              ,public                        ::    WriteNLevels          =>    WriteNLevels_BinsContainer
    procedure              ,public                        ::    WriteLevelsInBin      =>    WriteLevelsInBin_BinsContainer
    procedure              ,public                        ::    WriteQNsInBin         =>    WriteQNsInBin_BinsContainer
    procedure              ,public                        ::    WriteQNs_To_Bin       =>    WriteQNs_To_Bin_BinsContainer
    procedure              ,public                        ::    WriteQNsEnBin         =>    WriteQNsEnBin_BinsContainer
    procedure              ,public                        ::    WriteQNsFirst         =>    WriteQNsFirst_BinsContainer
    procedure              ,public                        ::    WriteBinsData         =>    WriteBinsData_BinsContainer
    procedure              ,public                        ::    ComputePartFunEnergy  =>    ComputePartFunEnergy_BinsContainer
    procedure              ,public                        ::    ReadPartFunEnergy     =>    ReadPartFunEnergy_BinsContainer
    procedure              ,public                        ::    ReadLevelToBin        =>    ReadLevelToBin_BinsContainer
  End Type
  
  logical   ,parameter                                      ::    i_Debug_Global = .False.

  contains


!________________________________________________________________________________________________________________________________!
Subroutine MainInitialize_BinsContainer( This, Input, iMol, DtbFldrPath, FldrPath, i_Debug )

  use Input_Class             ,only:  Input_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  integer                                     ,intent(in)     ::    iMol
  character(*)                                ,intent(in)     ::    DtbFldrPath
  character(*)                                ,intent(in)     ::    FldrPath
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  integer                                                     ::    Status
  logical                                                     ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "MainInitialize_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()
  

  This%NBins = Input%NBins(iMol)
  allocate( This%NBins_char, source=trim(adjustl( Input%NBins_char(iMol) )) )


  ! Allocate Bins
  if (.not. allocated(This%Bin)) then 
    allocate(This%Bin(This%NBins), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%Bin" )
    if (i_Debug_Loc)  call Logger%Write( "Allocated ", This%NBins, " Bins for the Molecule Nb", iMol )
  end if


  allocate( This%PathToBinMolDtbFldr, source = trim(adjustl( trim(adjustl(DtbFldrPath)) // 'Bins_' // This%NBins_char // '/' )) )
  allocate( This%PathToBinMolFldr,    source = trim(adjustl( trim(adjustl(FldrPath))    // 'Bins_' // This%NBins_char // '/' )) )
  if (i_Debug_Loc)  call Logger%Write( "Creating Folder ", This%PathToBinMolFldr )
  call system('mkdir -p ' // This%PathToBinMolFldr )


  This%BinDataFile(1) = trim(adjustl( This%PathToBinMolDtbFldr // trim(adjustl(Input%BinDataFile(iMol,1))) ))
  if (i_Debug_Loc) call Logger%Write( "I will read the Level-To-File Mapping from: This%BinDataFile(1) = ", This%BinDataFile(1) )

  This%BinDataFile(2) = trim(adjustl( This%PathToBinMolDtbFldr // trim(adjustl(Input%BinDataFile(iMol,2))) ))
  if (i_Debug_Loc) call Logger%Write( "I will read the Level-To-File Mapping from: This%BinDataFile(2) = ", This%BinDataFile(2) )
  

  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


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
Subroutine Compute_Level_To_Bin_BinsContainer( This, LevelsContainer, i_Debug )

  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    iLevels, iBins, jBins
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Compute_Level_To_Bin_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) call Logger%Exiting

end Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteNLevels_BinsContainer( This, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)  ::    This
  logical                           ,optional ,intent(in)  ::    i_Debug
  
  integer                                                  ::    Status, Unit
  character(:)                      ,allocatable           ::    FileName
  logical                                                  ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteNLevels_BinsContainer" )
  !i_Debug_Loc   =     Logger%On()
  

  allocate( FileName, source= adjustl(trim( This%PathToBinMolFldr // '/../NLevels.inp' )) )
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName ) 
    write(Unit,'(I6)') This%NBins
  close(Unit)


  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteLevelsInBin_BinsContainer( This, LevelsContainer, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)     ::    This
  Type(LevelsContainer_Type)                  ,intent(in)     ::    LevelsContainer
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
  if (i_Debug_Loc) call Logger%Entering( "WriteLevelsInBin_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()  


  if (i_Debug_Loc)  call Logger%Write( "Finding the first (v,J)s for each of the bins (vectors vqnFirst and jqnFirst) ... ")
  do iBins = 1,This%NBins
    write(iBins_char,'(I5)') iBins

    FileName = This%PathToBinMolFldr // '/levels_Bin' // trim(adjustl(iBins_char)) // '.inp'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName ) 
    
    !write(Unit,'(A150)') ('######################################################################################################################################################')
    !write(Unit,'(A150)') ('# vqn   : the vibrational q.n. of the i''th quantum state                                                                                             ')
    !write(Unit,'(A150)') ('# jqn   : the rotational q.n. of the i''th quantum state                                                                                              ')
    !write(Unit,'(A150)') ('# eint  : internal energy of i''th quantum state                                                                                                      ')
    !write(Unit,'(A150)') ('# egam  : Half width of i''th quantum state                                                                                                           ')
    !write(Unit,'(A150)') ('# rmin  : the position of the potential minimum (included centrifugal potential) for i''th quantum state                                              ')
    !write(Unit,'(A150)') ('# vmin  : the value of the potential minimun (inc. cent. pot.)                                                                                        ')
    !write(Unit,'(A150)') ('# vmax  : the value of the local potential maximum (inc. cent. pot.)                                                                                  ')
    !write(Unit,'(A150)') ('# tau   : the vibrational period of the i''th quantum state                                                                                           ')
    !write(Unit,'(A150)') ('# ri    : inner turning point                                                                                                                         ')
    !write(Unit,'(A150)') ('# ro    : outter turning point                                                                                                                        ')
    !write(Unit,'(A150)') ('# rmax  : location of maximum in centrifugal barrier                                                                                                  ')
    !write(Unit,'(A150)') ('######################################################################################################################################################')
    !write(Unit,'(A150)') ('#    vqn  jqn   eint            egam            rmin          rmax          vmin            vmax            tau           ri             ro           ')
    !write(Unit,'(A150)') ('######################################################################################################################################################')
    
    !write(Unit,'(A)') ('#================================================================================================================================================')
    !write(Unit,'(A)') ('#')                                                                                              
    !write(Unit,'(A)') ('#================================================================================================================================================')
    !write(Unit,'(A)') ('# vqn, jqn,         E[Eh],      EGam[au],      rMin[a0],      rMax[a0],      VMin[Eh],      VMax[Eh],       Tau[au],       rIn[a0],      rOut[a0]')

    do iLevels = 1,LevelsContainer%NStates 
         
      if (LevelsContainer%States(iLevels)%To_Bin == iBins) then

        write(Unit,1) int(LevelsContainer%States(iLevels)%vqn), ',', &
                      int(LevelsContainer%States(iLevels)%jqn), ',', &
                          LevelsContainer%States(iLevels)%eint, ',', &
                          LevelsContainer%States(iLevels)%egam, ',', &  
                          LevelsContainer%States(iLevels)%rmin, ',', &
                          LevelsContainer%States(iLevels)%rmax, ',', &
                          LevelsContainer%States(iLevels)%vmin, ',', &
                          LevelsContainer%States(iLevels)%vmax, ',', &
                          LevelsContainer%States(iLevels)%tau,  ',', &
                          LevelsContainer%States(iLevels)%ri,   ',', &
                          LevelsContainer%States(iLevels)%ro
      
      end if
       
    end do
   
    close(Unit) 
       
  end do

  !1 format(1X, 2I5, 9es15.7)
  1 format(1X, I4, A, I4, 9(A, es14.7) )

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteQNsInBin_BinsContainer( This, LevelsContainer, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)     ::    This
  Type(LevelsContainer_Type)                  ,intent(in)     ::    LevelsContainer
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
  if (i_Debug_Loc) call Logger%Entering( "WriteQNsInBin_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()  


  if (i_Debug_Loc)  call Logger%Write( "Finding the first (v,J)s for each of the bins (vectors vqnFirst and jqnFirst) ... ")
  do iBins = 1,This%NBins
    write(iBins_char,'(I5)') iBins

    
    FileName = This%PathToBinMolFldr // '/../QNsInBin' // trim(adjustl(iBins_char)) // '.csv'
    open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // FileName ) 
    write(Unit,*) "#      vqn       jqn     En [eV]"

    do iLevels = 1,LevelsContainer%NStates 
         
      if (LevelsContainer%States(iLevels)%To_Bin == iBins) then

        write(Unit,'(3X  2I10, d20.10)') LevelsContainer%States(iLevels)%vqn, LevelsContainer%States(iLevels)%jqn, LevelsContainer%States(iLevels)%einteV
      
      end if
       
    end do
   
    close(Unit) 
       
  end do


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteQNs_To_Bin_BinsContainer( This, LevelsContainer, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)     ::    This
  Type(LevelsContainer_Type)                  ,intent(in)     ::    LevelsContainer
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
  if (i_Debug_Loc) call Logger%Entering( "WriteQNs_To_Bin_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()  


  FileName = This%PathToBinMolFldr // '/QNs_To_Bin.csv'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )  
    do ijqn = 0,LevelsContainer%maxjqn
      write(Unit,*) This%qns_to_Bin(:,ijqn)
    end do   
  close(Unit) 


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteQNsFirst_BinsContainer( This, LevelsContainer, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)     ::    This
  Type(LevelsContainer_Type)                  ,intent(in)     ::    LevelsContainer
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
  if (i_Debug_Loc) call Logger%Entering( "WriteQNsFirst_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()  


  FileName = This%PathToBinMolFldr // '/QNsFirst.csv'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )  
    write(Unit,'(A)') ("# vqn, jqn, No Levels,  % Levels")
    do iBins = 1,This%NBins
      write(Unit,1) This%Bin(iBins)%vqnFirst, ',', This%Bin(iBins)%jqnFirst, ',', This%Bin(iBins)%NLevels, ',', real(This%Bin(iBins)%NLevels) * 1.e2_rkp / real(LevelsContainer%NStates)
    end do   
  close(Unit) 

  1 format (I5, A, I4, A, I10, A, es11.4)

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine WriteQNsEnBin_BinsContainer( This, LevelsContainer, i_Debug )

  use Input_Class             ,only:  Input_Type
  use LevelsContainer_Class   ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(in)     ::    This
  Type(LevelsContainer_Type)                  ,intent(in)     ::    LevelsContainer
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
  if (i_Debug_Loc) call Logger%Entering( "WriteQNsEnBin_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()  


  FileName = This%PathToBinMolFldr // '/QNsEnBin.csv'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // FileName )  
    write(Unit,'(A)') ("#Level, vqn, jqn,        E [Eh],    Degeneracy,   Bin")
    do iLevels = 1,LevelsContainer%NStates 
      write(Unit,1) iLevels, ',', int(LevelsContainer%States(iLevels)%vqn),  ',', int(LevelsContainer%States(iLevels)%jqn), ',', LevelsContainer%States(iLevels)%einteV, ',', LevelsContainer%States(iLevels)%g, ',', LevelsContainer%States(iLevels)%To_Bin
    end do     
  close(Unit) 

  1 format (I6, A, I4, A, I4, A, es14.7, A, es14.7, A, I6)

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


  ! Write the Levels Contained in Each Bin
  !call This%WriteLevelsInBin( LevelsContainer, i_Debug=i_Debug_Loc )


  ! ! Write the Mapping QNs Contained in Each Bin
  ! call This%WriteQNsInBin( LevelsContainer, i_Debug=i_Debug_Loc )


  ! ! Write the Mapping QNs -> Bin
  ! call This%WriteQNs_To_Bin( LevelsContainer, i_Debug=i_Debug_Loc )
  
  
  ! Write the Mapping Bins -> First Quantum Number
  call This%WriteQNsFirst( LevelsContainer, i_Debug=i_Debug_Loc )


  ! Write the Mapping qns and Energy [eV]-> Bin
  call This%WriteQNsEnBin( LevelsContainer, i_Debug=i_Debug_Loc )

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ComputePartFunEnergy_BinsContainer( This, Input, LevelsContainer, iMol, i_Debug )

  use Input_Class              ,only:  Input_Type
  use Degeneracy_Factory_Class ,only:  Degeneracy_Factory_Type
  use Degeneracy_Class         ,only:  Compute_PartFunc_State
  use LevelsContainer_Class    ,only:  LevelsContainer_Type

  class(BinsContainer_Type)                   ,intent(inout)  ::    This
  Type(Input_Type)                            ,intent(in)     ::    Input
  Type(LevelsContainer_Type)                  ,intent(inout)  ::    LevelsContainer
  integer                                     ,intent(in)     ::    iMol
  logical                           ,optional ,intent(in)     ::    i_Debug
  
  Type(Degeneracy_Factory_Type)                               ::    Degeneracy_Factory
  integer                                                     ::    iBin
  integer                                                     ::    iLevels
  integer                                                     ::    Status
  integer                                                     ::    Unit
  character(:)            ,allocatable                        ::    FileName
  character(10)                                               ::    T_char
  logical                                                     ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ComputePartFunEnergy_BinsContainer")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
      
  if (i_Debug_Loc)  call Logger%Write( "iMol = ", iMol )  

  ! Computing Bins Effective Energy Minima (the Minima specified by the User could differ from actual Molecule Levels Energies)
  LevelsContainer%States%einteV_scaled = LevelsContainer%States%einteV - LevelsContainer%MineinteV
  This%Bin%MineinteV = 1000.d0
  do iLevels = 1,LevelsContainer%NStates
    iBin = LevelsContainer%States(iLevels)%To_Bin
    if (LevelsContainer%States(iLevels)%einteV < This%Bin(iBin)%MineinteV) then 
      This%Bin(iBin)%MineinteV = LevelsContainer%States(iLevels)%einteV
    end if
  end do
  
  
  This%Bin(:)%TotEinteV     = Zero
  This%Bin(:)%Q             = Zero
  This%Bin(:)%TotEinteVInit = Zero
  This%Bin(:)%QInit         = Zero


  do iLevels = 1,LevelsContainer%NStates
    LevelsContainer%States(iLevels)%QInit = Compute_PartFunc_State( LevelsContainer%States(iLevels)%g, LevelsContainer%States(iLevels)%einteV_scaled, Input%TInit )
    LevelsContainer%States(iLevels)%Q     = Compute_PartFunc_State( LevelsContainer%States(iLevels)%g, LevelsContainer%States(iLevels)%einteV_scaled, Input%TInt(iMol)  )

    iBin = LevelsContainer%States(iLevels)%To_Bin

    This%Bin(iBin)%TotEinteVInit = This%Bin(iBin)%TotEinteVInit + LevelsContainer%States(iLevels)%QInit * LevelsContainer%States(iLevels)%einteV_scaled
    This%Bin(iBin)%QInit         = This%Bin(iBin)%QInit         + LevelsContainer%States(iLevels)%QInit

    This%Bin(iBin)%TotEinteV     = This%Bin(iBin)%TotEinteV     + LevelsContainer%States(iLevels)%Q     * LevelsContainer%States(iLevels)%einteV_scaled
    This%Bin(iBin)%Q             = This%Bin(iBin)%Q             + LevelsContainer%States(iLevels)%Q
  end do
  
  This%Bin%TotEinteVInit = This%Bin%TotEinteVInit / This%Bin%QInit
  This%Bin%TotEinteV     = This%Bin%TotEinteV     / This%Bin%Q
  
  
  ! Computing the Bins Ratios Partition Functions
  This%Q = Zero
  do iBin = 1,This%NBins
    This%Q = This%Q         + This%Bin(iBin)%Q 
  end do
  do iBin = 1,This%NBins
    This%Bin(iBin)%QRatio     = This%Bin(iBin)%Q     / This%Q
  end do
  
  
  ! Computing the Bins Ratios Partition Functions
  This%QInit = Zero
  do iBin = 1,This%NBins
    This%QInit = This%QInit + This%Bin(iBin)%QInit 
  end do
  do iBin = 1,This%NBins
    This%Bin(iBin)%QRatioInit = This%Bin(iBin)%QInit / This%QInit
  end do
  
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%Tint(iMol))
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.csv'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,'(A)')"# Level/Bin Ratio,       Part. Func,      Energy [eV]"
  if (Status/=0) call Error( "Error opening file: " // FileName )                      
    do iBin = 1,This%NBins
      write(Unit,1) This%Bin(iBin)%QRatio, ',', This%Bin(iBin)%Q, ',', This%Bin(iBin)%TotEinteV
    end do   
  close(Unit) 
  
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%TInit)
  if (i_Debug_Loc)  call Logger%Write( "Write the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.csv'
  open( File=FileName, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,'(A)') "# Level/Bin Ratio,       Part. Func,      Energy [eV]"
  if (Status/=0) call Error( "Error opening file: " // FileName )                      
    do iBin = 1,This%NBins
      write(Unit,1) This%Bin(iBin)%QRatioInit, ',', This%Bin(iBin)%QInit, ',', This%Bin(iBin)%TotEinteVInit
    end do   
  close(Unit) 

  1 format (es18.10E3,A,es18.10E3,A,es18.10E3)

  
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
      
  
  This%Bin(:)%TotEinteV     = Zero
  This%Bin(:)%Q             = Zero
  This%Bin(:)%TotEinteVInit = Zero
  This%Bin(:)%QInit         = Zero
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%Tint(iMol))
  if (i_Debug_Loc)  call Logger%Write( "Reading the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.csv'
  if (i_Debug_Loc)  call Logger%Write( "Open File ", FileName ) 
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error Opening file: " // FileName )                      
    read(Unit,*)
    do iBins = 1,This%NBins
      read(Unit, *, iostat=Status) This%Bin(iBins)%QRatio, This%Bin(iBins)%Q, This%Bin(iBins)%TotEinteV
      if (Status/=0) call Error( "Error Reading file: " // FileName )  
    end do   
  close(Unit) 
  
  ! Writing the Mapping Bins -> First Quantum Numbers
  write(T_char,"(I10)") floor(Input%TInit)
  if (i_Debug_Loc)  call Logger%Write( "Reading the Mapping Bins -> First Quantum Numbers" )  
  FileName = This%PathToBinMolFldr // '/T' // trim(adjustl(T_char))// '.csv'
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error Opening file: " // FileName )  
    read(Unit,*)       
    do iBins = 1,This%NBins
      read(Unit, *, iostat=Status) This%Bin(iBins)%QRatioInit, This%Bin(iBins)%QInit, This%Bin(iBins)%TotEinteVInit
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
             trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBinsCG_char(iMol))) // '/QNsEnBin.csv'
  if (i_Debug_Loc)  call Logger%Write( "Opening the File: ", trim(adjustl(FileName)) )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error Opening File: " // FileName )  
  read(Unit,*)
  do iLevels = 1,LevelsContainer%NStates 
    !read(Unit,'(I8, 2I7, 2es15.7, I10)')  iTemp1, iTemp2, iTemp3, xTemp1, xTemp2, LevelsContainer%States(iLevels)%ToBin
    read(Unit,'(I8, 2I7, 2es15.7, I10)', iostat=Status)  iTemp1, iTemp2, iTemp3, xTemp1, xTemp2, LevelsContainer%States(iLevels)%To_Bin
    if (Status/=0) call Error( "Error Reading File: " // FileName )
  end do     
  close(Unit) 

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module