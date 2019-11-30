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

Module KONIG_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error, CheckVariable
  implicit none

  private
  public    ::    KONIGInquiringVariables
  public    ::    KONIGOpeningKinetics
  public    ::    KONIGCreatingDirs
  public    ::    KONIGFinishingWritingKinetics
  public    ::    KONIGWritingThermo
  public    ::    KONIGComputingWritingPop0
  public    ::    KONIGRunningExtCode

  logical   ,parameter    ::    i_Debug_Global = .false.

  contains

!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   INQUIRING KONIG ENVIRONMENTAL VARIABLES
! ==============================================================================================================
Subroutine KONIGInquiringVariables( Input, i_Debug )
  
  use Input_Class           ,only:  Input_Type
  
  class(Input_Type)                           ,intent(inout)  ::    Input
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    len
  integer                                                     ::    status
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "KONIGInquiringVariables" )
  !i_Debug_Loc   =     Logger%On()


  if (i_Debug_Loc) call Logger%Write( "Inquiring for $KONIG_SOURCE_DIR Environment Variable" )
  call get_environment_variable ("KONIG_SOURCE_DIR", Input%KONIGOrigDir, len, status, .true.)
  if (status .ge. 2) then
    if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $KONIGOrigDir: status = ", status )
    stop
  elseif (status .eq. 1) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_SOURCE_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_SOURCE_DIR: length = ", len, " truncated to 200")
    len = 200
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_SOURCE_DIR exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "KONIG's Original Code Folder, Input%KONIGOrigDir = ", Input%KONIGOrigDir)
  end if



  if (i_Debug_Loc) call Logger%Write( "Inquiring for $KONIG_EXECUTABLE Environment Variable" )
  call get_environment_variable ("KONIG_EXECUTABLE", Input%KONIGRunCMD, len, status, .true.)
  if (status .ge. 2) then
    if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $KONIGRunCMD: status = ", status )
    stop
  elseif (status .eq. 1) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_EXECUTABLE does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_EXECUTABLE: length = ", len, " truncated to 200")
    len = 200
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_EXECUTABLE exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "KONIG's Executable, Input%KONIGRunCMD = ", Input%KONIGRunCMD)
  end if


  if (i_Debug_Loc) call Logger%Write( "Inquiring for $KONIG_DATABASE_DIR Environment Variable" )
  call get_environment_variable ("KONIG_DATABASE_DIR", Input%KONIGDtbPath, len, status, .true.)
  if (status .ge. 2) then
    if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $KONIGDtbPath: status = ", status )
    stop
  elseif (status .eq. 1) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_DATABASE_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_DATABASE_DIR: length = ", len, " truncated to 200")
    len = 200
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$KONIG_DATABASE_DIR exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "KONIG's Database Folder, Input%KONIGDtbPath = ", Input%KONIGDtbPath)
  end if


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   OPENING KINETIC FILE 
! ==============================================================================================================
Subroutine KONIGOpeningKinetics( Input, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(inout)  ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    exist_flag
  character(:)                    ,allocatable              ::    FileName
  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "KONIGOpeningKinetics" )
  !i_Debug_Loc   =     Logger%On()

  FileName = './database/kinetics/'// trim(adjustl(Input%System))
  if ( i_Debug ) call Logger%Write( "Writing File: ", FileName )
  inquire(file=trim(adjustl(FileName)), exist=exist_flag)
  if (exist_flag) then
    open( File=FileName, Unit=122, status='OLD', position="append", action="write", iostat=Status )
  else
    open( File=FileName, Unit=122, status='REPLACE', iostat=Status )
    write(122,*) 'Units=cm^3/s'
  end if 

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   CREATING OUTPUT FOLDERS
! ==============================================================================================================
Subroutine KONIGCreatingDirs( Input, iMol, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iMol
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iTtra
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "KONIGCreatingDirs" )
  !i_Debug_Loc   =     Logger%On()

  call system( 'mkdir -p ./output/' )
  
  do iTtra = 1,Input%NTtra
  
    call system( 'mkdir -p ./output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) )
    call system( 'mkdir -p ./output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '/output' )
    call system( 'scp ' // trim(adjustl(Input%OutputDir)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                 trim(adjustl(Input%Molecules_Name(iMol))) // '_' //  trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // &
                 '.dat ' // trim(adjustl(Input%OutputDir)) // '/../output' )
  
  end do

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   FINISHING TO WRITE ARRHENIUS COEFFICIENTS
! ==============================================================================================================
Subroutine KONIGFinishingWritingKinetics( Input, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iBins, jBins
  character(6)                                              ::    iBinsChar, jBinsChar
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "KONIGFinishingWritingKinetics" )
  !i_Debug_Loc   =     Logger%On()

  if (Input%Kinetics_O2DissFlg) then
  
    do iBins =  1,Input%NBins(2)
      write(iBinsChar,'(I6)') iBins
      write(122,'(A)') 'O2(' // trim(adjustl(iBinsChar)) // ')+O=O+O+O:+0.0166E+00,-1.5E+00,+5975E+01,2'
      write(122,'(A)') 'O2(' // trim(adjustl(iBinsChar)) // ')+C=O+O+C:+0.0166E+00,-1.5E+00,+5975E+01,2'
    end do
    
    do iBins = 1,Input%NBins(1)
      write(iBinsChar,'(I6)') iBins
      do jBins = iBins,Input%NBins(2)
        write(jBinsChar,'(I6)') jBins
        write(122,'(A)') 'O2(' // trim(adjustl(jBinsChar)) // ')+CO(' // trim(adjustl(iBinsChar)) // ')=O+O+CO(' // trim(adjustl(iBinsChar)) // '):+0.0033E+00,-1.5E+00,+5975E+01,2'
      end do
    end do
  
    do iBins = 1,Input%NBins(2)
      do jBins = iBins,Input%NBins(2)
        write(iBinsChar,'(I6)') iBins
        write(jBinsChar,'(I6)') jBins
        write(122,'(A)') 'O2(' // trim(adjustl(iBinsChar)) // ')+O2(' // trim(adjustl(jBinsChar)) // ')=O+O+O2:+0.0033E+00,-1.5E+00,+5975E+01,2'
      end do
    end do
  
  end if
        
  close(122)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   WRITING THERMO PROPERTIES
! ==============================================================================================================
Subroutine KONIGWritingThermo( Input, BinnedMolecule, i_Debug )

  use Input_Class           ,only:  Input_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  Type(BinsContainer_Type)   ,dimension(:)  ,intent(inout)  ::    BinnedMolecule
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    NEnLevels
  integer                                                   ::    iTtra, iBins, iMol
  character(6)                                              ::    NEnLevelsChar
  character(:)                    ,allocatable              ::    FileName
  character(:)                    ,allocatable              ::    ThermoFileOrig
  character(:)                    ,allocatable              ::    ThermoFileNew
  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "KONIGWritingThermo" )
  !i_Debug_Loc   =     Logger%On()

  do iMol = 1,Input%NMolecules
    if ( i_Debug ) call Logger%Write( "Molecule Nb", iMol )
    
    do iTtra = 1,Input%NTtra
    
      FileName = trim(adjustl(Input%OutputDir)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                 trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '.dat'
      open( File=FileName, Unit=127, status='OLD', iostat=Status )
      if (Status/=0) call Error( "Error opening file: " // FileName )
      if ( i_Debug ) call Logger%Write( "Opening file: " // FileName // " for reading." )   
      read(127,*)    
      
      ThermoFileOrig = trim(adjustl(Input%KONIGDtbPath)) // '/thermo/' // trim(adjustl(Input%Molecules_Name(iMol))) // '_Format'
      ThermoFileNew  = './database/thermo/' // trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra)))
      call system( 'rm -rf ' // trim(adjustl(ThermoFileNew)) )
      call system( 'scp ' // trim(adjustl(ThermoFileOrig )) // ' ' // trim(adjustl(ThermoFileNew)) )
      if ( i_Debug ) call Logger%Write( "scp " // trim(adjustl(ThermoFileOrig )) // " " // trim(adjustl(ThermoFileNew)) ) 
      
      open( File=ThermoFileNew, Unit=142, status='OLD', position="append", action="write", iostat=Status ) 
      if (Status/=0) call Error( "Error opening file: " // ThermoFileNew )
      if ( i_Debug ) call Logger%Write( "Opening file: " // trim(adjustl(ThermoFileNew)) // " for writing." )   
      
      NEnLevels = Input%BinFinal(iMol) - Input%BinStart(iMol) + 1
      write(NEnLevelsChar,'(I6)') NEnLevels
      write(142,'(A)') 'NB_ENERGY_LEVELS = ' // trim(adjustl(NEnLevelsChar))
      
      
      do iBins = 1,Input%NBins(iMol)
        read(127,'(3es20.10E3)') BinnedMolecule(iMol)%Bin(iBins)%QRatio, BinnedMolecule(iMol)%Bin(iBins)%Q,  BinnedMolecule(iMol)%Bin(iBins)%ToteinteV
      end do
      close(127)
          
          
      do iBins = Input%BinStart(iMol),Input%BinFinal(iMol)
        if (trim(Input%BSortMethod(iMol)) .eq. "State-Specific" ) then
          write(142,'(e14.8,3x,e15.8)') BinnedMolecule(iMol)%Bin(iBins)%Q,  BinnedMolecule(iMol)%Bin(iBins)%ToteinteV
        else
          write(142,'(e14.8,3x,e15.8)') BinnedMolecule(iMol)%Bin(iBins)%Q,  Zero !BinnedMolecule(iMol)%Bin(iBins)%ToteinteV
        end if
      end do
      close(142)
    
    
    end do   
    
  end do   ! iMol

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   COMPUTING AND WRITING INITIAL POPULATIONS
! ==============================================================================================================
Subroutine KONIGComputingWritingPop0( Input, BinnedMolecule, i_Debug )

  use Input_Class           ,only:  Input_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Parameters_Module     ,only:  Zero, One

  type(Input_Type)                          ,intent(in)     ::    Input
  Type(BinsContainer_Type)   ,dimension(:)  ,intent(inout)  ::    BinnedMolecule
  logical                         ,optional ,intent(in)     ::    i_Debug

  real(rkp)                                                 ::    QIniTot
  integer                                                   ::    NEnLevels
  integer                                                   ::    iComp, iBins, iMol
  character(:)                    ,allocatable              ::    FileName1
  character(:)                    ,allocatable              ::    FileName2
  logical                                                   ::    MolFlag
  integer                                                   ::    Status
  real(rkp)                                                 ::    MinEn
  integer                                                   ::    iBins0
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "KONIGComputingWritingPop0" )
  !i_Debug_Loc   =     Logger%On()

  !FileName1 = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/InitialMolFracs.dat'
  FileName1 = './InitialMolFracs.dat'
  open( File=FileName1, Unit=132, status='REPLACE', iostat=Status )
  
    write(132,*)     '# Nb of Components'
    write(132,'(I20)') Input%NComp
    
    write(132,*)     '# Nb of Bins per Component'
    do iComp = 1,Input%NComp
      MolFlag = .false.
      iMol=1
      do while ((iMol <= Input%NMolecules) .and. (MolFlag .eqv. .false.))
        if ( trim(adjustl(Input%CompName(iComp))) == trim(adjustl(Input%Molecules_Name(iMol))) ) then
          MolFlag = .true.
        else
          iMol=iMol+1
        end if
      end do
      if ( MolFlag ) then
        NEnLevels = Input%BinFinal(iMol) - Input%BinStart(iMol) + 1
        write(132,'(I20)') NEnLevels
      else 
        write(132,'(I20)') 0
      end if
    end do
    
    write(132,*)     '# Initial Mole Fractions at ', trim(adjustl(Input%TInit_char)), 'K for KONIG'
    do iComp = 1,Input%NComp   
      MolFlag = .false.
      iMol=1
      do while ((iMol <= Input%NMolecules) .and. (MolFlag .eqv. .false.))
        if ( trim(adjustl(Input%CompName(iComp))) == trim(adjustl(Input%Molecules_Name(iMol))) ) then
          MolFlag = .true.
        else
          iMol=iMol+1
        end if
      end do
      if ( MolFlag ) then

        MinEn  = 1.d20
        iBins0 = 1 
        do iBins = 1,Input%NBins(iMol)
          if (BinnedMolecule(iMol)%Bin(iBins)%ToteinteVInit < MinEn) then
            MinEn  = BinnedMolecule(iMol)%Bin(iBins)%ToteinteVInit
            iBins0 = iBins
          end if
        end do
        
        ! QIniTot = Zero
        ! FileName2 = trim(adjustl(Input%OutputDir)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
        !             trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(Input%TInit_char))// '.dat'
        ! open( File=FileName2, Unit=131, status='OLD', iostat=Status )
        ! read(131,*)
        ! do iBins = 1,Input%NBins(iMol)
        !   read(131,'(3es20.10E3)')  BinnedMolecule(iMol)%Bin(iBins)%QRatioInit, BinnedMolecule(iMol)%Bin(iBins)%QInit, BinnedMolecule(iMol)%Bin(iBins)%ToteinteVInit
        !   if ( (iBins >= Input%BinStart(iMol) ) .and. (iBins <= Input%BinFinal(iMol)) ) then
        !     QIniTot = QIniTot + BinnedMolecule(iMol)%Bin(iBins)%QInit
        !   end if
        ! end do
        ! close(131)
        
        ! do iBins = Input%BinStart(iMol),Input%BinFinal(iMol)
        !   write(132,'(d20.10)') max(BinnedMolecule(iMol)%Bin(iBins)%QInit / QIniTot,1.d-99)
        ! end do

        do iBins = Input%BinStart(iMol),Input%BinFinal(iMol)
          if (iBins == iBins0) then 
            write(132,'(d20.10)') 1.d0
          else
            write(132,'(d20.10)') 1.d-20
          end if
        end do

      else
        write(132,'(d20.10)') One
      end if
    end do
    
  close(132)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   RUNNING EXTERNAL CODE
! ==============================================================================================================
Subroutine KONIGRunningExtCode( Input, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iTtra
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "KONIGRunningExtCode" )
  !i_Debug_Loc   =     Logger%On()
        
  do iTtra = 1,Input%NTtra

    !call system( 'cd ./output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '/output && ' // trim(adjustl(Input%KONIGRunCMD)) // ' ' // &
    !             trim(adjustl(Input%KONIGOrigDir)) // '/input/' // trim(adjustl(Input%KONIGInputFileName)) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ' ' // trim(adjustl(Input%NProcForKONIGChar)) // ' ' // trim(adjustl(Input%NProcForKONIGChar)))
    
    !if ( i_Debug ) call Logger%Write( "cd ./output/T_" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "/output && " // trim(adjustl(Input%KONIGRunCMD)) // " " // &
    !             trim(adjustl(Input%KONIGOrigDir)) // "/input/" // trim(adjustl(Input%KONIGInputFileName)) // "_" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ' ' // trim(adjustl(Input%NProcForKONIGChar)) // ' ' // trim(adjustl(Input%NProcForKONIGChar)))
    
    
    call chdir( './output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '/output' )
    call system( trim(adjustl(Input%KONIGRunCMD)) // ' ' // trim(adjustl(Input%NProcExtCodeChar)) // ' ' // trim(adjustl(Input%NProcExtCodeChar)) )
    
    if ( i_Debug ) call Logger%Write( 'Running System Command: ' // trim(adjustl(Input%KONIGRunCMD)) // ' ' // trim(adjustl(Input%NProcExtCodeChar)) // ' ' // trim(adjustl(Input%NProcExtCodeChar)) )
                 
    !call system( 'scp ./output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '/output/box.dat    ../../output')
    !call system( 'scp ./output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '/output/pop_' // trim(adjustl(Input%Molecules_Name(1))) // '.dat ../../output')
    !call system( 'scp ./output/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '/output/Tint.dat   ../../output')

  end do   ! iTtra
  
  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
