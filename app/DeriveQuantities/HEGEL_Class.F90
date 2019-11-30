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

Module HEGEL_Class

  use Parameters_Module     ,only:  rkp, Zero, One, UKb, Ue, Una, Hartree_To_eV, au_To_kgm
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error, CheckVariable
  implicit none

  private
  public    :: HEGELInquiringVariables
  public    :: HEGELOpeningKinetics
  public    :: HEGELFinishingWritingKinetics
  public    :: HEGELComputeTables

  logical   ,parameter    ::    i_Debug_Global = .false.

  contains

!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   INQUIRING HEGEL ENVIRONMENTAL VARIABLES
! ==============================================================================================================
Subroutine HEGELInquiringVariables( Input, i_Debug )
  
  use Input_Class           ,only:  Input_Type
  
  class(Input_Type)                           ,intent(inout)  ::    Input
  logical                           ,optional ,intent(in)     ::    i_Debug

  integer                                                     ::    len
  integer                                                     ::    status
  logical                                                     ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "HEGELInquiringVariables" )
  !i_Debug_Loc   =     Logger%On()


  if (i_Debug_Loc) call Logger%Write( "Inquiring for $HEGEL_SOURCE_DIR Environment Variable" )
  call get_environment_variable ("HEGEL_SOURCE_DIR", Input%HEGELOrigDir, len, status, .true.)
  if (status .ge. 2) then
    if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $HEGELOrigDir: status = ", status )
    stop
  elseif (status .eq. 1) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_SOURCE_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_SOURCE_DIR: length = ", len, " truncated to 200")
    len = 200
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_SOURCE_DIR exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "HEGEL's Original Code Folder, Input%HEGELOrigDir = ", Input%HEGELOrigDir)
  end if



  if (i_Debug_Loc) call Logger%Write( "Inquiring for $HEGEL_EXECUTABLE Environment Variable" )
  call get_environment_variable ("HEGEL_EXECUTABLE", Input%HEGELRunCMD, len, status, .true.)
  if (status .ge. 2) then
    if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $HEGELRunCMD: status = ", status )
    stop
  elseif (status .eq. 1) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_EXECUTABLE does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_EXECUTABLE: length = ", len, " truncated to 200")
    len = 200
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_EXECUTABLE exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "HEGEL's Executable, Input%HEGELRunCMD = ", Input%HEGELRunCMD)
  end if


  if (i_Debug_Loc) call Logger%Write( "Inquiring for $HEGEL_DATABASE_DIR Environment Variable" )
  call get_environment_variable ("HEGEL_DATABASE_DIR", Input%HEGELDtbPath, len, status, .true.)
  if (status .ge. 2) then
    if (i_Debug_Loc) call Logger%Write( "get_environment_variable failed for $HEGELDtbPath: status = ", status )
    stop
  elseif (status .eq. 1) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_DATABASE_DIR does not exist" )
    stop
  elseif (status .eq. -1) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_DATABASE_DIR: length = ", len, " truncated to 200")
    len = 200
  end if
  if (len .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "$HEGEL_DATABASE_DIR exists, but has no value")
    stop
  end if
  if (status .eq. 0) then
    if (i_Debug_Loc) call Logger%Write( "HEGEL's Database Folder, Input%HEGELDtbPath = ", Input%HEGELDtbPath)
  end if


  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   OPENING KINETIC FILE 
! ==============================================================================================================
Subroutine HEGELOpeningKinetics( Input, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(inout)  ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    exist_flag
  character(:)                    ,allocatable              ::    FileName
  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "HEGELOpeningKinetics" )
  !i_Debug_Loc   =     Logger%On()

  FileName = './database/kinetics/'// trim(adjustl(Input%System))
  if ( i_Debug ) call Logger%Write( "Writing File: ", FileName )
  inquire(file=trim(adjustl(FileName)), exist=exist_flag)
  if (exist_flag) then
    open( File=FileName, Unit=132, status='OLD', position="append", action="write", iostat=Status )
  else
    open( File=FileName, Unit=132, status='REPLACE', iostat=Status )
    write(132,*) 'Units=cm^3/s'
  end if 

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!   FINISHING TO WRITE ARRHENIUS COEFFICIENTS
! ==============================================================================================================
Subroutine HEGELFinishingWritingKinetics( Input, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iBins, jBins
  character(6)                                              ::    iBinsChar, jBinsChar
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "HEGELFinishingWritingKinetics" )
  !i_Debug_Loc   =     Logger%On()

  if (Input%Kinetics_O2DissFlg) then
  
    do iBins =  1,Input%NBins(2)
      write(iBinsChar,'(I6)') iBins
      write(132,'(A)') 'O2(' // trim(adjustl(iBinsChar)) // ')+O=O+O+O:+0.0166E+00,-1.5E+00,+5975E+01,2'
      write(132,'(A)') 'O2(' // trim(adjustl(iBinsChar)) // ')+C=O+O+C:+0.0166E+00,-1.5E+00,+5975E+01,2'
    end do
    
    do iBins = 1,Input%NBins(1)
      write(iBinsChar,'(I6)') iBins
      do jBins = iBins,Input%NBins(2)
        write(jBinsChar,'(I6)') jBins
        write(132,'(A)') 'O2(' // trim(adjustl(jBinsChar)) // ')+CO(' // trim(adjustl(iBinsChar)) // ')=O+O+CO(' // trim(adjustl(iBinsChar)) // '):+0.0033E+00,-1.5E+00,+5975E+01,2'
      end do
    end do
  
    do iBins = 1,Input%NBins(2)
      do jBins = iBins,Input%NBins(2)
        write(iBinsChar,'(I6)') iBins
        write(jBinsChar,'(I6)') jBins
        write(132,'(A)') 'O2(' // trim(adjustl(iBinsChar)) // ')+O2(' // trim(adjustl(jBinsChar)) // ')=O+O+O2:+0.0033E+00,-1.5E+00,+5975E+01,2'
      end do
    end do
  
  end if
        
  close(132)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!> This function computes the partition function according to the definition from statistical mechanics
! ==============================================================================================================
Subroutine HEGELComputeTables( Input, LevelsContainer, BinnedMolecule, Species, i_Debug )

  use Input_Class           ,only:  Input_Type
  use LevelsContainer_Class ,only:  LevelsContainer_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Species_Class         ,only:  Species_Type

  Type(Input_Type)                            ,intent(in)     :: Input
  Type(LevelsContainer_Type)                  ,intent(inout)  :: LevelsContainer
  Type(BinsContainer_Type)                    ,intent(inout)  :: BinnedMolecule
  type(Species_Type)         ,dimension(:)    ,intent(in)     :: Species
  logical                           ,optional ,intent(in)     :: i_Debug

  real(rkp)          ,dimension(:),   allocatable             :: TGroupTable
  real(rkp)          ,dimension(:,:), allocatable             :: EMat, GMat
  real(rkp)          ,dimension(:,:), allocatable             :: EBinTable, CvBinTable, QBinTable

  integer                                                     :: iMol
  integer                                                     :: NLevels = 0
  integer                                                     :: iLevels
  integer            ,dimension(:),   allocatable             :: NLevelsPerGroup
  integer                                                     :: iBins
  character(len=10)                                           :: iBins_char
  character(len=150)                                          :: PathToTermo

  real(rkp)                                                   :: MolWeight 
  real(rkp)                                                   :: MolMass  
  character(len=4)                                            :: MoleculeName
  real(rkp)                                                   :: EMin
  real(rkp)                                       ,parameter  :: Tmin      = 250.0_rkp
  real(rkp)                                       ,parameter  :: dT        = 25.0_rkp
  integer                                         ,parameter  :: NTablePts = 800
  integer                                                     :: iPts
  integer                                                     :: FileUnit
  character(len=100)                                          :: DataFile
  integer                                                     :: astat, ios

  logical                                                     :: i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if ( i_Debug_Loc ) call Logger%Entering( "HEGELComputeTables" )
  !i_Debug_Loc   =     Logger%On()

  iMol = 1

  do iBins = 1,Input%NBins(1)
    NLevels = NLevels + BinnedMolecule%Bin(iBins)%NLevels
  end do

  MolWeight    = Species(1)%Mass * au_To_kgm
  MolMass      = MolWeight / Una 
  MoleculeName = Species(1)%Name

  !-----------------------!
  ! Allocate memory 
  allocate( NLevelsPerGroup(Input%NBins(iMol)), stat=astat )
  if (astat .ne. 0) stop "ERROR while allocating: 'NLevelsPerGroup'"
  NLevelsPerGroup = 0

  allocate( EMat(NLevels, Input%NBins(1)), stat=astat )
  if (astat .ne. 0) stop "ERROR while allocating: 'EMat'"
  EMat = Zero

  allocate( GMat(NLevels, Input%NBins(1)), stat=astat )
  if (astat .ne. 0) stop "ERROR while allocating: 'GMat'"
  GMat = Zero

  allocate( TGroupTable(NTablePts), stat=astat )
  if (astat .ne. 0) stop "ERROR while allocating: 'TGroupTable'"
  TGroupTable = Zero

  allocate( EBinTable(NTablePts, Input%NBins(1)), stat=astat )
  if (astat.ne.0) stop "ERROR while allocating: 'EBinTable'"
  EBinTable = Zero

  allocate( CvBinTable(NTablePts, Input%NBins(1)), stat=astat)
  if (astat .ne. 0) stop "ERROR while allocating: 'CvBinTable'"
  CvBinTable = Zero

  allocate( QBinTable(NTablePts, Input%NBins(1)), stat=astat )
  if (astat .ne. 0) stop "ERROR while allocating: 'QBinTable'"
  QBinTable = Zero


  TGroupTable(1) = Tmin
  do iPts = 2,NTablePts
     TGroupTable(iPts) = TGroupTable(iPts - 1) + dT
  end do

  do iLevels = 1,NLevels
    iBins                                 = LevelsContainer%States(iLevels)%ToBin
    NLevelsPerGroup(iBins)                = NLevelsPerGroup(iBins) + 1

    LevelsContainer%States(iLevels)%eintK = (LevelsContainer%States(iLevels)%eint*Hartree_To_eV - LevelsContainer%mineinteV) * Ue / Ukb
    EMat(NLevelsPerGroup(iBins),iBins)    = LevelsContainer%States(iLevels)%eintK

    GMat(NLevelsPerGroup(iBins),iBins)    = LevelsContainer%States(iLevels)%g
  end do   


  ! Subtract the energy of the first group (better conditioning)
  do iBins = 1,Input%NBins(1)

     EMin = EMat(1,iBins)

     EMat(1:NLevelsPerGroup(iBins),iBins) = EMat(1:NLevelsPerGroup(iBins),iBins) - EMin

     !print*, iBins, EMat(1,iBins) * Ukb * Una
  end do


  ! Formation energy of N (J/mol)
  !print*
  !print*,-0.5d0 * LevelsContainer%States(1)%eint * Hartree_To_eV * UE * Una

  !-----------------------!
  ! Compute thermo data
  do iBins = 1,Input%NBins(1)
     write(iBins_char,'(i4)') iBins
     !write(*,*) "Computing THERMO properties for GROUP " // trim(adjustl(iBins_char))

     do iPts = 1,NTablePts

        QBinTable(iPts,iBins) = BoltzPartFunction(NLevelsPerGroup(iBins), TGroupTable(iPts), GMat(1:NLevelsPerGroup(iBins),iBins), EMat(1:NLevelsPerGroup(iBins),iBins))
        
        call BoltzEnergyCv(NLevelsPerGroup(iBins), MolMass, TGroupTable(iPts), GMat(1:NLevelsPerGroup(iBins),iBins), EMat(1:NLevelsPerGroup(iBins),iBins), EBinTable(iPts,iBins), CvBinTable(iPts,iBins))
        
     end do
     
  end do  

  PathToTermo = './database/thermo/' // trim(adjustl(Input%NBins_char(1))) // 'g'
  call system( 'mkdir -p ' // trim(adjustl(PathToTermo)) )

  !-----------------------!
  ! Write temperature tables
  DataFile = trim(adjustl(PathToTermo)) // '/Tg_table.dat'
  open(NewUnit=FileUnit, file=trim(DataFile), status='unknown', form='formatted', action='write')
  do iPts = 1,NTablePts 
     write(FileUnit, 10) TGroupTable(iPts)
  enddo
  close(FileUnit)

  DataFile = trim(adjustl(PathToTermo)) // '/Tt_table.dat'
  open(NewUnit=FileUnit, file=trim(DataFile), status='unknown', form='formatted', action='write')
  do iPts = 1,NTablePts 
    write(FileUnit, 10) TGroupTable(iPts)
  enddo
  close(FileUnit)


  !-----------------------!
  ! Write group thermo data in binary format
  do iBins = 1,Input%NBins(1)
     write(iBins_char, '(i4)') iBins
     !write(*,*) "Writing THERMO properties for GROUP " // trim(adjustl(iBins_char))

     ! - partition function
     DataFile = trim(adjustl(PathToTermo)) // '/' // trim(adjustl(MoleculeName)) // '_' // trim(adjustl(iBins_char)) // '_qint.dat'
     open(NewUnit=FileUnit, file=trim(DataFile), status='replace', form='unformatted', access='direct',recl=rkp*NTablePts, action='write', CONVERT="BIG_ENDIAN")
     write(FileUnit,rec=1) ( QBinTable(iPts, iBins), iPts = 1,NTablePts )
     close(FileUnit)

     ! - energy per unit mass [J/kg]
     DataFile = trim(adjustl(PathToTermo)) // '/' // trim(adjustl(MoleculeName)) // '_' // trim(adjustl(iBins_char)) // '_eint.dat'
     open(NewUnit=FileUnit, file=trim(DataFile), status='replace', form='unformatted', access='direct', recl=rkp*NTablePts, action='write', CONVERT="BIG_ENDIAN")
     write(FileUnit,rec=1) ( EBinTable(iPts, iBins), iPts = 1,NTablePts )
     close(FileUnit)

     ! - specific heat per unit mass [J/(kg * K)]
     DataFile = trim(adjustl(PathToTermo)) // '/' // trim(adjustl(MoleculeName)) // '_' // trim(adjustl(iBins_char)) // '_cvint.dat'
     open(NewUnit=FileUnit, file=trim(DataFile), status='replace', form='unformatted', access='direct', recl=rkp*NTablePts, action='write', CONVERT="BIG_ENDIAN")
     write(FileUnit,rec=1) ( CvBinTable(iPts, iBins), iPts = 1,NTablePts )
     close(FileUnit)

  enddo
  ! write(*,*) 'Temperature = ', TGroupTable(391)
  ! write(*,*) 'QBinTable   = ', QBinTable(391, :)
  ! write(*,*) 'EBinTable   = ', EBinTable(391, :)
  ! write(*,*) 'CvBinTable  = ', CvBinTable(391, :)
  ! write(*,*)

  ! !-----------------------!
  ! ! Write kinetic file
  ! DataFile = 'output/nitrogen_15g'
  ! open(NewUnit=FileUnit,file=trim(DataFile),status='unknown',action='write')
  ! write(FileUnit,5)"Units=cm^3/s"
  ! ! N3 excitation
  ! kin_file = 'input/data_'//trim(Input%NBins(1)_char)//'groups/kf_Exh_N3_'//trim(Input%NBins(1)_char)//'.dat'
  ! open(NewUnit=FileUnit2,file=trim(kin_file),status='old',form='formatted',action='read',iostat=ios)
  ! if (ios.ne.0) then
  !    write(*,5)"ERROR:: file '"//trim(kin_file)//"' not found"     
  !    stop  
  ! endif

  ! read(FileUnit2,*)n_proc
  ! do p = 1,n_proc

  !    read(FileUnit2,*)i,j,afit,bfit,cfit

  !    write(i_char,'(i4)')i
  !    i_char = adjustl(i_char)

  !    write(j_char,'(i4)')j
  !    j_char = adjustl(j_char)

  !    write(af_char,'(e14.6)')afit
  !    af_char = adjustl(af_char)

  !    write(bf_char,'(e14.6)')bfit
  !    bf_char = adjustl(bf_char)

  !    write(cf_char,'(e14.6)')cfit
  !    cf_char = adjustl(cf_char)

  !    write(FileUnit,5)'N2_'//trim(i_char)//'+N=N2_'//trim(j_char)//'+N:'//trim(af_char)//','//trim(bf_char)//','//trim(cf_char)//',5'
  ! enddo
  ! close(FileUnit2)

  ! ! N3 dissociation
  ! kin_file = 'input/data_'//trim(Input%NBins(1)_char)//'groups/kf_Dh_N3_'//trim(Input%NBins(1)_char)//'.dat'
  ! open(NewUnit=FileUnit2,file=trim(kin_file),status='old',form='formatted',action='read',iostat=ios)
  ! if (ios.ne.0) then
  !    write(*,5)"ERROR:: file '"//trim(kin_file)//"' not found"     
  !    stop  
  ! endif
  ! read(FileUnit2,*)n_proc
  ! do p = 1,n_proc
  !    read(FileUnit2,*)i,afit,bfit,cfit

  !    write(i_char,'(i4)')i
  !    i_char = adjustl(i_char)

  !    write(af_char,'(e14.6)')afit
  !    af_char = adjustl(af_char)

  !    write(bf_char,'(e14.6)')bfit
  !    bf_char = adjustl(bf_char)

  !    write(cf_char,'(e14.6)')cfit
  !    cf_char = adjustl(cf_char)

  !    write(FileUnit,5)'N2_'//trim(i_char)//'+N=N+N+N:'//trim(af_char)//','//trim(bf_char)//','//trim(cf_char)//',2'
  ! enddo           

  ! ! N4 excitation
  ! kin_file = 'input/data_'//trim(Input%NBins(1)_char)//'groups/kf_Exh_N4_'//trim(Input%NBins(1)_char)//'.dat'
  ! open(NewUnit=FileUnit2,file=trim(kin_file),status='old',form='formatted',action='read',iostat=ios)
  ! if (ios.ne.0) then
  !    write(*,5)"ERROR:: file '"//trim(kin_file)//"' not found"     
  !    stop  
  ! endif
  ! read(FileUnit2,*)n_proc
  ! do p = 1,n_proc
  !    read(FileUnit2,*)i,j,k,l,afit,bfit,cfit

  !    write(i_char,'(i4)')i
  !    i_char = adjustl(i_char)

  !    write(j_char,'(i4)')j
  !    j_char = adjustl(j_char)

  !    write(k_char,'(i4)')k
  !    k_char = adjustl(k_char)

  !    write(l_char,'(i4)')l
  !    l_char = adjustl(l_char)

  !    write(af_char,'(e14.6)')afit
  !    af_char = adjustl(af_char)

  !    write(bf_char,'(e14.6)')bfit
  !    bf_char = adjustl(bf_char)

  !    write(cf_char,'(e14.6)')cfit
  !    cf_char = adjustl(cf_char)

  !    write(FileUnit,5)'N2_'//trim(i_char)//'+N2_'//trim(j_char)//'='// 'N2_'//trim(k_char)//'+N2_'//trim(l_char)//':'//trim(af_char)//','//trim(bf_char)//','//trim(cf_char)//',6'
  ! enddo        

  ! ! N4 combined excitation and dissociation
  ! kin_file = 'input/data_'//trim(Input%NBins(1)_char)//'groups/kf_DExh_N4_'//trim(Input%NBins(1)_char)//'.dat'
  ! open(NewUnit=FileUnit2,file=trim(kin_file),status='old',form='formatted',action='read',iostat=ios)
  ! if (ios.ne.0) then
  !    write(*,5)"ERROR:: file '"//trim(kin_file)//"' not found"     
  !    stop  
  ! endif
  ! read(FileUnit2,*)n_proc
  ! do p = 1,n_proc
  !    read(FileUnit2,*)i,j,k,afit,bfit,cfit

  !    write(i_char,'(i4)')i
  !    i_char = adjustl(i_char)

  !    write(j_char,'(i4)')j
  !    j_char = adjustl(j_char)

  !    write(k_char,'(i4)')k
  !    k_char = adjustl(k_char)

  !    write(af_char,'(e14.6)')afit
  !    af_char = adjustl(af_char)

  !    write(bf_char,'(e14.6)')bfit
  !    bf_char = adjustl(bf_char)

  !    write(cf_char,'(e14.6)')cfit
  !    cf_char = adjustl(cf_char)

  !    write(FileUnit,5)'N2_'//trim(i_char)//'+N2_'//trim(j_char)//'='// 'N+N+N2_'//trim(k_char)//':'//trim(af_char)//','//trim(bf_char)//','//trim(cf_char)//',2'
  ! enddo        
  ! close(FileUnit)

  !-----------------------!
  ! Deallocate memory
  deallocate(NLevelsPerGroup,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'NLevelsPerGroup'"

  deallocate(EMat,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'EMat'"

  deallocate(GMat,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'GMat'"

  deallocate(TGroupTable,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'TGroupTable'"

  deallocate(EBinTable,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'EBinTable'"

  deallocate(CvBinTable,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'CvBinTable'"

  deallocate(QBinTable,stat=astat)
  if (astat.ne.0) stop "ERROR while deallocating: 'QBinTable'"

  5  format(a)
  10 format(e14.6)

  if ( i_Debug_Loc ) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!> This function computes the partition function according to the definition from statistical mechanics
! ==============================================================================================================
Pure Function BoltzPartFunction( NLev, T, gi, Ei ) result( Q )

  integer                  ,intent(in) :: NLev    !< number of energy levels
  real(rkp)                ,intent(in) :: T       !< temperature [K]
  real(rkp) ,dimension(:)  ,intent(in) :: gi      !< degeneracies
  real(rkp) ,dimension(:)  ,intent(in) :: Ei      !< energy levels [K]
  real(rkp)                            :: Q       !< partition function 

  integer                              :: i
  real(rkp)                            :: ov_T

  ! Common factor
  ov_T = One / T
  
  ! Loop over the energy levels
  Q = Zero
  do i = 1,NLev
     Q = Q + gi(i) * exp( - Ei(i) * ov_T )
  end do

End Function 
!--------------------------------------------------------------------------------------------------------------------------------!


!______________________________________________________________________________________________________________!
! ==============================================================================================================
!> This function computes the energy and specific heat per unit mass for a Boltzmann distribution
! ==============================================================================================================
Pure Subroutine BoltzEnergyCv( NLev, mass, T, gi, Ei, e, cv )

  integer                  ,intent(in)  :: NLev    !< number of energy levels
  real(rkp)                ,intent(in)  :: mass    !< mass [kg]
  real(rkp)                ,intent(in)  :: T       !< temperature [K]
  real(rkp) ,dimension(:)  ,intent(in)  :: gi      !< degeneracies 
  real(rkp) ,dimension(:)  ,intent(in)  :: Ei      !< energies [K]
  real(rkp)                ,intent(out) :: e       !< energy [J/kg]
  real(rkp)                ,intent(out) :: cv      !< specific heat [J/(kg*K)]

  integer                               :: i
  real(rkp)                             :: dexp0, Edexp0, Elev, ov_T, ov_Ts
  real(rkp)                             :: sum1, sum2, sum3 

  ! Common factor
  ov_T  = One / T
  ov_Ts = ( ov_T / T ) * ( Ukb / mass )
  
  ! Loop over the energy levels
  sum1 = Zero
  sum2 = Zero
  sum3 = Zero
  do i = 1,NLev
     Elev   = Ei(i)
     dexp0  = gi(i) * exp( - Elev * ov_T )
     Edexp0 = Elev  * dexp0 
     sum1   = sum1  + dexp0
     sum2   = sum2  + Edexp0
     sum3   = sum3  + Elev * Edexp0
  enddo

  e  = ( sum2 / sum1 ) * ( Ukb / mass )
  cv = ov_Ts * ( sum3 / sum1 - ( sum2 / sum1 )**2 )
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
