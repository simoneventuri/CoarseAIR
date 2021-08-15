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

Module LevelsGenerator_Module

  use, intrinsic :: ieee_arithmetic      
  use Parameters_Module       ,only:  rkp, Half, Zero, One, Two, Three, Four, Five, Ten, Hartree_To_eV
  use Logger_Class            ,only:  Logger
  use Error_Class             ,only:  Error
  use Level_Class             ,only:  Level_Type

  implicit none

  private
  public    ::    ComputeEnergyLevels
  public    ::    ReadEnergyLevels

  logical   ::    i_Debug_Global = .False.

  contains
  
!________________________________________________________________________________________________________________________________!
Subroutine ComputeEnergyLevels(Input, Collision, iMol, i_Debug, i_Debug_Deep) 
! This procedure ...
  
  use Collision_Class               ,only:  Collision_Type
  use Level_Class                   ,only:  Level_Type
  use Input_Class                   ,only:  Input_Type

  type(Input_Type)                                       ,intent(in)     ::    Input
  type(Collision_Type)                                   ,intent(in)     ::    Collision
  integer                                                ,intent(in)     ::    iMol
  logical                                  ,optional     ,intent(in)     ::    i_Debug
  logical                                  ,optional     ,intent(in)     ::    i_Debug_Deep

  integer               ,dimension(:)   ,allocatable                     ::    vqnMax
  integer                                                                ::    Status
  type(Level_Type)                                                       ::    State
  logical                                                                ::    vNotEmpty
  integer                                                                ::    ix
  integer                                                                ::    ivqn, ijqn
  real(rkp)                                                              ::    EMin, EMax
  real(rkp)                                                              ::    ELevelCurrent
  real(rkp)                                                              ::    DeltaGridLocal
  integer                                                                ::    NGridLocal
  integer                                                                ::    NNodesMin, NNodesMax, NNodes
  real(rkp)                                                              ::    yLastMin, yLastMax, yLastNew
  real(rkp)                                                              ::    EintMin, EintMax
  real(rkp)                                                              ::    EintOld, EintNew
  real(rkp)                                                              ::    EErr
  real(rkp)                                                              ::    Vc_R2
  integer                                                                ::    ierr
  integer                                                                ::    Unit_Levels, Unit_vqn
  real(rkp)             ,dimension(:) ,allocatable                       ::    x 
  real(rkp)             ,dimension(:) ,allocatable                       ::    EDiat         
  real(rkp)             ,dimension(:) ,allocatable                       ::    y  
  real(rkp)                                                              ::    psiFinal  
  character(:)                                              ,allocatable ::    levels_file
  character(:)                                              ,allocatable ::    vqn_file
  real(rkp)                                                              ::    EintPrevious
  real(rkp)                                                              ::    GridLeft0
  real(rkp)                                                              ::    xmu, xmui, xmui2
  integer                                                                ::    iE
  real(rkp)                                                              ::    Em, E, Ep 
  real(rkp) ,dimension(2)                                                ::    EStart 
  real(rkp) ,dimension(500)                                              ::    Err
  integer   ,dimension(500)                                              ::    NChngsSign
  real(rkp)                                                 ,parameter   ::    MinValue = 9.999999e-99_rkp
  logical                                                                ::    i_Debug_Loc 
  logical                                                                ::    i_Debug_Loc_Deep

  integer   ,parameter                                                   ::    NBisect   =   40                 ! number of bisection steps to perform
  integer   ,parameter                                                   ::    NNewton   =   80                 ! number of newton-raphson steps to perform
  
  i_Debug_Loc      = i_Debug_Global; if ( present(i_Debug)      ) i_Debug_Loc      = i_Debug
  i_Debug_Loc_Deep = .False.;        if ( present(i_Debug_Deep) ) i_Debug_Loc_Deep = i_Debug_Deep
  if (i_Debug_Loc) call Logger%Entering( "ComputeEnergyLevels")
  !i_Debug_Loc   =     Logger%On()
  

  allocate( vqnMax(0:Input%ijqnMax), stat=Status )
  if (Status/=0) call Error( "Error allocating vqnMax" )
  vqnMax = 0
  if (i_Debug_Loc) call Logger%Write( "Allocated vqnMax with dimension (", Input%ijqnMax, "+1); Initialized to 0" )


  xmu   = Collision%MoleculesContainer(iMol)%Molecule%xmu
  xmui  = Collision%MoleculesContainer(iMol)%Molecule%xmui
  xmui2 = Collision%MoleculesContainer(iMol)%Molecule%xmui2
  if (i_Debug_Loc) call Logger%Write( "xmu   = ", xmu   )
  if (i_Debug_Loc) call Logger%Write( "xmui  = ", xmui  )
  if (i_Debug_Loc) call Logger%Write( "xmui2 = ", xmui2 )


  ! ! ==============================================================================================================
  ! !   OPENING FILE FOR WRITING PROCESSOR's ENERGY LIST
  ! ! ==============================================================================================================
  levels_file  = Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // '/Temp_' // trim(adjustl(Input%iNode_char)) // '_' // trim(adjustl(Input%iProc_char)) // '.inp'
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", levels_file )
  open( File=levels_file, NewUnit=Unit_Levels, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // levels_file ) 
  ! ! ==============================================================================================================
  
  
  ! ! ==============================================================================================================
  ! !   OPENING FILE FOR WRITING NUMBER OF vqns PER jqn
  ! ! ==============================================================================================================
  vqn_file  = Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // '/vqn_' // trim(adjustl(Input%iNode_char)) // '_' // trim(adjustl(Input%iProc_char)) // '.inp'
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", vqn_file )
  open( File=vqn_file, NewUnit=Unit_vqn, status='REPLACE', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // vqn_file ) 
  ! ! ==============================================================================================================
  

  ! ! ==============================================================================================================
  ! !   LOOP ON J QUANTUM NUMBER
  ! ! ==============================================================================================================
  vNotEmpty = .True.
  ijqn      = Input%ijqnMin
  do while ( (ijqn <= Input%ijqnMax) .and. (vNotEmpty) ) 
  EMin = Input%EStart
  
    if ( mod(ijqn, int(min(Input%NProc * Input%NNode,200))) == ((Input%iNode-1)*Input%NProc + Input%iProc - 1) ) then
  
 
      ! ! ==============================================================================================================
      ! !   COMPUTING DIATOMIC POTENTIAL MINIMAX POINTS
      ! ! ============================================================================================================== 
      Vc_R2 = xmui2 * ( ijqn + Half )**2  
      call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%FindMinimum( [Input%xExtremes(1), Input%xExtremes(2) / Two], Vc_R2, Input%EEpsilon, State%rmin, State%Vmin, ierr)

      call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%FindMaximum( [State%rmin, Input%xExtremes(2)], Vc_R2, Input%EEpsilon, State%rmax, State%Vmax, ierr )
      ! ! ==============================================================================================================
      

      ! ! ==============================================================================================================
      ! !   CONSTRUCTING x GRID
      ! ! ==============================================================================================================
      call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%TurningPoint( [1.e-1_rkp, State%rmin], Vc_R2, State%Vmax, GridLeft0, ierr, NBisection=NBisect, NNewton=NNewton )
      GridLeft0 = GridLeft0 - 0.5d0

      DeltaGridLocal = max( (State%rmax - GridLeft0) / (Input%NGrid-1), Input%DeltaGrid)
      NGridLocal     = int( (State%rmax - GridLeft0) / DeltaGridLocal ) + 1
      
      allocate(x(NGridLocal))
      allocate(EDiat(NGridLocal))
      allocate(y(NGridLocal))

      
      if (Input%BoundayConditions(2,3) <= 0.d0) then
        psiFinal = real(NGridLocal,rkp)
      else
        psiFinal = Input%BoundayConditions(2,3)
      end if
      
      x    = Zero
      x(1) = GridLeft0
      do ix = 2,NGridLocal
        x(ix) = x(ix-1) + DeltaGridLocal
      end do
      ! ! ==============================================================================================================
      
      
      ! ! ==============================================================================================================
      ! !   COMPUTING DIATOMIC POTENTIAL
      ! ! ============================================================================================================== 
      do ix = 1,NGridLocal
        EDiat(ix) = min(Collision%MoleculesContainer(iMol)%Molecule%DiatPot%EffectivePotential( x(ix), Vc_R2 ), 1.e300_rkp)
      end do
      EDiat = EDiat - State%Vmin
      ! ! ============================================================================================================== 


      ! ! ==============================================================================================================
      ! !   LOOP ON V QUANTUM NUMBER
      ! ! ============================================================================================================== 
      vNotEmpty     = .false.
      ivqn          = Input%ivqnMin
      ELevelCurrent = -1.e10_rkp
      EintPrevious  = -1.e100_rkp
      if (i_Debug_Loc) call Logger%Write( "ELevelCurrent = ", ELevelCurrent, "; State%Vmax = ", State%Vmax )
      do while ( (ivqn <= Input%ivqnMax) .and. (ELevelCurrent <= State%Vmax))! .and. (ELevelCurrent > EintPrevious + Input%DeltaE)) 
        EintPrevious = ELevelCurrent
        EStart       = Zero
        Err          = Zero
        NChngsSign   = 0 
            
        iE = 1;
        call Numerov(NGridLocal, Input%BoundayConditions, xmui, DeltaGridLocal, EDiat, EMin, y, yLastMin, NNodesMin)    
        Err(iE)        = yLastMin
        NChngsSign(iE) = int(NNodesMin)  
        !write(*,*) EMin + State%Vmin, NNodesMin, yLastMin


        iE        = iE + 1
        EMax      = EMin
        NNodesMax = NNodesMin
        do while (NNodesMin == NNodesMax)
          if (State%Vmax - ELevelCurrent > 1.d-2) then
            EMax = EMax + Input%DeltaE
          elseif (State%Vmax - ELevelCurrent > 1.d-3) then
            EMax = EMax + Input%DeltaE / 1.d1;
          elseif (State%Vmax - ELevelCurrent > 1.d-4) then
            EMax = EMax + Input%DeltaE / 1.d2;
          else
            EMax = EMax + Input%DeltaE / 1.d3;
          end if
          call Numerov(NGridLocal, Input%BoundayConditions, xmui, DeltaGridLocal, EDiat, EMax, y, yLastMax, NNodesMax)    
          Err(iE)        = yLastMax
          NChngsSign(iE) = int(NNodesMax)    
          !write(*,*) EMax + State%Vmin, NNodesMax, yLastMax
        end do
        !pause

        ! ! ==============================================================================================================
        ! !   BISECTION METHOD
        ! ! ============================================================================================================== 
        EStart  = [EMin, EMax]
        if ( product(Err) > Zero ) call Error( "ERRROR! Initial Guesses for Energies Bring to Psi(end) with Same Signs!" ) 
        if (yLastMin > 0.0) then
          Ep = EStart(1)
          Em = EStart(2)
        else
          Ep = EStart(2)
          Em = EStart(1)
        end if
        do while ( abs(Ep - Em) > Input%EEpsilon )
          iE      = iE+1

          E       = (Ep + Em) / Two
          call Numerov(NGridLocal, Input%BoundayConditions, xmui, DeltaGridLocal, EDiat, E, y, yLastNew, NNodes) 
          Err(iE)        = yLastNew;
          NChngsSign(iE) = int(NNodes);

          if (Err(iE) > Zero) then
            Ep = E
          else
            Em = E
          end if

        end do
        EintOld = (Ep + Em) / Two
        ! ! ============================================================================================================== 

        
        if (((NNodes == NNodesMin) .or. (NNodes == NNodesMax)) .and. (EintOld + State%Vmin <= State%Vmax)) then
          vNotEmpty     = .true.
          State%eint    = EintOld + State%Vmin
          ELevelCurrent = State%eint
          State%vqn     = ivqn
          State%jqn     = ijqn


          ! ! ==============================================================================================================
          ! !   COMPUTING TURNING POINTS
          ! ! ============================================================================================================== 
          call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%TurningPoint( [1.e-1_rkp, State%rmin], Vc_R2, State%eint, State%ri, ierr, NBisection=NBisect, NNewton=NNewton )
          ! ! ============================================================================================================== 
                  
          
          ! ! ==============================================================================================================
          ! !   COMPUTING TAU
          ! ! ============================================================================================================== 
          call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%Period( State%Eint, Vc_R2, .True., State%rMin, State%VMin, State%rMax, State%VMax, State%ri, State%ro, State%Tau, i_Debug_Loc_Deep )
          ! ! ============================================================================================================== 
          
          
          ! ! ==============================================================================================================
          ! !   COMPUTING RESONANCE WIDTH
          ! ! ============================================================================================================== 
          if (State%eint < 0.e0_rkp) then
            State%Egam = 0.e0_rkp
          else 
            call Collision%MoleculesContainer(iMol)%Molecule%DiatPot%ResonanceWidth( State%Eint, Vc_R2, State%rMin, State%rMax, State%Tau, State%Egam, i_Debug_Loc_Deep )
            !if (abs(State%egam) >= 1.e300_rkp) then
            if ( ieee_is_nan(abs(State%egam)) .or. (.not.ieee_is_finite(abs(State%egam))) ) then 
              write(*,'(a,i2,i4)')'WARNING:: State%egam is NAN or INF for (v,J) pair: ',ivqn,iJqn    
              print*,State%Egam 
              State%egam = Zero
            elseif (State%egam < MinValue) then
              State%egam = MinValue
            end if
          end if
          ! ! ============================================================================================================== 
          
          ! ! ==============================================================================================================
          ! !   WRITING ENERGY LEVEL
          ! ! ============================================================================================================== 
          write(Unit_Levels,'(X,2I5,*(es15.7))') State%vqn,  State%jqn, State%eint, State%egam, State%rmin, State%rmax, State%Vmin, State%Vmax, State%tau, State%ri, State%ro
          flush(Unit_Levels) 
          if (Input%PrintLevelsFlg) then
            write(*,'(X,2I5,*(es15.7))') State%vqn,  State%jqn, State%eint, State%egam, State%rmin, State%rmax, State%Vmin, State%Vmax, State%tau, State%ri, State%ro
          end if
          ! ! ============================================================================================================== 
          !pause


          ! ! ==============================================================================================================
          ! !   WRITING ENERGY LEVEL
          ! ! ============================================================================================================== 
          if (Input%WriteWF) then

            ! open( NewUnit=UnitWrite, File='./trajectories.bin', Action='WRITE', access="Stream", form="Unformatted", iostat=StatusWrite )
            ! if (StatusWrite/=0) call Error( "Error writing the binary data file for wave function: " // './trajectories.bin'  ) 
            ! write(UnitWrite) int(This%NTraj, rkp)

            ! if (i_Debug_Loc) call Logger%Write( "Reading the trajectory data: bMax, bSampled, Qini, Qfin" )
            ! rewind(DataFile%Unit)                                                                                           
            ! read(DataFile%Unit,*)                                                                                         
            ! do iTraj = 1,This%NTraj                                                                                       
            !   read(DataFile%Unit,*,iostat=DataFile%Status) Idx, iPES, This%bMax(iTraj), This%bSampled(iTraj), This%Qini(:,iTraj), This%Qfin(:,iTraj)
            !   if (DataFile%Status/=0) call Error( "Error reading the data file for statistics: " // DataFile%Name  )  

            !   if (This%StatWritesBinaryFlg) then
            !     write(UnitWrite) int(Idx,  rkp)
            !     write(UnitWrite) int(iPES, rkp)
            !     write(UnitWrite) This%bMax(iTraj)
            !     write(UnitWrite) This%bSampled(iTraj)
            !     do iCond=1,This%NCond
            !       write(UnitWrite) This%Qini(iCond,iTraj)
            !     end do
            !     do iCond=1,This%NCond
            !       write(UnitWrite) This%Qfin(iCond,iTraj)
            !     end do
            !   end if

            ! end do                                                                                                        
            ! if (i_Debug_Loc) then
            !   call Logger%Write( "-> Done reading the trajectory data" )
            !   call Logger%Write( "-> Last line: iTraj = ", "This%bMax(iTraj) = ", This%bMax(This%NTraj), "This%bSampled(iTraj) = ", This%bSampled(This%NTraj), Fi="i9", Fr="es15.8")
            ! end if
            ! call DataFile%Close()

            !if (This%StatWritesBinaryFlg) close(UnitWrite)


          end if
          ! ! ============================================================================================================== 
          !pause
          
          ivqn = ivqn + 1
        elseif (EintOld + State%Vmin > State%Vmax) then
          ELevelCurrent = State%Vmax + 1.d0
        end if
        EMin = EintOld + Input%EEpsilon*1.e1_rkp


      end do    
      vqnMax(ijqn) = ivqn
      ! ! ============================================================================================================== 
      
      
      deallocate(x)
      deallocate(EDiat)
      deallocate(y)
      
    end if

    
    ! ! ==============================================================================================================
    ! !   WRITING NB of vqns PER jqn
    ! ! ============================================================================================================== 
    write(Unit_vqn,'(I5)') vqnMax(ijqn)
    flush(Unit_vqn)
    ! ! ============================================================================================================== 


    ijqn = ijqn + 1
  end do
  ! ! ==============================================================================================================


  close(Unit_Levels)
  
  close(Unit_vqn)

  deallocate( vqnMax )

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadEnergyLevels( Input, Collision, iMol, i_Debug ) 
! This procedure ...
  
  use Input_Class                   ,only:  Input_Type
  use Collision_Class               ,only:  Collision_Type
  use LevelsContainer_Class         ,only:  LevelsContainer_Type

  type(Input_Type)                                       ,intent(in)     ::    Input
  type(Collision_Type)                                   ,intent(inout)  ::    Collision
  integer                                                ,intent(in)     ::    iMol
  logical                                  ,optional     ,intent(in)     ::    i_Debug

  integer                                                                ::    NProc, NNode
  integer                                                                ::    ijqnMax
  integer                                                                ::    Status
  integer                                                                ::    iState
  character(:)                                              ,allocatable ::    levels_file
  character(:)                                              ,allocatable ::    vqn_file
  integer               ,dimension(:,:)                     ,allocatable ::    vqnMax
  integer               ,dimension(:)                       ,allocatable ::    vqnMaxOverall  
  integer                                                                ::    ijqn, ivqn
  integer                                                                ::    i
  integer                                                                ::    iNode
  character(2)                                                           ::    iNode_char
  integer                                                                ::    iProc
  integer                                                                ::    iRun
  character(2)                                                           ::    iProc_char
  integer                                                                ::    Unit_Levels, Unit_vqn
  integer                                                                ::    NStates
  Type(LevelsContainer_Type)                                             ::    LevelsContainer
  character(150)                                                         ::    FileName
  logical                                                                ::    i_Debug_Loc 

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadEnergyLevels")
  !i_Debug_Loc   =     Logger%On()


  NProc   = Input%NProc
  NNode   = Input%NNode
  ijqnMax = Input%ijqnMax

  allocate( vqnMax(0:ijqnMax,NProc*NNode), stat=Status )
  if (Status/=0) call Error( "Error allocating vqnMax" )
  vqnMax = 0
  if (i_Debug_Loc) call Logger%Write( "Allocated vqnMax with dimension (", ijqnMax, "+1,", NProc*NNode, "); Initialized to 0" )


  allocate( vqnMaxOverall(-1:ijqnMax), stat=Status )
  if (Status/=0) call Error( "Error allocating vqnMaxOverall" )
  vqnMaxOverall = 0
  if (i_Debug_Loc) call Logger%Write( "Allocated vqnMaxOverall with dimension (", ijqnMax, "+2); Initialized to 0" )  


  iRun = 1
  do iNode = 1,NNode
    write( iNode_char, '(I2)' ) iNode
  
    do iProc = 1,NProc
      write( iProc_char, '(I2)' ) iProc
      
      if (iRun <= 200) then

        vqn_file  = Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // '/vqn_' // trim(adjustl(iNode_char)) // '_' // trim(adjustl(iProc_char)) // '.inp'
        if (i_Debug_Loc) call Logger%Write( "Opening File vqn_file", vqn_file )
        open(NewUnit=Unit_vqn, File=vqn_file, status='OLD', iostat=Status)
        if (Status/=0) call Error( "Impossible to open vqn_file" )
        
          ijqn = Input%ijqnMin
          do while ((Status == 0) .and. (ijqn <= ijqnMax))
          
            read(Unit_vqn,*,iostat=Status) vqnMax(ijqn,iRun) 
          
            ijqn = ijqn + 1
          end do
          
        close(Unit_vqn)      
        
        iRun = iRun + 1
      end if   
    
    end do
    
  end do


  vqnMaxOverall(0) = sum(vqnMax(0,:))
  do ijqn = 1,ijqnMax
    vqnMaxOverall(ijqn) = vqnMaxOverall(ijqn-1) + sum(vqnMax(ijqn,:))
  end do
  NStates = vqnMaxOverall(ijqnMax)


  ! ==============================================================================================================
  !   7.4. CONSTRUCTING THE LEVEL CONTAINER
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Calling LevelsContainer%Initialize" )
  call LevelsContainer%Initialize( Input, Collision%MoleculesContainer(iMol)%Molecule%DiatPot, iMol, 0, NStates=NStates, ReCheckFlg=.True., i_Debug=i_Debug_Loc )
  if (i_Debug_Loc) call Logger%Write( "-> Done with LevelsContainer%Initialize" )
  ! ==============================================================================================================  

  
  iRun = 1
  do iNode = 1,NNode
    write( iNode_char, '(I2)' ) iNode
    
    do iProc = 1,NProc
      write( iProc_char, '(I2)' ) iProc
      
      if (iRun <= 200) then
    
        levels_file  = Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // '/Temp_' // trim(adjustl(iNode_char)) // '_' // trim(adjustl(iProc_char)) // '.inp'
        if (i_Debug_Loc) call Logger%Write( "Reading Levels from ", levels_file )
        open( NewUnit=Unit_Levels, File=levels_file, status='OLD', iostat=Status)
      
          ijqn = Input%ijqnMin
          do while ((Status == 0) .and. (ijqn <= ijqnMax))
          
            do ivqn=0,vqnMax(ijqn,iRun)-1
              read(Unit_Levels,'(X,2I5,*(es15.7))',iostat=Status) LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%vqn,   & 
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%jqn,   &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%eint,  &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%egam,  &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%rmin,  &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%rmax,  &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%Vmin,  &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%Vmax,  &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%tau,   &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%ri,    &
                                                                  LevelsContainer%States(vqnMaxOverall(ijqn-1)+ivqn+1)%ro
            end do
          
            ijqn = ijqn + 1
          end do
          
        close(Unit_Levels)   
        
        iRun = iRun + 1
      end if
    
    end do
    
  end do


  FileName = Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // trim(adjustl(Input%GeneratedLevelsFile(iMol)))
  if (i_Debug_Loc) call Logger%Write( "Writing Levels in ", FileName )
  call LevelsContainer%WriteList( FileName, SortLevelsFlg=Input%SortLevelsFlg, i_Debug=i_Debug_Loc ) 

  
  call system('rm -rf ' // Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // '/Temp_*' )
                   
  call system('rm -rf ' // Collision%MoleculesContainer(iMol)%Molecule%PathToMolFldr // '/vqn_*'  )

  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine Numerov(NGrid, psi, xmui, DeltaGrid, EDiat, E, y, yLast, NNodes) 
! This procedure ...

  integer                                                ,intent(in)     ::    NGrid
  real(rkp)             ,dimension(2,3)                  ,intent(in)     ::    psi
  real(rkp)                                              ,intent(in)     ::    xmui
  real(rkp)                                              ,intent(in)     ::    DeltaGrid
  real(rkp)             ,dimension(NGrid)                ,intent(in)     ::    EDiat
  real(rkp)                                              ,intent(in)     ::    E
  real(rkp)             ,dimension(NGrid)                ,intent(out)    ::    y
  real(rkp)                                              ,intent(out)    ::    yLast    
  integer                                                ,intent(out)    ::    NNodes
  
  real(rkp)             ,dimension(NGrid)                                ::    g
  integer                                                                ::    ix
  integer                                                                ::    iBC

  NNodes           = 0
  g                = Two / xmui * (EDiat - E) 
  iBC              = 1

  do ix = 1,NGrid
  
    if (ix == int(psi(2,iBC))) then
    
      y(int(psi(2,iBC))) = psi(1,iBC)
      iBC = iBC + 1
      
    else 

      y(ix) = (DeltaGrid**2 / 12.e0_rkp * ( Ten*y(ix-1)*g(ix-1) + y(ix-2)*g(ix-2) ) + Two*y(ix-1) - y(ix-2) ) / ( One - DeltaGrid**2 * g(ix) / 12.e0_rkp )
      
      if ( (y(ix)*y(ix-1) <= Zero) .and. (y(ix-1) /= 0) ) then
        NNodes = NNodes + 1
      end if
      
    end if
    
  end do

  yLast = y(size(y,1))
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
