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

Module KonigPostprocessing_Class

  use Parameters_Module          ,only:  rkp, Zero, Two
  use Logger_Class               ,only:  Logger, LogLevel_INFO, LogLevel_DEBUG
  use Error_Class                ,only:  Error
  use KonigPostprocessingT_Class ,only:  KonigPostprocessingT_Type

  implicit none

  !private
  public    ::    KonigPostprocessing_Type

  Type      ::    KonigPostprocessing_Type
    integer                                                              ::    NComponents
    character(5)                     ,dimension(:)   ,allocatable        ::    ComponentsName
    type(KonigPostprocessingT_Type)  ,dimension(:)   ,allocatable        ::    Ttra

  contains
    private
    
    procedure ,public   ::    ReadComponents
    procedure ,public   ::    ReadBox
    procedure ,public   ::    ReadPop
    procedure ,public   ::    ReadTemperature
    procedure ,public   ::    DeallocateKonigPost
    
    procedure ,public   ::    EvenlySpacedVector
    procedure ,public   ::    AllocateOutputToAnalyze
    procedure ,public   ::    Allocate2DOutputToAnalyze
    procedure ,public   ::    LinearInterpAtAbscissaGrid
    procedure ,public   ::    Histogram
    procedure ,public   ::    WriteTimeGrid
    procedure ,public   ::    WriteOutputSums
    procedure ,public   ::    WriteOutputHist
    
  End Type

  logical   ,parameter    ::    i_Debug_Global = .True.

  contains


Pure Subroutine EvenlySpacedVector(This, VectorMin, VectorMax, NIntervals, VectorScale, Vector)
    
  class(KonigPostprocessing_Type)      ,intent(in)  ::    This
  real(rkp)                            ,intent(in)  ::    VectorMin
  real(rkp)                            ,intent(in)  ::    VectorMax
  integer                              ,intent(in)  ::    NIntervals
  character(3)                         ,intent(in)  ::    VectorScale
  real(rkp) ,dimension(:) ,allocatable ,intent(out) ::    Vector
  
  real(rkp)                                         ::    h
  integer                                           ::    iIntervals

  allocate(Vector(NIntervals+1) )
  
  if ( (VectorScale == 'lin') .or. (VectorScale == 'LIN') ) then
  
    h = dabs(VectorMax - VectorMin) / real(NIntervals)
  
    Vector(1) = VectorMin
    do iIntervals = 2,NIntervals+1
       Vector(iIntervals) = Vector(1) + (iIntervals-1) * h
    end do
  
  elseif ( (VectorScale == 'log') .or. (VectorScale == 'LOG') ) then
    
    h = dabs((dlog10(VectorMax) - dlog10(VectorMin))) / real(NIntervals)
  
    Vector(1) = VectorMin
    do iIntervals = 2,NIntervals+1
       Vector(iIntervals) = 10.0_rkp**(dlog10(Vector(1)) + (iIntervals-1) * h)
    end do

  else
  
    !call Error( "Vector Scale NOT Specified. Please, choose lin / log." )

  end if
  
End Subroutine


Subroutine AllocateOutputToAnalyze( This, Input, NIntervals, NBins, OutputAtNodes, OutputSum, OutputSqSum, OutputHist, i_Debug )
    
  use Input_Class ,only:  Input_Type
    
  class(KonigPostprocessing_Type)        ,intent(in)  ::    This
  type(Input_Type)                       ,intent(in)  ::    Input
  integer                                ,intent(in)  ::    NIntervals
  integer                                ,intent(in)  ::    NBins
  real(rkp) ,dimension(:,:) ,allocatable ,intent(out) ::    OutputAtNodes
  real(rkp) ,dimension(:)   ,allocatable ,intent(out) ::    OutputSum
  real(rkp) ,dimension(:)   ,allocatable ,intent(out) ::    OutputSqSum
  integer   ,dimension(:,:) ,allocatable ,intent(out) ::    OutputHist
  logical                      ,optional ,intent(in)  ::    i_Debug
  
  integer                                             ::    Status
  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AllocateOutputToAnalyze")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  

  allocate(OutputAtNodes(NIntervals+1,Input%NFwdProp), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputAtNodes" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputAtNodes with dimension = (", NIntervals+1, ',', Input%NFwdProp, ")" )
  OutputAtNodes = Zero
  
  allocate(OutputSum(NIntervals+1), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputSum" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputSum with dimension = (", NIntervals+1, ")" )
  OutputSum = Zero
  
  allocate(OutputSqSum(NIntervals+1), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputSqSum" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputSqSum with dimension = (", NIntervals+1, ")" )
  OutputSqSum = Zero
  
  allocate(OutputHist(NBins, NIntervals+1), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputHist" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputHist with dimension = (", NBins, ",", NIntervals+1, ")" )
  OutputHist = 0
  
  
  call Logger%Exiting
  
End Subroutine


Subroutine Allocate2DOutputToAnalyze( This, SecondDim, Input, NIntervals, NBins, OutputAtNodes, OutputSum, OutputSqSum, OutputHist, i_Debug )
    
  use Input_Class ,only:  Input_Type
    
  class(KonigPostprocessing_Type)          ,intent(in)  ::    This
  integer                                  ,intent(in)  ::    SecondDim
  type(Input_Type)                         ,intent(in)  ::    Input
  integer                                  ,intent(in)  ::    NIntervals
  integer                                  ,intent(in)  ::    NBins
  real(rkp) ,dimension(:,:,:) ,allocatable ,intent(out) ::    OutputAtNodes
  real(rkp) ,dimension(:,:)   ,allocatable ,intent(out) ::    OutputSum
  real(rkp) ,dimension(:,:)   ,allocatable ,intent(out) ::    OutputSqSum
  integer   ,dimension(:,:,:) ,allocatable ,intent(out) ::    OutputHist
  logical                        ,optional ,intent(in)  ::    i_Debug
  
  integer                                               ::    Status
  logical                                               ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "AllocateOutputToAnalyze")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()

  allocate(OutputAtNodes(NIntervals+1,SecondDim,Input%NFwdProp), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputAtNodes" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputAtNodes with dimension = (", NIntervals+1, ',', SecondDim, ',', Input%NFwdProp, ")" )
  OutputAtNodes = Zero
  
  allocate(OutputSum(NIntervals+1,SecondDim), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputSum" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputSum with dimension = (", NIntervals+1, ',', SecondDim, ")" )
  OutputSum = Zero
  
  allocate(OutputSqSum(NIntervals+1,SecondDim), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputSqSum" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputSqSum with dimension = (", NIntervals+1, ',', SecondDim, ")" )
  OutputSqSum = Zero
  
  allocate(OutputHist(NBins, SecondDim, NIntervals+1), Stat=Status )
  if (Status/=0) call Error( "Error allocating OutputHist" )
  if (i_Debug_Loc) call Logger%Write( "Allocated OutputHist with dimension = (", NBins, ',', SecondDim, ',', NIntervals+1, ")" )
  OutputHist = 0
  
  
  call Logger%Exiting
  
End Subroutine


Subroutine ReadComponents( This, OutputPath, i_Debug )

  class(KonigPostprocessing_Type)           ,intent(inout)  ::    This
  character(:)                 ,allocatable ,intent(in)     ::    OutputPath
  logical                         ,optional ,intent(in)     ::    i_Debug

  character(:)     ,allocatable                             ::    FileName
  integer                                                   ::    Status
  integer                                                   ::    Unit
  integer                                                   ::    iComp
  integer                                                   ::    temp
  character(100)                                            ::    line
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadComponents")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()
  
  
  FileName = trim(adjustl(OutputPath)) // 'components'
  if (i_Debug_Loc) call Logger%Write( "Reading File: ", FileName )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  
    read(Unit,'(A100)') line
    read(line(34:100),'(I10)') This%NComponents
    if (i_Debug_Loc) call Logger%Write( "There are ", This%NComponents, " Components in the System" )
    
    allocate(This%ComponentsName(This%NComponents), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%ComponentsName" )
    if (i_Debug_Loc) call Logger%Write( "Allocated This%ComponentsName with dimension = (", This%NComponents, ")" )
    
    do iComp=1,This%NComponents
      read(Unit,'(I2,A)') temp, This%ComponentsName(iComp)
      if (i_Debug_Loc) call Logger%Write( "Component Nb", iComp, " = ", This%ComponentsName(iComp) )
    end do
  
  close(Unit)
  
  
  call Logger%Exiting

End Subroutine


Subroutine ReadBox( This, OutputPath, Input, iTtra, i_Debug )

  use Input_Class ,only:  Input_Type

  class(KonigPostprocessing_Type)           ,intent(inout)  ::    This
  character(:)                 ,allocatable ,intent(in)     ::    OutputPath
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iTtra
  logical                         ,optional ,intent(in)     ::    i_Debug

  character(:)                              ,allocatable    ::    FileName
  integer                                                   ::    Status
  integer                                                   ::    Unit
  integer                                                   ::    iLines
  integer                                                   ::    iComp
  integer                                                   ::    iTemp
  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadBox")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()


  FileName = trim(adjustl(OutputPath)) // 'box.dat'
  
  
  if (i_Debug_Loc) call Logger%Write( "Reading File: ", FileName )
  open( File=FileName, Unit=127, status='OLD', action='READ', iostat=Status )
  read(127,*,iostat=Status)
  iLines = 0
  do while (Status==0) 
    read(127,*,iostat=Status)
    iLines = iLines + 1
  end do
  This%Ttra(iTtra)%NTimeSteps = iLines-1
  
  allocate(This%Ttra(iTtra)%Time(This%Ttra(iTtra)%NTimeSteps), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%Time" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%Time with dimension = (", This%Ttra(iTtra)%NTimeSteps , ")" )
  
  allocate(This%Ttra(iTtra)%X(This%Ttra(iTtra)%NTimeSteps, This%NComponents), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%X" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%X with dimension = (", This%Ttra(iTtra)%NTimeSteps, ',', This%NComponents, ")" )
  
  allocate(This%Ttra(iTtra)%TTranslational(This%Ttra(iTtra)%NTimeSteps), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%TTranslational" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%TTranslational with dimension = (", This%Ttra(iTtra)%NTimeSteps, ")" )
  
  allocate(This%Ttra(iTtra)%Rho(This%Ttra(iTtra)%NTimeSteps), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%Rho" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%Rho with dimension = (", This%Ttra(iTtra)%NTimeSteps, ")" )
  
  allocate(This%Ttra(iTtra)%P(This%Ttra(iTtra)%NTimeSteps), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%P" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%P with dimension = (", This%Ttra(iTtra)%NTimeSteps, ")" )
  
  allocate(This%Ttra(iTtra)%Nd(This%Ttra(iTtra)%NTimeSteps), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%Nd" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%Nd with dimension = (", This%Ttra(iTtra)%NTimeSteps, ")" )
  
  allocate(This%Ttra(iTtra)%Energy(This%Ttra(iTtra)%NTimeSteps), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Ttra(iTtra)%Energy" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Ttra(iTtra)%Energy with dimension = (", This%Ttra(iTtra)%NTimeSteps, ")" )

  rewind(127)
  read(127,*)
  do iLines = 1,This%Ttra(iTtra)%NTimeSteps
    read(127,'(*(d20.10))',iostat=Status) This%Ttra(iTtra)%Time(iLines), (This%Ttra(iTtra)%X(iLines,iComp),iComp=1,This%NComponents), This%Ttra(iTtra)%TTranslational(iLines), &
                              This%Ttra(iTtra)%Rho(iLines), This%Ttra(iTtra)%P(iLines), This%Ttra(iTtra)%Nd(iLines), This%Ttra(iTtra)%Energy(iLines)
  end do
  close(127)


  call Logger%Exiting

End Subroutine


Subroutine ReadPop( This, OutputPath, Input, iTtra, i_Debug )

  use Input_Class         ,only:  Input_Type

  class(KonigPostprocessing_Type)           ,intent(inout)  ::    This
  character(:)                 ,allocatable ,intent(in)     ::    OutputPath
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iTtra
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  character(:)     ,allocatable                             ::    FileName
  integer                                                   ::    Status
  integer                                                   ::    Unit
  integer                                                   ::    iMol
  integer                                                   ::    iLines
  integer                                                   ::    iBins
  real(rkp)                                                 ::    temp
  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadPop")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()


  iMol = 1
  allocate(This%Ttra(iTtra)%Pop(This%Ttra(iTtra)%NTimeSteps, Input%NBins(iMol) ), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Pop" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Pop with dimension = (", This%Ttra(iTtra)%NTimeSteps, ',', Input%NBins(iMol), ")" )


  do iMol = 1,1
  
    FileName = trim(adjustl(OutputPath)) // 'pop_' // trim(adjustl(Input%Molecules_Name(iMol))) // '.dat'
    if (i_Debug_Loc) call Logger%Write( "Reading File: ", FileName )
    open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
    read(Unit,*)
    read(Unit,*)
    do iLines = 1,This%Ttra(iTtra)%NTimeSteps
      do iBins = 1,Input%NBins(iMol)
        read(Unit,*) temp, This%Ttra(iTtra)%Pop(iLines,iBins)
        if (Status/=0) call Error( "Error Reading This%Ttra(iTtra)%Pop." )
      end do
      read(Unit,*)
    end do
    close(Unit)
    
  end do
  
  call Logger%Exiting
  
End Subroutine

  
Subroutine ReadTemperature( This, OutputPath, Input, iTtra, i_Debug )

  use Input_Class         ,only:  Input_Type

  class(KonigPostprocessing_Type)           ,intent(inout)  ::    This
  character(:)                 ,allocatable ,intent(in)     ::    OutputPath
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iTtra
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  character(:)     ,allocatable                             ::    FileName
  integer                                                   ::    Status
  integer                                                   ::    Unit
  integer                                                   ::    iMol
  integer                                                   ::    iLines
  integer                                                   ::    iBins
  real(rkp)                                                 ::    temp
  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadTemperature")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()
  

  allocate(This%Ttra(iTtra)%Tint(This%Ttra(iTtra)%NTimeSteps, Input%NMolecules), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%Tint" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%Tint with dimension = (", This%Ttra(iTtra)%NTimeSteps, "," , Input%NMolecules, " )" )

  FileName = trim(adjustl(OutputPath)) // '/Tint.dat'
  
  if (i_Debug_Loc) call Logger%Write( "Reading File: ", FileName )
  open( File=FileName, NewUnit=Unit, status='OLD', iostat=Status )
  read(Unit,*)
  read(Unit,*)
  do iLines = 1,This%Ttra(iTtra)%NTimeSteps
    read(Unit,*) temp, This%Ttra(iTtra)%Tint(iLines,:)
  end do
  close(Unit)


  call Logger%Exiting

End Subroutine
 

Pure Subroutine LinearInterpAtAbscissaGrid(This, x, y, xNodes, yNodes)
         
  class(KonigPostprocessing_Type)               ,intent(in)    ::    This
  real(rkp)                       ,dimension(:) ,intent(in)    :: x
  real(rkp)                       ,dimension(:) ,intent(in)    :: y
  real(rkp)                       ,dimension(:) ,intent(in)    :: xNodes
  real(rkp)                       ,dimension(:) ,intent(inout) :: yNodes
  
  integer                                                      :: iNodes
  integer                                                      :: ix

  do iNodes = 1,size(xNodes,1)
  
    ix=1
    do 
      if ( (x(ix) >= xNodes(iNodes)) .or. (ix>=size(x,1)) ) exit
      ix = ix + 1
    end do
     
    if (x(ix) >= xNodes(iNodes)) then
      yNodes(iNodes) = y(ix)
    else
      yNodes(iNodes) = (y(ix)-y(ix-1))/(x(ix)-x(ix-1)) * (xNodes(iNodes)-x(ix-1)) + y(ix-1)
    end if

  end do
 
End Subroutine

 
Pure Subroutine Histogram(This, BinsExtremesVec, Vector, HistVector)
         
  class(KonigPostprocessing_Type)  ,intent(in)    ::    This
  real(rkp)  ,dimension(:), intent(in)    :: BinsExtremesVec
  real(rkp)  ,dimension(:), intent(in)    :: Vector
  integer    ,dimension(:), intent(inout) :: HistVector
    
  integer                                 :: i, j, j_target, j_target_max, j_target_min
  logical                                 :: FLAG_def
    
  do i=1,size(Vector,1)
    
    if (Vector(i) < BinsExtremesVec(2)) then

      HistVector(1) = HistVector(1) + 1
      
    elseif (Vector(i) >= BinsExtremesVec(size(BinsExtremesVec,1)-1)) then
   
      HistVector(size(HistVector,1)) = HistVector(size(HistVector,1)) + 1

    else
      
      j_target_min = 2
      j_target_max = size(BinsExtremesVec,1)-1
      j_target = floor(real(j_target_min+j_target_max)/Two)
    
      FLAG_def = .FALSE.
      do

        if (Vector(i) < BinsExtremesVec(j_target)) then

          j_target_min = j_target_min
          j_target_max = j_target
          j_target     = floor(real(j_target_min+j_target_max)/Two)
          
        elseif (Vector(i) >= BinsExtremesVec(j_target+1)) then
       
          j_target_min = j_target+1
          j_target_max = j_target_max
          j_target     = floor(real(j_target_min+j_target_max)/Two)
          
        elseif ((Vector(i) >= BinsExtremesVec(j_target)) .and. (Vector(i) < BinsExtremesVec(j_target+1))) then
       
          HistVector(j_target) = HistVector(j_target) + 1
          FLAG_def       = .TRUE.
          
        end if                     
        if (FLAG_def) exit

      end do

    end if      
 
  end do

End Subroutine


Subroutine DeallocateKonigPost( This, Input, iTtra, i_Debug )

  use Input_Class         ,only:  Input_Type
  
  class(KonigPostprocessing_Type)           ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iTtra
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) ) i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "DeallocateKonigPost")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()
  
  deallocate(This%Ttra(iTtra)%Time, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%Time" )
  
  deallocate(This%Ttra(iTtra)%X, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%X" )
  
  deallocate(This%Ttra(iTtra)%TTranslational, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%TTranslational" )

  deallocate(This%Ttra(iTtra)%Rho, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%Rho" )
  
  deallocate(This%Ttra(iTtra)%P, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%P" )
  
  deallocate(This%Ttra(iTtra)%Nd, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%Nd" )
  
  deallocate(This%Ttra(iTtra)%Energy, Stat=Status )
  if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%Energy" )
    

  if ( trim(adjustl(Input%FWD_Temperatures)) == 'yes') then
  
    deallocate(This%Ttra(iTtra)%Tint, Stat=Status )
    if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%Tint" )
    
  end if

  if (trim(adjustl(Input%FWD_Populations)) == 'yes') then
  
    deallocate(This%Ttra(iTtra)%Pop, Stat=Status )
    if (Status/=0) call Error( "Error deallocating This%Ttra(iTtra)%Pop" )
    
  end if
  
  call Logger%Exiting

End Subroutine


Subroutine WriteTimeGrid( This, OutputPath, FileName, Input, Grid, i_Debug )

  use Input_Class         ,only:  Input_Type
  
  class(KonigPostprocessing_Type)           ,intent(in)     ::    This
  character(:) ,allocatable                 ,intent(in)     ::    OutputPath
  character(:) ,allocatable                 ,intent(in)     ::    FileName
  type(Input_Type)                          ,intent(in)     ::    Input
  real(rkp), dimension(:)                   ,intent(in)     ::    Grid
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  character(:) ,allocatable                                 ::    FileNameFinal
  integer                                                   ::    Unit
  integer                                                   ::    iGrid
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteTimeGrid")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()
  
  
  FileNameFinal = trim(adjustl(OutputPath)) // trim(adjustl(FileName))
  if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileNameFinal )
  open( File=FileNameFinal, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,*) '# NTimeNodes = ', Input%NTimeNodes
  write(Unit,*) '#         Time  Nodes       '
  do iGrid = 1,size(Grid,1)
    write(Unit,'(d20.10,3X)') Grid(iGrid)
  end do
  close(Unit)


  call Logger%Exiting

End Subroutine


Subroutine WriteOutputSums( This, OutputPath, FileName, Input, Grid, OutputSum, OutputSqSum, i_Debug )

  use Input_Class         ,only:  Input_Type
  
  class(KonigPostprocessing_Type)           ,intent(in)     ::    This
  character(:) ,allocatable                 ,intent(in)     ::    OutputPath
  character(:) ,allocatable                 ,intent(in)     ::    FileName
  type(Input_Type)                          ,intent(in)     ::    Input
  real(rkp), dimension(:)                   ,intent(in)     ::    Grid
  real(rkp), dimension(:)                   ,intent(in)     ::    OutputSum
  real(rkp), dimension(:)                   ,intent(in)     ::    OutputSqSum
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  character(:) ,allocatable                                 ::    FileNameFinal
  integer                                                   ::    Unit
  integer                                                   ::    iGrid
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteOutputSums")  !, Active = i_Debug_Loc )
  ! i_Debug_Loc   =     Logger%On()
  
  
  FileNameFinal = trim(adjustl(OutputPath)) // trim(adjustl(FileName))
  if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileNameFinal )
  open( File=FileNameFinal, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,*) '# NFwdProp = ', Input%NFwdProp
  write(Unit,*) '#              Time            Output Sums        Output Sq. Sums'
  do iGrid = 1,size(Grid,1)
    write(Unit,'(3(d20.10,3X))') Grid(iGrid), OutputSum(iGrid), OutputSqSum(iGrid)
  end do
  close(Unit)


  call Logger%Exiting

End Subroutine


Subroutine WriteOutputHist( This, OutputPath, FileName, Input, BinsExtremes, Hist, i_Debug )

  use Input_Class         ,only:  Input_Type
  
  class(KonigPostprocessing_Type)           ,intent(in)     ::    This
  character(:) ,allocatable                 ,intent(in)     ::    OutputPath
  character(:) ,allocatable                 ,intent(in)     ::    FileName
  type(Input_Type)                          ,intent(in)     ::    Input
  real(rkp), dimension(:)                   ,intent(in)     ::    BinsExtremes
  integer,   dimension(:,:)                 ,intent(in)     ::    Hist
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  character(:) ,allocatable                                 ::    FileNameFinal
  integer                                                   ::    Unit
  integer                                                   ::    iBins
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "WriteOutputHist")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  FileNameFinal = trim(adjustl(OutputPath)) // trim(adjustl(FileName))
  if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileNameFinal )
  open( File=FileNameFinal, NewUnit=Unit, status='REPLACE', iostat=Status )
  write(Unit,*) '# NFwdProp = ', Input%NFwdProp
  write(Unit,*) '#    Bin Inf Extreme           Nb of Values at Abscissa Location'
  do iBins = 1,size(Hist,1)
    write(Unit,'(d20.10,*(3X,I10))') BinsExtremes(iBins), Hist(iBins,:)
  end do
  close(Unit)


  call Logger%Exiting

End Subroutine


End Module
