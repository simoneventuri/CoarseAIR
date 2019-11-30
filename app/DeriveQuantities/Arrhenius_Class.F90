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

Module Arrhenius_Class

  use Parameters_Module     ,only:  rkp, Zero
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error, CheckVariable

  implicit none

  private
  public    ::    Arrhenius_Type

  Type    ::    Arrhenius_Type
  
    character(20) ,dimension(:) ,allocatable ::    ReactionIn
    character(20) ,dimension(:) ,allocatable ::    ReactionFn
    integer       ,dimension(:) ,allocatable ::    TypeVec
    integer       ,dimension(:) ,allocatable ::    ProcToPair
    integer       ,dimension(:), allocatable ::    MolOK
    
  contains
    private
    procedure ,public   ::      CreatingFormat
    procedure ,public   ::      ComputingWriting

  End Type
  
  public    ::    ArrOpeningFiles
  public    ::    ArrAllocating
  public    ::    ArrDeAllocating
  public    ::    ArrCreatingVectors
  public    ::    ArrMergeIntExch
  logical   ,parameter    ::    i_Debug_Global = .false.

  contains


! ==============================================================================================================
!   CREATING FORMAT for ARRHENIUS COEFFICIENTS OUTPUT FILE
! ==============================================================================================================
Subroutine CreatingFormat( This, Input, Collision, Rates, iMol, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use Rates_Class           ,only:  Rates_Type
  use Parameters_Module     ,only:  Zero

  class(Arrhenius_Type)                     ,intent(inout)  ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  type(Rates_Type)                          ,intent(in)     ::    Rates
  integer                                   ,intent(in)     ::    iMol
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer   ,dimension(:,:) ,allocatable                    ::    AtomsToPairs
  integer   ,dimension(:)   ,allocatable                    ::    OtherAtoms
  integer                                                   ::    iP, iPTemp, jP
  integer                                                   ::    iArr
  integer                                                   ::    iBins, jBins
  character(6)                                              ::    iBinsChar, jBinsChar 
  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "CreatingFormat" )
  !i_Debug_Loc   =     Logger%On()

  if (.NOT. allocated(This%MolOK)) then
    if (i_Debug_Loc) call Logger%Write( "Allocating a counter for the Molecules" )
    allocate(This%MolOK(Input%NMolecules), Stat=Status )
    if (Status/=0) call Error( "Error allocating This%MolOK" )
    if (i_Debug_Loc) call Logger%Write( "Allocated ", Input%NMolecules, " This%MolOK" )
    if (i_Debug_Loc) call Logger%Write( "Collision%Pairs(:)%To_Molecule ", Collision%Pairs(:)%To_Molecule )
    This%MolOK = 0
  end if
  
  select case (Input%NAtoms)       ! Setting the index of the two atoms associated to each pair
  case (3)
    allocate(AtomsToPairs(2,3))
    AtomsToPairs(:,1) = [1,2]
    AtomsToPairs(:,2) = [1,3]
    AtomsToPairs(:,3) = [2,3]
  end select
  
  select case (Input%NAtoms)       ! Setting the index of the two atoms associated to each pair
  case (3)
    allocate(OtherAtoms(3))
    OtherAtoms(1) = 3
    OtherAtoms(2) = 2
    OtherAtoms(3) = 1
  end select
  
  do iPTemp = 1,Collision%NPairs
    if (Collision%Pairs(iPTemp)%To_Molecule == iMol) then
      iP=iPTemp
      exit
    end if
  end do
  allocate(This%ReactionIn(Input%NBins(Collision%Pairs(iP)%To_Molecule)))
  iArr = 1
  do iBins = 1,Input%NBins(Collision%Pairs(iP)%To_Molecule)
    write(iBinsChar,'(I6)') iBins
    if (Input%NTtra == 1) then
      This%ReactionIn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(iP)%To_Molecule))) // "(" // trim(adjustl(iBinsChar)) // ")" // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(iP))))
    else
      This%ReactionIn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(iP)%To_Molecule))) // "_" // trim(adjustl(iBinsChar)) // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(iP))))
    end if
    iArr = iArr+1
  end do

  allocate(This%ReactionFn(Rates%FinStatesVec(Collision%NPairs+2)), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%ReactionFn" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%ReactionFn with dimension = (", Rates%FinStatesVec(Collision%NPairs+2),  ")" )
  
  allocate(This%TypeVec(Rates%FinStatesVec(Collision%NPairs+2)), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%TypeVec" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%TypeVec with dimension = (", Rates%FinStatesVec(Collision%NPairs+2),  ")" )
  This%TypeVec = 0
  
  allocate(This%ProcToPair(Rates%FinStatesVec(Collision%NPairs+2)), Stat=Status )
  if (Status/=0) call Error( "Error allocating This%ProcToPair" )
  if (i_Debug_Loc) call Logger%Write( "Allocated This%ProcToPair with dimension = (", Rates%FinStatesVec(Collision%NPairs+2),  ")" )
  This%ProcToPair = 0
  
  iArr = 1
  This%ReactionFn(iArr) = trim(adjustl(Input%AtomsName(1))) // "+" // trim(adjustl(Input%AtomsName(2))) // "+" // trim(adjustl(Input%AtomsName(3)))
  iArr = 2
  
  do jP = 1,Collision%NPairs 
    if (i_Debug_Loc) call Logger%Write( "Pair = ", jP )
  
    if (Collision%Pairs(jP)%To_Molecule == 0) then
      if (i_Debug_Loc) call Logger%Write( "Collision%Pairs(jP)%To_Molecule == 0" )

      if (trim(adjustl(Input%AtomsName(AtomsToPairs(1,jP)))) == trim(adjustl(Input%AtomsName(AtomsToPairs(2,jP))))) then
        This%ReactionFn(iArr) = trim(adjustl(Input%AtomsName(AtomsToPairs(1,jP)))) // "2+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
      else
        This%ReactionFn(iArr) = trim(adjustl(Input%AtomsName(AtomsToPairs(1,jP)))) // trim(adjustl(Input%AtomsName(AtomsToPairs(2,jP)))) // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
      end if
      This%TypeVec(iArr) = 2
    
    else
    
      do jBins = 1,Input%NBins(Collision%Pairs(jP)%To_Molecule)
        write(jBinsChar,'(I6)') jBins

        if (Input%NTtra == 1) then
          This%ReactionFn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(jP)%To_Molecule))) // "(" // trim(adjustl(jBinsChar)) // ")" // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
        else
          This%ReactionFn(iArr) = trim(adjustl(Input%Molecules_Name(Collision%Pairs(jP)%To_Molecule))) // "_" // trim(adjustl(jBinsChar)) // "+" // trim(adjustl(Input%AtomsName(OtherAtoms(jP))))
        end if

        if (This%MolOK(Collision%Pairs(jP)%To_Molecule) == 0) then
          if (Collision%Pairs(jP)%To_Molecule == Collision%Pairs(iP)%To_Molecule ) then
            This%TypeVec(iArr)    = 1
            This%ProcToPair(iArr) = 1 
          else 
            This%TypeVec(iArr)    = 2
            This%ProcToPair(iArr) = jP
          end if
        else
          This%TypeVec(iArr)    = -1
        end if
        !write(*,*) "Pair = ", jP, "; iBins = ", iBins, "; iArr = ", iArr, "; This%TypeVec(iArr) = ", This%TypeVec(iArr), "; This%ProcToPair(iArr) = ", This%ProcToPair(iArr)
        
        iArr = iArr+1
      end do

      if (i_Debug_Loc) call Logger%Write( "This%ProcToPair = ", This%ProcToPair )
      
      This%MolOK(Collision%Pairs(jP)%To_Molecule) = 1
      
    end if
  end do
  deallocate(AtomsToPairs)

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   OPENING FILE FOR WRITING ARRHENIUS COEFFICIENTS
! ==============================================================================================================
!______________________________________________________________________________________________________________!
Subroutine ArrOpeningFiles( Input, iMol, iBins, iTint, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  integer                                   ,intent(in)     ::    iMol
  integer                                   ,intent(in)     ::    iBins
  integer                                   ,intent(in)     ::    iTint
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  character(2)                                              ::    iTintChar
  character(:)                             ,allocatable     ::    iTintCharTemp
  character(:)                             ,allocatable     ::    FileName
  logical                                                   ::    exist_flag
  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ArrOpeningFiles" )
  !i_Debug_Loc   =     Logger%On()
  
  call system('mkdir -p ' // trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/')
  
  if (Input%NTint==1) then
    iTintCharTemp = 'Therm' 
    
    if (iBins == Input%BinStart(iMol)) then
    
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '.dat'
      if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
      inquire(file=trim(adjustl(FileName)), exist=exist_flag)
      if (exist_flag) then
        open( File=FileName, Unit=102, status='OLD', position="append", action="write", iostat=Status )
      else
        open( File=FileName, Unit=102, status='REPLACE', iostat=Status )
        write(102,*) 'Units=cm^3/s'
      end if 
      
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-Diss.dat'
      if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
      inquire(file=trim(adjustl(FileName)), exist=exist_flag)
      if (exist_flag) then
        open( File=FileName, Unit=103, status='OLD', position="append", action="write", iostat=Status )
      else
        open( File=FileName, Unit=103, status='REPLACE', iostat=Status )
        write(103,*) 'Units=cm^3/s'
      end if 
      
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-Exo.dat'
      if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
      inquire(file=trim(adjustl(FileName)), exist=exist_flag)
      if (exist_flag) then
        open( File=FileName, Unit=104, status='OLD', position="append", action="write", iostat=Status )
      else
        open( File=FileName, Unit=104, status='REPLACE', iostat=Status )
        write(104,*) 'Units=cm^3/s'
      end if 
      
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-Endo.dat'
      if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
      inquire(file=trim(adjustl(FileName)), exist=exist_flag)
      if (exist_flag) then
        open( File=FileName, Unit=105, status='OLD', position="append", action="write", iostat=Status )
      else
        open( File=FileName, Unit=105, status='REPLACE', iostat=Status )
        write(105,*) 'Units=cm^3/s'
      end if 
      
      FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-Reac.dat'
      if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
      inquire(file=trim(adjustl(FileName)), exist=exist_flag)
      if (exist_flag) then
        open( File=FileName, Unit=106, status='OLD', position="append", action="write", iostat=Status )
      else
        open( File=FileName, Unit=106, status='REPLACE', iostat=Status )
        write(106,*) 'Units=cm^3/s'
      end if 
    
    end if
    
  else
    write(iTintChar,'(I2)') iTint
    iTintCharTemp = "Tint" // trim(adjustl(iTintChar))

    FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '.dat'
    if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
    inquire(file=trim(adjustl(FileName)), exist=exist_flag)
    if (exist_flag) then
      open( File=FileName, Unit=102, status='OLD', position="append", action="write", iostat=Status )
    else
      open( File=FileName, Unit=102, status='REPLACE', iostat=Status )
      write(102,*) 'Units=cm^3/s'
    end if 
    
    FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-Diss.dat'
    if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
    inquire(file=trim(adjustl(FileName)), exist=exist_flag)
    if (exist_flag) then
      open( File=FileName, Unit=103, status='OLD', position="append", action="write", iostat=Status )
    else
      open( File=FileName, Unit=103, status='REPLACE', iostat=Status )
      write(103,*) 'Units=cm^3/s'
    end if 
    
    FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-InelExo.dat'
    if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
    inquire(file=trim(adjustl(FileName)), exist=exist_flag)
    if (exist_flag) then
      open( File=FileName, Unit=104, status='OLD', position="append", action="write", iostat=Status )
    else
      open( File=FileName, Unit=104, status='REPLACE', iostat=Status )
      write(104,*) 'Units=cm^3/s'
    end if 
    
    FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-InelEndo.dat'
    if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
    inquire(file=trim(adjustl(FileName)), exist=exist_flag)
    if (exist_flag) then
      open( File=FileName, Unit=105, status='OLD', position="append", action="write", iostat=Status )
    else
      open( File=FileName, Unit=105, status='REPLACE', iostat=Status )
      write(105,*) 'Units=cm^3/s'
    end if 
    
    FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-ExchExo.dat'
    if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
    inquire(file=trim(adjustl(FileName)), exist=exist_flag)
    if (exist_flag) then
      open( File=FileName, Unit=106, status='OLD', position="append", action="write", iostat=Status )
    else
      open( File=FileName, Unit=106, status='REPLACE', iostat=Status )
      write(106,*) 'Units=cm^3/s'
    end if 
    
    FileName = trim(adjustl(Input%OutputDir)) // '/'// trim(adjustl(Input%System)) // '/Arrhenius-Coeffs/' // trim(adjustl(iTintCharTemp)) // '-ExchEndo.dat'
    if (i_Debug_Loc) call Logger%Write( "Writing File: ", FileName )
    inquire(file=trim(adjustl(FileName)), exist=exist_flag)
    if (exist_flag) then
      open( File=FileName, Unit=107, status='OLD', position="append", action="write", iostat=Status )
    else
      open( File=FileName, Unit=107, status='REPLACE', iostat=Status )
      write(107,*) 'Units=cm^3/s'
    end if 
     
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   CREATING FORMAT for ARRHENIUS COEFFICIENTS OUTPUT FILE
! ==============================================================================================================
Subroutine ArrAllocating( Input, Collision, BinnedMolecule, Rates, iBins, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Rates_Class           ,only:  Rates_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  Type(BinsContainer_Type)                  ,intent(inout)  ::    BinnedMolecule
  type(Rates_Type)                          ,intent(in)     ::    Rates
  integer                                   ,intent(in)     ::    iBins
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ArrAllocating" )
  !i_Debug_Loc   =     Logger%On()
  
   if (.NOT. allocated(BinnedMolecule%Bin(iBins)%RateConst_Arr) ) then
            
    allocate(BinnedMolecule%Bin(iBins)%RateConst_Arr(Rates%FinStatesVec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
    if (Status/=0) call Error( "Error allocating BinnedMolecule%Bin(iBins)%RateConst_Arr" )
    if (i_Debug_Loc) call Logger%Write( "Allocated BinnedMolecule%Bin(iBins)%RateConst_Arr with dimension = (", Rates%FinStatesVec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
    BinnedMolecule%Bin(iBins)%RateConst_Arr = Zero
    
    allocate(BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(Rates%FinStatesVec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
    if (Status/=0) call Error( "Error allocating BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2" )
    if (i_Debug_Loc) call Logger%Write( "Allocated BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2 with dimension = (", Rates%FinStatesVec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
    BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2 = Zero
  
    allocate(BinnedMolecule%Bin(iBins)%CArr(Rates%FinStatesVec(Collision%NPairs+2),max(Input%NTtra,3)), Stat=Status )
    if (Status/=0) call Error( "Error allocating BinnedMolecule%Bin(iBins)%CArr" )
    if (i_Debug_Loc) call Logger%Write( "Allocated BinnedMolecule%Bin(iBins)%CArr with dimension = (", Rates%FinStatesVec(Collision%NPairs+2), ",", max(Input%NTtra,3), ")" )
    BinnedMolecule%Bin(iBins)%CArr = Zero
    
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   CREATING FORMAT for ARRHENIUS COEFFICIENTS OUTPUT FILE
! ==============================================================================================================
Subroutine ArrDeAllocating( Input, BinnedMolecule, iBins, i_Debug )

  use Input_Class           ,only:  Input_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Parameters_Module     ,only:  Zero

  type(Input_Type)                          ,intent(in)     ::    Input
  Type(BinsContainer_Type)                  ,intent(inout)  ::    BinnedMolecule
  integer                                   ,intent(in)     ::    iBins
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ArrDeAllocating" )
  !i_Debug_Loc   =     Logger%On()
  
  if ( allocated(BinnedMolecule%Bin(iBins)%RateConst_Arr) ) then
            
    deallocate(BinnedMolecule%Bin(iBins)%RateConst_Arr, Stat=Status )
    if (Status/=0) call Error( "Error deallocating BinnedMolecule%Bin(iBins)%RateConst_Arr" )
    
    deallocate(BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2, Stat=Status )
    if (Status/=0) call Error( "Error deallocating BinnedMolecule%Bin(iBins)%RateConst_Arr" )
  
    deallocate(BinnedMolecule%Bin(iBins)%CArr, Stat=Status )
    if (Status/=0) call Error( "Error deallocating BinnedMolecule%Bin(iBins)%CArr" )            
  
  end if
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   CREATING VECTORS OF TEMPERATURES and RATES
! ==============================================================================================================
Subroutine ArrCreatingVectors( Input, BinnedMolecule, iBins, Ttra_Vec, iTtra, i_Debug )

  use Input_Class           ,only:  Input_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Parameters_Module     ,only:  Zero, Two, Ten

  type(Input_Type)                          ,intent(in)     ::    Input
  Type(BinsContainer_Type)                  ,intent(inout)  ::    BinnedMolecule
  integer                                   ,intent(in)     ::    iBins
  real(rkp)                  ,dimension(:)  ,intent(inout)  ::    Ttra_Vec
  integer                                   ,intent(in)     ::    iTtra
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ArrCreatingVectors" )
  !i_Debug_Loc   =     Logger%On()
  

  if (Input%NTtra == 1) then
  
    Ttra_Vec(iTtra) = Input%Ttra
    BinnedMolecule%Bin(iBins)%RateConst_Arr(:,iTtra)       = BinnedMolecule%Bin(iBins)%RateConst_Final(:)
    BinnedMolecule%Bin(iBins)%RateConst_Arr(:,2:3)         = Zero
    
    BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(:,iTtra) = BinnedMolecule%Bin(iBins)%RateConst_Sigma2(:)
    BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(:,2:3)   = Zero
  
  elseif (Input%NTtra == 2) then
  
    Ttra_Vec(iTtra) = Input%Ttra
    
    BinnedMolecule%Bin(iBins)%RateConst_Arr(:,iTtra)       = BinnedMolecule%Bin(iBins)%RateConst_Final(:)
    BinnedMolecule%Bin(iBins)%RateConst_Arr(:,3)           = BinnedMolecule%Bin(iBins)%RateConst_Arr(:,2)  * Ten
    
    BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(:,iTtra) = BinnedMolecule%Bin(iBins)%RateConst_Sigma2(:)
    BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(:,3)     = BinnedMolecule%Bin(iBins)%RateConst_Sigma2(:) * Ten**2
    
    Ttra_Vec(3) = Ttra_Vec(2) * Two
  
  elseif (Input%NTtra > 2) then
  
    Ttra_Vec(iTtra) = Input%Ttra

    BinnedMolecule%Bin(iBins)%RateConst_Arr(:,iTtra)       = BinnedMolecule%Bin(iBins)%RateConst_Final(:)
    BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(:,iTtra) = BinnedMolecule%Bin(iBins)%RateConst_Sigma2(:)
    
  end if
  
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   CREATING VECTORS OF TEMPERATURES and RATES
! ==============================================================================================================
Subroutine ArrMergeIntExch( Input, Collision, BinnedMolecule, Rates, iBins, iTtra, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use Rates_Class           ,only:  Rates_Type
  use Parameters_Module     ,only:  Zero, Two, Ten

  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  Type(BinsContainer_Type)                  ,intent(inout)  ::    BinnedMolecule
  type(Rates_Type)                          ,intent(in)     ::    Rates
  integer                                   ,intent(in)     ::    iBins
  integer                                   ,intent(in)     ::    iTtra
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    Status
  integer                                                   ::    StartProc, EndProc
  integer                                                   ::    iP
  logical                                                   ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ArrMergeIntExch" )
  !i_Debug_Loc   =     Logger%On()
  
  StartProc = Rates%FinStatesVec(2) + 1
  EndProc   = Rates%FinStatesVec(3) 
  
  do iP = 2,Collision%NPairs
  
    if (Collision%Pairs(iP)%To_Molecule == Collision%Pairs(1)%To_Molecule) then
    
      if (i_Debug_Loc) call Logger%Write( "Considering Internal Exchange as Inelastic Process: Merging Rates for Pair iP = ", iP )
  
      BinnedMolecule%Bin(iBins)%RateConst_Arr(StartProc:EndProc,iTtra)       = BinnedMolecule%Bin(iBins)%RateConst_Arr(StartProc:EndProc,iTtra)    + &
                                                                               BinnedMolecule%Bin(iBins)%RateConst_Arr(Rates%FinStatesVec(iP+1)+1:Rates%FinStatesVec(iP+2),iTtra)

      BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(StartProc:EndProc,iTtra) = BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(StartProc:EndProc,iTtra) + &
                                                                               BinnedMolecule%Bin(iBins)%RateConst_ArrSigma2(Rates%FinStatesVec(iP+1)+1:Rates%FinStatesVec(iP+2),iTtra)
      
    end if
    
  end do 
  
  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
Subroutine ComputingWriting( This, Input, Collision, BinnedMolecule, iMol, LevelsContainer, Rates, iBins, iBinsChar, Ttra_Vec, i_Debug )

  use Input_Class           ,only:  Input_Type
  use Collision_Class       ,only:  Collision_Type
  use BinsContainer_Class   ,only:  BinsContainer_Type
  use LevelsContainer_Class ,only:  LevelsContainer_Type
  use Rates_Class           ,only:  Rates_Type
  use Parameters_Module     ,only:  Zero, One, Two, Ten
  use C_interface           ,only:  fit 


  class(Arrhenius_Type)                     ,intent(in)     ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Collision_Type)                      ,intent(in)     ::    Collision
  Type(BinsContainer_Type) ,dimension(:)    ,intent(inout)  ::    BinnedMolecule
  integer                                   ,intent(in)     ::    iMol
  Type(LevelsContainer_Type) ,dimension(:)  ,intent(in)     ::    LevelsContainer
  type(Rates_Type)                          ,intent(in)     ::    Rates
  integer                                   ,intent(in)     ::    iBins
  character(6)                              ,intent(in)     ::    iBinsChar
  real(rkp)      ,dimension(:)              ,intent(in)     ::    Ttra_Vec
  logical                         ,optional ,intent(in)     ::    i_Debug

  integer                                                   ::    iP, iPTemp
  integer                                                   ::    iTtra
  real(rkp)                                                 ::    CorrFactor
  integer                                                   ::    iBinsFn
  character(:)                ,allocatable                  ::    format_char
  integer                                                   ::    errInt
  real(rkp)                                                 ::    err
  integer                                                   ::    Status
  logical                                   ,parameter      ::    CheckFlg     = .true.
  real(rkp)                                 ,parameter      ::    Tolerance    = 1.d0
  real(rkp)      ,dimension(Input%NTtra)                    ::    Trial
  real(rkp)      ,dimension(Input%NTtra)                    ::    ErrorTemp
  integer                                                   ::    NTtraFound
  logical        ,dimension(max(Input%NTtra,3))             ::    TFound       
  real(rkp)      ,dimension(max(Input%NTtra,3))             ::    RatesTempVec, RatesTempVecBis
  real(rkp)      ,dimension(max(Input%NTtra,3))             ::    TempTempVec  
  logical                                                   ::    i_Debug_Loc


  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ComputingWriting" )
  !i_Debug_Loc   =     Logger%On()
  
  do iPTemp = 1,Collision%NPairs
    if (Collision%Pairs(iPTemp)%To_Molecule == iMol) then
      iP=iPTemp
      exit
    end if
  end do

  do iBinsFn = 1,Rates%FinStatesVec(Collision%NPairs+2)
    
    CorrFactor = One
    if ( (iBinsFn == 1) ) then
      CorrFactor = Input%DissCorrFactor
    end if 
    BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,:) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,:) * CorrFactor
    
    if (sum(BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,:)) /= Zero) then
        
              
      ! ==============================================================================================================
      !   COMPUTING ARRHENIUS COEFFICIENTS
      ! ==============================================================================================================
      if (Input%NTtra == 1) then
        
        if (BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,1) < Input%MinForRates) then 
          BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,1) = Input%MinForRates
          errInt = 1
        else
          errInt = 0
        end if
!        if (BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(iBinsFn) < Input%MinForRates) then 
!          BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(iBinsFn) = Input%MinForRates
!        end if
        
!        BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Final(iBinsFn)
        BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,1)
        BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,2) = Zero
        BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,3) = Zero
        
        
      else
        

        ! NTtraFound   = Input%NTtra
        ! TempTempVec  = Ttra_Vec
        ! RatesTempVec = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,:)

        NTtraFound   = 0
        TFound       = .False.
        RatesTempVec = Zero
        TempTempVec  = Zero
        do iTtra = 1,Input%NTtra
          if (BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,iTtra) < Input%MinForRates) then 
            BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,iTtra) = Input%MinForRates
          else
            TFound(iTtra)            = .True.
            NTtraFound               = NTtraFound + 1
            TempTempVec(NTtraFound)  = Ttra_Vec(iTtra)
            RatesTempVec(NTtraFound) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,iTtra)
          end if
        end do

        if ( (.not. TFound(1)) .and. (.not. TFound(2)) ) then
          iTtra                    = 1
          TFound(iTtra)            = .True.
          NTtraFound               = NTtraFound + 1
          TempTempVec(NTtraFound)  = Ttra_Vec(iTtra)
          RatesTempVec(NTtraFound) = BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,iTtra)
        end if
        RatesTempVecBis = RatesTempVec
        
        !write(*,*) ' '
        if (NTtraFound > 2) then
          !write(*,*) 'Orig Rates   ', BinnedMolecule(iMol)%Bin(iBins)%RateConst_Arr(iBinsFn,:)
          !write(*,*) 'Temperatures ', TempTempVec(1:NTtraFound)
          !write(*,*) 'Acc Rates    ', RatesTempVec(1:NTtraFound)
          !write(*,*) 'Found Vec    ', TFound
          call fit(NTtraFound, TempTempVec(1:NTtraFound), RatesTempVec(1:NTtraFound), BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1), BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,2), BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,3), err )
          !write(*,*) 'Arr Coeff    ', BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
          errInt = int(err) - 1
          !write(*,*) err
        else
          errInt = 1
          !write(*,*) 'Less than 3 Temperatures'
          !write(*,*) 'Found Vec    ', TFound
        end if

        if (errInt == 0) then
        
          if ((BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1) >= Input%MaxForRates) .or. (BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1) == Input%MinForRates) .or. (BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1) == Zero)) then
            
            errInt = 1
        
          elseif (CheckFlg) then
            
            Trial(1:NTtraFound) = exp( log(BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1)) + BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,2) * log(TempTempVec(1:NTtraFound)) - BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,3) / TempTempVec(1:NTtraFound) )
            !write(*,*) 'Trial        ', Trial 

            !write(*,*) 'RatesTempVecBis ', RatesTempVecBis
            ErrorTemp(1:NTtraFound) = abs( Trial(1:NTtraFound) - RatesTempVecBis(1:NTtraFound) ) / RatesTempVecBis(1:NTtraFound)
            !write(*,*) 'ErrorTemp    ', ErrorTemp 
            
            if ((maxval(ErrorTemp) /= maxval(ErrorTemp)) .or. (maxval(ErrorTemp(1:NTtraFound)) > Tolerance)) then
              errInt = 1
              !write(*,*) 'errInt set to ', errInt 
            end if
          
          end if  
            
        end if
        
      end if         
      ! ==============================================================================================================
        
        
      ! ==============================================================================================================
      !   WRITING ARRHENIUS COEFFICIENTS
      ! ==============================================================================================================
      if (iBinsFn == 1) then
        format_char =  "'" // trim(adjustl(This%ReactionIn(iBins))) // "=" // trim(adjustl(This%ReactionFn(iBinsFn))) // ":', SP,  ES11.4, ',', SP,  ES11.4, ',', SP,  ES11.4, ',2'"
      else
        format_char =  "'" // trim(adjustl(This%ReactionIn(iBins))) // "=" // trim(adjustl(This%ReactionFn(iBinsFn))) // ":', SP,  ES11.4, ',', SP,  ES11.4, ',', SP,  ES11.4, ',5'"
      end if
      format_char = '(' // trim(adjustl(format_char)) // ')'
      
      if (errInt == 0) then
      
        !if ((Input%WriteArrFlg) .and. (BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1) >= Input%MinForRates)) then
        if ( (Input%WriteArrFlg) .and. (This%TypeVec(iBinsFn) /= -1) ) then

          if ( (This%TypeVec(iBinsFn) == 1) .and. (iBins > iBinsFn - Rates%FinStatesVec(iP+1) ) .and. (Input%Kinetics_InelFlg) .and. (Input%Kinetics_ExoFlg) )  then
            write(104,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
            write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)

          elseif ( (This%TypeVec(iBinsFn) == 1) .and. (iBins < iBinsFn - Rates%FinStatesVec(iP+1) ) .and. (Input%Kinetics_InelFlg) .and. (Input%Kinetics_EndoFlg) ) then
            write(105,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
          
          elseif ( (This%TypeVec(iBinsFn) == 0)  .and. (Input%Kinetics_DissFlg) ) then
            write(103,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
            write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)

          elseif (This%TypeVec(iBinsFn) == 2) then

            if ( (BinnedMolecule(iMol)%Bin(iBins)%ToteinteV + LevelsContainer(iMol)%MineinteV >= BinnedMolecule(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%Bin(iBinsFn-Rates%FinStatesVec(This%ProcToPair(iBinsFn)+1))%ToteinteV + LevelsContainer(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%MineinteV) .and. (Input%Kinetics_ExchFlg) .and. (Input%Kinetics_ExoFlg) ) then
              write(106,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
              write(102,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
            
            elseif ( (BinnedMolecule(iMol)%Bin(iBins)%ToteinteV + LevelsContainer(iMol)%MineinteV < BinnedMolecule(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%Bin(iBinsFn-Rates%FinStatesVec(This%ProcToPair(iBinsFn)+1))%ToteinteV + LevelsContainer(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%MineinteV) .and. (Input%Kinetics_ExchFlg) .and. (Input%Kinetics_EndoFlg) ) then
              write(107,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(iBins)%CArr(iBinsFn,1:3)
            
            end if

          end if

        end if
        


        if ( (Input%WriteKonigFlg) .and. (This%TypeVec(iBinsFn) /= -1) ) then

          if ( (This%TypeVec(iBinsFn) == 1) .and. (IBins > iBinsFn - Rates%FinStatesVec(iP+1)) .and. (Input%Kinetics_InelFlg) .and. (Input%Kinetics_ExoFlg) ) then
            write(122,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)

          elseif ( (This%TypeVec(iBinsFn) == 1) .and. (IBins < iBinsFn - Rates%FinStatesVec(iP+1) ) .and. (Input%Kinetics_InelFlg) .and. (Input%Kinetics_EndoFlg) ) then
            write(122,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)

          elseif ( (This%TypeVec(iBinsFn) == 0) .and. (Input%Kinetics_DissFlg) ) then
            write(122,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)

          elseif (This%TypeVec(iBinsFn) == 2) then

            if ( (BinnedMolecule(iMol)%Bin(iBins)%ToteinteV + LevelsContainer(iMol)%MineinteV >=  BinnedMolecule(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%Bin(iBinsFn-Rates%FinStatesVec(This%ProcToPair(iBinsFn)+1))%ToteinteV +  LevelsContainer(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%MineinteV) .and.  (Input%Kinetics_ExchFlg) .and. (Input%Kinetics_ExoFlg) ) then
              write(122,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
            
            elseif ( (BinnedMolecule(iMol)%Bin(iBins)%ToteinteV + LevelsContainer(iMol)%MineinteV <  BinnedMolecule(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%Bin(iBinsFn-Rates%FinStatesVec(This%ProcToPair(iBinsFn)+1))%ToteinteV + LevelsContainer(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%MineinteV) .and. (Input%Kinetics_ExchFlg) .and. (Input%Kinetics_EndoFlg) ) then
              write(122,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
            end if

          end if

          if (iBins == Input%BinFinal(iMol)) then

            if ( (Input%Kinetics_COCDissFlg) ) then
              write(122,'(A)') 'CO(' // trim(adjustl(iBinsChar)) // ')+C=C+O+C:+5.6458E-04,-1.0E+00,+129E+03,2'
            end if

          end if

        end if



        if ( (Input%WriteHegelFlg) .and. (This%TypeVec(iBinsFn) /= -1) ) then
        
          if ( (This%TypeVec(iBinsFn) == 1) .and. (IBins > iBinsFn - Rates%FinStatesVec(iP+1)) .and. (Input%Kinetics_InelFlg) .and. (Input%Kinetics_ExoFlg) ) then
            write(132,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
          
          elseif ( (This%TypeVec(iBinsFn) == 1) .and. (IBins < iBinsFn - Rates%FinStatesVec(iP+1) ) .and. (Input%Kinetics_InelFlg) .and. (Input%Kinetics_EndoFlg) ) then
            write(132,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
            
          elseif ( (This%TypeVec(iBinsFn) == 0) .and. (Input%Kinetics_DissFlg) ) then
            write(132,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
           
          elseif (This%TypeVec(iBinsFn) == 2) then
            
            if ( (BinnedMolecule(iMol)%Bin(iBins)%ToteinteV + LevelsContainer(iMol)%MineinteV >= BinnedMolecule(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%Bin(iBinsFn-Rates%FinStatesVec(This%ProcToPair(iBinsFn)+1))%ToteinteV + LevelsContainer(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%MineinteV) .and. (Input%Kinetics_ExchFlg) .and. (Input%Kinetics_ExoFlg) ) then
              write(132,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
            
            elseif ( (This%TypeVec(iBinsFn) == 2) .and. (BinnedMolecule(iMol)%Bin(iBins)%ToteinteV + LevelsContainer(iMol)%MineinteV < BinnedMolecule(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%Bin(iBinsFn-Rates%FinStatesVec(This%ProcToPair(iBinsFn)+1))%ToteinteV + LevelsContainer(Collision%Pairs(This%ProcToPair(iBinsFn))%To_Molecule)%MineinteV) .and. (Input%Kinetics_ExchFlg) .and. (Input%Kinetics_EndoFlg) ) then
              write(132,trim(adjustl(format_char))) BinnedMolecule(iMol)%Bin(IBins)%CArr(iBinsFn,1:3)
            end if

          end if 
          
          if (iBins == Input%BinFinal(iMol)) then
          
            if ( (Input%Kinetics_COCDissFlg) ) then
              write(132,'(A)') 'CO(' // trim(adjustl(iBinsChar)) // ')+C=C+O+C:+5.6458E-04,-1.0E+00,+129E+03,2'
            end if
          
          end if          
                    
        end if

      
      end if
      ! ==============================================================================================================
      
    end if
  
  end do   ! iBinsFn
  
  if ( (Input%WriteArrFlg) .and. ((iBins == Input%BinFinal(iMol)) .or. (Input%NTint /= 1)) ) then      
    close(102)  
    close(103)  
    close(104)  
    close(105)  
    close(106) 
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!

End Module
