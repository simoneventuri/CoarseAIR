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

! ==============================================================================================================
!            Regression Testing for CoarseAIR (Coarse-Grained QCT for Atmospheric Mixtures)
! ==============================================================================================================
! Driver program for classical trajectory program for a general number of atoms.
! The code uses a bulirsch-stoer ODE integrator.
!
! Mandatory Input Arguments: - Total Number of Processors
!
! ==============================================================================================================
Program RegressionTesting

  use Parameters_Module     ,only:  rkp, Zero
  use Input_Class           ,only:  Input_Type
  use Error_Class           ,only:  Error
  use RegressionTesting_Module     
  
  implicit none

  type(Input_Type)                        ::    Input
  
  character(:)  ,allocatable              ::    BenchFolder
  character(:)  ,allocatable              ::    TestFolder
  character(150)                          ::    FileName_Bench
  character(150)                          ::    FileName_Test
  
  integer ,dimension(:,:)   ,allocatable  ::    IntMatrix_Bench
  integer ,dimension(:,:)   ,allocatable  ::    IntMatrix_Test
  
  real(rkp) ,dimension(:,:) ,allocatable  ::    RealMatrix_Bench
  real(rkp) ,dimension(:,:) ,allocatable  ::    RealMatrix_Test
  
  real(rkp) ,dimension(:)   ,allocatable  ::    IntError
  real(rkp) ,dimension(:)   ,allocatable  ::    RealError
  
  integer                                 ::    IntFlg
  
  integer                                 ::    iMol
  integer                                 ::    iTtra
  character(10)                           ::    T_char
  integer                                 ::    iLevels
  character(5)                            ::    iLevelsChar
  integer                                 ::    iTemp                                     
  integer                                 ::    Unit
 
  
!  call CPU_Time( StartTime )
  

!! ==============================================================================================================
!!  READING THE PROGRAM ARGUMENT
!! ==============================================================================================================            
!  call getarg( 1, Input%NProc_char )
!  if (trim(adjustl(Input%NProc_char)) .eq. '') then
    Input%NProc_char  = '1'
    Input%NProc       =  1
!  else 
!    read( Input%NProc_char, '(I3)' ) Input%NProc
!  end if
!! ==============================================================================================================


! ==============================================================================================================
!   INITIALIZING INPUT
! ==============================================================================================================
  call Input%Initialize( i_Debug=.false. )
! ==============================================================================================================


  if (trim(adjustl(Input%BSortMethod(1))) == 'State-Specific') then
    BenchFolder = './BenchMarks_StS'
    TestFolder  = './Test_StS'
    Input%BinStart(1) = Input%BSortInfo(1,1)
    Input%BinFinal(1) = Input%BSortInfo(1,2)
  elseif (trim(adjustl(Input%BSortMethod(1))) == 'Vib-Specific') then
    BenchFolder = './BenchMarks_VS'
    TestFolder  = './Test_VS'
    Input%BinStart(1) = Input%BSortInfo(1,1)
    Input%BinFinal(1) = Input%BSortInfo(1,2)
  elseif (trim(adjustl(Input%BSortMethod(1))) == 'RoVib-CG') then
    BenchFolder = './BenchMarks_CG'
    TestFolder  = './Test_CG'
    Input%BinStart(1) = 1
    Input%BinFinal(1) = Input%BSortInfo(1,1)
  end if
  
  
! ==============================================================================================================
!   READING qnsEnBin.dat
! ==============================================================================================================
  do iMol = 1,Input%NMolecules
  
    FileName_Bench = trim(adjustl(BenchFolder)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qnsEnBin.dat'
    FileName_Test  = trim(adjustl(TestFolder))  // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qnsEnBin.dat'
    
    iTemp = 0   
    call ReadingRegTest( Input, FileName_Bench, iTemp, 1, IntMatrix_Bench, RealMatrix_Bench)
    call ReadingRegTest( Input, FileName_Test,  iTemp, 1, IntMatrix_Test,  RealMatrix_Test)
    
    call ComputingErrorsRegTest( 1, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
    
    if (IntError(2)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Levels' vqns! --> Check qnsEnBin.dat!"
      write(*,'(A)')    "      Error in Levels' vqns! --> Check qnsEnBin.dat!"
      close(Unit)
    end if
    if (IntError(3)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Levels' jqns! --> Check qnsEnBin.dat!"
      write(*,'(A)')    "      Error in Levels' jqns! --> Check qnsEnBin.dat!"
      close(Unit)
    end if
    if (RealError(1) > 1.d-100) then 
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Levels' Energies! --> Check qnsEnBin.dat!"
      write(*,'(A)')    "      Error in Levels' Energies! --> Check qnsEnBin.dat!"
      close(Unit)
    end if
    if (RealError(2) > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Levels' Degeneracies! --> Check qnsEnBin.dat!"
      write(*,'(A)')    "      Error in Levels' Degeneracies! --> Check qnsEnBin.dat!"
      close(Unit)
    end if
    if (IntError(4)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in LevelsToBins! --> Check qnsEnBin.dat!"
      write(*,'(A)')    "      Error in LevelsToBins! --> Check qnsEnBin.dat!"
      close(Unit)
    end if
    
    call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
 
    write(*,'(A)')    "      Checked qnsEnBin.dat!"     
  end do
! ==============================================================================================================  


!! ==============================================================================================================
!!   READING qns_Bin.dat
!! ==============================================================================================================
!  do iMol = 1,Input%NMolecules
!  
!    FileName_Bench = trim(adjustl(BenchFolder)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
!                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qns_Bin.dat'
!    FileName_Test  = trim(adjustl(TestFolder))  // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
!                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qns_Bin.dat'
!    
!    iTemp = 0           
!    call ReadingRegTest( FileName_Bench, iTemp, 2, IntMatrix_Bench, RealMatrix_Bench)
!    call ReadingRegTest( FileName_Test,  iTemp, 2, IntMatrix_Test,  RealMatrix_Test)
!    
!    call ComputingErrorsRegTest( 2, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
!    
!    if (IntError(1)  > 1.d-100) then
!      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
!      write(Unit,'(A)') "      Error in Levels' vqns! --> Check qns_Bin.dat!"
!      write(*,'(A)')    "      Error in Levels' vqns! --> Check qns_Bin.dat!"
!      close(Unit)
!    end if
!    if (IntError(2)  > 1.d-100) then
!      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
!      write(Unit,'(A)') "      Error in Levels' jqns! --> Check qns_Bin.dat!"
!      write(*,'(A)')    "      Error in Levels' jqns! --> Check qns_Bin.dat!"
!      close(Unit)
!    end if
!    if (RealError(1) > 1.d-100) then 
!      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
!      write(Unit,'(A)') "      Error in Levels' Energies! --> Check qns_Bin.dat!"
!      write(*,'(A)')    "      Error in Levels' Energies! --> Check qns_Bin.dat!"
!      close(Unit)
!    end if
!    
!    call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
! 
!     write(*,'(A)')    "      Checked qns_Bin.dat!"   
!  end do
!! ==============================================================================================================  


! ==============================================================================================================
!   READING qnsFirst.dat
! ==============================================================================================================
  do iMol = 1,Input%NMolecules
  
    FileName_Bench = trim(adjustl(BenchFolder)) // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qnsFirst.dat'
    FileName_Test  = trim(adjustl(TestFolder))  // '/'// trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/qnsFirst.dat'
    
    iTemp = 0           
    call ReadingRegTest( Input, FileName_Bench, iTemp, 3, IntMatrix_Bench, RealMatrix_Bench)
    call ReadingRegTest( Input, FileName_Test,  iTemp, 3, IntMatrix_Test,  RealMatrix_Test)
    
    call ComputingErrorsRegTest( 3, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
    
    if (IntError(1)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Levels/Bins' First vqns! --> Check qnsFirst.dat!"
      write(*,'(A)')    "      Error in Levels/Bins' First vqns! --> Check qnsFirst.dat!"
      close(Unit)
    end if
    if (IntError(2)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Levels/Bins' First jqns! --> Check qnsFirst.dat!"
      write(*,'(A)')    "      Error in Levels/Bins' First jqns! --> Check qnsFirst.dat!"
      close(Unit)
    end if
    if (IntError(3)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Nb of Levels/Bins' per Bin! --> Check qnsFirst.dat!"
      write(*,'(A)')    "      Error in Nb of Levels/Bins' per Bin! --> Check qnsFirst.dat!"
      close(Unit)
    end if
    if (RealError(1) > 1.d-100) then 
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Percentage of Levels per Bins! --> Check qnsFirst.dat!"
      write(*,'(A)')    "      Error in Percentage of Levels per Bins! --> Check qnsFirst.dat!"
      close(Unit)
    end if
    
    call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
    
    write(*,'(A)')    "      Checked qnsFirst.dat!"
  end do
! ==============================================================================================================


! ==============================================================================================================
!   READING T*.dat
! ==============================================================================================================
  write(T_char,"(I10)") floor(Input%TInit)
  
  do iMol = 1,Input%NMolecules
    
    FileName_Bench = trim(adjustl(BenchFolder)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(T_char))// '.dat'
    FileName_Test  = trim(adjustl(TestFolder)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                     trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(T_char))// '.dat'
    
    iTemp = 0           
    call ReadingRegTest( Input, FileName_Bench, iTemp, 4, IntMatrix_Bench, RealMatrix_Bench)
    call ReadingRegTest( Input, FileName_Test,  iTemp, 4, IntMatrix_Test,  RealMatrix_Test)
    
    call ComputingErrorsRegTest( 4, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
    
    if (RealError(1)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Level/Bin's Ratios! --> Check T" // trim(adjustl(T_char))// ".dat!"
      write(*,'(A)')    "      Error in Level/Bin's Ratios! --> Check T" // trim(adjustl(T_char))// ".dat!"
      close(Unit)
    end if
    if (RealError(2)  > 1.d-100) then
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Level/Bin's Partition Functions! --> Check T" // trim(adjustl(T_char))// ".dat!"
      write(*,'(A)')    "      Error in Level/Bin's Partition Functions! --> Check T" // trim(adjustl(T_char))// ".dat!"
      close(Unit)
    end if
    if (RealError(3) > 1.d-100) then 
      open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
      write(Unit,'(A)') "      Error in Level/Bin's Energies! --> Check T" // trim(adjustl(T_char))// ".dat!"
      write(*,'(A)')    "      Error in Level/Bin's Energies! --> Check T" // trim(adjustl(T_char))// ".dat!"
      close(Unit)
    end if
    
    call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
    
    write(*,'(A)')    "      Checked T'TInit'.dat!"
  end do
    
  do iTtra = 1,Input%NTtra
  
    do iMol = 1,Input%NMolecules
    
      FileName_Bench = trim(adjustl(BenchFolder)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                       trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(Input%TtraVecIntChar(iTtra)))// '.dat'
      FileName_Test  = trim(adjustl(TestFolder)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(iMol))) // '/' // &
                       trim(adjustl(Input%Molecules_Name(iMol))) // '_' // trim(adjustl(Input%NBins_char(iMol))) // '/T' // trim(adjustl(Input%TtraVecIntChar(iTtra)))// '.dat'
      
      iTemp = 0
      call ReadingRegTest( Input, FileName_Bench, iTemp, 4, IntMatrix_Bench, RealMatrix_Bench)
      call ReadingRegTest( Input, FileName_Test,  iTemp, 4, IntMatrix_Test,  RealMatrix_Test)
      
      call ComputingErrorsRegTest( 4, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      if (RealError(1)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Level/Bin's Ratios! --> Check T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!"
        write(*,'(A)')    "      Error in Level/Bin's Ratios! --> Check T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!"
        close(Unit)
      end if
      if (RealError(2)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Level/Bin's Partition Functions! --> Check T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!"
        write(*,'(A)')    "      Error in Level/Bin's Partition Functions! --> Check T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!"
        close(Unit)
      end if
      if (RealError(3) > 1.d-100) then 
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Level/Bin's Energies! --> Check T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!"
        write(*,'(A)')    "      Error in Level/Bin's Energies! --> Check T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!"
        close(Unit)
      end if
      
      call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      write(*,'(A)')    "      Checked T" // trim(adjustl(Input%TtraVecIntChar(iTtra))) // ".dat!" 
    end do
      
  end do
! ==============================================================================================================  


! ==============================================================================================================
!   READING trajectories-old.out
! ==============================================================================================================
  do iTtra = 1,Input%NTtra
  
    do iLevels = Input%BinStart(1),Input%BinFinal(1)
      write(iLevelsChar,'(I5)') iLevels
    
      FileName_Bench = trim(adjustl(BenchFolder)) // '/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // &
                       '/Bins_' // trim(adjustl(iLevelsChar)) // '_0/trajectories-old.out'
      FileName_Test  = trim(adjustl(TestFolder))  // '/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // &
                       '/Bins_' // trim(adjustl(iLevelsChar)) // '_0/trajectories-old.out'
      
      iTemp = 0           
      call ReadingRegTest( Input, FileName_Bench, iTemp, 5, IntMatrix_Bench, RealMatrix_Bench)
      call ReadingRegTest( Input, FileName_Test,  iTemp, 5, IntMatrix_Test,  RealMatrix_Test)
      
      call ComputingErrorsRegTest( 5, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      if (IntError(1)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in iTraj! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in iTraj! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(1)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Impact Parameters Strata! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Impact Parameters Strata! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(2)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Sampled Impact Parameters! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Sampled Impact Parameters! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(3)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Initial Rotational QNs! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Initial Rotational QNs! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(4)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Initial Vibrational QNs! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Initial Vibrational QNs! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(5)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Initial Arrangements! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Initial Arrangements! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(6)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Final Rotational QNs! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Final Rotational QNs! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(7)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Final Vibrational QNs! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Final Vibrational QNs! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      if (RealError(8)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Final Arrangements! --> Check trajectories-old.dat!"
        write(*,'(A)')    "      Error in Final Arrangements! --> Check trajectories-old.dat!"
        close(Unit)
      end if
      
      call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      write(*,'(A)')    "      Checked trajectories-old.dat for Level/Bin " // trim(adjustl(iLevelsChar)) // " @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "K"
    end do
    
  end do
! ==============================================================================================================


! ==============================================================================================================
!   READING statistics.out
! ==============================================================================================================
  do iTtra = 1,Input%NTtra
  
    do iLevels = Input%BinStart(1),Input%BinFinal(1)
      write(iLevelsChar,'(I5)') iLevels
    
      FileName_Bench = trim(adjustl(BenchFolder)) // '/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // &
                       '/Bins_' // trim(adjustl(iLevelsChar)) // '_0/statistics.out'
      FileName_Test  = trim(adjustl(TestFolder))  // '/T_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) // &
                       '/Bins_' // trim(adjustl(iLevelsChar)) // '_0/statistics.out'
             
      iTemp = 0         
      call ReadingRegTest( Input, FileName_Bench, iTemp, 6, IntMatrix_Bench, RealMatrix_Bench)
      call ReadingRegTest( Input, FileName_Test,  iTemp, 6, IntMatrix_Test,  RealMatrix_Test)
      
      call ComputingErrorsRegTest( 6, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      if (IntError(4)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Initial Rotational QNs! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Initial Rotational QNs! --> Check statistics.dat!"
        close(Unit)
      end if
      if (IntError(5)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Initial Vibrational QNs! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Initial Vibrational QNs! --> Check statistics.dat!"
        close(Unit)
      end if
      if (IntError(6)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Initial Arrangements! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Initial Arrangements! --> Check statistics.dat!"
        close(Unit)
      end if
      if (IntError(1)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Final Rotational QNs! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Final Rotational QNs! --> Check statistics.dat!"
        close(Unit)
      end if
      if (IntError(2)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Final Vibrational QNs! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Final Vibrational QNs! --> Check statistics.dat!"
        close(Unit)
      end if
      if (IntError(3)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Final Arrangements! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Final Arrangements! --> Check statistics.dat!"
        close(Unit)
      end if
      if (RealError(1)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Cross Sections! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Cross Sections! --> Check statistics.dat!"
        close(Unit)
      end if
      if (RealError(2)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Cross Sections' Standard Deviations! --> Check statistics.dat!"
        write(*,'(A)')    "      Error in Cross Sections' Standard Deviations! --> Check statistics.dat!"
        close(Unit)
      end if
      
      call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      write(*,'(A)')    "      Checked statistics.dat for Level/Bin " // trim(adjustl(iLevelsChar)) // " @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "K"
    end do
    
  end do
! ==============================================================================================================


! ==============================================================================================================
!   READING Bin*.dat (Rates)
! ==============================================================================================================
  if ((((Input%CompQuant == 'yes') .or. (Input%CompQuant == 'YES')) .and. (Input%WriteFormRatesFlg)) .or. ((Input%RunPost == 'yes') .or. (Input%RunPost == 'YES'))) then
  
    do iTtra = 1,Input%NTtra
    
      do iLevels = Input%BinStart(1),Input%BinFinal(1)
        write(iLevelsChar,'(I5)') iLevels
      
        FileName_Bench = trim(adjustl(BenchFolder)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(1))) // '/Rates/T_' // &
                         trim(adjustl(Input%TtraVecIntChar(iTtra))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) //'/Bin' // trim(adjustl(iLevelsChar)) // '.dat'
        
        FileName_Test  = trim(adjustl(TestFolder))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(1))) // '/Rates/T_' // &
                         trim(adjustl(Input%TtraVecIntChar(iTtra))) // '_' // trim(adjustl(Input%TtraVecIntChar(iTtra))) //'/Bin' // trim(adjustl(iLevelsChar)) // '.dat'
               
        iTemp = 0
        call ReadingRegTest( Input, FileName_Bench, iTemp, 7, IntMatrix_Bench, RealMatrix_Bench)
        call ReadingRegTest( Input, FileName_Test,  iTemp, 7, IntMatrix_Test,  RealMatrix_Test)
        
        call ComputingErrorsRegTest( 7, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
        
        if (IntError(1)  > 1.d-100) then
          open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
          write(Unit,'(A)') "      Error in Processes! --> Check Bin" // trim(adjustl(iLevelsChar)) // ".dat!"
          write(*,'(A)')    "      Error in Processes! --> Check Bin" // trim(adjustl(iLevelsChar)) // ".dat!"
          close(Unit)
        end if
        if (RealError(1)  > 1.d-100) then
          open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
          write(Unit,'(A)') "      Error in Rates! --> Check Bin" // trim(adjustl(iLevelsChar)) // ".dat!"
          write(*,'(A)')    "      Error in Rates! --> Check Bin" // trim(adjustl(iLevelsChar)) // ".dat!"
          close(Unit)
        end if
        if (RealError(2)  > 1.d-100) then
          open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
          write(Unit,'(A)') "      Error in Rates' Standard Deviations! --> Check Bin" // trim(adjustl(iLevelsChar)) // ".dat!"
          write(*,'(A)')    "      Error in Rates' Standard Deviations! --> Check Bin" // trim(adjustl(iLevelsChar)) // ".dat!"
          close(Unit)
        end if
        
        call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
        
        write(*,'(A)')    "      Checked Bin" // trim(adjustl(iLevelsChar)) // ".dat @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "K"
      end do
      
    end do
    
  end if
! ==============================================================================================================


! ==============================================================================================================
!   READING Rates.dat
! ==============================================================================================================
  if ((Input%CompQuant == 'yes') .or. (Input%CompQuant == 'YES')) then
  
    if (Input%WriteAllTsFlg) then
  
      FileName_Bench = trim(adjustl(BenchFolder)) // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(1))) // '/Rates/Rates.dat'
      FileName_Test  = trim(adjustl(TestFolder))  // '/' // trim(adjustl(Input%System)) // '/' // trim(adjustl(Input%Molecules_Name(1))) // '/Rates/Rates.dat'
             
      iTemp = 0
      call ReadingRegTest( Input, FileName_Bench, iTemp, 8, IntMatrix_Bench, RealMatrix_Bench)
      call ReadingRegTest( Input, FileName_Test,  iTemp, 8, IntMatrix_Test,  RealMatrix_Test)
      
      call ComputingErrorsRegTest( 8, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
      
      if (IntError(1)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Levels/Bins! --> Check Rates.dat!"
        write(*,'(A)')    "      Error in Levels/Bins! --> Check Rates.dat!"
        close(Unit)
      end if
      if (IntError(2)  > 1.d-100) then
        open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
        write(Unit,'(A)') "      Error in Processes! --> Check Rates.dat!"
        write(*,'(A)')    "      Error in Processes! --> Check Rates.dat!"
        close(Unit)
      end if
      do iTtra = 1,Input%NTtra
        if (RealError(1 + int(2*(iTtra-1)))  > 1.d-100) then
          open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
          write(Unit,'(A)') "      Error in Rates @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "! --> Check Rates.dat!"
          write(*,'(A)')    "      Error in Rates @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "! --> Check Rates.dat!"
          close(Unit)
        end if
        if (RealError(2 + int(2*(iTtra-1)))  > 1.d-100) then
          open( File='./ERROR', NewUnit=Unit, Status='unknown', Access='append')
          write(Unit,'(A)') "      Error in Rates' Standard Deviations @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "! --> Check Rates.dat!"
          write(*,'(A)')    "      Error in Rates' Standard Deviations @ " // trim(adjustl(Input%TtraVecIntChar(iTtra))) // "! --> Check Rates.dat!"
          close(Unit)
        end if
      end do
      
      call DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )
    
    end if
    
  end if
! ==============================================================================================================


End Program
