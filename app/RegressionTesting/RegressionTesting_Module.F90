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

Module RegressionTesting_Module

  use Parameters_Module     ,only:  rkp, Zero
  use Error_Class           ,only:  Error
  
  implicit none
  
  contains


! ==============================================================================================================
!   READING OUTPUT FILE
! ==============================================================================================================
Subroutine ReadingRegTest( Input, FileName, NLines, IntFlg, IntMatrix, RealMatrix)
  
  use Input_Class              ,only:  Input_Type
  
  type(Input_Type)                          ,intent(in)     ::    Input
  character(150)                            ,intent(in)     ::    FileName
  integer                                   ,intent(inout)  ::    NLines
  integer                                   ,intent(in)     ::    IntFlg
  integer   ,dimension(:,:) ,allocatable    ,intent(out)    ::    IntMatrix
  real(rkp) ,dimension(:,:) ,allocatable    ,intent(out)    ::    RealMatrix
  
  integer                                                   ::    iLine
  integer                                                   ::    iTemp
  integer                                                   ::    Status, Unit
  
  open( File=trim(adjustl(FileName)), NewUnit=Unit, Action='READ', iostat=Status )
  if (Status/=0)then
    open( File='./ERROR', Unit=133, Status='unknown', Access='append')
      write(133,'(A)')  "Impossible to read " // trim(adjustl(FileName))
      write(*,'(A)')    "Impossible to read " // trim(adjustl(FileName))
    close(133)
  end if
    
    if (IntFlg == 1) then
      if (NLines==0) then
        read(Unit, *, iostat=Status)
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      read(Unit, *, iostat=Status)
      allocate(IntMatrix(NLines,4))
      IntMatrix=0
      allocate(RealMatrix(NLines,2))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit, *, iostat=Status) IntMatrix(iLine,1:3), RealMatrix(iLine,1:2), IntMatrix(iLine,4)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')    "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 2) then
      if (NLines==0) then
        read(Unit, *, iostat=Status)
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      read(Unit, *, iostat=Status)
      allocate(IntMatrix(NLines,4))
      IntMatrix=0
      allocate(RealMatrix(NLines,2))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit, *, iostat=Status) IntMatrix(iLine,1:2), RealMatrix(iLine,1)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')    "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 3) then
      if (NLines==0) then
        read(Unit, *, iostat=Status)
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      read(Unit, *, iostat=Status)
      allocate(IntMatrix(NLines,4))
      IntMatrix=0
      allocate(RealMatrix(NLines,2))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit, *, iostat=Status) IntMatrix(iLine,1:3), RealMatrix(iLine,1)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')    "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 4) then
      if (NLines==0) then
        read(Unit, *, iostat=Status)
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      read(Unit, *, iostat=Status)
      allocate(IntMatrix(NLines,0))
      IntMatrix=0
      allocate(RealMatrix(NLines,3))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit, *, iostat=Status) RealMatrix(iLine,1:3)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')    "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 5) then
      if (NLines==0) then
        read(Unit, *, iostat=Status)
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      read(Unit, *, iostat=Status)
      allocate(IntMatrix(NLines,1))
      IntMatrix=0
      allocate(RealMatrix(NLines,8))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit, *, iostat=Status) IntMatrix(iLine,1), RealMatrix(iLine,1:8)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')    "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 6) then
      if (NLines==0) then
        read(Unit, *, iostat=Status)
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      read(Unit, *, iostat=Status)
      allocate(IntMatrix(NLines,6))
      IntMatrix=0
      allocate(RealMatrix(NLines,2))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit, *, iostat=Status) IntMatrix(iLine,1:6), RealMatrix(iLine,1:2)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')   "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 7) then
      if (NLines==0) then
        do iTemp=1,5
          read(Unit, *, iostat=Status)
        end do
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      do iTemp=1,5
        read(Unit, *, iostat=Status)
      end do
      allocate(IntMatrix(NLines,1))
      IntMatrix=0
      allocate(RealMatrix(NLines,2))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit,'(20X,I20,2(d20.10))', iostat=Status) IntMatrix(iLine,1), RealMatrix(iLine,1:2)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')   "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    elseif (IntFlg == 8) then
      if (NLines==0) then
        do iTemp=1,5
          read(Unit, *, iostat=Status)
        end do
        do
          read(Unit, *, iostat=Status)
          if (Status/=0) exit
          NLines=NLines+1 
        end do
        rewind(Unit)
      end if
      do iTemp=1,5
        read(Unit, *, iostat=Status)
      end do
      allocate(IntMatrix(NLines,2))
      IntMatrix=0
      allocate(RealMatrix(NLines,2*Input%NTtra))
      RealMatrix=Zero
      do iLine = 1,NLines
        read(Unit,'(20X,I20,2(d20.10))', iostat=Status) IntMatrix(iLine,1), RealMatrix(iLine,1:2*Input%NTtra)
        if (Status/=0)then
          open( File='./ERROR', Unit=133, Status='unknown', Access='append')
            write(133,'(A)') "Impossible to read " // trim(adjustl(FileName))
            write(*,'(A)')   "Impossible to read " // trim(adjustl(FileName))
          close(133)
        end if
      end do
    end if
    
  close(Unit)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   COMPARING TEST MATRIXES TO BENCHMARK MATRIXES, COMPUTING RMS
! ==============================================================================================================
Pure Subroutine ComputingErrorsRegTest( IntFlg, IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )

  integer                                   ,intent(in)       ::    IntFlg
  integer   ,dimension(:,:)                 ,intent(inout)    ::    IntMatrix_Bench
  real(rkp) ,dimension(:,:)                 ,intent(inout)    ::    RealMatrix_Bench
  integer   ,dimension(:,:)                 ,intent(inout)    ::    IntMatrix_Test
  real(rkp) ,dimension(:,:)                 ,intent(inout)    ::    RealMatrix_Test
  real(rkp) ,dimension(:)   ,allocatable    ,intent(out)      ::    IntError
  real(rkp) ,dimension(:)   ,allocatable    ,intent(out)      ::    RealError
  
  integer                                                     ::    iCol, iLine
  real(rkp) ,dimension(:)   ,allocatable                      ::    IntMean_Bench, RealMean_Bench
  
  allocate(IntError(size(IntMatrix_Bench,2)))
  IntError  = Zero
  allocate(RealError(size(RealMatrix_Bench,2)))
  RealError = Zero
  
  allocate(IntMean_Bench(size(IntMatrix_Bench,2)))
  IntMean_Bench  = Zero
  allocate(RealMean_Bench(size(RealMatrix_Bench,2)))
  RealMean_Bench = Zero
  
  if (IntFlg /= 4) then
    do iCol = 1,size(IntMatrix_Bench,2)
      do iLine = 1,size(IntMatrix_Bench,1)
        IntError(iCol)      = IntError(iCol) + (IntMatrix_Bench(iLine,iCol) - IntMatrix_Test(iLine,iCol))**2
        IntMean_Bench(iCol) = IntMean_Bench(iCol) + real(IntMatrix_Bench(iLine,iCol),rkp)
      end do
    end do
    IntMean_Bench = IntMean_Bench / size(IntMatrix_Bench,1)
    IntError      = sqrt( IntError / size(IntMatrix_Bench,1) ) / IntMean_Bench
  end if
  
  do iCol = 1,size(RealMatrix_Bench,2)
    do iLine = 1,size(RealMatrix_Bench,1)
      RealError(iCol)      = RealError(iCol) + (RealMatrix_Bench(iLine,iCol) - RealMatrix_Test(iLine,iCol))**2
      RealMean_Bench(iCol) = RealMean_Bench(iCol) + RealMatrix_Bench(iLine,iCol)
    end do
  end do
  RealMean_Bench = RealMean_Bench / size(RealMatrix_Bench,1)
  RealError      = sqrt( RealError / size(RealMatrix_Bench,1) ) / RealMean_Bench
  
  deallocate(IntMean_Bench)
  deallocate(RealMean_Bench)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ==============================================================================================================
!   COMPARING TEST MATRIXES TO BENCHMARK MATRIXES, COMPUTING RMS
! ==============================================================================================================
Subroutine DeallocatingRegTest( IntMatrix_Bench, RealMatrix_Bench, IntMatrix_Test, RealMatrix_Test, IntError, RealError )

  integer   ,dimension(:,:) ,allocatable ,intent(inout)    ::    IntMatrix_Bench
  real(rkp) ,dimension(:,:) ,allocatable ,intent(inout)    ::    RealMatrix_Bench
  integer   ,dimension(:,:) ,allocatable ,intent(inout)    ::    IntMatrix_Test
  real(rkp) ,dimension(:,:) ,allocatable ,intent(inout)    ::    RealMatrix_Test
  real(rkp) ,dimension(:)   ,allocatable ,intent(inout)    ::    IntError
  real(rkp) ,dimension(:)   ,allocatable ,intent(inout)    ::    RealError
  
  deallocate(IntMatrix_Bench)
  deallocate(RealMatrix_Bench)
  deallocate(IntMatrix_Test)
  deallocate(RealMatrix_Test)
  
  deallocate(IntError)
  deallocate(RealError)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module
