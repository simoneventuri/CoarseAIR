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

Module Matrix_Inversion_Module

  implicit none

  private
  public    ::    Inverse_Matrix
!   public    ::    SolveLinearSystem
!
!   Interface   SolveLinearSystem
!     Module Procedure  SolveLinearSystem
!     Module Procedure  SolveLinearSystem_Several_RHS
!   End Interface

  contains

! Returns the inverse of a matrix calculated by finding the LU decomposition.  Depends on LAPACK.
! Downloaded at: http://fortranwiki.org/fortran/show/inv
Function Inverse_Matrix( A ) result( Ainv )

  real(8)   ,dimension(:,:)                             ,intent(in)     ::    A                               !< Matrix to be inverted
  real(8)   ,dimension(size(A,1),size(A,2))                             ::    Ainv                            !< Invserse matrix

  real(8)   ,dimension(size(A,1))                                       ::    Work                            ! Work array for LAPACK
  integer   ,dimension(size(A,1))                                       ::    ipiv                            ! Pivot indices
  integer                                                               ::    N
  integer                                                               ::    Info

! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  Ainv    =     A                                                                                               ! Store A in Ainv to prevent it from being overwritten by LAPACK
  N       =     size(A,1)

! ==============================================================================================================
!    COMPUTING AN LU FACTORIZATION USING PARTIAL PIVOTING
! ==============================================================================================================
! This section computes an LU factorization of a general M-by-N matrix A using partial pivoting with row
! interchanges. It uses the 'DGETRF' procedure from Lapack.
! ==============================================================================================================
  call DGETRF( n, n, Ainv, n, ipiv, Info )
  if (Info /= 0) then
     stop 'Matrix is numerically singular!'
  end if
! ==============================================================================================================


! ==============================================================================================================
!    COMPUTING THE INVERSE OF A MATRIX USING THE LU FACTORIZATION
! ==============================================================================================================
! This section computes the inverse of a matrix using the LU factorization computed by DGETRF.
! It uses the 'DGETRI' procedure from Lapack.
! ==============================================================================================================
  call DGETRI(n, Ainv, n, ipiv, Work, n, Info)
  if (Info /= 0) then
     stop 'Matrix inversion failed!'
  end if
! ==============================================================================================================

End Function
!
!
! ! Solve the linear system A x = b
! Function SolveLinearSystem( A, b ) result( x )
!
!   real(8)   ,dimension(:,:)     ,intent(in)     ::    A
!   real(8)   ,dimension(:)       ,intent(in)     ::    b
!   real(8)   ,dimension( size(b) )               ::    x
!
!
! !   character(1)    ::    TRANS
! !   integer         ::    M, N, NRHS
! !   integer                 ,parameter     ::    NB    =   128                   ! The optimum block size
! !   real(rkp)   ,dimension(:)     ,allocatable                            ::    WORK                            !
! !   integer                                                               ::    LWORK                           !
!
!   integer                                                               ::    INFO                            ! Indicator
!   integer                                                               ::    M                               ! The number of rows of the matrix A.  M >= 0.
!   integer                                                               ::    N                               ! The number of columns of the matrix A.  N >= 0.
!   integer                                                               ::    NRHS                            ! The number of right hand sides, i.e., the number of columns of the matrices B and X. NRHS >=0.
!   integer                                                               ::    LDA                             ! The leading dimension of the array A.  LDA >= max(1,M).
!   integer                                                               ::    LDB                             ! The leading dimension of the array B. LDB >= MAX(1,M,N).
!   integer                                                               ::    LDC
!   integer                                                               ::    MN                              !
!   integer                                                               ::    LWORK                           !
!   integer                 ,parameter                                    ::    NB    =   128                   ! The optimum block size
!   character(len=1)        ,parameter                                    ::    TRANS =   'N'                   ! 'N': the linear system involves A; 'T': the linear system involves A**T.
! !   real(rkp)   ,dimension(:)     ,allocatable                            ::    b                               !
! !   real(rkp)   ,dimension(:)     ,allocatable                            ::    d                               !
! !   real(rkp)   ,dimension(:,:)   ,allocatable                            ::    A                               !
! !   real(rkp)   ,dimension(:,:)   ,allocatable                            ::    C                               ! On entry, the P-by-N constrain matrix C. On exit, the upper triangle of the subarray (1:P,N-P+1:N) contains the P-by-P upper triangular matrix R.
!   real(8)   ,dimension(:)     ,allocatable                            ::    WORK                            !
!
!
!   real(8)   ,dimension( size(A,1), size(A,2) )      ::    A_
!   real(8)   ,dimension( size(b,1),1 )               ::    b_
!
!   A_      =   A
!   b_(:,1) =   b
!
!   M       =     size(A,1)
!   N       =     size(A,2)
! !   NRHS    =     1
!
!   LDA       =       max(1,M)                                                                                    ! The leading dimension of the array A.  LDA >= max(1,M).
!   LDB       =       max(1,M,N)                                                                                  ! The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!
!     NRHS      =       1                                                                                           ! The number of right hand sides, i.e., the number of columns of the matrices B and X. NRHS >=0.
!     MN        =       min(M,N)                                                                                    ! Minimum value of M and N
!     LWORK     =       max( 1, MN + max(MN,NRHS) * NB )                                                            ! The dimension of the array WORK
!     allocate( WORK(LWORK) )                                                                                       ! Allocating array
!
!     call DGELS( TRANS, M, N, NRHS, A, LDA, b, LDB, WORK, -1, INFO )                                               !
!
!     LWORK     =       int( WORK(1) )
!
!     if ( LWORK > 0 ) then
!       deallocate( WORK ); allocate( WORK(LWORK) )                                                                 ! Allocating array
!     else
!       LWORK   =       size(WORK)
!     end if
!
!
!     call DGELS( TRANS, M, N, NRHS, A_, LDA, b_, LDB, WORK, LWORK, INFO )                                            ! Fitting the data: Solve Ax = b using LAPACK
! !     if (i_Debug_Loc) then
! !       call Logger%Write( " -> INFO = ", Info )
! !       select case (INFO)
! !         case(0);    call Logger%Write( " -> Successful exit." )
! !         case default;  call Logger%Write( " -> The Interpolate fit matrix is singular!" )
! !       end select
! !     end if
! !     if ( INFO /= 0 ) call Error%Raise( "Error in DGELS", ProcName=ProcName )
!     x  =       b_(1:N, 1)                                                                                      ! Getting the solution vector
! ! ==============================================================================================================
!
! End Function
!
!
!
! ! Solve the linear system A x = b
! Function SolveLinearSystem_Several_RHS( A, b ) result( x )
!
!   real(8)   ,dimension(:,:)     ,intent(in)     ::    A
!   real(8)   ,dimension(:,:)     ,intent(in)     ::    b
!   real(8)   ,dimension( size(b,1), size(b,2) )  ::    x
!
!
! !   character(1)    ::    TRANS
! !   integer         ::    M, N, NRHS
! !   integer                 ,parameter     ::    NB    =   128                   ! The optimum block size
! !   real(rkp)   ,dimension(:)     ,allocatable                            ::    WORK                            !
! !   integer                                                               ::    LWORK                           !
!
!   integer                                                               ::    INFO                            ! Indicator
!   integer                                                               ::    M                               ! The number of rows of the matrix A.  M >= 0.
!   integer                                                               ::    N                               ! The number of columns of the matrix A.  N >= 0.
!   integer                                                               ::    NRHS                            ! The number of right hand sides, i.e., the number of columns of the matrices B and X. NRHS >=0.
!   integer                                                               ::    LDA                             ! The leading dimension of the array A.  LDA >= max(1,M).
!   integer                                                               ::    LDB                             ! The leading dimension of the array B. LDB >= MAX(1,M,N).
!   integer                                                               ::    MN                              !
!   integer                                                               ::    LWORK                           !
!   integer                 ,parameter                                    ::    NB    =   128                   ! The optimum block size
!   character(len=1)        ,parameter                                    ::    TRANS =   'N'                   ! 'N': the linear system involves A; 'T': the linear system involves A**T.
! !   real(rkp)   ,dimension(:)     ,allocatable                            ::    b                               !
! !   real(rkp)   ,dimension(:)     ,allocatable                            ::    d                               !
! !   real(rkp)   ,dimension(:,:)   ,allocatable                            ::    A                               !
! !   real(rkp)   ,dimension(:,:)   ,allocatable                            ::    C                               ! On entry, the P-by-N constrain matrix C. On exit, the upper triangle of the subarray (1:P,N-P+1:N) contains the P-by-P upper triangular matrix R.
!   real(8)   ,dimension(:)     ,allocatable                            ::    WORK                            !
!
!
!   real(8)   ,dimension( size(A,1), size(A,2) )      ::    A_
!   real(8)   ,dimension( size(b,1), size(b,2) )      ::    b_
!
!   A_      =   A
!   b_      =   b
!
!   M       =     size(A,1)
!   N       =     size(A,2)
!   NRHS    =     size(b,2)
!
!   LDA       =       max(1,M)                                                                                    ! The leading dimension of the array A.  LDA >= max(1,M).
!   LDB       =       max(1,M,N)                                                                                  ! The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!
!     MN        =       min(M,N)                                                                                    ! Minimum value of M and N
!     LWORK     =       max( 1, MN + max(MN,NRHS) * NB )                                                            ! The dimension of the array WORK
!     allocate( WORK(LWORK) )                                                                                       ! Allocating array
!
!     call DGELS( TRANS, M, N, NRHS, A, LDA, b, LDB, WORK, -1, INFO )                                               !
!
!     LWORK     =       int( WORK(1) )
!
!     if ( LWORK > 0 ) then
!       deallocate( WORK ); allocate( WORK(LWORK) )                                                                 ! Allocating array
!     else
!       LWORK   =       size(WORK)
!     end if
!
!
!     call DGELS( TRANS, M, N, NRHS, A_, LDA, b_, LDB, WORK, LWORK, INFO )                                            ! Fitting the data: Solve Ax = b using LAPACK
! !     if (i_Debug_Loc) then
! !       call Logger%Write( " -> INFO = ", Info )
! !       select case (INFO)
! !         case(0);    call Logger%Write( " -> Successful exit." )
! !         case default;  call Logger%Write( " -> The Interpolate fit matrix is singular!" )
! !       end select
! !     end if
! !     if ( INFO /= 0 ) call Error%Raise( "Error in DGELS", ProcName=ProcName )
!     x  =       b_                                                                                      ! Getting the solution vector
! ! ==============================================================================================================
!
! End Function

End Module
