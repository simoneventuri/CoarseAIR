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

Module RandomVector_Module
  
  use Parameters_Module       ,only:  rkp, pi, Two
  use Error_Class             ,only:  Error
  use Logger_Class            ,only:  Logger

  implicit none

  private
  public    ::    RandsetVec
  public    ::    RanddwsVec
  public    ::    RandNorm
  
  integer   ,dimension(:,:) ,allocatable ::    iSeedStoreVec ! ~ iseed3
  real(rkp) ,dimension(:,:) ,allocatable ::    xMat          ! ~ x

  integer   ,parameter                   ::    ib             =    256
  logical   ,parameter                   ::    i_Debug_Global = .False.
  
  contains


Subroutine RandsetVec( iseed, i_Debug )
!The subroutine sets up data for random number generator.

!Written by david schwenke (nasa ames) for machines with a maximum integer of (2**32 - 1) or larger and has precision of 1 in 2**64,
!amended by kieran f. lim (stanford university) for machines with a maximum integer of (2**16 - 1) or larger
!note:  dec microvax has a maximum integer of (2**31 - 1)
!latest update: 12 april 1990
!this routine hits integer overflow only if the maximum integer is (2**16 -1) = 65535 or smaller
!either read in existing sequence of random numbers from unit 50
!or start a new sequence:
! 1.  convert iseed to base ib  (ib=256=2**8)
!     0  <=  iseed  <= (2**31 - 1)
! 2.  take iseed3=2*iseed+1
!     1  <=  iseed3  <= (2**32 - 1)
!     with iseed3 odd
! 3.  exercise random number generator in case
!     iseed not choosen between 0 and (2**31 - 1)
! 4.  put 100 random numbers into array x: this is
!     new sequence
!random number generator is function rand1 which uses the multipicative congruential method
!reference: d.e. knuth, "the art of computer programming", vol. 2 of "seminumerical algorithms" (addison-wesley, 1981), 2nd edition.
!calculations in base (ib=256=2**8)

  integer ,dimension(:) ,intent(inout)  ::    iSeed
  logical     ,optional ,intent(in)     ::    i_Debug

  integer                               ::    Status
  integer                               ::    is, i, iLine, VecSize
  logical                               ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "RandsetVec")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()
  
  
  VecSize = size(iSeed,1)
  
  allocate( iSeedStoreVec(VecSize,8), stat=Status )
  if (Status/=0) call Error( "Error allocating iSeedStoreVec" )
  if (i_Debug_Loc) call Logger%Write( "Allocated iSeedStoreVec with dimensions = (", VecSize, ",8)")

  allocate( xMat(VecSize,100), stat=Status )
  if (Status/=0) call Error( "Error allocating xMat" )
  if (i_Debug_Loc) call Logger%Write( "Allocated xMat with dimensions = (", VecSize, ",100)")


  if (i_Debug_Loc) write(Logger%Unit,"(6x,'[RandsetVec]: Initiating a Random Nb Sequence')")

  do iLine = 1,size(iSeed)
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[RandsetVec]: -> Input iseed in base 10 = ',g0)") iSeed(iLine)

    ! convert iSeed to base ib. answer is stored in iSeed3
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[randset]: Convert iseed to base 256')")
    is = iSeed(iLine)
    if ( mod(is,2) == 0 ) is = is-1
    do i = 1,8
      iSeedStoreVec(iLine,i) =   mod(is,ib)
      is                     =   is / ib
    end do
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[RandsetVec]: iseed3 = 2 * iseed + 1   in base 256 (i.e. base 2**8) = ',8i4)") (iSeedStoreVec(iLine,9-i), i=1,8)

    ! exercise random number generator in case iSeed not choosen between 0 and (2**31 - 1). put 100 random numbers into matrix xMat: this is new sequence
    do i = 1,size(xMat,2)
      xMat(iLine,i) = Rand12(iSeedStoreVec(iLine,:))
    end do
    
    iSeed(iLine) = 0
  end do
  

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Function RanddwsVec(iLine) result(RandNum)
! function rand(idum)
! gets the next random number from a sequence of random numbers
! by shuffling the sequence
!
! reference: d.e. knuth, "the art of computer programming",
! vol. 2 of "seminumerical algorithms" (addison-wesley, 1981),
! 2nd edition. (pg. 32)
!
! written by david schwenke (nasa ames)
! amended by kieran f. lim (stanford university)
!
! idum is a dummy number

  integer                  ,intent(in)  ::    iLine
  
  real(8)                               ::    RandNum
  integer                               ::    i
  logical                               ::    i_Debug_Loc


!   real(8) ,save :: isave = 0.01_8
!   isave   =  isave + 0.1_8
!   if (isave>=0.999_8) isave=0.05_8
!   RandNum = isave! Bruno
!   RandNum = 0.5_8! Bruno
!   return
!   RandNum = 0.5d0! Bruno
!   return

  !i_Debug_Loc = .True.!i_Debug_Global; 
  !!if ( present(i_Debug) )i_Debug_Loc = i_Debug
  !if (i_Debug_Loc) call Logger%Entering( "RanddwsVec")  !, Active = i_Debug_Loc )
  !!i_Debug_Loc   =     Logger%On()

  i               =   int( 99d0 * xMat(iLine,100) ) + 1
  RandNum         =   xMat(iLine,100)
  xMat(iLine,100) =   xMat(iLine,i)
  xMat(iLine,i)   =   Rand12(iSeedStoreVec(iLine,:))

  !if (i_Debug_Loc) call Logger%Write( " RandNum = ", RandNum )
  !if (i_Debug_Loc) call Logger%Exiting
  
End Function


Function rand12( iseed )

! Function rand1(iseed) random number generatation using the multipicative congruential method.
! Written by david schwenke (NASA Ames) for machines with a maximum integer of (2**32 - 1) or larger and has precision of 1 in 2**64
! Amended by kieran f. lim (stanford university) for machines with a maximum integer of (2**16 - 1) or larger
! NOTE:  dec microvax has a maximum integer of (2**31 - 1)
!
! This routine hits integer overflow only if the maximum integer is (2**16 -1) = 65535 or smaller
!
! If x(n) is the current seed, then the next one x(n+1) is given by x(n+1)=(a*x(n)+c)mod(m) and the random number is x(n+1)/m.
!
! The arithmatic is carried out using 8 "digit" numbers of base 2**8 (=256), thus iseed is 8 element array with each element a "digit". 
! The units "digit" is the first element of iseed. 
! A ninth "digit" is provided for overflows in the intermediate steps for arithmetic operations
!
! parameters:
!  m = 2**64 = (1) (0) (0) (0) (0) (0) (0) (0) (0) base 2**8
!  a = 6,364,136,223,846,793,005 base 10
!    = (22609) (62509) (19605) (322557) base 2**16
!    = (88) (81) (244) (45) (76) (149) (127) (45) base 256
!  c = 0
!
! These parameters are from d.e. knuth, "the art of computer programming", vol. 2 of "seminumerical algorithms" (addison-wesley
! reading, mass. 1981) 2nd edition, section 3.3.4, table 1 at pages 102-104 and are atributed to c.e. haynes

  integer ,dimension(8)                 ,intent(inout)  ::    iseed
  real(rkp)                                             ::    rand12

  integer ,dimension(8)                     ,parameter  ::    ia = [45 , 127, 149, 76 , 45 , 244, 81 , 88 ]   ! base 2**8 information
  integer ,dimension(8)                     ,parameter  ::    ic = [0  , 0  , 0  , 0  , 0  , 0  , 0   , 0 ]   ! base 2**8 information
  real(rkp)                                 ,parameter  ::    bi = 3.90625E-3_rkp
  integer                                               ::    i, j, k, ip
  integer ,dimension(16)                                ::    id

  do i = 1,8                              ! !     id will equal iseed*ia+ic. set it equal to ic
    id(i)     =   ic(i)
    id(i+8)   =   0
  end do

! Form ia*iseed+ic
  do j = 1,8
    do i = 1,9-j
      k       =   j + i - 1
      ip      =   ia(j) * iseed(i)        ! ip is unnormalized product of k digit of result so far. It should never exceed ib*(ib-1).
      do                                  ! Loop to normalize the result (carry forword to next digits if necessary)
        ip    =   ip + id(k)
        id(k) =   mod(ip,ib)
        ip    =   ip / ib
        if ( (ip==0) .or. (k==8) ) exit
        k     =   k + 1
      end do
    end do
  end do
  iseed(1:8)  =   id(1:8)

! now determine dfloating point random number
  rand12      =   dfloat(iseed(1))
  do i = 2,8
    rand12    =   dfloat(iseed(i)) + rand12 * bi
  end do

  rand12      =   rand12 * bi

End Function


Function RandNorm(mean, sd) result( RandN )

  real(rkp)                   ,intent(in) :: mean
  real(rkp)                   ,intent(in) :: sd
  real(rkp)                               :: RandN
 
  real(rkp)    ,dimension(2)              :: array
  integer      ,dimension(64)             :: seed
  integer                                 :: i,j
  
  !call random_init(.true., .false.)
  !call random_seed(put = seed)
  call RANDOM_NUMBER(array)
  
  i = 1
  RandN      = mean + sd * sqrt(-Two * log(array(i))) * sin(Two * pi * array(i+1))
  !write(*,*) 'mean = ', mean, ', sd = ', sd, ', randn = ', randn 

End Function


End Module
