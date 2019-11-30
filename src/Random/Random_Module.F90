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

Module Random_Module
  
  use Parameters_Module       ,only:  rkp, pi, Two
  use Logger_Class            ,only:  Logger

  implicit none

  private
  public    ::    randset
  public    ::    randget
  public    ::    randdws
  public    ::    randdwsBIS
  public    ::    randnorm


  logical   ,parameter      ::    i_Debug_Global = .False.
  integer ,parameter        ::    ib    =    256
  integer ,dimension(8)     ::    iseed3
  real(rkp) ,dimension(100) ::    x

  contains

Subroutine randset( iseed, i_Debug )
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

  integer                                   ,intent(inout)  ::    iseed
  logical                         ,optional ,intent(in)     ::    i_Debug

  logical                                                   ::    i_Debug_Loc
  integer                                                   ::    is, i
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "randset")  !, Active = i_Debug_Loc )
  !i_Debug_Loc   =     Logger%On()

  if (i_Debug_Loc) then
    write(Logger%Unit,"(6x,'[randset]: Initiating a random number sequence')")
    write(Logger%Unit,"(6x,'[randset]: -> input iseed in base 10 = ',g0)") iseed
  end if

  ! To continue a sequence, set iseed to a negitive number. In this case the program will try to read in data from unit 50.
  if ( iseed < 0 ) then
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[randset]: Continuing random number sequence')")
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[randset]: -> reading from unit 50 for iseed3/x')")
    read(50) iseed3 , x
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[randset]: -> iseed3 (in base 2**8) = ',5x(1x,8i4))") (iseed3(9-i), i=1,8)

  ! convert iseed to base ib. answer is stored in iseed3
  else
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[randset]: Convert iseed to base 256')")
    is    =   iseed
    if ( mod(is,2) == 0 ) is = is-1
    do i = 1,8
      iseed3(i) =   mod(is,ib)
      is        =   is / ib
    end do
    if (i_Debug_Loc) write(Logger%Unit,"(6x,'[randset]: iseed3 = 2 * iseed + 1   in base 256 (i.e. base 2**8) = ',8i4)") (iseed3(9-i), i=1,8)

    ! exercise random number generator in case iseed not choosen between 0 and (2**31 - 1). put 100 random numbers into array x: this is new sequence
    do i = 1,size(x)
      x(i)    =   rand12(iseed3)
    end do
    iseed=0
  end if

  if (i_Debug_Loc) call Logger%Exiting

End Subroutine


Subroutine randget()

  integer     ::      i

  rewind 50
  write(50) iseed3, x

  if (i_Debug_Global) write(Logger%Unit,"(4x,'[randget]: save random no. seq. on unit 50 - next seed is',5x,8i4)") (iseed3(9-i), i=1,8)

End Subroutine



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
Function randdws() result(RandNum)
  real(8)                   ::    RandNum
  integer                   ::    i
!   real(8) ,save :: isave = 0.01_8
!   isave   =  isave + 0.1_8
!   if (isave>=0.999_8) isave=0.05_8
!   RandNum = isave! Bruno
!   RandNum = 0.5_8! Bruno
!   return
!   RandNum = 0.5d0! Bruno
!   return
  i         =   int( 99d0 * x(100) ) + 1
  RandNum   =   x(100)
  x(100)    =   x(i)
  x(i)      =   rand12(iseed3)
End Function


! function rand1(iseed)
! random number generatation using the multipicative congruential
! method.

! written by david schwenke (nasa ames)
! for machines with a maximum integer of (2**32 - 1) or larger
! and has precision of 1 in 2**64

! amended by kieran f. lim (stanford university)
! for machines with a maximum integer of (2**16 - 1) or larger
! note:  dec microvax has a maximum integer of (2**31 - 1)

! this routine hits integer overflow only if the maximum
! integer is (2**16 -1) = 65535 or smaller

! if x(n) is the current seed, then the next one x(n+1) is given by

! x(n+1)=(a*x(n)+c)mod(m)
!                         and the random number is x(n+1)/m.

! the arithmatic is carried out using 8 "digit" numbers
! of base 2**8 (=256), thus iseed is 8 element array with
! each element a "digit". the units "digit" is
! the first element of iseed.
! a ninth "digit" is provided for overflows in the intermediate
! steps for arithmetic operations

! parameters:
!  m = 2**64 = (1) (0) (0) (0) (0) (0) (0) (0) (0) base 2**8
!  a = 6,364,136,223,846,793,005 base 10
!    = (22609) (62509) (19605) (322557) base 2**16
!    = (88) (81) (244) (45) (76) (149) (127) (45) base 256
!  c = 0

! these parameters are from d.e. knuth, "the art of computer
! programming", vol. 2 of "seminumerical algorithms" (addison-wesley
! reading, mass. 1981) 2nd edition,
! section 3.3.4, table 1 at pages 102-104
! and are atributed to c.e. haynes
Function rand12( iseed )

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


End Module