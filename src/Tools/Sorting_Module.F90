Module Sorting_Module

  use Parameters_Module           ,only:  rkp

  implicit none

  private
  public   ::   hpsort

  Interface             hpsort
    Module procedure    hpsort_R8
    Module procedure    hpsort_I4
    Module procedure    hpsort_I8
  End Interface

  contains

! Heapsort algorithm

Subroutine hpsort_R8( ra, Idx )

  real(rkp)   ,dimension(:)             ,intent(inout)  ::    ra
  integer     ,dimension(:)   ,optional ,intent(out)    ::    Idx

  real(rkp)                                             ::    rra
  integer                                               ::    i, ir, j, k, irra, n
  integer     ,dimension(size(ra))                      ::    Idx_

  n         =   size(ra)
  do i = 1,n
    Idx_(i)  =   i
  end do

  if (n < 2) then
    if ( present(Idx) ) Idx = Idx_
    return
  end if

  k     =   n/2+1
  ir    =   n
  do
    if ( k > 1 ) then
      k         =   k-1
      rra       =   ra(k)
      irra      =   Idx_(k)
    else
      rra       =   ra(ir)
      irra      =   Idx_(ir)
      ra(ir)    =   ra(1)
      Idx_(ir)  =   Idx_(1)
      ir        =   ir-1
      if ( ir == 1 ) then
        ra(1)   =   rra
        Idx_(1) =   irra
        exit
      end if
    end if

    i           =   k
    j           =   k+k
    do
      if ( j > ir ) exit
      if ( j < ir ) then
        if ( ra(j) < ra(j+1))j=j+1
      end if
      if ( rra < ra(j) ) then
        ra(i)   =   ra(j)
        Idx_(i) =   Idx_(j)
        i=j
        j=j+j
      else
        j=ir+1
      end if
    end do

    ra(i)       =   rra
    Idx_(i)     =   irra

  end do

  if ( present(Idx) ) Idx = Idx_

End Subroutine


! Subroutine hpsort_R8I8( ra, Idx )

!   real(rkp)      ,dimension(:)             ,intent(inout)  ::    ra
!   integer(rkp)   ,dimension(:)   ,optional ,intent(out)    ::    Idx

!   real(rkp)                                                ::    rra
!   integer(rkp)                                             ::    i, ir, j, k, irra, n
!   integer(rkp)   ,dimension(size(ra))                      ::    Idx_

!   n         =   size(ra)
!   do i = 1,n
!     Idx_(i)  =   i
!   end do

!   if (n < 2) then
!     if ( present(Idx) ) Idx = Idx_
!     return
!   end if

!   k     =   n/2+1
!   ir    =   n
!   do
!     if ( k > 1 ) then
!       k         =   k-1
!       rra       =   ra(k)
!       irra      =   Idx_(k)
!     else
!       rra       =   ra(ir)
!       irra      =   Idx_(ir)
!       ra(ir)    =   ra(1)
!       Idx_(ir)  =   Idx_(1)
!       ir        =   ir-1
!       if ( ir == 1 ) then
!         ra(1)   =   rra
!         Idx_(1) =   irra
!         exit
!       end if
!     end if

!     i           =   k
!     j           =   k+k
!     do
!       if ( j > ir ) exit
!       if ( j < ir ) then
!         if ( ra(j) < ra(j+1))j=j+1
!       end if
!       if ( rra < ra(j) ) then
!         ra(i)   =   ra(j)
!         Idx_(i) =   Idx_(j)
!         i=j
!         j=j+j
!       else
!         j=ir+1
!       end if
!     end do

!     ra(i)       =   rra
!     Idx_(i)     =   irra

!   end do

!   if ( present(Idx) ) Idx = Idx_

! End Subroutine


Subroutine hpsort_I4( ra, Idx )

  integer   ,dimension(:)             ,intent(inout)  ::    ra
  integer   ,dimension(:)   ,optional ,intent(out)    ::    Idx

  integer                                             ::    rra
  integer                                             ::    i, ir, j, k, irra, n
  integer   ,dimension(size(ra))                      ::    Idx_

  n         =   size(ra)
  do i = 1,n
    Idx_(i)  =   i
  end do

  if (n < 2) then
    if ( present(Idx) ) Idx = Idx_
    return
  end if
  n         =   size(ra)
  do i = 1,n
    Idx_(i)  =   i
  end do

  if (n < 2) return

  k     =   n/2+1
  ir    =   n
  do
    if ( k > 1 ) then
      k         =   k-1
      rra       =   ra(k)
      irra      =   Idx_(k)
    else
      rra       =   ra(ir)
      irra      =   Idx_(ir)
      ra(ir)    =   ra(1)
      Idx_(ir)  =   Idx_(1)
      ir        =   ir-1
      if ( ir == 1 ) then
        ra(1)   =   rra
        Idx_(1) =   irra
        exit
      end if
    end if

    i           =   k
    j           =   k+k
    do
      if ( j > ir ) exit
      if ( j < ir ) then
        if ( ra(j) < ra(j+1))j=j+1
      end if
      if ( rra < ra(j) ) then
        ra(i)   =   ra(j)
        Idx_(i) =   Idx_(j)
        i=j
        j=j+j
      else
        j=ir+1
      end if
    end do

    ra(i)       =   rra
    Idx_(i)     =   irra

  end do

  if ( present(Idx) ) Idx = Idx_

End Subroutine


Subroutine hpsort_I8( ra, Idx )

  integer(rkp)   ,dimension(:)             ,intent(inout)  ::    ra
  integer(rkp)   ,dimension(:)   ,optional ,intent(out)    ::    Idx

  integer(rkp)                                             ::    rra
  integer(rkp)                                             ::    i, ir, j, k, irra, n
  integer(rkp)   ,dimension(size(ra))                      ::    Idx_

  n         =   size(ra)
  do i = 1,n
    Idx_(i)  =   i
  end do

  if (n < 2) then
    if ( present(Idx) ) Idx = Idx_
    return
  end if
  n         =   size(ra)
  do i = 1,n
    Idx_(i)  =   i
  end do

  if (n < 2) return

  k     =   n/2+1
  ir    =   n
  do
    if ( k > 1 ) then
      k         =   k-1
      rra       =   ra(k)
      irra      =   Idx_(k)
    else
      rra       =   ra(ir)
      irra      =   Idx_(ir)
      ra(ir)    =   ra(1)
      Idx_(ir)  =   Idx_(1)
      ir        =   ir-1
      if ( ir == 1 ) then
        ra(1)   =   rra
        Idx_(1) =   irra
        exit
      end if
    end if

    i           =   k
    j           =   k+k
    do
      if ( j > ir ) exit
      if ( j < ir ) then
        if ( ra(j) < ra(j+1))j=j+1
      end if
      if ( rra < ra(j) ) then
        ra(i)   =   ra(j)
        Idx_(i) =   Idx_(j)
        i=j
        j=j+j
      else
        j=ir+1
      end if
    end do

    ra(i)       =   rra
    Idx_(i)     =   irra

  end do

  if ( present(Idx) ) Idx = Idx_

End Subroutine


End Module