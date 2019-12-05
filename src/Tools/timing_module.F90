Module Timing_Module

  use Parameters_Module     ,only:    rkp, zero

  implicit none

  public

  real(rkp)   ::    ti1   =   zero
  real(rkp)   ::    ti2   =   zero
  real(rkp)   ::    ti3   =   zero
  real(rkp)   ::    ti4   =   zero
  real(rkp)   ::    ti5   =   zero
  real(rkp)   ::    ti6   =   zero
  real(rkp)   ::    ti7   =   zero
  real(rkp)   ::    ti8   =   zero
  real(rkp)   ::    ti9   =   zero
  real(rkp)   ::    ti10  =   zero

  integer     ::    ic1   =   0
  integer     ::    ic2   =   0
  integer     ::    ic3   =   0   ! vphase
  integer     ::    ic4   =   0
  integer     ::    ic5   =   0
  integer     ::    ic6   =   0
  integer     ::    ic7   =   0
  integer     ::    ic8   =   0
  integer     ::    ic9   =   0
  integer     ::    ic10  =   0

!   real(rkp)     ::    t_Start     =     Zero
!   real(rkp)     ::    t_End       =     Zero
!   real(rkp)     ::    t_Start2    =     Zero
!   real(rkp)     ::    t_End2      =     Zero

  real(rkp)     ::    t_vphase    =     Zero  ! 3
  integer       ::    i_vphase    =     0

  real(rkp)     ::    t_inicon    =     Zero  ! 5
  integer       ::    i_inicon    =     0

  real(rkp)     ::    t_potcl     =     Zero  ! 6
  integer       ::    i_potcl     =     0

  real(rkp)     ::    t_pot       =     Zero  ! 7
  integer       ::    i_pot       =     0

  real(rkp)     ::    t_total     =     Zero

End Module