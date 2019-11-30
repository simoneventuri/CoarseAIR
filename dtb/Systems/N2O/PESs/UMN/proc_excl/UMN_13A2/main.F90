
!-----------------------------------------------------------------------------!
!> This program prints out the needed code lines for the triplet A'' N2-O  
!! based on those for the complete surface N2-O2 PES 
  program main 

    implicit none

    integer, parameter                    :: hpc_4    = 4
    integer(hpc_4), parameter             :: in_file  = 10
    integer(hpc_4), parameter             :: out_file = in_file + 1
    integer(hpc_4)                        :: c, i, j, k, kp, ios
    integer(hpc_4), dimension(0:2253)     :: pos_eq
    integer(hpc_4), dimension(3)          :: n_dpdr = [150, 150, 125]
    character(len=1)                      :: kp_char
    character(len=120), dimension(1:161)  :: p_new
    character(len=120), dimension(1:150)  :: dpdr_new
    character(len=120), dimension(0:2253) :: p_orig
    character(len=120), dimension(0:2253) :: dpdr_orig
    logical, dimension(1:161)             :: found 

    !-----------------------!
    ! Read original code listing 
    print*
    write(*,5)'Reading original coefficients for POTENTIAL...'

    open(unit=in_file,file='input/orig.dat',status='old',action='read',iostat=ios)

    if (ios.ne.0) stop "ERROR: file 'orig.dat' not found!"

    do i = 0,2253
       read(in_file,5)p_orig(i)
       pos_eq(i) = index(trim(p_orig(i)), '=')
    enddo

    close(in_file)

    write(*,5)'Done!'

    !-----------------------!
    ! Read modified code listing (formulas are attached later)
    print*
    write(*,5)'Reading remaining coefficients (formula is appended later)..'

    open(unit=in_file,file='input/mod.dat',status='old',action='read',iostat=ios)

    if (ios.ne.0) stop "ERROR: file 'red.dat' not found!"

    do i = 1,161
       read(in_file,5)p_new(i)
    enddo

    close(in_file)

    write(*,5)'Done!'

    found = .false.

    open(unit=out_file,file='output/new.dat',status='unknown',action='write')

    c = 0
    do i = 1,161
       found(i) = .false.
       search_loop : do j = 0,2253
          if (trim(p_orig(j)(1:pos_eq(j) - 1)).eq.trim(p_new(i))) then
             write(out_file,5)trim(p_orig(j))
             c = c + 1
             exit search_loop
          endif
       enddo search_loop
    enddo

    close(out_file)

    !-----------------------!
    ! Re-initialization 
    print*
    write(*,5)'Reading original coefficients for FORCE...'

    open(unit=in_file,file='input/orig_der.dat',status='old',action='read',iostat=ios)

    if (ios.ne.0) stop "ERROR: file 'orig_der.dat' not found!"

    pos_eq = -1
    do i = 0,2253
       read(in_file,5)dpdr_orig(i)
       pos_eq(i) = index(trim(dpdr_orig(i)), '=')
    enddo

    close(in_file)

    ! Read non-zero entries i = 4
    do k = 1,3

      kp = k + 3
      write(kp_char,'(i1)')kp
      kp_char = adjustl(kp_char)

      open(unit=in_file,file='input/mod_der'//trim(kp_char)//'.dat',status='old',action='read',iostat=ios)

      if (ios.ne.0) stop "ERROR: file 'mod_der4.dat' not found!"

      do i = 1,n_dpdr(k)
         read(in_file,5)dpdr_new(i)
         !print*,trim(dpdr_new(i))
      enddo

      close(in_file)

      write(*,5)'Done!'

      found = .false.
    
      open(unit=out_file,file='output/new_der'//trim(kp_char)//'.dat',status='unknown',action='write')

      c = 0
      do i = 1,n_dpdr(k)
         found(i) = .false.
         search_loop2 : do j = 0,2253
           if (trim(dpdr_orig(j)(1:pos_eq(j) - 1)).eq.trim(dpdr_new(i))) then
               write(out_file,5)trim(dpdr_orig(j))
               c = c + 1
               exit search_loop2
           endif
         enddo search_loop2
      enddo

      close(out_file)

    enddo

    write(*,5)'Done!'

5 format(a)

  end program main
!-----------------------------------------------------------------------------!
