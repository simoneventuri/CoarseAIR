! -*-F90-*-
!===============================================================================================================
! 
! Coarse-Grained QCT for Atmospheric Mixtures (CoarseAIR) 
! 
! Copyright (C) 2018 Berkan Bolkan and Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
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

Module O2_Basel_DiatomicPotential_Class

! The LeRoy model is used.

  use Parameters_Module         ,only:  rkp, Pi, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten
  use Logger_Class              ,only:  Logger
  use Error_Class               ,only:  Error
  use DiatomicPotential_Class   ,only:  DiatomicPotential_Type
  use RKHS            ! This module needs to be used by your code

  implicit none

  private
  public  ::    O2_Basel_DiatomicPotential_Type

  Type  ,extends(DiatomicPotential_Type)  ::    O2_Basel_DiatomicPotential_Type

    type(kernel)                           :: pes1              ! The kernel type is needed to set up and evaluate a RKHS model
    type(kernel)                           :: pes2              ! The kernel type is needed to set up and evaluate a RKHS model
    !type(kernel)                           :: pes3              ! The kernel type is needed to set up and evaluate a RKHS model

    logical                                :: stored   = .false.
    logical                                :: kread    = .false.
    logical                                :: ker1     = .false.
    logical                                :: ker2     = .false.
    
    integer                                :: na1
    integer                                :: na2
    integer                                :: nda1
    integer                                :: nda2
    real(rkp), allocatable, dimension(:,:) :: asy_array1
    real(rkp), allocatable, dimension(:,:) :: asy_array2
    !real(rkp), allocatable, dimension(:,:) :: asy_array3
    real(rkp), allocatable, dimension(:,:) :: darray1
    real(rkp), allocatable, dimension(:,:) :: darray2
  contains
    procedure          ::    Initialize        =>    Initialize_O2_Basel_DiatomicPotential
    procedure          ::    Compute_Vd_dVd    =>    Compute_Vd_dVd_O2_Basel
    procedure          ::    DiatomicPotential =>    DiatomicPotential_O2_Basel
    procedure ,private ::    ReadData_O2_Basel
  End Type

  character(*)                    ,parameter    :: Name_DiatPot    = 'Basel'

  real(rkp)                       ,parameter    :: dk26f1 = One / 14.0_rkp
  real(rkp)                       ,parameter    :: dk26f2 = One / 18.0_rkp
  real(rkp)                       ,parameter    :: dk24f1 = One / 15.0_rkp
  real(rkp)                       ,parameter    :: dk24f2 = Two / 21.0_rkp
  real(rkp)                       ,parameter    :: dk25f1 = Two / 21.0_rkp
  real(rkp)                       ,parameter    :: dk25f2 = One / 14.0_rkp
  real(rkp)                       ,parameter    :: akf1   = Two / Three
  real(rkp)                       ,parameter    :: lambda = 1e-10_rkp

  logical                         ,parameter    :: i_Debug_Global = .False.
 
  contains

!________________________________________________________________________________________________________________________________!
Subroutine Initialize_O2_Basel_DiatomicPotential( This, Input, SpeciesName, iMol, Mass1, Mass2, i_Debug )

  use Input_Class               ,only:  Input_Type

  class(O2_Basel_DiatomicPotential_Type),intent(out)  ::    This
  type(Input_Type)                      ,intent(in)   ::    Input
  character(:) ,allocatable             ,intent(in)   ::    SpeciesName
  integer                               ,intent(in)   ::    iMol
  real(rkp)                             ,intent(in)   ::    Mass1
  real(rkp)                             ,intent(in)   ::    Mass2
  logical ,optional                     ,intent(in)   ::    i_Debug
  
  logical                                             ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_O2_Basel_DiatomicPotential" )
  !i_Debug_Loc   =     Logger%On()
  
  allocate( This%Name        ,source = trim(Name_DiatPot) )
  allocate( This%SpeciesName ,source = trim(adjustl(SpeciesName)) )
  This%iMol         =    iMol
  This%Initialized  =   .True.

  This%RedMass = Mass1 * Mass2 / (Mass1 + Mass2)
  This%xmui    = One  / This%RedMass            ! Computing the inverse of the target reduced mass [1/a.u.]
  This%xmui2   = Half * This%xmui


  ! ==============================================================================================================
  !   READING VARIABLES
  ! ==============================================================================================================
  call This%ReadData_O2_Basel( Input, i_Debug=i_Debug_Loc )
  ! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadData_O2_Basel( This, Input, i_Debug )

  use Input_Class                         ,only:  Input_Type

  class(O2_Basel_DiatomicPotential_Type)  ,intent(inout)  ::    This
  type(Input_Type)                        ,intent(in)     ::    Input
  logical                       ,optional ,intent(in)     ::    i_Debug

  character(150)                                          ::    O2_Basel_Fldr
  integer                                                 ::    ii
  integer                                                 ::    Unit
  integer                                                 ::    Status
  logical                                                 ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadData_NO2_Basel_PES" )
  !i_Debug_Loc   =     Logger%On()

  O2_Basel_Fldr = trim(adjustl(Input%DtbPath))  // '/Molecules/O2/Basel/Params/'
  if (i_Debug_Loc) call Logger%Write( "Reading from Folder: ", trim(O2_Basel_Fldr) )


  if (.not. This%ker1) then
    inquire(file=trim(O2_Basel_Fldr)//"pes1.kernel", exist=This%ker1)   ! file_exists will be true if the file exists and false otherwise
  end if

  if (.not. This%ker2) then
    inquire(file=trim(O2_Basel_Fldr)//"pes2.kernel", exist=This%ker2)   ! file_exists will be true if the file exists and false otherwise
  end if

  !if (.not. ker3) then
  !  inquire(file=trim(O2_Basel_Fldr)//"pes23.kernel", exist=ker3)   ! file_exists will be true if the file exists and false otherwise
  !end if


  if (.not. This%stored ) then
    if (i_Debug_Loc) call Logger%Write( ".not. This%stored")

    if (i_Debug_Loc) call Logger%Write( "Reading from File: ", trim(O2_Basel_Fldr)//"asymp.dat" )
    open(NewUnit=Unit, file=trim(O2_Basel_Fldr)//"asymp.dat", status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // trim(O2_Basel_Fldr)//"asymp.dat" )

      read(Unit,*) This%na1
      if (i_Debug_Loc) call Logger%Write( "This%na1 = ", This%na1 )
      allocate( This%asy_array1(This%na1,2) )
      do ii = 1,This%na1
        read(Unit,*) This%asy_array1(ii,1), This%asy_array1(ii,2)
      end do

      read(Unit,*) This%na2
      if (i_Debug_Loc) call Logger%Write( "This%na2 = ", This%na2 )
      allocate( This%asy_array2(This%na2,2) )
      do ii = 1, This%na2
        read(Unit,*) This%asy_array2(ii,1), This%asy_array2(ii,2)
      end do

      !read(1001,*) This%na3
      !allocate( This%asy_array3(This%na3,2) )
      !do ii = 1, This%na3
      !  read(1001,*) This%asy_array3(ii,1), This%asy_array3(ii,2)
      !end do
      !This%na3 = This%na1
      !allocate( This%asy_array3(This%na3,2) )
      !This%asy_array3 = This%asy_array1

      read(Unit,*) This%nda1
      if (i_Debug_Loc) call Logger%Write( "This%nda1 = ", This%nda1 )
      allocate( This%darray1(This%nda1,2))
      do ii = 1,This%nda1
        read(Unit,*) This%darray1(ii,1), This%darray1(ii,2)
      end do

      read(Unit,*) This%nda2
      if (i_Debug_Loc) call Logger%Write( "This%nda2 = ", This%nda2 )
      allocate( This%darray2(This%nda2,2))
      do ii = 1,This%nda2
        read(Unit,*) This%darray2(ii,1), This%darray2(ii,2)
      end do

      This%stored = .true.

    close(Unit)

  end if


  if (.not. This%kread) then
  !if (This%ker1 .and. This%ker2 .and. This%ker3 ) then
    if (This%ker1 .and. This%ker2 ) then
      if (i_Debug_Loc) call Logger%Write( "Reading pes*.kernel Files" )

      call This%pes1%load_from_file(trim(O2_Basel_Fldr)//"pes1.kernel")
      call This%pes2%load_from_file(trim(O2_Basel_Fldr)//"pes2.kernel")
    !  call This%pes3%load_from_file(trim(O2_Basel_Fldr)//"pes23.kernel")
      This%kread = .true.
    else
      if (i_Debug_Loc) call Logger%Write( "Computing pes*.kernel Files" )

      call This%pes1%read_grid(trim(O2_Basel_Fldr)//"pes1.csv")
    !  print*,"IAMHERE"
      call This%pes1%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
      call This%pes1%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
      call This%pes1%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL) 
    !
      call This%pes1%calculate_coefficients_slow(lambda)
    !  call This%pes1%calculate_coefficients_fast()
    !
      call This%pes1%calculate_sums()
    !
      call This%pes1%save_to_file(trim(O2_Basel_Fldr)//"pes1.kernel")
    

      call This%pes2%read_grid(trim(O2_Basel_Fldr)//"pes2.csv")
    !  print*,"IAMHERE"
      call This%pes2%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
      call This%pes2%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
      call This%pes2%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
    !
      call This%pes2%calculate_coefficients_slow(lambda)
    !  call This%pes2%calculate_coefficients_fast()
    !
      call This%pes2%calculate_sums()
    !
      call This%pes2%save_to_file(trim(O2_Basel_Fldr)//"pes2.kernel")


    !  call This%pes3%read_grid(trim(O2_Basel_Fldr)//"pes3.csv")
    !  print*,"IAMHERE"
    !  call This%pes3%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
    !  call This%pes3%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
    !  call This%pes3%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
    !
    !  call This%pes3%calculate_coefficients_slow(lambda)
    !  call This%pes3%calculate_coefficients_fast()
    !
    !  call This%pes3%calculate_sums()
    !
    !  call This%pes3%save_to_file(trim(O2_Basel_Fldr)//"pes3.kernel")
      

      This%kread = .true.
    end if
  end if


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Elemental Function DiatomicPotential_O2_Basel( This, R ) result( V )

  class(O2_Basel_DiatomicPotential_Type),intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                                          :: V
    
  real(rkp)                                          :: VDiat
  real(rkp)                                          :: dVDiat

  call diato2(R, This%nda1, This%darray1, VDiat, dVDiat)
  
  V = VDiat                                                                                       
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!  


!________________________________________________________________________________________________________________________________!
Elemental Subroutine Compute_Vd_dVd_O2_Basel( This, R, V, dV )

  class(O2_Basel_DiatomicPotential_Type),intent(in)  :: This
  real(rkp)                             ,intent(in)  :: R
  real(rkp)                             ,intent(out) :: V
  real(rkp)                             ,intent(out) :: dV
  
  real(rkp)                                          :: VDiat
  real(rkp)                                          :: dVDiat
  
  call diato2(R, This%nda1, This%darray1, VDiat, dVDiat)
  
  V  =  VDiat                                                                                       
  dV = dVDiat 

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine diato2(r, nda1, darray1, ener, der)
  
  ! use surface2
  ! use rep_ker
  
  implicit none
  
  real(rkp)                 ,intent(in)  :: r
  integer                   ,intent(in)  :: nda1
  real(rkp) ,dimension(:,:) ,intent(in)  :: darray1
  real(rkp)                 ,intent(out) :: ener
  real(rkp)                 ,intent(out) :: der

  integer                                :: kk

  ener = Zero 
  der  = Zero
  do kk = 1,nda1
      ener = ener +  drker26(r,darray1(kk,1)) * darray1(kk,2)
      der  = der  + ddrker26(r,darray1(kk,1)) * darray1(kk,2)
  end do

  return

End Subroutine diato2
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function drker26(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: drker26, xl, xs

  xl = x
  xs = xi
  if (x .lt. xi) then
    xl = xi
    xs = x
  end if

  drker26 = dk26f1/xl**7 - dk26f2*xs/xl**8

End Function drker26
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function ddrker26(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: ddrker26, xl, xs

  if (x .lt. xi) then
    ddrker26 = -dk26f2/xi**8
  else
    ddrker26 = -Seven*dk26f1/x**8 + Eight*dk26f2*xi/x**9
  end if

End Function ddrker26
!--------------------------------------------------------------------------------------------------------------------------------!

      
End Module