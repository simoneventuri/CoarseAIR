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

Module NO2_Basel_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Pi, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use RKHS            ! This module needs to be used by your code

  implicit none

  private
  public    ::    NO2_Basel_PES_Type


  Type    ,extends(PES_Type)          ::    NO2_Basel_PES_Type

    type(kernel)                           :: pes1              ! The kernel type is needed to set up and evaluate a RKHS model
    type(kernel)                           :: pes2              ! The kernel type is needed to set up and evaluate a RKHS model
    !type(kernel)                           :: pes3              ! The kernel type is needed to set up and evaluate a RKHS model

    integer                                :: iO2
    integer                                :: iNO
    integer                                :: jNO

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
    !real(rkp), allocatable, dimension(:,:) :: darray3

    real(rkp)                              :: dr1
    real(rkp)                              :: dr2
    real(rkp)                              :: dr3
  
  contains
    procedure          ::  Initialize     =>    Initialize_NO2_Basel_PES
    procedure          ::  Output         =>    Output_NO2_Basel_PES
    procedure          ::  Compute        =>    Compute_NO2_Basel_PES_1d
    procedure          ::  Potential      =>    NO2_Basel_Potential_From_R
    procedure          ::  TriatPotential =>    NO2_Basel_Potential_From_R_OnlyTriat
    procedure ,private ::  ReadData_NO2_Basel_PES
    procedure ,private ::  NO_PES
    procedure ,private ::  pes3d
    procedure ,private ::  calcener
  End Type

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
  

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NASA PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_NO2_Basel_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                         ,only:  Input_Type
  use Atom_Class                          ,only:  Atom_Type
  use DiatomicPotential_Factory_Class    ,only:  DiatomicPotential_Factory_Type

  class(NO2_Basel_PES_Type)                 ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms  
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    ii, iP
  character(*)                    ,parameter                ::    Name_PES = 'NO2_Basel'
  integer         ,dimension(3,2)                           ::    iA
  integer         ,dimension(6)                             ::    MatrixTemp = (/ 1, 2, 1, 3, 2, 3 /)
  type(DiatomicPotential_Factory_Type)                      ::    DiatPotFactory
  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_NO2_Basel_PES" )
  !i_Debug_Loc   =     Logger%On()
    
  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair
  
  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]


  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( 'Atom 1 Name = ', Atoms(1)%Name ) 
  if (i_Debug_Loc) call Logger%Write( 'Atom 2 Name = ', Atoms(2)%Name ) 
  if (i_Debug_Loc) call Logger%Write( 'Atom 3 Name = ', Atoms(3)%Name ) 

  ii = 1
  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )

    if ( ( trim(adjustl(Atoms(MatrixTemp((iP-1)*2+1))%Name) ) .eq. "O") .and. ( trim(adjustl(Atoms(MatrixTemp((iP-1)*2+2))%Name)) .eq. "O" ) ) then
      !allocate( O2_Basel_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      This%iO2 = iP
      if (i_Debug_Loc) call Logger%Write( 'This%iO2 pair =',iP )
    else
     ! allocate( NO_Basel_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
      if (ii == 1) then
        This%jNO = iP
        if (i_Debug_Loc) call Logger%Write( 'This%jNO pair =',iP )
      else
        This%iNO = iP
        if (i_Debug_Loc) call Logger%Write( 'This%iNO pair =',iP )
      end if
      ii=ii+1
    end if
  end do
 ! ==============================================================================================================


  ! ==============================================================================================================
  !   READING VARIABLES
  ! ==============================================================================================================
  call This%ReadData_NO2_Basel_PES( Input, iPES, i_Debug=i_Debug_Loc )
  ! ==============================================================================================================


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadData_NO2_Basel_PES( This, Input, iPES, i_Debug )

  use Input_Class                         ,only:  Input_Type

  class(NO2_Basel_PES_Type)               ,intent(inout)  ::    This
  type(Input_Type)                        ,intent(in)     ::    Input
  integer                                 ,intent(in)     ::    iPES
  logical                       ,optional ,intent(in)     ::    i_Debug

  character(150)                                          ::    NO2_Basel_Fldr
  integer                                                 ::    ii
  integer                                                 ::    Unit
  integer                                                 ::    Status
  logical                                                 ::    i_Debug_Loc

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadData_NO2_Basel_PES" )
  !i_Debug_Loc   =     Logger%On()


  NO2_Basel_Fldr = trim(adjustl(Input%DtbPath))  // '/Systems/NO2/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '/'


  if (.not. This%ker1) then
    inquire(file=trim(NO2_Basel_Fldr)//"pes1.kernel", exist=This%ker1)   ! file_exists will be true if the file exists and false otherwise
  end if

  if (.not. This%ker2) then
    inquire(file=trim(NO2_Basel_Fldr)//"pes2.kernel", exist=This%ker2)   ! file_exists will be true if the file exists and false otherwise
  end if

  !if (.not. ker3) then
  !  inquire(file=trim(NO2_Basel_Fldr)//"pes23.kernel", exist=ker3)   ! file_exists will be true if the file exists and false otherwise
  !end if


  if (i_Debug_Loc) call Logger%Write( "Reading dr1, dr2, dr3 from File: ", trim(NO2_Basel_Fldr)//"Params.csv" )
  open(NewUnit=Unit, file=trim(NO2_Basel_Fldr)//"Params.csv", status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // trim(NO2_Basel_Fldr)//"Params.csv" )
    read(Unit,*) This%dr1, This%dr2, This%dr3
  close(Unit)
  if (i_Debug_Loc) call Logger%Write( "dr1 = ", This%dr1 )
  if (i_Debug_Loc) call Logger%Write( "dr2 = ", This%dr2 )
  if (i_Debug_Loc) call Logger%Write( "dr3 = ", This%dr3 )



  if (.not. This%stored ) then
  
    if (i_Debug_Loc) call Logger%Write( "Reading from File: ", trim(NO2_Basel_Fldr)//"asymp.dat" )
    open(NewUnit=Unit, file=trim(NO2_Basel_Fldr)//"asymp.dat", status='OLD', iostat=Status )
    if (Status/=0) call Error( "Error opening file: " // trim(NO2_Basel_Fldr)//"asymp.dat" )

      read(Unit,*) This%na1
      allocate( This%asy_array1(This%na1,2) )
      do ii = 1,This%na1
        read(Unit,*) This%asy_array1(ii,1), This%asy_array1(ii,2)
      end do

      read(Unit,*) This%na2
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
      allocate( This%darray1(This%nda1,2))
      do ii = 1,This%nda1
        read(Unit,*) This%darray1(ii,1), This%darray1(ii,2)
      end do

      read(Unit,*) This%nda2
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

      call This%pes1%load_from_file(trim(NO2_Basel_Fldr)//"pes1.kernel")
      call This%pes2%load_from_file(trim(NO2_Basel_Fldr)//"pes2.kernel")
    !  call This%pes3%load_from_file(trim(NO2_Basel_Fldr)//"pes23.kernel")
      This%kread = .true.
    else
      if (i_Debug_Loc) call Logger%Write( "Computing pes*.kernel Files" )

      call This%pes1%read_grid(trim(NO2_Basel_Fldr)//"pes1.csv")
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
      call This%pes1%save_to_file(trim(NO2_Basel_Fldr)//"pes1.kernel")
    

      call This%pes2%read_grid(trim(NO2_Basel_Fldr)//"pes2.csv")
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
      call This%pes2%save_to_file(trim(NO2_Basel_Fldr)//"pes2.kernel")


    !  call This%pes3%read_grid(trim(NO2_Basel_Fldr)//"pes3.csv")
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
    !  call This%pes3%save_to_file(trim(NO2_Basel_Fldr)//"pes3.kernel")
      

      This%kread = .true.
    end if
  end if


  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_NO2_Basel_PES( This, Unit )

  class(NO2_Basel_PES_Type)               ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function NO2_Basel_Potential_From_R( This, R, Q ) result( V )

  class(NO2_Basel_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  integer                                                    ::    iP
  real(rkp)                                                  ::    VAll
  real(rkp) ,dimension(3)                                    ::    dVAll
  real(rkp)                                                  ::    t1, t2
  
  ! Computing the diatomic potential energies associated to the 3 internuclear distances
  call This%NO_PES(R(This%iO2), R(This%iNO), R(This%jNO), VAll, dVAll)
  
  V  = VAll

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function NO2_Basel_Potential_From_R_OnlyTriat( This, R, Q ) result( V )

  class(NO2_Basel_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  integer                                                    ::    iP
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp) ,dimension(3)                                    ::    dVDiat
  real(rkp)                                                  ::    VAll
  real(rkp) ,dimension(3)                                    ::    dVAll
  real(rkp)                                                  ::    t1, t2

  ! Computing the diatomic potential energies associated to the 3 internuclear distances
  do iP=1,3
    call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do     

  call This%NO_PES(R(This%iO2), R(This%iNO), R(This%jNO), VAll, dVAll)

  V = VAll - sum(VDiat)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_NO2_Basel_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(NO2_Basel_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  integer                                                    ::    iP
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp) ,dimension(3)                                    ::    dVDiat
  real(rkp)                                                  ::    VAll
  real(rkp) ,dimension(3)                                    ::    dVAll
  real(rkp)                                                  ::    t1, t2

  !call cpu_time ( t1 )

  dVdQ  = Zero


  call This%NO_PES(R(This%iO2), R(This%iNO), R(This%jNO), VAll, dVAll)


  ! ! Computing the diatomic potential energies associated to the 3 internuclear distances
  ! call This%Pairs(This%iO2)%Vd%Compute_Vd_dVd( R(This%iO2), VDiat(This%iO2), dVDiat(This%iO2) )
  ! call This%Pairs(This%iNO)%Vd%Compute_Vd_dVd( R(This%iNO), VDiat(This%iNO), dVDiat(This%iNO) )
  ! call This%Pairs(This%jNO)%Vd%Compute_Vd_dVd( R(This%jNO), VDiat(This%jNO), dVDiat(This%jNO) )

  V       = VAll           ! - sum(VDiat)
  dVdR(1) = dVAll(This%iO2)! - dVDiat(1)
  dVdR(2) = dVAll(This%iNO)! - dVDiat(2)
  dVdR(3) = dVAll(This%jNO)! - dVDiat(3)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!



! **************************************************************************************************************
! **************************************************************************************************************
!                                   PRIVATE PROCEDURES for BASEL PESs
! **************************************************************************************************************
! **************************************************************************************************************

!________________________________________________________________________________________________________________________________!
Subroutine NO_PES(This, tmpr1, tmpr2, tmpr3, totener, dvdr)
  
  class(NO2_Basel_PES_Type) ,intent(in)  :: This
  real(rkp)                 ,intent(in)  :: tmpr1
  real(rkp)                 ,intent(in)  :: tmpr2
  real(rkp)                 ,intent(in)  :: tmpr3
  real(rkp)                 ,intent(out) :: totener
  real(rkp), dimension(3)   ,intent(out) :: dvdr

  real(rkp), dimension(3)                :: xp
  real(rkp), dimension(3)                :: tmpdvdr
  real(rkp), dimension(3)                :: r
  real(rkp)                              :: ener
  integer                                :: ii

  real(rkp), parameter                   :: dx = 0.0005e0_rkp

  !      3
  !      N
  !     / \
  !  r3/   \r2
  !   /     \
  !  /       \
  ! O---------O
  ! 1   r1    2

  !O+ON
  !OO ON NO

  r(1) = tmpr1  !OO
  r(2) = tmpr2  !ON
  r(3) = tmpr3  !NO
  dvdr = Zero

  call This%pes3d(r(1),r(2),r(3),totener,dvdr)

  if ( any(isNaN(dvdr)) ) go to 100
  if ((dvdr(1)-1) .eq. dvdr(1) .or. (dvdr(3)-1) .eq. dvdr(3) .or. (dvdr(2)-1) .eq. dvdr(2) ) go to 100
  return
  !100 print*,"NaN or Infinity"
  
  100 if ((r(1)+r(2)-r(3))<1d-10) then
  !==========================================================
  !   Forward difference to calculate first derivative      =
  !==========================================================
      xp=r
      xp(1)=r(1)+dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(1)=(ener-totener)/dx
      xp=r
      xp(2)=r(2)+dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(2)=(ener-totener)/dx
      xp=r
      xp(3)=r(3)-dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(3)=(totener-ener)/dx
      return
  end if
  
  if ((r(2)+r(3)-r(1))<1d-10) then
  !==========================================================
  !   Forward difference to calculate first derivative      =
  !==========================================================
      xp=r
      xp(1)=r(1)-dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(1)=(totener-ener)/dx
      xp=r
      xp(2)=r(2)+dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(2)=(ener-totener)/dx
      xp=r
      xp(3)=r(3)+dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(3)=(ener-totener)/dx
      return
  end if
  
  if ((r(1)+r(3)-r(2))<1d-10) then
  !==========================================================
  !   Forward difference to calculate first derivative      =
  !==========================================================
      xp=r
      xp(1)=r(1)+dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(1)=(ener-totener)/dx
      xp=r
      xp(2)=r(2)-dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(2)=(totener-ener)/dx
      xp=r
      xp(3)=r(3)+dx
      call This%pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
      dvdr(3)=(ener-totener)/dx
      return
  end if

End Subroutine NO_PES
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine pes3d(This, r12, r23, r31, totener, dvdr)
  
  class(NO2_Basel_PES_Type) ,intent(in)  :: This  
  real(rkp)                 ,intent(in)  :: r12
  real(rkp)                 ,intent(in)  :: r23
  real(rkp)                 ,intent(in)  :: r31
  real(rkp)                 ,intent(out) :: totener
  real(rkp) ,dimension(3)   ,intent(out) :: dvdr

  real(rkp)                              :: capr
  real(rkp)                              :: theta
  real(rkp) ,dimension(3,3)              :: derv
  real(rkp) ,dimension(3,3)              :: deri
  real(rkp) ,dimension(3,3)              :: dw

  real(rkp) ,dimension(3)                :: derj
  real(rkp) ,dimension(3)                :: w
  real(rkp) ,dimension(3)                :: ener
  real(rkp) ,dimension(3)                :: r
  real(rkp) ,dimension(3)                :: m

  real(rkp) ,parameter                   :: m1 = 15.99491462e0_rkp
  real(rkp) ,parameter                   :: m2 = m1
  real(rkp) ,parameter                   :: m3 = 14.003074e0_rkp

  real(rkp) ,dimension(3)                :: aa
  real(rkp)                              :: r1, r2, r3

  !      3
  !      N
  !     / \
  !  r3/   \r2
  !   /     \
  !  /       \
  ! O---------O
  ! 1   r1    2

  !=====================================================
  !
  !PES1                     PES2            PES3
  !        3             3                 3
  !        N             N                 N
  !  r2   /               \               /
  !      /  r3             \r1        r1 /
  !theta/             theta/\           /\ theta r2
  !O___/_____O         r2 /  \         /  \
  !1   r1    2           /    O      1O    \
  !                     / r3  2             \
  !                    O1                r3  O2

  r1 = r12
  r2 = r31
  r3 = r23

  ! aa = [0.15000000000000000000d+01,    0.18919597989949759054d+01,    0.10863257659814278266d+01]
  ! r1 = aa(1)
  ! r2 = aa(2)
  ! r3 = aa(3)

  r(1) = r1
  r(2) = r2
  r(3) = r3
  m(1) = m2
  m(2) = m1
  m(3) = m3
  call crdtrf(r,m,capr,theta,derv)

  call This%calcener(capr,r1,theta, ener(1), derj, 1)
 
  deri(1,1) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
  deri(1,3) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
  deri(1,2) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

  r(1) = r3
  r(2) = r2
  r(3) = r1
  m(1) = m2
  m(2) = m3
  m(3) = m1
  call crdtrf(r,m,capr,theta,derv)
  call This%calcener(capr,r3,theta, ener(2), derj, 2)

  deri(2,2) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
  deri(2,3) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
  deri(2,1) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

  r(1) = r2
  r(2) = r3
  r(3) = r1
  m(1) = m1
  m(2) = m3
  m(3) = m2
  call crdtrf(r,m,capr,theta,derv)
  call This%calcener(capr,r2,theta, ener(3), derj, 3)

  deri(3,3) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
  deri(3,2) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
  deri(3,1) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

  call weights(r1, r3, r2, This%dr1, This%dr2, This%dr3, w, dw)
  totener = sum(ener*w)
 
  dvdr(1) = sum(deri(:,1)*w)+sum(dw(:,1)*ener)
  dvdr(2) = sum(deri(:,2)*w)+sum(dw(:,2)*ener)
  dvdr(3) = sum(deri(:,3)*w)+sum(dw(:,3)*ener)

  ! write(*,*) 'deri(1,:) = ', deri(1,:)
  ! write(*,*) 'deri(2,:) = ', deri(2,:)
  ! write(*,*) 'deri(3,:) = ', deri(3,:)
  ! write(*,*) 'totener   = ', totener
  ! write(*,*) 'dvdr      = ', dvdr
  ! write(*,*) 
  ! pause

End Subroutine pes3d
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine calcener(This, capr, smlr, theta, ener, der, sno)

  class(NO2_Basel_PES_Type) ,intent(in)    :: This    
  real(rkp)                 ,intent(in)    :: capr
  real(rkp)                 ,intent(in)    :: smlr
  real(rkp)                 ,intent(in)    :: theta
  real(rkp)                 ,intent(out)   :: ener
  real(rkp) ,dimension(3)   ,intent(out)   :: der
  integer                   ,intent(in)    :: sno

  real(rkp)                                 :: asener, anener, asder, z1, z2
  real(rkp) ,dimension(3)                   :: ander, x
  integer                                   :: kk, ii

  x(1)=(One-cos(theta))/Two
  x(2)=capr
  x(3)=smlr

  if (sno==1) then
    asener = Zero
    asder  = Zero
    do kk = 1,This%na1
      asener = asener + drker26(smlr,  This%asy_array1(kk,1)) * This%asy_array1(kk,2)
      asder  = asder  + ddrker26(smlr, This%asy_array1(kk,1)) * This%asy_array1(kk,2)
    end do

    anener = Zero
    ander  = Zero
    call This%pes1%evaluate_fast(x,anener,ander)

    ener   = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/Two
  end if

  if (sno==2) then
    asener = Zero
    asder  = Zero
    do kk = 1,This%na2
      asener = asener + drker26(smlr,  This%asy_array2(kk,1)) * This%asy_array2(kk,2)
      asder  = asder  + ddrker26(smlr, This%asy_array2(kk,1)) * This%asy_array2(kk,2)
    end do

    anener = Zero
    ander  = Zero
    call This%pes2%evaluate_fast(x,anener,ander)

    ener   = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/Two
  end if

  if (sno==3) then
    asener = Zero
    asder  = Zero
    do kk = 1,This%na2
      asener = asener + drker26(smlr,  This%asy_array2(kk,1)) * This%asy_array2(kk,2)
      asder  = asder  + ddrker26(smlr, This%asy_array2(kk,1)) * This%asy_array2(kk,2)
    end do

    anener = Zero
    ander  = Zero
    call This%pes2%evaluate_fast(x,anener,ander)

    ener   = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/Two
  end if

  return

End Subroutine calcener
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine crdtrf(r, m, capr, theta, derv)

  implicit none

  real(rkp) ,dimension(:)   ,intent(in)  :: r(3)
  real(rkp) ,dimension(:)   ,intent(in)  :: m(3)
  real(rkp)                 ,intent(out) :: capr
  real(rkp)                 ,intent(out) :: theta
  real(rkp) ,dimension(:,:) ,intent(out) :: derv(3,3)

  real(rkp)                              :: cmu1, cmu2, rcm1, rcm2

  real(rkp) ,parameter                   :: pc = sqrt(epsilon(One))

!     #       3       #
!     #      /|\      #
!     #     / | \     #
!     #  r3/ R|  \r2  #
!     #   /   |th \   #
!     #  /____|____\  #
!     # 1    r1    2  #

  cmu1 = m(2) / (m(1)+m(2))
  cmu2 = m(1) / (m(1)+m(2))

  rcm1 = r(1) * cmu1
  rcm2 = r(1) * cmu2

  capr = sqrt (r(2)**2/r(1)*rcm1 + r(3)**2/r(1)*rcm2 - rcm1*rcm2 )
  if (abs(capr) < pc) capr = pc

  theta = (rcm2**2+capr**2-r(2)**2)/Two/rcm2/capr
  theta = min(One, max(-One,theta))
  theta = acos(theta)

  derv(1,1) = -cmu1*cmu2*r(1)/capr

  derv(1,2) = r(2)*cmu1/capr

  derv(1,3) = r(3)*cmu2/capr

  derv(2,1) = One 
  derv(2,2) = Zero
  derv(2,3) = Zero

  derv(3,1) = (derv(1,1)/capr*cos(theta)+cos(theta)/r(1)-(capr*derv(1,1)+rcm2*cmu2)/rcm2/capr) /sqrt(One-cos(theta)**2)

  derv(3,2) = (r(2)/rcm2/capr-derv(1,2)/rcm2+cos(theta)/capr*derv(1,2))/sqrt(One-cos(theta)**2)

  derv(3,3) = (cos(theta)/capr*derv(1,3)-derv(1,3)/rcm2)/sqrt(One-cos(theta)**2)

  return

End Subroutine crdtrf
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Subroutine weights(d1, d2, d3, dr1, dr2, dr3, w, dw)
  
  implicit none
  
  real(rkp)                 ,intent(in)  :: d1
  real(rkp)                 ,intent(in)  :: d2
  real(rkp)                 ,intent(in)  :: d3
  real(rkp)                 ,intent(in)  :: dr1
  real(rkp)                 ,intent(in)  :: dr2
  real(rkp)                 ,intent(in)  :: dr3
  real(rkp) ,dimension(:)   ,intent(out) :: w(3)
  real(rkp) ,dimension(:,:) ,intent(out) :: dw(3,3)
  
  real(rkp)                              :: wsum, pw1, pw2, pw3, r1, r2, r3
  
  integer   ,parameter                   :: power = 2
  real(rkp) ,parameter                   :: wtol  = epsilon(One)

  r1=d1
  r2=d2
  r3=d3

  wsum = Zero
  do while(wsum < wtol)
    pw1 = exp(-(r1/dr1)**power)
    pw2 = exp(-(r2/dr2)**power)
    pw3 = exp(-(r3/dr3)**power)
    wsum = pw1+pw2+pw3
    if(wsum < wtol) then
  !    r1=d1/(d1+d2+d3)*3.0d0
  !    r2=d2/(d1+d2+d3)*3.0d0
  !    r3=d3/(d1+d2+d3)*3.0d0
      r1=r1-0.1d0
      r2=r2-0.1d0
      r3=r3-0.1d0
    end if
  end do

    w(1) = pw1 / wsum
    w(2) = pw2 / wsum
    w(3) = pw3 / wsum

    dw(1,1) =-power*((r1/dr1)**(power-1))*(pw1*pw2+pw1*pw3)/dr1/wsum**2
    dw(1,2) = power*((r2/dr2)**(power-1))*pw1*pw2/dr2/wsum**2
    dw(1,3) = power*((r3/dr3)**(power-1))*pw1*pw3/dr3/wsum**2

    dw(2,1) = power*((r1/dr1)**(power-1))*pw2*pw1/dr1/wsum**2
    dw(2,2) =-power*((r2/dr2)**(power-1))*(pw2*pw1+pw2*pw3)/dr2/wsum**2
    dw(2,3) = power*((r3/dr3)**(power-1))*pw2*pw3/dr3/wsum**2


    dw(3,1) = power*((r1/dr1)**(power-1))*pw3*pw1/dr1/wsum**2
    dw(3,2) = power*((r2/dr2)**(power-1))*pw3*pw2/dr2/wsum**2
    dw(3,3) =-power*((r3/dr3)**(power-1))*(pw3*pw1+pw3*pw2)/dr3/wsum**2

  return

End Subroutine weights
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function drker24(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: drker24, xl, xs

  xl = x
  xs = xi
  if (x .lt. xi) then
    xl = xi
    xs = x
  end if

  drker24 = dk24f1/xl**5 - dk24f2*xs/xl**6

End Function drker24
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function ddrker24(x,xi)
 
  implicit none
 
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: ddrker24, xl, xs

  if (x .lt. xi) then
    ddrker24 = -dk24f2/xi**6
  else
    ddrker24 = -Five*dk24f1/x**6 + Six*dk24f2*xi/x**7
  end if

End Function ddrker24
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function drker25(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: drker25, xl, xs

  xl = x
  xs = xi
  if (x .lt. xi) then
    xl = xi
    xs = x
  end if

  drker25 = dk25f1/xl**6 - dk25f2*xs/xl**7

End Function drker25
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function ddrker25(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: ddrker25, xl, xs

  if (x .lt. xi) then
    ddrker25 = -dk25f2/xi**7
  else
    ddrker25 = -Six*dk25f1/x**7 + Seven*dk25f2*xi/x**8
  end if

End Function ddrker25
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


!________________________________________________________________________________________________________________________________!
Pure Function atker23(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: atker23, xl, xs

  xl = x
  xs = xi
  if (x .lt. xi) then
    xl = xi
    xs = x
  end if

  atker23 = One + xs*xl + Two*xs**2*xl  - akf1*xs**3

End Function atker23
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function datker23(x,xi)
  
  implicit none
  
  real(rkp) ,intent(in) :: x, xi
  real(rkp)             :: datker23, xl, xs

  if (x .lt. xi) then
    datker23 = xi + Four*x*xi  - Three*akf1*x**2
  else
    datker23 = xi + Two*xi**2
  end if

End Function datker23
!________________________________________________________________________________________________________________________________!
  

End Module NO2_Basel_PES_Class