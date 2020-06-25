! -*-F90-*-
!===============================================================================================================
! 
! Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
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

Module CNO_Basel_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Kcm_To_Hartree, B_To_Ang
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use Input_Class           ,only:  Input_Type
  use RKHS
  
  implicit none

  private
  public    :: CNO_Basel_PES_Type 
  
  Type    ,extends(PES_Type)                :: CNO_Basel_PES_Type

    character(10)                           :: model_output
    real(rkp)                               :: Rmin
    real(rkp)                               :: Vmin
    real(rkp)                               :: dr1, dr2, dr3
    real(rkp), allocatable, dimension(:,:)  :: asy_array1, asy_array2, asy_array3, &
                                               darray1, darray2, darray3  
    integer                                 :: na1, na2, na3, nda1, nda2, nda3
    integer                                 :: iCN
    integer                                 :: iCO
    integer                                 :: iNO
    
    type(kernel)                            :: pes1         ! The kernel type is needed to set up and evaluate a RKHS model
    type(kernel)                            :: pes2         ! The kernel type is needed to set up and evaluate a RKHS model
    type(kernel)                            :: pes3         ! The kernel type is needed to set up and evaluate a RKHS model
    logical                                 :: ComputeDiatPotFlg = .False.
    logical                                 :: stored   = .false.
    logical                                 :: kread    = .false.
    logical                                 :: ker1     = .false.
    logical                                 :: ker2     = .false.
    logical                                 :: ker3     = .false.

  contains
    procedure                               ::  Initialize        =>    Initialize_CNO_Basel_PES
    procedure                               ::  Output            =>    Output_CNO_Basel_PES
    procedure                               ::  Compute           =>    Compute_CNO_Basel_PES_1d
    procedure                               ::  Potential         =>    CNO_Basel_Potential_From_R
    procedure                               ::  TriatPotential    =>    CNO_Basel_Potential_From_R_OnlyTriat
    procedure                               ::  DiatPotential     =>    CNO_Basel_Potential_From_R_OnlyDiat
    procedure, private                      ::  ReadData_CNO_Basel_PES
    procedure, private                      ::  calcener
    procedure, private                      ::  OCN_PES
    procedure, private                      ::  pes3d
  End Type
  
  type(Input_Type)                          ::    Input
  
  logical                     ,parameter    :: i_Debug_Global = .False.
  real(rkp)                   ,parameter    :: dk26f1  = One/14.0_rkp
  real(rkp)                   ,parameter    :: dk26f2  = One/18.0_rkp
  real(rkp)                   ,parameter    :: wtol    = epsilon(1.0d0)
  real(rkp)                   ,parameter    :: pc      = sqrt(epsilon(1.0d0))
  real(rkp)                   ,parameter    :: m1 = 12.0d0, m2 = 14.003074d0, m3 = 15.99491462d0
  real(rkp)                   ,parameter    :: dr1_2ap  = 0.71d0, dr2_2ap  = 0.69d0, dr3_2ap  = 0.68d0
  real(rkp)                   ,parameter    :: dr1_2app = 0.48d0, dr2_2app = 0.40d0, dr3_2app = 0.28d0
  real(rkp)                   ,parameter    :: dr1_4app = 1.06d0, dr2_4app = 0.80d0, dr3_4app = 1.0d0

  contains

! **************************************************************************************************************
! **************************************************************************************************************
!                                       DEFERRED PROCEDURES for CNO Basel PESs
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_CNO_Basel_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(CNO_Basel_PES_Type)                 ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms  
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  type(DiatomicPotential_Factory_Type)                       ::    DiatPotFactory   
  integer                                                   ::    iP
  character(*)                    ,parameter                ::    Name_PES = 'CNO_Basel'
  integer         ,dimension(3,2)                           ::    iA
  integer         ,dimension(6)                             ::    MatrixTemp = (/ 1, 2, 1, 3, 2, 3 /) 

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_CNO_Basel_PES" )
  !i_Debug_Loc   =     Logger%On()


  This%Name         =   Name_PES
  This%Model        =   Input%PES_Model(iPES)
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  ! allocate( This%mMiMn(3) )
  ! This%mMiMn(1:2) = - Atoms(1:2)%Mass / Atoms(3)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )
  
  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( 'Atom 1 Name = ', Atoms(1)%Name ) 
  if (i_Debug_Loc) call Logger%Write( 'Atom 2 Name = ', Atoms(2)%Name ) 
  if (i_Debug_Loc) call Logger%Write( 'Atom 3 Name = ', Atoms(3)%Name ) 

  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )

    if ( ( trim(adjustl(Atoms(MatrixTemp((iP-1)*2+1))%Name) ) .eq. "C") .and. ( trim(adjustl(Atoms(MatrixTemp((iP-1)*2+2))%Name)) .eq. "N" ) ) then
      This%iCN = iP
      if (i_Debug_Loc) call Logger%Write(  'This%iCN pair = ',iP )
    elseif ( ( trim(adjustl(Atoms(MatrixTemp((iP-1)*2+1))%Name) ) .eq. "C") .and. ( trim(adjustl(Atoms(MatrixTemp((iP-1)*2+2))%Name)) .eq. "O" ) ) then 
       This%iCO = iP
       if (i_Debug_Loc) call Logger%Write( 'This%iCO pair = ',iP )
    else
       This%iNO = iP
       if (i_Debug_Loc) call Logger%Write( 'This%iNO pair = ',iP )
    end if

  end do
  ! ==============================================================================================================

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiaPotFactory%Construct" )
  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )
  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
  ! ==============================================================================================================

  if    (adjustl(trim(This%Model))  == 'Basel_2A1') then
    This%dr1   =  dr1_2ap
    This%dr2   =  dr2_2ap
    This%dr3   =  dr3_2ap
 elseif (adjustl(trim(This%Model))  == 'Basel_2A2') then
    This%dr1   =  dr1_2app
    This%dr2   =  dr2_2app
    This%dr3   =  dr3_2app
  elseif (adjustl(trim(This%Model)) == 'Basel_4A2') then
    This%dr1   =  dr1_4app
    This%dr2   =  dr2_4app
    This%dr3   =  dr3_4app
  else
    write(*,*) 'PES_Model different from Basel_2A1, Basel_2A2 or Basel_4A2'
    stop
  end if

  ! ==============================================================================================================
  ! READING PES DATA
  ! ==============================================================================================================
  call This%ReadData_CNO_Basel_PES( Input, iPES, i_Debug=i_Debug_Loc )
  ! ==============================================================================================================
  
  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine ReadData_CNO_Basel_PES( This, Input, iPES, i_Debug )

  use Input_Class                         ,only:  Input_Type

  class(CNO_Basel_PES_Type)               ,intent(inout)  ::    This
  type(Input_Type)                        ,intent(in)     ::    Input
  integer                                 ,intent(in)     ::    iPES
  logical                       ,optional ,intent(in)     ::    i_Debug

  integer                                                 ::    ii
  integer                                                 ::    Unit
  integer                                                 ::    Status
  logical                                                 ::    i_Debug_Loc
  character(150)                                          ::    Basel_CNO_folder

  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "ReadData_CNO_Basel_PES" )
  !i_Debug_Loc   =     Logger%On()

  Basel_CNO_folder = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '/'
  if (i_Debug_Loc) call Logger%Write( "Folder for Reading Basel CNO Parameters: ", Basel_CNO_folder )

  if (i_Debug_Loc) call Logger%Write( "Calling External Module for Initializing PES and Reading Basel CNO Parameters: ", Basel_CNO_folder )
  
  if (.not. This%ker1) then
     inquire(file=Basel_CNO_folder//"/pes1.kernel", exist=This%ker1)   ! file_exists will be true if the file exists and false otherwise
  end if

  if (.not. This%ker2) then
     inquire(file=Basel_CNO_folder//"/pes2.kernel", exist=This%ker2)   ! file_exists will be true if the file exists and false otherwise
  end if

  if (.not. This%ker3) then
     inquire(file=Basel_CNO_folder//"/pes3.kernel", exist=This%ker3)   ! file_exists will be true if the file exists and false otherwise
  end if

  if (.not. This%stored ) then
     if (i_Debug_Loc) call Logger%Write( "Reading from File: ", trim(Basel_CNO_Folder)//"asymp.dat" )
     open(NewUnit=Unit,file=Basel_CNO_folder//"/asymp.dat", status = 'OLD', iostat=Status)
     if (Status/=0) call Error( "Error opening file: " // trim(Basel_CNO_Folder)//"asymp.dat" )
  
     read(Unit,*)This%na1
     allocate(This%asy_array1(This%na1,2))
     do ii = 1, This%na1
        read(Unit,*)This%asy_array1(ii,1), This%asy_array1(ii,2)
     end do
     
     read(Unit,*)This%na2
     allocate(This%asy_array2(This%na2,2))
     do ii = 1, This%na2
        read(Unit,*)This%asy_array2(ii,1), This%asy_array2(ii,2)
     end do
  
     read(Unit,*)This%na3
     allocate(This%asy_array3(This%na3,2))
     do ii = 1, This%na3
        read(Unit,*)This%asy_array3(ii,1), This%asy_array3(ii,2)
     end do
     
     read(Unit,*)This%nda1
     allocate(This%darray1(This%nda1,2))
     do ii = 1, This%nda1
        read(Unit,*)This%darray1(ii,1), This%darray1(ii,2)
     end do
     
     read(Unit,*)This%nda2
     allocate(This%darray2(This%nda2,2))
     do ii = 1, This%nda2
        read(Unit,*)This%darray2(ii,1), This%darray2(ii,2)
     end do
     
     read(Unit,*)This%nda3
     allocate(This%darray3(This%nda3,2))
     do ii = 1, This%nda3
        read(Unit,*)This%darray3(ii,1), This%darray3(ii,2)
     end do
     
     This%stored = .true.

     close(Unit)
     
  end if
  
  if (.not. This%kread) then
     if (This%ker1 .and. This%ker2 .and. This%ker3 ) then
        if (i_Debug_Loc) call Logger%Write( "Reading pes*.kernel Files" )  

        call This%pes1%load_from_file(Basel_CNO_folder//"/pes1.kernel")
        call This%pes2%load_from_file(Basel_CNO_folder//"/pes2.kernel")
        call This%pes3%load_from_file(Basel_CNO_folder//"/pes3.kernel")
        This%kread = .true.
     else
        if (i_Debug_Loc) call Logger%Write( "Computing pes*.kernel Files" )
        
        call This%pes1%read_grid(Basel_CNO_folder//"/pes1.csv")
  !    print*,"IAMHERE"
        call This%pes1%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
        call This%pes1%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
        call This%pes1%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL) 
        !    call pes11%calculate_coefficients_slow(lambda)
        call This%pes1%calculate_coefficients_fast()
        
        call This%pes1%calculate_sums()
        
        call This%pes1%save_to_file(Basel_CNO_folder//"/pes1.kernel")
        !
        call This%pes2%read_grid(Basel_CNO_folder//"/pes2.csv")
        !    print*,"IAMHERE"
        call This%pes2%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
        call This%pes2%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
        call This%pes2%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
        
        !    call pes12%calculate_coefficients_slow(lambda)
        call This%pes2%calculate_coefficients_fast()
        
        call This%pes2%calculate_sums()
        
        call This%pes2%save_to_file(Basel_CNO_folder//"/pes2.kernel")
        
        call This%pes3%read_grid(Basel_CNO_folder//"/pes3.csv")
        
        call This%pes3%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
        call This%pes3%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
        call This%pes3%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
        
        !    call pes13%calculate_coefficients_slow(lambda)
        call This%pes3%calculate_coefficients_fast()
        
        call This%pes3%calculate_sums()
        
        call This%pes3%save_to_file(Basel_CNO_folder//"/pes3.kernel")
        
        This%kread = .true.
     end if
  end if
  
  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_CNO_Basel_PES( This, Unit )

  class(CNO_Basel_PES_Type)               ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function CNO_Basel_Potential_From_R( This, R, Q ) result( V )
   
  class(CNO_Basel_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp)                                                  ::    totener
  real(rkp) ,dimension(3)                                    ::    dvdr, RTemp

  RTemp(1) = R(This%iCN)
  RTemp(2) = R(This%iNO)
  RTemp(3) = R(This%iCO)

  call This%OCN_PES(RTemp,totener,dvdr)
  V = totener

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function CNO_Basel_Potential_From_R_OnlyDiat( This, R, Q ) result( V )

  class(CNO_Basel_PES_Type)                     ,intent(in)  ::    This   !< PES class
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree]. 
  integer                                                    ::    iP
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp) ,dimension(3)                                    ::    dVDiat

  do iP=1,3
     call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do
  
  !-------------------------!
  ! Form and sum contributions of 2-body interactions 
  V = sum(VDiat)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function CNO_Basel_Potential_From_R_OnlyTriat( This, R, Q ) result( V )
   
  class(CNO_Basel_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    RTemp

  real(rkp)                                                  ::    totener
  real(rkp) ,dimension(3)                                    ::    dvdr, VDiat, dVDiat
  integer                                                    ::    iP

  RTemp(1) = R(This%iCN)
  RTemp(2) = R(This%iNO)
  RTemp(3) = R(This%iCO)

  do iP=1,3
     call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
  end do

  call This%OCN_PES(RTemp,totener,dvdr)

  V = totener - sum(VDiat)

End Function CNO_Basel_Potential_From_R_OnlyTriat
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_CNO_Basel_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(CNO_Basel_PES_Type)                     ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat, dVdrTemp
  real(rkp) ,dimension(3)                                    ::    X, Y, Z
  real(rkp) ,dimension(3)                                    ::    RTemp

  real(rkp)                                                  ::    totener

  RTemp(1) = R(This%iCN)
  RTemp(2) = R(This%iNO)
  RTemp(3) = R(This%iCO)
  
  dvdr = Zero

  call This%OCN_PES(RTemp,totener,dVdrTemp)

  V = totener
  dVdR(1) = dVdrTemp(This%iCN)! - dVDiat(1)
  dVdR(2) = dVdrTemp(This%iNO)! - dVDiat(1)
  dVdR(3) = dVdrTemp(This%iCO)! - dVDiat(1)

  dVdQ = Zero
  call This%TransToCart_3Atoms( R, Q, dVdR, dVdQ)

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
! **************************************************************************************************************
!                                         PRIVATE PROCEDURES for LEPS
! **************************************************************************************************************

!________________________________________________________________________________________________________________________________
pure function drker26(x,xi)
  
  implicit none
  
  real(rkp), intent(in) :: x, xi
  real(rkp)             :: drker26, xl, xs

  xl = x
  xs = xi
  if (x .lt. xi) then
    xl = xi
    xs = x
  end if

  drker26 = dk26f1/xl**7 - dk26f2*xs/xl**8

end function drker26
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!  
pure function ddrker26(x,xi)
  
  implicit none
  
  real(rkp), intent(in) :: x, xi
  real(rkp)             :: ddrker26, xl, xs

  if (x .lt. xi) then
    ddrker26 = -dk26f2/xi**8
  else
    ddrker26 = -7.0d0*dk26f1/x**8 + 8.0d0*dk26f2*xi/x**9
  end if

end function ddrker26
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine pes3d(This,r12,r23,r31,totener,dvdr)
  
  implicit none
  
  class(CNO_Basel_PES_Type) ,intent(in)  ::    This  
  real(rkp)                 ,intent(in)  ::    r12, r23, r31
  real(rkp), dimension(:)   ,intent(out) ::    dvdr(3)
  real(rkp)                 ,intent(out) ::    totener
  real(rkp)                              ::    capr, theta
  real(rkp), dimension(:,:)              ::    derv(3,3), deri(3,3), dw(3,3)
  real(rkp), dimension(:)                ::    derj(3), w(3), ener(3), r(3), m(3)

  !      3
  !      O
  !     / \
  !  r3/   \r2
  !   /     \
  !  /       \
  ! C---------N
  ! 1   r1    2

  !=====================================================
  !
  !PES1                     PES2            PES3
  !        3             3                 3
  !        O             O                 O
  !  r3   /               \               /
  !      /  r2             \r1        r1 /
  !     /th                /\           /\  r3
  !C___/_____N         r3 /th\         /  \
  !1   r1    2           /    N      1C    \
  !                     / r2  2        theta\
  !                   1C                r2  N2
  !

      r(1)=r12
      r(2)=r23
      r(3)=r31
      m(1) = m1
      m(2) = m2
      m(3) = m3
      call crdtrf(r,m,capr,theta,derv)
      call This%calcener(capr,r12,theta, ener(1), derj, 1)
          
      deri(1,1) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
      deri(1,2) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
      deri(1,3) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)
    
      r(1)=r23
      r(2)=r12
      r(3)=r31
      m(1) = m3
      m(2) = m2
      m(3) = m1
      call crdtrf(r,m,capr,theta,derv)
      call This%calcener(capr,r23,theta, ener(2), derj, 2)

      deri(2,2) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
      deri(2,1) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
      deri(2,3) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

      r(1)=r31
      r(2)=r12
      r(3)=r23
      m(1) = m3
      m(2) = m1
      m(3) = m2
      call crdtrf(r,m,capr,theta,derv)
      call This%calcener(capr,r31,theta, ener(3), derj, 3)
      ener(3) = ener(3) + 0.090024267768d0

      deri(3,3) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
      deri(3,1) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
      deri(3,2) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

      call weights(This,r12, r23, r31, w, dw)
      totener = sum(ener*w)
     
      dvdr(1) = sum(deri(:,1)*w)+sum(dw(:,1)*ener)
      dvdr(2) = sum(deri(:,2)*w)+sum(dw(:,2)*ener)
      dvdr(3) = sum(deri(:,3)*w)+sum(dw(:,3)*ener)

end subroutine pes3d
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine crdtrf(r,m,capr,theta, derv)
  
  implicit none
  
  real(rkp), dimension(:)    ,intent(in)  ::   r(3), m(3)
  real(rkp), dimension(:,:)  ,intent(out) ::   derv(3,3)
  real(rkp)                  ,intent(out) ::   capr, theta
  real(rkp)                               :: cmu1, cmu2, rcm1, rcm2

  cmu1 = m(2) / (m(1)+m(2))
  cmu2 = m(1) / (m(1)+m(2))

  rcm1 = r(1) * cmu1
  rcm2 = r(1) * cmu2

  capr = sqrt (r(2)**2/r(1)*rcm1 + r(3)**2/r(1)*rcm2 - rcm1*rcm2 )
  if (abs(capr) < pc) capr = pc

  theta = (rcm2**2+capr**2-r(2)**2)/2.0d0/rcm2/capr
  theta=min(1.0d0,max(-1.0d0,theta))
  theta = acos(theta)

  derv(1,1) = -cmu1*cmu2*r(1)/capr

  derv(1,2) = r(2)*cmu1/capr

  derv(1,3) = r(3)*cmu2/capr

  derv(2,1) = 1.0d0
  derv(2,2) = 0.0d0
  derv(2,3) = 0.0d0

  derv(3,1) = (derv(1,1)/capr*cos(theta)+cos(theta)/r(1)-(capr*derv(1,1)+rcm2*cmu2)/rcm2/capr)&
              /sqrt(1.0d0-cos(theta)**2)

  derv(3,2) = (r(2)/rcm2/capr-derv(1,2)/rcm2+cos(theta)/capr*derv(1,2))/sqrt(1.0d0-cos(theta)**2)

  derv(3,3) = (cos(theta)/capr*derv(1,3)-derv(1,3)/rcm2)/sqrt(1.0d0-cos(theta)**2)

  return

end subroutine crdtrf
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine weights(This,d1,d2,d3,w,dw)
  
  implicit none

  class(CNO_Basel_PES_Type) ,intent(in)  ::    This
  real(rkp)                 ,intent(in)  ::    d1, d2, d3
  real(rkp), dimension(:)   ,intent(out) ::    w(3)
  real(rkp), dimension(:,:) ,intent(out) ::    dw(3,3)
  real(rkp)                              ::    wsum, pw1, pw2, pw3, r1, r2, r3
  integer, parameter                     ::    power = 2


  r1=d1
  r2=d2
  r3=d3

  wsum = 0.0d0
  do while(wsum < wtol)
    pw1 = exp(-(r1/This%dr1)**power)
    pw2 = exp(-(r2/This%dr2)**power)
    pw3 = exp(-(r3/This%dr3)**power)
    wsum = pw1+pw2+pw3
    if(wsum < wtol) then
      r1=r1-0.1d0
      r2=r2-0.1d0
      r3=r3-0.1d0
    end if
  end do

    w(1) = pw1 / wsum
    w(2) = pw2 / wsum
    w(3) = pw3 / wsum

    dw(1,1) =-power*((r1/This%dr1)**(power-1))*(pw1*pw2+pw1*pw3)/This%dr1/wsum**2
    dw(1,2) = power*((r2/This%dr2)**(power-1))*pw1*pw2/This%dr2/wsum**2
    dw(1,3) = power*((r3/This%dr3)**(power-1))*pw1*pw3/This%dr3/wsum**2

    dw(2,1) = power*((r1/This%dr1)**(power-1))*pw2*pw1/This%dr1/wsum**2
    dw(2,2) =-power*((r2/This%dr2)**(power-1))*(pw2*pw1+pw2*pw3)/This%dr2/wsum**2
    dw(2,3) = power*((r3/This%dr3)**(power-1))*pw2*pw3/This%dr3/wsum**2


    dw(3,1) = power*((r1/This%dr1)**(power-1))*pw3*pw1/This%dr1/wsum**2
    dw(3,2) = power*((r2/This%dr2)**(power-1))*pw3*pw2/This%dr2/wsum**2
    dw(3,3) =-power*((r3/This%dr3)**(power-1))*(pw3*pw1+pw3*pw2)/This%dr3/wsum**2

  return

end subroutine weights
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine calcener(This,capr,smlr,theta, ener, der, sno)
!use rep_ker

  implicit none
  
  class(CNO_Basel_PES_Type) ,intent(in)  :: This
  integer, intent(in)                    :: sno
  real(rkp), intent(in)                  :: capr, smlr, theta
  real(rkp), intent(out)                 :: ener
  real(rkp), dimension(:), intent(out)   :: der(3)

  real(rkp)                              :: lambda
  real(rkp),parameter                    :: pi = acos(-1.0d0), piby180 = pi/180.0d0
  real(rkp)                              :: asener, anener, asder, z1, z2
  real(rkp), dimension(:)                :: ander(3), x(3)
  integer                                :: kk, ii

  x(1)=(1.0d0-cos(theta))/2.0d0
  x(2)=capr
  x(3)=smlr

  !print*, 'here 1'
  if (sno==1) then
    asener = 0.0d0
    asder=0.0d0
    do kk = 1,This%na1
      asener = asener + drker26(smlr,This%asy_array1(kk,1))*This%asy_array1(kk,2)
      asder = asder + ddrker26(smlr,This%asy_array1(kk,1))*This%asy_array1(kk,2)
    end do

    anener = 0.0d0
    ander=0.0d0
    call This%pes1%evaluate_fast(x,anener,ander)

      ener = asener+anener
      der(1) = ander(2)
      der(2) = asder+ander(3)
      der(3) = ander(1)*sin(theta)/2.0d0
  end if

  if (sno==2) then
    asener = 0.0d0
    asder=0.0d0
    do kk = 1,This%na2
      asener = asener + drker26(smlr,This%asy_array2(kk,1))*This%asy_array2(kk,2)
      asder = asder + ddrker26(smlr,This%asy_array2(kk,1))*This%asy_array2(kk,2)
    end do

    anener = 0.0d0
    ander=0.0d0
    call This%pes2%evaluate_fast(x,anener,ander)

      ener = asener+anener
      der(1) = ander(2)
      der(2) = asder+ander(3)
      der(3) = ander(1)*sin(theta)/2.0d0
  end if

  if (sno==3) then
    asener = 0.0d0
    asder=0.0d0
    do kk = 1,This%na3
      asener = asener + drker26(smlr,This%asy_array3(kk,1))*This%asy_array3(kk,2)
      asder = asder + ddrker26(smlr,This%asy_array3(kk,1))*This%asy_array3(kk,2)
    end do

    anener = 0.0d0
    ander=0.0d0
    call This%pes3%evaluate_fast(x,anener,ander)

      ener = asener+anener
      der(1) = ander(2)
      der(2) = asder+ander(3)
      der(3) = ander(1)*sin(theta)/2.0d0
  end if

  return

end subroutine calcener
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
subroutine OCN_PES(This,tmpr, totener, dvdr)
  
  implicit none
  
  class(CNO_Basel_PES_Type) ,intent(in)   ::    This
  real(rkp), dimension (:)  ,intent (in)  ::    tmpr(3)
  real(rkp), dimension (:)  ,intent (out) ::    dvdr(3)
  real(rkp)                 ,intent (out) ::    totener
  real(rkp), dimension (:)                :: xp(3), tmpdvdr(3), r(3)
  real(rkp), parameter                    :: dx = 0.0005d0
  real(rkp)                               :: ener
  integer                                 :: ii

  !      3
  !      O
  !     / \
  !  r3/   \r2
  !   /     \
  !  /       \
  ! C---------N
  ! 1   r1    2

  r(1)=tmpr(1) !CN
  r(2)=tmpr(2) !NO
  r(3)=tmpr(3) !CO

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

end subroutine OCN_PES

end module