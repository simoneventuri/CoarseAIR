! -*-F90-*
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


!!! - Find Diatomic Potential
!!! - Remove Diatomic Contribution
!!! - Orient Derivative of Diatomic in Cartesian Coordinates
!!! - Test it 



Module N2O_Basel_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only: rkp, Zero, One, Two, B_To_Ang, Kcm_To_Hartree, KcmAng_To_HartB, cmm1_To_Hartree, Pi
  use PES_Class             ,only: PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only: Logger
  use Error_Class           ,only: Error

  implicit none

  !-------------------------!
  ! Default: make everything private
  private
  public :: N2O_Basel_PES_Type

  !-------------------------!
  type, extends(PES_Type) :: N2O_Basel_PES_Type
  integer                   :: iN2
  integer                   :: iNO
  integer                   :: jNO
  integer                   :: OIdx
  integer                   :: NIdx1
  integer                   :: NIdx2

  !-------------------------------------------------------------------------------------      
  ! The dimensions of the PES grid is specified here (number of points for the three coordinates)
  integer                 :: nAlpha  !number of gridpoints for alpha
  integer                 :: nVdW    !25,&! !number of gridpoints for bigR
  integer                 :: nR      !number of gridpoints for smallR
  !-------------------------------------------------------------------------------------  
                 
  !Array to store the grid data
  real(rkp) ,dimension(:,:) ,allocatable       :: dataArray1, dataArray2, dataArray3
  real(rkp) ,dimension(:)   ,allocatable       :: r_grid1, r_grid2, r_grid3                   
    
  !Arrays for the weights (coefficients)
  real(rkp) ,dimension(:) ,allocatable       :: w_tot1, w_R1
  real(rkp) ,dimension(:) ,allocatable       :: w_tot2, w_R2
  real(rkp) ,dimension(:) ,allocatable       :: w_tot3, w_R3

  !####CHANGE CUTOFF FOR SWITCHING FUNCTION HERE####
  real(rkp)  :: dxSwitch  !determines the "cutoff" of the switching
  real(rkp)  :: aSwitch   !determines plateau of switching function
  contains
  procedure :: Initialize     => Initialize_N2O_Basel_PES
  procedure :: Output         => Output_N2O_Basel_PES
  procedure :: Compute        => Compute_N2O_Basel_PES_1d
  procedure :: Potential      => N2O_Basel_Potential_From_R
  procedure :: DiatPotential  => N2O_Basel_Potential_From_R_OnlyDiat
  procedure :: TriatPotential => N2O_Basel_Potential_From_R_OnlyTriat
  procedure :: Switchfac      => Switchfac_Basel
  procedure :: asymp          => asymp_Basel
  procedure :: dasymp         => dasymp_Basel
  procedure :: Epot           => Epot_Basel
  procedure :: adEdVdW        => adEdVdW_Basel
  procedure :: adEdR          => adEdR_Basel
  procedure :: adEdalpha      => adEdalpha_Basel
  procedure :: dEdalpha       => dEdalpha_Basel
  procedure :: dEdVdW         => dEdVdW_Basel
  procedure :: CalcWeights    => CalcWeights_Basel
  procedure :: dEdR           => dEdR_Basel
  procedure :: USERE          => USERE_Basel
  end type
  
  logical                         ,parameter    ::    i_Debug_Global = .False.
  
  !-------------------------!
  ! Local (private) parameters

  real(rkp) ,dimension(:)   ,allocatable :: w_tot
  real(rkp) ,dimension(:,:) ,allocatable :: dataArray
  real(rkp) ,dimension(:)   ,allocatable :: r_grid
  real(rkp) ,dimension(:)   ,allocatable :: w_R

  logical :: printe_flag = .false.
  logical :: printe_off  = .false.

  contains

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for UMN PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_N2O_Basel_PES( This, Input, Atoms, iPES, i_Debug )
  
  use Input_Class                     ,only:  Input_Type
  use Atom_Class                      ,only:  Atom_Type
  use DiatomicPotential_Factory_Class  ,only:  DiatomicPotential_Factory_Type


  class(N2O_Basel_PES_Type), intent(out)    :: This
  type(Input_Type), intent(in)              :: Input
  type(Atom_Type) ,dimension(:), intent(in) :: Atoms
  integer, intent(in)                       :: iPES
  logical, optional ,intent(in)             :: i_Debug
  
  integer                                   :: iP
  integer                                   :: k
  integer                                   :: ii, jA
  integer                                   :: Status
  integer                                   :: Unit
  integer, dimension(6)                     :: MatrixTemp = (/ 1, 2, 1, 3, 2, 3 /)
  character(:), allocatable                 :: Basel_N2O_folder
  character(*), parameter                   :: Name_PES = 'N2O_Basel'
  character(80)                             :: line
  integer         ,dimension(3,2)           :: iA
  type(DiatomicPotential_Factory_Type)       :: DiatPotFactory
  logical                                   :: i_Debug_Loc
  integer                                   :: i
  integer                                   :: IOS, IOS1, IOS2 !used to check file status

  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N2O_Basel_PES" )
  !i_Debug_Loc   =     Logger%On()
   
  This%Name         = Name_PES
  This%Model        = Input%PES_Model(iPES)
  This%Initialized  = .True.
  This%CartCoordFlg = .True.
  This%NPairs       = 3                 ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  ! allocate( This%mMiMn(3) )
  ! This%mMiMn(1:2) = - Atoms(1:2)%Mass / Atoms(3)%Mass 
  ! if (i_Debug_Loc) call Logger%Write( "This%mMiMn = ", This%mMiMn )
  
  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]

  ii = 1
  do jA = 1,3
    if (trim(adjustl(Atoms(jA)%Name)) == 'O') then
      This%OIdx  = jA
      if (i_Debug_Loc) call Logger%Write( 'O  Atom Nb = ', This%OIdx )
    elseif (trim(adjustl(Atoms(jA)%Name)) == 'N') then
      if (ii==1) then
        This%NIdx1 = jA
        if (i_Debug_Loc) call Logger%Write( 'N1 Atom Nb = ', This%NIdx1 )
        ii = ii + 1
      else
        This%NIdx2 = jA
        if (i_Debug_Loc) call Logger%Write( 'N2 Atom Nb = ', This%NIdx2 )
      end if
    end if
  end do

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiatPotFactory%Construct" )
  ii=1
  do iP = 1,This%NPairs
    call DiatPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )

    if ( ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+1))) ) .ne. "O") .and. ( trim(adjustl(Input%AtomsName(MatrixTemp((iP-1)*2+2)))) .ne. "O" ) ) then
      This%iN2 = iP
      if (i_Debug_Loc) call Logger%Write( 'N2 pair =',iP )
    else
      if (ii == 1) then
      This%iNO=iP
      if (i_Debug_Loc) call Logger%Write( 'First NO pair =',iP )
      else
      This%jNO=iP
      if (i_Debug_Loc) call Logger%Write( 'Second NO pair =',iP )
      end if
      ii=ii+1
    end if

  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
 ! ==============================================================================================================

  
  if (adjustl(trim(This%Model)) == 'Basel_13A1') then
    This%nAlpha   = 11    !number of gridpoints for alpha
    This%nVdW     = 19    !25,&! !number of gridpoints for bigR
    This%nR       = 13    !number of gridpoints for smallR
    This%dxSwitch = 1.d0  !determines the "cutoff" of the switching
    This%aSwitch  = 4.d0  !determines plateau of switching function
  elseif (adjustl(trim(This%Model)) == 'Basel_13A2') then
    This%nAlpha   = 11    !number of gridpoints for alpha
    This%nVdW     = 18    !25,&! !number of gridpoints for bigR
    This%nR       = 12    !number of gridpoints for smallR
    This%dxSwitch = 1.2d0 !determines the "cutoff" of the switching
    This%aSwitch  = 4.d0  !determines plateau of switching function
  else
    write(*,*) 'PES_Model different from Basel or Basel_13A2'
    stop
  end if

  if (.not. allocated(dataArray)) allocate(dataArray( This%nAlpha*This%nVdW*This%nR,4))
  allocate(This%dataArray1(This%nAlpha*This%nVdW*This%nR,4))
  allocate(This%dataArray2(This%nAlpha*This%nVdW*This%nR,4))
  allocate(This%dataArray3(This%nAlpha*This%nVdW*This%nR,4))

  if (.not. allocated(r_grid)) allocate(r_grid( This%nR))                 
  allocate(This%r_grid1(This%nR))
  allocate(This%r_grid2(This%nR))
  allocate(This%r_grid3(This%nR))              
    
  if (.not. allocated(w_tot)) allocate(w_tot( This%nAlpha*This%nR*This%nVdW))
  allocate(This%w_tot1(This%nAlpha*This%nR*This%nVdW))
  allocate(This%w_tot2(This%nAlpha*This%nR*This%nVdW))
  allocate(This%w_tot3(This%nAlpha*This%nR*This%nVdW))

  if (.not. allocated(w_R)) allocate(w_R( This%nR))
  allocate(This%w_R1(This%nR))
  allocate(This%w_R2(This%nR))
  allocate(This%w_R3(This%nR))


  Basel_N2O_folder = trim(adjustl(Input%DtbPath))  // '/' // trim(adjustl(Input%System)) // '/PESs/' // trim(adjustl(Input%PES_Model(iPES))) // '/data/'
  if (i_Debug_Loc) call Logger%Write( "Folder for Reading Basel N2O Parameters: ", Basel_N2O_folder )

  if (i_Debug_Loc) call Logger%Write( "Calling External Module for Initializing PES and Reading Basel N2O Parameters: ", Basel_N2O_folder )


  call readData2(20,IOS ,This%dataArray1, Basel_N2O_folder // "/data1.txt",This%nAlpha*This%nVdW*This%nR,4)
  call readData2(20,IOS1,This%dataArray2, Basel_N2O_folder // "/data2.txt",This%nAlpha*This%nVdW*This%nR,4)
  call readData2(20,IOS2,This%dataArray3, Basel_N2O_folder // "/data3.txt",This%nAlpha*This%nVdW*This%nR,4)
  if ((IOS /= 0).or.(IOS1 /= 0).or.(IOS2 /= 0)) then
    print*, "One of the data files that contains the grid in the data folder could not be read properly."
    print*, "program terminated."
    stop
  end if 
      
  !change alpha to new coordinate y=(1-cos(alpha))/2
  do i=1,This%nR*This%nVdW*This%nAlpha
    This%dataArray1(i,1)=y(This%dataArray1(i,1))
    This%dataArray2(i,1)=y(This%dataArray2(i,1))
    This%dataArray3(i,1)=y(This%dataArray3(i,1))
  end do
  !read in r_grid (needed for asymptote correction)
  do i=1,This%nR
    This%r_grid1(i)=This%dataArray1(1+(i-1)*This%nVdW,2)
    This%r_grid2(i)=This%dataArray2(1+(i-1)*This%nVdW,2)
    This%r_grid3(i)=This%dataArray3(1+(i-1)*This%nVdW,2)
  end do
  
  !read weight/coefficient data for surface 1 from files
  call readData(21,IOS1,This%w_R1,   Basel_N2O_folder // "/w_r1.dat",This%nR)
  call readData(22,IOS2,This%w_tot1, Basel_N2O_folder // "/w_tot1.dat",This%nAlpha*This%nR*This%nVdW)    
  !Error check
  if ((IOS1 /= 0).or.(IOS2 /= 0).or.(IOS /= 0)) then
    print*, "One of the weight-files in the data folder for surface 1 could not be read properly."
    print*, "program terminated."
    stop
  end if    

  !read weight/coefficient data for surface 2 from files
  call readData(21,IOS1,This%w_R2,   Basel_N2O_folder // "/w_r2.dat",This%nR)
  call readData(22,IOS2,This%w_tot2, Basel_N2O_folder // "/w_tot2.dat",This%nAlpha*This%nR*This%nVdW)
  !Error check
  if ((IOS1 /= 0).or.(IOS2 /= 0).or.(IOS /= 0)) then
    print*, "One of the weight-files in the data folder for surface 2 could not be read properly."
    print*, "program terminated."
    stop
  end if 
  
  !read weight/coefficient data for surface 3 from files
  call readData(21,IOS1,This%w_R3,   Basel_N2O_folder // "/w_r3.dat",This%nR)
  call readData(22,IOS2,This%w_tot3, Basel_N2O_folder // "/w_tot3.dat",This%nAlpha*This%nR*This%nVdW)
  !Error check
  if ((IOS1 /= 0).or.(IOS2 /= 0).or.(IOS /= 0)) then
    print*, "One of the weight-files in the data folder for surface 3 could not be read properly."
    print*, "program terminated."
    stop
  end if  



  if (i_Debug_Loc) call Logger%Exiting()
 
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_N2O_Basel_PES( This, Unit )

  class(N2O_Basel_PES_Type)               ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)")  This%Name
  write(Unit,"('PES Model: ',g0)") This%Model
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N2O_Basel_Potential_From_R( This, R, Q ) result( V )
  
  class(N2O_Basel_PES_Type)                     ,intent(in)  ::    This   !< PES class       
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    RTemp
  real(rkp) ,dimension(3)                                    ::    X, Y, Z
  real(rkp)                                                  ::    EU
  real(rkp) ,dimension(3)                                    ::    DX, DY, DZ
  logical                                                    ::    QECONT
  real(rkp) ,dimension(3)                                    ::    ECONT
  real(rkp)                                                  ::    VDiat, VTriat

  X = [Q((This%NIdx1-1)*3+1), Q((This%OIdx-1)*3+1), Q((This%NIdx2-1)*3+1)] 
  Y = [Q((This%NIdx1-1)*3+2), Q((This%OIdx-1)*3+2), Q((This%NIdx2-1)*3+2)]
  Z = [Q((This%NIdx1-1)*3+3), Q((This%OIdx-1)*3+3), Q((This%NIdx2-1)*3+3)]

  ! RTemp(1) = sqrt( (Q(1)  - Q(4))**2 + (Q(2)  - Q(5))**2 + (Q(3)  - Q(6))**2 )
  ! RTemp(2) = sqrt( (Q(1)  - Q(7))**2 + (Q(2)  - Q(8))**2 + (Q(3)  - Q(9))**2 )
  ! RTemp(3) = sqrt( (Q(7)  - Q(4))**2 + (Q(8)  - Q(5))**2 + (Q(9)  - Q(6))**2 )

  ! !-------------------------!
  ! ! Evaluate 2-body interactions
  ! VDiat =         This%Pairs(1)%Vd%DiatomicPotential(RTemp(1))
  ! VDiat = VDiat + This%Pairs(2)%Vd%DiatomicPotential(RTemp(2))
  ! VDiat = VDiat + This%Pairs(3)%Vd%DiatomicPotential(RTemp(3))
  ! V     = VDiat

  call This%USERE(EU, X*B_To_Ang, Y*B_To_Ang, Z*B_To_Ang, DX, DY, DZ, QECONT, ECONT, 3)
  V = EU * Kcm_To_Hartree
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N2O_Basel_Potential_From_R_OnlyDiat( This, R, Q ) result( V )

  class(N2O_Basel_PES_Type)                     ,intent(in)  ::    This   !< PES class
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree]. 

  real(rkp) ,dimension(3)                                    ::    X, Y, Z
  real(rkp) ,dimension(3)                                    ::    RTemp


  X = [Q((This%NIdx1-1)*3+1), Q((This%OIdx-1)*3+1), Q((This%NIdx2-1)*3+1)]
  Y = [Q((This%NIdx1-1)*3+2), Q((This%OIdx-1)*3+2), Q((This%NIdx2-1)*3+2)]
  Z = [Q((This%NIdx1-1)*3+3), Q((This%OIdx-1)*3+3), Q((This%NIdx2-1)*3+3)]

  RTemp(1) = sqrt( (X(This%OIdx)  - X(This%NIdx1))**2 + (Y(This%OIdx)  - Y(This%NIdx1))**2 + (Z(This%OIdx)  - Z(This%NIdx1))**2 )
  RTemp(2) = sqrt( (X(This%OIdx)  - X(This%NIdx2))**2 + (Y(This%OIdx)  - Y(This%NIdx2))**2 + (Z(This%OIdx)  - Z(This%NIdx2))**2 )
  RTemp(3) = sqrt( (X(This%NIdx2) - X(This%NIdx1))**2 + (Y(This%NIdx2) - Y(This%NIdx1))**2 + (Z(This%NIdx2) - Z(This%NIdx1))**2 )

  !-------------------------!
  ! Form and sum contributions of 2-body interactions 
  V = This%Pairs(1)%Vd%DiatomicPotential(RTemp(1))
  V = V + This%Pairs(2)%Vd%DiatomicPotential(RTemp(2))
  V = V + This%Pairs(3)%Vd%DiatomicPotential(RTemp(3))
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N2O_Basel_Potential_From_R_OnlyTriat( This, R, Q ) result( V )

  class(N2O_Basel_PES_Type)                     ,intent(in)  ::    This   !< PES class
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    X, Y, Z
  real(rkp)                                                  ::    EU, VDiat
  real(rkp) ,dimension(3)                                    ::    DX, DY, DZ
  logical                                                    ::    QECONT
  real(rkp) ,dimension(3)                                    ::    ECONT    
  real(rkp) ,dimension(3)                                    ::    RTemp


  X = [Q((This%NIdx1-1)*3+1), Q((This%OIdx-1)*3+1), Q((This%NIdx2-1)*3+1)]
  Y = [Q((This%NIdx1-1)*3+2), Q((This%OIdx-1)*3+2), Q((This%NIdx2-1)*3+2)]
  Z = [Q((This%NIdx1-1)*3+3), Q((This%OIdx-1)*3+3), Q((This%NIdx2-1)*3+3)]


  RTemp(1) = sqrt( (X(This%OIdx)  - X(This%NIdx1))**2 + (Y(This%OIdx)  - Y(This%NIdx1))**2 + (Z(This%OIdx)  - Z(This%NIdx1))**2 )
  RTemp(2) = sqrt( (X(This%OIdx)  - X(This%NIdx2))**2 + (Y(This%OIdx)  - Y(This%NIdx2))**2 + (Z(This%OIdx)  - Z(This%NIdx2))**2 )
  RTemp(3) = sqrt( (X(This%NIdx2) - X(This%NIdx1))**2 + (Y(This%NIdx2) - Y(This%NIdx1))**2 + (Z(This%NIdx2) - Z(This%NIdx1))**2 )

  !-------------------------!
  ! Form and sum contributions of 2-body interactions 
  VDiat = This%Pairs(This%iNO)%Vd%DiatomicPotential(RTemp(1))
  VDiat = VDiat + This%Pairs(This%jNO)%Vd%DiatomicPotential(RTemp(2))
  VDiat = VDiat + This%Pairs(This%iN2)%Vd%DiatomicPotential(RTemp(3))

  call This%USERE(EU, X*B_To_Ang, Y*B_To_Ang, Z*B_To_Ang, DX, DY, DZ, QECONT, ECONT, 3)
  V = EU * Kcm_To_Hartree - VDiat

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_N2O_Basel_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(N2O_Basel_PES_Type)                     ,intent(in)  ::    This   !< PES class
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp) ,dimension(3)                                    ::    RTemp
  real(rkp) ,dimension(3)                                    ::    X, Y, Z
  real(rkp)                                                  ::    EU
  real(rkp) ,dimension(3)                                    ::    DX, DY, DZ
  logical                                                    ::    QECONT
  real(rkp) ,dimension(3)                                    ::    ECONT
  real(rkp)                                                  ::    VDiat, VTriat

  X = [Q((This%NIdx1-1)*3+1), Q((This%OIdx-1)*3+1), Q((This%NIdx2-1)*3+1)]
  Y = [Q((This%NIdx1-1)*3+2), Q((This%OIdx-1)*3+2), Q((This%NIdx2-1)*3+2)]
  Z = [Q((This%NIdx1-1)*3+3), Q((This%OIdx-1)*3+3), Q((This%NIdx2-1)*3+3)]

  RTemp(1) = sqrt( (X(This%OIdx)  - X(This%NIdx1))**2 + (Y(This%OIdx)  - Y(This%NIdx1))**2 + (Z(This%OIdx)  - Z(This%NIdx1))**2 )
  RTemp(2) = sqrt( (X(This%OIdx)  - X(This%NIdx2))**2 + (Y(This%OIdx)  - Y(This%NIdx2))**2 + (Z(This%OIdx)  - Z(This%NIdx2))**2 )
  RTemp(3) = sqrt( (X(This%NIdx2) - X(This%NIdx1))**2 + (Y(This%NIdx2) - Y(This%NIdx1))**2 + (Z(This%NIdx2) - Z(This%NIdx1))**2 )

  !-------------------------!
  ! Form and sum contributions of 2-body interactions 
  VDiat = This%Pairs(This%iNO)%Vd%DiatomicPotential(RTemp(1))
  VDiat = VDiat + This%Pairs(This%jNO)%Vd%DiatomicPotential(RTemp(2))
  VDiat = VDiat + This%Pairs(This%iN2)%Vd%DiatomicPotential(RTemp(3))

  V     = VDiat * Kcm_To_Hartree

  ! call This%USERE(EU, X*B_To_Ang, Y*B_To_Ang, Z*B_To_Ang, DX, DY, DZ, QECONT, ECONT, 3)
  ! V       = EU * Kcm_To_Hartree
  ! dVdQ(1) = DX(1)
  ! dVdQ(4) = DX(2)
  ! dVdQ(7) = DX(3)
  ! dVdQ(2) = DY(1)
  ! dVdQ(5) = DY(2)
  ! dVdQ(8) = DY(3)
  ! dVdQ(3) = DZ(1)
  ! dVdQ(6) = DZ(2)
  ! dVdQ(9) = DZ(3)
  ! dVdQ    = dVdQ * Kcm_To_Hartree / B_To_Ang

  dVdR = Zero

  dVdQ = Zero
  call This%TransToCart_3Atoms( R, Q, dVdR, dVdQ)
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!   atom 1: N
!   atom 2: O
!   atom 3: N
!
! The PES for the NO+N in the ^3A' state is call 
!
!     call  call initialize_PES()
!     call USERE(EU,X,Y,Z,DX,DY,DZ,QECONT,ECONT,NATOMX)
!
!     EU   - User energy to be returned
!     X,Y,Z - Coordinates for energy evaluation
!     DX,DY,DZ - Forces. Must be modified. Append (dE/dX,...)
!     QECONT - Flag for analysis (0-no analysis,>0=analysis)
!     ECONT(NATOMX) - Analysis array to be filled if QECONT>0.
!     NATOMX - Number of atoms = 3
!
!     
  
Subroutine Switchfac_Basel(This, x,f,dfdx)
  !This is the subroutine to calculate the switching factors
  
  class(N2O_Basel_PES_Type) ,intent(in)  :: This
  real(rkp), intent(in)                  :: x
  real(rkp), intent(out)                 :: f
  real(rkp), intent(out), optional       :: dfdx


  !HERE THE SWITCHING FUNCTION IS SPECIFIED. BE SURE TO ALSO CHANGE THE DERIVATIVE
  !SHOULD YOU CHANGE IT!!!

  !new function
  f= dexp(-(x/This%dxSwitch)**This%aSwitch)

  if (present(dfdx)) dfdx= -(This%aSwitch/x)*f*(x/This%dxSwitch)**This%aSwitch
     
End Subroutine 


!calculates the new asymptote depending for r'
Function asymp_Basel(This, R, surf) result( asymp )
  
  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: R
  integer,   intent(in) :: surf
  integer               :: counter
  real(rkp)             :: asymp

  asymp = 0.d0
  
  !read appropriate weights, depending on which surface is selected
  if      (surf == 1) then
    w_R    = This%w_R1
    r_grid = This%r_grid1
  else if (surf == 2) then
    w_R    = This%w_R2
    r_grid = This%r_grid2
  else if (surf == 3) then
    w_R    = This%w_R3
    r_grid = This%r_grid3
  else
    print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
    stop
  end if 
  do counter=1,This%nR
    asymp = asymp + w_R(counter)*q_ker(R, r_grid(counter))
  end do        
end function



Function dasymp_Basel(This, R, surf) result( dasymp )

  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: R
  !which surface shall be calculated
  integer, intent(in)   :: surf
  integer               :: counter
  real(rkp)             :: dasymp

  dasymp = 0.d0
  
      !read appropriate weights, depending on which surface is selected
  if      (surf == 1) then
    w_R    = This%w_R1
    r_grid = This%r_grid1
  else if (surf == 2) then
    w_R    = This%w_R2
    r_grid = This%r_grid2
  else if (surf == 3) then
    w_R    = This%w_R3
    r_grid = This%r_grid3
  else
    print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
    stop
  end if 
  do counter=1,This%nR
    dasymp = dasymp + w_R(counter)*dq_ker(R, r_grid(counter))
  end do   

end function



!calculates the energy V(alpha',r',R') for an off-grid point
Function Epot_Basel(This, R,VdW, alpha, surf) result( Epot )
  
  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  !input coordinates
  real(rkp), intent(in) :: R,VdW,alpha
  integer,   intent(in) :: surf
  real(rkp)             :: Epot

  real(rkp)             :: y_val !coordinate transformation
  !which surface shall be calculated?
  integer               :: counter

  y_val = y(alpha)    !transforms alpha to y
  Epot = 0d0

  !read appropriate weights, depending on which surface is selected
  if      (surf == 1) then
    w_tot     = This%w_tot1
    dataArray = This%dataArray1
  else if (surf == 2) then
    w_tot     = This%w_tot2
    dataArray = This%dataArray2
  else if (surf == 3) then
    w_tot     = This%w_tot3
    dataArray = This%dataArray3
  else
    print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
    stop
  end if 
  
  
  do counter=1,This%nAlpha*This%nR*This%nVdW
    Epot = Epot + w_tot(counter)*k_ker(y_val,dataArray(counter,1)) * q_ker(R,dataArray(counter,2))*q_ker(VdW,dataArray(counter,3))
  end do
  Epot = Epot + This%asymp(R,surf)

end function



!calculates the partial derivative of V(alpha,r,R) with respect to VdW analytically
Function adEdVdW_Basel(This, R,VdW,alpha,surf) result(adEdVdW)
  
  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: alpha,R,VdW
  real(rkp) :: y_val
  integer, intent(in) :: surf 
  integer :: counter
  real(rkp)             :: adEdVdW

  y_val = y(alpha)    !transforms alpha to y
  adEdVdW = 0d0
  !read appropriate weights, depending on which surface is selected
  if      (surf == 1) then
    w_tot = This%w_tot1
    dataArray = This%dataArray1
  else if (surf == 2) then
    w_tot = This%w_tot2
    dataArray = This%dataArray2
  else if (surf == 3) then
    w_tot = This%w_tot3
    dataArray = This%dataArray3
  else
    print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
    stop
  end if 
  
  do counter=1,This%nAlpha*This%nR*This%nVdW
    adEdVdW = adEdVdW + w_tot(counter)*k_ker(y_val, dataArray(counter,1)) * q_ker(R, dataArray(counter,2))*dq_ker(VdW, dataArray(counter,3))
  end do

end function



!calculates the partial derivative of V(alpha,r,R) with respect to r analytically
Function adEdR_Basel(This, R,VdW,alpha,surf) result(adEdR)
  
  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: alpha,R,VdW
  real(rkp) :: y_val
  integer, intent(in) :: surf 
  integer :: counter
  real(rkp)             :: adEdR

  y_val = y(alpha)    !transforms alpha to y
  adEdR = 0d0
  !read appropriate weights, depending on which surface is selected
  if      (surf == 1) then
    w_tot = This%w_tot1
    dataArray = This%dataArray1
  else if (surf == 2) then
    w_tot = This%w_tot2
    dataArray = This%dataArray2
  else if (surf == 3) then
    w_tot = This%w_tot3
    dataArray = This%dataArray3
  else
    print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
    stop
  end if 

  do counter=1,This%nAlpha*This%nR*This%nVdW
    adEdR = adEdR + w_tot(counter)*k_ker(y_val, dataArray(counter,1)) * dq_ker(R, dataArray(counter,2))*q_ker(VdW, dataArray(counter,3))
  end do
  adEdR = adEdR + This%dasymp(R,surf)

end function



!calculates the partial derivative of V(alpha,r,R) with respect to alpha analytically
Function adEdalpha_Basel(This, R,VdW,alpha,surf) result(adEdalpha)
  
  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: alpha,R,VdW
  real(rkp) :: y_val,dy_val
  integer, intent(in) :: surf 
  integer :: counter
  real(rkp) :: adEdalpha

  y_val = y(alpha)    !transforms alpha to y
  dy_val = dy(alpha)
  adEdalpha = 0d0
  
  !read appropriate weights, depending on which surface is selected
  if      (surf == 1) then
    w_tot     = This%w_tot1
    dataArray = This%dataArray1
  else if (surf == 2) then
    w_tot     = This%w_tot2
    dataArray = This%dataArray2
  else if (surf == 3) then
    w_tot     = This%w_tot3
    dataArray = This%dataArray3
  else
    print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
    stop
  end if 
  
  do counter=1,This%nAlpha*This%nR*This%nVdW
    adEdalpha = adEdalpha + w_tot(counter)*dk_ker(y_val, dataArray(counter,1))*dy_val * q_ker(R, dataArray(counter,2))*q_ker(VdW, dataArray(counter,3))
  end do

end function


!calculates the partial derivative of V(alpha,r,R) with respect to VdW numerically
Function dEdVdW_Basel(This, R,VdW,alpha,surf) result(dEdVdW)

  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: alpha,R,VdW
  integer, intent(in) :: surf    
  real(rkp) :: eps
  real(rkp), dimension(4) :: stencil
  real(rkp) :: dEdVdW
  !construct epsilon
  !eps = dabs(dsqrt(epsilon(0.d0))*VdW)
  !if(eps < dsqrt(epsilon(0.d0))) eps = dsqrt(epsilon(0.d0))
  eps = 1d-5
  !build stencil
  stencil(1) = This%Epot(R, VdW +2*eps, alpha, surf)
  stencil(2) = This%Epot(R, VdW +1*eps, alpha, surf)
  stencil(3) = This%Epot(R, VdW -1*eps, alpha, surf)
  stencil(4) = This%Epot(R, VdW -2*eps, alpha, surf)
  dEdVdW = derivative_FivePointStencil(stencil,eps)   

end function



!calculates the partial derivative of V(alpha,r,R) with respect to r numerically
Function dEdR_Basel(This, R,VdW,alpha,surf) result(dEdR)

  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: alpha,R,VdW  
  integer, intent(in) :: surf 
  real(rkp) :: eps
  real(rkp), dimension(4) :: stencil
  real(rkp) :: dEdR
  !construct epsilon
  !eps = dabs(dsqrt(epsilon(0.d0))*R)
  !if(eps < dsqrt(epsilon(0.d0))) eps = dsqrt(epsilon(0.d0))
  eps = 1d-5
  !build stencil
  stencil(1) = This%Epot(R +2*eps, VdW, alpha, surf)
  stencil(2) = This%Epot(R +1*eps, VdW, alpha, surf)
  stencil(3) = This%Epot(R -1*eps, VdW, alpha, surf)
  stencil(4) = This%Epot(R -2*eps, VdW, alpha, surf)
  dEdR = derivative_FivePointStencil(stencil,eps)   

end function



!calculates the partial derivative of V(alpha,r,R) with respect to alpha numerically
Function dEdalpha_Basel(This, R,VdW,alpha,surf) result(dEdalpha)

  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp), intent(in) :: alpha,R,VdW
  integer, intent(in) :: surf 
  real(rkp) :: eps, tempy, yplusdy, tempalpha
  real(rkp), dimension(4) :: stencil
  real(rkp) :: dEdalpha
  !construct epsilon
  !eps = dabs(dsqrt(epsilon(0.d0))*alpha)
  !if(eps < dsqrt(epsilon(0.d0))) eps = dsqrt(epsilon(0.d0))
  eps = 1d-5        
  !This is a long block to construct an appropriate displacement for alpha
  tempy = y(alpha)
  yplusdy = tempy + eps
  if (yplusdy > 1.d0) yplusdy = tempy - eps 
  tempalpha = dacos(1-2*yplusdy)
  eps = dabs(alpha-tempalpha)
  !build stencil            
  stencil(1) = This%Epot(R, VdW, alpha +2*eps, surf)
  stencil(2) = This%Epot(R, VdW, alpha +1*eps, surf)
  stencil(3) = This%Epot(R, VdW, alpha -1*eps, surf)
  stencil(4) = This%Epot(R, VdW, alpha -2*eps, surf)
  dEdalpha = derivative_FivePointStencil(stencil,eps)

end function


SUBROUTINE USERE_Basel(This, EU,X,Y,Z,DX,DY,DZ,QECONT,ECONT,NATOMX)
  !
  !     THE USER ENERGY ROUTINE. DOES NOTHING IN THE NORMAL SYSTEM
  !     EXCEPT SET EU TO ZERO. QECONT IS A FLAG WHICH WILL BE SET TO 1
  !     WHEN THE ROUTINE IS CALLED IN THE QECONTIS SECTION. INDIVIDUAL
  !     ATOM USER ENERGIES ARE THEN RETURNED IN ECONT, AND THE DERIVATIVE
  !     ARRAYS SHOULD NOT BE ACCESSED.
  !
  !     EU   - User energy to be returned
  !     X,Y,Z - Coordinates for energy evaluation
  !     DX,DY,DZ - Forces. Must be modified. Append (dE/dX,...)
  !     QECONT - Flag for analysis (0-no analysis,>0=analysis)
  !     ECONT(NATOMX) - Analysis array to be filled if QECONT>0.
  !     NATOMX - Number of atoms
  !
  !     Author: Robert Bruccoleri
  !
  
!*****************************************************************************************
! END OF MODIFICATIONS
!*****************************************************************************************  

  class(N2O_Basel_PES_Type) ,intent(in)  ::    This

  real(rkp) EU
  INTEGER NATOMX
  real(rkp) X(NATOMX),Y(NATOMX),Z(NATOMX)
  real(rkp) DX(NATOMX),DY(NATOMX),DZ(NATOMX)
  LOGICAL QECONT
  real(rkp) ECONT(NATOMX)
  !
  INTEGER I
  !
  !EU=0.0  !***THIS WAS ORIGINALLY IN CHARMM, BUT I COMMENTED OUT TO CONTINUE DECLARATIONS
  
!*****************************************************************************************
! BEGIN OF MODIFICATIONS
!*****************************************************************************************

  !loop variable (integer i is already defined by standard charmm code)
  integer j
  
  !conversion from eV to kcal/mol
  real(rkp), parameter :: ev2kcalmol = 23.06035d0 
  
  !Array that stores masses. This is needed for some things like calculating the centre of
  !mass, needed to calculate Jacobi coordinates. YOU NEED TO CHANGE THIS WHERE IT'S DEFINED!
  real(rkp), dimension(NATOMX)  :: mass
  
  !Output coordinates (Jacobi)
  !1st index: surface
  !2nd index: Jacobi coordinate (r, bigR, alpha)
  real(rkp), dimension(3,3) :: jacCoord
  
  !Array for the derivatives of Jacobi coordinates with respect to the cartesian coordinates
  !1st index: surface
  !2nd index: Jacobi coordinate (r,bigR,alpha)
  !3rd index: cartesian coordinate (x1,y1,z1,x2,y2,z2,x3,y3,z3)
  real(rkp), dimension(3,3,9) :: dJacdCart
  
  !Array for the derivatives of the potential with respect to the Jacobi coordinates 
  !1st index: surface
  !2nd index: Jacobi coordinate (r,bigR,alpha)
  real(rkp), dimension(3,3) :: dEdJac
  
  !All this is just needed for the PES switching routine
  !weights for the 3 surfaces, w0: unnormalized, w: normalized
  real(rkp), dimension(3)   :: w0,w  
  !derivatives of raw (unnormalized) weights with respect to r
  real(rkp), dimension(3)   :: dw0dr 
  !derivatives of weight j with respect to r_i (1st index j, 2nd index i)
  real(rkp), dimension(3,3) :: dwdr
  !derivatives of weight j with respect to raw weight i (1st index j, 2nd index i)
  real(rkp), dimension(3,3) :: dwdw0
  !derivatives of weight j with respect to cart. coord. i (1st index j, 2nd index i)
  real(rkp), dimension(3,9) :: dwjdxi 
  !derivatives of potential j with respect to cart. coord. i (1st index j, 2nd index i)
  real(rkp), dimension(3,9) :: dEjdxi
  !temporary array that stores either dwjdxi or dEjdxi, because I need to swap back coordinates
  !for surface 2 and 3!
  real(rkp), dimension(9) :: ddxi
  !potential energies on the "pure" surfaces
  real(rkp), dimension(3)   :: Esurf
                      
!-------------------------------------CHANGE HERE-------------------------------------  
! If this value is set to .true., reactive MD simulations are performed (the PESs are switched) 
! dynamically. If it is set to .false., the simulations will be performed on PES 1 only (bond 
! dissociation of the diatomic is then not possible).
  logical :: use_surface_switching = .true.
!-------------------------------------------------------------------------------------  
  
  !####CHOOSE HERE HOW THE DERIVATIVES OF THE PES SHOULD BE CALCULATED!####
  !####Analytical derivatives were found to be superior in this implementation####
  logical :: use_analytical_derivatives = .true.
  
  !THIS IS AN ERROR CHECK. IF WE DON'T HAVE 3 ATOMS, THIS CODE MAKES NO SENSE AT ALL!
  if(NATOMX /= 3) then
    write(*,*) "WARNING: Number of atoms is not 3."
    write(*,*) "This code is designed for 3-atomic systems ONLY."
    write(*,*) "program was terminated."
    stop
  end if
  
!-------------------------------------CHANGE HERE-------------------------------------        
! Masses of the atoms need to be specified here, because we need to calculate the centre of mass.
! The unit is not important, as long as all masses are given in the same unit.
! The reordering of the masses for the other surfaces is performed automatically.
  mass(1) = 14.007d0  !this corresponds to atom 1 
  mass(2) = 15.9994d0 !this corresponds to atom 2 
  mass(3) = 14.007d0 !this corresponds to atom 3 
!-------------------------------------------------------------------------------------     
  
  call calc_Jacobi_coordinates(X,Y,Z,jacCoord,mass,dJacdCart)
!DEBUG: Jacobi coordinates for all 3 surfaces are calculated fine (was tested)
!and deemed bug-free
!  print*, "Jacobi coordinates: ", "    r    ", "    R    ", "    alpha    "
!  print*, "(surface 1)       : ", jacCoord(1,1), jacCoord(1,2), jacCoord(1,3)
!  print*, "(surface 2)       : ", jacCoord(2,1), jacCoord(2,2), jacCoord(2,3)
!  print*, "(surface 3)       : ", jacCoord(3,1), jacCoord(3,2), jacCoord(3,3)
  
  !print*, "dJacdCart:", dJacdCart
 

if(use_surface_switching) then
!calculate weights and their respective derivatives etc.
  call This%CalcWeights(w,w0,dw0dr,dwdr,dwdw0,jacCoord(:,1))
else
!give surface 1 the full weight
  w    = 0.d0
  w(1) = 1.d0
end if
 
!calculate the derivatives of the potential surfaces 1-3 with respect to the Jacobi
!coordinates (r,R,alpha).
 do i=1,3
   !if a weight is 0, we can just set all the derivatives to 0 and skip expensive
   !calculations
   if(w(i) == 0.d0) then
     dEdJac(i,:) = 0.d0
     cycle 
   end if
   
   if(.not.use_analytical_derivatives) then
   !do it numerically
     dEdJac(i,1) = This%dEdR(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
     dEdJac(i,2) = This%dEdVdW(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
     dEdJac(i,3) = This%dEdalpha(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)    
   else
   !do it analytically
     dEdJac(i,1) = This%adEdR(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
     dEdJac(i,2) = This%adEdVdW(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
     dEdJac(i,3) = This%adEdalpha(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)   
   end if 
 end do

 
  !if we need to switch surfaces, we do it like this
  if(use_surface_switching) then
    !calculate energies on "pure" surfaces
    do j=1,3
      !if a weight is 0, we can just set the energy to 0 and skip expensive calculations
      if(w(j) == 0.d0) then
        Esurf(j) = 0.d0
        cycle 
      end if
      Esurf(j) = This%Epot(jacCoord(j,1),jacCoord(j,2),jacCoord(j,3),j)
    end do
    !Assign mixed energy value:
    EU = real((Esurf(1)*w(1) + Esurf(2)*w(2) + Esurf(3)*w(3))*ev2kcalmol,rkp)
     
    !calculate all 18 dwjdxi and dEjdxi
    do j=1,3
      do i=1,9             
        dEjdxi(j,i) = sum(dEdJac(j,:)*dJacdCart(j,:,i))
        
        
        !calculate dwjdxi by using the raw weight derivatives
        dwjdxi(j,i) = dwdw0(j,1)*dw0dr(1)*dJacdCart(1,1,i) + dwdw0(j,2)*dw0dr(2)*dJacdCart(2,1,i) &
               +dwdw0(j,3)*dw0dr(3)*dJacdCart(3,1,i)
        !print*, "dEjdxi: ", dEjdxi(j,i)
      end do
    end do
      
   !Assign derivatives:
   do i=1,3
     DX(i) = real(sum(w*dEjdxi(:,1+(i-1)*3) + dwjdxi(:,1+(i-1)*3)*Esurf)*ev2kcalmol,rkp)
     DY(i) = real(sum(w*dEjdxi(:,2+(i-1)*3) + dwjdxi(:,2+(i-1)*3)*Esurf)*ev2kcalmol,rkp)
     DZ(i) = real(sum(w*dEjdxi(:,3+(i-1)*3) + dwjdxi(:,3+(i-1)*3)*Esurf)*ev2kcalmol,rkp)
   end do
            
  !if not, we just consider the first surface and that's it  
  else
    !Assign energy:
    EU = real(This%Epot(jacCoord(1,1),jacCoord(1,2),jacCoord(1,3),1)*ev2kcalmol,rkp)
    
    !Assign derivatives:
    do i=1,3
      DX(i) = real((dEdJac(1,1)*dJacdCart(1,1,1+(i-1)*3) + dEdJac(1,2)*dJacdCart(1,2,1+(i-1)*3) &
           + dEdJac(1,3)*dJacdCart(1,3,1+(i-1)*3))*ev2kcalmol,rkp)
      DY(i) = real((dEdJac(1,1)*dJacdCart(1,1,2+(i-1)*3) + dEdJac(1,2)*dJacdCart(1,2,2+(i-1)*3) &
           + dEdJac(1,3)*dJacdCart(1,3,2+(i-1)*3))*ev2kcalmol,rkp)
      DZ(i) = real((dEdJac(1,1)*dJacdCart(1,1,3+(i-1)*3) + dEdJac(1,2)*dJacdCart(1,2,3+(i-1)*3) &
           + dEdJac(1,3)*dJacdCart(1,3,3+(i-1)*3))*ev2kcalmol,rkp)
    end do 
  end if
   
!*****************************************************************************************
! END OF MODIFICATIONS
!*****************************************************************************************  

  IF(QECONT) THEN
   DO I=1,NATOMX
    ECONT(I)=0.0
   ENDDO
  ENDIF
  RETURN
END SUBROUTINE


SUBROUTINE CalcWeights_Basel(This, w,w0,dw0dr,dwdr,dwdw0,r)

  class(N2O_Basel_PES_Type) ,intent(in)  :: This

  real(rkp), dimension(3),   intent(in)  :: r
  real(rkp), dimension(3),   intent(out) :: w,w0,dw0dr
  real(rkp), dimension(3,3), intent(out) :: dwdr
  real(rkp), dimension(3,3), intent(out) :: dwdw0

  real(rkp), dimension(3), save          :: w_save = (/1.d0,0.d0,0.d0/) !this is a safety measure, so that we never get w=0 for all weights!

  real(rkp)                              :: sumw0

  !for loops
  integer i

  do i=1,3
    call This%switchfac(r(i),w0(i),dw0dr(i))
  end do
  sumw0 = sum(w0)
  !this is implemented to prevent division by zero in case we get full atomization
  if(sumw0 < epsilon(0.d0)) sumw0 = 10*epsilon(0.d0)

  !calculate the normalized weights
  w = w0/sumw0

  !check if we have the unfortunate case of all weights being 0. If yes, we load the weights from the previous step.
  if((w(1) == 0.d0).and.(w(2) == 0.d0).and.(w(3) == 0.d0)) then
    w = w_save
  end if

  !store old weights
  w_save = w

  !calculate derivatives of the normalized weights with respect to all different r
  dwdr(1,1)  = (dw0dr(1)*sumw0 - w0(1)*dw0dr(1))/(sumw0**2)
  dwdr(1,2)  = (-w0(1)*dw0dr(2))/(sumw0**2)
  dwdr(1,3)  = (-w0(1)*dw0dr(3))/(sumw0**2)

  dwdr(2,1)  = (-w0(2)*dw0dr(1))/(sumw0**2)
  dwdr(2,2)  = (dw0dr(2)*sumw0 - w0(2)*dw0dr(2))/(sumw0**2)
  dwdr(2,3)  = (-w0(2)*dw0dr(3))/(sumw0**2)

  dwdr(3,1)  = (-w0(3)*dw0dr(1))/(sumw0**2)
  dwdr(3,2)  = (-w0(3)*dw0dr(2))/(sumw0**2)
  dwdr(3,3)  = (dw0dr(3)*sumw0 - w0(3)*dw0dr(3))/(sumw0**2)

  dwdw0(1,1) = (sumw0-w0(1))/(sumw0**2)
  dwdw0(1,2) = (-w0(1))/(sumw0**2)
  dwdw0(1,3) = (-w0(1))/(sumw0**2)

  dwdw0(2,1) = (-w0(2))/(sumw0**2)
  dwdw0(2,2) = (sumw0-w0(2))/(sumw0**2)
  dwdw0(2,3) = (-w0(2))/(sumw0**2)

  dwdw0(3,1) = (-w0(3))/(sumw0**2)
  dwdw0(3,2) = (-w0(3))/(sumw0**2)
  dwdw0(3,3) = (sumw0-w0(3))/(sumw0**2)

End SUBROUTINE 



!***** For Charmm just remove until here and the comments with 3 exclamations (!!!) ***********
!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains the 1-dimensional kernel functions and their derivatives
!
!///////////////////////////////////////////////////////////////////////////////

!the reproducing kernel q_ker for distance like variables
real(rkp) Pure Function q_ker(x,r)
  implicit none
  real(rkp) ,intent(in) :: x,r
  real(rkp)             :: x_l,x_s
  x_l=r
     x_s=x
  if(r < x) then
    x_s=r
    x_l=x
  endif
  q_ker = (1d0/(14d0*x_l**7))*(1d0-7*x_s/(9*x_l))
end function q_ker


!derivative of q_ker
real(rkp) Pure Function dq_ker(x,r)
  implicit none
  real(rkp) ,intent(in) :: x,r
  if(x < r) then
    dq_ker = -1d0/(18d0*r**8)
  else
    dq_ker = -(1d0/18d0)*(9*x-8*r)/(x**9)
  end if
end function dq_ker


!the reproducing kernel k_ker for angle like variables
real(rkp) Pure Function k_ker(x,r)      !k=2 Anglelike
  implicit none
  real(rkp) ,intent(in) :: x,r
  real(rkp)             :: x_l,x_s
  x_l=r !Large value of the coordinate
  x_s=x !Small value of the coordinate
  if(r < x) then !Here the smaller and larger value are chosen
    x_s=r
    x_l=x
  endif
  !Only two terms (up to n=1) were taken of the Gauss´ hypergeometric function expansion  
  k_ker=1+x_s*x_l+2*x_s**2*x_l*(1.d0-x_s/(3*x_l))
end function k_ker


!derivative of k_ker
real(rkp) Pure Function dk_ker(x,r)
  implicit none
  real(rkp) ,intent(in) :: x,r
  if(x < r) then
    dk_ker = r+4*x*r*(1.d0-x/(3*r))-(2*x**2)/3d0
  else
    dk_ker = r+2*r**2*(1.d0-r/(3*x))+(2*r**3)/(3*x)
  end if
end function dk_ker



!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains routines to read files and store them in arrays
!
!///////////////////////////////////////////////////////////////////////////////

!subroutine to read in files
subroutine readData(unitnr,ios,array, inputfile,dim1)
  implicit none
  !declaration of variables
  integer :: unitnr,dim1
  integer :: i,ios
  real(rkp) :: array(dim1)
  character(len=*) :: inputfile

  open(unitnr,file=inputfile, status='old', action='read', &
       iostat=ios)
  !error check
  if (ios /= 0) then
    return
  end if 
  do i=1,dim1
    read(unitnr,*) array(i)
  end do
  close(unit=unitnr)    
end subroutine readData



subroutine readData2(unitnr,ios,array, inputfile,dim1,dim2)
  implicit none
  !declaration of variables
  integer :: unitnr,dim1,dim2
  integer :: i,ios
  real(rkp) :: array(dim1,dim2)
  character(len=*) :: inputfile

  open(unitnr,file=inputfile, status='old', action='read', &
       iostat=ios)
  !error check
  if (ios /= 0) then
    return
  end if 
  do i=1,dim1
    read(unitnr,*) array(i,:)
  end do
  close(unit=unitnr)    
end subroutine readData2



!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains all the routines necessary to return the custom energies
!    and forces from the RKHS interpolation
!
!///////////////////////////////////////////////////////////////////////////////
!Reads in the weights for the RKHS interpolation. This is invoked at the startup of CHARMM


!transforms alpha to the new variable y=(1-cos(alpha))/2
real(rkp) Pure Function y(alpha) 
  implicit none
  real(rkp),intent(in) :: alpha
  y = (1d0-DCOS(alpha*PI/180d0))/2d0
end function y

real(rkp) Pure Function dy(alpha) !used for calculating the analytical derivative with respect to alpha (chainrule)
  implicit none
  real(rkp),intent(in) :: alpha
  dy = (PI/360d0)*DSIN(alpha*PI/180d0)
end function dy


!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains functions and subroutines needed for unit or coordinate
!    conversion or similar things. 
!
!///////////////////////////////////////////////////////////////////////////////

Pure Subroutine calc_Jacobi_coordinates(X,Y,Z,jacCoord,mass,dJacdCart)
!This subroutine calculates the Jacobi coordinates from the cartesian coordinates
!The indices of the Jacobi coordinate arrays coorespond to different surfaces:
!(the numbers here are the atom labels)
!index 1 (surface 1) : 1=2---3 
!index 2 (surface 2) : 1=3---2
!index 3 (surface 3) : 2=3---1
  implicit none
  
  !Input coordinates (cartesian)
  real(rkp), dimension(3), intent(in) :: X, Y, Z
  
  !Input masses (needed for some centre of mass)
  real(rkp), dimension(3), intent(in) :: mass
  
  !Output coordinates (Jacobi)
  !1st index: surface
  !2nd index: Jacobi coordinate (r, bigR, alpha)
  real(rkp), dimension(3,3), intent(out) :: jacCoord
  
  !Output derivatives of Jacobi coordinates with respect to cartesian coordinates
  !1st index: surface
  !2nd index: Jacobi coordinate (r, bigR, alpha)
  !3rd index: cartesian coordinate (x1,y1,z1,x2,y2,z2,x3,y3,z3)
  real(rkp), intent(out) :: dJacdCart(3,3,9)
  
  !temporary array to store Jacobi coordinates. We need this for swapping back!
  real(rkp) :: dJdC_tmp(3,9)
  
  !Utility variable for convenience that stores all the coordinates in one array
  real(rkp), dimension(9) :: coord
  
  !temporary mass array (gets reordered depending on which surface is calculated)
  real(rkp), dimension(3) :: mass_tmp
  
  !loop variables
  integer i
  
  
  !calculate the Jacobi coordinates for surface 1 (1=2- -3)
    !store coordinates appropriately
  do i=1,3
    coord(1+3*(i-1)) = real(X(i),8)
    coord(2+3*(i-1)) = real(Y(i),8)
    coord(3+3*(i-1)) = real(Z(i),8)
  end do
  call cartesian2jacobi(coord,jacCoord(1,:),dJacdCart(1,:,:),mass)
  
  !calculate the Jacobi coordinates for surface 2 (1=3- -2)
    !store coordinates appropriately
  coord(1) = real(X(1),8)
  coord(2) = real(Y(1),8)
  coord(3) = real(Z(1),8)
  coord(4) = real(X(3),8)
  coord(5) = real(Y(3),8)
  coord(6) = real(Z(3),8)
  coord(7) = real(X(2),8)
  coord(8) = real(Y(2),8)
  coord(9) = real(Z(2),8)
  !build a temporary mass array with swapped around masses
  mass_tmp(1) = mass(1)
  mass_tmp(2) = mass(3)
  mass_tmp(3) = mass(2)
  call cartesian2jacobi(coord,jacCoord(2,:),dJacdCart(2,:,:),mass_tmp)
  
  !calculate the Jacobi coordinates for surface 3 (2=3- -1)
    !store coordinates appropriately
  coord(1) = real(X(2),8)
  coord(2) = real(Y(2),8)
  coord(3) = real(Z(2),8)
  coord(4) = real(X(3),8)
  coord(5) = real(Y(3),8)
  coord(6) = real(Z(3),8)
  coord(7) = real(X(1),8)
  coord(8) = real(Y(1),8)
  coord(9) = real(Z(1),8)
  !build a temporary mass array with swapped around masses
  mass_tmp(1) = mass(2)
  mass_tmp(2) = mass(3)
  mass_tmp(3) = mass(1)
  call cartesian2jacobi(coord,jacCoord(3,:),dJacdCart(3,:,:),mass_tmp) 
  
  !Swap back coordinates for surface 2
  dJdC_tmp = dJacdCart(2,:,:)
  dJacdCart(2,:,4:6) = dJdC_tmp(:,7:9)
  dJacdCart(2,:,7:9) = dJdC_tmp(:,4:6)
  
  !Swap back coordinates for surface 3
  dJdC_tmp = dJacdCart(3,:,:)
  dJacdCart(3,:,1:3) = dJdC_tmp(:,7:9)
  dJacdCart(3,:,4:6) = dJdC_tmp(:,1:3)
  dJacdCart(3,:,7:9) = dJdC_tmp(:,4:6)

  return
end subroutine calc_Jacobi_coordinates


Pure Subroutine cartesian2jacobi(x,xi,dxidx,m)
  !This subroutine expects a 9-vector as input that contains the cartesian coordinates,        
  !x -> cartesian coordinates (x1,y1,z1,x2,y2,z2,x3,y3,z3)

  !a 3-vector as output that contains the Jacobi coordinates,
  !xi -> 3 Jacobi coordinates (r,R,alpha)

  !an array to store the derivatives of the Jacobi coordinates with respect to the
  !cartesian coordinates. The 1st index corresponds to the Jacobi coordinate (r,R,alpha),
  !while the 2nd index corresponds to the cartesian coordinate (x1,y1,z1,x2,y2,z2,x3,y3,z3).
  !dxidx -> derivatives of the Jacobi coordinates with respect to cartesians

  !Also, a 3-vector containing the masses of the atoms is required. This is because the Jacobi
  !coordinates reference the centre of mass between atom1 and atom2.
  !m -> masses of the 3 atoms
        
  implicit none
  real(rkp), parameter   :: PI = 3.141592653589793238462643383279502884197d0
  real(rkp), intent(in)  :: x(9)
  real(rkp), intent(in)  :: m(3)
  real(rkp), intent(out) :: xi(3)
  real(rkp), intent(out) :: dxidx(3,9)
  real(rkp)::  x12cm(3),dalphadx2(3),dalphadx12cm(3),dalphadx3(3)
           
  !distance of atom1 and atom2 (Jacobi coordinate r)
  xi(1)  = distance(x(1:3),x(4:6))
      
  !centre of mass of atom1 and atom2 (needed to calculate bigR)
  x12cm  = (m(1)*x(1:3)+m(2)*x(4:6))/(m(1)+m(2))
           
  !distance of atom3 and centre of mass of atom1 and atom2 (Jacobi coordinate bigR)
  xi(2)  = distance(x12cm,x(7:9))

  !Derivatives of the Jacobi coordinates with respect to the cartesian coordinates
  !dr/dx
  dxidx(1,1:3) =  (x(1:3)-x(4:6))/xi(1)
  dxidx(1,4:6) = -dxidx(1,1:3)
  dxidx(1,7:9) =  0.d0
  
  !dR/dx       
  dxidx(2,1:3) = (x12cm(1:3)-x(7:9))/xi(2)*m(1)/(m(1)+m(2)) 
  dxidx(2,4:6) = (x12cm(1:3)-x(7:9))/xi(2)*m(2)/(m(1)+m(2))   
  dxidx(2,7:9) = (x(7:9)-x12cm(1:3))/xi(2)     

  !calculate Jacobi coordinate alpha and derivatives needed to get the derivatives of
  !alpha with respect to cartesian coordinates (all in radians!!!)
  call angle(x(4:6),x12cm(1:3),x(7:9),xi(3),dalphadx2,dalphadx12cm,dalphadx3)
      
  !calculate the derivatives of alpha with respect to cartesian coordinates
  dxidx(3,1:3) = dalphadx12cm*m(1)/(m(1)+m(2))
  dxidx(3,4:6) = dalphadx12cm*m(2)/(m(1)+m(2)) + dalphadx2
  dxidx(3,7:9) = dalphadx3  
      
  !transform the derivatives of alpha and alpha itself from radians 
  !to degrees
  dxidx(3,1:9) = dxidx(3,1:9)*180.d0/PI
  xi(3)        = xi(3)*180.d0/PI
      
  return     
end subroutine cartesian2jacobi
    


Pure Subroutine angle(a,b,c,phi,dphida,dphidb,dphidc)
!This subroutine calculates the defined by the points a,b,c and the derivative of the angle
!with respect to these three points (all in radians)
  implicit none
  real(rkp),intent(in) ::a(3),b(3),c(3)
  real(rkp),intent(out)::phi
  real(rkp),intent(out),optional::dphida(3),dphidb(3),dphidc(3)
  real(rkp) :: y(3),z(3),y2,z2,ny,nz,yz,nynz,y0z0
  real(rkp) :: dcosphidy(3),dcosphidz(3),dcosdphi,dphidcos
  real(rkp) :: dphidy(3),dphidz(3)    
  ! y=a-b, z=c-b
  y   =  a-b
  z   =  c-b
  ! y2, z2
  y2  = sum(y**2) 
  z2  = sum(z**2) 
  ! abs(y), abs(z)
  ny  = dsqrt(y2)
  nz  = dsqrt(z2)
  ! yz
  yz  = sum(y*z)
  ! abs(y)abs(z)
  nynz =ny*nz   
  ! y0*z0

  y0z0=   yz/nynz
  
  y0z0=   min(1d0,max(-1d0,y0z0))         !=cosphi
  phi =   dacos(y0z0)
      
  if (.not.present(dphida).or..not.present(dphidb).or..not.present(dphidc)) return
      
  ! calculate derivatives
  ! y=a-b, z=c-b
  dcosphidy(1:3)=(z(1:3)-y(1:3)/y2*yz)/nynz 
  dcosphidz(1:3)=(y(1:3)-z(1:3)/z2*yz)/nynz 
  ! dcosdphi
  dcosdphi = -dsin(phi)                    
  ! dphidcos 

  !THIS WAS ADDED BY OLIVER TO CIRCUMVENT DIVISION BY ZERO!!! IMPORTANT!!!
  if(dabs(dcosdphi) < epsilon(0.d0)) dcosdphi = epsilon(0.d0)
  
  dphidcos =1d0/dcosdphi 
  !print*, "dphidcos", dphidcos
  
  ! dphidy dphidz
  dphidy(1:3) = dphidcos*dcosphidy
  dphidz(1:3) = dphidcos*dcosphidz   
  ! dphida= dphidy*dyda = dphidy
  dphida =  dphidy
  ! dphidc= dphidz*dzdc = dphidz
  dphidc =  dphidz
  ! dphida= dphidy*dydb + dphidz*dzdb = - dphidy - dphidz 
  dphidb = -dphida-dphidc
  
  return
end subroutine angle  
    

real(rkp) Pure function distance(a,b)
  implicit none 
  !real(rkp) ::distance
  real(rkp),intent(in) ::a(3),b(3)
  distance=dsqrt(sum((a-b)*(a-b)))
end function distance
  

!*****************************************************************************************
!    EVERYTHING BETWEEN THESE *** MARKERS WAS ADDED! EVERYTHING ELSE IN THIS SUBROUTINE IS
!    VANILLA CHARMM CODE! YOU SHOULD BE ABLE TO RETURN TO VANILLA CODE BY COMMENTING ALL
!    THE STUFF BETWEEN THESE MARKERS OUT. IMPORTANT: UNCOMMENT EU=0.0 LINE RIGHT AFTER
!    THE VARIABLE DECLARATIONS IF YOU WANT VANILLA CHARMM CODE (THIS WAS THERE ORIGINALLY)
!*****************************************************************************************
real(rkp) Pure function derivative_FivePointStencil(stencil,dx)
  implicit none
  real(rkp), dimension(4), intent(in) :: stencil
  real(rkp),               intent(in) :: dx
  derivative_FivePointStencil = (-stencil(1)+8*stencil(2)-8*stencil(3)+stencil(4))/(12*dx)
  return
end function derivative_FivePointStencil


End Module 