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

Module LEPS_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Kcm_To_Hartree, B_To_Ang
  use PES_Class             ,only:  PES_Type, DiaPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error
  use Input_Class           ,only:  Input_Type

  implicit none

  private
  public    :: LEPS_PES_Type 
  
  Type    ,extends(PES_Type)                :: LEPS_PES_Type
    real(rkp), dimension(3)                 :: re
    real(rkp), dimension(3)                 :: De
    real(rkp), dimension(3)                 :: DeA
    real(rkp), dimension(3)                 :: Alpha
    real(rkp), dimension(3)                 :: Beta
    real(rkp)                               :: S
    real(rkp)                               :: k
    real(rkp)                               :: A0
    real(rkp)                               :: A2
    real(rkp)                               :: A4
    character(10)                           :: model_output = 'LEPS'
    real(rkp)                               :: Rmin
    real(rkp)                               :: Vmin
    logical                                 :: ComputeDiatPotFlg = .False.
  contains
    procedure                               ::  Initialize        =>    Initialize_LEPS_PES
    procedure                               ::  Output            =>    Output_LEPS_PES
    procedure                               ::  Compute           =>    Compute_LEPS_PES_1d
    procedure                               ::  Potential         =>    LEPS_Potential_From_R
    procedure                               ::  TriatPotential    =>    LEPS_Potential_From_R_OnlyTriat
  End Type
  
  type(Input_Type)                          ::    Input
  
  logical                     ,parameter    :: i_Debug_Global = .False.

  contains


! **************************************************************************************************************
! **************************************************************************************************************
!                                       DEFERRED PROCEDURES for LEPS PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_LEPS_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use DiatomicPotential_Factory_Class     ,only:  DiatomicPotential_Factory_Type
  
  class(LEPS_PES_Type)                      ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms 
  integer                                   ,intent(in)     ::    iPES
  type(DiatomicPotential_Factory_Type)                       ::    DiaPotFactory
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP, jP
  character(*)                    ,parameter                ::    Name_PES = 'LEPS'
  
  character(:)                    ,allocatable              ::    LEPS_file
  integer                                                   ::    Status
  integer                                                   ::    Unit, idum
  
  character(:)                    ,allocatable              ::    LEPS_input_temp
  integer(8)                                                ::    iostat_SPES_input
  character(:)                    ,allocatable              ::    i_case
  character(150)                                            ::    line_input
  integer(8)                                                ::    i_char
  integer                                                   ::    ierr
  integer(8)                                                ::    N_diat_points = 1000
  integer(8)                                                ::    i
  integer         ,dimension(3,2)                           ::    iA
  real(rkp)                       ,parameter                ::    pi = 3.141592653589793_rkp
  character(1)    ,dimension(3)                             ::    AtomName
  integer         ,dimension(3)                             ::    ToIp
  integer                                                   ::    DiffAtomsFlg = 0
  ! integer         ,dimension(3,3)                           ::    iA_to_Pair

  logical                                                   ::    i_Debug_Loc

  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_LEPS_PES" )
  !i_Debug_Loc   =     Logger%On()


  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .False.
  This%NPairs       =   3               ! Setting the number of atom-atom pairs
  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair

  iA(1,:)           =   [1,2]
  iA(2,:)           =   [1,3]
  iA(3,:)           =   [2,3]
  ! iA_to_Pair(1,:)   =   [0,1,2]
  ! iA_to_Pair(2,:)   =   [1,0,3]
  ! iA_to_Pair(3,:)   =   [2,3,0]

  ! ==============================================================================================================
  !   CONSTRUCTING THE DIATOMIC POTENTIAL OBJECT
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "Constructing the diatomic potential object" )
  if (i_Debug_Loc) call Logger%Write( "-> Calling DiaPotFactory%Construct" )
  do iP = 1,This%NPairs
    call DiaPotFactory%Construct( Atoms, iA(iP,:), Input, This%Pairs(iP)%Vd, i_Debug=i_Debug_Loc )
  end do
  if (i_Debug_Loc) call Logger%Write( "-> Done constructing the diatomic potential" )
  ! ==============================================================================================================


  
  ! ==============================================================================================================
  !     READING LEPS INPUT FILE
  ! ==============================================================================================================
  if (i_Debug_Loc) call Logger%Write( "PES_ParamsFile = ", Input%PES_ParamsFile )
  if (adjustl(trim(Input%PES_ParamsFile)) == 'NONE') then
    LEPS_file = trim(adjustl(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) // '/PESs/LEPS/LEPS.dat'
  elseif (adjustl(trim(Input%PES_ParamsFile)) == 'Local') then
    LEPS_file = adjustl(trim(Input%OutputDir))  // '/' // trim(adjustl(Input%System)) //  '/LEPS.dat'    
  else 
    LEPS_file = adjustl(trim(Input%DtbPath))  // '/Systems/' // trim(adjustl(Input%System)) // '/PESs/LEPS/' // trim(adjustl(Input%PES_ParamsFile))
  end if
  if (i_Debug_Loc) call Logger%Write( "Reading the 'LEPS' Parameters file" )
  if (i_Debug_Loc) call Logger%Write( "-> Opening file: ", LEPS_file )
  open( File=LEPS_file, NewUnit=Unit, status='OLD', iostat=Status )
  if (Status/=0) call Error( "Error opening file: " // LEPS_file )
  
  do 
    
    read(Unit,'(A)',iostat=Status) line_input
         
    if (Status /= 0) then
      exit
    else
    
      if ( (line_input(1:1) == '#') .or. (line_input(1:10) == '          ') ) then
        continue
      else
      
        i_char = 1
        do
          if (line_input(i_char:i_char) == '=') exit
          i_char = i_char + 1
        end do
        i_case = adjustl(trim(line_input(1:(i_char-2))))

        select case (adjustl((TRIM(i_case))))
      
          case("Computing Diatomic Component?")
            line_input = line_input(i_char+2:150)
            if (trim(adjustl(line_input)) == 'yes') then
              This%ComputeDiatPotFlg = .True.
            end if 
            if (i_Debug_Loc) call Logger%Write( "Computing Diatomic Component?:      This%ComputeDiatPotFlg  = ", This%ComputeDiatPotFlg )

          case("Model for LEPS")
            This%model_output = line_input(i_char+2:150)
            if (i_Debug_Loc) call Logger%Write( "Model for LEPS:      This%model_output  = ", This%model_output )


          case("Atom(1)")
            line_input = line_input(i_char+2:150)
            AtomName(1) = (trim(adjustl(line_input)))
            DiffAtomsFlg = 1

          case("Atom(2)")
            line_input = line_input(i_char+2:150)
            AtomName(2) = (trim(adjustl(line_input)))
            DiffAtomsFlg = 1

          case("Atom(3)")
            line_input = line_input(i_char+2:150)
            AtomName(3) = (trim(adjustl(line_input)))
            DiffAtomsFlg = 1

          if (DiffAtomsFlg > 0) then
            do iP=1,3
              do jP=1,3
                if ( ( ( AtomName(iA(iP,1)) == Atoms(iA(jP,1))%Name ) .and. ( AtomName(iA(iP,2)) == Atoms(iA(jP,2))%Name ) ) .or. ( ( AtomName(iA(iP,1)) == Atoms(iA(jP,2))%Name ) .and. ( AtomName(iA(iP,2)) == Atoms(iA(jP,1))%Name ) ) ) then
                  ToIp(iP) = jP
                end if
              end do
            end do
            if (i_Debug_Loc) call Logger%Write( "Pair 1:      ToIp(1)  = ", ToIp(1) )
            if (i_Debug_Loc) call Logger%Write( "Pair 2:      ToIp(2)  = ", ToIp(2) )
            if (i_Debug_Loc) call Logger%Write( "Pair 3:      ToIp(3)  = ", ToIp(3) )
          else
            ToIp = [1,2,3]
          end if


          case("re(1)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%re(ToIp(1))
            
          case("re(2)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%re(ToIp(2))
            
          case("re(3)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%re(ToIp(3))

          case("re")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%re(1)
            this%re = this%re(1)
            
         
          case("De(1)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%De(ToIp(1))
            
          case("De(2)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%De(ToIp(2))
            
          case("De(3)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%De(ToIp(3))

          case("De")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%De(1)
            this%De = this%De(1)
            
            
          case("DeA(1)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%DeA(ToIp(1))
            
          case("DeA(2)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%DeA(ToIp(2))
            
          case("DeA(3)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%DeA(ToIp(3))

          case("DeA")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.4)') this%DeA(1)
            this%DeA = this%DeA(1)
            
         
          case("Alpha(1)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Alpha(ToIp(1)) 
            
          case("Alpha(2)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Alpha(ToIp(2)) 
            
          case("Alpha(3)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Alpha(ToIp(3)) 

          case("Alpha")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Alpha(1) 
            this%Alpha = this%Alpha(1) 
             
         
          case("Beta(1)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Beta(ToIp(1))
            
          case("Beta(2)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Beta(ToIp(2))
            
          case("Beta(3)")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Beta(ToIp(3))

          case("Beta")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%Beta(1)
            this%Beta = this%Beta(1)

         
          case("S")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%S
            
            
          case("A0")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%A0
            
          case("A2")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%A2
            
          case("A4")
            line_input = line_input(i_char+2:150)
            READ(line_input, '(d20.10)') this%A4
            
            
        end select
      
      end if
   
    end if
   
  end do

  close(Unit)
  
  this%S     = -this%S
  this%k     = this%S**2
  if (i_Debug_Loc) call Logger%Write( "LEPS re(1):       this%re(1)    = ", this%re(1) )
  if (i_Debug_Loc) call Logger%Write( "LEPS re(2):       this%re(2)    = ", this%re(2) )
  if (i_Debug_Loc) call Logger%Write( "LEPS re(3):       this%re(3)    = ", this%re(3) )
  if (i_Debug_Loc) call Logger%Write( "LEPS De(1):      this%De(1)   = ", this%De(1) )
  if (i_Debug_Loc) call Logger%Write( "LEPS De(2):      this%De(2)   = ", this%De(2) )
  if (i_Debug_Loc) call Logger%Write( "LEPS De(3):      this%De(3)   = ", this%De(3) )
  if (i_Debug_Loc) call Logger%Write( "LEPS DeA(1):      this%DeA(1)   = ", this%DeA(1) )
  if (i_Debug_Loc) call Logger%Write( "LEPS DeA(2):      this%DeA(2)   = ", this%DeA(2) )
  if (i_Debug_Loc) call Logger%Write( "LEPS DeA(3):      this%DeA(3)   = ", this%DeA(3) )
  if (i_Debug_Loc) call Logger%Write( "LEPS Alpha(1):    this%Alpha(1) = ", this%Alpha(1) )
  if (i_Debug_Loc) call Logger%Write( "LEPS Alpha(2):    this%Alpha(2) = ", this%Alpha(2) )
  if (i_Debug_Loc) call Logger%Write( "LEPS Alpha(3):    this%Alpha(3) = ", this%Alpha(3) ) 
  if (i_Debug_Loc) call Logger%Write( "LEPS Beta(1):     this%Beta(1)  = ", this%Beta(1) )    
  if (i_Debug_Loc) call Logger%Write( "LEPS Beta(2):     this%Beta(2)  = ", this%Beta(2) )
  if (i_Debug_Loc) call Logger%Write( "LEPS Beta(3):     this%Beta(3)  = ", this%Beta(3) )
  if (i_Debug_Loc) call Logger%Write( "LEPS S:           this%S        = ", this%S )
  if (i_Debug_Loc) call Logger%Write( "LEPS A0:          this%A0       = ", this%A0 )
  if (i_Debug_Loc) call Logger%Write( "LEPS A2:          this%A2       = ", this%A2 )
  if (i_Debug_Loc) call Logger%Write( "LEPS A4:          this%A4       = ", this%A4 )
  
  if (i_Debug_Loc) call Logger%Exiting
    
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_LEPS_PES( This, Unit )

  class(LEPS_PES_Type)                 ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function LEPS_Potential_From_R( This, R, Q ) result( V )
   
  class(LEPS_PES_Type)                          ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    EB, EA
  real(rkp) ,dimension(3)                                    ::    DedR, dEAdR
  real(rkp) ,dimension(3)                                    ::    QQ, A
  real(rkp) ,dimension(3)                                    ::    dQQdR, dAdR
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat
  real(rkp)                                                  ::    VTheta  
  integer                                                    ::    iP 
  
  !!! LEPS Distance are expressed in Angstrom. Therefore R NASA -> R LEPS : 1 -> B_To_Ang
  !!! LEPS Energies are expressed in Kcal/mol. Therefore V LEPS -> V NASA : 1 -> Kcm_To_Hartree

  !call computing_LEPS(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, R(:,iTraj) * B_To_Ang, V(iTraj))
  !call computing_LEPS_all(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, This%A0, This%A2, This%A4, R(:) * B_To_Ang, This%model_output,  V)

  call Compute_EB_EA( R(:)*B_To_Ang,  This%De,  This%DeA,  This%Beta,  This%Alpha,  This%re,  This%k, EB, EA, DedR, dEAdR )
  call compute_Q_A( EB, EA, DedR, dEAdR, QQ, A, dQQdR, dAdR )

  if (.not. This%ComputeDiatPotFlg) then 
    do iP=1,3
      call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
    end do
  else 
    call ComputeDiatomic( QQ, dQQdR, A, dAdR, VDiat, dVDiat )
    VDiat  = VDiat  * Kcm_To_Hartree
  end if

  VTheta = Zero
  !if ( trim(This%model_output) == trim('LEPS-MOD') ) then
  !  call computing_LEPS_Vtheta(This%De, This%Alpha, This%re, R(:) * B_To_Ang, This%A0, This%A2, This%A4, VTheta)    
  !  VTheta = VTheta * Kcm_To_Hartree
  !end if

  call ComputeTriatomic( A, dAdR, VTriat, dVTriat )
  VTriat = VTriat  * Kcm_To_Hartree
  !call ComputeAll( QQ, dQQdR, A, dAdR, V, dVdR )
  !V    = V    * Kcm_To_Hartree

  V = sum(VDiat) + VTriat + VTheta

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Pure Function LEPS_Potential_From_R_OnlyTriat( This, R, Q ) result( V )
   
  class(LEPS_PES_Type)                          ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(3)                                    ::    EB, EA
  real(rkp) ,dimension(3)                                    ::    DedR, dEAdR
  real(rkp) ,dimension(3)                                    ::    QQ, A
  real(rkp) ,dimension(3)                                    ::    dQQdR, dAdR
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat
  real(rkp)                                                  ::    VTheta  
  integer                                                    ::    iP 
  
  !!! LEPS Distance are expressed in Angstrom. Therefore R NASA -> R LEPS : 1 -> B_To_Ang
  !!! LEPS Energies are expressed in Kcal/mol. Therefore V LEPS -> V NASA : 1 -> Kcm_To_Hartree

  !call computing_LEPS(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, R(:,iTraj) * B_To_Ang, V(iTraj))
  !call computing_LEPS_all(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, This%A0, This%A2, This%A4, R(:) * B_To_Ang, This%model_output,  V)
  

  call Compute_EB_EA( R(:)*B_To_Ang,  This%De,  This%DeA,  This%Beta,  This%Alpha,  This%re,  This%k, EB, EA, DedR, dEAdR )
  call compute_Q_A( EB, EA, DedR, dEAdR, QQ, A, dQQdR, dAdR )

  VTheta = Zero
  !if ( trim(This%model_output) == trim('LEPS-MOD') ) then
  !  call computing_LEPS_Vtheta(This%De, This%Alpha, This%re, R(:) * B_To_Ang, This%A0, This%A2, This%A4, VTheta)    
  !  VTheta = VTheta * Kcm_To_Hartree
  !end if

  call ComputeTriatomic( A, dAdR, VTriat, dVTriat )
  VTriat = VTriat * Kcm_To_Hartree
  !call ComputeAll( QQ, dQQdR, A, dAdR, V, dVdR )
  !V    = V    * Kcm_To_Hartree

  V = VTriat + VTheta

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_LEPS_PES_1d( This, R, Q, V, dVdR, dVdQ )

  class(LEPS_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  real(rkp) ,dimension(3)                                    ::    EB, EA
  real(rkp) ,dimension(3)                                    ::    DedR, dEAdR
  real(rkp) ,dimension(3)                                    ::    QQ, A
  real(rkp) ,dimension(3)                                    ::    dQQdR, dAdR
  real(rkp) ,dimension(3)                                    ::    VDiat
  real(rkp)                                                  ::    VTriat
  real(rkp) ,dimension(3)                                    ::    dVDiat, dVTriat
  real(rkp)                                                  ::    VTheta  
  integer                                                    ::    iP       
  
  !!! LEPS Distance are expressed in Angstrom. Therefore R NASA -> R LEPS : 1 -> B_To_Ang
  !!! LEPS Energies are expressed in Kcal/mol. Therefore V LEPS -> V NASA : 1 -> Kcm_To_Hartree
  !!! LEPS Forces are expressed in (Kcal/mol) / Angstrom. Therefore dV LEPS -> dV NASA : 1 -> Kcm_To_Hartree * B_To_Ang

  !call computing_LEPS(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, R(:) * B_To_Ang, V)
  
  !call computing_LEPS_all(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, This%A0, This%A2, This%A4, R(:) * B_To_Ang, This%model_output,  V)                  
  !call computing_LEPS_gradient(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, This%A0, This%A2, This%A4, R(:) * B_To_Ang, This%model_output, dVdR(:))    
  ! call computing_LEPS_num_gradient(This%S, This%De, This%DeA, This%Alpha, This%Beta, This%re, This%A0, This%A2, This%A4, R(:) * B_To_Ang, This%model_output, 1.d-8, int(8,8), dVdR(:))

  !V    = V    * Kcm_To_Hartree
  !dVdR = dVdR * Kcm_To_Hartree * B_To_Ang

  !rite(*,*) This%De,  This%DeA,  This%Beta,  This%Alpha,  This%re,  This%k
  !pause    
  call Compute_EB_EA( R(:)*B_To_Ang,  This%De,  This%DeA,  This%Beta,  This%Alpha,  This%re,  This%k, EB, EA, DedR, dEAdR )
  call compute_Q_A( EB, EA, DedR, dEAdR, QQ, A, dQQdR, dAdR )

  if (.not. This%ComputeDiatPotFlg) then 
    do iP=1,3
      call This%Pairs(iP)%Vd%Compute_Vd_dVd( R(iP), VDiat(iP), dVDiat(iP) )
    end do
  else 
    call ComputeDiatomic( QQ, dQQdR, A, dAdR, VDiat, dVDiat )
    VDiat  = VDiat  * Kcm_To_Hartree
    dVDiat = dVDiat * Kcm_To_Hartree * B_To_Ang
  end if

  VTheta = Zero
  !if ( trim(This%model_output) == trim('LEPS-MOD') ) then
  !  call computing_LEPS_Vtheta(This%De, This%Alpha, This%re, R(:) * B_To_Ang, This%A0, This%A2, This%A4, VTheta)    
  !  VTheta = VTheta * Kcm_To_Hartree    
  !end if

  call ComputeTriatomic( A, dAdR, VTriat, dVTriat )
  VTriat  = VTriat  * Kcm_To_Hartree
  dVTriat = dVTriat * Kcm_To_Hartree * B_To_Ang
  !call ComputeAll( QQ, dQQdR, A, dAdR, V, dVdR )
  !V    = V    * Kcm_To_Hartree
  !dVdR = dVdR * Kcm_To_Hartree * B_To_Ang

  V    = sum(VDiat) + VTriat + VTheta
  dVdR = dVDiat     + dVTriat

  dVdQ = Zero
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
! **************************************************************************************************************
!                                         PRIVATE PROCEDURES for LEPS
! **************************************************************************************************************
!________________________________________________________________________________________________________________________________!     
Pure Subroutine compute_EB_EA( r, De, DeA, Beta, Alpha, re, k, EB, EA, DedR, dEAdR )

  real(rkp)     ,dimension(3) ,intent(in)  :: r  
  real(rkp)     ,dimension(3) ,intent(in)  :: De  
  real(rkp)     ,dimension(3) ,intent(in)  :: DeA
  real(rkp)     ,dimension(3) ,intent(in)  :: Beta  
  real(rkp)     ,dimension(3) ,intent(in)  :: Alpha
  real(rkp)     ,dimension(3) ,intent(in)  :: re
  real(rkp)                   ,intent(in)  :: k
  real(rkp)     ,dimension(3) ,intent(out) :: EB
  real(rkp)     ,dimension(3) ,intent(out) :: EA
  real(rkp)     ,dimension(3) ,intent(out) :: DedR
  real(rkp)     ,dimension(3) ,intent(out) :: dEAdR

  real(rkp)     ,dimension(3)              :: XpntA, XpntB

  XpntB = - Beta  * (r - re)
  EB    = De/Two * ( exp(Two*XpntB) - Two*exp(XpntB) ) 
  DedR = - Two * De/Two * Beta * exp(XpntB) * (exp(XpntB) - One)

  XpntA = - Alpha * (r - re)
  EA    = DeA/Two * ( exp(Two*XpntA) + Two*exp(XpntA) ) * (One-k)/(One+k) 
  dEAdR = - Two * DeA/Two * Alpha * exp(XpntA) * (exp(XpntA) + One) * (One-k)/(One+k) 

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!     
Pure Subroutine compute_Q_A( EB, EA, DedR, dEAdR, Q, A, dQdR, dAdR )
 
  real(rkp)     ,dimension(3) ,intent(in)  :: EA
  real(rkp)     ,dimension(3) ,intent(in)  :: EB
  real(rkp)     ,dimension(3) ,intent(in)  :: DedR
  real(rkp)     ,dimension(3) ,intent(in)  :: dEAdR
  real(rkp)     ,dimension(3) ,intent(out) :: Q
  real(rkp)     ,dimension(3) ,intent(out) :: A
  real(rkp)     ,dimension(3) ,intent(out) :: dQdR
  real(rkp)     ,dimension(3) ,intent(out) :: dAdR

  Q    = EB + EA 
  dQdR = DedR + dEAdR

  A    = EB - EA 
  dAdR = DedR - dEAdR

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!     
Pure Subroutine ComputeDiatomic( Q, dQdR, A, dAdR, V, dVdR )
 
  real(rkp)     ,dimension(3) ,intent(in)  :: Q
  real(rkp)     ,dimension(3) ,intent(in)  :: dQdR 
  real(rkp)     ,dimension(3) ,intent(in)  :: A
  real(rkp)     ,dimension(3) ,intent(in)  :: dAdR
  real(rkp)     ,dimension(3) ,intent(out) :: V
  real(rkp)     ,dimension(3) ,intent(out) :: dVdR

  V    = Q + A
  dVdR = dQdR + dAdR

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!     
Pure Subroutine ComputeTriatomic( A, dAdR, V, dVdR )
 
  real(rkp)     ,dimension(3) ,intent(in)  :: A
  real(rkp)     ,dimension(3) ,intent(in)  :: dAdR
  real(rkp)                   ,intent(out) :: V
  real(rkp)     ,dimension(3) ,intent(out) :: dVdR

  real(rkp)                                :: Temp
  real(rkp)     ,dimension(3)              :: dVdA

  Temp = ( sum(A**2) - A(1)*A(2) - A(3)*A(2) - A(1)*A(3) )**Half

  V    = - (Temp + sum(A))
  dVdA = (sum(A) -  (Three * A) ) /(Two*Temp) - One

  dVdR = dVdA * dAdR

End Subroutine
!-------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!     
Pure Subroutine ComputeAll( Q, dQdR, A, dAdR, V, dVdR )

  real(rkp)     ,dimension(3) ,intent(in)  :: Q
  real(rkp)     ,dimension(3) ,intent(in)  :: dQdR 
  real(rkp)     ,dimension(3) ,intent(in)  :: A
  real(rkp)     ,dimension(3) ,intent(in)  :: dAdR
  real(rkp)                   ,intent(out) :: V
  real(rkp)     ,dimension(3) ,intent(out) :: dVdR

  real(rkp)                                :: Temp
  real(rkp)     ,dimension(3)              :: dVdA

  V    =   sum(Q) - ( sum(A**2) - A(1)*A(2) - A(3)*A(2) - A(1)*A(3) )**Half
  dVdA = - sum(A)/(Two*Temp) - One + (Three * A)

  dVdR = dVdA * dAdR

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!    
Pure Subroutine computing_LEPS_W( De, Alpha, re, r, W )

   real(rkp) ,dimension(3) ,intent(in)  :: De
   real(rkp) ,dimension(3) ,intent(in)  :: Alpha
   real(rkp) ,dimension(3) ,intent(in)  :: re
   real(rkp) ,dimension(3) ,intent(in)  :: r
   
   real(rkp)               ,intent(out) :: W
   
   real(rkp) ,dimension(3)              :: Vrep
   integer                              :: i, j
   
   W = Zero
   
   Vrep = De*dexp(-Two*Alpha*(r-re))
   
   do i=1,3
   
      if (abs(Vrep(i))<1.0e-100_rkp) Vrep(i)=Zero
   
   end do
   
   do i=1,2
   
      do j=i+1,3
         
         W = W + Vrep(i)*Vrep(j)
      
      end do
   
   end do
          
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!     
Pure Subroutine computing_LEPS_Vtheta( De, Alpha, re, r, A0, A2, A4, Vtheta )
      
   real(rkp) ,dimension(3) ,intent(in)  :: De
   real(rkp) ,dimension(3) ,intent(in)  :: Alpha
   real(rkp) ,dimension(3) ,intent(in)  :: re
   real(rkp) ,dimension(3) ,intent(in)  :: r
   real(rkp)               ,intent(in)  :: A0
   real(rkp)               ,intent(in)  :: A2
   real(rkp)               ,intent(in)  :: A4
   
   real(rkp)               ,intent(out) :: Vtheta
   
   real(rkp)                            :: W
   real(rkp)                            :: angle
   
   call computing_LEPS_W(De, Alpha, re, r, W)
   
   angle = (r(1)**2 + r(2)**2 - r(3)**2) / (Two*r(1)*r(2))

   Vtheta = (A0 + A2*dcos(angle)**2 + A4*dcos(angle)**4) * Half * dtanh(0.4_rkp*(dsqrt(W) - Seven) + One)
          
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! **************************************************************************************************************
!Pure Subroutine computing_LEPS( S, De, DeA, Alpha, Beta, re, r, V_LEPS )
!      
!   real(rkp)               ,intent(in)  :: S
!   real(rkp) ,dimension(3) ,intent(in)  :: De
!   real(rkp) ,dimension(3) ,intent(in)  :: DeA
!   real(rkp) ,dimension(3) ,intent(in)  :: Alpha
!   real(rkp) ,dimension(3) ,intent(in)  :: Beta
!   real(rkp) ,dimension(3) ,intent(in)  :: re
!   real(rkp) ,dimension(3) ,intent(in)  :: r
!   
!   real(rkp)               ,intent(out) :: V_LEPS
!   
!   real(rkp) ,dimension(6)              :: E
!   real(rkp) ,dimension(6)              :: I
!   
!   call computing_singlet(S, De, Alpha, re, r, E(1:3))
!   
!   call computing_triplet(S, DeA, Beta, re, r, E(4:6))
!   
!   call computing_LEPS_integrals(E, I)

!   call computing_LEPS_potential(I, S, V_LEPS)
!        
!End Subroutine
!!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________!
!Pure Subroutine computing_singlet( S, De, Alpha, re, r, E )
!      
!   real(rkp)                ,intent(in)  :: S
!   real(rkp) ,dimension(3)  ,intent(in)  :: De
!   real(rkp) ,dimension(3)  ,intent(in)  :: Alpha
!   real(rkp) ,dimension(3)  ,intent(in)  :: re
!   real(rkp) ,dimension(3)  ,intent(in)  :: r
!   
!   real(rkp) ,dimension(3)  ,intent(out) :: E
!   
!   E=Zero

!   E = (One+S)*De*(dexp(-Two*Alpha*(r-re)) - Two*dexp(-Alpha*(r-re)))
!         
!End Subroutine
!!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________! 
!Pure Subroutine computing_triplet( S, DeA, Beta, re, r, E )
!      
!   real(rkp),                 intent(in)  :: S
!   real(rkp), dimension(3),   intent(in)  :: DeA
!   real(rkp), dimension(3),   intent(in)  :: Beta
!   real(rkp), dimension(3),   intent(in)  :: re
!   real(rkp), dimension(3),   intent(in)  :: r
!   
!   real(rkp), dimension(3),   intent(out) :: E
!   
!   E=Zero
!   
!   E = (One-S)*DeA*(dexp(-Two*Beta*(r-re)) + Two*dexp(-Beta*(r-re)))
!          
!End Subroutine 
!!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________!
!Pure Subroutine computing_LEPS_integrals(E, I)
!      
!   real(rkp) ,dimension(6)   ,intent(in)  :: E
!   
!   real(rkp) ,dimension(6)   ,intent(out) :: I

!   real(rkp) ,dimension(6)                :: EE
!   real(rkp) ,dimension(6,6)              :: A, AA
!   integer                                :: INFO
!   integer   ,dimension(6)                :: IPIV

!   real(rkp) ,dimension(6)   ,parameter   :: A1 = (/  One,  Zero,  Zero,   One,  Zero,  Zero /)
!   real(rkp) ,dimension(6)   ,parameter   :: A2 = (/ Zero,   One,  Zero,  Zero,   One,  Zero /)
!   real(rkp) ,dimension(6)   ,parameter   :: A3 = (/ Zero,  Zero,   One,  Zero,  Zero,   One /)
!   real(rkp) ,dimension(6)   ,parameter   :: A4 = (/  One,  Zero,  Zero,  -One,  Zero,  Zero /)
!   real(rkp) ,dimension(6)   ,parameter   :: A5 = (/ Zero,   One,  Zero,  Zero,  -One,  Zero /)
!   real(rkp) ,dimension(6)   ,parameter   :: A6 = (/ Zero,  Zero,   One,  Zero,  Zero,  -One /)
!   
!   A(1,:) = A1
!   A(2,:) = A2
!   A(3,:) = A3
!   A(4,:) = A4
!   A(5,:) = A5
!   A(6,:) = A6
!   
!   AA = A
!   EE = E
!   I  = Zero 

!!   call dgesv(6, 1, AA, 6, IPIV, EE, 6, INFO)
!   
!   if (INFO == 0) then
!   
!      I = EE
!      
!   end if
!    
!End Subroutine
!!--------------------------------------------------------------------------------------------------------------------------------!


!!________________________________________________________________________________________________________________________________!     
!Pure Subroutine computing_LEPS_potential( I, S, V_LEPS )
!      
!      real(rkp)               ,intent(in)  :: S
!      real(rkp) ,dimension(6) ,intent(in)  :: I
!      
!      real(rkp)               ,intent(out) :: V_LEPS

!      V_LEPS = (sum(I(1:3)) - dsqrt(Half*( (I(4)-I(5))**2 + (I(4)-I(6))**2 + (I(5)-I(6))**2 )))/(One+S)
!         
!End Subroutine
!!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!    
Pure Subroutine computing_LEPS_gradient( S, De, DeA, Alpha, Beta, re, A0, A2, A4, r, model_output, dv )
      
   real(rkp)                    ,intent(in)  :: S
   real(rkp)                    ,intent(in)  :: A0
   real(rkp)                    ,intent(in)  :: A2
   real(rkp)                    ,intent(in)  :: A4
   real(rkp)     ,dimension(3)  ,intent(in)  :: De
   real(rkp)     ,dimension(3)  ,intent(in)  :: DeA
   real(rkp)     ,dimension(3)  ,intent(in)  :: Alpha
   real(rkp)     ,dimension(3)  ,intent(in)  :: Beta
   real(rkp)     ,dimension(3)  ,intent(in)  :: re
   real(rkp)     ,dimension(3)  ,intent(in)  :: r
   character(10)                ,intent(in)  :: model_output
   
   real(rkp)     ,dimension(3)  ,intent(out) :: dv
   
   dv = Zero

   dv(1) = (One/(Two*(One+S)))*                                                                                         &                                          
     (De(1)*(One+S)*(-Two*Alpha(1)*dexp(-Two*Alpha(1)*(r(1)-re(1))) + Two*Alpha(1)*dexp(-Alpha(1)*(r(1)-re(1))))   + &
      DeA(1)*(One-S)*(-Two*Beta(1) *dexp(-Two*Beta(1) *(r(1)-re(1))) - Two*Beta(1) *dexp(-Beta(1) *(r(1)-re(1)))))  - &          
     ((0.353553_rkp*((Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))      - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + &
                   Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2))))))     * & 
      (De(1)*(One+S)*(-Two*Alpha(1)*dexp(-Two*Alpha(1)*(r(1)-re(1)))+Two*Alpha(1)*dexp(-Alpha(1)*(r(1)-re(1))))    - &
       DeA(1)*(One-S)*(-Two*Beta(1) *dexp(-Two*Beta(1) *(r(1)-re(1)))-Two*Beta(1) *dexp(-Beta(1) *(r(1)-re(1)))))   + &
                   (Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))        - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + & 
                   Half*(-De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))        + & 
                           DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))     * &
      (De(1)*(One+S)*(-Two*Alpha(1)*dexp(-Two*Alpha(1)*(r(1)-re(1)))+Two*Alpha(1)*dexp(-Alpha(1)*(r(1)-re(1))))    - &
       DeA(1)*(One-S)*(-Two*Beta(1) *dexp(-Two*Beta(1) *(r(1)-re(1)))-Two*Beta(1) *dexp(-Beta(1) *(r(1)-re(1))))))) / &
                  ((Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))        - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + &
                   Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2))))))**2  + &
                  (Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2)))))      + &
                    Half*(De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))        - &
                           DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))**2  + &
                   (Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))        - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + &
                   Half*(-De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))        + &
            DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))**2)**Half) / (One+S)

   dv(2) = (One/(Two*(One+S)))*                                                                                         &
     (De(2)*(One+S)*(-Two*Alpha(2)*dexp(-Two*Alpha(2)*(r(2)-re(2))) + Two*Alpha(2)*dexp(-Alpha(2)*(r(2)-re(2))))   + &
      DeA(2)*(One-S)*(-Two*Beta(2) *dexp(-Two*Beta(2) *(r(2)-re(2))) - Two*Beta(2) *dexp(-Beta(2) *(r(2)-re(2)))))  - &
     ((0.353553_rkp*((Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))      - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + &
                   Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2))))))     * &
   (-De(2)*(One+S)*(-Two*Alpha(2)*dexp(-Two*Alpha(2)*(r(2)-re(2)))+Two*Alpha(2)*dexp(-Alpha(2)*(r(2)-re(2))))      + &
     DeA(2)*(One-S)*(-Two*Beta(2)* dexp(-Two*Beta(2) *(r(2)-re(2)))-Two*Beta(2) *dexp(-Beta(2) *(r(2)-re(2)))))     + &
                  (Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2)))))      + &
                    Half*(De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))        - &
                           DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))     * &
    (-De(2)*(One+S)*(-Two*Alpha(2)*dexp(-Two*Alpha(2)*(r(2)-re(2)))+Two*Alpha(2)*dexp(-Alpha(2)*(r(2)-re(2))))     + &
      DeA(2)*(One-S)*(-Two*Beta(2) *dexp(-Two*Beta(2) *(r(2)-re(2)))-Two*Beta(2) *dexp(-Beta(2) *(r(2)-re(2)))))))  / &
                  ((Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))        - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + &
                   Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2))))))**2  + &
                  (Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))        + &
                           DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2)))))      + &
                    Half*(De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))        - &
                           DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))**2  + &
                   (Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))        - &
                           DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))      + &
                   Half*(-De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))        + &
            DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))**2)**Half) / (One+S)

   dv(3) = (One/(Two*(One+S)))*                                                                                         &
     (De(3)*(One+S)*(-Two*Alpha(3)*dexp(-Two*Alpha(3)*(r(3)-re(3))) + Two*Alpha(3)*dexp(-Alpha(3)*(r(3)-re(3))))   + &
      DeA(3)*(One-S)*(-Two*Beta(3) *dexp(-Two*Beta(3) *(r(3)-re(3))) - Two*Beta(3) *dexp(-Beta(3) *(r(3)-re(3)))))  - &   
     ((0.353553_rkp*((Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))     + &
                            DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2)))))     + &
                     Half*(De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))       - &
                            DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))    * &
      (De(3)*(One+S)*(-Two*Alpha(3)*dexp(-Two*Alpha(3)*(r(3)-re(3)))+Two*Alpha(3)*dexp(-Alpha(3)*(r(3)-re(3))))    - &
      DeA(3)*(One-S)*(-Two*Beta(3) *dexp(-Two*Beta(3) *(r(3)-re(3)))-Two*Beta(3) *dexp(-Beta(3) *(r(3)-re(3)))))    + & 
                    (Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))       - & 
                            DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))     + &
                    Half*(-De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))       + & 
                            DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))    * &
    (-De(3)*(One+S)*(-Two*Alpha(3)*dexp(-Two*Alpha(3)*(r(3)-re(3)))+Two*Alpha(3)*dexp(-Alpha(3)*(r(3)-re(3))))     + &
      DeA(3)*(One-S)*(-Two*Beta(3) *dexp(-Two*Beta(3) *(r(3)-re(3)))-Two*Beta(3) *dexp(-Beta(3) *(r(3)-re(3)))))))  / & 
                   ((Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))       - & 
                            DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))     + & 
                    Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))       + &
                            DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2))))))**2 + & 
                   (Half*(-De(2)*(One+S)*(dexp(-Two*Alpha(2)*(r(2)-re(2)))-Two*dexp(-Alpha(2)*(r(2)-re(2))))       + & 
                            DeA(2)*(One-S)*(dexp(-Two*Beta(2) *(r(2)-re(2)))+Two*dexp(-Beta(2) *(r(2)-re(2)))))     + &
                     Half*(De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))       - &
                            DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))**2 + & 
                    (Half*(De(1)*(One+S)*(dexp(-Two*Alpha(1)*(r(1)-re(1)))-Two*dexp(-Alpha(1)*(r(1)-re(1))))       - & 
                            DeA(1)*(One-S)*(dexp(-Two*Beta(1) *(r(1)-re(1)))+Two*dexp(-Beta(1) *(r(1)-re(1)))))     + & 
                    Half*(-De(3)*(One+S)*(dexp(-Two*Alpha(3)*(r(3)-re(3)))-Two*dexp(-Alpha(3)*(r(3)-re(3))))       + & 
            DeA(3)*(One-S)*(dexp(-Two*Beta(3) *(r(3)-re(3)))+Two*dexp(-Beta(3) *(r(3)-re(3))))))**2)**Half) / (One+S)

   if ( trim(model_output) == trim('LEPS-MOD') ) then  
   
   dv(1) = dv(1) + Half*((A2*(r(1)**2+r(2)**2-r(3)**2))/(r(1)*r(2)**2)-(A2*(r(1)**2+r(2)**2-r(3)**2)**2)/(Two*r(1)**3*r(2)**2)+& 
                    (A4*(r(1)**2+r(2)**2-r(3)**2)**3)/(Two*r(1)**3*r(2)**4)                                           - &
                    (A4*(r(1)**2+r(2)**2-r(3)**2)**4)/(Four*r(1)**5*r(2)**4))*dtanh(One+0.4_rkp*(-Seven               + &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                       De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half))     + &
                   (0.1_rkp*(A0+(A2*(r(1)**2+r(2)**2-r(3)**2)**2)/(Four*r(1)**2*r(2)**2)+(A4*(r(1)**2+r(2)**2-r(3)**2)**4) / &
                   (16.0_rkp*r(1)**4*r(2)**4))*dcosh(min(One+0.4_rkp*(-Seven                                          + &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + & 
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + & 
                 De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half),3.0e2_rkp))**(-2)*&
       (-Two*Alpha(1)*De(1)*De(2)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(1)*(r(1)-re(1)))               - & 
         Two*Alpha(1)*De(1)*De(3)*dexp(-Two*Alpha(3)*(r(3)-re(3)))*dexp(-Two*Alpha(1)*(r(1)-re(1)))))             / & 
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + & 
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + & 
                       De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half

   dv(2) = dv(2) + Half*((A2*(r(1)**2+r(2)**2-r(3)**2))/(r(1)**2*r(2))-(A2*(r(1)**2+r(2)**2-r(3)**2)**2)/(Two*r(1)**2*r(2)**3)+&
                    (A4*(r(1)**2+r(2)**2-r(3)**2)**3)/(Two*r(1)**4*r(2)**3)                                           - &
                    (A4*(r(1)**2+r(2)**2-r(3)**2)**4)/(Four*r(1)**4*r(2)**5))*dtanh(One+0.4_rkp*(-Seven               + &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                       De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half))     + &
                   (0.1_rkp*(A0+(A2*(r(1)**2+r(2)**2-r(3)**2)**2)/(Four*r(1)**2*r(2)**2)+(A4*(r(1)**2+r(2)**2-r(3)**2)**4) / &
                   (16.0_rkp*r(1)**4*r(2)**4))*dcosh(min(One+0.4_rkp*(-Seven                                          + &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                 De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half),3.0e2_rkp))**(-2)*&
       (-Two*Alpha(2)*De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))               - &
         Two*Alpha(2)*De(2)*De(3)*dexp(-Two*Alpha(3)*(r(3)-re(3)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))))             / &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                       De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half

   dv(3) = dv(3) + Half*(-((A2*r(3)*(r(1)**2+r(2)**2-r(3)**2))/(r(1)**2*r(2)**2))-(A4*r(3)*(r(1)**2+r(2)**2-r(3)**2)**3)/ &
                       (Two*r(1)**4*r(2)**4))*dtanh(One+0.4_rkp*(-Seven                                               + &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                       De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half))     + &
                       (0.1_rkp*(A0+(A2*(r(1)**2+r(2)**2-r(3)**2)**2)/(Four*r(1)**2*r(2)**2)                          + &
                      (A4*(r(1)**2+r(2)**2-r(3)**2)**4)/(16.0_rkp*r(1)**4*r(2)**4))*dcosh(min(One+0.4_rkp*(-Seven     + &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                 De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half),Three))**(-2)*&
       (-Two*Alpha(3)*De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))               - &
         Two*Alpha(3)*De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))))             / &
                      (De(1)*De(2)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(2)*(r(2)-re(2)))              + &
                       De(1)*De(3)*dexp(-Two*Alpha(1)*(r(1)-re(1)))*dexp(-Two*Alpha(3)*(r(3)-re(3)))              + &
                       De(2)*De(3)*dexp(-Two*Alpha(2)*(r(2)-re(2)))*dexp(-Two*Alpha(3)*(r(3)-re(3))))**Half
                                      
   end if

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! !________________________________________________________________________________________________________________________________!
! Subroutine computing_LEPS_num_gradient( S, De, DeA, Alpha, Beta, re, A0, A2, A4, r, model_output, h, order, dv )
      
!    real(rkp)                     ,intent(in)  :: S
!    real(rkp)                     ,intent(in)  :: A0
!    real(rkp)                     ,intent(in)  :: A2
!    real(rkp)                     ,intent(in)  :: A4
!    real(rkp)     ,dimension(3)   ,intent(in)  :: De
!    real(rkp)     ,dimension(3)   ,intent(in)  :: DeA
!    real(rkp)     ,dimension(3)   ,intent(in)  :: Alpha
!    real(rkp)     ,dimension(3)   ,intent(in)  :: Beta
!    real(rkp)     ,dimension(3)   ,intent(in)  :: re
!    real(rkp)     ,dimension(3)   ,intent(in)  :: r
!    character(10)                 ,intent(in)  :: model_output
!    real(rkp)                     ,intent(in)  :: h
!    integer                       ,intent(in)  :: order
   
!    real(rkp)     ,dimension(3)   ,intent(out) :: dv
   
!    real(rkp)     ,dimension(3,order,3)        :: r_temp
!    real(rkp)     ,dimension(order,3)          :: v 
!    integer                                    :: i, j, k
!    real(rkp)     ,dimension(order)            :: coeff
!    integer       ,dimension(order)            :: grid
    
!    if (order == 2) then
   
!       coeff = (/ -Half, Half /)
!       grid  = (/ -1, 1 /)
      
!    elseif (order == 4) then
   
!       coeff = (/ One/12.0_rkp, -Two/Three, Two/Three, -One/12.0_rkp /)
!       grid  = (/ -2, -1, 1, 2 /)
   
!    elseif (order == 6) then
   
!       coeff = (/ -One/60.0_rkp, Three/20.0_rkp, -Three/Four, Three/Four, -Three/20.0_rkp, One/60.0_rkp /)
!       grid  = (/ -3, -2, -1, 1, 2, 3 /)
   
!    elseif (order == 8) then
   
!       coeff = (/ One/280.0_rkp, -Four/105.0_rkp, One/5.0_rkp, -Four/Five, Four/Five, -One/Five, Four/105.0_rkp, -One/280.0_rkp /)
!       grid  = (/ -4, -3, -2, -1, 1, 2, 3, 4 /)
   
!    end if
   
!    do i=1,3
   
!       do k=1,order
      
!          r_temp(:,k,i) = r
         
!          r_temp(i,:,i) = r_temp(i,:,i) + h*grid
         
!          call computing_LEPS_all(S, De, DeA, Alpha, Beta, re, A0, A2, A4, r_temp(:,k,i), model_output, v(k,i))
      
!       end do
      
!       dv(i) = sum(v(:,i) * coeff(:)/h)
      
!    end do
         
! End Subroutine
! !--------------------------------------------------------------------------------------------------------------------------------!


End Module
