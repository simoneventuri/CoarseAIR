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

Module N4_NASA_PES_Class

#include "../qct.inc"

  use Parameters_Module     ,only:  rkp, Zero, Half, One, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten, Pi
  use PES_Class             ,only:  PES_Type, DiatPotContainer_Type
  use Logger_Class          ,only:  Logger
  use Error_Class           ,only:  Error

  implicit none

  ! Pre-processor conditional compilation to suppress verbosity
!#undef N4_PES_VERBOSE   
!!#define N4_EXP_FIX

  private
  public    ::    N4_NASA_PES_Type


  Type, extends(PES_Type) :: N4_NASA_PES_Type
    real(rkp) ::  shiftof0
  contains
    procedure ::  Initialize     =>    Initialize_N4_NASA_PES
    procedure ::  Output         =>    Output_N4_NASA_PES
    procedure ::  Compute        =>    Compute_N4_NASA_PES_1d
    procedure ::  Potential      =>    N4_NASA_Potential_From_R
    procedure ::  TriatPotential =>    N4_NASA_Potential_From_R_OnlyTriat
  End Type

  logical                                 ,parameter      :: i_Debug_Global = .False.

  integer                                 ,parameter      :: lwork  = 50

  real(rkp)                               ,parameter      :: retd   = 2.76149524d0
  real(rkp)                               ,parameter      :: betatd = 1.075722195d0
  real(rkp)                               ,parameter      :: detd   = 0.0955135476d0
                                                
  real(rkp)                               ,parameter      :: x6     =  75.63d0
  real(rkp)                               ,parameter      :: f2     =  0.129334171818588d0
  real(rkp)                               ,parameter      :: f220   =  5.220578416136951d-2
  real(rkp)                               ,parameter      :: f221   = -1.159697407223882d-2
  real(rkp)                               ,parameter      :: f222   =  1.450724604023408d-3

  real(rkp)                               ,parameter      :: x8     =  2489d0
  real(rkp)                               ,parameter      :: g2     =  0.624578507475241d0
  real(rkp)                               ,parameter      :: g220   =  0.157105504988713d0
  real(rkp)                               ,parameter      :: g221   = -2.188839577146523d-2
  real(rkp)                               ,parameter      :: g222   =  1.003716064096833d-3
  real(rkp)                               ,parameter      :: bswth  =  Four
  real(rkp)                               ,parameter      :: eps    =  1.d-6 ! eps added to rbig in case rbig=0   RLJ 8/09
  real(rkp)                               ,parameter      :: dpmn   = -14.9681846800000d0
  real(rkp)                               ,parameter      :: vibdispt = Two*dsqrt(dabs(dpmn))/(dpmn*dpmn)
  real(rkp)                               ,parameter      :: damp   =  Six
  real(rkp)                               ,parameter      :: damp4  =  damp*damp*damp*damp   
  real(rkp)                               ,parameter      :: damp5  =  damp4*damp            
  real(rkp)                               ,parameter      :: damp6  =  damp5*damp            
  real(rkp)                               ,parameter      :: beta   = -2.18d0
  real(rkp)                               ,parameter      :: beta70 = beta*70.d0          
  real(rkp)                               ,parameter      :: re     =  2.1d0
  real(rkp)                               ,parameter      :: re2    =  Two*re     
  real(rkp) ,dimension(20)                ,parameter      :: coef   = [1.0345018d0,9.7550153d0,2.394727d-1,-4.4270639d-1,4.2206130d0,    &
                                                                       -9.6318443d0,-1.7070046d1,1.5595772d0,9.0354414d-1,-3.7645959d0,  &
                                                                       6.2330259d-1,1.5681477d0,8.1289901d-1,-2.7320998d-1,8.4521948d0,  &
                                                                       5.7363853d0,-1.3138353d0,-7.5160389d-1,1.3283815d-1,0.0038d0]
  integer                                 ,parameter      :: idp  = 10
  integer                                 ,parameter      :: idq  = 5
  integer                                 ,parameter      :: idd  = 5
  integer                                 ,parameter      :: npow = 6
  real(rkp)                               ,parameter      :: ao   = Two
  real(rkp)                               ,parameter      :: a1   = Six
  real(rkp)                               ,parameter      :: a0   = Two
  real(rkp)                               ,parameter      :: einf = -109.036558873442d0
  real(rkp)                               ,parameter      :: c6   = -58.893858362275d0  
  real(rkp)                               ,parameter      :: d    = 2.5d0
  real(rkp)                               ,parameter      :: co   = 49.d0
  real(rkp) ,dimension(idp)               ,parameter      :: pes  = [146.736533087545d0,101.160591312956d0,-251.323812509419d0,-499.020252425371d0,-210.074594675137d0, &
                                                                     135.174682961986d0,214.7693593595d0,425.317131479472d0,570.993394092895d0,266.199319038817d0]
  real(rkp)                               ,parameter      :: a2   = Two
  real(rkp) ,dimension(idq)               ,parameter      :: qm   = [17.3398942423333d0,10.2017829193229d0,-4.7194296899574d0,-9.61643152091987d0,-2.12839249419926d0]
  real(rkp)                               ,parameter      :: a3   = Two
  real(rkp) ,dimension(idd)               ,parameter      :: dip  = [-42.6175278814367d0,873.328188249918d0,984.492872053711d0,518.057930455541d0,124.185667297002d0]
  real(rkp)                               ,parameter      :: dp1   = (873.238d0*dexp(-2.2d0)*2.1d0 - 42.6175d0)
  real(rkp)                               ,parameter      :: dp2   = dp1
  real(rkp)                               ,parameter      :: preqq = 1.5d0*dsqrt(Two)
  real(rkp)                               ,parameter      :: db2td = -Two*detd*betatd !-0.205492086153018d0 
 
  contains

! **************************************************************************************************************
! **************************************************************************************************************
!                                      DEFERRED PROCEDURES for NASA PES
! **************************************************************************************************************
! **************************************************************************************************************
Subroutine Initialize_N4_NASA_PES( This, Input, Atoms, iPES, i_Debug )

  use Input_Class                        ,only:  Input_Type
  use Atom_Class                         ,only:  Atom_Type
  use N2_LeRoy_DiatomicPotential_Class   ,only:  N2_LeRoy_DiatomicPotential_Type
  
  class(N4_NASA_PES_Type)                   ,intent(out)    ::    This
  type(Input_Type)                          ,intent(in)     ::    Input
  type(Atom_Type) ,dimension(:)             ,intent(in)     ::    Atoms  
  integer                                   ,intent(in)     ::    iPES
  logical                         ,optional ,intent(in)     ::    i_Debug
  
  integer                                                   ::    iP
  character(*)                    ,parameter                ::    Name_PES = 'N4_NASA'
  real(rkp)                                                 ::    Temp
  real(rkp)                                                 ::    veccsd
  real(rkp)                                                 ::    veuse
  integer         ,dimension(6,2)                           ::    iA

  logical                                                   ::    i_Debug_Loc
  
  i_Debug_Loc = i_Debug_Global; if ( present(i_Debug) )i_Debug_Loc = i_Debug
  if (i_Debug_Loc) call Logger%Entering( "Initialize_N4_NASA_PES" )
  !i_Debug_Loc   =     Logger%On()
    

  This%Name         =   Name_PES
  This%Initialized  =   .True.
  This%CartCoordFlg =   .True.
  This%NPairs       =   6               ! Setting the number of atom-atom pairs

  iA(1,:) = [1, 2]
  iA(2,:) = [1, 3]
  iA(3,:) = [1, 4]
  iA(4,:) = [2, 3]
  iA(5,:) = [2, 4]
  iA(6,:) = [3, 4] 

  allocate( This%Pairs(This%NPairs) )   ! Allocating the Pairs array which contains the polymorphic Diatomi-Potential associated to each pair
  do iP = 1,This%NPairs
    allocate( N2_LeRoy_DiatomicPotential_Type :: This%Pairs(iP)%Vd  )
  end do

  call diatccsdd(re, veccsd, Temp, 1)    
  if (i_Debug_Loc) call Logger%Write( " ccsd potential at re:     veccsd = ", veccsd )

  call This%Pairs(1)%Vd%Compute_Vd_dVd( re, veuse, Temp )       
  if (i_Debug_Loc) call Logger%Write( " diatomic potential at re: veuse  = ", veuse )

  This%shiftof0 = Two * (veccsd - veuse)
  if (i_Debug_Loc) call Logger%Write( " shift of diatomic zero of energy = ", This%shiftof0/Two, "; shift of VTd = ", This%shiftof0 )                                        

  if (i_Debug_Loc) call Logger%Exiting()

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Output_N4_NASA_PES( This, Unit )

  class(N4_NASA_PES_Type)                 ,intent(in)     ::    This
  integer                                 ,intent(in)     ::    Unit
  
  write(Unit,"('PES Name: ',g0)") This%Name
  write(Unit,"('N4 PES with cta bug fixed')")
  write(Unit,"('acpf energies, 1s2s core and leroy n2 with hyper radius repulsion')")
  
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N4_NASA_Potential_From_R( This, R, Q ) result( V )

  use N2_LeRoy_DiatomicPotential_Class  ,only:  N2_LeRoy_DiatomicPotential_Type

  class(N4_NASA_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp) ,dimension(12)                                   ::    dVdQin
  real(rkp)                                                  ::    t1, t2
  
  call n4fitd(This%shiftof0, Q(1:3), Q(4:6), Q(7:9), Q(10:12), V, dVdQin)

End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Function N4_NASA_Potential_From_R_OnlyTriat( This, R, Q ) result( V )

  class(N4_NASA_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R           !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q           !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                                  ::    V           !< Potential energy in [hartree].

  real(rkp)                                                  ::    t1, t2
  
  V = 0.d0
  
End Function
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine Compute_N4_NASA_PES_1d( This, R, Q, V, dVdR, dVdQ )

  use N2_LeRoy_DiatomicPotential_Class  ,only:  N2_LeRoy_DiatomicPotential_Type

  class(N4_NASA_PES_Type)                       ,intent(in)  ::    This
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    R            !< Distances of atom-atom pairs [bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(in)  ::    Q            !< Atom Cartisian Coordinates [bohr]. Dim=(NAtoms*3) 
  real(rkp)                                     ,intent(out) ::    V            !< Potential energy in [hartree].
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdR         !< Derivative of the potential wrt pair distances [hartree/bohr]. Dim=(NPairs)
  real(rkp) ,dimension(:)          CONTIGUOUS   ,intent(out) ::    dVdQ         !< Derivative of the potential wrt atom coordinates [hartree/bohr]. Dim=(NAtoms*3)

  integer                                                    ::    at, dir
  real(rkp)                                                  ::    t1, t2
  real(rkp), dimension(12)                                   ::    dVdQin
  real(rkp), dimension(3)                                    ::    dVdR4

  dVdR = Zero
  
  call n4fitd(This%shiftof0, Q(1:3), Q(4:6), Q(7:9), Q(10:12), V, dVdQin)

  ! Extract derivatives with respect to Cartesian coordinates of fourth atom
  dVdR4 = dVdQin(10:12)

  ! Form derivatives with respect to Cartesian coordinates of the first three atoms
  ! by applying the chain-rule and by exploiting the fact that all atoms have the 
  ! same mass. Hence dx_4/dx_i = -1, i = 1,2,3
  atom_loop : do at = 1,3 
     dir_loop : do dir = 1,3
        dVdQ(3*(at - 1) + dir) = dVdQin(3*(at - 1) + dir) - dVdR4(dir)     
     enddo dir_loop
  enddo atom_loop

End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


! ! **************************************************************************************************************
! ! **************************************************************************************************************
! !                                   PRIVATE PROCEDURES for NASA PES
! ! **************************************************************************************************************
! ! **************************************************************************************************************
Subroutine n4fitd(shiftof0, xyz1, xyz2, xyz3, xyz4, fit, dfit)
!
! N4 PES including dissociation and N4 Td
! DWS 8/4/09
! xyzi: cartesian coordinates of atom i in a.u.
!
  real(rkp), intent(in)                                   ::    shiftof0
  real(rkp) ,dimension(3)                 ,intent(in)     ::    xyz1, xyz2, xyz3, xyz4
  real(rkp)                               ,intent(out)    ::    fit
  real(rkp) ,dimension(12)                ,intent(out)    ::    dfit

  real(rkp) ,dimension(3,4)                               ::    cart
  real(rkp) ,dimension(4,4)                               ::    vmat
  real(rkp) ,dimension(4,4)                               ::    vmatc
  real(rkp) ,dimension(12,4)                              ::    dpart
  real(rkp)                                               ::    sv
  real(rkp) ,dimension(4)                                 ::    vec
  real(rkp) ,dimension(4)                                 ::    eig
  integer   ,dimension(20)                                ::    iwork
  real(rkp) ,dimension(lwork)                             ::    work  
  integer   ,dimension(4)                                 ::    ifail
  integer                                                 ::    neig
  real(rkp)                                               ::    vl
  real(rkp)                                               ::    ul
  integer                                                 ::    info
  real(rkp) ,dimension(4)                                 ::    fv1, fv2
  real(rkp) ,dimension(4,4)                               ::    cpy
  integer                                                 ::    i, j, k, iok

  vmat  = 1.d-5
  dpart = Zero
                                              
  cart(:,1) = xyz1(:)                                                
  cart(:,2) = xyz2(:)                                                
  cart(:,3) = xyz3(:)                                                
  cart(:,4) = xyz4(:)                                                

  !write(*,*) 'Atom1, ', cart(:,1), '; Atom2, ', cart(:,2), '; Atom3, ', cart(:,3), '; Atom4, ', cart(:,4)


!#ifdef N4_PES_VERBOSE  
  !write(*,*) '** n4fitd **** cart(:,1) = ', cart(:,1)
  !write(*,*) '** n4fitd **** cart(:,2) = ', cart(:,2)
  !write(*,*) '** n4fitd **** cart(:,3) = ', cart(:,3)
  !write(*,*) '** n4fitd **** cart(:,4) = ', cart(:,4)

  !write(*,*) '** n4fitd **** Entering vtdd'
!#endif
  !call vtdd(  xyz1, xyz2, xyz3, xyz4, vmat(1,1), dpart)
  call vtdd(  xyz1, xyz2, xyz3, xyz4, vmat(1,1), dpart(:,1))  
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** Done with vtdd'
  !write(*,*) '** n4fitd **** Entering vn2n2d for the 1st time'
!#endif
  call vn2n2d(xyz1, xyz2, xyz3, xyz4, vmat(2,2), dpart(:,2), 0)  
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** Done with vn2n2d'
  !write(*,*) '** n4fitd **** Entering vn2n2d for the 2nd time'
!#endif
  call vn2n2d(xyz1, xyz3, xyz2, xyz4, vmat(3,3), dpart(:,3), 0)  
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** Done with vn2n2d'
  !write(*,*) '** n4fitd **** Entering vn2n2d for the 3rd time'
!#endif
  call vn2n2d(xyz1, xyz4, xyz2, xyz3, vmat(4,4), dpart(:,4), 0)  
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** Done with vn2n2d'
!#endif

  do i=1,3
    sv           = dpart(i+3,3)
    dpart(i+3,3) = dpart(i+6,3)
    dpart(i+6,3) = sv
    sv           = dpart(i+3,4)
    dpart(i+3,4) = dpart(i+6,4)
    dpart(i+6,4) = dpart(i+9,4)
    dpart(i+9,4) = sv
  end do
!#ifdef N4_PES_VERBOSE  
  !write(*,*) '** n4fitd **** dpart = ', dpart
!#endif
  !stop
  vmat(1,1) = vmat(1,1) - shiftof0                                      
  do i=2,4
    vmat(1,i) = 6.0847866d-02
    vmat(i,1) = 6.0847866d-02
  end do
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** vmat = ', vmat
!#endif   
  do i=1,4
    do j=1,4
      vmatc(j,i) = vmat(j,i)
    end do
  end do
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** vmatc = ', vmatc
!#endif  
  !     
  ! Call to MKL (Find the eigenvalues and eigenvectors)
  ! ---------------------------------------------------
  ! Compute only the 1st eigenvalue   ! MARCO
  ! call dsyevx('V','I','L',4,vmat,4,vl,ul,1,1,zero,neig, eig,vec,4, work,lwork,iwork,ifail,info)
  ! LAPACK  dsyevx (JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO)      
  ! ---------------------------------------------------------------
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** Entering dsyevx'
!#endif
  call dsyevx('V', 'I', 'L', 4, vmat, 4, vl, ul, 1, 1, Zero, neig, eig, vec, 4, work, lwork, iwork, ifail, info)
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** Done with dsyevx '
!#endif

  do i=1,4
    vec(i) = vec(i)**2
  end do
!#ifdef N4_PES_VERBOSE
  !write(*,*) '** n4fitd **** vec = ', vec
!#endif

  if (info .ne. 0) then
    !write(*,*) '** n4fitd **** info = ', info, '; STOPPING!'
    stop
  end if

  fit = eig(1)
!#ifdef N4_PES_VERBOSE  
  !write(*,*) '** n4fitd **** fit = ', fit
!#endif

  do i=1,12

    dfit(i) = dpart(i,1)*vec(1) + dpart(i,2)*vec(2) + dpart(i,3)*vec(3) + dpart(i,4)*vec(4)
!#ifdef N4_PES_VERBOSE
    !write(*,*) '** n4fitd **** dfit(i) = ', dfit(i), ' for i=', i
!#endif

    if ( (dfit(i) .ne. dfit(i)) .or. (dfit(i) .gt. 1.d20) .or. (dfit(i) .lt. -1.d20) ) then
      
!#ifdef N4_PES_VERBOSE
      !write(*,*) '** n4fitd **** got nan for dfit ', i
      !write(*,*) '** n4fitd **** dparts ', dpart(i,1), dpart(i,2), dpart(i,3), dpart(i,4)
      !write(*,*) '** n4fitd **** vecs ', vec(1), vec(2), vec(3), vec(4)
      !write(*,*) '** n4fitd **** vmat: '

      do k=1,4
        !write(*,*) (vmat(k,j), j=1,4)
      end do
!#endif

      if ( (dpart(i,2) .ne. dpart(i,2)) .or. (dpart(i,2) .gt. 1d20) .or. (dpart(i,2) .le. -1d20) ) then
        call vn2n2d( xyz1, xyz2, xyz3, xyz4, vmat(2,2), dpart(1,2), 1)
      else if ( (dpart(i,3) .ne. dpart(i,3)) .or. (dpart(i,3) .gt. 1d20) .or. (dpart(i,3) .le. -1d20) ) then
        call vn2n2d( xyz1, xyz3, xyz2, xyz4, vmat(3,3), dpart(1,3), 1)
      else if ( (dpart(i,4) .ne. dpart(i,4)) .or. (dpart(i,4) .gt. 1d20) .or. (dpart(i,4) .lt. -1d20) ) then
        call vn2n2d( xyz1, xyz4, xyz2, xyz3, vmat(4,4), dpart(1,4), 1)
      end if

      stop
    end if

  end do

 
  if (dfit(1) .ge. Zero) then
    iok = 1
  else
    if (dfit(1) .lt. Zero) then
      iok = 1
    else
      !write(*,*) '** n4fitd **** have not a number for first der: '
      do i=1,4
        !write(*,*) '** n4fitd **** i = ', i, '; vec(i) = ', vec(i), '; dpart(1,i)', dpart(1,i)
        if ( (dpart(1,i) .le. Zero) .or. (dpart(1,i) .gt. Zero) )then
          iok=iok+1
        else
          if      (i .eq. 2) then
            !write(*,*) '** n4fitd **** i=2; Entering vn2n2d'
            call vn2n2d(xyz1, xyz2, xyz3, xyz4, vmat(2,2), dpart(1,2), 1)
          else if (i .eq. 3) then
            !write(*,*) '** n4fitd **** i=3; Entering vn2n2d'
            call vn2n2d(xyz1, xyz3, xyz2, xyz4, vmat(3,3), dpart(1,3), 1)
          else if (i .eq. 4) then
            !write(*,*) '** n4fitd **** i=4; Entering vn2n2d'
            call vn2n2d(xyz1, xyz4, xyz2, xyz3, vmat(4,4), dpart(1,4), 1)
          end if
        end if
      end do
      stop
    end if
  end if
  
  return
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine vtdd(x1, x2, x3, x4, fit, dfit)
! compute N4 Td potential and derivatives

  real(rkp) ,dimension(3)                 ,intent(in)     ::    x1, x2, x3, x4
  real(rkp)                               ,intent(out)    ::    fit
  real(rkp) ,dimension(12)                ,intent(out)    ::    dfit
           
  real(rkp) ,dimension(12)                                ::    dr12, dr13, dr14, dr23, dr24, dr34
  real(rkp) ,dimension(12)                                ::    dvx12, dvx13, dvx14, dvx23, dvx24, dvx34
  real(rkp)                                               ::    dvx12dr, dvx13dr, dvx14dr, dvx23dr, dvx24dr, dvx34dr
  real(rkp)                                               ::    ex12, ex13, ex14, ex23, ex24, ex34
  real(rkp)                                               ::    vx12, vx13, vx14, vx23, vx24, vx34
  real(rkp)                                               ::    r12, r13, r14, r23, r24, r34
  integer                                                 ::    i
  
  dr12 = Zero 
  dr13 = Zero 
  dr14 = Zero 
  dr23 = Zero 
  dr24 = Zero 
  dr34 = Zero 
  
  dvx12 = Zero 
  dvx13 = Zero 
  dvx14 = Zero 
  dvx23 = Zero 
  dvx24 = Zero 
  dvx34 = Zero 
           
  r12 = Zero
  r13 = Zero
  r14 = Zero
  r23 = Zero
  r24 = Zero
  r34 = Zero
  do i=1,3
    r12 = r12 + ( x1(i) - x2(i) )**2
    r13 = r13 + ( x1(i) - x3(i) )**2
    r14 = r14 + ( x1(i) - x4(i) )**2
    r23 = r23 + ( x2(i) - x3(i) )**2
    r24 = r24 + ( x2(i) - x4(i) )**2
    r34 = r34 + ( x3(i) - x4(i) )**2
  end do
  r12 = sqrt(r12)
  r13 = sqrt(r13)
  r14 = sqrt(r14)
  r23 = sqrt(r23)
  r24 = sqrt(r24)
  r34 = sqrt(r34)
!#ifdef N4_PES_VERBOSE  
  !write(*,*) '**** vtdd **** r12 = ', r12, '; r13 = ', r13, '; r14 = ', r14, '; r23 = ', r23, '; r24 = ', r24, '; r34 = ', r34
!#endif

  do i=1,3
    dr12(i)   = (x1(i) - x2(i)) / r12
    dr13(i)   = (x1(i) - x3(i)) / r13
    dr14(i)   = (x1(i) - x4(i)) / r14
    dr23(i+3) = (x2(i) - x3(i)) / r23
    dr24(i+3) = (x2(i) - x4(i)) / r24
    dr34(i+6) = (x3(i) - x4(i)) / r34
    dr12(i+3) = -dr12(i)
    dr13(i+6) = -dr13(i)
    dr14(i+9) = -dr14(i)
    dr23(i+6) = -dr23(i+3)
    dr24(i+9) = -dr24(i+3)
    dr34(i+9) = -dr34(i+6)
  end do
!#ifdef N4_PES_VERBOSE  
  !write(*,*) '**** vtdd **** dr12 = ', dr12, '; dr13 = ', dr13, '; dr14 = ', dr14, '; dr23 = ', dr23, '; dr24 = ', dr24, '; dr34 = ', dr34
!#endif

  ex12 = exp( -betatd * (r12 - retd) )
  vx12 = (ex12 - One)**2 * detd
  ex13 = exp( -betatd * (r13 - retd) )
  vx13 = (ex13 - One)**2 * detd
  ex14 = exp( -betatd * (r14 - retd) )
  vx14 = (ex14 - One)**2 * detd
  ex23 = exp( -betatd * (r23 - retd) )
  vx23 = (ex23 - One)**2 * detd
  ex24 = exp( -betatd * (r24 - retd) )
  vx24 = (ex24 - One)**2 * detd
  ex34 = exp( -betatd * (r34 - retd) )
  vx34 = (ex34 - One)**2 * detd
  fit  = vx12 + vx13 + vx14 + vx23 + vx24 + vx34 - 218.47601249d0 + 0.10d0 ! from my ccsd(t) caln and empirical shift to get X peak height about right
!#ifdef N4_PES_VERBOSE
  !write(*,*) '**** vtdd **** fit = ', fit
!#endif

! compute derivatives wrt atomic cartesian coordinates
  dvx12dr = db2td * ex12 * (ex12 - One)
  dvx13dr = db2td * ex13 * (ex13 - One)      
  dvx14dr = db2td * ex14 * (ex14 - One)      
  dvx23dr = db2td * ex23 * (ex23 - One)      
  dvx24dr = db2td * ex24 * (ex24 - One)      
  dvx34dr = db2td * ex34 * (ex34 - One) 
  
  do i=1,3
    dvx12(i)   = dvx12dr * dr12(i)
    dvx13(i)   = dvx13dr * dr13(i)
    dvx14(i)   = dvx14dr * dr14(i)
    dvx23(i+3) = dvx23dr * dr23(i+3)
    dvx24(i+3) = dvx24dr * dr24(i+3)
    dvx34(i+6) = dvx34dr * dr34(i+6)   
    dvx12(i+3) = dvx12dr * dr12(i+3)
    dvx13(i+6) = dvx13dr * dr13(i+6)
    dvx14(i+9) = dvx14dr * dr14(i+9)
    dvx23(i+6) = dvx23dr * dr23(i+6)
    dvx24(i+9) = dvx24dr * dr24(i+9)
    dvx34(i+9) = dvx34dr * dr34(i+9)   
  end do

  do i=1,12
    dfit(i) = dvx12(i) + dvx13(i) + dvx14(i) + dvx23(i) + dvx24(i) + dvx34(i)
  end do     
!#ifdef N4_PES_VERBOSE  
  !write(*,*) '**** vtdd **** dfit = ', dfit
!#endif

  return
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine vn2n2d(x1, x2, x3, x4, fit, dfit, iprt)

  use N2_LeRoy_DiatomicPotential_Class  ,only:  N2_LeRoy_DiatomicPotential_Type

  real(rkp) ,dimension(3)                 ,intent(in)     :: x1, x2, x3, x4
  real(rkp)                               ,intent(out)    :: fit
  real(rkp) ,dimension(12)                ,intent(out)    :: dfit
  integer                                 ,intent(in)     :: iprt

  integer                                                 :: i, j, ipst
  type(N2_LeRoy_DiatomicPotential_Type)                   :: DiatPot
  real(rkp) ,dimension(12)                                :: dr1, dr2, dr3, dr4, dr5, dr6
  real(rkp) ,dimension(12)                                :: dxa1, dxa2, dxa3, dxb1, dxb2, dxb3
  real(rkp) ,dimension(12)                                :: dxa, dxb  
  real(rkp) ,dimension(3)                                 :: x12, x13, x14, x23, x24, x34
  real(rkp) ,dimension(12)                                :: dtheta1, dtheta2
  real(rkp) ,dimension(12)                                :: drbig
  real(rkp) ,dimension(12)                                :: dphi
  real(rkp) ,dimension(12)                                :: ddot3
  real(rkp) ,dimension(3)                                 :: xcm
  real(rkp) ,dimension(3,2)                               :: rv
  real(rkp) ,dimension(3)                                 :: rab, rcd
  real(rkp)                                               :: dp1, dp2 
  real(rkp)                                               :: dot1, dot2, dot3
  real(rkp)                                               :: rbig
  real(rkp)                                               :: theta1arg, theta2arg, theta1, theta2
  real(rkp)                                               :: xa1, xa2, xa3, xb1, xb2, xb3
  real(rkp)                                               :: xa, xb
  real(rkp)                                               :: rxaxb
  real(rkp)                                               :: phi, tphi
  real(rkp)                                               :: trat1, tpre1, xpre1, ypre1, zpre1, trat2, tpre2, xpre2, ypre2, zpre2
  real(rkp)                                               :: dws1, dws2
  real(rkp)                                               :: rxa1, rxa2, rxa3, rxb1, rxb2, rxb3
  real(rkp)                                               :: v12, v34
  real(rkp)                                               :: dv12dr1, dv34dr2
  real(rkp)                                               :: ddp1dr1, ddp2dr2
  real(rkp)                                               :: q1, q2
  real(rkp)                                               :: dq1dr1, dq2dr2
  real(rkp)                                               :: vibdispa, vibdispb
  real(rkp)                                               :: dvibdisp1, dvibdisp2
  real(rkp)                                               :: q12
  real(rkp)                                               :: ct1, ct2
  real(rkp)                                               :: st1, st2
  real(rkp)                                               :: cph, sph
  real(rkp)                                               :: cp, sp
  real(rkp)                                               :: cst1, cst2
  real(rkp)                                               :: y2a, y2b, z2a, z2b, y2, y3, y4, y5a, y5b, y5, y6a, y6b, y6
  real(rkp)                                               :: p2a, p2b, p21a, p21b, p22a, p22b
  real(rkp)                                               :: c6000p, c8000p, c6000, c8000
  real(rkp)                                               :: dc6000r1, dc6000r2, dc6000t1, dc6000t2, dc8000r1, dc8000r2, dc8000t1, dc8000t2
  real(rkp)                                               :: qq
  real(rkp)                                               :: dqqdr1, dqqdr2, dqqdph, dqqdt1, dqqdt2
  real(rkp)                                               :: denom1, denom2a, denom2, denom3
  real(rkp)                                               :: vlr, vlro
  real(rkp)                                               :: ex1
  real(rkp)                                               :: pairwise
  real(rkp)                                               :: sargtyp
  real(rkp)                                               :: swth, dswth, dswthab, dswthc, swth1, dswthdr1, dswthdr2, dswthdrb
  real(rkp)                                               :: dvlrdr1, dvlrdr2, dvlrdt1, dvlrdt2, dvlrdrb, dvlrph
  real(rkp)                                               :: r1re, r2re, r1re2, r2re2
  real(rkp)                                               :: ddr1, ddr2, dr12, dr1p2
  real(rkp)                                               :: bra, brb
  real(rkp)                                               :: dpairdr3, dpairdr4, dpairdr5, dpairdr6
  real(rkp)                                               :: dbradr1, dbradr2, dbradr3, dbradr4, dbradr5, dbradr6
  real(rkp)                                               :: dbrbdrb, dbrbdr1, dbrbdr2, dbrbdr3, dbrbdr4, dbrbdr5, dbrbdr6
  real(rkp)                                               :: func1a, func2a, func1, func2
  real(rkp)                                               :: dfunc1dr1, dfunc1dr2, dfunc1dt1, dfunc1dt2, dfunc1dph, dfunc2dr1, dfunc2dr2, dfunc2dt1, dfunc2dt2, dfunc2dph
  real(rkp)                                               :: dfitr1, dfitr2, dfitt1, dfitt2, dfitph, dfitrb, dfitr3, dfitr4, dfitr5, dfitr6
  real(rkp)                                               :: bot
  real(rkp)                                               :: r1, r2, r3, r4, r5, r6


  !write(*,*) "**** vn2n2d *********************************************************************************************************"


  dr1 = Zero 
  dr2 = Zero 
  dr3 = Zero 
  dr4 = Zero 
  dr5 = Zero 
  dr6 = Zero 

  !!! drab = Zero 
  !!! drcd = Zero 
  !!! dxcm = Zero 
  !!! drv  = Zero 

  dxa1 = Zero 
  dxa2 = Zero 
  dxa3 = Zero 
  dxb1 = Zero 
  dxb2 = Zero 
  dxb3 = Zero 
  
  dxa  = Zero 
  dxb  = Zero 

  if (iprt .ne. 0) then
    !write(*,*) 'input carts: '
    do i=1,3
      !write(*,*) x1(i), x2(i), x3(i), x4(i)
    end do
  end if

  r1 = Zero
  r5 = Zero
  r3 = Zero
  r6 = Zero
  r4 = Zero
  r2 = Zero
  do i=1,3
    x12(i) = x1(i)-x2(i)
    x13(i) = x1(i)-x3(i)
    x14(i) = x1(i)-x4(i)
    x23(i) = x2(i)-x3(i)
    x24(i) = x2(i)-x4(i)
    x34(i) = x3(i)-x4(i)

    r1 = r1 + x12(i)**2
    r5 = r5 + x13(i)**2
    r3 = r3 + x14(i)**2
    r6 = r6 + x23(i)**2
    r4 = r4 + x24(i)**2
    r2 = r2 + x34(i)**2
  end do
  r1 = sqrt(r1)
  r5 = sqrt(r5)
  r3 = sqrt(r3)
  r6 = sqrt(r6)
  r4 = sqrt(r4)
  r2 = sqrt(r2)
!#ifdef N4_PES_VERBOSE  
   !write(*,*) '**** vn2n2d ** r1 = ', r1, '; r2 = ', r2, '; r3 = ', r3, '; r4 = ', r4, '; r5 = ', r5, '; r6 = ', r6 
!#endif

  dot1 = Zero
  dot2 = Zero
  rbig = Zero
  do j=1,3
    rab(j)  = Half * ( x1(j) + x2(j) )
    rcd(j)  = Half * ( x3(j) + x4(j) )
    xcm(j)  = rcd(j) - rab(j)
    rv(j,1) =  x1(j) - rab(j)
    rv(j,2) =  x3(j) - rcd(j)
    rbig    = rbig + ( rcd(j) - rab(j) )**2
    dot1    = dot1 + ( x1(j)  - x2(j)  ) * ( rcd(j) - rab(j) )
    dot2    = dot2 + ( x3(j)  - x4(j)  ) * ( rcd(j) - rab(j) )
  end do
  rbig = sqrt(rbig)
!#ifdef N4_PES_VERBOSE
  !write(*,*) '**** vn2n2d ** dot1 = ', dot1, '; dot2 = ', dot2
  !write(*,*) '**** vn2n2d ** rbig = ', rbig
!#endif

  theta1arg = dot1 / (r1 * rbig)                                          
  if ( abs(theta1arg) .gt. One ) then                                     
    if ( theta1arg .gt. Zero ) then                                         
      theta1arg = One                                                  
    else                                                             
      theta1arg = -One                                                  
    end if                                                           
  end if    
!#ifdef N4_PES_VERBOSE
  !write(*,*) '**** vn2n2d ** theta1arg = ', theta1arg
!#endif

  theta1 = acos(theta1arg)                                            
  if (theta1 .ne. theta1) then
    !write(*,*) '**** vn2n2d ** theta1 nan ', dot1, r1, rbig, theta1arg
  end if
  theta2arg = dot2 / (r2 * rbig)                                          
  if ( abs(theta2arg) .gt. One ) then                                     
    if ( theta2arg .gt. Zero ) then                                         
      theta2arg = One                                                   
    else                                                             
      theta2arg = -One                                                  
    end if                                                           
  end if  
  theta2 = acos(theta2arg)  

  !write(*,*) "**** vn2n2d ** rv  = ", rv
  !write(*,*) "**** vn2n2d ** xcm = ", xcm

  xa1    = rv(2,1) * xcm(3) - rv(3,1) * xcm(2)
  xa2    = rv(3,1) * xcm(1) - rv(1,1) * xcm(3)
  xa3    = rv(1,1) * xcm(2) - rv(2,1) * xcm(1)
  xb1    = rv(2,2) * xcm(3) - rv(3,2) * xcm(2)
  xb2    = rv(3,2) * xcm(1) - rv(1,2) * xcm(3)
  xb3    = rv(1,2) * xcm(2) - rv(2,2) * xcm(1)
  xa     = sqrt(xa1**2 + xa2**2 + xa3**2)
  xb     = sqrt(xb1**2 + xb2**2 + xb3**2)
  rxaxb  = One / (xa * xb)
  
  !write(*,*) "**** vn2n2d ** xa = ", xa
  !write(*,*) "**** vn2n2d ** xb = ", xb

  if ( (abs(xa) .le. 1.d-6) .or. (abs(xb) .le. 1.d-6) ) then
    phi  = Zero
    tphi = Zero
    dot3 = Zero
  else
    dot3 = ( xa1 * xb1 + xa2 * xb2 + xa3 * xb3 ) * rxaxb
    tphi = One / ( sqrt( abs(One - dot3**2) ) + 1.d-14)                         
    if (iprt .ne. 0) then
      !write(*,*) '**** vn2n2d ** dot3 = ',   dot3, xa1, xb1, xa2, xb2, xa3, xb3, rxaxb
      !write(*,*) '**** vn2n2d ** dot3-1d0 ', dot3-One, abs(abs(dot3)-One)
    end if
    ipst = 0
    if ( abs(abs(dot3) - One) .lt. 1.d-10) then                              
      if (dot3 .gt. Zero) then                                             
        phi  = Zero                                                        
        ipst = 1
      else                                                            
        phi  = pi                                                         
        ipst = 2
      end if                                                          
    else                                                             
      phi  = acos(dot3)
      ipst = 3
    end if                                                           
    if (iprt .gt. 0) then
      !write(*,*) '**** vn2n2d ** phi, ipst ', phi, ipst, pi
    end if
  end if

  ! derivatives of r's and Jacobi coordinates wrt atomic cartesians    (RLJ)
  do i=1,3
    dr1(i)     = x12(i) / r1
    dr1(i+3)   = -dr1(i)
    dr2(i+6)   = x34(i) / r2                                              
    dr2(i+9)   = -dr2(i+6)
    dr3(i)     = (x1(i) - x4(i)) / r3
    dr3(i+9)   = -dr3(i)
    dr4(i+3)   = (x2(i) - x4(i)) / r4
    dr4(i+9)   = -dr4(i+3)
    dr5(i)     = (x1(i) - x3(i)) / r5
    dr5(i+6)   = -dr5(i)
    dr6(i+3)   = (x2(i) - x3(i)) / r6
    dr6(i+6)   = -dr6(i+3)
    drbig(i)   = -Half * xcm(i)  / rbig
    drbig(i+3) =  drbig(i)
    drbig(i+6) = -drbig(i)
    drbig(i+9) = -drbig(i)
  end do
  
  trat1 =  dot1 / (r1 * rbig)
  trat2 =  dot2 / (r2 * rbig)
  tpre1 =  One  / sqrt(One - trat1**2)
  xpre1 =  One  / (r1 * rbig)
  ypre1 =  Half * dot1 / (r1 * rbig**2)
  zpre1 =  dot1 / (rbig * r1**2)
  tpre2 =  One  / sqrt(One - trat2**2)
  xpre2 =  One  / (r2 * rbig)
  ypre2 =  Half * dot2 / (r2 * rbig**2)
  zpre2 =  dot2 / (rbig * r2**2)
  dws1  = -One  / sqrt(abs(r1 * r1 * rbig * rbig - dot1 * dot1) + 1.d-14)              
  dws2  = -One  / sqrt(abs(r2 * r2 * rbig * rbig - dot2 * dot2) + 1.d-14)         

  if (iprt .ne. 0)then
    !write(*,*) '**** vn2n2d ** dws1, dws2 ', dws1, dws2
    !write(*,*) '**** vn2n2d ** sqrt arg ',   r1*r1*rbig*rbig, dot1*dot1, r1*r1*rbig*rbig-dot1*dot1
    !write(*,*) '**** vn2n2d ** r1,rbig ',    r1, r2, rbig
    !write(*,*) '**** vn2n2d ** dot1 ',  dot1
    !write(*,*) '**** vn2n2d ** rab ',   rab
    !write(*,*) '**** vn2n2d ** rcd ',   rcd
    !write(*,*) '**** vn2n2d ** x12 ',   x12
    !write(*,*) '**** vn2n2d ** dr1 ',   dr1
    !write(*,*) '**** vn2n2d ** drbig ', drbig
  end if

  do i=1,3
    dtheta1(i)   = dws1 * ( rcd(i) - rab(i) - Half * x12(i) - dot1 * ( (dr1(i)   / r1) + (drbig(i)   / rbig) ) )                                           
    dtheta1(i+3) = dws1 * ( rab(i) - rcd(i) - Half * x12(i) - dot1 * ( (dr1(i+3) / r1) + (drbig(i+3) / rbig) ) )                           
    dtheta1(i+6) = dws1 * (  Half * x12(i) - dot1 * ( (dr1(i+6) / r1) + (drbig(i+6) / rbig) ) )                                         
    dtheta1(i+9) = dws1 * (  Half * x12(i) - dot1 * ( (dr1(i+9) / r1) + (drbig(i+9) / rbig) ) )                                         

    dtheta2(i)   = dws2 * ( -Half * x34(i) - dot2 * ( (dr2(i)   / r2) + (drbig(i)   / rbig) ) )                                           
    dtheta2(i+3) = dws2 * ( -Half * x34(i) - dot2 * ( (dr2(i+3) / r2) + (drbig(i+3) / rbig) ) )                           
    dtheta2(i+6) = dws2 * ( rcd(i) - rab(i) + Half * x34(i) - dot2 * ( (dr2(i+6) / r2) + (drbig(i+6) / rbig) ) )                                         
    dtheta2(i+9) = dws2 * ( rab(i) - rcd(i) + Half * x34(i) - dot2 * ( (dr2(i+9) / r2) + (drbig(i+9) / rbig) ) )                                         
  end do

  dxa1(2)  =  Half*(rv(3,1) + xcm(3))
  dxa1(3)  = -Half*(rv(2,1) + xcm(2)) 
  dxa1(5)  =  Half*(rv(3,1) - xcm(3)) 
  dxa1(6)  = -Half*(rv(2,1) - xcm(2))    
  dxa1(8)  = -Half*rv(3,1)
  dxa1(9)  =  Half*rv(2,1)
  dxa1(11) = -Half*rv(3,1)
  dxa1(12) =  Half*rv(2,1)
  dxb1(2)  =  Half*rv(3,2)
  dxb1(3)  = -Half*rv(2,2)
  dxb1(5)  =  Half*rv(3,2)
  dxb1(6)  = -Half*rv(2,2)   
  dxb1(8)  = -Half*(rv(3,2) - xcm(3))
  dxb1(9)  =  Half*(rv(2,2) - xcm(2)) 
  dxb1(11) = -Half*(rv(3,2) + xcm(3)) 
  dxb1(12) =  Half*(rv(2,2) + xcm(2))    

  dxa2(1)  = -Half*(rv(3,1) + xcm(3))
  dxa2(3)  =  Half*(rv(1,1) + xcm(1)) 
  dxa2(4)  = -Half*(rv(3,1) - xcm(3)) 
  dxa2(6)  =  Half*(rv(1,1) - xcm(1))    
  dxa2(7)  =  Half*rv(3,1)
  dxa2(9)  = -Half*rv(1,1)
  dxa2(10) =  Half*rv(3,1)
  dxa2(12) = -Half*rv(1,1)
  dxb2(1)  = -Half*rv(3,2)
  dxb2(3)  =  Half*rv(1,2)
  dxb2(4)  = -Half*rv(3,2)
  dxb2(6)  =  Half*rv(1,2)   
  dxb2(7)  =  Half*(rv(3,2) - xcm(3))
  dxb2(9)  = -Half*(rv(1,2) - xcm(1)) 
  dxb2(10) =  Half*(rv(3,2) + xcm(3)) 
  dxb2(12) = -Half*(rv(1,2) + xcm(1))    

  dxa3(1)  =  Half*(rv(2,1) + xcm(2))
  dxa3(2)  = -Half*(rv(1,1) + xcm(1)) 
  dxa3(4)  =  Half*(rv(2,1) - xcm(2)) 
  dxa3(5)  = -Half*(rv(1,1) - xcm(1))    
  dxa3(7)  = -Half*rv(2,1)
  dxa3(8)  =  Half*rv(1,1)
  dxa3(10) = -Half*rv(2,1)
  dxa3(11) =  Half*rv(1,1)
  dxb3(1)  =  Half*rv(2,2)
  dxb3(2)  = -Half*rv(1,2)
  dxb3(4)  =  Half*rv(2,2)
  dxb3(5)  = -Half*rv(1,2)   
  dxb3(7)  = -Half*(rv(2,2) - xcm(2))
  dxb3(8)  =  Half*(rv(1,2) - xcm(1)) 
  dxb3(10) = -Half*(rv(2,2) + xcm(2)) 
  dxb3(11) =  Half*(rv(1,2) + xcm(1))    

  rxa1 = xa1 / xa
  rxa2 = xa2 / xa
  rxa3 = xa3 / xa
  rxb1 = xb1 / xb
  rxb2 = xb2 / xb
  rxb3 = xb3 / xb

  dxa(1)  = rxa2*dxa2(1)  + rxa3*dxa3(1)
  dxa(4)  = rxa2*dxa2(4)  + rxa3*dxa3(4)
  dxa(7)  = rxa2*dxa2(7)  + rxa3*dxa3(7)
  dxa(10) = rxa2*dxa2(10) + rxa3*dxa3(10)
  dxa(2)  = rxa1*dxa1(2)  + rxa3*dxa3(2)
  dxa(5)  = rxa1*dxa1(5)  + rxa3*dxa3(5)
  dxa(8)  = rxa1*dxa1(8)  + rxa3*dxa3(8)
  dxa(11) = rxa1*dxa1(11) + rxa3*dxa3(11)
  dxa(3)  = rxa1*dxa1(3)  + rxa2*dxa2(3)
  dxa(6)  = rxa1*dxa1(6)  + rxa2*dxa2(6)
  dxa(9)  = rxa1*dxa1(9)  + rxa2*dxa2(9)
  dxa(12) = rxa1*dxa1(12) + rxa2*dxa2(12)
  dxb(1)  = rxb2*dxb2(1)  + rxb3*dxb3(1)
  dxb(4)  = rxb2*dxb2(4)  + rxb3*dxb3(4)
  dxb(7)  = rxb2*dxb2(7)  + rxb3*dxb3(7)
  dxb(10) = rxb2*dxb2(10) + rxb3*dxb3(10)
  dxb(2)  = rxb1*dxb1(2)  + rxb3*dxb3(2)
  dxb(5)  = rxb1*dxb1(5)  + rxb3*dxb3(5)
  dxb(8)  = rxb1*dxb1(8)  + rxb3*dxb3(8)
  dxb(11) = rxb1*dxb1(11) + rxb3*dxb3(11)
  dxb(3)  = rxb1*dxb1(3)  + rxb2*dxb2(3)
  dxb(6)  = rxb1*dxb1(6)  + rxb2*dxb2(6)
  dxb(9)  = rxb1*dxb1(9)  + rxb2*dxb2(9)
  dxb(12) = rxb1*dxb1(12) + rxb2*dxb2(12)

  !write(*,*) "**** vn2n2d ** dxa = ", dxa
  !write(*,*) "**** vn2n2d ** dxb = ", dxb

  do i=1,12
    ddot3(i) = (xa1 * dxb1(i) + xb1 * dxa1(i) + xa2 * dxb2(i) + xb2 * dxa2(i) + xa3 * dxb3(i) + xb3 * dxa3(i)) * rxaxb - dot3 * rxaxb * (xa * dxb(i) + xb * dxa(i))
    dphi(i)  = -tphi * ddot3(i)                                          
  end do
  if (iprt .ne. 0) then
    !write(*,*) '**** vn2n2d ** tphi = ', tphi
    do i=1,12
      !write(*,*) ddot3(i)
    end do
  end if
  
  ! get diatomic potential curves
  if (iprt .ne. 0) then                                                 
    !write(*,*) r1, r2, rbig, theta1, theta2, phi                          
    !write(*,*)('**** vn2n2d **  jacobi ders ')
    do i=1,12
      !write(*,*) dr1(i), dr2(i), drbig(i), dtheta1(i), dtheta2(i), dphi(i)
    end do
    !write(*,*)('**** vn2n2d ** other rij ders ')
    do i=1,12
      !write(*,*) dr3(i), dr4(i), dr5(i), dr6(i)
    end do
  end if                                                            

  call DiatPot%Compute_Vd_dVd( r1, v12, dv12dr1 )
  call DiatPot%Compute_Vd_dVd( r2, v34, dv34dr2 )
  ! dpi: dipole polarizability
  call diatccsdd(r1, dp1, ddp1dr1, 3)
  call diatccsdd(r2, dp2, ddp2dr2, 3)
  ! qi: quadrupole moment
  call diatccsdd(r1, q1,  dq1dr1,  2)
  call diatccsdd(r2, q2,  dq2dr2,  2)
 
  vibdispa = dp1*dp2/(dpmn*dpmn)
  bot      = One/(sqrt(abs(dp1)) + sqrt(abs(dp2)))   
  vibdispb = vibdispa*sqrt(abs(dpmn))*Two*bot   

  dvibdisp1 = ddp1dr1 * vibdispb * (One - Half * sqrt(abs(dp1)) * bot) / dp1    
  dvibdisp2 = ddp2dr2 * vibdispb * (One - Half * sqrt(abs(dp2)) * bot) / dp2    
  c6000p    = x6 * vibdispb
  c8000p    = x8 * vibdispb
  q12       = q1 * q2
  ct1 = cos(theta1)
  ct2 = cos(theta2)
  st1 = sin(theta1)
  st2 = sin(theta2)
  cph = cos(phi)
  sph = sin(phi)
 
  cp   = cos(Two * phi)
  sp   = sin(Two * phi)
  cst1 = ct1 * st1
  cst2 = ct2 * st2
  y2a  = ct1 * ct1
  y2b  = ct2 * ct2
  z2a  = st1 * st1
  z2b  = st2 * st2

  y2  = y2a + y2b
  y3  = y2a * y2b
  y4  = z2a * z2b * cp 
  y5a = y2a * y2a
  y5b = y2b * y2b
  y5  = y5a + y5b
  y6a = z2a * cp                                                            
  y6b = z2b * cp                                                            
  y6  = y6a + y6b

  p2a  = Half * (Three * y2a - One)
  p2b  = Half * (Three * y2b - One)
  p21a = -Three * cst1 
  p21b = -Three * cst2 
  p22a =  Three * z2a
  p22b =  Three * z2b

  c6000    = -c6000p         * (One + f2 * (p2a + p2b) + f220 * p2a * p2b + f221 * p21a * p21b + f222 * p22a * p22b)
  dc6000r1 = -dvibdisp1 * x6 * (One + f2 * (p2a + p2b) + f220 * p2a * p2b + f221 * p21a * p21b + f222 * p22a * p22b)                              
  dc6000r2 = -dvibdisp2 * x6 * (One + f2 * (p2a + p2b) + f220 * p2a * p2b + f221 * p21a * p21b + f222 * p22a * p22b)                              

  dc6000t1 = -c6000p * (cst1 * ( -Three * f2 - Three * f220 * p2b + 18d0 * f222 * z2b) + Nine * f221 * cst2 * (y2a - z2a) )
  dc6000t2 = -c6000p * (cst2 * ( -Three * f2 - Three * f220 * p2a + 18d0 * f222 * z2a) + Nine * f221 * cst1 * (y2b - z2b) )
  c8000    = -c8000p * (One + g2 * (p2a + p2b) + g220 * p2a * p2b + g221 * p21a * p21b + g222 * p22a * p22b)

  dc8000r1 = -dvibdisp1 * x8 * (One + g2 * (p2a + p2b) + g220 * p2a * p2b + g221 * p21a * p21b + g222 * p22a * p22b)                              
  dc8000r2 = -dvibdisp2 * x8 * (One + g2 * (p2a + p2b) + g220 * p2a * p2b + g221 * p21a * p21b + g222 * p22a * p22b)

  dc8000t1 = -c8000p * (cst1 * (-Three * g2 - Three * g220 * p2b + 18d0 * g222 * z2b) + Nine * g221 * cst2 * (y2a - z2a) )
  dc8000t2 = -c8000p * (cst2 * (-Three * g2 - Three * g220 * p2a + 18d0 * g222 * z2a) + Nine * g221 * cst1 * (y2b - z2b) )

  qq = preqq * (Half * p2a * p2b -(p21a * p21b / Nine) * cph + (p22a * p22b / 72d0) * cp) * q12
  
  if (iprt .ne. 0) then
    !write(*,*) '**** vn2n2d ** q1,q2 ', qq, dq1dr1, dq2dr2, q1, q2
  end if

  if (qq .ne. Zero) then                                                 
    dqqdr1 = qq * (dq1dr1 / q1)                                          
    dqqdr2 = qq * (dq2dr2 / q2)                                          
  else                                                              
    dqqdr1 = Zero                                                      
    dqqdr2 = Zero                                                       
  end if                                                            
  dqqdph = preqq * ( (p21a * p21b / Nine) * sph - (p22a * p22b / 36.d0) * sp) * q12
  dqqdt1 = preqq * ( -1.5d0 * cst1 * p2b - cst2 * (y2a - z2a) * cph + 0.25d0 * cst1 * z2b * cp) * q12
  dqqdt2 = preqq * ( -1.5d0 * cst2 * p2a - cst1 * (y2b - z2b) * cph + 0.25d0 * cst2 * z2a * cp) * q12

  denom1  = One / ( rbig**6 + damp6)
  denom2a = One / ( rbig**4 + damp4)
  denom2  =  denom2a * denom2a
  denom3  = One / ( rbig**5 + damp5)

  vlr = (c6000 * denom1) + (c8000 * denom2) + (qq * denom3)
  ex1 = exp(beta * r3) + exp(beta * r4) + exp(beta * r5) + exp(beta * r6)
  pairwise = ex1 * 7.d1
  sargtyp  = (r1 + r2) / (rbig + eps)  
!#ifndef N4_EXP_FIX
  swth     = One / (One + exp(bswth * (sargtyp - 1.2d0))) 
!#else
  swth     = One / (One + exp(min(bswth * (sargtyp - 1.2d0),300.d0))) 
!#endif
  if (swth .gt. 1.d-12) then                                             
    dswth = -swth * swth * bswth * exp(bswth * (sargtyp - 1.2d0))                 
  else                                                              
     dswth=Zero                                                        
  end if                                                            
  if (iprt .ne. 0) then
    !write(*,*) '**** vn2n2d ** dswth ', swth, bswth, sargtyp
  end if
  dswthab = dswth / (rbig + eps)                                          
  if (iprt .ne. 0) then
    !write(*,*) '**** vn2n2d ** dswthab ', dswthab, dswth, rbig, eps
  end if
  swth1    = bswth * swth * (One - swth)
  dswthc   = -dswthab * (r1 + r2) / (rbig + eps)                                
  dswthdr1 =   -swth1             / (rbig + eps)
  dswthdr2 =   -swth1             / (rbig + eps)
  dswthdrb =    swth1 * (r1 + r2) / (rbig + eps)**2
  vlro = vlr
  vlr  = vlr * swth                                                      

  dvlrdr1 = vlro * dswthab                                              
  dvlrdr2 = vlro * dswthab                                              
  dvlrdt1 = swth * (dc6000t1 * denom1 + dc8000t1 * denom2 + dqqdt1 * denom3)
  dvlrdt2 = swth * (dc6000t2 * denom1 + dc8000t2 * denom2 + dqqdt2 * denom3)
  dvlrdrb = vlro * dswthc - swth * (Six * c6000 * (denom1**2) * (rbig**5) + Eight * c8000 * denom2 * denom2a * (rbig**3) + Five * qq * (denom3**2) * (rbig**4))
  if (iprt .ne. 0) then
    !write(*,*) '**** vn2n2d ** dvlrdr1 ', dvlrdr1, swth, dqqdr1, denom3
    !write(*,*) '**** vn2n2d ** dvlrdr2 ', dvlrdr2, dqqdr2, dc6000r2, denom1, dc8000r2, denom2
  end if
  dvlrdr1 = dvlrdr1 + swth * (dqqdr1 * denom3 + dc6000r1 * denom1 + dc8000r1 * denom2)                                            
  dvlrdr2 = dvlrdr2 + swth * (dqqdr2 * denom3 + dc6000r2 * denom1 + dc8000r2 * denom2)                                            
  dvlrph  = denom3  * dqqdph * swth  

  r1re  = One / (r1 + re)
  r2re  = One / (r2 + re)
  r1re2 = r1re * r1re
  r2re2 = r2re * r2re
  ddr1  =(r1 - 2.1d0) * r1re
  ddr2  =(r2 - 2.1d0) * r2re
  dr12  = ddr1 * ddr2
  dr1p2 = ddr1 + ddr2

  bra = pairwise * swth
  brb = bra * rbig
  if (iprt .ne. 0) then
    !write(*,*) '**** vn2n2d ** pairwise, bra, brb ', pairwise, bra, brb
  end if

  dpairdr3 = beta70 * exp(beta * r3)
  dpairdr4 = beta70 * exp(beta * r4)
  dpairdr5 = beta70 * exp(beta * r5)
  dpairdr6 = beta70 * exp(beta * r6)

  dbradr1 = pairwise * dswthab                                      
  dbradr2 = pairwise * dswthab                                      
  dbradr3 = swth * dpairdr3
  dbradr4 = swth * dpairdr4
  dbradr5 = swth * dpairdr5
  dbradr6 = swth * dpairdr6

  !dbrbdrb = bra
  dbrbdr1 = rbig * dbradr1
  dbrbdr2 = rbig * dbradr2
  dbrbdr3 = rbig * dbradr3
  dbrbdr4 = rbig * dbradr4
  dbrbdr5 = rbig * dbradr5
  dbrbdr6 = rbig * dbradr6

  func1a = (coef(2) + coef(7)*y2 + coef(13)*y4 + coef(16)*y5)
  func2a = (coef(4) + coef(9)*y2)  

  func1  = coef(1)  + coef(6) * y2 + coef(10) * y3 + coef(12) * y4 + coef(15) * y5 + coef(18) * y6 + dr1p2 * func1a
  func2  = coef(3)  + coef(8) * y2 + coef(11) * y3 + coef(14) * y4 + coef(17) * y5 + coef(19) * y6 + dr1p2 * func2a + dr12 * coef(5)

  fit    = bra * func1 + brb * func2 + pairwise + vlr + (v12 + v34) + coef(20)
  if (iprt .gt. 0) then
    !write(*,*) '**** vn2n2d ** fit: ', fit, bra, func1, func2, pairwise, vlr, v12, v34, coef(20)
  end if

  dfunc1dr1 = func1a * re2 * r1re2
  dfunc1dr2 = func1a * re2 * r2re2
  dfunc1dt1 = -Two * cst1 * (coef(6) + coef(10) * y2b - coef(12) * z2b * cp + Two * coef(15) * y2a - coef(18) * cp + dr1p2 * (coef(7) - coef(13) * z2b * cp + Two * coef(16) * y2a))
  dfunc1dt2 = -Two * cst2 * (coef(6) + coef(10) * y2a - coef(12) * z2a * cp + Two * coef(15) * y2b - coef(18) * cp + dr1p2 * (coef(7) - coef(13) * z2a * cp + Two * coef(16) * y2b))
  dfunc1dph = -Two * sp * ( (coef(12) + coef(13) * dr1p2) * z2a * z2b + coef(18) * (z2a + z2b))
  dfunc2dr1 = (re2 * r1re2) * (func2a + coef(5) * ddr2)
  dfunc2dr2 = (re2 * r2re2) * (func2a + coef(5) * ddr1)
  dfunc2dt1 = -Two * cst1 * (coef(8) + coef(11) * y2b - coef(14) * z2b * cp + Two * coef(17) * y2a - coef(19) * cp + dr1p2 * coef(9))
  dfunc2dt2 = -Two * cst2 * (coef(8) + coef(11) * y2a - coef(14) * z2a * cp + Two * coef(17) * y2b - coef(19) * cp + dr1p2 * coef(9))
  dfunc2dph = -Two * sp * (coef(14) * z2a * z2b + coef(19) * (z2a + z2b))
 
  dfitr1 = dbradr1 * func1 + bra * dfunc1dr1 + dbrbdr1 * func2 + brb * dfunc2dr1 + dvlrdr1 + dv12dr1 
  dfitr2 = dbradr2 * func1 + bra * dfunc1dr2 + dbrbdr2 * func2 + brb * dfunc2dr2 + dvlrdr2 + dv34dr2                                                              
  dfitr3 = dbradr3 * func1 + dbrbdr3 * func2 + dpairdr3
  dfitr4 = dbradr4 * func1 + dbrbdr4 * func2 + dpairdr4
  dfitr5 = dbradr5 * func1 + dbrbdr5 * func2 + dpairdr5
  dfitr6 = dbradr6 * func1 + dbrbdr6 * func2 + dpairdr6
  dfitt1 = bra * dfunc1dt1 + brb * dfunc2dt1 + dvlrdt1 
  dfitt2 = bra * dfunc1dt2 + brb * dfunc2dt2 + dvlrdt2 
  dfitph = bra * dfunc1dph + brb * dfunc2dph + dvlrph  
  dfitrb = dvlrdrb + pairwise * dswthc * (func1 + rbig * func2) + bra * func2 
  if (iprt .gt. 0) then
    !write(*,*) '**** vn2n2d ** dfitr1: ', dfitr1, dbradr1, func1, bra, dfunc1dr1, dbrbdr1, func2, brb, dfunc2dr1, dvlrdr1, dv12dr1
    !write(*,*) '**** vn2n2d ** dfitr2: ', dfitr2, dbradr2, func1, bra, dfunc1dr2, dbrbdr2, func2, brb, dfunc2dr2, dvlrdr2, dv34dr2
    !write(*,*) '**** vn2n2d ** dfitt1: ', dfitt1
    !write(*,*) '**** vn2n2d ** dfitt2: ', dfitt2
    !write(*,*) '**** vn2n2d ** dfitph: ', dfitph
    !write(*,*) '**** vn2n2d ** dfitrb: ', dfitrb
    !write(*,*) '**** vn2n2d ** dfitr3: ', dfitr3
    !write(*,*) '**** vn2n2d ** dfitr4: ', dfitr4
    !write(*,*) '**** vn2n2d ** dfitr5: ', dfitr5
    !write(*,*) '**** vn2n2d ** dfitr6: ', dfitr6
  end if
  
  do i=1,12
    dfit(i) = dfitr1 * dr1(i)     + & 
              dfitr2 * dr2(i)     + &
              dfitr3 * dr3(i)     + &
              dfitr4 * dr4(i)     + &
              dfitr5 * dr5(i)     + &
              dfitr6 * dr6(i)     + &
              dfitt1 * dtheta1(i) + &
              dfitt2 * dtheta2(i) + &
              dfitph * dphi(i)    + &
              dfitrb * drbig(i)
  end do

  return
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


!________________________________________________________________________________________________________________________________!
Subroutine diatccsdd(r, ans, dansdr, iwant)
! iwant=1 for pes
! iwant=2 for qm
! iwant=3 for dipp

  real(rkp)                               ,intent(in)     :: r
  real(rkp)                               ,intent(out)    :: ans
  real(rkp)                               ,intent(out)    :: dansdr
  integer                                 ,intent(in)     :: iwant

  real(rkp)                                               :: r2, r4, r6
  real(rkp)                                               :: d2, d4, d6
  real(rkp)                                               :: c8, c10
  real(rkp)                                               :: vlr
  real(rkp)                                               :: vrep
  real(rkp)                                               :: dvlrdr
  real(rkp)                                               :: dvrepdr
  real(rkp)                                               :: dr
  real(rkp)                                               :: ex
  real(rkp)                                               :: dexdr
  real(rkp)                                               :: poly, dpoly
  real(rkp)                                               :: sum, dsum
  integer                                                 :: i

  if (iwant .eq. 1) then

    ! csd diatomic potential
    r2 = r  * r
    r4 = r2 * r2
    r6 = r2 * r4
    
    d2 = d  * d
    d4 = d2 * d2
    d6 = d2 * d4
    
    c8  = c6 * 20.d0
    c10 = c6 * 500.d0
    vlr = einf + c6 / (r6+d6) + c8 / (r4+d4)**2 + c10 / (r2+d2)**5
    
    vrep    = co * exp( -ao * r ) / r
    dvlrdr  = -Six * c6 * (r4*r) / (r6+d6)**2 - Eight * c8 * (r2*r) / (r4+d4)**3 - Ten * c10 * r / (r2+d2)**6
    dvrepdr = -vrep * ao - vrep / r
    dr    = r-re
    ex    = exp( -a1 * r ) * r**npow
    dexdr = -a1 * ex + ex * npow / r
    poly  = pes(idp)
    do i=idp-1,1,-1
      poly = pes(i) + poly * dr
    end do

    dpoly  = (((((((( Nine  * pes(10) * dr + &
                      Eight * pes(9)) * dr + & 
                      Seven * pes(8)) * dr + &
                      Six   * pes(7)) * dr + &
                      Five  * pes(6)) * dr + &
                      Four  * pes(5)) * dr + &
                      Three * pes(4)) * dr + &
                      Two   * pes(3)) * dr + &
                              pes(2))
    ans    = vrep + poly * ex + vlr
    dansdr = dvrepdr + dvlrdr + poly * dexdr + ex * dpoly

  else if (iwant .eq. 2) then
    ! quadrupole
    ex  = exp( -a2 * r )
    dr  = r-re
    sum = qm(idq)
    do i=idq-1,1,-1
      sum = qm(i) + sum * dr
    end do
    ans    = sum * ex * r
    dansdr = sum * ex - a2 * ans + r * ex * ( ( (Four * qm(5) * dr + Three * qm(4) ) * dr + Two * qm(3) ) * dr + qm(2) )

  else if (iwant .eq. 3)then
    ! dipole polarizability
    ex   = exp( -a3 * r )
    dr   = r-re
    sum  = dip(idd)
    dsum = Zero                                                      
    do i=idd-1,2,-1
      dsum =    sum + dsum * dr                                                
      sum  = dip(i)  + sum * dr
    end do
    ans    = sum * ex * r + dip(1)
    dansdr = ex * r * dsum + sum * ex - a3 * r * ex * sum      

  else
    !write(*,*) '**** diatccsdd * unknown iwant: ', iwant
    stop
  end if

  return
End Subroutine
!--------------------------------------------------------------------------------------------------------------------------------!


End Module